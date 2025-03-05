#SPectral ANalysis software (SPAN)
#Written by Daniele Gasparri#
#SPectral ANalysis software (SPAN).
#Written by Daniele Gasparri#

"""
    Copyright (C) 2020-2025, Daniele Gasparri

    E-mail: daniele.gasparri@gmail.com

    SPAN is a GUI interface that allows to modify and analyse 1D astronomical spectra.

    1. This software is licensed **for non-commercial use only**.
    2. The source code may be **freely redistributed**, but this license notice must always be included.
    3. Any user who redistributes or uses this software **must properly attribute the original author**.
    4. The source code **may be modified** for non-commercial purposes, but any modifications must be clearly documented.
    5. **Commercial use is strictly prohibited** without prior written permission from the author.

    DISCLAIMER:
    THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES, OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT, OR OTHERWISE, ARISING FROM, OUT OF, OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""

#******************************************************************************************
#******************************************************************************************
#*************************** SPECTRA ANALYSIS FUNCTIONS FOR SPAN **************************
#******************************************************************************************
#******************************************************************************************


try:#Local imports
    from span_functions import spec_manipul as spman
    from span_functions import system_span as stm
    from span_functions import utilities as uti
    from span_functions import build_templates as template

except ModuleNotFoundError: #local import if executed as package
    from . import spec_manipul as spman
    from . import system_span as stm
    from . import utilities as uti
    from . import build_templates as template

#pPXF import
from ppxf.ppxf import ppxf
import ppxf.ppxf_util as util
import ppxf.sps_util as lib
from urllib import request
from pathlib import Path

#Python imports
import numpy as np
import math as mt
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.gridspec as gridspec

from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.constants import h,k,c
from scipy.integrate import quad
from scipy.optimize import leastsq
import scipy.stats
from scipy.interpolate import griddata
from scipy.interpolate import CubicSpline
from scipy.signal import savgol_filter

from time import perf_counter as clock
from os import path
import os
import joblib

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C, Matern
from sklearn.model_selection import GridSearchCV

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(CURRENT_DIR)

# 1) BLACK BODY FITTING FUNCTION. Adapted from: https://github.com/Professor-G/BlackbodyFit
def blackbody(wavelength, T):

    """
    Planck's law, which describes the black body radiation
    of a source in thermal equilibrium at a given temperature T.
    """

    return 2*h*c**2 / (wavelength**5 * (np.e**(h*c / (wavelength*k*T)) - 1))


#*****************************************************************************************************
# 2) blackbody fitting
def blackbody_fit(wavelength, flux, initial_wave, final_wave, t_guess, with_plots, save_plot, result_plot_dir, spec_name):

    """
    This function fits the Planck function to a spectrum with the user defined wavelength
    range and initial temperature guess
    Input: wavelength and flux arrays of the spectrum, float initial and final
    wavelength range to fit, int temperature initial guess, bool show or not
    the plot of the fit
    Output: effective temperature of the fitted model, array of the flux residuals
            of the fitted wavelength range.

    """

    #extracting the wavelength and the flux within the selected band
    wave_bb = wavelength[(wavelength >= initial_wave) & (wavelength <= final_wave)]
    flux_bb = flux[(wavelength >= initial_wave) & (wavelength <= final_wave)]

    wavelength = wave_bb
    flux = flux_bb

    #normalising the flux to the median wavelength
    median_wave = np.median(wavelength)
    epsilon_wave = 5.

    norm_flux = spman.norm_spec(wavelength, flux, median_wave, epsilon_wave, flux)
    flux = norm_flux

    wavelength = wavelength*1e-9 #convert to meters for SI consistency

    def blackbody_flux(wave, flux, T):
        blackbody_flux = blackbody(wave, T)

        flux_integral = np.abs(np.trapz(flux, wave))
        planck_integral = np.abs(quad(blackbody, np.min(wave), np.max(wave), args = T ))[0]
        scale_factor = flux_integral / planck_integral

        return scale_factor*blackbody_flux

    def residuals(T, y, lam):
        return y - blackbody_flux(wavelength, flux, T)

    t0 = np.array([t_guess]) #the initial temperature guess for the optimization
    T = leastsq(residuals, t0, args = (flux, wavelength))[0].astype('float32')
    bbody_fit = blackbody(wavelength, T)

    instrumental_integral = np.abs(np.trapz(flux, wavelength))
    planck_integral = np.abs(quad(blackbody, np.min(wavelength), np.max(wavelength), args = T))[0]
    scale_factor = instrumental_integral / planck_integral
    y = scale_factor*bbody_fit
    wavelength = wavelength*1e9
    temperature = int(round(T[0]))
    residual_bb = flux-y

    if with_plots or save_plot:
        fig = plt.figure()
        ax = fig.add_subplot(211)
        plt.plot(wavelength, flux, label = "DATA")
        plt.plot(wavelength, y, 'r-',  label = "FIT")
        plt.title("Blackbody fitting, T = " + str(temperature) +' K' )
        plt.xlabel(r"$\lambda$ (nm)")
        plt.ylabel("Relative flux")
        plt.legend(loc = 1)
        plt.tight_layout()

        ax2 = fig.add_subplot(212)
        plt.title("Residuals" )
        plt.ylabel("Relative flux")
        plt.xlabel(r"$\lambda$ (nm)")
        plt.plot(wavelength, residual_bb, 'g.', label = "residuals")
        plt.tight_layout()

        if with_plots:
            plt.show()
        else:
            plt.savefig(result_plot_dir + '/'+ 'blackbody_'+ spec_name + '.png', format='png', dpi=300)
        plt.close()

    return temperature, residual_bb


#*****************************************************************************************************
# 3) Cross-correlation, adapted from  pyasl, but working
def crosscorrRV(w, f, tw, tf, rvmin, rvmax, drv, skipedge=0, edgeTapering=None):

    """
    This function calculates the cross-correlation function between the input spectrum
    and the provided template spectrum, within an interval of velocity dispersion, with a
    delta v step of velocity dispersion.
    Input: wavelength and flux arrays of the spectrum,wavelength and flux arrays of the
           template, float minimum and maximum radial velocity interval to compute the
           cross-correlation, float step of velocity dispersion.
    Output: array of the velocities, array containing the cross-correlation function

    """

  # Copy and cut wavelength and flux arrays
    w, f = w.copy(), f.copy()
    if skipedge > 0:
        w, f = w[skipedge:-skipedge], f[skipedge:-skipedge]

    if edgeTapering is not None:
        # Smooth the edges using a sine
        if isinstance(edgeTapering, float):
            edgeTapering = [edgeTapering, edgeTapering]

        indi = np.where(w < w[0]+edgeTapering[0])[0]
        f[indi] *= np.sin((w[indi] - w[0])/edgeTapering[0]*np.pi/2.0)
        # Carry out edge tapering (right edge)
        indi = np.where(w > (w[-1]-edgeTapering[1]))[0]
        f[indi] *= np.sin((w[indi] - w[indi[0]])/edgeTapering[1]*np.pi/2.0 + np.pi/2.0)

    # Speed of light in km/s
    c = 299792.458

    # Calculate the cross correlation
    drvs = np.arange(rvmin, rvmax, drv)
    cc = np.zeros(len(drvs))
    for i, rv in enumerate(drvs):

        # Apply the Doppler shift
        fi = interpolate.interp1d(tw*(1.0 + rv/c), tf)
        # Shifted template evaluated at location of spectrum
        cc[i] = np.sum(f * fi(w))

    return drvs, cc




# 3B) Considering the cosmologica redshift instead
def crosscorrZ(w, f, tw, tf, zmin, zmax, dz, skipedge=0, edgeTapering=None):

    """
    Calculate the cross-correlation function between the input spectrum
    and the provided template over a range of cosmological redshifts, with a step dz.
    Input:
        - w, f: wavelength and flux arrays of the spectrum.
        - tw, tf: wavelength and flux arrays of the template.
        - zmin, zmax: minimum and maximum redshift interval to compute the cross-correlation.
        - dz: redshift step.
    Output:
        - array of redshifts, array containing the cross-correlation function.
    """

    # Copy and cut wavelength and flux arrays
    w, f = w.copy(), f.copy()
    if skipedge > 0:
        w, f = w[skipedge:-skipedge], f[skipedge:-skipedge]

    # Edge tapering to smooth boundaries
    if edgeTapering is not None:
        if isinstance(edgeTapering, float):
            edgeTapering = [edgeTapering, edgeTapering]

        indi = np.where(w < w[0]+edgeTapering[0])[0]
        f[indi] *= np.sin((w[indi] - w[0])/edgeTapering[0]*np.pi/2.0)
        indi = np.where(w > (w[-1]-edgeTapering[1]))[0]
        f[indi] *= np.sin((w[indi] - w[indi[0]])/edgeTapering[1]*np.pi/2.0 + np.pi/2.0)

    # Calculate the cross correlation
    redshifts = np.arange(zmin, zmax, dz)
    cc = np.zeros(len(redshifts))
    for i, z in enumerate(redshifts):
        # Apply the cosmological redshift
        shifted_template_wavelengths = tw * (1 + z)
        # Use cubic spline interpolation
        fi = CubicSpline(shifted_template_wavelengths, tf, extrapolate=False)
        # fi = interpolate.interp1d(shifted_template_wavelengths, tf)
        # Evaluate cross-correlation
        cc[i] = np.sum(f * fi(w))

    return redshifts, cc




#*****************************************************************************************************
# 4) Apply Cross-correlation to span standard, use with cautions for galaxies. Works well for stellar spectra
def crosscorr (wavelength_spec, flux_spec, template, lambda_units_template, wave_interval, smooth_vel, vel_z_interval, is_vel):

    """
    This function applies the cross-correlation function to the selected spectrum
    with the selected template.
    Input: wavelength and flux arrays of the spectrum, path and name of the template
           string wavelength units of the templae (physical: 'nm', 'A', 'mu')
           Array of the wavelength interval to cross-correlate ([min_wave, max_wave]
           bool smoothing the template by adding a velocity dispersion broadening, int velocity
           dispersion value (sigma) to add to the template.
    Output: arrays or computed radial velocities and cross-correlation function, float
            radial velocity value corresponding at the maximum of the cross-correlation
            function, float maximum of the cross-correlation function, wavelength and flux arrays
            of the spectrum within the selected wavelength range, wavelength and flux arrays of the
            template within the selected wavelength range.

    """

    #reading the template
    wavelength_template_orig, flux_template_orig, step, name = stm.read_spec(template, lambda_units_template)

    #Sampling to the common, smallest, delta lambda
    step_spec = wavelength_spec[1]-wavelength_spec[0]
    step_template = step

    if step_spec > step_template:
        wavelength, flux, npoint_resampled = spman.resample(wavelength_spec, flux_spec, step_template)

        #resampling to linear also the template
        wavelength_template, flux_template, npoint_template = spman.resample(wavelength_template_orig, flux_template_orig, step_template)

    elif step_spec <= step_template:
        wavelength_template, flux_template, npoint_template = spman.resample(wavelength_template_orig, flux_template_orig, step_spec)
        #resampling to linear also the spec
        wavelength, flux, npoint_resampled = spman.resample(wavelength_spec, flux_spec, step_spec)


    if is_vel:
    #Extracting the arrays values, considering that they can be inverted, so I look for max and min values)
        low_lim_vel =  np.min(vel_z_interval)
        high_lim_vel = np.max(vel_z_interval)
        low_wave_interval = np.min(wave_interval)
        high_wave_interval = np.max(wave_interval)
        step_vel = 2.

    if not is_vel:
        low_lim_z =  np.min(vel_z_interval)
        high_lim_z = np.max(vel_z_interval)
        low_wave_interval = np.min(wave_interval)
        high_wave_interval = np.max(wave_interval)

        #putting some constraints on the z step as a function of the redshift interval to consider
        if high_lim_z - low_lim_z < 1:
            step_z = 0.0001
        if high_lim_z - low_lim_z >= 1 and high_lim_z - low_lim_z < 2:
            step_z = 0.001
        if high_lim_z - low_lim_z >= 2:
            step_z = 0.01


    #smoothing the template
    if(smooth_vel > 0):
        flux_gauss_template = spman.sigma_broad(wavelength_template, flux_template, smooth_vel)
        flux_template = flux_gauss_template


    #if the interval is zero:
    if low_wave_interval == 0:

        # Normalise the spectrum and the templates
        flux = flux / np.median(flux)
        flux_template = flux_template / np.median(flux_template)

        #continuum subtraction
        poly_degree = 3
        math_operation = 'divide'
        flux, continuum_model = spman.continuum(wavelength, flux, False, None, poly_degree, math_operation, False)
        flux_template, continuum_model = spman.continuum(wavelength_template, flux_template, False, None, poly_degree, math_operation, False)

        if is_vel:

            rv, cc = crosscorrRV(wavelength, flux, wavelength_template, flux_template, low_lim_vel, high_lim_vel, step_vel)

            #estrapolating the most probable rv value
            max_corr_fcn = np.argmax(cc)
            rv_at_max = rv[max_corr_fcn]
            cc_at_max = cc[max_corr_fcn]

            return rv, cc, rv_at_max, cc_at_max, wavelength, flux, wavelength_template, flux_template
        if not is_vel:
            z, cc = crosscorrZ(wavelength, flux, wavelength_template, flux_template, low_lim_z, high_lim_z, step_z)

            #estrapolating the most probable rv value
            max_corr_fcn = np.argmax(cc)
            z_at_max = z[max_corr_fcn]
            cc_at_max = cc[max_corr_fcn]
            return z, cc, z_at_max, cc_at_max, wavelength, flux, wavelength_template, flux_template


    #otherwise I need to extract lambda and flux both for the spectrum and the template.
    else:
        corr_wave = []
        corr_flux = []

        npoints_wave = len(wavelength)
        for i in range(npoints_wave):
            if wavelength[i]>= low_wave_interval and wavelength[i]<= high_wave_interval:
                corr_wave.append(wavelength[i])
                corr_flux.append(flux[i])

        #converting to numpy
        corr_wave_np = np.asarray(corr_wave, dtype = float)
        corr_flux_np = np.asarray(corr_flux, dtype = float)
        wavelength = corr_wave_np
        flux = corr_flux


        # Normalise the spectrum and the templates
        flux = flux / np.median(flux)
        flux_template = flux_template / np.median(flux_template)

        #continuum subtraction
        poly_degree = 3
        math_operation = 'divide'
        flux, continuum_model = spman.continuum(wavelength, flux, False, None, poly_degree, math_operation, False)
        flux_template, continuum_model = spman.continuum(wavelength_template, flux_template, False, None, poly_degree, math_operation, False)


        if is_vel:
            rv, cc = crosscorrRV(wavelength, flux, wavelength_template, flux_template, low_lim_vel, high_lim_vel, step_vel)

            #estrapolating the most probable rv value
            max_corr_fcn = np.argmax(cc)
            rv_at_max = rv[max_corr_fcn]
            cc_at_max = cc[max_corr_fcn]

            return rv, cc, rv_at_max, cc_at_max, wavelength, flux, wavelength_template, flux_template

        else:
            z, cc = crosscorrZ(wavelength, flux, wavelength_template, flux_template, low_lim_z, high_lim_z, step_z)

            #estrapolating the most probable rv value
            max_corr_fcn = np.argmax(cc)
            z_at_max = z[max_corr_fcn]
            cc_at_max = cc[max_corr_fcn]

            return z, cc, z_at_max, cc_at_max, wavelength, flux, wavelength_template, flux_template



#*************************************************************************************************
# 5) Velocity dispersion measurement
def sigma_measurement(wavelength, flux, spec_test_template, lambda_units_template, resolutionR_spec, resolution_template, banda1, banda1_cont, err_calc):
    c = 299792.458

    """
    This function performs a least square fit of the selected spectrum with a
    template in order to retrieve a fast estimation of the velocity dispersion.
    The results are a good estimate of the velocity dispersion, but for science
    results is is advisable to use the Kinematics tasks that perform the pPXF
    algorithm.
    Input: wavelength and flux arrays of the spectrum, path and name of the template,
            string wavelength units of the template (physical: 'nm', 'A', 'mu'), int spectral
            resolution of the template (R), spectral resolution of the template (0 for EMILES
            templates, SPAN will take care of the different resolution between optical and NIR bands),
            array containing the wavelength range of the spectral band to fit ([min_wave, max_wave])
            array containing the wavelength range of the continuum level within the band
            ([min_wave_continuum, max_wave_continuum]), bool value whether to calculate the
            uncertainties (MonteCarlo simulations) of not.
    Output: float velocity dispersion value, float uncertainty (if calculated), float
            chi2 of the fit, wavelength and flux arrays of the spectrum and the template
            within the fitted band, float mean velscale of the spectrum within the band.

    """

        #************************************ select the bands/lines ******************************
    wave_range = banda1
    snr_range = banda1_cont
    line_name = 'band1'
    band_cont = banda1_cont
    wave_norm = np.mean(banda1_cont)

    sigma_instrumental = (c/resolutionR_spec)/2.355
    resolution_lambda_fwhm_spec = (np.mean(wave_range)/resolutionR_spec)

    #if the resolution_template variable is zero, I assume they are (E)MILES spectra, otherwise I transform the resolution from R to sigma.
    if resolution_template != 0:
        resolution_lambda_fwhm_temp = (np.mean(wave_range)/resolution_template)
        resolution_temp = (c/resolution_template)/2.355

    #preparing variables for the determination of sigma
    initial_sigma = 0.
    final_sigma = 360. # Maximum meaningful sigma
    step_sigma = 2. #sigma step (km/s)
    number_values_sigma = round(final_sigma/step_sigma+1)
    sigma_values = np.zeros(number_values_sigma)

    #filling the sigma vector
    for i in range(1, number_values_sigma):
        sigma_values[i] = sigma_values[i-1]+step_sigma

    chisquare_fit = [] #vector for the fit chi squared values



    #********************************* Let's rock! *********************************

    #read the template
    wavelength_template, flux_template, step_template, name = stm.read_spec(spec_test_template, lambda_units_template)


    #rough continuum subtraction for more stable results
    flux, cont_flux = spman.sub_cont(wavelength, flux, 'divide')
    flux_template, cont_flux = spman.sub_cont(wavelength_template, flux_template, 'divide')

    #rebin to smallest # better to rebin all to a constant step, in case it's not
    optimal_step = 2*np.mean(wave_range)*step_sigma/c # I tested that this sampling is the optimal one: smaller values don't change the result, while greater values start to feel the quantization effect.

    wavelength_template, flux_template, points = spman.resample(wavelength_template, flux_template, optimal_step)
    wavelength, flux, point_spec = spman.resample(wavelength, flux, optimal_step)

    #uniform the wavelength grid
    interpfunc = interpolate.interp1d(wavelength_template, flux_template, kind = 'linear', fill_value='extrapolate')
    flux_template = (interpfunc(wavelength))

    #storing the original template flux
    flux_template_original = flux_template

    #extract the line flux and wavelength arrays
    line_wave = wavelength[(wavelength >= wave_range[0]) & (wavelength <= wave_range[1])]
    line_flux_spec = flux[(wavelength >= wave_range[0]) & (wavelength <= wave_range[1])]

    #resolution for the EMILES models: if lambda < 895, the resolution is constant with wavelength, which means that the resolution in sigma diminishes with the increasing wavelength.
    if resolution_template == 0:
        if np.mean(line_wave) < 895.:
            resolution_lambda_fwhm_temp = 0.251 #fwhm
            #converting to sigma
            resolution_temp = resolution_lambda_fwhm_temp/np.mean(line_wave)*c/2.3548 # in sigma velocity. it's an approximation, since the resolution in sigma changes with the wavelength, but if the band is small, the approximation is good.
            print (resolution_temp)
        else:
            resolution_vel_fwhm_temp = 60. #fwhm, in km/s
            resolution_lambda_fwhm_temp = (resolution_vel_fwhm_temp*np.mean(wave_range))/c
            resolution_temp = resolution_vel_fwhm_temp/2.3548
            print (resolution_temp)


    #normalise the spectrum
    epsilon_norm = optimal_step*200
    line_flux_spec_norm = spman.norm_spec(line_wave, line_flux_spec, wave_norm, epsilon_norm, line_flux_spec)

    #calculating the SNR
    snr_flux = line_flux_spec_norm[(line_wave >= snr_range[0]) & (line_wave <= snr_range[1])]
    snr = np.mean(snr_flux)/np.std(snr_flux)


    #TEST FITTING TAMPLATE
    #The idea is: select the template, broad to 200 km/s, normalise, extract the working band, extract the working continuum, fitting, finding the displacement with respect to the spectrum, then apply to the original, not broadened and normalized template.

    #preparing the template
    #broadening the template to 200 km/s
    flux_template_broad = spman.sigma_broad(wavelength, flux_template_original, 200.)

    #extract the line flux template
    line_flux_template_orig = flux_template_broad[(wavelength >= wave_range[0]) & (wavelength <= wave_range[1])]

    #normalize the template in the wavelength range
    line_flux_temp_test = spman.norm_spec(line_wave, line_flux_template_orig, wave_norm, epsilon_norm, line_flux_template_orig)

    #normalize the whole, original template
    flux_temp_test = spman.norm_spec(wavelength, flux_template_original, wave_norm, epsilon_norm, flux_template_original)

    #extract flux and wavelength of the continuum
    cont_spec_banda = line_flux_spec_norm[(line_wave >= band_cont[0]) & (line_wave <= band_cont[1])]
    cont_temp_banda = line_flux_temp_test[(line_wave >= band_cont[0]) & (line_wave <= band_cont[1])]
    wave_banda = line_wave[(line_wave >= band_cont[0]) & (line_wave <= band_cont[1])]

    #calculate the rms of the banda1
    rms_banda = np.std(cont_spec_banda)
    #the shift_step will be half or one rms
    shift_step = rms_banda/2.

    starting_flux_level = 0.95
    delta_starting_flux = 1-starting_flux_level
    cont_temp_banda = cont_temp_banda - delta_starting_flux
    end_flux_level = 1.1
    actual_level = starting_flux_level
    chisquare_test = []
    shift_array = []

    #fitting procedure
    while actual_level <= end_flux_level:
        #chi square
        chisqr_test = 0. #it is a chisquared
        for i in range(len(wave_banda)):
            chisqr_test = chisqr_test + (cont_spec_banda[i]-cont_temp_banda[i])**2/cont_temp_banda[i]

        chisquare_test.append(chisqr_test) #chi squared vector!
        shift_array.append(actual_level)
    #test to stop the cycle when reached the minimum chi square, without exploring other values
        actual_level = actual_level + shift_step
        cont_temp_banda = cont_temp_banda + shift_step

    min_index_test = np.argmin(chisquare_test)
    best_fit_level_temp = shift_array[min_index_test]
    print ('Adjusting the template to new level:', best_fit_level_temp)

    #applying the corrections
    flux_template_shifted = flux_temp_test + (best_fit_level_temp-1)
    line_flux_temp_test_shifted = line_flux_temp_test + (best_fit_level_temp-1)

    print ('line selected: ', line_name)
    print ('SNR line:', round(snr))
    print('Processing...')


    #********************** FITTING PROCEDURE, BOTH SPECTRAL REGION AND LINES **************************
    shape = len(line_wave)
    line_flux_template_norm = np.zeros(shape)

    for j in range(number_values_sigma):

        previous_line_flux_template_norm = line_flux_template_norm

        #broadening the template
        flux_template_broad = spman.sigma_broad(wavelength, flux_template_shifted, sigma_values[j])

        #extract the wavelength range for the template. I use it only to know the initial guesses
        line_flux_template_orig = flux_template_broad[(wavelength >= wave_range[0]) & (wavelength <= wave_range[1])]

        #normalise the template
        line_flux_template_norm = line_flux_template_orig


        #resudials
        chisqr = 0. #it is a chisquared
        for i in range(len(line_wave)):
            chisqr = chisqr + (line_flux_spec_norm[i]-line_flux_template_norm[i])**2/line_flux_template_norm[i]

        chisquare_fit.append(chisqr) #chi squared vector!

        #test to stop the cycle when reached the minimum chi square, without exploring other values
        min_index = np.argmin(chisquare_fit)
        if j > 0:
            if min_index < j:
                line_flux_template_norm = previous_line_flux_template_norm # if I enter this cycle, I use the previous stored value of the broadened template because its the best
                break

    if (j == 1): # that means that I found the minumum at 0, with no broadening of the template
        print ('Warning: the resolution of the template is greater than the sigma you want to measure. The value obtained will be just a superior limit')

    #finding the minimum chi square and the respective value of sigma
    min_value = np.argmin(chisquare_fit)

    sigma_vel = sigma_values[min_value] # this is the holy grail: the most probable velocity dispersion value
    min_residual = chisquare_fit[min_value] # this is the chi square of the best fit


    #***************************** Uncertainties with MonteCarlo simulations ************
    # I perturb the line_flux_temp_fit_norm by adding noise equal to the SNR of my line, then I do the Gaussian fit with the same non-noisy template and see how the sigma that I obtain from the Gaussian fit fluctuates.
    if err_calc:
        scale_noise = np.mean(snr_flux)/snr
        number_noisy_cont = 20 #how many syntethics? a lot, but you need to consider also computation time

        sigma_vel_err_tot_ = []

        print ('Calculating the error...')
        for k in range (number_noisy_cont):
            #generate the noisy line
            noise_array = np.random.standard_normal((len(line_wave),))
            noise_array_scaled = noise_array * scale_noise
            noisy_template = []
            for i in range(len(line_wave)):
                noisy_template.append(line_flux_template_norm[i] + noise_array_scaled[i])

            #maximum error on the sigma estimation
            max_error_sigma = int(round(700./snr)) #empirical value. Checked and seems ok
            if max_error_sigma > 70 and k == 0:
                max_error_sigma = 70
                print ('Warning: SNR < 10, large error and long processing time!')
            step_error = 2.
            initial_sigma_err = sigma_vel - max_error_sigma

            #if I reach negative values in the interval of possible sigma values:
            if initial_sigma_err < 0:
                initial_sigma_err = 0
            final_sigma_err = sigma_vel + max_error_sigma

            number_values_sigma_err = int(round((final_sigma_err-initial_sigma_err)/step_error+1))
            #filling the sigma vector
            sigma_values_err = []
            sigma_value_err = initial_sigma_err

            for i in range(1, number_values_sigma_err):
                sigma_values_err.append(sigma_value_err)
                sigma_value_err =  sigma_value_err + step_error

            residuals_err = []

            # fitting procedure
            for h in range(len(sigma_values_err)):
                #broadening the template
                flux_template_new = spman.sigma_broad(wavelength, flux_template_original, sigma_values_err[h])

                #extract the wavelength range for the template. I use it only to know the initial guesses
                line_flux_template_new = flux_template_new[(wavelength >= wave_range[0]) & (wavelength <= wave_range[1])]

                #normalize the template
                line_flux_template_new_norm = spman.norm_spec(line_wave, line_flux_template_new, wave_norm, epsilon_norm, line_flux_template_new)

                #resudials
                residual_err = 0. #it is a chisquared
                for i in range(len(line_wave)):
                    residual_err = residual_err + (noisy_template[i]-line_flux_template_new_norm[i])**2/line_flux_template_new_norm[i]

                residuals_err.append(residual_err) #chi squared vector!

                min_index = np.argmin(residuals_err)
                if h > 1:
                    if min_index < h:
                        break

            #finding the minimum
            min_value_err = np.argmin(residuals_err)
            sigma_vel_err = sigma_values_err[min_value_err]
            sigma_vel_err_tot_.append(sigma_vel_err)

        #storing the data
        sigma_vel_err_tot = np.asarray(sigma_vel_err_tot_, dtype = float)
        error_sigma_fit = np.std(sigma_vel_err_tot)

        #the total error is the quadratic sum of the error above + the quantization error due to the step of sigma values selected.
        total_error = mt.sqrt(error_sigma_fit**2+step_error**2 +step_sigma**2)
    else:
        total_error = 0.

    ####################################################################################################

    sigma_real_broadened = mt.sqrt(sigma_vel**2+resolution_temp**2)

    if (sigma_real_broadened == resolution_temp):
        sigma_real = sigma_real_broadened
    elif (sigma_real_broadened < sigma_instrumental):
        sigma_real = sigma_instrumental
        print ('WARNING: The real velocity dispersion is lower than the instrumental sigma of the spectrum. Do not trust the result!')
    else:
        sigma_real = mt.sqrt(sigma_real_broadened**2 - sigma_instrumental**2)

    print ('Resolution template in A (FWHM): ', resolution_lambda_fwhm_temp*10)
    print ('Resolution spectrum in A (FWHM): ', resolution_lambda_fwhm_spec*10)
    print ('Resolution sigma template (km/s): ', resolution_temp)
    print ('Resolution sigma spectrum: (km/s)', sigma_instrumental)
    print ('Best sigma gaussian broadening: ', sigma_vel)
    print ('Template best real total broadening: ', sigma_real_broadened , 'km/s')
    print ('Sigma spectrum = sqrt(best broadening^2- resolution sigma spectrum^2): ', sigma_real , 'km/s')

    return sigma_real, total_error, min_residual, line_wave, line_flux_spec_norm, line_flux_template_norm, sigma_instrumental


#*****************************************************************************************************
# 6) SIngle line fitting
def line_fitting (wavelength, flux, wave_interval, guess_param):

    """
    This function fits an emission or absorption line with a convolution of
    a gussian function with a line, in order to account for the continuum
    slope.
    Input: wavelength and flux arrays of the spectrum, array containing the
           wavelength range of the spectral that contains the line to fit ([min_wave, max_wave])
            array containing the initial guess for the fit [y_offset, line_wave, relative_intensity,
            sigma_gauss, slope_continuum_line, intercept_continuum_line].
    Output: wavelength and normalised flux arrays of the spectrum within the fitted band,
            array of the model fit, array of the best fit parameters found.

    """

    step = wavelength[1]-wavelength[0]
    wave1 = min(wave_interval)
    wave2 = max(wave_interval)

    #isolating the region of interest
    line_wave = wavelength[(wavelength >= wave1) & (wavelength <= wave2)]
    line_flux_spec = flux[(wavelength >= wave1) & (wavelength <= wave2)]

    #normalize the spectrum
    wave_norm = line_wave[10] # guessing for now. fix it later
    epsilon_norm = step*10
    line_flux_spec_norm = spman.norm_spec(line_wave, line_flux_spec, wave_norm, epsilon_norm, line_flux_spec)

    #fitting to the spectra
    popt_spec, pconv_spec = curve_fit(uti.Gauss_slope, line_wave, line_flux_spec_norm, p0=guess_param)

    fit = uti.Gauss_slope(line_wave, *popt_spec)
    return line_wave, line_flux_spec_norm, fit, popt_spec


#*****************************************************************************************************
# 7) fitting threee gaussians to the CaT lines
def cat_fitting (wavelength, flux):

    """
    This function fits the Calcium Trilet (CaT) lines in the NIR with
    a convolution of three gussians with a line, in order to account for the continuum
    slope.
    Input: wavelength and flux arrays of the spectrum.
    Output: wavelength and normalised flux arrays of the spectrum within the fitted band,
            array of the model fit, array of the best fit parameters found.

    """

    step = wavelength[1]-wavelength[0]

    wave1 = 844
    wave2 = 872

    #extract the line flux and wavelength arrays
    line_wave = wavelength[(wavelength >= wave1) & (wavelength <= wave2)]
    line_flux_spec = flux[(wavelength >= wave1) & (wavelength <= wave2)]

    #normalize the spectrum
    wave_norm = line_wave[10] # guessing for now. fix it later
    epsilon_norm = step*10
    line_flux_spec_norm = spman.norm_spec(line_wave, line_flux_spec, wave_norm, epsilon_norm, line_flux_spec)

    #initial guesses
    y0 = 1
    x0 = 850
    a = -0.8
    sigma = 0.1
    m = 0.1
    c = 1

    guess = [y0,x0,a,sigma,m,c]
    for i in range(3):
        if i == 0:
            guess = guess
        if i == 1:
            guess+= [y0,854,-0.6, 0.2, m, c]
        if i == 2:
            guess+= [y0,866,-0.6, 0.2, m, c]

    #fitting to the spectra
    popt_spec, pconv_spec = curve_fit(uti.multiple_gauss, line_wave, line_flux_spec_norm, p0=guess)
    fit = uti.multiple_gauss(line_wave, *popt_spec)

    return line_wave, line_flux_spec_norm, fit, popt_spec


#*****************************************************************************************************
# 8) kinematics with ppxf and EMILES SSP models
def ppxf_kinematics(wavelength, flux, wave1, wave2, FWHM_gal, is_resolution_gal_constant, R, z, sigma_guess, stellar_library, additive_degree, kin_moments, kin_noise, kin_fit_gas, kin_fit_stars, kin_best_noise, with_errors_kin, custom_lib, custom_lib_folder, custom_lib_suffix, dust_correction_gas, dust_correction_stars, tied_balmer, two_stellar_components, age_model1, met_model1, age_model2, met_model2, vel_guess1, sigma_guess1, vel_guess2, sigma_guess2, mask_lines, mc_sim):

    """
     This function uses the pPXF algorith to retrieve the n kinematics moments
     by fitting SPS and gas templates to the selected wavelength
     range of the selected spectrum.
     Input: wavelength and flux arrays of the spectrum, array containing the
            wavelength range to fit ([min_wave, max_wave]), delta lambda resolution
            of the spectrum in the wavelength range considered (FWHM value, in Angstrom),
            bool constant (True) or not (False) FWHM resolution, resolving power (R)
            redshift guess, velocity dispersion guess, stellar library to use
            additive degree polynomial to use for the fit, kin moments to fit (2-6)
            constant noise estimation of the spectrum, bool with (True) or without (False)
            gas component to include, bool auto noise estimation (True) or not (False),
            bool uncertainties estimation with MonteCarlo simulations.
     Output: array containing the kinematics moments fitted, array containing the formal errors
             array of the best fit template flux, wavelength array,
             array of the model fit, array of the best fit parameters found, components found,
             S/N of the spectrum measured in the fit range, array of uncertainties in the
             kinematics moments (zero if not estimated).
    """

    ppxf_default_lib = ["emiles", "fsps", "galaxev"]
    wavelength = wavelength*10
    wave1 = wave1*10
    wave2 = wave2*10
    galaxy = flux

    line_wave = wavelength[(wavelength >= wave1) & (wavelength <= wave2)]
    line_flux_spec = galaxy[(wavelength >= wave1) & (wavelength <= wave2)]

    #updating the variables
    galaxy = line_flux_spec
    wave = line_wave

    #normalise to unity
    galaxy = galaxy/np.median(galaxy)

    # In case I have High or low redshift
    redshift = z
    high_z = 0.01
    if redshift > high_z: #with high redshift I de-redshift the spectrum to the restframe, correct the FWHM of the galaxy and set the z to zero.
        lam_range_gal = np.array([np.min(wave), np.max(wave)])/(1 + z)
        FWHM_gal /= 1 + z
        z = 0
        redshift_0 = redshift
    else: #with low redshift (z<0.01) I just let pPXF to estimate the redshift
        lam_range_gal = np.array([np.min(wave), np.max(wave)])
        redshift_0 = 0

    print('Rebinning to log')
    galaxy, ln_lam1, velscale = util.log_rebin(lam_range_gal, galaxy)

    wave = np.exp(ln_lam1) #converting the ln wavelength to wavelength, but keeping the ln sampling
    noise = np.full_like(galaxy, kin_noise) #noise per pixel

    c = 299792.458

    lam_range_temp = [lam_range_gal[0]/1.02, lam_range_gal[1]*1.02]
    sps_name = stellar_library

    # Read SPS models file
    if not custom_lib:

        #requesting the pPXF preloaded templates, if needed
        if stellar_library in ppxf_default_lib:
            ppxf_dir = Path(util.__file__).parent
            basename = f"spectra_{sps_name}_9.0.npz"
            filename = ppxf_dir / 'sps_models' / basename
            if not filename.is_file():
                url = "https://raw.githubusercontent.com/micappe/ppxf_data/main/" + basename
                request.urlretrieve(url, filename)

        #loading the templates and convolve them with the FWHM of the galaxy spectrum
        if is_resolution_gal_constant:
            if stellar_library == 'xshooter':
                pathname_xsl = os.path.join(BASE_DIR, "spectralTemplates", "xsl_mod", "*XSL_SSP*.fits" )
                sps = template.xshooter(pathname_xsl, velscale, FWHM_gal, wave_range=lam_range_temp)
            else:
                sps = lib.sps_lib(filename, velscale, FWHM_gal, lam_range=lam_range_temp)

        if not is_resolution_gal_constant:
            FWHM_gal = wave/R
            if stellar_library == 'xshooter':
                pathname_xsl = os.path.join(BASE_DIR, "spectralTemplates", "xsl_mod", "*XSL_SSP*.fits" )
                sps = template.xshooter(pathname_xsl, velscale, FWHM_gal, wave_range=lam_range_temp, R = R)
            else:
                FWHM_gal = {"lam": wave, "fwhm": FWHM_gal}
                sps = lib.sps_lib(filename, velscale, FWHM_gal, lam_range=lam_range_temp)
                #extracting the FWHM gal in numpy array, needed for the gas template building
                FWHM_gal = FWHM_gal["fwhm"]


        stars_templates = sps.templates.reshape(sps.templates.shape[0], -1)


    if custom_lib:
        pathname = custom_lib_folder + '/' + custom_lib_suffix

        if is_resolution_gal_constant:
            sps = template.miles(pathname, velscale, FWHM_gal, wave_range=lam_range_temp)

        else:
            print('')
            print('Sorry, currently the custom template selection works only with constant FWHM spectral resolution')
            print('')

        #reshaping the templates
        stars_templates = sps.templates.reshape(sps.templates.shape[0], -1)


    #loading or not the mask emission, if activated an only for stars fitting
    if kin_fit_stars and mask_lines:
        # Compute a mask for gas emission lines
        goodpix = util.determine_goodpixels(ln_lam1, lam_range_temp, z)
    else:
        goodpix = None


    error_kinematics_mc = 0

    #considering one or two stellar components
    if two_stellar_components:

        if custom_lib or stellar_library == 'xshooter':
            #retrieving age and metallicity grids
            age_grid = sps.get_full_age_grid()
            met_grid = sps.get_full_metal_grid()
            age_bins = age_grid[:,0]
            age_values = age_bins[::-1]
            met_bins = met_grid[0,:]
            met_values = met_bins
        else:
            #using the wrapper to pPXF
            sps_data_ppxf = template.SPSLibWrapper(
                filename, velscale, FWHM_gal, lam_range=lam_range_temp
            )

            age_values = sps_data_ppxf.get_age_grid()[::-1]
            met_values = sps_data_ppxf.get_metal_grid()[::-1]

        model1, i_closest1, j_closest1 = pick_ssp_template(age_model1, met_model1, age_values, met_values, sps.templates)
        model2, i_closest2, j_closest2 = pick_ssp_template(age_model2, met_model2, age_values, met_values, sps.templates)

        model1 /= np.median(model1)
        model2 /= np.median(model2)
        stars_templates = np.column_stack([model1, model2, model1, model2])



###################### Only stellar ##################
    if kin_fit_stars:
        print ('Fitting only the stellar component')

        try:

            if not two_stellar_components:
                templates = stars_templates
                vel = c*np.log(1 + z)
                start = [vel, sigma_guess]
                n_temps = stars_templates.shape[1]
                component = [0]*n_temps
                gas_component = np.array(component) > 0
                moments = kin_moments
                global_search = False # No need for single component
                t = clock()
            else:

                templates = stars_templates
                vel = c*np.log(1 + z)
                vel1 = vel + vel_guess1
                vel2 = vel + vel_guess2
                start = [[vel1, sigma_guess1], [vel2, sigma_guess2]]
                component = [0, 0, 1, 1]
                moments = [kin_moments, kin_moments]
                global_search = True #in case of two stellar components, this keyword should be set to true, according to pPXF manual
                t = clock()


            #define the dust components, if activated
            if dust_correction_stars or dust_correction_gas:
                if (dust_correction_stars and dust_correction_gas):
                    print('You are fitting only stars, discarding the dust for gas')
                    dust_gas = None
                    dust_stars = {"start": [0.1, -0.1], "bounds": [[0, 4], [-1, 0.4]], "component": ~gas_component}
                    dust = [dust_stars]
                if not dust_correction_stars and dust_correction_gas:
                    print('You do not have gas to correct for dust. No dust correction applied')
                    dust = None
                if dust_correction_stars and not dust_correction_gas:
                    print('Considering dust for the stellar component')
                    dust_stars = {"start": [0.1, -0.1], "bounds": [[0, 4], [-1, 0.4]], "component": ~gas_component}
                    dust = [dust_stars]

            else:
                dust = None

            #routine to find automatically the best noise for ppxf
            if kin_best_noise:
                print('')
                print ('Running ppxf in silent mode to find the best noise level...')

                pp = ppxf(templates, galaxy, noise, velscale, start, goodpixels = goodpix,
                    moments=moments, degree= additive_degree,
                    lam=wave, lam_temp=sps.lam_temp, quiet = True, bias =0, dust = dust, component = component, global_search = global_search) #no penalty for estimation of the noise

                nonregul_deltachi_square = round((pp.chi2 - 1)*galaxy.size, 2)
                best_noise = np.full_like(galaxy, noise*mt.sqrt(pp.chi2))
                noise = best_noise

                print ('Best noise: ', round(best_noise[0],5))
                print ('Now fitting with this noise estimation')
                print ('')


            #do the fit!
            pp = ppxf(templates, galaxy, noise, velscale, start, goodpixels = goodpix,
                moments=moments, plot = True, degree= additive_degree,
                lam=wave, lam_temp=sps.lam_temp, dust = dust, component = component, global_search = global_search)


            if np.sum(component) == 0:
                errors = pp.error*np.sqrt(pp.chi2)  # Assume the fit is good chi2/DOF=1
            else:
                errors = [array * np.sqrt(pp.chi2) for array in pp.error]

            # errors = pp.error*np.sqrt(pp.chi2)  # Assume the fit is good chi2/DOF=1
            redshift_fit = (1 + redshift_0)*np.exp(pp.sol[0]/c) - 1  # eq. (5c) C22
            redshift_err = (1 + redshift_fit)*errors[0]/c            # eq. (5d) C22

            print("Formal errors in stellar component:")
            print("     dV    dsigma   dh3      dh4")
            if np.sum(component) == 0:
                print("".join("%8.2g" % f for f in errors))
                print('Elapsed time in pPXF: %.2f s' % (clock() - t))
                prec = int(1 - np.floor(np.log10(redshift_err)))  # two digits of uncertainty
                print(f"Best-fitting redshift z = {redshift_fit:#.{prec}f} "
                    f"+/- {redshift_err:#.{prec}f}")
            else:
                stellar_uncertainties = errors[0]
                print("".join("%8.2g" % f for f in stellar_uncertainties))
                print('Elapsed time in pPXF: %.2f s' % (clock() - t))
                prec = int(1 - np.floor(np.log10(redshift_err[0])))  # two digits of uncertainty
                print(f"Best-fitting redshift z = {redshift_fit[0]:#.{prec}f} "
                    f"+/- {redshift_err[0]:#.{prec}f}")


            #output kinematics parameters
            kinematics = pp.sol
            error_kinematics = errors

            #output fit_model
            bestfit_flux = pp.bestfit
            bestfit_wavelength = wave

            #adding the mock h3, h4, h5, h6 column to the kinematic array in order to not change the main code
            all_moments = 6
            if kin_moments < all_moments:
                missing_moments = all_moments - kin_moments
                moments_to_add = np.zeros(missing_moments)
                # kinematics = np.hstack((kinematics, moments_to_add))
                # error_kinematics = np.hstack((error_kinematics, moments_to_add))

                if np.sum(component) == 0:
                    kinematics = np.hstack((kinematics, moments_to_add))
                    error_kinematics = np.hstack((error_kinematics, moments_to_add))
                else:
                    components = np.max(component)
                    for k in range (components+1):

                        kinematics[k] = np.hstack((kinematics[k], moments_to_add))
                        error_kinematics[k] = np.hstack((error_kinematics[k], moments_to_add))


            residual = galaxy - bestfit_flux
            snr = 1/np.std(residual)
            print ('S/N of the spectrum:', round(snr))


        # Uncertainties estimation with MonteCarlo simulations
            if with_errors_kin: #calculating the errors of age and metallicity with MonteCarlo simulations
                print('Calculating the uncertainties with MonteCarlo simulations')

                n_sim = mc_sim #how many simulated templates I want to create. Watch out for the computation time!

                #if fitting only one stellar component
                if np.sum(component)==0:
                    #initialising the arrays containing the n_sim simulated kinematics
                    vel_dist = []
                    sigma_dist = []
                    h3_dist = []
                    h4_dist = []
                    h5_dist = []
                    h6_dist = []

                    for i in range(n_sim):
                        noisy_template = spman.add_noise(bestfit_wavelength, bestfit_flux, snr)

                        #no regularization!
                        pp = ppxf(templates, noisy_template, noise, velscale, start, goodpixels = goodpix,
                        moments=kin_moments, degree=additive_degree,
                        lam=bestfit_wavelength, lam_temp=sps.lam_temp,
                        quiet = True, dust = dust, component = component, global_search = global_search)

                        kinematics_mc = pp.sol

                        vel_mc = int(kinematics_mc[0])
                        sigma_mc = int(kinematics_mc[1])
                        vel_dist.append(vel_mc)
                        sigma_dist.append(sigma_mc)

                        if kin_moments > 2:
                            h3_mc = round(kinematics_mc[2],3)
                            h3_dist.append(h3_mc)

                        if kin_moments > 3:
                            h4_mc = round(kinematics_mc[3],3)
                            h4_dist.append(h4_mc)

                        if kin_moments > 4:
                            h5_mc = round(kinematics_mc[4],3)
                            h5_dist.append(h5_mc)

                        if kin_moments > 5:
                            h6_mc = round(kinematics_mc[5],3)
                            h6_dist.append(h6_mc)


                    error_vel = np.std(vel_dist)
                    error_sigma = np.std(sigma_dist)
                    error_h3 = 0
                    error_h4 = 0
                    error_h5 = 0
                    error_h6 = 0

                    if kin_moments > 1:
                        error_kinematics_mc = np.column_stack((error_vel, error_sigma, error_h3, error_h4, error_h5, error_h6))

                    if kin_moments > 2:
                        error_h3 = np.std(h3_dist)
                        error_kinematics_mc = np.column_stack((error_vel, error_sigma, error_h3, error_h4, error_h5, error_h6))

                    if kin_moments > 3:
                        error_h4 = np.std(h4_dist)
                        error_kinematics_mc = np.column_stack((error_vel, error_sigma, error_h3, error_h4, error_h5, error_h6))
                    if kin_moments > 4:
                        error_h5 = np.std(h5_dist)
                        error_kinematics_mc = np.column_stack((error_vel, error_sigma, error_h3, error_h4, error_h5, error_h6))
                    if kin_moments > 5:
                        error_h6 = np.std(h6_dist)
                        error_kinematics_mc = np.column_stack((error_vel, error_sigma, error_h3, error_h4, error_h5, error_h6))


                    print('Uncertainties with MonteCarlo simulations:')
                    print(error_kinematics_mc)


                # if fitting two stellar components
                else:
                    #initialising the arrays containing the n_sim simulated kinematics
                    vel_dist1 = []
                    sigma_dist1 = []
                    h3_dist1 = []
                    h4_dist1 = []
                    h5_dist1 = []
                    h6_dist1 = []

                    vel_dist2 = []
                    sigma_dist2 = []
                    h3_dist2 = []
                    h4_dist2 = []
                    h5_dist2 = []
                    h6_dist2 = []

                    for i in range(n_sim):
                        noisy_template = spman.add_noise(bestfit_wavelength, bestfit_flux, snr)

                        #fitting the noisy templates
                        pp = ppxf(templates, noisy_template, noise, velscale, start, goodpixels = goodpix,
                        moments=kin_moments, degree=additive_degree,
                        lam=bestfit_wavelength, lam_temp=sps.lam_temp,
                        quiet = True, dust = dust, component = component, global_search = global_search)

                        kinematics_mc = pp.sol

                        vel_mc1 = int(kinematics_mc[0][0])
                        sigma_mc1 = int(kinematics_mc[0][1])
                        vel_dist1.append(vel_mc1)
                        sigma_dist1.append(sigma_mc1)

                        vel_mc2 = int(kinematics_mc[1][0])
                        sigma_mc2 = int(kinematics_mc[1][1])
                        vel_dist2.append(vel_mc2)
                        sigma_dist2.append(sigma_mc2)

                        if kin_moments > 2:
                            h3_mc1 = round(kinematics_mc[0][2],3)
                            h3_dist1.append(h3_mc1)
                            h3_mc2 = round(kinematics_mc[1][2],3)
                            h3_dist2.append(h3_mc2)

                        if kin_moments > 3:
                            h4_mc1 = round(kinematics_mc[0][3],3)
                            h4_dist1.append(h4_mc1)
                            h4_mc2 = round(kinematics_mc[1][3],3)
                            h4_dist2.append(h4_mc2)

                        if kin_moments > 4:
                            h5_mc1 = round(kinematics_mc[0][4],3)
                            h5_dist1.append(h5_mc1)
                            h5_mc2 = round(kinematics_mc[1][4],3)
                            h5_dist2.append(h5_mc2)

                        if kin_moments > 5:
                            h6_mc1 = round(kinematics_mc[0][5],3)
                            h6_dist1.append(h6_mc1)
                            h6_mc2 = round(kinematics_mc[1][5],3)
                            h6_dist2.append(h6_mc2)


                    error_vel1 = np.std(vel_dist1)
                    error_sigma1 = np.std(sigma_dist1)
                    error_h31 = 0
                    error_h41 = 0
                    error_h51 = 0
                    error_h61 = 0

                    error_vel2 = np.std(vel_dist2)
                    error_sigma2 = np.std(sigma_dist2)
                    error_h32 = 0
                    error_h42 = 0
                    error_h52 = 0
                    error_h62 = 0

                    if kin_moments > 1:
                        error_kinematics_mc = np.column_stack((error_vel1, error_sigma1, error_h31, error_h41, error_h51, error_h61, error_vel2, error_sigma2, error_h32, error_h42, error_h52, error_h62))

                    if kin_moments > 2:
                        error_h31 = np.std(h3_dist1)
                        error_h32 = np.std(h3_dist2)
                        error_kinematics_mc = np.column_stack((error_vel1, error_sigma1, error_h31, error_h41, error_h51, error_h61, error_vel2, error_sigma2, error_h32, error_h42, error_h52, error_h62))

                    if kin_moments > 3:
                        error_h41 = np.std(h4_dist1)
                        error_h42 = np.std(h4_dist2)
                        error_kinematics_mc = np.column_stack((error_vel1, error_sigma1, error_h31, error_h41, error_h51, error_h61, error_vel2, error_sigma2, error_h32, error_h42, error_h52, error_h62))

                    if kin_moments > 4:
                        error_h51 = np.std(h5_dist1)
                        error_h52 = np.std(h5_dist2)
                        error_kinematics_mc = np.column_stack((error_vel1, error_sigma1, error_h31, error_h41, error_h51, error_h61, error_vel2, error_sigma2, error_h32, error_h42, error_h52, error_h62))

                    if kin_moments > 5:
                        error_h61 = np.std(h6_dist1)
                        error_h62 = np.std(h6_dist2)
                        error_kinematics_mc = np.column_stack((error_vel1, error_sigma1, error_h31, error_h41, error_h51, error_h61, error_vel2, error_sigma2, error_h32, error_h42, error_h52, error_h62))


                    print('Uncertainties with MonteCarlo simulations:')
                    print(error_kinematics_mc)

            components = component[0] #only to return the number of gas components, that is zero!
            return kinematics, error_kinematics, bestfit_flux, bestfit_wavelength, components, snr, error_kinematics_mc

        except Exception:
            print ('ERROR')
            kinematics = error_kinematics = bestfit_flux = bestfit_wavelength = component =  snr =  error_kinematics_mc = 0


#################### WITH GAS AND STARS #########################


    if kin_fit_gas:

        print ('Fitting the stars and at least one gas component')
        try:
            # tied_balmer = False
            tie_balmer=tied_balmer
            limit_doublets=False

            #retrieving the emission lines in the wavelength range
            gas_templates, gas_names, line_wave = emission_lines(
            sps.ln_lam_temp, lam_range_gal, FWHM_gal,
            tie_balmer=tie_balmer, limit_doublets=limit_doublets, wave_galaxy = wave)

            if tie_balmer:
                dust_correction_gas = True
                print ('With tied Balmer lines, I activate the gas dust correction for you')

            templates = np.column_stack([stars_templates, gas_templates])
            vel = c*np.log(1 + z)
            start = [vel, sigma_guess]
            n_temps = stars_templates.shape[1]

            # grouping the emission lines: 1) balmer, 2) forbidden, 3) others
            n_forbidden = np.sum(["[" in a for a in gas_names])
            if not tie_balmer:
                n_balmer = np.sum(["(" in a for a in gas_names])
            else:
                n_balmer = np.sum(["Balmer" in a for a in gas_names])
                print ('Tied Balmer lines')

            n_others = np.sum(["-" in a for a in gas_names])


            #looking for the existence of at least one line of each group in the selected spectral window
            if n_forbidden !=0 and n_balmer !=0 and n_others !=0:
                ##### THREE GAS COMPONETS
                print('Balmer, forbidden and other lines')
                component = [0]*n_temps + [1]*n_balmer + [2]*n_forbidden +[3]*n_others
                gas_component = np.array(component) > 0
                moments = [kin_moments, kin_moments, kin_moments, kin_moments]
                start = [start, start, start, start]

            if n_forbidden !=0 and n_balmer !=0 and n_others == 0:
                #####
                print ('Forbidden and Balmer lines')
                component = [0]*n_temps + [1]*n_balmer + [2]*n_forbidden
                gas_component = np.array(component) > 0
                moments = [kin_moments, kin_moments, kin_moments]
                start = [start, start, start]

            if n_forbidden !=0 and n_balmer == 0 and n_others !=0:
                #####
                print ('Forbidden and other lines')
                component = [0]*n_temps + [1]*n_others + [2]*n_forbidden
                gas_component = np.array(component) > 0
                moments = [kin_moments, kin_moments, kin_moments]
                start = [start, start, start]

            if n_forbidden !=0 and n_balmer == 0 and n_others ==0:
                #######
                print ('Only forbidden lines')
                component = [0]*n_temps + [1]*n_forbidden
                gas_component = np.array(component) > 0
                moments = [kin_moments, kin_moments]
                start = [start, start]

            if n_forbidden ==0 and n_balmer != 0 and n_others ==0:
                ######
                print('Only balmer lines')
                component = [0]*n_temps + [1]*n_balmer
                gas_component = np.array(component) > 0
                moments = [kin_moments, kin_moments]
                start = [start, start]

            if n_forbidden ==0 and n_balmer != 0 and n_others !=0:
                #######
                print ('Balmer and other lines')
                component = [0]*n_temps + [1]*n_balmer [2]*n_forbidden
                gas_component = np.array(component) > 0
                moments = [kin_moments, kin_moments, kin_moments]
                start = [start, start, start]

            if n_forbidden ==0 and n_balmer == 0 and n_others !=0:
                ########
                print ('Only other lines')
                component = [0]*n_temps + [1]*n_others
                gas_component = np.array(component) > 0
                moments = [kin_moments, kin_moments]
                start = [start, start]

            if n_forbidden ==0 and n_balmer == 0 and n_others ==0:
                ########### NO GAS COMPONENT
                print ('No gas lines found. Fitting only the stellar component')
                component = [0]*n_temps
                gas_component = np.array(component) > 0
                moments = kin_moments
                start = start
            t = clock()

            #define the dust components, if activated
            if dust_correction_stars or dust_correction_gas:
                if (dust_correction_stars and dust_correction_gas):
                    print('Considering dust for stars and gas')
                    dust_gas = {"start": [0.1], "bounds": [[0, 8]], "component": gas_component}
                    dust_stars = {"start": [0.1, -0.1], "bounds": [[0, 4], [-1, 0.4]], "component": ~gas_component}
                    dust = [dust_gas, dust_stars]
                if not dust_correction_stars and dust_correction_gas:
                    print ('Considering dust for gas')
                    dust_gas = {"start": [0.1], "bounds": [[0, 8]], "component": gas_component}
                    dust = [dust_gas]
                if dust_correction_stars and not dust_correction_gas:
                    print('Considering dust for the stellar component')
                    dust_stars = {"start": [0.1, -0.1], "bounds": [[0, 4], [-1, 0.4]], "component": ~gas_component}
                    dust = [dust_stars]

            else:
                dust = None

            #routine to find automatically the best noise for ppxf
            if kin_best_noise:
                print('')
                print ('Running ppxf in silent mode to find the best noise level...')

                pp = ppxf(templates, galaxy, noise, velscale, start, goodpixels = goodpix,
                    moments=moments, degree= additive_degree,
                    lam=wave, lam_temp=sps.lam_temp,component=component, gas_component=gas_component, gas_names=gas_names, quiet = True, bias = 0, dust = dust)

                nonregul_deltachi_square = round((pp.chi2 - 1)*galaxy.size, 2)
                best_noise = np.full_like(galaxy, noise*mt.sqrt(pp.chi2))
                noise = best_noise

                print ('Best noise: ', round(best_noise[0],5))
                print ('Now fitting with this noise estimation')
                print ('')

            #finally fitting
            if component != 0: #with gas
                pp = ppxf(templates, galaxy, noise, velscale, start, goodpixels = goodpix,
                    moments=moments, plot = True, degree= additive_degree,
                    lam=wave, lam_temp=sps.lam_temp,component=component, gas_component=gas_component, gas_names=gas_names, dust = dust)

            else: #without gas
                pp = ppxf(templates, galaxy, noise, velscale, start, goodpixels = goodpix,
                    moments=moments, plot = True, degree= additive_degree,
                    lam=wave, lam_temp=sps.lam_temp,component=component, dust = dust)

            # The updated best-fitting redshift is given by the following
            # lines (using equations 5 of Cappellari 2022, arXiv, C22)
            #Calculating the formal errors:
            if component == 0:
                errors = pp.error*np.sqrt(pp.chi2)  # Assume the fit is good chi2/DOF=1
            else:
                errors = [array * np.sqrt(pp.chi2) for array in pp.error]


            redshift_fit = (1 + redshift_0)*np.exp(pp.sol[0]/c) - 1  # eq. (5c) C22
            redshift_err = (1 + redshift_fit)*errors[0]/c            # eq. (5d) C22

            print("Formal errors in stellar component:")
            print("     dV    dsigma   dh3      dh4")
            if component == 0:
                print("".join("%8.2g" % f for f in errors))
                print('Elapsed time in pPXF: %.2f s' % (clock() - t))
                prec = int(1 - np.floor(np.log10(redshift_err)))  # two digits of uncertainty
                print(f"Best-fitting redshift z = {redshift_fit:#.{prec}f} "
                    f"+/- {redshift_err:#.{prec}f}")
            else:
                stellar_uncertainties = errors[0]
                print("".join("%8.2g" % f for f in stellar_uncertainties))
                print('Elapsed time in pPXF: %.2f s' % (clock() - t))
                prec = int(1 - np.floor(np.log10(redshift_err[0])))  # two digits of uncertainty
                print(f"Best-fitting redshift z = {redshift_fit[0]:#.{prec}f} "
                    f"+/- {redshift_err[0]:#.{prec}f}")


            #output kinematics parameters
            kinematics = pp.sol
            error_kinematics = errors

            #output fit_model
            bestfit_flux = pp.bestfit
            bestfit_wavelength = wave

            #adding the mock h3, h4, h5, and h6 column to the kinematic array in order to not change the main code
            all_moments = 6
            if kin_moments < all_moments:
                missing_moments = all_moments - kin_moments
                moments_to_add = np.zeros(missing_moments)
                if component == 0:
                    kinematics = np.hstack((kinematics, moments_to_add))
                    error_kinematics = np.hstack((error_kinematics, moments_to_add))
                else:
                    components = np.max(component)
                    for k in range (components+1):

                        kinematics[k] = np.hstack((kinematics[k], moments_to_add))
                        error_kinematics[k] = np.hstack((error_kinematics[k], moments_to_add))

            #extracting the output parameters
            bestfit_flux = pp.bestfit
            residual = galaxy - bestfit_flux
            snr = 1/np.std(residual)
            print ('S/N of the spectrum:', round(snr))


        # Uncertainties estimation on the stellar kinematics with MonteCarlo simulations
            if with_errors_kin: #calculating the errors of age and metallicity with MonteCarlo simulations

                start = [vel, sigma_guess]

                print('Calculating the uncertainties with MonteCarlo simulations')
                n_sim = mc_sim #how many simulated templates I want to create. Watch out for the computation time!

                #initialising the arrays containing the n_sim simulated kinematics
                vel_dist = []
                sigma_dist = []
                h3_dist = []
                h4_dist = []
                h5_dist = []
                h6_dist = []

                for i in range(n_sim):
                    noisy_template = spman.add_noise(bestfit_wavelength, bestfit_flux, snr)

                    #fitting!
                    pp = ppxf(templates, noisy_template, noise, velscale, start, goodpixels = goodpix,
                    moments=kin_moments, degree=additive_degree,
                    lam=bestfit_wavelength, lam_temp=sps.lam_temp,
                    component=0, quiet = True, dust = dust)

                    kinematics_mc = pp.sol

                    vel_mc = int(kinematics_mc[0])
                    sigma_mc = int(kinematics_mc[1])
                    vel_dist.append(vel_mc)
                    sigma_dist.append(sigma_mc)

                    if kin_moments > 2:
                        h3_mc = round(kinematics_mc[2],3)
                        h3_dist.append(h3_mc)

                    if kin_moments > 3:
                        h4_mc = round(kinematics_mc[3],3)
                        h4_dist.append(h4_mc)

                    if kin_moments > 4:
                        h5_mc = round(kinematics_mc[4],3)
                        h5_dist.append(h5_mc)

                    if kin_moments > 5:
                        h6_mc = round(kinematics_mc[5],3)
                        h6_dist.append(h6_mc)

                # calculating the uncertainties
                error_vel = np.std(vel_dist)
                error_sigma = np.std(sigma_dist)
                error_h3 = 0
                error_h4 = 0
                error_h5 = 0
                error_h6 = 0

                if kin_moments > 1:
                    error_kinematics_mc = np.column_stack((error_vel, error_sigma, error_h3, error_h4, error_h5, error_h6))

                if kin_moments > 2:
                    error_h3 = np.std(h3_dist)
                    error_kinematics_mc = np.column_stack((error_vel, error_sigma, error_h3, error_h4, error_h5, error_h6))

                if kin_moments > 3:
                    error_h4 = np.std(h4_dist)
                    error_kinematics_mc = np.column_stack((error_vel, error_sigma, error_h3, error_h4, error_h5, error_h6))
                if kin_moments > 4:
                    error_h5 = np.std(h5_dist)
                    error_kinematics_mc = np.column_stack((error_vel, error_sigma, error_h3, error_h4, error_h5, error_h6))
                if kin_moments > 5:
                    error_h6 = np.std(h6_dist)
                    error_kinematics_mc = np.column_stack((error_vel, error_sigma, error_h3, error_h4, error_h5, error_h6))


                print('Uncertainties with MonteCarlo simulations:')
                print(error_kinematics_mc)

            return kinematics, error_kinematics, bestfit_flux, bestfit_wavelength, component, snr, error_kinematics_mc

        except AssertionError:
            print ('The selected template does not cover the wavelength range you want to fit')
            kinematics = error_kinematics = bestfit_flux = bestfit_wavelength = component =  snr =  error_kinematics_mc = 0





#*****************************************************************************************************
# 9) stellar populations with ppxf
def ppxf_pop(wave, flux, wave1, wave2, FWHM_gal, z, sigma_guess, fit_components, with_plots, with_errors, save_plot, spec_name, regul_err, additive_degree, multiplicative_degree, tied_balmer, stellar_library, dust_correction_stars, dust_correction_gas, noise_per_pix, age_range, metal_range, custom_emiles, custom_emiles_folder, custom_npz, filename_npz, mask_emission, custom_temp_suffix, best_param, best_noise_estimate, frac_chi, convolve_temp, have_user_mask, mask_ranges, nrand, lg_age, result_plot_dir):

    """
     This function uses the pPXF algorith to retrieve the properties of the
     stellar populations (age, metallicity, alpha/Fe if available) and the non
     parametric Star Formation History (SFH) of a galaxy spectrum
     by fitting SPS and (eventually) gas templates to the selected wavelength
     range of the selected spectrum.
     Input: wavelength and flux arrays of the spectrum, min and max wavelength
            to fit, delta lambda resolution of the spectrum in the wavelength
            range considered (FWHM value, in Angstrom), float redshift guess, float
            velocity dispersion guess, string wether fit the gas ('with gas')
            or just the stars ()'whitout gas'), bool whether showing (True) or not
            (False) the plots, bool whether calculate (True) or not (False) the uncertainties,
            bool whether save (True) or not (False) the plots, string name of the spectrum,
            float regularization error, int degree of additive polynomials, int degree of
            multiplicative polynomials, bool wheter to tie (True) or not (False) the Balmer lines,
            string SPS library to use, bool whether to correct (True) or not (False) for the dust
            the stellar component, bool wheter to correct (True) or not (False) for the dust the
            gas component, floas neano noise per pixel, array with the minimum and maximum age
            to consider for the models, array with the minimum and maximum metallocity [M/H] to
            consider for the models, bool whether use (True) or not (False) custom (E)MILES models,
            path of the folder containing the custom (E)MILES models to use, bool whether to mask (True)
            or not (False) the emission lines, string with the common suffix of the custom templates to use,
            bool whether estimate (True) or not (False) automatically the best noise and regul. error,
            bool whether estimate (True) or not (False) only the noise level of the spectrum,
            float fraction of the delta chi2 to reach in case of auto determination of the noise and
            regul. error parameters, bool whether convolve (True) or not (False) the SPS templates
            to the resolution of the galaxy spectrum, bool whether to include (True) or not (False)
            a user defined mask, touple wavelength interval(s) to mask, int number of bootstrap
            simulations in case the 'with_errors' option is activated (True).
     Output: array containing the kinematics moments fitted, array containing the properties fo the
             stellar populations fitted weighted in luminosity, array containing the properties fo the
             stellar populations fitted weighted in mass, array with the formal uncertainties in the
             kinematics moments, array containing the flux of the best fit template found, array containing the
             wavelength grid of the best fit template found, array containing the flux of the best
             fit gas template found, float chi2 value of the fit, float lum age lower 1sigma uncertainties
             (if with_errors = True), float lum age upper 1sigma uncertainties
             (if with_errors = True), float lum met lower 1sigma uncertainties
             (if with_errors = True), float lum met upper 1sigma uncertainties
             (if with_errors = True), float lum alpha/Fe lower 1sigma uncertainties
             (if with_errors = True), float lum alpha/Fe upper 1sigma uncertainties
             (if with_errors = True), float mass age lower 1sigma uncertainties
             (if with_errors = True), float mass age upper 1sigma uncertainties
             (if with_errors = True), float mass met lower 1sigma uncertainties
             (if with_errors = True), float mass met upper 1sigma uncertainties
             (if with_errors = True), float mass alpha/Fe lower 1sigma uncertainties
             (if with_errors = True), float mass alpha/Fe upper 1sigma uncertainties
             (if with_errors = True), array of the emission corrected flux of the spectrum,
             array of the age bins of the templates, array of the mass fraction per age bin, array of the
             cumulative mass per age bin, floar S/N measured from the residuals, array of the
             light weights calculated by ppxf, array of the mass weights calculated by ppxf
    """

    ppxf_default_lib = ["emiles", "fsps", "galaxev"]
    #converting wavelength to angstrom
    wave = wave*10
    galaxy = flux


    #selecting the input range
    wave1 = wave1*10
    wave2 = wave2*10
    line_wave = wave[(wave >= wave1) & (wave <= wave2)]
    line_flux_spec = galaxy[(wave >= wave1) & (wave <= wave2)]

    #updating the variables
    galaxy = line_flux_spec
    wave = line_wave

    #normalise to unity
    galaxy = galaxy/np.median(galaxy)

    # Setting the new lambda ranges to the rest-frame
    high_z = 0.01
    lam_range_gal = np.array([np.min(wave), np.max(wave)])/(1 + z)
    if z > high_z:
        FWHM_gal /= 1 + z
    z = 0

    #Log rebin to the restframe wavelength
    print('Rebinning to log')
    galaxy, ln_lam1, velscale = util.log_rebin(lam_range_gal, galaxy)
    wave = np.exp(ln_lam1) #converting the ln wavelength to wavelength, but keeping the ln sampling

    noise = np.full_like(galaxy, noise_per_pix) #noise per pixel
    c = 299792.458

    #setting up the wavelength range of the templates with a little of margin (1.02)
    lam_range_temp = [np.min(lam_range_gal)/1.02, np.max(lam_range_gal)*1.02]

    min_age_range = np.min(age_range)
    max_age_range = np.max(age_range)
    min_met_range = np.min(metal_range)
    max_met_range = np.max(metal_range)

    sps_name = stellar_library
    #loading the ppxf templates...
    if not custom_emiles and not custom_npz: #Using the incorporated templates with SPAN

        #requesting the pPXF preloaded templates, only if needed
        if stellar_library in ppxf_default_lib:
            ppxf_dir = Path(util.__file__).parent
            basename = f"spectra_{sps_name}_9.0.npz"
            filename = ppxf_dir / 'sps_models' / basename
            if not filename.is_file():
                url = "https://raw.githubusercontent.com/micappe/ppxf_data/main/" + basename
                request.urlretrieve(url, filename)

        if stellar_library == 'xshooter': #Xshooter templates require the custom xshooter_ppxf module
            print(stellar_library)
            pathname_xsl = os.path.join(BASE_DIR, "spectralTemplates", "xsl_mod", "*XSL_SSP*.fits" )
            if convolve_temp:
                sps = template.xshooter(pathname_xsl, velscale, FWHM_gal, norm_range=[5070, 5950],wave_range=lam_range_temp, age_range = [min_age_range, max_age_range], metal_range = [min_met_range, max_met_range]) #normalization range
            else:
                sps = template.xshooter(pathname_xsl, velscale, norm_range=[5070, 5950],wave_range=lam_range_temp, age_range = [min_age_range, max_age_range], metal_range = [min_met_range, max_met_range]) #normalization range

        elif stellar_library == 'sMILES':
            print(stellar_library)
            pathname_smiles = os.path.join(BASE_DIR, "spectralTemplates", "sMILES_afeh", "M*.fits" ) #using only the M identified, so I do not give constrain on the IMF.
            if convolve_temp:
                sps = template.smiles(pathname_smiles, velscale, FWHM_gal, norm_range=[5070, 5950], wave_range=lam_range_temp, age_range = [min_age_range, max_age_range], metal_range = [min_met_range, max_met_range]) #normalization range
            else:
                sps = template.smiles(pathname_smiles, velscale, norm_range=[5070, 5950], wave_range=lam_range_temp, age_range = [min_age_range, max_age_range], metal_range = [min_met_range, max_met_range]) #normalization range

        else: #The other templates comes with pPXF distribution and require the sps_util module
            print(stellar_library)
            if convolve_temp:
                sps = lib.sps_lib(filename, velscale, FWHM_gal, norm_range=[5070, 5950], lam_range=lam_range_temp, age_range = [min_age_range, max_age_range], metal_range = [min_met_range, max_met_range])
            else:
                sps = lib.sps_lib(filename, velscale, norm_range=[5070, 5950], lam_range=lam_range_temp, age_range = [min_age_range, max_age_range], metal_range = [min_met_range, max_met_range])

        reg_dim = sps.templates.shape[1:]
        stars_templates = sps.templates.reshape(sps.templates.shape[0], -1)

    #Loading the custom emiles templates selected by the user
    if custom_emiles and not custom_npz:
        print('Custom EMILES')
        pathname = custom_emiles_folder + '/' + custom_temp_suffix
        if convolve_temp:
            sps = template.miles(pathname, velscale, FWHM_gal, norm_range=[5070, 5950],wave_range=lam_range_temp, age_range = [min_age_range, max_age_range], metal_range = [min_met_range, max_met_range]) #normalization range
        else:
            sps = template.miles(pathname, velscale, norm_range=[5070, 5950],wave_range=lam_range_temp, age_range = [min_age_range, max_age_range], metal_range = [min_met_range, max_met_range]) #normalization range

        reg_dim = sps.templates.shape[1:]
        stars_templates = sps.templates.reshape(sps.templates.shape[0], -1)


    # Loading the custom .npz templates
    if custom_npz:
        print ('Custon templates in .npz format')
        pathname_npz = filename_npz

        if convolve_temp:
            sps = lib.sps_lib(filename_npz, velscale, FWHM_gal, norm_range=[5070, 5950], lam_range=lam_range_temp, age_range = [min_age_range, max_age_range], metal_range = [min_met_range, max_met_range])
        else:
            sps = lib.sps_lib(filename_npz, velscale, norm_range=[5070, 5950], lam_range=lam_range_temp, age_range = [min_age_range, max_age_range], metal_range = [min_met_range, max_met_range])

        reg_dim = sps.templates.shape[1:]
        stars_templates = sps.templates.reshape(sps.templates.shape[0], -1)


    #loading or not the mask emission
    if mask_emission:
        # Compute a mask for gas emission lines
        goodpix = util.determine_goodpixels(ln_lam1, lam_range_temp, z)
    else:
        goodpix = None

    #now it is time to define the user mask, if masking is activated, and transform to Angstrom.
    if have_user_mask:
        if mask_emission:
            print('Goodpix and mask cannot be used together. Continuing with user mask and neglecting goodpix...')
            goodpix = None

        if z > high_z:
            corrected_mask_ranges = [(start / (1 + z)*10, end  /(1 + z)*10) for start, end in mask_ranges]
            mask_ranges = corrected_mask_ranges
            user_mask = spman.mask_spectrum(wave, mask_ranges)
        else:
            corrected_mask_ranges = [(start *10, end *10) for start, end in mask_ranges]
            mask_ranges = corrected_mask_ranges
            user_mask = spman.mask_spectrum(wave, mask_ranges)
        mask = user_mask
    else:
        mask = None

    #definying and check on regularization value
    if regul_err > 0:
        regularization = 1/regul_err
    else:
        regularization = 0
        print ('Non-regularized fit')

    age_err_lower = 0
    age_err_upper = 0
    met_err_lower = 0
    met_err_upper = 0
    alpha_err_lower = 0
    alpha_err_upper = 0
    mass_age_err_lower = 0
    mass_age_err_upper = 0
    mass_met_err_lower = 0
    mass_met_err_upper = 0
    mass_alpha_err_lower = 0
    mass_alpha_err_upper = 0


  ###################### Now without gas ##################
    if fit_components == 'without_gas':
        print ('Fitting without gas component')

        try:
            templates = stars_templates
            vel = c*np.log(1 + z)
            start = [vel, sigma_guess]
            n_temps = stars_templates.shape[1]
            component = [0]*n_temps
            gas_component = np.array(component) > 0
            moments = 4
            start = start
            gas = False
            t = clock()


            #define the dust components, if activated
            if dust_correction_stars or dust_correction_gas:
                if (dust_correction_stars and dust_correction_gas) and gas:
                    print('Considering dust for stars and gas')
                    dust_gas = {"start": [0.1], "bounds": [[0, 8]], "component": gas_component}
                    dust_stars = {"start": [0.1, -0.1], "bounds": [[0, 4], [-1, 0.4]], "component": ~gas_component}
                    dust = [dust_gas, dust_stars]
                if (dust_correction_stars and dust_correction_gas) and not gas:
                    print('You only have stars, considering only dust for the stellar component')
                    dust_stars = {"start": [0.1, -0.1], "bounds": [[0, 4], [-1, 0.4]], "component": ~gas_component}
                    dust = [dust_stars]
                if not dust_correction_stars and dust_correction_gas and not gas:
                    print('You do not have gas to correct for dust')
                    dust = None
                if not dust_correction_stars and dust_correction_gas and gas:
                    print ('Considering dust for gas')
                    dust_gas = {"start": [0.1], "bounds": [[0, 8]], "component": gas_component}
                    dust = [dust_gas]
                if dust_correction_stars and not dust_correction_gas:
                    print('Considering dust for the stellar component')
                    dust_stars = {"start": [0.1, -0.1], "bounds": [[0, 4], [-1, 0.4]], "component": ~gas_component}
                    dust = [dust_stars]

            else:
                dust = None


            #routine to find automatically the best parameters for ppxf (noise and regul_err)
            if best_param or best_noise_estimate:
                print('')
                print ('Running ppxf in silent mode to find the best noise level...')
                try_regularization = 0


                pp = ppxf(templates, galaxy, noise, velscale, start, goodpixels = goodpix,
                    moments=moments, degree=additive_degree, mdegree=multiplicative_degree,
                    lam=wave, lam_temp=sps.lam_temp,
                    regul=try_regularization, reg_dim=reg_dim,
                    component=component, dust = dust, quiet = True)
                nonregul_deltachi_square = round((pp.chi2 - 1)*galaxy.size, 2)
                best_noise = np.full_like(galaxy, noise*mt.sqrt(pp.chi2))
                noise = best_noise

                print ('Best noise: ', round(best_noise[0],5))

                if not best_param:
                    print ('Now fitting with this noise level...')
                    print('')

            if best_param:
                #now finding the best regul_err
                max_iter = 10 #maximum iteration in order to find the best regul err
                desired_deltachi_square = round(np.sqrt(2*galaxy.size),2)
                target_deltachi_square = round(desired_deltachi_square*frac_chi, 2) #the real delta chi is a fracion of the desired one
                epsilon_chi = 0.1*target_deltachi_square # if the deltachi2 found will be up to 10% smaller than the desired delta chi2, I will accept the parameters.
                min_meaningful_regul = 0.30/n_temps #empirical value from test and errors.
                print ('Maximum delta chi2: ',desired_deltachi_square)
                print ('Trying to reach target delta chi2: ',target_deltachi_square)
                current_deltachi_square = 0 #nonregul_deltachi_square


                min_regul_err, max_regul_err = 0, 0.2 # min regul err = 0 and max likely 0.05

                #starting from the regul_err guess, if it's reasonable
                if regul_err < max_regul_err:
                    max_regul_err = regul_err

                print('')
                print ('Running iteratively ppxf in silent mode to find the best regul err...')

                #this is the regul_err you entered in the GUI
                input_regul_err = regul_err

                #finding the best regul_err with the bisection algorithm
                for k in range(max_iter):
                    print('Trying regul error: ',regul_err)

                    if regul_err < min_meaningful_regul:
                        regul_err = round(min_meaningful_regul, 3)
                        print ('')
                        print ('WARNING: your spectra are too noisy for a proper regul err estimation')
                        print ('Minimum accettable regul err ', regul_err, ' reached. Using this regardless the delta chi2 value.')
                        print('')
                        regularization = 1/regul_err
                        break



                    pp = ppxf(templates, galaxy, best_noise, velscale, start, goodpixels = goodpix,
                        moments=moments, degree=additive_degree, mdegree=multiplicative_degree,
                        lam=wave, lam_temp=sps.lam_temp,
                        regul=1/regul_err, reg_dim=reg_dim,
                        component=component, dust = dust, quiet = True)
                    current_deltachi_square = round((pp.chi2 - 1)*galaxy.size, 2)
                    print(f"Current Delta Chi^2: {(pp.chi2 - 1)*galaxy.size:#.4g}")
                    print(f"Desired Delta Chi^2: {np.sqrt(2*galaxy.size):#.4g}")
                    print('')


                    #Checking if I reached the good value according th the tolerance epsilon_chi, and only if the current deltachi is smaller or equal to the derired, not greater.
                    if abs(target_deltachi_square - current_deltachi_square) < epsilon_chi:
                        print ('Best Regul. err found!', round(regul_err,3))
                        print('Now running ppxf with noise: ', round(best_noise[0],5), 'and Regul. err: ', round(regul_err,3))
                        print('')
                        regularization = 1/regul_err
                        break

                    #simple bisection method
                    elif current_deltachi_square > target_deltachi_square:
                        min_regul_err = regul_err
                    else:
                        max_regul_err = regul_err

                    #splitting the regul err interval and trying a new value
                    regul_err = round((min_regul_err + max_regul_err) / 2, 5)

                    if k == max_iter-1:
                        print ('Convergence not reached, using the input regul err')
                        regularization = 1/input_regul_err

                    #In case the regul err is too small, I adjust the search range to include greater values
                    if k == 1:
                        if regul_err == input_regul_err:
                            print ('The regul err you entered is too small. I will guess a better value for you')
                            max_regul_err = 0.2
                            min_regul_err = regul_err
                            regul_err = max_regul_err


            #do the fit!
            pp = ppxf(templates, galaxy, noise, velscale, start, goodpixels = goodpix,
                moments=moments, degree= additive_degree, mdegree=multiplicative_degree,
                lam=wave, lam_temp=sps.lam_temp,
                regul=regularization, reg_dim=reg_dim,
                mask = mask, dust = dust)


            #setting up the result parameters
            light_weights = pp.weights[~gas_component]
            light_weights = light_weights.reshape(reg_dim)
            mass_weights = light_weights/sps.flux #converting to mass weigths
            #Normalizing
            light_weights /= light_weights.sum() # Normalize to light fractions
            mass_weights /= mass_weights.sum()

# NOTE: Following what states Cappellari (in sps_util.py), please be aware that:
# "One can use the output attribute ``.flux`` to convert light-normalized
        # weights into mass weights, without repeating the ``ppxf`` fit.
        # However, when using regularization in ``ppxf`` the results will not
        # be identical. In fact, enforcing smoothness to the light-weights is
        # not quite the same as enforcing it to the mass-weights."


            # Retrieving the mean weighted age, metallicity, and alpha values (if available).
            # For the embedded pPXF libraries I need to extract the data from auxiliary functions, since I cannot modify the sps.util function.
            if custom_emiles or stellar_library in ['sMILES', 'xshooter']:
                print('\nLuminosity weighted stellar populations:')
                info_pop = sps.mean_age_metal(light_weights, lg_age)
                print('\nMass weighted stellar populations:')
                info_pop_mass = sps.mean_age_metal(mass_weights, lg_age)
                mass_light = 0  # No photometry info available
            else:
                sps_data_ppxf = template.SPSLibWrapper(
                    filename, velscale, fwhm_gal=FWHM_gal, age_range=[min_age_range, max_age_range],
                    lam_range=lam_range_temp, metal_range=[min_met_range, max_met_range],
                    norm_range=[5070, 5950], norm_type='mean'
                )
                print('\nLuminosity weighted stellar populations:')
                info_pop = sps_data_ppxf.mean_age_metal(light_weights, lg_age)
                print('\nMass weighted stellar populations:')
                info_pop_mass = sps_data_ppxf.mean_age_metal(mass_weights, lg_age)
                mass_light = sps.mass_to_light(mass_weights, band="v")

            # Printing output infos
            print(f"\nCurrent Delta Chi^2: {(pp.chi2 - 1) * galaxy.size:#.4g}")
            print(f"Desired Delta Chi^2: {np.sqrt(2 * galaxy.size):#.4g}")
            print(f"Chi^2: {pp.chi2:#.4g}")
            print(f"Elapsed time in pPXF: {clock() - t:.2f}")

            # Extracting the output parameters
            kinematics = pp.sol
            bestfit_flux = pp.bestfit
            bestfit_wave = wave
            bestfit_gas_flux = 0.
            chi_square = pp.chi2
            emission_corrected_flux = galaxy
            errors = pp.error * np.sqrt(pp.chi2)

            residual = galaxy - bestfit_flux
            snr = 1 / np.std(residual)
            print('S/N of the spectrum:', round(snr))

            # Adjusting weights and building the SFH plot
            if stellar_library == 'sMILES' and not custom_emiles:
                reduced_mass_weights = np.sum(mass_weights, axis=2)
                mass_weights_age_bin = np.sum(reduced_mass_weights, axis=1)[::-1]

                reduced_light_weights = np.sum(light_weights, axis=2)
                light_weights_age_bin = np.sum(reduced_light_weights, axis=1)[::-1]

                #retrieving age and metallicity grids
                age_grid = sps.get_full_age_grid()
                met_grid = sps.get_full_metal_grid()
                alpha_grid = sps.get_full_alpha_grid()
                age_bins = age_grid[:,0] #extracting
                age_bins = np.mean(age_bins, axis=1)[::-1] #inverting
                met_bins = met_grid[0,:] #extracting
                met_bins = np.mean(met_bins, axis=1)[::-1] #inverting
                alpha_bins = alpha_grid[0, 0, :] #extracting
                alpha_bins = alpha_bins[::-1] #inverting

            else:
                if custom_emiles or stellar_library == 'xshooter':
                    #retrieving age and metallicity grids
                    age_grid = sps.get_full_age_grid()
                    met_grid = sps.get_full_metal_grid()
                    age_bins = age_grid[:,0] #extracting
                    age_bins = age_bins[::-1] #inverting
                    met_bins = met_grid[0,:] #extracting
                    met_bins = met_bins[::-1] #inverting
                else:
                    age_bins = sps_data_ppxf.get_age_grid()[::-1] #extracting and inverting
                    met_bins = sps_data_ppxf.get_metal_grid()[::-1] #extracting and inverting

                mass_weights_age_bin = np.sum(mass_weights, axis=1)[::-1]
                light_weights_age_bin = np.sum(light_weights, axis=1)[::-1]

            if lg_age:
                age_bins = np.log10(age_bins) + 9

            cumulative_mass = np.cumsum(mass_weights_age_bin)


            #Doing plots
            if with_plots or save_plot:

                # Creating figure and grid
                fig = plt.figure(figsize=(13, 7))
                gs = gridspec.GridSpec(3, 2, height_ratios=[1.5, 1.7, 2])


                #*********** 1) First plot: fit of the spectrum across all the columns ***********
                ax1 = fig.add_subplot(gs[0, :])
                plt.sca(ax1)
                pp.plot()
                plt.tight_layout()


                #*********** 2) Second plot: light weights map ***********
                ax2 = fig.add_subplot(gs[1, 0])
                mean_lum_age = info_pop[0]
                mean_lum_met = info_pop[1]
                plt.sca(ax2)

                #For the embedded pPXF SSP, I need to call my external function, since I cannot modify any of the pPXF files
                if not custom_emiles and stellar_library != 'sMILES' and stellar_library != 'xshooter':
                    sps_data_ppxf.plot(light_weights, lg_age, cmap='BuPu')
                else:
                    template.plot_weights(light_weights, age_grid, met_grid, lg_age, cmap='BuPu')

                #Considering log or linear age grid
                if lg_age:
                    plt.title(f"Luminosity fraction   lg<Age> = {mean_lum_age:.3g} dex, <[M/H]> = {mean_lum_met:.2g} dex", fontsize=11)
                else:
                    plt.title(f"Luminosity fraction   <Age> = {mean_lum_age:.3g} Gyr, <[M/H]> = {mean_lum_met:.2g} dex", fontsize=11)
                plt.plot(mean_lum_age, mean_lum_met, 'ro')


                #*********** 3) Third plot: mass weights map ***********
                ax3 = fig.add_subplot(gs[1, 1])
                mean_mass_age = info_pop_mass[0]
                mean_mass_met = info_pop_mass[1]
                plt.sca(ax3)

                #For the embedded pPXF SSP, I need to call my external function, since I cannot modify any of the pPXF files
                if not custom_emiles and stellar_library != 'sMILES' and stellar_library != 'xshooter':
                    sps_data_ppxf.plot(mass_weights, lg_age, cmap='BuPu')
                else:
                    template.plot_weights(mass_weights, age_grid, met_grid, lg_age, cmap='BuPu')

                #Considering log or linear age grid
                if lg_age:
                    plt.title(f"Mass fraction   lg<Age> = {mean_mass_age:.3g} dex, <[M/H]> = {mean_mass_met:.2g} dex", fontsize=11)
                else:
                    plt.title(f"Mass fraction   <Age> = {mean_mass_age:.3g} Gyr, <[M/H]> = {mean_mass_met:.2g} dex", fontsize=11)

                plt.tight_layout()


                #*********** 4) Fourth plot: light weigthed SFH ***********
                ax4 = fig.add_subplot(gs[2, 0])
                plt.sca(ax4)
                plt.plot(age_bins, light_weights_age_bin, lw=2, color='black')
                ax4.set_ylim(bottom=0)
                ax4.set_xlim(left=np.min(age_bins))
                ax4.set_xlim(right=np.max(age_bins))

                #Considering log or linear age grid
                if lg_age:
                    plt.xlabel("lg Age (dex)", fontsize=11)
                else:
                    plt.xlabel("Age (Gyr)", fontsize=11)
                plt.ylabel("Fractional luminosity", fontsize=11)
                plt.title('Luminosity weighted SFH', fontsize=10)
                plt.tight_layout()


                #*********** 5) Fifth plot: mass weigthed SFH ***********
                ax5 = fig.add_subplot(gs[2, 1])
                plt.sca(ax5)
                plt.plot(age_bins, mass_weights_age_bin, lw=2, color='black')
                ax5.set_ylim(bottom=0)
                ax5.set_xlim(left=np.min(age_bins))
                ax5.set_xlim(right=np.max(age_bins))

                #Considering log or linear age grid
                if lg_age:
                    plt.xlabel("lg Age (dex)", fontsize=11)
                else:
                    plt.xlabel("Age (Gyr)", fontsize=11)

                plt.ylabel("Fractional mass", fontsize=11)
                plt.title('Mass weighted SFH', fontsize=10)
                plt.tight_layout()

                #In case I want to save the plots
                if save_plot:
                    # result_plot_dir = 'results/plots'
                    # os.makedirs(result_plot_dir, exist_ok=True)
                    plt.savefig(result_plot_dir + '/'+ 'SFH_plots_'+ spec_name + '.png', format='png', dpi=300)
                    plt.close()
                else:
                    plt.show()
                    plt.close()


                #In case of sMILES SSPs, I show also the [alpha/Fe]-[M/H] plot in another window
                if stellar_library == 'sMILES' and not custom_emiles:
                    plt.figure(figsize=(12, 4))
                    plt.subplot(121)
                    # sps.plot_alpha(light_weights,cmap='BuPu', title = 'Luminosity fraction')
                    template.plot_alpha_weights(light_weights, alpha_grid, met_grid, cmap='BuPu', title = 'Luminosity fraction')
                    plt.plot(info_pop[1], info_pop[2], 'ro')

                    plt.subplot(122)
                    # sps.plot_alpha(mass_weights,cmap='BuPu', title = 'Mass fraction' )
                    template.plot_alpha_weights(mass_weights, alpha_grid, met_grid, cmap='BuPu', title = 'Mass fraction')
                    if save_plot:
                        plt.savefig(result_plot_dir + '/'+ 'pop_pop_alpha_weights_'+ spec_name + '.png', format='png', dpi=300)
                        plt.close()
                    else:
                        plt.show()
                        plt.close()


            #UNCERTAINTIES WITH BOOTSTRAP!
            if with_errors:
                print('Estimating the uncertainties with bootstrapping')
                print('')
                bestfit0 = pp.bestfit.copy()
                resid = galaxy - bestfit0
                start = pp.sol.copy()

                np.random.seed(123)

                weights_array = np.empty((nrand, pp.weights.size))
                age_dist = []
                met_dist = []
                mass_age_dist = []
                mass_met_dist = []

                for j in range(nrand):
                    galaxy1 = bootstrap_residuals(bestfit0, resid)


                    t = clock()
                    pp = ppxf(templates, galaxy1, noise, velscale, start, goodpixels = goodpix,
                        moments=moments, degree= additive_degree, mdegree=multiplicative_degree,
                        lam=wave, lam_temp=sps.lam_temp,
                        reg_dim=reg_dim,
                        mask = mask, dust = dust, quiet =1)

                    print(f"{j + 1}/{nrand}: Elapsed time in pPXF: {clock() - t:.2f} s")

                    weights_array[j] = pp.weights

                    #setting up the result parameters
                    light_weights_err = pp.weights[~gas_component]
                    light_weights_err = light_weights_err.reshape(reg_dim)
                    light_weights_err /= light_weights_err.sum()
                    info_pop_err = sps.mean_age_metal(light_weights_err, lg_age)
                    mass_weights_err = light_weights_err/sps.flux
                    mass_weights_err /= mass_weights_err.sum()              # Normalize to mass fractions

                    if custom_emiles or stellar_library == 'xshooter' or stellar_library =='sMILES':
                        mass_light_err = 0
                    else:
                        mass_light_err = sps.mass_to_light(mass_weights_err, band="v")

                    info_pop_mass_err = sps.mean_age_metal(mass_weights_err, lg_age)


                    age_dist.append(info_pop_err[0])
                    met_dist.append(info_pop_err[1])

                    if stellar_library == 'sMILES' and not custom_emiles:
                        alpha_dist.append(info_pop_err[2])
                        mass_alpha_dist.append(info_pop_mass_err[2])

                    age_dist.append(info_pop_err[0])
                    met_dist.append(info_pop_err[1])
                    mass_age_dist.append(info_pop_mass_err[0])
                    mass_met_dist.append(info_pop_mass_err[1])

                #sorting the values
                age_dist_sorted = np.sort(age_dist)
                met_dist_sorted = np.sort(met_dist)
                mass_age_dist_sorted = np.sort(mass_age_dist)
                mass_met_dist_sorted = np.sort(mass_met_dist)

                #calculating the 68% percentiles
                age_lower = np.percentile(age_dist_sorted, 16)
                age_median = np.percentile(age_dist_sorted, 50)
                age_upper = np.percentile(age_dist_sorted, 84)

                met_lower = np.percentile(met_dist_sorted, 16)
                met_median = np.percentile(met_dist_sorted, 50)
                met_upper = np.percentile(met_dist_sorted, 84)

                mass_age_lower = np.percentile(mass_age_dist_sorted, 16)
                mass_age_median = np.percentile(mass_age_dist_sorted, 50)
                mass_age_upper = np.percentile(mass_age_dist_sorted, 84)

                mass_met_lower = np.percentile(mass_met_dist_sorted, 16)
                mass_met_median = np.percentile(mass_met_dist_sorted, 50)
                mass_met_upper = np.percentile(mass_met_dist_sorted, 84)

                #calculating the lower and upper uncertainties
                age_err_lower = age_median - age_lower
                age_err_upper = age_upper - age_median
                met_err_lower = met_median - met_lower
                met_err_upper = met_upper - met_median

                #for mass weigthed values
                mass_age_err_lower = mass_age_median - mass_age_lower
                mass_age_err_upper = mass_age_upper - mass_age_median
                mass_met_err_lower = mass_met_median - mass_met_lower
                mass_met_err_upper = mass_met_upper - mass_met_median

                # for the sMILES models I also have the alpha/Fe to consider
                if stellar_library == 'sMILES' and not custom_emiles:
                    alpha_dist_sorted = np.sort(alpha_dist)
                    mass_alpha_dist_sorted = np.sort(mass_alpha_dist)
                    alpha_lower = np.percentile(alpha_dist_sorted, 16)
                    alpha_median = np.percentile(alpha_dist_sorted, 50)
                    alpha_upper = np.percentile(alpha_dist_sorted, 84)
                    alpha_err_lower = alpha_median - alpha_lower
                    alpha_err_upper = alpha_upper - alpha_median

                    mass_alpha_lower = np.percentile(mass_alpha_dist_sorted, 16)
                    mass_alpha_median = np.percentile(mass_alpha_dist_sorted, 50)
                    mass_alpha_upper = np.percentile(mass_alpha_dist_sorted, 84)
                    mass_alpha_err_lower = mass_alpha_median - mass_alpha_lower
                    mass_alpha_err_upper = mass_alpha_upper - mass_alpha_median
                else:
                    alpha_err_lower = 0
                    alpha_err_upper = 0
                    mass_alpha_err_lower = 0
                    mass_alpha_err_upper = 0

                print('')
                print(f"Error luminosity age: ({age_err_lower:.4g}, {age_err_upper:.4g})")
                print(f"Error luminosity met (dex)): ({met_err_lower:.4g}, {met_err_upper:.4g})")
                print(f"Error mass age: ({mass_age_err_lower:.4g}, {mass_age_err_upper:.4g})")
                print(f"Error mass met (dex)): ({mass_met_err_lower:.4g}, {mass_met_err_upper:.4g})")
                if stellar_library == 'sMILES' and not custom_emiles:
                    print(f"Error luminosity alpha/Fe (dex): ({alpha_err_lower:.4g}, {alpha_err_upper:.4g})")
                    print(f"Error mass alpha/Fe (dex): ({mass_alpha_err_lower:.4g}, {mass_alpha_err_upper:.4g})")


            return kinematics, info_pop, info_pop_mass, mass_light, errors, galaxy, bestfit_flux, bestfit_wave, bestfit_gas_flux, residual, chi_square, age_err_lower, age_err_upper, met_err_lower, met_err_upper, alpha_err_lower, alpha_err_upper, mass_age_err_lower, mass_age_err_upper, mass_met_err_lower, mass_met_err_upper, mass_alpha_err_lower, mass_alpha_err_upper, emission_corrected_flux, age_bins, light_weights_age_bin, mass_weights_age_bin, cumulative_mass, snr, light_weights, mass_weights

        except AssertionError:
            print ('The selected template does not cover the wavelength range you want to fit')
            kinematics=info_pop=info_pop_mass= mass_light= errors= galaxy= bestfit_flux= bestfit_wave= bestfit_gas_flux=residual= chi_square=age_err_lower=age_err_upper=met_err_lower=met_err_upper=alpha_err_lower=alpha_err_upper=mass_age_err_lower=mass_age_err_upper=mass_met_err_lower=mass_met_err_upper=mass_alpha_err_lower=mass_alpha_err_upper=emission_corrected_flux, age_bins, light_weights_age_bin, mass_weights_age_bin, cumulative_mass, snr, light_weights, mass_weights = 0

#################### WITH GAS #########################


    if fit_components == 'with_gas':

        print ('Fitting with at least one gas component')
        try:
            tie_balmer=tied_balmer
            limit_doublets=False

            #retrieving the emission lines in the wavelength range
            gas_templates, gas_names, line_wave = emission_lines(
            sps.ln_lam_temp, lam_range_gal, FWHM_gal,
            tie_balmer=tie_balmer, limit_doublets=limit_doublets)

            if tie_balmer:
                dust_correction_gas = True
                print ('With tied Balmer lines, I activate the gas dust correction for you')



            templates = np.column_stack([stars_templates, gas_templates])
            vel = c*np.log(1 + z)
            start = [vel, sigma_guess]
            n_temps = stars_templates.shape[1]

            # grouping the emission lines: 1) balmer, 2) forbidden, 3) others
            n_forbidden = np.sum(["[" in a for a in gas_names])
            if not tie_balmer:
                n_balmer = np.sum(["(" in a for a in gas_names])
            else:
                n_balmer = np.sum(["Balmer" in a for a in gas_names])
                print ('Tied Balmer lines')

            n_others = np.sum(["-" in a for a in gas_names])


            #looking for the existence of at least one line of each group in the selected spectral window
            if n_forbidden !=0 and n_balmer !=0 and n_others !=0:
                ##### THREE GAS COMPONETS
                gas = True
                print('Balmer, forbidden and other lines')

                component = [0]*n_temps + [1]*n_balmer + [2]*n_forbidden +[3]*n_others
                gas_component = np.array(component) > 0
                moments = [4, 2, 2, 2]
                start = [start, start, start, start]

            if n_forbidden !=0 and n_balmer !=0 and n_others == 0:
                #####
                gas = True
                print ('Forbidden and Balmer lines')
                component = [0]*n_temps + [1]*n_balmer + [2]*n_forbidden
                gas_component = np.array(component) > 0
                moments = [4, 2, 2]
                start = [start, start, start]

            if n_forbidden !=0 and n_balmer == 0 and n_others !=0:
                #####
                gas = True
                print ('Forbidden and other lines')
                component = [0]*n_temps + [1]*n_others + [2]*n_forbidden
                gas_component = np.array(component) > 0
                moments = [4, 2, 2]
                start = [start, start, start]

            if n_forbidden !=0 and n_balmer == 0 and n_others ==0:
                #######
                gas = True
                print ('Only forbidden lines')
                component = [0]*n_temps + [1]*n_forbidden
                gas_component = np.array(component) > 0
                moments = [4, 2]
                start = [start, start]

            if n_forbidden ==0 and n_balmer != 0 and n_others ==0:
                ######
                gas = True
                print('Only balmer lines')
                component = [0]*n_temps + [1]*n_balmer
                gas_component = np.array(component) > 0
                moments = [4, 2]
                start = [start, start]

            if n_forbidden ==0 and n_balmer != 0 and n_others !=0:
                #######
                gas = True
                print ('Balmer and other lines')
                component = [0]*n_temps + [1]*n_balmer [2]*n_forbidden
                gas_component = np.array(component) > 0
                moments = [4, 2, 2]
                start = [start, start, start]

            if n_forbidden ==0 and n_balmer == 0 and n_others !=0:
                ########
                gas = True
                print ('Only other lines')
                component = [0]*n_temps + [1]*n_others
                gas_component = np.array(component) > 0
                moments = [4, 2]
                start = [start, start]

            if n_forbidden ==0 and n_balmer == 0 and n_others ==0:
                ########### NO GAS COMPONENT
                gas = False
                print ('No gas lines found')
                # check_gas_cond = 0
                component = [0]*n_temps
                gas_component = np.array(component) > 0
                moments = 4
                start = start

            t = clock()


            #define the dust components, if activated
            if dust_correction_stars or dust_correction_gas:
                if (dust_correction_stars and dust_correction_gas) and gas:
                    print('Considering dust for stars and gas')
                    dust_gas = {"start": [0.1], "bounds": [[0, 8]], "component": gas_component}
                    dust_stars = {"start": [0.1, -0.1], "bounds": [[0, 4], [-1, 0.4]], "component": ~gas_component}
                    dust = [dust_gas, dust_stars]
                if (dust_correction_stars and dust_correction_gas) and not gas:
                    print('You only have stars, considering only dust for the stellar component')
                    dust_stars = {"start": [0.1, -0.1], "bounds": [[0, 4], [-1, 0.4]], "component": ~gas_component}
                    dust = [dust_stars]
                if not dust_correction_stars and dust_correction_gas and not gas:
                    print('You do not have gas to correct for dust')
                    dust = None
                if not dust_correction_stars and dust_correction_gas and gas:
                    print ('Considering dust for gas')
                    dust_gas = {"start": [0.1], "bounds": [[0, 8]], "component": gas_component}
                    dust = [dust_gas]
                if dust_correction_stars and not dust_correction_gas:
                    print('Considering dust for the stellar component')
                    dust_stars = {"start": [0.1, -0.1], "bounds": [[0, 4], [-1, 0.4]], "component": ~gas_component}
                    dust = [dust_stars]

            else:
                dust = None


            #routine to find automatically the best parameters for ppxf (noise and regul_err)
            if best_param or best_noise_estimate:
                print('')
                print ('Running ppxf in silent mode to find the best noise level...')
                try_regularization = 0

                if component != 0:
                    pp = ppxf(templates, galaxy, noise, velscale, start, goodpixels = goodpix,
                        moments=moments, degree=additive_degree, mdegree=multiplicative_degree,
                        lam=wave, lam_temp=sps.lam_temp,
                        regul=try_regularization, reg_dim=reg_dim,
                        component=component, gas_component=gas_component,
                        gas_names=gas_names, dust=dust, quiet = True)

                else:
                    pp = ppxf(templates, galaxy, noise, velscale, start, goodpixels = goodpix,
                        moments=moments, degree=additive_degree, mdegree=multiplicative_degree,
                        lam=wave, lam_temp=sps.lam_temp,
                        regul=regularization, reg_dim=reg_dim,
                        component=component, dust = dust, quiet = True)
                nonregul_deltachi_square = round((pp.chi2 - 1)*galaxy.size, 2)
                best_noise = np.full_like(galaxy, noise*mt.sqrt(pp.chi2))
                noise = best_noise

                print ('Best noise: ', round(best_noise[0],5))

                if not best_param:
                    print ('Now fitting with this noise level...')
                    print('')

            if best_param:
                #now finding the best regul_err
                max_iter = 10 #maximum iteration in order to find the best regul err
                desired_deltachi_square = round(np.sqrt(2*galaxy.size),2)
                target_deltachi_square = round(desired_deltachi_square*frac_chi, 2) #the real delta chi is a fracion of the desired one
                epsilon_chi = 0.1*target_deltachi_square # if the deltachi2 found will be up to 10% smaller than the desired delta chi2, I will accept the parameters.
                min_meaningful_regul = 0.30/n_temps #empirical value from test and errors.
                print ('Maximum delta chi2: ',desired_deltachi_square)
                print ('Trying to reach target delta chi2: ',target_deltachi_square)
                current_deltachi_square = 0 #nonregul_deltachi_square


                min_regul_err, max_regul_err = 0, 0.2 # min regul err = 0 and max likely 0.05

                #starting from the regul_err guess, if it's reasonable
                if regul_err < max_regul_err:
                    max_regul_err = regul_err

                print('')
                print ('Running iteratively ppxf in silent mode to find the best regul err...')

                #this is the regul_err you entered in the GUI
                input_regul_err = regul_err

                #finding the best regul_err with the bisection algorithm
                for k in range(max_iter):
                    print('Trying regul error: ',regul_err)

                    if regul_err < min_meaningful_regul:
                        regul_err = round(min_meaningful_regul, 3)
                        print ('')
                        print ('WARNING: your spectra are too noisy for a proper regul err estimation')
                        print ('Minimum accettable regul err ', regul_err, ' reached. Using this regardless the delta chi2 value.')
                        print('')
                        regularization = 1/regul_err
                        break

                    #with gas
                    if gas:
                        pp = ppxf(templates, galaxy, best_noise, velscale, start, goodpixels = goodpix,
                            moments=moments, degree=additive_degree, mdegree=multiplicative_degree,
                            lam=wave, lam_temp=sps.lam_temp,
                            regul=1/regul_err, reg_dim=reg_dim,
                            component=component, gas_component=gas_component,
                            gas_names=gas_names, dust=dust, mask=mask, quiet = True)

                        current_deltachi_square = round((pp.chi2 - 1)*galaxy.size, 2)
                        print(f"Current Delta Chi^2: {(pp.chi2 - 1)*galaxy.size:#.4g}")
                        print(f"Desired Delta Chi^2: {np.sqrt(2*galaxy.size):#.4g}")
                        print('Target delta Chi^2: ', target_deltachi_square)
                        print('')

                    #in case I did not find gas lines
                    else:
                        pp = ppxf(templates, galaxy, best_noise, velscale, start, goodpixels = goodpix,
                            moments=moments, degree=additive_degree, mdegree=multiplicative_degree,
                            lam=wave, lam_temp=sps.lam_temp,
                            regul=1/regul_err, reg_dim=reg_dim,
                            component=component, dust=dust, mask = mask, quiet = True)
                        current_deltachi_square = round((pp.chi2 - 1)*galaxy.size, 2)
                        print(f"Current Delta Chi^2: {(pp.chi2 - 1)*galaxy.size:#.4g}")
                        print(f"Desired Delta Chi^2: {np.sqrt(2*galaxy.size):#.4g}")
                        print('')

                    #Checking if I reached the good value according th the tolerance epsilon_chi, and only if the current deltachi is smaller or equal to the derired, not greater.
                    if abs(target_deltachi_square - current_deltachi_square) < epsilon_chi:
                        print ('Best Regul. err found!', round(regul_err,3))
                        print('Now running ppxf with noise: ', round(best_noise[0],5), 'and Regul. err: ', round(regul_err,3))
                        print('')
                        regularization = 1/regul_err
                        break

                    #simple bisection method
                    elif current_deltachi_square > target_deltachi_square:
                        min_regul_err = regul_err
                    else:
                        max_regul_err = regul_err

                    #splitting the regul err interval and trying a new value
                    regul_err = round((min_regul_err + max_regul_err) / 2, 5)

                    if k == max_iter-1:
                        print ('Convergence not reached, using the input regul err')
                        regularization = 1/input_regul_err

                    #In case the regul err is too small, I adjust the search range to include greater values
                    if k == 1:
                        if regul_err == input_regul_err:
                            print ('The regul err you entered is too small. I will guess a better value for you')
                            max_regul_err = 0.2
                            min_regul_err = regul_err
                            regul_err = max_regul_err


            #finally fitting
            if gas:
                pp = ppxf(templates, galaxy, noise, velscale, start, goodpixels = goodpix,
                    moments=moments, degree=additive_degree, mdegree=multiplicative_degree,
                    lam=wave, lam_temp=sps.lam_temp,
                    regul=regularization, reg_dim=reg_dim,
                    component=component, gas_component=gas_component,
                    gas_names=gas_names, dust=dust, mask = mask)
            else:
                pp = ppxf(templates, galaxy, noise, velscale, start, goodpixels = goodpix,
                    moments=moments, degree=additive_degree, mdegree=multiplicative_degree,
                    lam=wave, lam_temp=sps.lam_temp,
                    regul=regularization, reg_dim=reg_dim,
                    component=component, mask = mask, dust = dust)

                    #setting up the result parameters
            light_weights = pp.weights[~gas_component]
            light_weights = light_weights.reshape(reg_dim)
            mass_weights = light_weights/sps.flux #converting to mass weigths
            #Normalizing
            light_weights /= light_weights.sum() # Normalize to light fractions
            mass_weights /= mass_weights.sum()              # Normalize to mass fractions

# NOTE: Following what states Cappellari (in sps_util.py), please be aware that:
# "One can use the output attribute ``.flux`` to convert light-normalized
        # weights into mass weights, without repeating the ``ppxf`` fit.
        # However, when using regularization in ``ppxf`` the results will not
        # be identical. In fact, enforcing smoothness to the light-weights is
        # not quite the same as enforcing it to the mass-weights."


            # Retrieving the mean weighted age, metallicity, and alpha values (if available).
            # For the embedded pPXF libraries I need to extract the data from auxiliary functions, since I cannot modify the sps.util function.
            if custom_emiles or stellar_library in ['sMILES', 'xshooter']:
                print('\nLuminosity weighted stellar populations:')
                info_pop = sps.mean_age_metal(light_weights, lg_age)
                print('\nMass weighted stellar populations:')
                info_pop_mass = sps.mean_age_metal(mass_weights, lg_age)
                mass_light = 0  # No photometry info available
            else:
                sps_data_ppxf = template.SPSLibWrapper(
                    filename, velscale, fwhm_gal=FWHM_gal, age_range=[min_age_range, max_age_range],
                    lam_range=lam_range_temp, metal_range=[min_met_range, max_met_range],
                    norm_range=[5070, 5950], norm_type='mean'
                )
                print('\nLuminosity weighted stellar populations:')
                info_pop = sps_data_ppxf.mean_age_metal(light_weights, lg_age)
                print('\nMass weighted stellar populations:')
                info_pop_mass = sps_data_ppxf.mean_age_metal(mass_weights, lg_age)
                mass_light = sps.mass_to_light(mass_weights, band="v")

            # Printing output infos
            print(f"\nCurrent Delta Chi^2: {(pp.chi2 - 1) * galaxy.size:#.4g}")
            print(f"Desired Delta Chi^2: {np.sqrt(2 * galaxy.size):#.4g}")
            print(f"Chi^2: {pp.chi2:#.4g}")
            print(f"Elapsed time in pPXF: {clock() - t:.2f}")

            # Extracting the output parameters
            kinematics = pp.sol
            bestfit_flux = pp.bestfit
            bestfit_wave = wave
            bestfit_gas_flux = pp.gas_bestfit
            chi_square = pp.chi2
            residual = galaxy - bestfit_flux
            snr = 1/np.std(residual)
            errors = pp.error[0]*np.sqrt(pp.chi2)
            try:
                emission_corrected_flux = galaxy - pp.gas_bestfit
            except TypeError:
                emission_corrected_flux = galaxy

            print ('S/N of the spectrum:', round(snr))

            # Adjusting weights and building the SFH plot
            if stellar_library == 'sMILES' and not custom_emiles:
                reduced_mass_weights = np.sum(mass_weights, axis=2)
                mass_weights_age_bin = np.sum(reduced_mass_weights, axis=1)[::-1]

                reduced_light_weights = np.sum(light_weights, axis=2)
                light_weights_age_bin = np.sum(reduced_light_weights, axis=1)[::-1]

                #retrieving age and metallicity grids
                age_grid = sps.get_full_age_grid()
                met_grid = sps.get_full_metal_grid()
                alpha_grid = sps.get_full_alpha_grid()
                age_bins = age_grid[:,0] #extracting
                age_bins = np.mean(age_bins, axis=1)[::-1] #inverting
                met_bins = met_grid[0,:] #extracting
                met_bins = np.mean(met_bins, axis=1)[::-1] #inverting
                alpha_bins = alpha_grid[0, 0, :] #extracting
                alpha_bins = alpha_bins[::-1] #inverting


            else:
                if custom_emiles or stellar_library == 'xshooter':
                    #retrieving age and metallicity grids
                    age_grid = sps.get_full_age_grid()
                    met_grid = sps.get_full_metal_grid()
                    age_bins = age_grid[:,0] #extracting
                    age_bins = age_bins[::-1] #inverting
                    met_bins = met_grid[0,:] #extracting
                    met_bins = met_bins[::-1] #inverting
                else:
                    age_bins = sps_data_ppxf.get_age_grid()[::-1] #extracting and inverting
                    met_bins = sps_data_ppxf.get_metal_grid()[::-1] #extracting and inverting

                mass_weights_age_bin = np.sum(mass_weights, axis=1)[::-1]
                light_weights_age_bin = np.sum(light_weights, axis=1)[::-1]


            if lg_age:
                age_bins = np.log10(age_bins) + 9

            cumulative_mass = np.cumsum(mass_weights_age_bin)


            #Doing plots
            if with_plots or save_plot:

                # Creating figure and grid
                fig = plt.figure(figsize=(13, 7))
                gs = gridspec.GridSpec(3, 2, height_ratios=[1.5, 1.7, 2])


                #*********** 1) First plot: fit of the spectrum across all the columns ***********
                ax1 = fig.add_subplot(gs[0, :])
                plt.sca(ax1)
                pp.plot()
                plt.tight_layout()


                #*********** 2) Second plot: light weights map ***********
                ax2 = fig.add_subplot(gs[1, 0])
                mean_lum_age = info_pop[0]
                mean_lum_met = info_pop[1]
                plt.sca(ax2)

                #For the embedded pPXF SSP, I need to call my external function, since I cannot modify any of the pPXF files
                if not custom_emiles and stellar_library != 'sMILES' and stellar_library != 'xshooter':
                    sps_data_ppxf.plot(light_weights, lg_age, cmap='BuPu')
                else:
                    template.plot_weights(light_weights, age_grid, met_grid, lg_age, cmap='BuPu')

                #Considering log or linear age grid
                if lg_age:
                    plt.title(f"Luminosity fraction   lg<Age> = {mean_lum_age:.3g} dex, <[M/H]> = {mean_lum_met:.2g} dex", fontsize=11)
                else:
                    plt.title(f"Luminosity fraction   <Age> = {mean_lum_age:.3g} Gyr, <[M/H]> = {mean_lum_met:.2g} dex", fontsize=11)
                plt.plot(mean_lum_age, mean_lum_met, 'ro')


                #*********** 3) Third plot: mass weights map ***********
                ax3 = fig.add_subplot(gs[1, 1])
                mean_mass_age = info_pop_mass[0]
                mean_mass_met = info_pop_mass[1]
                plt.sca(ax3)

                #For the embedded pPXF SSP, I need to call my external function, since I cannot modify any of the pPXF files
                if not custom_emiles and stellar_library != 'sMILES' and stellar_library != 'xshooter':
                    sps_data_ppxf.plot(mass_weights, lg_age, cmap='BuPu')
                else:
                    template.plot_weights(mass_weights, age_grid, met_grid, lg_age, cmap='BuPu')

                #Considering log or linear age grid
                if lg_age:
                    plt.title(f"Mass fraction   lg<Age> = {mean_mass_age:.3g} dex, <[M/H]> = {mean_mass_met:.2g} dex", fontsize=11)
                else:
                    plt.title(f"Mass fraction   <Age> = {mean_mass_age:.3g} Gyr, <[M/H]> = {mean_mass_met:.2g} dex", fontsize=11)

                plt.tight_layout()


                #*********** 4) Fourth plot: light weigthed SFH ***********
                ax4 = fig.add_subplot(gs[2, 0])
                plt.sca(ax4)
                plt.plot(age_bins, light_weights_age_bin, lw=2, color='black')
                ax4.set_ylim(bottom=0)
                ax4.set_xlim(left=np.min(age_bins))
                ax4.set_xlim(right=np.max(age_bins))

                #Considering log or linear age grid
                if lg_age:
                    plt.xlabel("lg Age (dex)", fontsize=11)
                else:
                    plt.xlabel("Age (Gyr)", fontsize=11)
                plt.ylabel("Fractional luminosity", fontsize=11)
                plt.title('Luminosity weighted SFH', fontsize=10)
                plt.tight_layout()


                #*********** 5) Fifth plot: mass weigthed SFH ***********
                ax5 = fig.add_subplot(gs[2, 1])
                plt.sca(ax5)
                plt.plot(age_bins, mass_weights_age_bin, lw=2, color='black')
                ax5.set_ylim(bottom=0)
                ax5.set_xlim(left=np.min(age_bins))
                ax5.set_xlim(right=np.max(age_bins))

                #Considering log or linear age grid
                if lg_age:
                    plt.xlabel("lg Age (dex)", fontsize=11)
                else:
                    plt.xlabel("Age (Gyr)", fontsize=11)

                plt.ylabel("Fractional mass", fontsize=11)
                plt.title('Mass weighted SFH', fontsize=10)
                plt.tight_layout()

                #In case I want to save the plots
                if save_plot:
                    plt.savefig(result_plot_dir + '/'+ 'SFH_plots_'+ spec_name + '.png', format='png', dpi=300)
                    plt.close()
                else:
                    plt.show()
                    plt.close()


                #In case of sMILES SSPs, I show also the [alpha/Fe]-[M/H] plot in another window
                if stellar_library == 'sMILES' and not custom_emiles:
                    plt.figure(figsize=(12, 4))
                    plt.subplot(121)
                    template.plot_alpha_weights(light_weights, alpha_grid, met_grid, cmap='BuPu', title = 'Luminosity fraction')
                    plt.plot(info_pop[1], info_pop[2], 'ro')

                    plt.subplot(122)
                    template.plot_alpha_weights(mass_weights, alpha_grid, met_grid, cmap='BuPu', title = 'Mass fraction')

                    if save_plot:
                        plt.savefig(result_plot_dir + '/'+ 'pop_pop_alpha_weights_'+ spec_name + '.png', format='png', dpi=300)
                        plt.close()
                    else:
                        plt.show()
                        plt.close()


            #UNCERTAINTIES WITH BOOTSTRAP!
            if with_errors:
                print('Estimating the uncertainties with bootstrapping')
                print('')
                bestfit0 = pp.bestfit.copy()
                resid = galaxy - bestfit0
                start = pp.sol.copy()
                np.random.seed(123)

                weights_array = np.empty((nrand, pp.weights.size))
                age_dist = []
                met_dist = []
                alpha_dist = []
                mass_age_dist = []
                mass_met_dist = []
                mass_alpha_dist = []

                for j in range(nrand):
                    galaxy1 = bootstrap_residuals(bestfit0, resid)

                    t = clock() #starting the clock

                    #finally fitting
                    if gas:
                        pp = ppxf(templates, galaxy1, noise, velscale, start, goodpixels = goodpix,
                            moments=moments, degree=additive_degree, mdegree=multiplicative_degree,
                            lam=wave, lam_temp=sps.lam_temp,
                            reg_dim=reg_dim,
                            component=component, gas_component=gas_component,
                            gas_names=gas_names, dust=dust, mask = mask, quiet =1)
                    else:
                        pp = ppxf(templates, galaxy1, noise, velscale, start, goodpixels = goodpix,
                            moments=moments, degree=additive_degree, mdegree=multiplicative_degree,
                            lam=wave, lam_temp=sps.lam_temp,
                            reg_dim=reg_dim,
                            component=component, mask = mask, dust = dust, quiet =1)


                    print(f"{j + 1}/{nrand}: Elapsed time in pPXF: {clock() - t:.2f} s")

                    weights_array[j] = pp.weights


                    #setting up the result parameters
                    light_weights_err = pp.weights[~gas_component]
                    light_weights_err = light_weights_err.reshape(reg_dim)
                    light_weights_err /= light_weights_err.sum()
                    info_pop_err = sps.mean_age_metal(light_weights_err, lg_age)
                    mass_weights_err = light_weights_err/sps.flux
                    mass_weights_err /= mass_weights_err.sum()              # Normalize to mass fractions


                    if custom_emiles or stellar_library == 'xshooter' or stellar_library =='sMILES':
                        mass_light_err = 0
                    else:
                        mass_light_err = sps.mass_to_light(mass_weights_err, band="v")

                    info_pop_mass_err = sps.mean_age_metal(mass_weights_err, lg_age)


                    age_dist.append(info_pop_err[0])
                    met_dist.append(info_pop_err[1])

                    if stellar_library == 'sMILES' and not custom_emiles:
                        alpha_dist.append(info_pop_err[2])
                        mass_alpha_dist.append(info_pop_mass_err[2])

                    age_dist.append(info_pop_err[0])
                    met_dist.append(info_pop_err[1])
                    mass_age_dist.append(info_pop_mass_err[0])
                    mass_met_dist.append(info_pop_mass_err[1])

                #sorting the values
                age_dist_sorted = np.sort(age_dist)
                met_dist_sorted = np.sort(met_dist)
                mass_age_dist_sorted = np.sort(mass_age_dist)
                mass_met_dist_sorted = np.sort(mass_met_dist)

                #calculating the 68% percentiles
                age_lower = np.percentile(age_dist_sorted, 16)
                age_median = np.percentile(age_dist_sorted, 50)
                age_upper = np.percentile(age_dist_sorted, 84)

                met_lower = np.percentile(met_dist_sorted, 16)
                met_median = np.percentile(met_dist_sorted, 50)
                met_upper = np.percentile(met_dist_sorted, 84)

                mass_age_lower = np.percentile(mass_age_dist_sorted, 16)
                mass_age_median = np.percentile(mass_age_dist_sorted, 50)
                mass_age_upper = np.percentile(mass_age_dist_sorted, 84)

                mass_met_lower = np.percentile(mass_met_dist_sorted, 16)
                mass_met_median = np.percentile(mass_met_dist_sorted, 50)
                mass_met_upper = np.percentile(mass_met_dist_sorted, 84)

                #calculating the lower and upper uncertainties
                age_err_lower = age_median - age_lower
                age_err_upper = age_upper - age_median
                met_err_lower = met_median - met_lower
                met_err_upper = met_upper - met_median

                #for mass weigthed values
                mass_age_err_lower = mass_age_median - mass_age_lower
                mass_age_err_upper = mass_age_upper - mass_age_median
                mass_met_err_lower = mass_met_median - mass_met_lower
                mass_met_err_upper = mass_met_upper - mass_met_median

                # for the sMILES models I also have the alpha/Fe to consider
                if stellar_library == 'sMILES' and not custom_emiles:
                    alpha_dist_sorted = np.sort(alpha_dist)
                    mass_alpha_dist_sorted = np.sort(mass_alpha_dist)
                    alpha_lower = np.percentile(alpha_dist_sorted, 16)
                    alpha_median = np.percentile(alpha_dist_sorted, 50)
                    alpha_upper = np.percentile(alpha_dist_sorted, 84)
                    alpha_err_lower = alpha_median - alpha_lower
                    alpha_err_upper = alpha_upper - alpha_median

                    mass_alpha_lower = np.percentile(mass_alpha_dist_sorted, 16)
                    mass_alpha_median = np.percentile(mass_alpha_dist_sorted, 50)
                    mass_alpha_upper = np.percentile(mass_alpha_dist_sorted, 84)
                    mass_alpha_err_lower = mass_alpha_median - mass_alpha_lower
                    mass_alpha_err_upper = mass_alpha_upper - mass_alpha_median
                else:
                    alpha_err_lower = 0
                    alpha_err_upper = 0
                    mass_alpha_err_lower = 0
                    mass_alpha_err_upper = 0

                print('')
                print(f"Error luminosity age: ({age_err_lower:.4g}, {age_err_upper:.4g})")
                print(f"Error luminosity met (dex)): ({met_err_lower:.4g}, {met_err_upper:.4g})")
                print(f"Error mass age: ({mass_age_err_lower:.4g}, {mass_age_err_upper:.4g})")
                print(f"Error mass met (dex)): ({mass_met_err_lower:.4g}, {mass_met_err_upper:.4g})")
                if stellar_library == 'sMILES' and not custom_emiles:
                    print(f"Error luminosity alpha/Fe (dex): ({alpha_err_lower:.4g}, {alpha_err_upper:.4g})")
                    print(f"Error mass alpha/Fe (dex): ({mass_alpha_err_lower:.4g}, {mass_alpha_err_upper:.4g})")



            plt.close()

            return kinematics, info_pop, info_pop_mass, mass_light, errors, galaxy, bestfit_flux, bestfit_wave, bestfit_gas_flux, residual, chi_square, age_err_lower, age_err_upper, met_err_lower, met_err_upper, alpha_err_lower, alpha_err_upper, mass_age_err_lower, mass_age_err_upper, mass_met_err_lower, mass_met_err_upper, mass_alpha_err_lower, mass_alpha_err_upper, emission_corrected_flux, age_bins, light_weights_age_bin, mass_weights_age_bin, cumulative_mass, snr, light_weights, mass_weights

        except AssertionError:
            print ('The selected template does not cover the wavelength range you want to fit')
            kinematics=info_pop=info_pop_mass= mass_light= errors= bestfit_flux= bestfit_wave= bestfit_gas_flux=residual= chi_square=age_err_lower=age_err_upper=met_err_lower=met_err_upper=alpha_err_lower=alpha_err_upper=mass_age_err_lower=mass_age_err_upper=mass_met_err_lower=mass_met_err_upper=mass_alpha_err_lower=mass_alpha_err_upper=emission_corrected_flux, age_bins, light_weights_age_bin, mass_weights_age_bin, cumulative_mass, snr, light_weights, mass_weights = 0


#*****************************************************************************************************
# 10) Modified emission line list from the ppxf_util.py of the ppxf package.
def emission_lines(ln_lam_temp, lam_range_gal, FWHM_gal, pixel=True,
                   tie_balmer=False, limit_doublets=False, vacuum=False, wave_galaxy = None):
    """
    Generates an array of Gaussian emission lines to be used as gas templates in PPXF.

    Daniele Gasparri:
    Added the 'wave_galaxy' array, which is the galaxy wavelength array, needed to compute the FWHM_gal values in the gas
    emission lines when the FWHM_gal is not constant (i.e. when we are working with spectra with
    a fixed resolving power R. Needed for stars and gas kinematics task when selecting 'Spec. constant R resolution:'.)

    ****************************************************************************
    ADDITIONAL LINES CAN BE ADDED BY EDITING THE CODE OF THIS PROCEDURE, WHICH
    IS MEANT AS A TEMPLATE TO BE COPIED AND MODIFIED BY THE USERS AS NEEDED.
    ****************************************************************************


    Output Parameters
    -----------------

    emission_lines: ndarray
        Array of dimensions ``[ln_lam_temp.size, line_wave.size]`` containing
        the gas templates, one per array column.

    line_names: ndarray
        Array of strings with the name of each line, or group of lines'

    line_wave: ndarray
        Central wavelength of the lines, one for each gas template'

    """
    #        Balmer:     H10       H9         H8        Heps    Hdelta    Hgamma    Hbeta     Halpha
    balmer = np.array([3798.983, 3836.479, 3890.158, 3971.202, 4102.899, 4341.691, 4862.691, 6564.632])  # vacuum wavelengths

    if tie_balmer:

        # Balmer decrement for Case B recombination (T=1e4 K, ne=100 cm^-3)
        # from Storey & Hummer (1995) https://ui.adsabs.harvard.edu/abs/1995MNRAS.272...41S
        # In electronic form https://cdsarc.u-strasbg.fr/viz-bin/Cat?VI/64
        # See Table B.7 of Dopita & Sutherland 2003 https://www.amazon.com/dp/3540433627
        # Also see Table 4.2 of Osterbrock & Ferland 2006 https://www.amazon.co.uk/dp/1891389343/
        wave = balmer
        if not vacuum:
            wave = util.vac_to_air(wave)

        #if FWHM_gal is an array, I need to extract the FWHM values corresponding to the emission lines of the gas template
        if isinstance(FWHM_gal, np.ndarray):
            FWHM_gal_line = np.array([FWHM_gal[find_nearest(wave_galaxy, lw)] for lw in wave])
            gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal_line, pixel)
        else:
            gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel)

        ratios = np.array([0.0530, 0.0731, 0.105, 0.159, 0.259, 0.468, 1, 2.86])
        ratios *= wave[-2]/wave  # Account for varying log-sampled pixel size in Angstrom
        emission_lines = gauss @ ratios
        line_names = ['Balmer']
        w = (lam_range_gal[0] < wave) & (wave < lam_range_gal[1])
        line_wave = np.mean(wave[w]) if np.any(w) else np.mean(wave)

    else:

        line_wave = balmer
        if not vacuum:
            line_wave = util.vac_to_air(line_wave)
        line_names = ['(H10)', '(H9)', '(H8)', '(Heps)', '(Hdelta)', '(Hgamma)', '(Hbeta)', '(Halpha)']

        #if FWHM_gal is an array, I need to extract the FWHM values corresponding to the emission lines of the gas template
        if isinstance(FWHM_gal, np.ndarray):
            FWHM_gal_line = np.array([FWHM_gal[find_nearest(wave_galaxy, lw)] for lw in line_wave])
            emission_lines = util.gaussian(ln_lam_temp, line_wave, FWHM_gal_line, pixel)
        else:
            emission_lines = util.gaussian(ln_lam_temp, line_wave, FWHM_gal, pixel)


    if limit_doublets:

        # The line ratio of this doublet lam3727/lam3729 is constrained by
        # atomic physics to lie in the range 0.28--1.47 (e.g. fig.5.8 of
        # Osterbrock & Ferland 2006 https://www.amazon.co.uk/dp/1891389343/).
        # We model this doublet as a linear combination of two doublets with the
        # maximum and minimum ratios, to limit the ratio to the desired range.
        #       -----[OII]-----
        wave = [3727.092, 3729.875]    # vacuum wavelengths
        if not vacuum:
            wave = util.vac_to_air(wave)
        names = ['[OII]3726_d1', '[OII]3726_d2']

        if isinstance(FWHM_gal, np.ndarray):
            FWHM_gal_line = np.array([FWHM_gal[find_nearest(wave_galaxy, lw)] for lw in wave])
            gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal_line, pixel)
        else:
            gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel)

        doublets = gauss @ [[1, 1], [0.28, 1.47]]  # produces *two* doublets
        emission_lines = np.column_stack([emission_lines, doublets])
        line_names = np.append(line_names, names)
        line_wave = np.append(line_wave, wave)

        # The line ratio of this doublet lam6717/lam6731 is constrained by
        # atomic physics to lie in the range 0.44--1.43 (e.g. fig.5.8 of
        # Osterbrock & Ferland 2006 https://www.amazon.co.uk/dp/1891389343/).
        # We model this doublet as a linear combination of two doublets with the
        # maximum and minimum ratios, to limit the ratio to the desired range.
        #        -----[SII]-----
        wave = [6718.294, 6732.674]    # vacuum wavelengths
        if not vacuum:
            wave = util.vac_to_air(wave)
        names = ['[SII]6731_d1', '[SII]6731_d2']

        #if FWHM_gal is an array, I need to extract the FWHM values corresponding to the emission lines of the gas template
        if isinstance(FWHM_gal, np.ndarray):
            FWHM_gal_line = np.array([FWHM_gal[find_nearest(wave_galaxy, lw)] for lw in wave])
            gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal_line, pixel)
        else:
            gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel)

        doublets = gauss @ [[0.44, 1.43], [1, 1]]  # produces *two* doublets
        emission_lines = np.column_stack([emission_lines, doublets])
        line_names = np.append(line_names, names)
        line_wave = np.append(line_wave, wave)

    else:

        # Here the two doublets are free to have any ratio
        #         -----[OII]-----     -----[SII]-----
        # wave = [3727.092, 3729.875, 6718.294, 6732.674]  # vacuum wavelengths
        wave = [3727.092, 3729.875, 5198.4, 5201.35, 6718.294, 6732.674]  # vacuum wavelengths with NI "empirical"
        # wave = [3727.092, 3729.875, 5196.45, 5198.94, 6718.294, 6732.674] #right NI wavelengths from the emission file of GIST
        if not vacuum:
            wave = util.vac_to_air(wave)
        # names = ['[OII]3726', '[OII]3729', '[SII]6716', '[SII]6731']
        names = ['[OII]3726', '[OII]3729', '[NI]5196', '[NI]5198', '[SII]6716', '[SII]6731']

        #if FWHM_gal is an array, I need to extract the FWHM values corresponding to the emission lines of the gas template
        if isinstance(FWHM_gal, np.ndarray):
            FWHM_gal_line = np.array([FWHM_gal[find_nearest(wave_galaxy, lw)] for lw in wave])
            gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal_line, pixel)
        else:
            gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel)

        emission_lines = np.column_stack([emission_lines, gauss])
        line_names = np.append(line_names, names)
        line_wave = np.append(line_wave, wave)

    # Here the lines are free to have any ratio
    #       -----[NeIII]-----    HeII      HeI
    wave = [3968.59, 3869.86, 4687.015, 5877.243]  # vacuum wavelengths
    if not vacuum:
        wave = util.vac_to_air(wave)
    names = ['[NeIII]3968', '[NeIII]3869', '-HeII4687-', '-HeI5876-']

    if isinstance(FWHM_gal, np.ndarray):
        FWHM_gal_line = np.array([FWHM_gal[find_nearest(wave_galaxy, lw)] for lw in wave])
        gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal_line, pixel)
    else:
        gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel)

    emission_lines = np.column_stack([emission_lines, gauss])
    line_names = np.append(line_names, names)
    line_wave = np.append(line_wave, wave)

    # NIR H lines
    #       paeps      pad      pab
    wave = [10052.1, 10941.1, 12821.6]  # vacuum wavelengths
    if not vacuum:
        wave = util.vac_to_air(wave)
    names = ['-PaEps-', '-Pad-', '-Pab-']

    if isinstance(FWHM_gal, np.ndarray):
        FWHM_gal_line = np.array([FWHM_gal[find_nearest(wave_galaxy, lw)] for lw in wave])
        gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal_line, pixel)
    else:
        gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel)

    emission_lines = np.column_stack([emission_lines, gauss])
    line_names = np.append(line_names, names)
    line_wave = np.append(line_wave, wave)



    ######### Doublets with fixed ratios #########

    # To keep the flux ratio of a doublet fixed, we place the two lines in a single template
    #        -----[OIII]-----
    wave = [4960.295, 5008.240]    # vacuum wavelengths
    if not vacuum:
        wave = util.vac_to_air(wave)


    if isinstance(FWHM_gal, np.ndarray):
        FWHM_gal_line = np.array([FWHM_gal[find_nearest(wave_galaxy, lw)] for lw in wave])
        doublet = util.gaussian(ln_lam_temp, wave, FWHM_gal_line, pixel) @ [0.33, 1]
    else:
        doublet = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel) @ [0.33, 1]



    # doublet = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel) @ [0.33, 1]
    emission_lines = np.column_stack([emission_lines, doublet])
    line_names = np.append(line_names, '[OIII]5007_d')  # single template for this doublet
    line_wave = np.append(line_wave, wave[1])

    # To keep the flux ratio of a doublet fixed, we place the two lines in a single template
    #        -----[OI]-----
    wave = [6302.040, 6365.535]    # vacuum wavelengths
    if not vacuum:
        wave = util.vac_to_air(wave)



    if isinstance(FWHM_gal, np.ndarray):
        FWHM_gal_line = np.array([FWHM_gal[find_nearest(wave_galaxy, lw)] for lw in wave])
        doublet = util.gaussian(ln_lam_temp, wave, FWHM_gal_line, pixel) @ [1, 0.33]
    else:
        doublet = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel) @ [1, 0.33]



    # doublet = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel) @ [1, 0.33]
    emission_lines = np.column_stack([emission_lines, doublet])
    line_names = np.append(line_names, '[OI]6300_d')  # single template for this doublet
    line_wave = np.append(line_wave, wave[0])

    # To keep the flux ratio of a doublet fixed, we place the two lines in a single template
    #       -----[NII]-----
    wave = [6549.860, 6585.271]    # air wavelengths
    if not vacuum:
        wave = util.vac_to_air(wave)



    if isinstance(FWHM_gal, np.ndarray):
        FWHM_gal_line = np.array([FWHM_gal[find_nearest(wave_galaxy, lw)] for lw in wave])
        doublet = util.gaussian(ln_lam_temp, wave, FWHM_gal_line, pixel) @ [0.33, 1]
    else:
        doublet = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel) @ [0.33, 1]


    # doublet = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel) @ [0.33, 1]
    emission_lines = np.column_stack([emission_lines, doublet])
    line_names = np.append(line_names, '[NII]6583_d')  # single template for this doublet
    line_wave = np.append(line_wave, wave[1])

    # Only include lines falling within the estimated fitted wavelength range.
    #
    w = (lam_range_gal[0] < line_wave) & (line_wave < lam_range_gal[1])
    emission_lines = emission_lines[:, w]
    line_names = line_names[w]
    line_wave = line_wave[w]

    print('Emission lines included in gas templates:')
    print(line_names)

    return emission_lines, line_names, line_wave


#*****************************************************************************************************
# 11) stellar populations with LICK and ssp models

def lick_pop(ssp_lick_indices, ssp_lick_indices_err, ssp_model_name, interp_model):


    """
     This function calculates the properties of the stellar populations of
     galaxies via interpolation of the Lick/IDS indices Hbeta-[MgFe]' and Fem - Mgb
     using the Lick/IDS indices measured in a galaxy spectrum and the Lick/IDS indices
     measured for the following SSP models: Thomas2010 (published in Thomas et al. 2011),
     Xshooter SSP library (XSL) of Verro et al. 2022, MILES and sMILES models, all with Salpeter IMF.
     The interpolation between the measured and model Lick/IDS indices is carried out in two ways:
     1) Linear n-dimensional interpolation using the griddata function and 2) via
     Gaussian Process Regression (GPR) machine learning based model. This latter model gives
     generally better and much faster results. The learning models as well as the Lick/IDS indices measured
     for the SSP models used have been calculated with SPAN and are stored in the system_files subdirectory.
     Input: array of the Lick/IDS indices used for interpolation and measured from the galaxy spectrum,
            array of the uncertainties of the same indices, string SSP model name, string interpolation
            model to use for interpolation.
     Output: float interpolated age (in Gyr), float interpolated metallicity (dex), float
             interpolated alpha/Fe (where available, dex), float error in age, float error
             in metallicity, float error in alpha/Fe (if available).
    """


    #loading the Lick/IDS indices for the selected SSP model
    ssp_models_folder = 'system_files/'
    ssp_model_file = ssp_models_folder + ssp_model_name + '_lick.txt'
    ssp_model_file = os.path.join(BASE_DIR, ssp_model_file )

    #loading the model
    ssp_lick = np.loadtxt(ssp_model_file, delimiter=' ')

    #loading the model data in sigle arrays
    age_teo, met_teo, alpha_teo, hb_teo, mg2_teo, mgb_teo, fe5270_teo, fe5335_teo, fem_teo, mgfe_teo = ssp_lick.T

    #extracting single indices and the uncertanties
    Hbeta = ssp_lick_indices[0]
    MgFe = ssp_lick_indices[1]
    Fem = ssp_lick_indices[2]
    Mgb = ssp_lick_indices[3]

    Hbetae = ssp_lick_indices_err[0]
    MgFee = ssp_lick_indices_err[1]
    Feme = ssp_lick_indices_err[2]
    Mgbe = ssp_lick_indices_err[3]


    if ssp_model_name == 'Thomas2010' or ssp_model_name == 'miles' or ssp_model_name == 'smiles':

        #interpolate
        print('')
        print('Interpolating the values...')
        values = np.column_stack((age_teo, met_teo, alpha_teo))
        lick_indices_ml = np.array([[Hbeta, MgFe, Fem, Mgb]])
        if interp_model == 'griddata':
            #interpolation with griddata
            print('With griddata function')

            results = griddata((hb_teo, mgfe_teo, fem_teo, mgb_teo), values, lick_indices_ml, method='linear')
            age_oss = results[:, 0]
            met_oss = results[:, 1]
            alpha_oss = results[:, 2]

            age_interp = age_oss[0]
            met_interp = met_oss[0]
            alpha_interp = alpha_oss[0]

            #Uncertainties
            print('')
            print('Calculating the uncertainties...')
            sim_number = 10

            #define the arrays containing normal fluctuations of the EW of the indices with respect to their errors
            Hbeta_sim_array = np.random.normal(loc=Hbeta, scale=Hbetae, size=sim_number)
            MgFe_sim_array = np.random.normal(loc=MgFe, scale=MgFee, size=sim_number)
            Fem_sim_array = np.random.normal(loc=Fem, scale=Feme, size=sim_number)
            Mgb_sim_array = np.random.normal(loc=Mgb, scale=Mgbe, size=sim_number)

            #preparing the array of the stellar population simulated parameters
            age_sim = []
            met_sim = []
            alpha_sim = []

            #doing the simulation
            for g in range (sim_number):

                #points to interpolate
                points_for_param_sim = np.column_stack((Hbeta_sim_array[g], MgFe_sim_array[g], Fem_sim_array[g], Mgb_sim_array[g]))

                # Interpolate
                age_sim.append(griddata((hb_teo, mgfe_teo, fem_teo, mgb_teo), age_teo, points_for_param_sim, method='linear'))
                met_sim.append(griddata((hb_teo, mgfe_teo, fem_teo, mgb_teo), met_teo, points_for_param_sim, method='linear'))
                alpha_sim.append(griddata((hb_teo, mgfe_teo, fem_teo, mgb_teo), alpha_teo, points_for_param_sim, method='linear'))

            #remove the nan
            age_sim = [value for value in age_sim if not np.isnan(value)]
            met_sim = [value for value in met_sim if not np.isnan(value)]
            alpha_sim = [value for value in alpha_sim if not np.isnan(value)]

            #finally calculating the std and associate that to the error in age, met and alpha
            err_age = np.std(age_sim)
            err_met = np.std(met_sim)
            err_alpha = np.std(alpha_sim)


        if interp_model == 'GPR':
            # TEST INTERPOLATION WITH MACHINE LEARNING
            # File names for the trained models
            model_age_file = os.path.join(BASE_DIR, "system_files", ssp_model_name + "_gpr_age_model.pkl" )
            model_met_file = os.path.join(BASE_DIR, "system_files", ssp_model_name + "_gpr_met_model.pkl" )
            model_alpha_file = os.path.join(BASE_DIR, "system_files", ssp_model_name + "_gpr_alpha_model.pkl" )

            # Kernel with initial parameters
            print('With Gaussian Process Regression (GPR)')
            kernel = C(1.0, (1e-4, 4e1)) * RBF(1.0, (1e-4, 4e1))

            #better kernel for age
            kernel_age = C(1.0, (1e-4, 1e2)) * Matern(length_scale=1.0, length_scale_bounds=(1e-4, 1e2), nu=1.5)

            # Function to traing and save the model
            def train_and_save_model(X_train, y_train, kernel, filename, alpha):
                gpr = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=10, alpha=alpha)
                gpr.fit(X_train, y_train)
                joblib.dump(gpr, filename)
                return gpr

            # Load of train the model, if not already saved to the disc
            if os.path.exists(model_age_file) and os.path.exists(model_met_file) and os.path.exists(model_alpha_file):
                print('Loading trained models...')
                gpr_age = joblib.load(model_age_file)
                gpr_met = joblib.load(model_met_file)
                gpr_alpha = joblib.load(model_alpha_file)
            else:
                print('Training models...')
                X_train = np.column_stack((hb_teo, mgfe_teo, fem_teo, mgb_teo))
                gpr_age = train_and_save_model(X_train, age_teo, kernel_age, model_age_file, Hbetae)
                gpr_met = train_and_save_model(X_train, met_teo, kernel, model_met_file, MgFee)
                gpr_alpha = train_and_save_model(X_train, alpha_teo, kernel, model_alpha_file, Feme)

            # Interpolation and uncertainties with the trained model
            print('Now interpolating...')
            age_interp, err_age = gpr_age.predict(lick_indices_ml, return_std=True)
            met_interp, err_met = gpr_met.predict(lick_indices_ml, return_std=True)
            alpha_interp, err_alpha = gpr_alpha.predict(lick_indices_ml, return_std=True)

            age_interp = age_interp[0]
            met_interp = met_interp[0]
            alpha_interp = alpha_interp[0]
            err_age = err_age[0]
            err_met = err_met[0]
            err_alpha = err_alpha[0]


    #the same thing for xshooter or any model without alpha enhancment
    if ssp_model_name == 'xshooter':
        print('')
        print('Interpolating the values...')

        #interpolate
        print('')
        print('Interpolating the values...')
        values = np.column_stack((age_teo, met_teo))

        #interpolation with griddata
        print('Only with griddata function')
        lick_indices_ml = np.array([[Hbeta, MgFe]])
        results = griddata((hb_teo, mgfe_teo), values, lick_indices_ml, method='linear')
        age_oss = results[:, 0]
        met_oss = results[:, 1]
        age_interp = age_oss[0]
        met_interp = met_oss[0]
        alpha_interp = 0

        #Uncertainties
        print('')
        print('Calculating the uncertainties...')

        sim_number = 30

        #define the arrays containing normal fluctuations of the EW of the indices with respect to their errors
        Hbeta_sim_array = np.random.normal(loc=Hbeta, scale=Hbetae, size=sim_number)
        MgFe_sim_array = np.random.normal(loc=MgFe, scale=MgFee, size=sim_number)


        #preparing the array of the stellar population simulated parameters
        age_sim = []
        met_sim = []

        #doing the simulation
        for g in range (sim_number):

            #points to interpolate
            points_for_param_sim = np.column_stack((Hbeta_sim_array[g], MgFe_sim_array[g]))

            # Interpolate
            age_sim.append(griddata((hb_teo, mgfe_teo), age_teo, points_for_param_sim, method='linear'))
            met_sim.append(griddata((hb_teo, mgfe_teo), met_teo, points_for_param_sim, method='linear'))


        #remove the nan
        age_sim = [value for value in age_sim if not np.isnan(value)]
        met_sim = [value for value in met_sim if not np.isnan(value)]


        #finally calculating the std and associates it to the error in age, met and alpha
        err_age = np.std(age_sim)
        err_met = np.std(met_sim)

        err_alpha = 0

    return age_interp, met_interp, alpha_interp, err_age, err_met, err_alpha




    #**********************plotting**************************

def lick_grids(ssp_model_name, ssp_lick_indices, ssp_lick_indices_err, age, show_plot, save_plot, spectra_list_name, result_plot_dir):

    """
     This function plots the measured values of the Hbeta-[MgFe]' and Fem-Mgb indices for the galaxy spectrum
     into the index grid of the selected SSP models.
     Input: string SSP model name, array of the Lick/IDS indices used for interpolation and measured from the galaxy spectrum,
            array of the uncertainties of the same indices, float mean luminosity age estimated via interpolation
            with the SSP models for the n spectra or value for the single spectrum,
            bool whether show (True) or not (False) the plot, bool whether to ssave (True) or not (False) the plot
            in a png high resolution image, string name of the spectra to plot, or the single spectrum.
     Output: A matplot window or a PNG image with the model grids and the data point(s).

    """


    # if with_plots or save_plots:
    ssp_models_folder = 'system_files/'
    ssp_model_file = ssp_models_folder + ssp_model_name + '_lick.txt'
    ssp_model_file = os.path.join(BASE_DIR, ssp_model_file )

    #extracting single indices and the uncertanties
    Hbeta = ssp_lick_indices[:,0]
    MgFe = ssp_lick_indices[:,1]
    Fem = ssp_lick_indices[:,2]
    Mgb = ssp_lick_indices[:,3]

    Hbetae = ssp_lick_indices_err[:,0]
    MgFee = ssp_lick_indices_err[:,1]
    Feme = ssp_lick_indices_err[:,2]
    Mgbe = ssp_lick_indices_err[:,3]

    if ssp_model_name == 'Thomas2010':
        data = np.genfromtxt(os.path.join(BASE_DIR, ssp_model_file), delimiter=' ', skip_header=True)
        met_values = [-1.35, -0.33, 0, 0.35, 0.67]
        age_values = [0.6, 0.8, 1, 2, 4, 10, 15]
        alpha_values = [-0.3, 0, 0.3, 0.5]

    if ssp_model_name == 'xshooter':
        data = np.genfromtxt(os.path.join(BASE_DIR, "system_files", "xshooter_lick_plot.txt"), delimiter=' ', skip_header=True)
        met_values = [-1.2, -0.8, -0.4, 0, 0.1, 0.2]
        age_values = [0.79, 1, 2, 3.98, 6.31, 10, 15.85]
        alpha_values = [0]

    if ssp_model_name == 'miles':
        data = np.genfromtxt(os.path.join(BASE_DIR, "system_files", "miles_lick_plot.txt"), delimiter=' ', skip_header=True)
        met_values = [-0.96, -0.66, -0.35, 0.06, 0.26, 0.4]
        age_values = [0.6, 0.8, 1., 2., 4., 10., 14.]
        alpha_values = [0, 0.4]

    if ssp_model_name == 'smiles':
        data = np.genfromtxt(os.path.join(BASE_DIR, "system_files", "smiles_lick_plot.txt"), delimiter=' ', skip_header=True)
        met_values = [-0.96, -0.66, -0.35, 0.06, 0.26]
        age_values = [0.6, 0.8, 1., 2., 4., 10., 14.]
        alpha_values = [-0.2, 0, 0.2, 0.4, 0.6]


    age_values = np.array(age_values)
    met_values = np.array(met_values)

    #round to the closest value. Useful for the mgb-fem grid which strongly depends on age, so I need to fix that
    age_closest_alpha_plot = min(age_values, key=lambda x: abs(x - age))

    plt.figure(figsize=(10, 6))

    # age
    plt.subplot(1, 2, 1)

    alpha_value_grid = 0. # fixing a solar alpha enhancment for the mgfe-hbeta plot


    if ssp_model_name == 'Thomas2010':

        for i in range(len(age_values)):
            age_idx = np.where((data[:, 0] == age_values[i]) & (data[:, 2] == alpha_value_grid))[0]
            plt.plot(data[age_idx, 9], data[age_idx, 3], color='black', linewidth=-1)
            plt.text(data[age_idx[4], 9] + 0.1, data[age_idx[4], 3] + 0., f'{age_values[i]} Gyr', fontsize=10, color='darkgrey')


        # met
        for h in range(len(met_values)):
            met_idx = np.where((data[:, 1] == met_values[h]) & (data[:, 2] == alpha_value_grid))[0]
            plt.plot(data[met_idx, 9], data[met_idx, 3], color='black', linewidth=1)
            plt.text(data[met_idx[15], 9] - 0.1, data[met_idx[15], 3] - 0.25, f'{met_values[h]}', fontsize=10, color='darkgray')
            if h == 0:
                plt.text(data[met_idx[15], 9] - 0.7, data[met_idx[15], 3] - 0.25, f'[Fe/H]=', fontsize=10, color='darkgray')

        plt.xlabel("[MgFe]' (\u00c5)", fontsize=14)
        plt.ylabel('H\u03B2 (\u00c5)', fontsize=14)
        plt.tick_params(axis='both', which='both', direction='in', left=True, right=True, top=True, bottom=True)
        plt.minorticks_on()
        plt.gca().yaxis.set_minor_locator(MultipleLocator(0.2))
        plt.gca().xaxis.set_minor_locator(MultipleLocator(0.25))
        plt.xlim(0.4, 5.80)
        plt.ylim(0.7, 5.3)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)

        #plotting the value
        plt.scatter(MgFe, Hbeta, color='red', s = 16)
        plt.errorbar(MgFe, Hbeta, xerr=MgFee, yerr=Hbetae, linestyle='None', ecolor = 'black', capsize=2)

        # alpha values, second plot
        plt.subplot(1, 2, 2)
        for i in range(len(alpha_values)):
            alpha_idx = np.where((data[:, 2] == alpha_values[i]) & (data[:, 0] == age_closest_alpha_plot))[0]
            plt.plot(data[alpha_idx, 5], data[alpha_idx, 8], color='black', linewidth=-1)
            plt.text(data[alpha_idx[4], 5] - 0.5, data[alpha_idx[4], 8] + 0.05, f' [\u03B1/Fe]={alpha_values[i]}', fontsize=10,
                    color='darkgray')

        # #minor ticks
        plt.minorticks_on()
        plt.gca().yaxis.set_minor_locator(MultipleLocator(0.2))
        plt.gca().xaxis.set_minor_locator(MultipleLocator(0.25))

        plt.yticks(np.arange(1, 6), fontsize=14)
        plt.xlabel('Mgb (\u00c5)', fontsize=14)
        plt.ylabel('<Fe> (\u00c5)', fontsize=14)
        plt.tick_params(axis='both', which='both', direction='in', left=True, right=True, top=True, bottom=True)
        plt.xlim(0.4, 6.5)
        plt.ylim(0.7, 5.3)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)

        #plotting the value
        plt.scatter(Mgb, Fem, color='red', s = 16)
        plt.errorbar(Mgb, Fem, xerr=Mgbe, yerr=Feme, linestyle='None', ecolor = 'black', capsize=2)

        plt.tight_layout()

        if show_plot:
            plt.show()
            plt.close()
        if save_plot:
            # result_plot_dir = 'results/plots'
            # os.makedirs(result_plot_dir, exist_ok=True)
            model_grids_file = result_plot_dir + '/'+ 'index_grids_' + spectra_list_name + '.png'
            plt.savefig(model_grids_file, format='png', dpi=300) #I must save a png image because the eps file does not reproduce well the grids. Don't know why...
            print ('Index-index diagrams for SSP models and data points saved in: ', model_grids_file)
            plt.close()






    if ssp_model_name == 'xshooter':
        for i in range(len(age_values)):
            age_idx = np.where((data[:, 0] == age_values[i]) & (data[:, 2] == alpha_value_grid))[0]
            plt.plot(data[age_idx, 9], data[age_idx, 3], color='black', linewidth=-1)
            plt.text(data[age_idx[4], 9] + 0.25, data[age_idx[4], 3] - 0., f'{age_values[i]} Gyr', fontsize=10, color='darkgrey')


        # met
        for h in range(len(met_values)):
            met_idx = np.where((data[:, 1] == met_values[h]) & (data[:, 2] == alpha_value_grid))[0]
            plt.plot(data[met_idx, 9], data[met_idx, 3], color='black', linewidth=1)

            if h == 0:
                plt.text(data[met_idx[25], 9] - 0.6, data[met_idx[25], 3] - 0.2, f'[Fe/H]=', fontsize=10, color='darkgray')
                plt.text(data[met_idx[25], 9] -0.1 , data[met_idx[25], 3] - 0.2, f'{met_values[h]}', fontsize=10, color='darkgray')

            elif h == (len(met_values)-2):
                plt.text(data[met_idx[25], 9] +0 , data[met_idx[25], 3] - 0.2, f'{met_values[h]}', fontsize=10, color='darkgray')

            elif h == (len(met_values)-1):
                plt.text(data[met_idx[25], 9] +0.1 , data[met_idx[25], 3] - 0.2, f'{met_values[h]}', fontsize=10, color='darkgray')
            else:
                plt.text(data[met_idx[25], 9] - 0.1, data[met_idx[25], 3] - 0.2, f'{met_values[h]}', fontsize=10, color='darkgray')


        plt.xlabel("[MgFe]' (\u00c5)", fontsize=14)
        plt.ylabel('H\u03B2 (\u00c5)', fontsize=14)
        plt.tick_params(axis='both', which='both', direction='in', left=True, right=True, top=True, bottom=True)
        plt.minorticks_on()
        plt.gca().yaxis.set_minor_locator(MultipleLocator(0.2))
        plt.gca().xaxis.set_minor_locator(MultipleLocator(0.25))
        plt.xlim(0.5, 4.8)
        plt.ylim(1.2, 8)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)

        #plotting the value
        plt.scatter(MgFe, Hbeta, color='red', s = 16)
        plt.errorbar(MgFe, Hbeta, xerr=MgFee, yerr=Hbetae, linestyle='None', ecolor = 'black', capsize=2)

        # alpha values, second plot
        plt.subplot(1, 2, 2)
        for i in range(len(alpha_values)):
            alpha_idx = np.where((data[:, 2] == alpha_values[i]) & (data[:, 0] == age_closest_alpha_plot))[0]
            plt.plot(data[alpha_idx, 5], data[alpha_idx, 8], color='black', linewidth=-1)
            plt.text(data[alpha_idx[4], 5] - 0.5, data[alpha_idx[4], 8] + 0.05, f' [\u03B1/Fe]={alpha_values[i]}', fontsize=10,
                    color='darkgray')

        #minor ticks
        plt.minorticks_on()
        plt.gca().yaxis.set_minor_locator(MultipleLocator(0.2))
        plt.gca().xaxis.set_minor_locator(MultipleLocator(0.25))

        plt.yticks(np.arange(1, 6), fontsize=14)
        plt.xlabel('Mgb (\u00c5)', fontsize=14)
        plt.ylabel('<Fe> (\u00c5)', fontsize=14)
        plt.tick_params(axis='both', which='both', direction='in', left=True, right=True, top=True, bottom=True)
        plt.xlim(0.4, 6.5)
        plt.ylim(0.7, 5.3)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)

        #plotting the value
        plt.scatter(Mgb, Fem, color='black', s = 16)
        plt.errorbar(Mgb, Fem, xerr=Mgbe, yerr=Feme, linestyle='None', ecolor = 'black', capsize=2)

        plt.tight_layout()

        if show_plot:
            plt.show()
            plt.close()
        if save_plot:
            # result_plot_dir = 'results/plots'
            # os.makedirs(result_plot_dir, exist_ok=True)
            model_grids_file = result_plot_dir + '/'+ 'index_grids_' + spectra_list_name + '.png'
            plt.savefig(model_grids_file, format='png', dpi=300) #I must save a png image because the eps file does not reproduce well the grids. Don't know why...
            print ('Index-index diagrams for SSP models and data points saved in: ', model_grids_file)
            plt.close()



    if ssp_model_name == 'miles':
        for i in range(len(age_values)):
            age_idx = np.where((data[:, 0] == age_values[i]) & (data[:, 2] == alpha_value_grid))[0]
            plt.plot(data[age_idx, 9], data[age_idx, 3], color='black', linewidth=-1)
            plt.text(data[age_idx[4], 9] + 0.4, data[age_idx[4], 3] - 0.1, f'{age_values[i]} Gyr', fontsize=10, color='darkgrey')


        # met
        for h in range(len(met_values)):
            met_idx = np.where((data[:, 1] == met_values[h]) & (data[:, 2] == alpha_value_grid))[0]
            plt.plot(data[met_idx, 9], data[met_idx, 3], color='black', linewidth=1)

            if h == 0:
                plt.text(data[met_idx[35], 9] - 0.7, data[met_idx[35], 3] - 0.7, f'[Fe/H]=', fontsize=10, color='darkgray')
                plt.text(data[met_idx[35], 9] -0.2 , data[met_idx[35], 3] - 0.7, f'{met_values[h]}', fontsize=10, color='darkgray')
            elif h == (len(met_values)-2):
                plt.text(data[met_idx[35], 9] -0.1 , data[met_idx[35], 3] - 0.25, f'{met_values[h]}', fontsize=10, color='darkgray')

            elif h == (len(met_values)-1):
                plt.text(data[met_idx[35], 9] +0.1 , data[met_idx[35], 3] - 0.2, f'{met_values[h]}', fontsize=10, color='darkgray')
            else:
                plt.text(data[met_idx[35], 9] - 0.1, data[met_idx[35], 3] - 0.2, f'{met_values[h]}', fontsize=10, color='darkgray')


        plt.xlabel("[MgFe]' (\u00c5)", fontsize=14)
        plt.ylabel('H\u03B2 (\u00c5)', fontsize=14)
        plt.tick_params(axis='both', which='both', direction='in', left=True, right=True, top=True, bottom=True)
        plt.minorticks_on()
        plt.gca().yaxis.set_minor_locator(MultipleLocator(0.2))
        plt.gca().xaxis.set_minor_locator(MultipleLocator(0.25))
        plt.xlim(0.5, 5)
        plt.ylim(1.2, 6.3)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)

        #plotting the value
        plt.scatter(MgFe, Hbeta, color='red', s = 16)
        plt.errorbar(MgFe, Hbeta, xerr=MgFee, yerr=Hbetae, linestyle='None', ecolor = 'black', capsize=2)

        # alpha values, second plot
        plt.subplot(1, 2, 2)
        for i in range(len(alpha_values)):
            alpha_idx = np.where((data[:, 2] == alpha_values[i]) & (data[:, 0] == age_closest_alpha_plot))[0]
            plt.plot(data[alpha_idx, 5], data[alpha_idx, 8], color='black', linewidth=-1)
            plt.text(data[alpha_idx[4], 5] - 0.5, data[alpha_idx[4], 8] + 0.05, f' [\u03B1/Fe]={alpha_values[i]}', fontsize=10,
                    color='darkgray')

        # minor ticks
        plt.minorticks_on()
        plt.gca().yaxis.set_minor_locator(MultipleLocator(0.2))
        plt.gca().xaxis.set_minor_locator(MultipleLocator(0.25))

        plt.yticks(np.arange(1, 6), fontsize=14)
        plt.xlabel('Mgb (\u00c5)', fontsize=14)
        plt.ylabel('<Fe> (\u00c5)', fontsize=14)
        plt.tick_params(axis='both', which='both', direction='in', left=True, right=True, top=True, bottom=True)
        plt.xlim(1, 6.5)
        plt.ylim(1.1, 5.3)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)

        #plotting the value
        plt.scatter(Mgb, Fem, color='black', s = 16)
        plt.errorbar(Mgb, Fem, xerr=Mgbe, yerr=Feme, linestyle='None', ecolor = 'black', capsize=2)

        plt.tight_layout()

        if show_plot:
            plt.show()
            plt.close()
        if save_plot:
            # result_plot_dir = 'results/plots'
            # os.makedirs(result_plot_dir, exist_ok=True)
            model_grids_file = result_plot_dir + '/'+ 'index_grids_' + spectra_list_name + '.png'
            plt.savefig(model_grids_file, format='png', dpi=300) #I must save a png image because the eps file does not reproduce well the grids. Don't know why...
            print ('Index-index diagrams for SSP models and data points saved in: ', model_grids_file)
            plt.close()



    if ssp_model_name == 'smiles':

        for i in range(len(age_values)):
            age_idx = np.where((data[:, 0] == age_values[i]) & (data[:, 2] == alpha_value_grid))[0]
            plt.plot(data[age_idx, 9], data[age_idx, 3], color='black', linewidth=-1)
            plt.text(data[age_idx[4], 9] + 0.1, data[age_idx[4], 3] + 0., f'{age_values[i]} Gyr', fontsize=10, color='darkgrey')


        # met
        for h in range(len(met_values)):
            met_idx = np.where((data[:, 1] == met_values[h]) & (data[:, 2] == alpha_value_grid))[0]
            plt.plot(data[met_idx, 9], data[met_idx, 3], color='black', linewidth=1)

            if h == 0:
                plt.text(data[met_idx[36], 9] - 0.7, data[met_idx[36], 3] - 0.75, f'[Fe/H]=', fontsize=10, color='darkgray')
                plt.text(data[met_idx[36], 9] -0.2 , data[met_idx[35], 3] - 0.75, f'{met_values[h]}', fontsize=10, color='darkgray')
            else:
                plt.text(data[met_idx[36], 9] - 0.1, data[met_idx[36], 3] - 0.25, f'{met_values[h]}', fontsize=10, color='darkgray')

        plt.xlabel("[MgFe]' (\u00c5)", fontsize=14)
        plt.ylabel('H\u03B2 (\u00c5)', fontsize=14)
        plt.tick_params(axis='both', which='both', direction='in', left=True, right=True, top=True, bottom=True)
        plt.minorticks_on()
        plt.gca().yaxis.set_minor_locator(MultipleLocator(0.2))
        plt.gca().xaxis.set_minor_locator(MultipleLocator(0.25))
        plt.xlim(0.6, 4.6)
        plt.ylim(1.2, 6.2)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)

        #plotting the value
        plt.scatter(MgFe, Hbeta, color='red', s = 16)
        plt.errorbar(MgFe, Hbeta, xerr=MgFee, yerr=Hbetae, linestyle='None', ecolor = 'black', capsize=2)

        # alpha values, second plot
        plt.subplot(1, 2, 2)
        for i in range(len(alpha_values)):
            alpha_idx = np.where((data[:, 2] == alpha_values[i]) & (data[:, 0] == age_closest_alpha_plot))[0]
            plt.plot(data[alpha_idx, 5], data[alpha_idx, 8], color='black', linewidth=-1)
            plt.text(data[alpha_idx[4], 5] - 0.5, data[alpha_idx[4], 8] + 0.05, f' [\u03B1/Fe]={alpha_values[i]}', fontsize=10,
                    color='darkgray')

        # #minor ticks
        plt.minorticks_on()
        plt.gca().yaxis.set_minor_locator(MultipleLocator(0.2))
        plt.gca().xaxis.set_minor_locator(MultipleLocator(0.25))

        plt.yticks(np.arange(1, 6), fontsize=14)
        plt.xlabel('Mgb (\u00c5)', fontsize=14)
        plt.ylabel('<Fe> (\u00c5)', fontsize=14)
        plt.tick_params(axis='both', which='both', direction='in', left=True, right=True, top=True, bottom=True)
        plt.xlim(0.4, 6.5)
        plt.ylim(0.7, 5.3)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)

        #plotting the value
        plt.scatter(Mgb, Fem, color='black', s = 16)
        plt.errorbar(Mgb, Fem, xerr=Mgbe, yerr=Feme, linestyle='None', ecolor = 'black', capsize=2)

        plt.tight_layout()

        if show_plot:
            plt.show()
            plt.close()
        if save_plot:
            # result_plot_dir = 'results/plots'
            # os.makedirs(result_plot_dir, exist_ok=True)
            model_grids_file = result_plot_dir + '/'+ 'index_grids_' + spectra_list_name + '.png'
            plt.savefig(model_grids_file, format='png', dpi=300) #I must save a png image because the eps file does not reproduce well the grids. Don't know why...
            print ('Index-index diagrams for SSP models and data points saved in: ', model_grids_file)
            plt.close()


#wild bootstrap function for uncertainties for ppxf populations by Cappellari
def bootstrap_residuals(model, resid, wild=True):
    """
    https://en.wikipedia.org/wiki/Bootstrapping_(statistics)#Resampling_residuals
    https://en.wikipedia.org/wiki/Bootstrapping_(statistics)#Wild_bootstrap

    Davidson & Flachaire (2008) eq.(12) gives the recommended form
    of the wild bootstrapping probability used here.

    https://doi.org/10.1016/j.jeconom.2008.08.003

    :param spec: model (e.g. best fitting spectrum)
    :param res: residuals (best_fit - observed)
    :param wild: use wild bootstrap to allow for variable errors
    :return: new model with bootstrapped residuals

    """
    if wild:    # Wild Bootstrapping: generates -resid or resid with prob=1/2
        eps = resid*(2*np.random.randint(2, size=resid.size) - 1)
    else:       # Standard Bootstrapping: random selection with repetition
        eps = np.random.choice(resid, size=resid.size)

    return model + eps



# simple function to find the nearest FWHM_gal value corresponding to the emission line funcion, in case of kinematics with gas and variable FWHM_gal (i.e. constant resolving power R).
def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx



# Find a specific SSP template for ppxf kinematics module when you use the two fit stellar components
def pick_ssp_template(desired_age, desired_metal, ages, metals, templates):
    """
    Find the nearest SSP template to the age and metallicity values provided by
    the user for the stellar and gas kinematics task when the two stellar component fit is
    activated.
    Return (model, i_closest, j_closest)
        where model is the spectral template in pPXF standard (n_wave,)
    """
    # If the age grid is in log10:
    # i_closest = np.argmin(np.abs(np.log10(ages) - np.log10(desired_age)))

    # If the age are linear (Gyr) values:
    i_closest = np.argmin(np.abs(ages - desired_age))
    j_closest = np.argmin(np.abs(metals - desired_metal))

    # Retrieving the nearest age and metallicity values for the SSP models
    best_age   = ages[i_closest]
    best_metal = metals[j_closest]

    # Retrieving the corresponding template
    model = templates[:, i_closest, j_closest]

    # Checking if the difference between the desired age and metallicity values is large:
    age_diff   = abs(best_age - desired_age)
    metal_diff = abs(best_metal - desired_metal)
    age_threshold   = 0.1   # adjust the value if needed
    metal_threshold = 0.06  # adjust the value if needed

    if age_diff > age_threshold or metal_diff > metal_threshold:
        msg = (
            f"WARNING: Not found a template for age={desired_age}Gyr, [M/H]={desired_metal}. "
            f"Selected the nearest template with age={best_age:.2f}Gyr, [M/H]={best_metal:.2f}."
        )
        print(msg)

    return model, i_closest, j_closest

#********************** END OF SPECTRA ANALYSIS FUNCTIONS *********************************
#******************************************************************************************
