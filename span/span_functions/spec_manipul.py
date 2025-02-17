#SPectral ANalysis software (SPAN)
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
#*************************** SPECTRA MANIPULATION FUNCTIONS FOR SPAN **********************
#******************************************************************************************
#******************************************************************************************


#Imposting ppxf_util to use the varsmooth function
import ppxf.ppxf_util as util

#Python imports
import numpy as np
import math as mt
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from astropy.convolution import Gaussian1DKernel
from astropy.convolution import convolve
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation

import scipy.stats
from scipy import interpolate
from scipy import ndimage, misc
from scipy.ndimage import gaussian_filter1d
from scipy.signal import butter, filtfilt
from scipy.signal import correlate2d
from scipy.interpolate import interp1d

from skimage.restoration import denoise_wavelet


#*****************************************************************************************************
# 1) RESAMPLE
def resample(wavelength, flux, new_step):
    """
    Resamples a non-linearly spaced spectrum to a linearly spaced wavelength array.

    Parameters:
        wavelength (array-like): Original non-linear wavelength array.
        flux (array-like): Original flux array corresponding to the wavelengths.
        new_step (float): Desired linear step for the resampled wavelength array.

    Returns:
        res_wave (ndarray): Resampled linearly spaced wavelength array.
        res_flux (ndarray): Resampled flux array.
        npoint_resampled (int): Number of points in the resampled spectrum.
    """

    initial_wave = np.min(wavelength)
    final_wave = np.max(wavelength)

    # Create a linearly spaced wavelength array from initial_wave to final_wave
    try:
        res_wave = np.arange(initial_wave, final_wave, new_step)

        # Interpolate the flux values to match the new linearly spaced wavelength grid
        interpfunc = interpolate.interp1d(wavelength, flux, kind='linear', bounds_error=False, fill_value="extrapolate")
        res_flux = interpfunc(res_wave)

    except Exception as e:
        print('Rebinning step too small. Skipping')
        res_wave = wavelength
        res_flux = flux

    npoint_resampled = len(res_wave)

    return res_wave, res_flux, npoint_resampled


#*****************************************************************************************************
# 2) Normalise spectrum
def norm_spec(wavelength, flux, wave_norm, epsilon_avg, flux_to_normalize):

    """

    This function normalises the selected spectrum to a reference wavelength
    defined by the user. In order to avoid to normalise the flux to a peak
    generated by a dead pixel or a sky line falling in the wavelength chosen,
    the normalisation flux is averaged within a wavelength range defined by
    'epsilon_avg', usually taken to be 5 to 10 pix. The normalisation flux
    will be then the average within the flux range included in the interval
    (wave_norm-epsilon_avg) -- (wave_norm+epsilon_avg).
    Input: wavelength array, flux array, reference wavelength to normalise the
    flux, real flux array to normalise.
    Output: Normalised flux array.

    """

    # Masking for selecting the wavelength interval within wave_norm
    mask = (wavelength >= (wave_norm - epsilon_avg)) & (wavelength <= (wave_norm + epsilon_avg))

    # Mean flux values within the window around wave_norm
    avg_flux_ref = np.mean(flux[mask])

    # Normalize the flux
    norm_flux = flux_to_normalize / avg_flux_ref

    return norm_flux


#*****************************************************************************************************
# 3) Degradation of the spectrum
def degrade(wavelength, flux, original_resolution, final_resolution, verbose):

    """

    This function reduces the resolution of the selected spectrum from the
    original resolution R to a user defined resolution R. In order to achieve
    this, the spectrum will be log rebinned, reduced in resolution R and then
    linear rebinned again.
    Input: wavelength array, flux array, original resolution in R, final
    resolution in R, verbose mode (True or False).
    Output: wavelength array, new flux array.

    """

    # Speed of light in km/s
    c = 299792.458

    # Set up the initial parameters
    initial_wave = wavelength[0]
    final_wave = wavelength[-1]
    initial_step = wavelength[1] - wavelength[0]

    # Calculate the velocity per pixel for the initial wavelength
    vel_pix = c / initial_wave * initial_step

    # 1) Convert to log sampling
    log_wave = [initial_wave]
    while log_wave[-1] < final_wave:
        new_step = (vel_pix / c) * log_wave[-1]
        log_wave.append(log_wave[-1] + new_step)

    log_wave = np.array(log_wave)  # Convert to a numpy array for efficient computation

    # Verbose mode logging
    if verbose:
        print('Resampling to log...')
        print(f'New step at lambda initial {log_wave[0]} : {log_wave[1] - log_wave[0]} nm')
        print(f'New step at lambda final {log_wave[-1]} : {log_wave[-1] - log_wave[-2]} nm')

    # Interpolate the flux on the new log wavelength scale
    interpfunc = interpolate.interp1d(wavelength, flux, kind='linear', fill_value='extrapolate')
    log_flux = interpfunc(log_wave)

    # Calculate original and final resolutions in terms of sigma
    fwhm_to_sigma = 2.3548
    sigma_original = (c / original_resolution) / fwhm_to_sigma
    sigma_final = (c / final_resolution) / fwhm_to_sigma

    # Verbose mode logging
    if verbose:
        print(f'Original resolution in sigma vel: {sigma_original}')
        print(f'Final resolution in sigma vel: {sigma_final}')
        print(f'Velocity per pixel: {vel_pix}')

    # 2) Degrade resolution if the final resolution is lower
    if sigma_final > sigma_original:
        sigma_to_broad = mt.sqrt(sigma_final**2 - sigma_original**2)
        gauss_stdev_pix = sigma_to_broad / vel_pix

        if verbose:
            print(f'Sigma to broaden the spectrum: {sigma_to_broad}')
            print(f'Gaussian kernel sigma to apply (in pixels): {gauss_stdev_pix}')

        # Convolve with a Gaussian kernel
        kernel = Gaussian1DKernel(gauss_stdev_pix)
        log_degraded_flux = convolve(log_flux, kernel)
    else:
        print('Final resolution is greater than or equal to the original; no degradation applied.')
        return wavelength, flux

    # 3) Resample to the original wavelength grid
    interpfunc1 = interpolate.interp1d(log_wave, log_degraded_flux, kind='linear', fill_value='extrapolate')
    degraded_flux = interpfunc1(wavelength)

    return wavelength, degraded_flux


#*****************************************************************************************************
# 4) BIS: DEGRADATION OF THE SPECTRA TO DELTA LAMBDA
def degrade_lambda(wavelength, flux, original_resolution_lambda, final_resolution_lambda):

    """

    This function reduces the resolution of the selected spectrum from
    the original resolution delta lambda (FWHM) to a user defined
    resolution delta lambda (FWHM).
    Input: wavelength array, flux array, original resolution in FWHM
    (Angstrom), final resolution in FWHM (in Angstrom).
    Output: wavelength array, new flux array.

    """

    fwhm_to_sigma = 2.3548
    wave_components = len(wavelength)
    step = wavelength[1]-wavelength[0] # STEP SUPPOSED LINEAR!

    original_resolution_lambda_nm = original_resolution_lambda/10. #converting to nm!
    final_resolution_lambda_nm = final_resolution_lambda/10. #converting to nm!


    print ('Original resolution: ', original_resolution_lambda, 'A')
    print ('Final resolution in sigma vel: ', final_resolution_lambda , 'A')


    #2) Degrading the resolution only if the resolution selected is smaller than the initial
    if ((final_resolution_lambda_nm - original_resolution_lambda_nm)> 0):
        real_value_to_broad = mt.sqrt(final_resolution_lambda_nm**2-original_resolution_lambda_nm**2)

        gauss_fwhm_pix = real_value_to_broad/step
        gauss_stdev_pix = gauss_fwhm_pix/fwhm_to_sigma

        print ('Sigma to broad the spectra: ', real_value_to_broad*10, 'A')
        print ('Gaussian kernel sigma to apply (in pixels): ', gauss_stdev_pix)

        kernel = Gaussian1DKernel(gauss_stdev_pix)
        degraded_flux = convolve(flux, kernel)

    else:
        print('Final resolution is less than the original: cannot degrade the spectrum. Doing nothing.')
        return wavelength, flux

    return wavelength, degraded_flux


#*************************************************************************************************
# 5) degrade from R to delta lambda
def degradeRtoFWHM(wavelength, flux, R_resolution, FWHM_resolution):

    """

    This function reduces the resolution of the selected spectrum from
    the constant R resolution to a user defined resolution delta lambda (FWHM).
    Input: wavelength array, flux array, original resolution in R,
    final resolution in delta lambda FWHM (in Angstrom).
    Output: wavelength array, new flux array.

    """

    fwhm_to_sigma = 2.3548

    # Convert FWHM resolution to nanometers
    FWHM_resolution_nm = FWHM_resolution / 10.0

    # Calculate FWHM of original resolution in nm
    R_resolution_FWHM_nm = wavelength / R_resolution

    # Set the target FWHM for the final resolution as a constant array
    final_resolution_FWHM_nm = np.full_like(R_resolution_FWHM_nm, FWHM_resolution_nm, dtype=float)

    # Compute the required broadening values
    try:
        # Calculate broadening in FWHM
        real_value_to_broad = np.sqrt(final_resolution_FWHM_nm**2 - R_resolution_FWHM_nm**2)
        gauss_sigma = real_value_to_broad / fwhm_to_sigma

        # Perform variable sigma convolution using ppxf.util varsmooth
        degraded_flux = util.varsmooth(wavelength, flux, gauss_sigma, xout=None, oversample=1)
        return wavelength, degraded_flux
    except ValueError:
        print('You want to improve the resolution? That''s impossible! Skypping...')
        return wavelength, flux


#*************************************************************************************************
# 6) TRIS: DEGRADATION OF THE SPECTRA FOR LICK INDICES MEASUREMENT
def degrade_to_lick(wavelength, flux, original_resolution, res_delta_lambda):

    """

    This function reduces the resolution of the selected spectrum to match
    that of the Lick/IDS line-strength indices, both for delta lambda constant
    resolution spectrum or to R constant resolution spectrum.
    Input: wavelength array, flux array, original resolution, in R or
    FWHM (Angstrom), bool type of resolution of the input spectrum (True
    for delta lambda constant, False to R constant).
    Output: wavelength array, new flux array.

    """

    fwhm_to_sigma = 2.3548
    wave_components = len(wavelength)
    step = wavelength[1]-wavelength[0] # STEP SUPPOSED LINEAR!

    #IF THE RESOLUTION GIVEN IS IN FWHM!
    if res_delta_lambda:
        original_resolution_lambda_nm = original_resolution/10. #converting to nm!
        final_resolution_lambda_nm_lick = 8.4/10. #converting to nm!


        print ('Original resolution: ', original_resolution, 'A')

        #2) Degrading the resolution only if the resolution selected is smaller than the initial
        if ((final_resolution_lambda_nm_lick - original_resolution_lambda_nm)> 0):

            real_value_to_broad = mt.sqrt(final_resolution_lambda_nm_lick**2-original_resolution_lambda_nm**2)

            gauss_fwhm_pix = real_value_to_broad/step
            gauss_stdev_pix = gauss_fwhm_pix/fwhm_to_sigma

            print ('Sigma to broad the spectra: ', real_value_to_broad*10, 'A')
            print ('Gaussian kernel sigma to apply (in pixels): ', gauss_stdev_pix)

            kernel = Gaussian1DKernel(gauss_stdev_pix)
            degraded_flux = convolve(flux, kernel)

        else:
            print('Final resolution is less than the original: cannot degrade the spectrum. Doing nothing.')
            degraded_flux = flux
            return wavelength, degraded_flux

        return wavelength, degraded_flux

    #IF THE RESOLUTION GIVEN IS IN R!
    if not res_delta_lambda:
        final_resolution_lambda_nm_lick = 8.4/10. #converting to nm!

        original_resolution_lambda_nm = wavelength/original_resolution #array contenente le risoluzioni in FWHM
        final_resolution_lambda_nm = np.full_like(original_resolution_lambda_nm, final_resolution_lambda_nm_lick, dtype=float)#Array with the same size containing the Lick/IDS rresolution, in FEHM (8.4 A)
        real_value_to_broad = np.zeros_like(wavelength)
        degraded_flux = np.zeros_like(flux)

        for i in range (len(wavelength)):

            real_value_to_broad[i] = mt.sqrt(final_resolution_lambda_nm[i]**2-original_resolution_lambda_nm[i]**2)

        gauss_sigma = real_value_to_broad/fwhm_to_sigma

        #using the varsmooth function of ppxf.util for convolution with variable sigma. Works great!
        degraded_flux = util.varsmooth(wavelength, flux, gauss_sigma, xout=None, oversample=1)
        return wavelength, degraded_flux


#*****************************************************************************************************
# 7) Rough continuum subtraction
def sub_cont(wavelength, flux, operation):

    """

    This function computes a rough continuum estimation of the selected
    spectrum by reducing it to a very low resolution. This rough continuum
    model can be subtracted or divided to the original spectrum.
    Input: wavelength array, flux array, operation to perform to the
    spectrum ('subtract' or 'divide').
    Output: continuum subtracted or divided flux array, continuum model flux array.

    """

    cont_wave, cont_flux = degrade(wavelength, flux, 5000., 50., False) # degrade the spectra to a fixed value
    if operation == 'divide':
        norm_flux = flux/cont_flux #I divide the continuum model to the flux
    if operation == 'subtract':
        norm_flux = flux - cont_flux #I subtract the continuum model to the flux
    return norm_flux, cont_flux


#*************************************************************************************************
# 8) Continuum fitting with polynomials and using the mask function
def continuum(wavelength, flux, want_to_maks, mask_ranges, poly_degree, math_operation, with_plots):

    """

    This function computes a more accurate polynomial continuum model
    of the selected spectrum. The continuum model can be subtracted or
    divided to the original spectrum.
    Input: wavelength array, flux array, bool whether to apply or not
    the mask, array of masked pixels, polynomial degree to fit, operation
    to perform to the spectrum ('subtract' or 'divide'), bool whether show
    or not the plot.
    Output: continuum subtracted or divided flux array, continuum model flux array.

    """

    if want_to_maks:
        mask = mask_spectrum(wavelength, mask_ranges)
        wavelength_masked = wavelength[mask]
        flux_masked = flux[mask]
    else:
        wavelength_masked = wavelength
        flux_masked = flux

    coefficients = np.polyfit(wavelength_masked, flux_masked, poly_degree)
    continuum_model = np.polyval(coefficients, wavelength)

    if math_operation == 'subtract':
        new_flux = flux - continuum_model
    elif math_operation == 'divide':
        new_flux = flux/continuum_model

    if with_plots:
        # Show the results with masked regions
        fig, ax = plt.subplots(figsize=(10, 6))

        # Plotting the original spectrum
        ax.plot(wavelength, flux, label='Original spectrum')

        # Plotting the continuum model
        ax.plot(wavelength, continuum_model, label='Continuum model')

        # Plotting the continuum corrected spectrum
        ax.plot(wavelength, new_flux, label='Corrected spectrum')

        # Showing the masked regions with axvspan
        if want_to_maks:
            for mask_range in mask_ranges:
                ax.axvspan(mask_range[0], mask_range[1], color='gray', alpha=0.5)

        ax.legend()
        ax.set_xlabel('Wavelength')
        ax.set_ylabel('Flux')
        plt.title('Continuum fitting')
        plt.show()
        plt.close()
    return new_flux, continuum_model


#*************************************************************************************************
# 9) Mask the spectrum
def mask_spectrum(wavelength, mask_ranges):

    """

    This function computes a pixel mask to apply to the 'continuum' function above.
    Input: wavelength array, array of masked pixels.
    Output: array of masked wavelength regions.

    """

    mask = np.ones_like(wavelength, dtype=bool)

    for mask_range in mask_ranges:
        mask = mask & ((wavelength < mask_range[0]) | (wavelength > mask_range[1]))

    return mask


#*************************************************************************************************
# 10) Sigma clipping
def sigma_clip(wavelength, flux, clipping, resolution, sigma_vel):

    """

    This function computes a dynamic sigma clipping of the flux of the
    selected spectrum in order to delete spikes. The algorithm requires
    the the mean real velocity dispersion which sets an accurate frequency
    cut-off between spike-spurious and real features. For velocity dispersions
    greater than 100 km/s, the functions will performs a pre-cleaning of the highest
    frequency flux values. If your resolution in sigma (velscale) is comparable to
    the real velocity dispersion of the spectrum, you should input a sigma_vel value
    smaller than 100 km/s in order to deactivate the pre-cleaning and perform only the cleaning,.
    Input: wavelength array, flux array, sigma clipping factor, approximate
    resolution in R (even with constant delta lambda, just input a mean
    resolution value in R) and the velocity dispersion (in km/s).
    Output: wavelength array of the spectrum, cleaned flux array of the spectrum.

    """

    c = 299792.458  # Speed of light in km/s
    npoints = len(wavelength)
    base_smooth_width = 41  # Base width for dynamic smoothing
    original_flux = flux  # Store original flux for reference
    clipped_flux = np.zeros(npoints)

    # Define parameters for dynamic smoothing based on velocity dispersion
    if sigma_vel == 0:
        dynamic_size = 41
        sigma_coeff = 0
    else:
        sigma_instrum = c / resolution
        sigma_coeff = sigma_vel / sigma_instrum
        dynamic_size = round(base_smooth_width * sigma_coeff)

        # Adjust window size for small or large sigma_coeff
        if sigma_coeff < 1:
            dynamic_size = int(41)
            window_width_small = 1
        else:
            small_frequency_size = 51
            window_width_small = small_frequency_size

    window_width_dyn = int(dynamic_size)

    # Ensure window width is odd for proper smoothing
    if window_width_small % 2 == 0:
        window_width_small += 1

    clipping_dyn = clipping  # Dynamic clipping threshold
    clipping_small = 2.5  # Clipping for initial small frequency filtering

    elements_flux = len(flux)
    n_clipped = 0  # Counter for clipped points
    sigma_flux_ext_small = np.zeros(window_width_small)
    sigma_flux_ext_dyn = np.zeros(window_width_dyn)

    print('Velocity dispersion:', sigma_vel)
    print('Smooth windows:', window_width_small, window_width_dyn)

    # Step size between wavelength points
    step1 = wavelength[1] - wavelength[0]
    step2 = wavelength[-1] - wavelength[-2]
    epsilon = 1e-4  # Small threshold for comparisons

    # Pre-cleaning stage for high velocity dispersion values
    if sigma_coeff > 1.5:
        threshold = 50
        print('Performing pre-cleaning with smooth size:', window_width_small)
        counter = 1
        while threshold > 0:
            n_clipped = 0  # Reset clipped counter
            smooth_flux = scipy.ndimage.filters.uniform_filter(flux, size=window_width_small)

            # Process each flux point
            for k in range(elements_flux):

                # Apply clipping in central window, ignoring edges
                if window_width_small <= k <= (elements_flux - 1 - window_width_small):
                    sigma_flux_ext_small[:] = flux[int(k - (window_width_small - 1) / 2): int(k + (window_width_small + 1) / 2)]
                    sigma = np.std(sigma_flux_ext_small)

                    # Apply clipping threshold
                    if abs(smooth_flux[k] - flux[k]) <= sigma * clipping_small:
                        clipped_flux[k] = flux[k]
                    else:
                        clipped_flux[k] = smooth_flux[k]
                        n_clipped += 1

                # Preserve flux at edges
                else:
                    clipped_flux[k] = flux[k]

            print('Iteration N.', counter, 'Clipped points:', n_clipped)
            flux = clipped_flux  # Update flux with clipped values
            counter += 1
            threshold = n_clipped  # Continue until no more points are clipped

    # Dynamic sigma clipping using adjusted window size
    print('\nPerforming cleaning with dynamic smooth size:', window_width_dyn)
    counter = 1
    threshold = 50

    while threshold > 0:
        n_clipped = 0
        smooth_flux = scipy.ndimage.filters.uniform_filter(flux, size=window_width_dyn)

        # Process each flux point for dynamic smoothing
        for k in range(elements_flux):

            # Ignore edges and apply clipping in central window
            if window_width_dyn <= k <= (elements_flux - window_width_dyn - 1):
                sigma_flux_ext_dyn[:] = flux[int(k - (window_width_dyn - 1) / 2): int(k + (window_width_dyn + 1) / 2)]
                sigma = np.std(sigma_flux_ext_dyn)

                # Apply dynamic clipping threshold
                if abs(smooth_flux[k] - flux[k]) <= sigma * clipping_dyn:
                    clipped_flux[k] = flux[k]
                else:
                    clipped_flux[k] = smooth_flux[k]
                    n_clipped += 1

            # Preserve flux at edges
            else:
                clipped_flux[k] = flux[k]

        print('Iteration N.', counter, 'Clipped points:', n_clipped)
        flux = clipped_flux  # Update flux with clipped values
        counter += 1
        threshold = n_clipped  # Continue until no more points are clipped

    # Return cleaned spectrum
    new_flux = flux
    return wavelength, new_flux


#*************************************************************************************************
# 11) Log rebin
def log_rebin(wavelength, flux, velscale):

    """

    This function computes a log rebinning of the selected spectrum based
    on the user input velscale.
    IMPORTANT: SPAN automatically linearise the input spectra once they are
    loaded into the listbox. However, this function is not required to perform
    tasks that require the log rebinned spectra (e.g. ppxf). For all the situations
    where a log rebinned spectrum is required, SPAN will automatically convert to the log
    space, will perform the task and then will reconvert to the linear pixel space.
    Input: wavelength array, flux array, velocity scale to rebin (in km/s).
    Output: wavelength array of the log-spaced spectrum, flux array of the log-spaced spectrum.
           IMPORTANT: the output wavelength scale will be expressed in nm,
                      only the step will be log-spaced.

    """

    c = 299792.458
    wave_components = len(wavelength)
    log_wave = []
    initial_wave = wavelength[0]
    final_wave = wavelength[wave_components-1]
    initial_step = wavelength[1]- wavelength[0]
    log_wave.append(wavelength[0])
    tmp_wave = wavelength[0]
    i = 0

    #if velscale input value == 0, I use the initial step
    if velscale == 0:
        velscale = (c*initial_step)/wavelength[0]


    # 1) converting to log resample so sigma = c/R = const
    while tmp_wave <= final_wave:
        new_step = (velscale/c)*log_wave[i]
        log_wave.append(log_wave[i]+new_step)
        tmp_wave = tmp_wave +new_step
        i = i+1

     #Now calculate the interpolated flux at the res_wave points
    log_wave = np.array(log_wave)
    interpfunc = interpolate.interp1d(wavelength, flux, kind = 'linear',fill_value='extrapolate')
    log_flux = (interpfunc(log_wave))

    return log_wave, log_flux


#*************************************************************************************************
# 12) Doppler correction
def dopcor(wavelength, flux, value, is_velocity):

    """

    This function computes the doppler correction of the selected
    spectrum and gives the rest-frame spectrum.
    Input: wavelength array, flux array, recession velocity of the spectrum (in km/s) or redshift (z).
    Output: wavelength array of the rest-frame spectrum.

    """
    if is_velocity:
        print('Doppler correction with velocity = ', value)
        c = 299792.458
        z = value/c
        corr_wave = wavelength/(1+z)
    else:
        print('Redshift correction with z = ', value)
        z = value
        corr_wave = wavelength/(1+z)

    return corr_wave, flux


#*************************************************************************************************
# 13) Sigma broadening: broad the spectra to a user defined sigma value (km/s)
def sigma_broad(wavelength, flux, sigma_to_broad):

    """

    This function computes the broadening of the spectrum to a user
    specified velocity dispersion. Remember that the final real broadening
    will be the square root of the quadrature sum of the instrumental broadening
    and the broadening value selected.
    Input: wavelength array of the spectrum, flux array of the spectrum, velocity scale to
    add to the spectrum.
    Output: velocity broadened flux array of the spectrum.

    """

    c = 299792.458
    step = wavelength[1]-wavelength[0]
    sigma_pix = c/wavelength[0]*step

    if sigma_to_broad == 0:
        flux_gauss = flux
        return flux_gauss

    #rebinning to sigma = cost
    log_wave, log_flux = log_rebin(wavelength,flux, sigma_pix)

    gauss_stdev_pix = sigma_to_broad/sigma_pix
    kernel = Gaussian1DKernel(gauss_stdev_pix)
    log_flux_gauss = convolve(log_flux, kernel)

    #back to the original, linear sample
    interpfunc1 = interpolate.interp1d(log_wave, log_flux_gauss, kind = 'linear')
    flux_gauss = (interpfunc1(wavelength))

    return flux_gauss


#*************************************************************************************************
# 14) Add noise
def add_noise(wavelength, flux, snr):

    """

    This function adds a poissonian noise to the spectrum corresponding to
    the user input sigmal-to-noise.
    Remember that the final real signal-to-noise will be the square root of
    the quadrature sum of the original signal-to-noise and the signal-to-noise
    value selected.
    Input: wavelength array of the spectrum, flux array of the spectrum,
           signal-to-noise to add to the spectrum.
    Output: noise added flux array of the spectrum.

    """

    #normalise the continuum. This gives me the possibility to have the new spectrum at level = 1, usefull for adding noise
    norm_flux, cont_flux = sub_cont(wavelength, flux, 'divide')
    npoints = len(wavelength)

    #generate the noise array
    noise_array = np.random.standard_normal((npoints,))
    scale = 1./snr
    noise_array_scaled = noise_array * scale
    noisy_norm_flux = norm_flux + noise_array_scaled

    #give back the continuum shape
    noisy_flux = noisy_norm_flux * cont_flux
    return noisy_flux


#*************************************************************************************************
# 15) Simple box window moving average
def mov_avg(flux, window_size):

    """

    This function computes a simple box-window moving average of the
    selected spectrum.
    Input: flux array of the spectrum, size (in pix) of the box-window
    Output: moving averaged flux array of the spectrum.

    """

    avg_flux = np.convolve(flux, np.ones((window_size,))/window_size, mode='same') # for now I decided to not handle the edges with the keyword 'same'
    return avg_flux


#*************************************************************************************************
# 16) Gaussian box moving average
def mov_avg_gauss(wavelength, flux, sigma):

    """

    This function computes a gaussian-kernel moving average of the
    selected spectrum.
    Input: flux array of the spectrum, size (in pix) of the box-window
    Output: moving averaged flux array of the spectrum.

    """

    try:
        denoised_flux = gaussian_filter1d(flux, sigma=sigma)
        return denoised_flux
    except:
        print ('Error applying the gaussian moving average. Skipping...')
        return flux


#*************************************************************************************************
# 17) Heliocentric calculation and correction on the spectrum
def helio_corr(wavelength, flux, epoch, where, ra_obj, dec_obj):

    """

    This function computes the heliocentric correction for the selected
    spectrum and then apply the correction. It uses the 'EarthLocation'
    function of the Astropy library in order to resolve the name of the
    location insterted by the user and to retrieve the informations of
    the heliocentric correction to apply.
    IMPORTANT: The first time that SPAN calls this function, it requires
    an internet connection, otherwise it raises an error.
    Input: wavelength array of the spectrum, flux array of the spectrum,
    date of the observation, location, RA and Dec. of the observation.
    Output: float heliocentric calculated velocity, heliocentric corrected
            wavelength array, heliocentric corrected flux array of the spectrum.

    """

    location = EarthLocation.of_site(where)
    sc = SkyCoord(ra=ra_obj*u.deg, dec=dec_obj*u.deg)
    heliocorr = sc.radial_velocity_correction('heliocentric', obstime=Time(epoch), location = location)
    heliocorr = heliocorr.to(u.km/u.s)
    location_list = EarthLocation.get_site_names()

    heliocorr = (str(heliocorr))

    #extract the value from heliocorr
    helio_number = 0.
    for heliocorr in heliocorr.split():
        try:
            helio_number = (float(heliocorr))
        except ValueError:
            pass

    #correct_the_spectrum
    new_wavelength, new_flux = dopcor(wavelength, flux, helio_number, True)

    return helio_number, new_wavelength, new_flux


#*************************************************************************************************
# 18) simple cropping function
def crop_spec(wavelength, flux, wave_interval):

    """

    This function performs a simple cropping of the spectrum within
    the wavelength range defined by the user.
    Input: wavelength array of the spectrum, flux array of the spectrum,
    wavelength interval to crop the spectrum (2 component array containing
    the minimum and maximum values of the wavelength range considered).
    Output: cropped wavelength array, cropped flux array of the spectrum.

    """

    wave1 = np.min(wave_interval)
    wave2 = np.max(wave_interval)
    cropped_wave = wavelength[(wavelength >= wave1) & (wavelength <= wave2)]
    cropped_flux = flux[(wavelength >= wave1) & (wavelength <= wave2)]

    return cropped_wave, cropped_flux


#*************************************************************************************************
# 19) wavelets denoising
def wavelet_cleaning(wavelength, flux, sigma, wavelets_layers):

    """

    This function applies a wavelet based denoise to the selected spectrum
    from the ski-image module, using the pywavelets library.
    Input: wavelength array of the spectrum, flux array of the spectrum,
    sigma, number of wavelet layers.
    Output: wavelet corrected flux array of the spectrum.

    """

    #normalising the spectrum
    epsilon_norm = 2
    denoised_flux = denoise_wavelet(flux, sigma=sigma, wavelet='sym5', mode='soft', wavelet_levels=wavelets_layers)

    return denoised_flux


#*************************************************************************************************
# 20) lowpass filter
def lowpass(wavelength, flux, cut_off, order):

    """

    This function applies a Butterworth low-pass filter to filter-out
    the spatially higher frequencies of the selected spectrum, using the
    ski-image module.
    Input: wavelength array, flux array, float frequency cut-off, int order
    of the filter
    Output: filtered flux array of the spectrum.

    """

    try:
        b, a = butter(order, cut_off, btype='lowpass', analog=False)
        denoised_flux = filtfilt(b, a, flux)
        return denoised_flux
    except:
        print ('Error applying the lowpass filter. Skipping...')
        return flux


#*************************************************************************************************
# 21) Bandpass filter
def bandpass(wavelength, flux, lower_cut_off, upper_cut_off, order):

    """

    This function applies a Butterworth band-pass filter to filter-out
    the highest and lowest spatial frequencies of the selected spectrum, using the
    ski-image module.
    Input: wavelength array, flux array, float low frequency cut-off, float high frequency cut-off,
    int order of the filter.
    Output: filtered flux array of the spectrum.

    """

    try:
        b, a = butter(order, [lower_cut_off, upper_cut_off], btype='band', analog=False)
        denoised_flux = filtfilt(b, a, flux)
        return denoised_flux
    except:
        print ('Error applying the bandpass filter. Skipping...')
        return flux


#********************** END OF SPECTRA MANIPULATION FUNCTIONS *****************************
#******************************************************************************************
