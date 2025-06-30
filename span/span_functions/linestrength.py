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
#*************************** LINE-STRENGTH FUNCTIONS FOR SPAN *****************************
#******************************************************************************************
#******************************************************************************************

try:#Local imports
    from FreeSimpleGUI_local import FreeSimpleGUI as sg
    from span_functions import spec_manipul as spman
    from span_functions import system_span as stm

except ModuleNotFoundError: #local import if executed as package
    from ..FreeSimpleGUI_local import FreeSimpleGUI as sg
    #SPAN functions import
    from . import spec_manipul as spman
    from . import system_span as stm

#Python imports
import numpy as np
import math as mt
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from scipy import interpolate
from scipy.interpolate import interp1d
import os


CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(CURRENT_DIR)


#*****************************************************************************************************
# 1) EXTRACT INDEX limits in order to measure the line-strength
def extract_index(wavelength, flux, index):

    """
    This function extracts the wavelength and the flux within the line-strength
    index definition bands, regardless the position of the continuum bands
    (simmetrycal to the line or just to one side)
    Input: wavelength array, flux array, index definition array.
    Output: index wavelength array, index flux array.

    """

    left_band = index[0]
    right_band = index[3]
    left_index = index[4]
    right_index = index[5]
    epsilon = 50 #arbitrary values in A to add to the index window
    index_flux = []
    index_wave = []
    wave_components = len(wavelength)

    #;case 1): Line band in the middle
    if (left_band <= left_index and left_index <= right_band):
        left_band_tmp = left_band-epsilon
        right_band_tmp = right_band+epsilon
        left_index_tmp = left_index

    #;case 2): Continuum bands on the left
    if (left_band <= left_index and left_index >= right_band):
        left_band_tmp = left_band -epsilon
        right_band_tmp = right_index+epsilon
        right_index_tmp = right_band

#;case 3): Continuum bands on the right
    if (left_band >= left_index and left_band <= right_band):
        left_band_tmp = left_index-epsilon
        right_band_tmp = right_band+epsilon
        left_index_tmp = left_band

    for i in range (wave_components):
        if (wavelength[i] >= left_band_tmp and wavelength[i] <= right_band_tmp):
            index_flux.append(flux[i])
            index_wave.append(wavelength[i])
    return index_wave, index_flux


#*****************************************************************************************************
# 2) Generating a pseudo continuum, useful for the Equivalent Width (EW) calculation
def idx_cont(index, wavelength, flux):

    """
    This function computes a syntethic psedudo-continuum level (straigth line) over the line defininition
    of a line-strength index by iterpolating between the mean level of the two pseudo-continuum
    definition bands. THe results is a straigth line between the mean wavelength blue and red bands
    of the index definition.
    Input: array index definintion wavelength, wavelength array, flux array of the spectrum.
    Output: array of the syntethic flux generated between the line definintion of the considered index,
            array of the real flux between the line definintion of the considered index,
            array of the wavelength grid between the line definintion of the considered index,
            array of the mean fluxes of the blue and red index definition band,
            array of the mean wavelength of the blue and red index definition band,
            float of the central wavelength value of the blue index definition band,
            float of the central wavelength value of the red index definition band,
            float mean flux value of the blue index definition band,
            float mean flux value of the red index definition band,



    """

    left_wave_a = index[0]
    left_wave_b = index[1]
    right_wave_a = index[2]
    right_wave_b = index[3]
    index_left_band = index[4]
    index_right_band = index[5]
    nx = len(wavelength)

    cont_flux_left = 0.
    points_left = 0.
    cont_flux_right = 0.
    points_right = 0.

    #find the total flux on the blue and red continuum band
    for i in range(nx):
        if (wavelength[i] >= left_wave_a and wavelength[i] <= left_wave_b):
            cont_flux_left = cont_flux_left + flux[i]
            points_left += 1
        if (wavelength[i] >= right_wave_a and wavelength[i] <= right_wave_b):
            cont_flux_right = cont_flux_right + flux[i]
            points_right += 1


    #finding the average blue and red flux
    avg_left_flux = cont_flux_left/points_left
    avg_right_flux = cont_flux_right/points_right

    #finding the central blue and red band of the continuum
    central_left_wave = (left_wave_a + left_wave_b) /2.
    central_right_wave = (right_wave_a + right_wave_b ) /2.

    #add the data to arrays
    flux_ref_cont = [avg_left_flux, avg_right_flux]
    wave_ref_cont = [central_left_wave, central_right_wave]

    #finding the x points (lambda) where interpolate the continuum
    wave_pseudo_cont = []
    line_flux = []

    for i in range(nx):
        if (wavelength[i] >= index_left_band and wavelength[i] <= index_right_band):
            wave_pseudo_cont.append(wavelength[i])
            line_flux.append(flux[i])

    #doing the interpolation here
    interpfunc = interpolate.interp1d(wave_ref_cont, flux_ref_cont, kind = 'linear',fill_value='extrapolate')
    flux_pseudo_cont = (interpfunc(wave_pseudo_cont))

    return flux_pseudo_cont, line_flux, wave_pseudo_cont, flux_ref_cont, wave_ref_cont, central_left_wave, central_right_wave, avg_left_flux, avg_right_flux


#*****************************************************************************************************
#3) Calculating the EW of an index in A
def eq_width(flux_pseudo_cont, line_flux, step):

    """
    This function calculates the equivalent width in A of the selected spectrum
    for a selected index.
    Input: flux array of the two pseudo-continuum bands,flux array of the line band,
    lambda step of the seledted spectrum, supposed linear.
    Output: float equivalent width value, in Angstrom.

    """

    ew = 0.
    points = len(line_flux)
    for i in range(points):
        ew = ew + (1- line_flux[i]/flux_pseudo_cont[i])*step
    return ew


#*****************************************************************************************************
#4) Calculating the EW of an index in mag
def eq_width_mag(flux_pseudo_cont, line_flux, step,lambda_blue_line, lambda_red_line):

    """
    This function calculates the equivalent width in magnitudes of the selected spectrum
    for a selected index.
    Input: flux array of the two pseudo-continuum bands,flux array of the line band,
    lambda step of the seledted spectrum, supposed linear.
    Output: float equivalent width value, in magnitudes.

    """

    ew_mag = 0.
    points = len(line_flux)
    for i in range(points):
        ew_mag = ew_mag + (line_flux[i]/flux_pseudo_cont[i]*step)
    #Check for negative argument of the log!
    if ew_mag <= 0:
        ew_mag = 999
    else:
        ew_mag = -2.5*mt.log10(1/(lambda_red_line-lambda_blue_line)*ew_mag)
    return ew_mag


#*****************************************************************************************************
#5) Uncertainties of the EW values
def ew_err (index, wavelength, flux, step, flux_pseudo_cont, wave_pseudo_cont, flux_ref_cont, wave_ref_cont):

    """
    This function calculates uncertainties, in Angstrom, of a
    line-strength index of the selected spectrum, by performing MonteCarlo simulations.
    The function generates a straigth line pseudo-continuum covering all the index band
    definition range, creates n pseudo-continuum with the signal-to-noise of the
    real pseudo-continuum regions of the index and measures the equivalent width.
    The standard deviation of the values around the zero is taken as representative of
    the uncertainties of the equivalent width measured in the spectrum, for the same
    signal-to-noise.
    Input: array of index definition wavelength, wavelength array of the selected
    spectrum, flux array of the selected spectrum, delta lambda step supposed linear,
    array of pseudo-continuum flux, 2 component array of the mean flux in each of the
    pseudo-continuum bands, 2 component array of the mean wavelength of the pseudo-continuum bands.
    Output: uncentainty on the equivalent width measurement, flux standard deviation of
    thepseudo cintinuum bands.

    """

    nx = len(wavelength)
    nx_check = len(flux)
    if (nx != nx_check):
        print ('ERROR: wavelength componens different from flux. Something went wrong!')

    #;definying useful arrays and variables
    flux_red_band = []
    flux_blue_band = []
    red_wave = []
    blue_wave = []
    h = 0

    #;extract the flux and lambda arrays from the continuum blue and red bands.
    for i in range(nx):
        if (wavelength[i] >= index[0] and wavelength[i] <= index[1]):
            flux_blue_band.append(flux[i])
            blue_wave.append(wavelength[i])

        if (wavelength[i] >= index[2] and wavelength[i] <= index[3]):
            flux_red_band.append(flux[i])
            red_wave.append(wavelength[i])

    wave_all = []
    flux_all = []

    for i in range(len(blue_wave)):
        wave_all.append(blue_wave[i])
        flux_all.append(flux_blue_band[i])

    for j in range(len(red_wave)):
        wave_all.append(red_wave[j])
        flux_all.append(flux_red_band[j])

    #interpolate to create the a fake continuum where measure the EW
    interpfunc = interpolate.interp1d(wave_ref_cont, flux_ref_cont, kind = 'linear', fill_value='extrapolate')
    interp_all = (interpfunc(wave_all))

    #calculating the residual continuum
    components = len(flux_all)
    residuals = []
    for i in range(components):
        residuals.append(flux_all[i] - interp_all[i])

    sigma_cont_real = np.std(residuals)

    scale = sigma_cont_real
    line_components = len(flux_pseudo_cont)

    number_noisy_cont = 100 #how many syntethics? a lot!
    ews= []


    for k in range(number_noisy_cont):
        noise_array = np.random.standard_normal((line_components,))
        noise_array_scaled = noise_array * scale
        noisy_cont = []

        for i in range(line_components):
            noisy_cont.append(flux_pseudo_cont[i] + noise_array_scaled[i])

        #;measuring the EW
        ew = eq_width(flux_pseudo_cont, noisy_cont, step)
        ews.append(ew)

    error = np.std(ews)
    return error, sigma_cont_real


#*****************************************************************************************************
# 6) reading index file
def read_idx(index_file):

    """
    This function reads the formatted index ASCII file containing the names and the
    definitions of the line-strength indices. See the 'readme_span.txt' file to find out
    how it has to be formatted.
    Input: path and name of the index ASCII file,
    Output: string containing the names of the indices found and nXm array containing the
            wavelength definitions of the indices.

    """

    indices = []
    idx_names = np.loadtxt(index_file, dtype = 'str', delimiter = ' ', max_rows = 1) #reading the first line only, with the header
    indices = np.loadtxt(index_file, comments = '#', skiprows = 1) #reading the other lines
    return idx_names, indices


#*****************************************************************************************************
# 7) Equivalenth width measurements
def ew_measurement(wavelength, flux, index_file, is_usr_idx, want_plot, verbose, calculate_error, save_plot, spec_name, normalize_spec, result_plot_dir):

    """
    This function measures the equivalent width, in Angstrom and magnitudes,
    the signal-to-noise and plots the line-strength index of the selected spectrum.
    Input: wavelength array, flux array, string index file name, bool whether is a
           single used defined index (True) or not (False, bool whether display (True)
           or not (False) the plots, bool whether display the verbose (True) or not (False)
           output in the terminal, bool whether calculate (True) or not (False) the
           uncertainties with MonteCarlo simulations, bool whether save (True) or not (False)
           the plots in PGN images, string name of the spectra or the spectrum, bool whether
           normalise (True) or not (False) the spectrum.
    Output: string containing the name of the index, float equivalent width value measured (in Angstrom),
            float uncertainty (if calculated) of the measured index, float S/N estimated within the
            index definition blue and red bands, float equivalent width value measured (in magnitudes),
            float uncertainty (if calculated, in magnitudes).

    """

    new_step = wavelength[1]-wavelength[0]

    #reading the index file
    if not is_usr_idx: #if index_file is a file containing the index definitions
        id_array, index = read_idx(index_file)
        num_indices = len(id_array)

    else: # if index_file is a numpy array already containing the definition of a usr_idx
        id_array = 'usr_idx'
        num_indices = 1
        index = index_file

    if verbose:
        print ('Number of indices to measure = ', num_indices)

        if not is_usr_idx:
            index_transpose = np.transpose(index) #only for the following visualization
            for i in range(num_indices):
                print (id_array[i], index_transpose[i,:]) #printing the indices names and their wavelengths
            print ('')
        else:
            print(id_array)
            print(index)
            print ('')

    #defining the arrays containing the data
    ew_array = np.zeros(num_indices)
    ew_array_mag = np.zeros(num_indices)
    snr_ew_array = np.zeros(num_indices)
    err_array = np.zeros(num_indices)
    err_array_mag = np.zeros(num_indices)
    wave_limits_spec = np.array([wavelength[0], wavelength[len(wavelength)-1]])
    #and the name of the spectrum, without the path
    spec_name_no_path = os.path.basename(spec_name)

    #************************************* Cycling all over the indices **********************************

        #if I have a list of indices
    if not is_usr_idx and normalize_spec:

        # Arranging the arrays and variables for plotting
        num_good_indices = 0
        num_bad_indices = 0
        good_id_array = []
        good_index = []
        good_ew = []
        good_err = []
        good_snr_per_pix = []
        good_index_wave = []
        good_index_flux = []
        good_band = []
        good_interp_flux_bands = []
        good_interp_lambda = []
        good_interp_cont_flux = []

        #now Cycling for the indices, discarding those not valid because not included in the wavelength range of the spectrum
        for t in range (num_indices):
            lambda_ref_norm = index[4,t]
            min_idx = np.min(index[:,t])
            max_idx = np.max(index[:, t])
            if (min_idx < wave_limits_spec[0] or max_idx > wave_limits_spec[1] or lambda_ref_norm > np.max(wave_limits_spec) or lambda_ref_norm < np.min(wave_limits_spec)):
                print ('Skipping the index: ', id_array[t])
                num_bad_indices += 1
                continue
            else:
                num_good_indices +=1

            # 1)Normalizing the spectra
            epsilon_wave = new_step*10. #just an epsilon value to average the flux for the normalization value
            if normalize_spec:
                norm_flux = spman.norm_spec(wavelength, flux, lambda_ref_norm, epsilon_wave, flux)
            elif not normalize_spec:
                norm_flux = flux

            # 2) Extract the spectral region around the index (+/- 50 A) to speed up the process
            index_wave, index_flux = extract_index(wavelength, norm_flux, index[:,t])

            #3) Extract the pseudo continuum
            interp_cont_flux, line_flux, interp_lambda, flux_ref_cont, lambda_ref_cont, central_left_lambda, central_right_lambda, avg_left_flux, avg_right_flux = idx_cont(index[:,t], index_wave, index_flux)

            #4) Determining the EW of the index
            ew = eq_width(interp_cont_flux, line_flux, new_step)

            #4) BIS Determining the EW of the index in mag
            ew_mag = eq_width_mag(interp_cont_flux, line_flux, new_step, index[4,t],index[5,t])

            #5) Calculate the error via MonteCarlo simulation, if you want
            if calculate_error:
                err, sigma_cont = ew_err(index[:,t], index_wave, index_flux, new_step, interp_cont_flux, interp_lambda, flux_ref_cont, lambda_ref_cont)
            else:
                err = 0.

            #trasform in A
            ew = ew
            err = err

            #calculating the errors in magnitudes
            err_mag= 0.434*abs(err/ew)

            #fill the arrays
            ew_array[t] = ew
            err_array[t] = err

            ew_array_mag[t] = ew_mag
            err_array_mag[t] = err_mag

            # Calculate the snr
            if calculate_error:
                snr_per_pix = (avg_left_flux + avg_right_flux)/(2*sigma_cont)
                pix_per_a = 1/(new_step)
                snr_per_a = snr_per_pix*mt.sqrt(pix_per_a)

                #fill the array
                snr_ew_array[t] = snr_per_pix
            else:
                snr_ew_array[t] = 0

            band = []
            for i in range (len(wavelength)):
                if (wavelength[i] >= central_left_lambda and wavelength[i] <= central_right_lambda):
                    band.append(wavelength[i])
            interpfunc = interpolate.interp1d(lambda_ref_cont, flux_ref_cont, kind = 'linear')
            interp_flux_bands = (interpfunc(band))

            #storing the values of only the valid measured indices and only if I want to display the plots in the 'Preview result' mode
            if(want_plot or save_plot):
                good_id_array.append(id_array[t])
                good_ew.append(ew)
                good_err.append(err)
                good_index.append(index[:,t])
                good_snr_per_pix.append(snr_per_pix)
                good_index_wave.append(index_wave)
                good_index_flux.append(index_flux)

                good_band.append(band)
                good_interp_flux_bands.append(interp_flux_bands)
                good_interp_lambda.append(interp_lambda)
                good_interp_cont_flux.append(interp_cont_flux)


        # In the 'Preview result' mode I want to arrange all the plots of the valid indices in one matplot window. In order to do this, I need to perform some make-up.
        if(want_plot or save_plot):

            #converting to numpy array and use only the good indices
            good_index = np.array(good_index).T

            try:
                cols = int(np.ceil(np.sqrt(num_good_indices)))
                rows = int(np.ceil(num_good_indices / cols))

                #create the figure and the subplots
                fig, axes = plt.subplots(rows, cols, figsize=(13, 7))

                # Flatten the axes array
                axes = axes.flatten()
                #dynamic font size for the titles
                font_size = max(12 - num_good_indices // 3, 8)
            except Exception as e:
                print('')
                print ('No valid indices to plot')
                print('')
                want_plot = False

            #Cycling for all the good indices measured before
            for t in range (num_good_indices):

                #set the y limits for the plots
                ylim_low = np.mean(interp_flux_bands)-0.6
                ylim_high = np.mean(interp_flux_bands)+0.8

                yeps = 0.1
                ew_string = str(round(good_ew[t],2))
                err_string = str(round(good_err[t],2))
                snr_ew_string = str(round(good_snr_per_pix[t],0))

                #if the bands are symmetric with respect to the line band or not:
                if (good_index[3, t] < good_index[5, t]):
                    axes[t].set_title(good_id_array[t] + ' ' + ew_string + r'$\pm$' + err_string + r' $\AA$.  SNR =' + snr_ew_string, fontsize = font_size)
                    axes[t].plot(good_index_wave[t], good_index_flux[t], linewidth=0.5, color = 'green')
                    axes[t].set_xlim(good_index[0, t]-5., good_index[5, t]+5.)
                    axes[t].set_ylim (ylim_low, ylim_high)
                    axes[t].set_xlabel('Wavelength A', fontsize = 9)
                    axes[t].set_ylabel('Flux', fontsize = 9)
                    axes[t].tick_params(axis = 'both', labelsize = 9)

                else:
                    axes[t].set_title(good_id_array[t] + ' ' + ew_string + r'$\pm$' + err_string + r' $\AA$.  SNR =' + snr_ew_string, fontsize = font_size)
                    axes[t].plot(good_index_wave[t], good_index_flux[t], linewidth=0.5, color = 'green')
                    axes[t].set_xlim(good_index[0, t]-5., good_index[3, t]+5.)
                    axes[t].set_ylim(ylim_low, ylim_high)
                    axes[t].set_xlabel('Wavelength A', fontsize = 9)
                    axes[t].set_ylabel('Flux', fontsize = 9)
                    axes[t].tick_params(axis = 'both', labelsize = 9)

                #pseudo continuum
                axes[t].plot(good_band[t], good_interp_flux_bands[t], linewidth = 0.5, color = 'black')
                axes[t].plot(good_interp_lambda[t], good_interp_cont_flux[t], linewidth = 1, color = 'blue')

                #polygons
                x_polygon_bband = [good_index[0, t], good_index[1,t], good_index[1,t], good_index[0,t], good_index[0,t]]
                y_polygon_bband = [0.4+yeps, 0.4+yeps, 1.8-yeps, 1.8-yeps, 0.4+yeps]
                axes[t].fill(x_polygon_bband,y_polygon_bband, 'blue')

                x_polygon_rband = [good_index[2,t], good_index[3,t], good_index[3,t], good_index[2,t], good_index[2,t]]
                y_polygon_rband = [0.4+yeps, 0.4+yeps, 1.8-yeps, 1.8-yeps, 0.4+yeps]
                axes[t].fill(x_polygon_rband,y_polygon_rband, 'red')

                x_polygon_line = [good_index[4,t], good_index[5,t], good_index[5,t], good_index[4,t], good_index[4,t]]
                y_polygon_line = [0.4+yeps, 0.4+yeps, 1.8-yeps, 1.8-yeps, 0.4+yeps]
                axes[t].fill(x_polygon_line,y_polygon_line, 'gray')

                if t == num_good_indices - 1:
                # removing the empty plots at the end of the grid, if any
                    for ax in axes:
                        if not ax.has_data():
                            fig.delaxes(ax)

                    # finally plotting
                    plt.tight_layout()

                    if want_plot:
                        plt.show()
                        plt.close()
                    else:
                        plt.savefig(result_plot_dir + '/'+spec_name_no_path + '_LS.png', format='png', dpi=300)
                        plt.close()

        return id_array, ew_array, err_array, snr_ew_array, ew_array_mag, err_array_mag


        # if I have just a single index
    else:
        id_array = 'usr_idx'
        if normalize_spec:

            #normalization wavelength
            lambda_ref_norm = index[4]

            # 1)Normalising the spectra
            epsilon_wave = new_step*10. #just an epsilon value to average the flux for the normalization value
            norm_flux = spman.norm_spec(wavelength, flux, lambda_ref_norm, epsilon_wave, flux)

        if not normalize_spec:
            norm_flux = flux


        # 2) Extract the spectral region around the index (+/- 50 A) to speed up the process
        index_wave, index_flux = extract_index(wavelength, norm_flux, index)

        #3) Extract the pseudo continuum
        interp_cont_flux, line_flux, interp_lambda, flux_ref_cont, lambda_ref_cont, central_left_lambda, central_right_lambda, avg_left_flux, avg_right_flux = idx_cont(index, index_wave, index_flux)

        #4) Determining the EW of the index
        ew = eq_width(interp_cont_flux, line_flux, new_step)


        #4) BIS Determining the EW of the index in MAG
        ew_mag = eq_width_mag(interp_cont_flux, line_flux, new_step, index[4],index[5])

        #5) Calculate the error via MonteCarlo simulation
        if calculate_error:
            err, sigma_cont = ew_err(index, index_wave, index_flux, new_step, interp_cont_flux, interp_lambda, flux_ref_cont, lambda_ref_cont)
        else:
            err = 0.

        #fill the vectors
        ew = ew
        err = err
        err_mag= 0.434*abs(err/ew)

        # Calculate the snr
        if calculate_error:
            snr_per_pix = (avg_left_flux + avg_right_flux)/(2*sigma_cont)
            pix_per_a = 1/(new_step)
            snr_per_a = snr_per_pix*mt.sqrt(pix_per_a)
        else:
            snr_per_pix = 0
            snr_per_a = 0

        # doing the plots for the single index and only if I want
        if(want_plot or save_plot):
            band = []
            for i in range (len(wavelength)):
                if (wavelength[i] >= central_left_lambda and wavelength[i] <= central_right_lambda):
                    band.append(wavelength[i])
            interpfunc = interpolate.interp1d(lambda_ref_cont, flux_ref_cont, kind = 'linear')
            interp_flux_bands = (interpfunc(band))

            #set the y limits for the plots
            ylim_low = np.mean(interp_flux_bands)-0.6
            ylim_high = np.mean(interp_flux_bands)+0.8

            yeps = 0.1
            ew_string = str(round(ew,3))
            err_string = str(round(err,3))
            snr_ew_string = str(round(snr_per_pix,0))

            if (index[3] < index[5]):
                plt.title('EW usr index ' + ew_string + r'$\pm$' + err_string + r' $\AA$.  SNR =' + snr_ew_string)
                plt.plot(index_wave, index_flux, linewidth=0.5, color = 'green')
                plt.xlim(index[0]-5., index[5]+5.)
                plt.ylim(ylim_low, ylim_high)
                plt.xlabel('Wavelength A', fontsize = 9)
                plt.ylabel('Flux', fontsize = 9)
                plt.tick_params(axis = 'both', labelsize = 9)

            else:
                plt.title('EW usr index ' + ew_string + r'$\pm$' + err_string + r' $\AA$.  SNR =' + snr_ew_string)
                plt.plot(index_wave, index_flux, linewidth=0.5, color = 'green')
                plt.xlim(index[0]-5., index[3]+5.)
                plt.ylim(ylim_low, ylim_high)
                plt.xlabel('Wavelength A', fontsize = 9)
                plt.ylabel('Flux', fontsize = 9)
                plt.tick_params(axis = 'both', labelsize = 9)

            #pseudo continuum
            plt.plot(band, interp_flux_bands, linewidth = 0.5, color = 'black')
            plt.plot(interp_lambda, interp_cont_flux, linewidth = 1, color = 'blue')

            #polygons
            x_polygon_bband = [index[0], index[1], index[1], index[0], index[0]]
            y_polygon_bband = [ylim_low+yeps, ylim_low+yeps, ylim_high-yeps, ylim_high-yeps, ylim_low+yeps]
            plt.fill(x_polygon_bband,y_polygon_bband, 'blue')

            x_polygon_rband = [index[2], index[3], index[3], index[2], index[2]]
            y_polygon_rband = [ylim_low+yeps, ylim_low+yeps, ylim_high-yeps, ylim_high-yeps, ylim_low+yeps]
            plt.fill(x_polygon_rband,y_polygon_rband, 'red')

            x_polygon_line = [index[4], index[5], index[5], index[4], index[4]]
            y_polygon_line = [ylim_low+yeps, ylim_low+yeps, ylim_high-yeps, ylim_high-yeps, ylim_low+yeps]
            plt.fill(x_polygon_line,y_polygon_line, 'gray')

            if want_plot:
                plt.show()
                plt.close()
            else:
                # result_plot_dir = 'results/plots'
                # os.makedirs(result_plot_dir, exist_ok=True)
                plt.savefig(result_plot_dir + '/'+spec_name_no_path + '_' + id_array +'.png', format='png', dpi=300)
                plt.close()


    return id_array, ew, err, snr_per_pix, ew_mag, err_mag


#*****************************************************************************************************
#8) Calculating the velocity dispersion coefficients to correct the EW due to velocity dispersion broadening
def sigma_coeff (spectra_file, index_file, lambda_units, is_usr_idx, want_plot, smooth_value, save_plot, result_plot_dir):

    """
    This function calculates the velocity dispersion correction coefficients as a
    function of the velocity dispersion broadening from a sample of templates or
    stellar spectra with zero velocity dispersion. It takes the sample spectra, broaden up to 400 km/s
    with 50 km/s steps, calculates the EW of the selected index for any spectra and any broadening
    value and returns the spline coefficients that best interpolate the mean EW of the
    broadened semple for a fixed index and the relative uncertainties.
    The correction coefficients from the spline interpolation of the broadened index values
    measured from the sample spectra are: C0, C1, C2, C3, following the method and
    nomenclature of Trager et al. 1998 for the Lick/IDS indices.
    Input: string of the spectra list containing the path and names of template (or stellar) spectra
           from which determine the correction coefficients, string ASCII file containing the
           names and definitions of the line-strenght indices for which calculate the broadening
           coefficients, string wavelength units scale (A, nm, mu), bool whether calculate the correction
           coefficients from an index ASCII external file (False) or a user defined single index (True),
           bool, whether display (True) or not (False) the plots, bool whether reduce the resolution of the
           template spectra to a user defined velscale, whether to save (true) or not (False) the plots
           in a PNG high resolution image(s).
    Output: string array containing the names of the line-strength indices considered, float array of the
            spline coefficients that best interpolate the mean EW of the broadened templates as a function
            of the velocity dispersion, float array of the
            spline coefficients that best interpolate the standard deviation from the mean EW of the broadened templates as a function
            of the velocity dispersion, float array with the mean EW of the broadened template spectra
            as a function of the velocity dispersion broadening, float array with the standard deviation
            from the mean EW of the broadened template spectra as a function of the velocity dispersion broadening,
            int binary (0-1) value used to cancel the task from the main SPAN interface in the 'Process all' mode.

    """


    #check to add the absolute path in case the spectra list is given in relative path
    with open(spectra_file, 'r') as f:
        spec_names = []
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            # If I gave the relative path, try to solve the absolute path
            if not os.path.isabs(line):
                spec_names.append(os.path.join(BASE_DIR, line))

            #If I have the absolute path
            else:
                # Loading the spectra list
                spec_names = np.loadtxt(spectra_file, dtype='str', delimiter=' ', usecols=[0])

                # Normalisation: spec_names must always be a list, even with only one element (one spectrum)
                if np.ndim(spec_names) == 0:
                    spec_names = [spec_names]  # trasnforming to a list if I have one spectrum and loadtxt sees it as a scalar
                else:
                    spec_names = spec_names.tolist()  #Converting the 1D array in a list
                break #assuming I have all absolute or relative paths, braking the cycle


    number_spec = len(spec_names)
    print ('Number of spectra', number_spec)
    #array of sigmas to be tested
    sigma_array = np.array([0, 50, 100, 150, 200, 250, 300, 350, 400])
    sigma_values = len (sigma_array)

    if is_usr_idx:
        print ('User index')
        ew_sigma = np.zeros((sigma_values, number_spec))
        ew_mean = np.zeros(sigma_values)
        ew_std = np.zeros(sigma_values)
        idx_array = 'usr_idx'
        # for every spectrum in my list:
        stop_condition = 0
        for i in range(number_spec):

            #reading the spectrum
            wavelength, flux, original_step, obj_name = stm.read_spec(spec_names[i], lambda_units)

            #preparing the plotting variables and titles
            if want_plot or save_plot:
                if i == 0:
                    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (9.5,4.5))
                    fig.suptitle('Broadening coefficients for ' + idx_array)


            wave_limits = np.array([np.min(wavelength), np.max(wavelength)])
            if (np.min(index_file) < wave_limits[0] or np.max(index_file > wave_limits[1])):
                if i == 0:
                    sg.popup ('The index definition wavelength is out of at least one spectrum. Trying other spectra...')
                else:
                    print ('The index definition wavelength is out of the spectrum')
                continue

            if not sg.OneLineProgressMeter('Task progress', i+1, number_spec,  'single', 'Processing spectra:', orientation='h',button_color=('white','red')):
                print ('***CANCELLED***')
                print ('')
                stop_condition = 1
                break

            #resample if step not constant
            step1 = wavelength[1]-wavelength[0]
            step2 = wavelength[len(wavelength)-1]- wavelength[len(wavelength)-2]
            epsilon = 1e-4
            if abs(step1-step2) > epsilon:
                wavelength, flux, npoint_resampled = spman.resample(wavelength, flux, original_step)
                print('Spectrum resampled to a linear step')

            if smooth_value != 0:
                flux = spman.sigma_broad(wavelength, flux, smooth_value)

            #for every sigma value:
            for j in range(sigma_values):

                # broadening the spectra
                flux_broadened = spman.sigma_broad(wavelength, flux, sigma_array[j])

                #measuring the index/indices, for one index
                if j == 0:
                    id_array, ew_orig, err, snr_ew_array, ew_array_mag, err_array_mag = ew_measurement(wavelength, flux_broadened, index_file, True, False, False, False, False, 'fake_name', True, result_plot_dir)
                    ew_sigma[j,i] = 0.

                #measuring the EW
                else:
                    id_array, ew, err, snr_ew_array, ew_array_mag, err_array_mag = ew_measurement(wavelength, flux_broadened, index_file, True, False, False, False, False, 'fake_name', True, result_plot_dir)

                    #storing the ew in the array containing row = sigma values; columns = spectra number
                    ew_sigma[j,i] = (ew-ew_orig)/ew_orig

            if want_plot:
                if i == 0:
                    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (9.5,4.5))
                    fig.suptitle('Broadening coefficients for ' + idx_array)

                    ax1.plot(sigma_array, ew_sigma[:,i], ls = 'none', marker = 'o', color = 'black', markersize = 1, label = 'Single values')
                else:
                    ax1.plot(sigma_array, ew_sigma[:,i], ls = 'none', marker = 'o', color = 'black', markersize = 1)

            #plotting
            if save_plot:
                if i == 0:
                    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (9.5,4.5))
                    fig.suptitle('Broadening coefficients for ' + idx_array)

                    ax1.plot(sigma_array, ew_sigma[:,i], ls = 'none', marker = 'o', color = 'black', markersize = 1, label = 'Single values')
                else:
                    ax1.plot(sigma_array, ew_sigma[:,i], ls = 'none', marker = 'o', color = 'black', markersize = 1)

        #filling the vectors with mean and std
        for h in range (sigma_values):
            ew_mean[h] = np.mean(ew_sigma[h,:])
            ew_std[h] = np.std(ew_sigma[h,:])

        ew_coeff = np.polyfit(sigma_array, ew_mean, 3)
        err_coeff = np.polyfit(sigma_array, ew_std, 3)

        #plotting
        if want_plot:
            sigma_for_fit = np.linspace(np.min(sigma_array), np.max(sigma_array), 500)
            ew_fit = np.polyval(ew_coeff, sigma_for_fit)
            ax1.errorbar(sigma_array, ew_mean, yerr = ew_std, color = 'red', ls='none', marker='o',markersize=5., label = 'Mean values')
            ax1.plot(sigma_for_fit, ew_fit, '-', color = 'red', label = 'Fit')
            ax1.set_xlabel('Broadening (km/s)')
            ax1.set_ylabel('Relative variation')
            ax1.legend(fontsize = 10)

            ax2.plot(sigma_array, ew_std, color = 'red')
            ax2.set_xlabel('Broadening (km/s)')
            ax2.set_ylabel('Error (1 sigma)')
            plt.show()
            plt.close()

        #plotting
        if save_plot:
            sigma_for_fit = np.linspace(np.min(sigma_array), np.max(sigma_array), 500)
            ew_fit = np.polyval(ew_coeff, sigma_for_fit)
            ax1.errorbar(sigma_array, ew_mean, yerr = ew_std, color = 'red', ls='none', marker='o',markersize=5., label = 'Mean values')
            ax1.plot(sigma_for_fit, ew_fit, '-', color = 'red', label = 'Fit')
            ax1.set_xlabel('Broadening (km/s)')
            ax1.set_ylabel('Relative variation')
            ax1.legend(fontsize = 10)

            ax2.plot(sigma_array, ew_std, color = 'red')
            ax2.set_xlabel('Broadening (km/s)')
            ax2.set_ylabel('Error (1 sigma)')

            plt.savefig(result_plot_dir + '/'+ 'sigma_coeff_' + idx_array +'.png', format='png', dpi=300)
            plt.close() # remember to always close the plots!

        return idx_array, ew_coeff, err_coeff, ew_mean, ew_std, stop_condition


    #if I have an index file
    if not is_usr_idx:
        idx_array, indices = read_idx(index_file)
        num_indices = len(idx_array)
        ew_sigma = np.zeros((sigma_values, number_spec))
        ew_mean = np.zeros((sigma_values, num_indices))
        ew_std = np.zeros((sigma_values, num_indices))
        num_coeff = 4
        ew_coeff_array = np.zeros((num_coeff, num_indices))
        err_coeff_array = np.zeros((num_coeff, num_indices))
        print ('Number of indices: ', num_indices)

        # for every index in my list:
        stop_condition = 0
        for k in range (num_indices):

            if not sg.OneLineProgressMeter('Task progress', k+1, num_indices,  'single', 'Processing index:', orientation='h',button_color=('white','red')):
                print ('***CANCELLED***')
                print ('')
                stop_condition = 1
                break

            # for every spectrum in my list:
            for i in range(number_spec):

                #reading the spectrum
                wavelength, flux, original_step, obj_name = stm.read_spec(spec_names[i], lambda_units)

                #preparing the plotting variables and titles
                if want_plot or save_plot:
                    if i == 0:
                        fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (9.5,4.5))
                        fig.suptitle('Broadening coefficients for ' + idx_array[k])

                #checking limits
                wave_limits_spec = np.array([np.min(wavelength), np.max(wavelength)])
                min_idx = np.min(indices[:,k])
                max_idx = np.max(indices[:, k])
                if (min_idx < wave_limits_spec[0] or max_idx > wave_limits_spec[1]):
                    print ('Warning: index not in the spectrum limits. Skipping')
                    continue

                #resample if step not constant
                step1 = wavelength[1]-wavelength[0]
                step2 = wavelength[len(wavelength)-1]- wavelength[len(wavelength)-2]
                epsilon = 1e-4
                if abs(step1-step2) > epsilon:
                    wavelength, flux, npoint_resampled = spman.resample(wavelength, flux, original_step)
                    print('Spectrum resampled to a linear step')

                if smooth_value != 0:
                    flux = spman.sigma_broad(wavelength, flux, smooth_value)

                #for every sigma value:
                for j in range(sigma_values):

                    # broadening the spectra
                    flux_broadened = spman.sigma_broad(wavelength, flux, sigma_array[j])

                    #measuring the index/indices, for one index
                    if j == 0:
                        id_array, ew_orig, err, snr_ew_array, ew_array_mag, err_array_mag = ew_measurement(wavelength, flux_broadened, indices[:,k], True, False, False, False, False, 'fake_name', True, result_plot_dir)
                        ew_sigma[j,i] = 0.
                    #measuring the EW
                    else:
                        id_array, ew, err, snr_ew_array, ew_array_mag, err_array_mag = ew_measurement(wavelength, flux_broadened, indices[:,k], True, False, False, False, False, 'fake_name', True, result_plot_dir)

                        #storing the ew in the array containing row = sigma values; columns = spectra number
                        ew_sigma[j,i] = (ew-ew_orig)/ew_orig

                # plotting the single values
                if want_plot:
                    if i == 0:
                        ax1.plot(sigma_array, ew_sigma[:,i], ls = 'none', marker = 'o', color = 'black', markersize = 1, label = 'Single values')
                    else:
                        ax1.plot(sigma_array, ew_sigma[:,i], ls = 'none', marker = 'o', color = 'black', markersize = 1)

                #Preparing the plots to be saved
                if save_plot:
                    if i == 0:
                        fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (9.5,4.5))
                        fig.suptitle('Broadening coefficients for ' + idx_array[k])

                        ax1.plot(sigma_array, ew_sigma[:,i], ls = 'none', marker = 'o', color = 'black', markersize = 1, label = 'Single values')
                    else:
                        ax1.plot(sigma_array, ew_sigma[:,i], ls = 'none', marker = 'o', color = 'black', markersize = 1)

            #filling the vectors with mean and std
            for h in range (sigma_values):
                ew_mean[h,k] = np.mean(ew_sigma[h,:])
                ew_std[h,k] = np.std(ew_sigma[h,:])

            ew_coeff = np.polyfit(sigma_array, ew_mean[:,k], 3)
            err_coeff = np.polyfit(sigma_array, ew_std[:,k], 3)

            ew_coeff_array[:,k] = ew_coeff
            err_coeff_array[:,k] = err_coeff

            #plotting
            if want_plot:
                sigma_for_fit = np.linspace(np.min(sigma_array), np.max(sigma_array), 500)
                ew_fit = np.polyval(ew_coeff, sigma_for_fit)
                ax1.errorbar(sigma_array,ew_mean[:,k], yerr = ew_std[:,k], color = 'red', ls='none', marker='o',markersize=5., label = 'Mean values')
                ax1.plot(sigma_for_fit, ew_fit, '-', color = 'red', label = 'Fit')
                ax1.set_xlabel('Broadening (km/s)')
                ax1.set_ylabel('Relative variation')
                ax1.legend(fontsize = 10)

                ax2.plot(sigma_array, ew_std[:,k], color = 'red')
                ax2.set_xlabel('Broadening (km/s)')
                ax2.set_ylabel('Relative error')
                plt.show()
                plt.close()

            #plotting
            if save_plot:
                sigma_for_fit = np.linspace(np.min(sigma_array), np.max(sigma_array), 500)
                ew_fit = np.polyval(ew_coeff, sigma_for_fit)
                ax1.errorbar(sigma_array,ew_mean[:,k], yerr = ew_std[:,k], color = 'red', ls='none', marker='o',markersize=5., label = 'Mean values')
                ax1.plot(sigma_for_fit, ew_fit, '-', color = 'red', label = 'Fit')
                ax1.set_xlabel('Broadening (km/s)')
                ax1.set_ylabel('Relative variation')
                ax1.legend(fontsize = 10)

                ax2.plot(sigma_array, ew_std[:,k], color = 'red')
                ax2.set_xlabel('Broadening (km/s)')
                ax2.set_ylabel('Relative error')

                plt.savefig(result_plot_dir + '/'+ 'sigma_coeff_' + idx_array[k] +'.png', format='png', dpi=300)
                plt.close() #Remember to always close the plots!

            plt.close()
        return idx_array, ew_coeff_array, err_coeff_array, ew_mean, ew_std, stop_condition


#*****************************************************************************************************
#9) Correction EW for the sigma broadening
def corr_ew(ew_file, corr_file, sigma_file):

    """
    This function corrects the raw EW values for
    a series of spectra for the zero velocity dispersion frame, using the
    correction coefficients estimated with the corr_ew function.
    Input: array containing the raw EWs and their uncertainties (in Angstrom) measured by
           the ew_measurement function for the n indices and m spectra,
           ASCII file containing the spectra names and the 4 spline valocity dispersion correction
           coefficients and their uncerteinties (C0, C1, C2, C3) for the n indices generated by
           the 'calculate velocity dispersion coefficients' option in the 'line-strength parametes'
           window of SPAN,
           ASCII file containing the spectra names and the velocity dispersion measured for the m spectra
           generated by the 'Velocity dispersion' task of SPAN.
    Output: string array of the index names corrected, array of the corrected EW and corrected uncertainties
            (in Angstrom).

    """

    #reading the files
    data_ew = pd.read_csv(ew_file, header=None, sep = ' ')
    data_corr = pd.read_csv(corr_file, header=None, sep = ' ')

    data_number = len(data_ew.index)

    index_number = round((len(data_ew.columns)-1)/2)
    total_column = len(data_ew.columns)
    ew_values_starting_at = 1
    ew_values_end_at = index_number+1
    err_ew_starting_at = ew_values_end_at
    err_ew_ends_at = total_column
    n_column_ew_to_correct = index_number

    column_corr_file = len(data_corr.columns)
    corr_ew_ends_at = round(column_corr_file/2)
    coeff_number = len(data_corr.index)

    ew_values = pd.DataFrame(data_ew.iloc[1:data_number, ew_values_starting_at:ew_values_end_at])
    err_values = pd.DataFrame(data_ew.iloc[1:data_number, err_ew_starting_at:err_ew_ends_at ])
    index_names = pd.DataFrame(data_ew.iloc[0:1, ew_values_starting_at:err_ew_ends_at ])
    all_idx = pd.DataFrame(data_ew.iloc[0:1, 0:err_ew_ends_at ])
    sigma_values = np.loadtxt(sigma_file, usecols = [1])

    spec_names = pd.DataFrame(data_ew.iloc[1:data_number, 0:1])

    corr_ew_values = pd.DataFrame(data_corr.iloc[1:coeff_number, 0:corr_ew_ends_at])
    corr_err_values = pd.DataFrame(data_corr.iloc[1:coeff_number, corr_ew_ends_at:column_corr_file])
    corr_index_names = pd.DataFrame(data_corr.iloc[0:1, 0:column_corr_file])

    ew_values_np = ew_values.to_numpy(dtype = float)
    err_values_np = err_values.to_numpy(dtype = float)
    sigma_values_np = sigma_values
    corr_ew_values_np = corr_ew_values.to_numpy(dtype = float)
    corr_err_values_np = corr_err_values.to_numpy(dtype = float)
    corr_index_names_np = corr_index_names.to_numpy(dtype = str)
    spec_names_np = spec_names.to_numpy()
    index_names_np = index_names.to_numpy(dtype = str)
    all_idx_np = np.loadtxt(ew_file, dtype = 'str', delimiter = ' ' , max_rows = 1, comments = '##')

    spectra_number = data_number -1
    new_ew = np.zeros((spectra_number, n_column_ew_to_correct))
    ew_correction = np.zeros((spectra_number, n_column_ew_to_correct))
    sigma_correction = np.zeros((spectra_number, n_column_ew_to_correct))
    new_err_tot = np.zeros((spectra_number, n_column_ew_to_correct))

    #wiping the zeros
    epsilon = 1e-5
    ew_values_np[ew_values_np == 0] = epsilon

    #checking correspondence between the files
    is_equal = np.char.equal(corr_index_names_np, index_names_np)
    if False in is_equal:
        print ('Ops: seems that theere is no correspondence between the ew file and the correction one!')
        print ('Doing nothing')
        return 0

    for i in range  (n_column_ew_to_correct):
        for j in range (spectra_number):

            #correcting ews
            ew_correction[j,i] = (corr_ew_values_np[3,i] + corr_ew_values_np[2,i]*sigma_values_np[j]+corr_ew_values_np[1,i]*sigma_values_np[j]**2+corr_ew_values_np[0,i]*sigma_values_np[j]**3)
            new_ew[j,i] = ew_values_np[j,i]/ (ew_correction[j,i]+1)

            #correcting uncertainties
            sigma_correction[j,i] = abs(corr_err_values_np[3,i] + corr_err_values_np[2,i]*sigma_values_np[j]+corr_err_values_np[1,i]*sigma_values_np[j]**2+corr_err_values_np[0,i]*sigma_values_np[j]**3)

            ##total uncertainties
            new_err_tot[j,i] =  mt.sqrt(sigma_correction[j,i]**2 + (err_values_np[j,i]/abs(ew_values_np[j,i]))**2)*abs(new_ew[j,i]) #c'è uno zero o qualcosa che non va nella righa 9

    #Ok, putting together
    new_ew_data = np.column_stack((spec_names_np, new_ew, new_err_tot))

    return all_idx_np, new_ew_data


#*****************************************************************************************************
#10) Correction EW for the sigma broadening of Lick/IDS indices
def corr_ew_lick(ew_values, err_values, ew_mag_values, coeff_file, sigma_value):

    """
    This function correct the raw EW values of the Lick/IDS indices from the
    broadening coefficients of Trager et al. 1998 and returns the corrected
    Lick/IDS EWs and the corrected uncertainties.
    Input: array of the raw EW of the Lick/IDS indices for the m spectra,
           array of the uncertainties of the raw Lick/IDS EW for the n indices and m spectra,
           array of the raw EW in magnitudes of the raw Lick/IDS indices for the n indices and m spectra,
           ASCII file containing the correction coefficients for the n indices stored
           in the system_files subdirectory of SPAN,
           array of the velocity dispersion measured for the m spectra and the n Lick/IDS indices.
    Output: arrays of: corrected EW, corrected uncertainties, corrected EW in magnitudes,
            corrected uncertainties in magnitudes.

    """

    #Converting the EW and EW_err to numpy
    ew_values_np = ew_values
    err_values_np = err_values

    #Open and reading the file containing the coefficients
    data_corr = pd.read_csv(coeff_file, header=None, sep = ' ')
    data_number = len(ew_values)

    corr_ew_values = data_corr.iloc[1:, :].values
    corr_ew_values_np = corr_ew_values.astype(float)

    #creating the array of the corrected ews and errors
    ew_correction = np.zeros(data_number)
    new_ew_err = np.zeros(data_number)
    new_ew = np.zeros(data_number)


    # CORRECTING THE EWS AND ERRORS, CONSIDERING BOTH IN ANGSTROM AND MAGNITUDES
    for i in range(data_number):
        #correcting ews
        ew_correction[i] = (corr_ew_values_np[3,i] + corr_ew_values_np[2,i]*sigma_value+corr_ew_values_np[1,i]*sigma_value**2+corr_ew_values_np[0,i]*sigma_value**3)
        new_ew[i] = ew_values_np[i]*ew_correction[i]

    new_ew_err = err_values_np*ew_correction
    new_ew_mag = ew_mag_values*ew_correction
    new_ew_err_mag = np.zeros(len(new_ew_err))
    for i in range(data_number):
        if new_ew[i] == 0:
            new_ew[i] = 0
        else:
            new_ew_err_mag[i]= 0.434*abs(new_ew_err[i]/new_ew[i])

    return new_ew, new_ew_err, new_ew_mag, new_ew_err_mag

#********************** END OF LNE-STRENGTH FUNCTIONS *************************************
#******************************************************************************************
