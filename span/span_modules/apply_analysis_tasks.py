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

# Functions to apply the tasks of the Spectral analysis frame

try: #try local import if executed as script
    #GUI import
    from FreeSimpleGUI_local import FreeSimpleGUI as sg

    #SPAN functions import
    from span_functions import system_span as stm
    from span_functions import utilities as uti
    from span_functions import spec_manipul as spman
    from span_functions import spec_math as spmt
    from span_functions import linestrength as ls
    from span_functions import spec_analysis as span
    from span_functions import cube_extract as cubextr
    from params import SpectraParams

except ModuleNotFoundError: #local import if executed as package
    #GUI import
    from span.FreeSimpleGUI_local import FreeSimpleGUI as sg

    #SPAN functions import
    from span.span_functions import system_span as stm
    from span.span_functions import utilities as uti
    from span.span_functions import spec_manipul as spman
    from span.span_functions import spec_math as spmt
    from span.span_functions import linestrength as ls
    from span.span_functions import spec_analysis as span
    from span.span_functions import cube_extract as cubextr
    from .params import SpectraParams

import numpy as np
import os
import glob
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.backends.backend_tkagg
from scipy.signal import argrelextrema
from dataclasses import replace


CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(CURRENT_DIR)



def apply_blackbody_fitting(event, save_plot, params):

    """
    Applies blackbody fitting to a spectrum within a specified wavelength range,
    using the blackbody_fit function of the spectral analysis module.

    Returns:
    - float: Best-fit blackbody temperature.
    - numpy array: Residuals of the blackbody fit.
    """

    wavelength, flux, wave1_bb, wave2_bb, t_guess, prev_spec_nopath, result_plot_dir, task_done, task_done2, task_analysis = params.wavelength, params.flux, params.wave1_bb, params.wave2_bb, params.t_guess, params.prev_spec_nopath, params.result_plot_dir, params.task_done, params.task_done2, params.task_analysis


    if event == 'Process all':
        task_done2 = 1
    else:
        task_done = 1
        task_analysis = 1

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_done2=task_done2, task_analysis=task_analysis)

    if wave1_bb >= wave2_bb:
        if event == "Process all":
            print("Invalid wavelength range: the first wavelength cannot be greater than the second.")
        else:
            sg.popup("Invalid wavelength range: the first wavelength cannot be greater than the second.")
        return None, None, None, None, params

    if (wave1_bb < wavelength[0] or wave2_bb > wavelength[-1]):
        if event == "Process all":
            print("Wavelength interval for blackbody fitting exceeds the spectrum's range!")
        else:
            sg.popup("Wavelength interval for blackbody fitting exceeds the spectrum's range!")
        return None, None, None, None, params

    try:
        if event == "Process all":
            preview = False
            temperature_bb, residual_bb, T_err, chi2_red = span.blackbody_fit(wavelength, flux, wave1_bb, wave2_bb, t_guess, preview, save_plot, result_plot_dir, prev_spec_nopath)
        if event == "Preview result":
            preview = True
            save_plot = False
            temperature_bb, residual_bb, T_err, chi2_red = span.blackbody_fit(wavelength, flux, wave1_bb, wave2_bb, t_guess, preview, save_plot, result_plot_dir, prev_spec_nopath)
        if event == "Process selected":
            preview = False
            save_plot = False
            temperature_bb, residual_bb, T_err, chi2_red = span.blackbody_fit(wavelength, flux, wave1_bb, wave2_bb, t_guess, preview, save_plot, result_plot_dir, prev_spec_nopath)

        print(f"Best Blackbody temperature: {int(temperature_bb)} K ± {int(T_err)}")
        print(f"Chi2: {chi2_red}")
        print('')
        return temperature_bb, residual_bb, T_err, chi2_red, params

    except Exception:
        if event == "Process all":
            print("Black-body fitting failed")
        else:
            sg.popup("Black-body fitting failed")
        return None, None, None, None, params



def apply_cross_correlation(event, save_plot, params):

    """
    Performs cross-correlation between a spectrum and a template,
    using the crosscorr function of the spectral analysis module.

    Returns:
    - Measured velocity or z value, error

    """

    # assigning params to local variables
    wavelength = params.wavelength
    flux = params.flux
    template_crosscorr = params.template_crosscorr
    lambda_units_template_crosscorr = params.lambda_units_template_crosscorr
    low_wave_corr = params.low_wave_corr
    high_wave_corr = params.high_wave_corr
    wave_interval_corr = params.wave_interval_corr
    vel_interval_corr = params.vel_interval_corr
    low_vel_corr = params.low_vel_corr
    high_vel_corr = params.high_vel_corr
    z_interval_corr = params.z_interval_corr
    is_z_xcorr = params.is_z_xcorr
    is_vel_xcorr = params.is_vel_xcorr
    low_z_corr = params.low_z_corr
    high_z_corr = params.high_z_corr
    xcorr_limit_wave_range = params.xcorr_limit_wave_range
    xcorr_vel_step = params.xcorr_vel_step
    xcorr_z_step = params.xcorr_z_step
    smooth_value_crosscorr = params.smooth_value_crosscorr
    smooth_template_crosscorr = params.smooth_template_crosscorr
    interval_corr = params.interval_corr
    prev_spec_nopath = params.prev_spec_nopath
    result_plot_dir = params.result_plot_dir
    task_done = params.task_done
    task_done2 = params.task_done2
    task_analysis = params.task_analysis

    if event == 'Process all':
        task_done2 = 1
    else:
        task_done = 1
        task_analysis = 1
        # print('*** Cross-correlation ***')

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_done2=task_done2, task_analysis=task_analysis)

    # Ensure the template file exists
    if not os.path.isfile(template_crosscorr):
        if event == "Process all":
            print("The template file does not exist. Skipping...")
        else:
            sg.popup("The template file does not exist. Skipping...")
        return None, None, params

    mode = 'z' if is_z_xcorr else 'v'
    if mode == 'z':
        grid = np.arange(low_z_corr, high_z_corr, xcorr_z_step)
    if mode == 'v':
        grid = np.arange(low_vel_corr, high_vel_corr, xcorr_vel_step)

    if mode == 'v' and (low_vel_corr < -50000 or high_vel_corr > 50000):
        print("WARNING: You have selected a velocity range beyond ±50000 km/s. Consider using redshift mode.")

    if xcorr_limit_wave_range:
        real_low_wave_corr = np.min(wave_interval_corr)
        real_high_wave_corr = np.max(wave_interval_corr)
        mask = (wavelength >= real_low_wave_corr) & (wavelength <= real_high_wave_corr)
        wave_obs = wavelength[mask]
        flux_obs = flux[mask]
    else:
        real_low_wave_corr = np.min(wavelength)
        real_high_wave_corr = np.max(wavelength)

        wave_obs = wavelength
        flux_obs = flux

    try:
        wave_temp_xcorr, flux_temp_xcorr, step_temp_xcorr, name_temp_xcorr = stm.read_spec(template_crosscorr, lambda_units_template_crosscorr)
        wave_limits_template_xcorr = np.array([wave_temp_xcorr[0], wave_temp_xcorr[len(wave_temp_xcorr)-1]])
    except Exception:
        if event == "Process all":
            print('Cannot read the template.')
        else:
            sg.popup('Cannot read the template.')
        return None, None, params

    #Fix the template name for aesthetic purposes
    template_crosscorr_nopath = os.path.splitext(os.path.basename(template_crosscorr))[0] #no path for showing and saving things

    # smooth template
    if smooth_value_crosscorr > 0 and smooth_template_crosscorr:
        flux_temp_xcorr = spman.sigma_broad(wave_temp_xcorr, flux_temp_xcorr, smooth_value_crosscorr)

    # Check wavelength limits
    if (real_low_wave_corr < np.min(wavelength) or real_high_wave_corr > np.max(wavelength) or (wave_limits_template_xcorr[1] < real_low_wave_corr or wave_limits_template_xcorr[0] > real_high_wave_corr)):
        if event == "Process all":
            print("The template does not cover the wavelength range you want to cross-correlate!")
        else:
            sg.popup("The template does not cover the wavelength range you want to cross-correlate!")
        return None, None, params

    #Cross-correlate
    try:
        best_shift, cc_values, grid_used, sigma = span.estimate_from_template(
            wave_obs, flux_obs, wave_temp_xcorr, flux_temp_xcorr, grid, mode=mode)

        unit = 'z' if mode == 'z' else 'km/s'
        if mode == 'v':
            print(f"Best {unit}: {best_shift:.2f} ± {sigma:.2f}")
            print('')
        if mode == 'z':
            print(f"Best {unit}: {best_shift:.5f} ± {sigma:.5f}")
            print('')


        if event == 'Preview result':
            save_crosscorr_plot = False
            span.plot_cross_correlation(grid_used, cc_values, best_shift, wave_obs, flux_obs, wave_temp_xcorr, flux_temp_xcorr, best_shift, save_crosscorr_plot, prev_spec_nopath, result_plot_dir, mode=mode)
        if event == 'Process all' and save_plot:
            span.plot_cross_correlation(grid_used, cc_values, best_shift, wave_obs, flux_obs, wave_temp_xcorr, flux_temp_xcorr, best_shift, save_plot, prev_spec_nopath, result_plot_dir, mode=mode)

        return best_shift, sigma, params

    except Exception:
        if event == "Process all":
            print('Cannot find cross-correlation within the ranges you inserted. Try again with different (smaller?) ranges')
        else:
            sg.Popup ('Cannot find cross-correlation within the ranges you inserted. Try again with different (smaller?) ranges')

        return None, None, params



def apply_velocity_dispersion(event, save_plot, params):

    """
    Performs a simple fit of the spectrum to retrieve the velocity dispersion
    between a spectrum and a template, using the sigma_measurement function of the spectral analysis module.

    Returns:
    - Velocity dispersion
    - Uncertainty of the velocity dispersion
    - Chi2 of the fit

    """

    # params to local variables
    wavelength = params.wavelength
    flux = params.flux
    template_sigma = params.template_sigma
    lambda_units_template_sigma = params.lambda_units_template_sigma
    resolution_spec = params.resolution_spec
    resolution_template = params.resolution_template
    band_sigma = params.band_sigma
    cont_sigma = params.cont_sigma
    band_custom = params.band_custom
    prev_spec_nopath = params.prev_spec_nopath
    task_done = params.task_done
    task_done2 = params.task_done2
    task_analysis = params.task_analysis
    result_plot_dir = params.result_plot_dir


    # 1) HEADER
    if event == 'Process all':
        task_done2 = 1
    else:
        task_done = 1
        task_analysis = 1

        #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_done2=task_done2, task_analysis=task_analysis)

    # 2) FILES AND PARAMETER CHECK
    try:
        # Check if template file exists
        if not os.path.isfile(template_sigma):
            if event == "Process all":
                print("The template file does not exist. Skipping...")
            else:
                sg.popup("The template file does not exist. Skipping...")
            return None, None, None, params

        # Read template spectrum
        wave_temp_sigma, flux_temp_sigma, _, _ = stm.read_spec(template_sigma, lambda_units_template_sigma)
        wave_limits_template_sigma = np.array([wave_temp_sigma[0], wave_temp_sigma[-1]])

        # Check if the band limits are within both the spectrum and template
        if ((band_sigma[0] < wavelength[0] or cont_sigma[0] < wavelength[0] or
            band_sigma[1] > wavelength[-1] or cont_sigma[1] > wavelength[-1] or
            band_sigma[0] < wave_limits_template_sigma[0] or cont_sigma[0] < wave_limits_template_sigma[0] or
            band_sigma[1] > wave_limits_template_sigma[-1] or cont_sigma[1] > wave_limits_template_sigma[-1])):
            if event == "Process all":
                print("Wavelength interval for the band is out of the range of the spectrum or the template!")
            else:
                sg.popup("Wavelength interval for the band is out of the range of the spectrum or the template!")
            return None, None, None, params

        # Ensure band_sigma and cont_sigma are properly ordered
        if band_custom and (band_sigma[0] > band_sigma[1] or cont_sigma[0] > cont_sigma[1]):
            if event == "Process all":
                print("It seems we have a problem. Did you invert the wavelength range?")
            else:
                sg.popup("It seems we have a problem. Did you invert the wavelength range?")
            return None, None, None, params


        # 3) CALCULATION
        # Process based on the event type
        if event == "Process selected" or event == "Process all":
            compute_errors = True
        else:
            compute_errors = False

        sigma, error, chisqr, band_wave, band_flux, band_flux_template_fitted, sigma_instrumental = span.sigma_measurement(
            wavelength, flux, template_sigma, lambda_units_template_sigma, resolution_spec, resolution_template,
            band_sigma, cont_sigma, compute_errors
        )

        print(prev_spec_nopath)
        print(f"Sigma: {sigma} +/- {round(error,2)} km/s    Chi-Square: {chisqr}\n")
        print('')

        # 4) PLOTTING
        # If previewing, show the plot
        if event == "Preview result" or (event == 'Process all' and save_plot):
            if sigma == sigma_instrumental:
                sg.popup("WARNING: The real velocity dispersion is lower than the instrumental sigma of the spectrum. Do not trust the result!")

            # Plot results
            sigma_title = str(round(sigma, 1))

            fig, (ax1, ax2) = plt.subplots(2, figsize=(8.5, 7), gridspec_kw={"height_ratios": [3, 1]})
            fig.suptitle(f"{prev_spec_nopath}  Sigma: {sigma_title} km/s")

            # Spectrum and fitted template
            ax1.plot(band_wave, band_flux, label="Spectrum")
            ax1.plot(band_wave, band_flux_template_fitted, label="Fitted template")
            ax1.set_xlabel("Wavelength (A)")
            ax1.set_ylabel("Norm flux")
            ax1.legend(fontsize=10)

            # Residuals
            ax2.plot(band_wave, band_flux - band_flux_template_fitted, linewidth=0.5, label="Residuals")
            ax2.hlines(y=0, xmin=min(band_wave), xmax=max(band_wave), linestyles="--", lw=2, color="r")
            ax2.set_xlabel("Wavelength (A)")
            ax2.set_ylabel("Residuals")
            ax2.legend(fontsize=10)

            if event == "Preview result":
                plt.show()
                plt.close()
            else:
                plt.savefig(result_plot_dir + '/'+ 'sigma_vel_' + prev_spec_nopath + '.png', format='png', dpi=300)
                plt.close()


        return sigma, error, chisqr, params

    except Exception as e:
        if event == "Process all":
            print(f"Velocity dispersion measurement failed: {str(e)}")
        else:
            sg.popup(f"Velocity dispersion measurement failed: {str(e)}")
        return None, None, None, params  # Ensures consistent return type



def apply_ew_measurement_single(event, save_plot, params):

    """
    Measures the Equivalent Width (EW) for a single spectral index,
    using the ew_measurement function of the linestrength module.

    Returns:
    - ew_array: float, Measured equivalent width
    - err_array: float, Measurement error
    - ew_array_mag: float, EW in magnitude scale
    - err_array_mag: float, EW magnitude error
    - snr_ew_array: float, Signal-to-noise ratio of EW

    """

    # params to local variables
    wavelength = params.wavelength
    flux = params.flux
    index_usr = params.index_usr
    wave_limits = params.wave_limits
    prev_spec =  params.prev_spec
    prev_spec_nopath = params.prev_spec_nopath
    result_plot_dir = params.result_plot_dir
    task_done = params.task_done
    task_done2 = params.task_done2
    task_analysis = params.task_analysis


    # 1) HEADER
    if event == 'Process all':
        task_done2 = 1
    else:
        task_done = 1
        task_analysis = 1
        # print('*** Equivalent width measurement for a single index ***')

        #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_done2=task_done2, task_analysis=task_analysis)

    # Check if the index definition is within spectrum limits
    if np.min(index_usr) < wave_limits[0] or np.max(index_usr) > wave_limits[1]:
        if event == "Process all":
            print("The index definition wavelength exceeds the wavelength range of the spectrum")
        else:
            sg.popup("The index definition wavelength exceeds the wavelength range of the spectrum")
        return None, None, None, None, None, None, params

    # Check if the index definition is in the correct order
    if index_usr[0] > index_usr[1] or index_usr[2] > index_usr[3] or index_usr[4] > index_usr[5]:
        if event == "Process all":
            print("It seems we have a problem. Did you invert the wavelengths of the indices?")
        else:
            sg.popup("It seems we have a problem. Did you invert the wavelengths of the indices?")
        return None, None, None, None, None, None, params

    if event == "Preview result":
        # Parameters for EW measurement
        plot = True
        verbose = True
        with_uncertainties = True
        save_plot = False
        normalise_spec = True
    elif event == 'Process all' and save_plot:
        plot = False
        verbose = True
        with_uncertainties = True
        save_plot = True
        normalise_spec = True
    else:
        plot = False
        verbose = True
        with_uncertainties = True
        save_plot = False
        normalise_spec = True

    # Perform EW measurement
    try:
        idx, ew, err, snr_ew, ew_mag, err_mag = ls.ew_measurement(
            wavelength, flux, index_usr, True, plot, verbose, with_uncertainties,
            save_plot, prev_spec, normalise_spec, result_plot_dir
        )

        print("EW:", round(ew, 3), "+/-", round(err, 3))
        print("EW Mag:", round(ew_mag, 3), "+/-", round(err_mag, 3))
        print("SNR:", int(snr_ew), "per pix")
        print("")

        return idx, ew, err, snr_ew, ew_mag, err_mag, params

    except Exception as e:
        if event == "Process all":
            print(f"EW measurement failed. Try adjusting the parameters: {e}")
        else:
            sg.popup(f"EW measurement failed. Try adjusting the parameters: {e}")
        return None, None, None, None, None, None, params



def apply_ew_measurement_list(event, save_plot, params):

    """
    Measures the Equivalent Width (EW) for a list of spectral indices,
    using the ew_measurement function of the linestrength module.

    Returns:
    - id_array: list, Index IDs
    - ew_array: array, Measured equivalent width values
    - err_array: array, Measurement errors
    - snr_ew_array: array, Signal-to-noise ratios
    - ew_array_mag: array, EW in magnitude scale
    - err_array_mag: array, EW magnitude errors

    """

    wavelength = params.wavelength
    flux = params.flux
    index_file = params.index_file
    single_index = params.single_index
    prev_spec = params.prev_spec
    prev_spec_nopath = params.prev_spec_nopath
    result_plot_dir = params.result_plot_dir
    task_done = params.task_done
    task_done2 = params.task_done2
    task_analysis = params.task_analysis


    # 1) HEADER
    if event == 'Process all':
        task_done2 = 1
    else:
        task_done = 1
        task_analysis = 1

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_done2=task_done2, task_analysis=task_analysis)

    # Check if the index file exists
    if not os.path.isfile(index_file):
        if event == "Process all":
            print("The index file does not exist. Skipping...")
        else:
            sg.popup("The index file does not exist. Skipping...")
        return None, None, None, None, None, None, params

    # Try reading the index file
    try:
        idx_names, indices = ls.read_idx(index_file)
    except ValueError:
        if event == "Process all":
            print("At least one index in the file is not valid")
        else:
            sg.popup("At least one index in the file is not valid")
        return None, None, None, None, None, None, params

    # Validate index length
    if len(indices[:, 0]) < 6:
        if event == "Process all":
            print("The length of at least one index is not correct")
        else:
            sg.popup("The length of at least one index is not correct")
        return None, None, None, None, None, None, params

    # Check for wavelength order issues
    bad_idx = []
    for t in range(len(idx_names)):
        if indices[0, t] > indices[1, t] or indices[2, t] > indices[3, t] or indices[4, t] > indices[5, t]:
            bad_idx.append(idx_names[t])

    if bad_idx:
        if event == "Process all":
            print("It seems we have a problem. Did you invert the wavelengths of these indices?", bad_idx)
        else:
            sg.popup("It seems we have a problem. Did you invert the wavelengths of these indices?", bad_idx)
        return None, None, None, None, None, None, params

    if event == "Preview result":
        # Parameters for EW measurement
        plot = True
        verbose = True
        with_uncertainties = True
        save_plot = False
        normalise_spec = True
    elif event == 'Process all' and save_plot:
        plot = False
        verbose = True
        with_uncertainties = True
        save_plot = True
        normalise_spec = True
    else:
        plot = False
        verbose = True
        with_uncertainties = True
        save_plot = False
        normalise_spec = True

    # Perform EW measurement
    try:
        id_array, ew_array, err_array, snr_ew_array, ew_array_mag, err_array_mag = ls.ew_measurement(wavelength, flux, index_file, single_index, plot, verbose, with_uncertainties, save_plot, prev_spec, normalise_spec, result_plot_dir)

        print(id_array)
        print("EW in A:", np.round(ew_array, 3))
        print("EW in Mag:", np.round(ew_array_mag, 3))
        print("SNR:", np.round(snr_ew_array), "per pix")
        print('')
        return id_array, ew_array, err_array, snr_ew_array, ew_array_mag, err_array_mag, params

    except Exception as e:
        if event == "Process all":
            print(f"EW measurement failed. Try adjusting the parameters: {e}")
        else:
            sg.popup(f"EW measurement failed. Try adjusting the parameters: {e}")
        return None, None, None, None, None, None, params



def apply_lick_indices_ew_measurement(event, save_plot, i, params):

    """
    Measures the Equivalent Width (EW) of Lick/IDS indices,
    with optional emission correction, Doppler correction,
    and velocity dispersion corrections, using the ew_measurement
    function of the linestrength module amd the ppxf_pop function
    of the spectral analysis module for gas correction, doppler
    correction and auto sigma correction options.

    Returns:
    - tuple: (lick_id_array, lick_ew_array, lick_err_array, lick_snr_ew_array, lick_ew_array_mag, lick_err_array_mag)

    """


    wavelength = params.wavelength
    flux = params.flux
    lick_index_file = params.lick_index_file
    lick_correct_emission = params.lick_correct_emission
    dop_correction_lick = params.dop_correction_lick
    correct_ew_sigma = params.correct_ew_sigma
    radio_lick_sigma_auto = params.radio_lick_sigma_auto
    radio_lick_sigma_single = params.radio_lick_sigma_single
    radio_lick_sigma_list = params.radio_lick_sigma_list
    sigma_lick_coeff_file = params.sigma_lick_coeff_file
    lick_constant_fwhm = params.lick_constant_fwhm
    spec_lick_res_fwhm = params.spec_lick_res_fwhm
    spec_lick_res_r = params.spec_lick_res_r
    z_guess_lick_emission = params.z_guess_lick_emission
    stellar_parameters_lick = params.stellar_parameters_lick
    ssp_model = params.ssp_model
    interp_model = params.interp_model
    sigma_lick_file = params.sigma_lick_file
    spectra_number_to_process = params.spectra_number_to_process
    prev_spec = params.prev_spec
    prev_spec_nopath = params.prev_spec_nopath
    result_plot_dir = params.result_plot_dir
    task_done = params.task_done
    task_done2 = params.task_done2
    task_analysis = params.task_analysis


    # 1) HEADER
    if event == 'Process all':
        task_done2 = 1
    else:
        task_done = 1
        task_analysis = 1
        # print('*** Equivalent Width Measurement of Lick/IDS Indices ***')

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_done2=task_done2, task_analysis=task_analysis)

    # 2) FILES AND PARAMETER CHECK
    age = 0
    met = 0
    alpha = 0
    err_age = 0
    err_met = 0
    err_alpha = 0
    lick_for_ssp = np.zeros(14)
    ssp_lick_indices_list = None
    ssp_lick_indices_err_list = None

    try:
        # 1) Read Lick index file
        lick_idx_names, lick_indices = ls.read_idx(lick_index_file)
        lick_wavelength = wavelength
        lick_flux = flux

        # Spectrum wavelength range checks
        lick_wave_limits = np.array([lick_wavelength[0], lick_wavelength[-1]])
        lick_wave_lower_limit, lick_wave_upper_limit = 4200, 6300

        if (lick_wave_limits[0] < lick_wave_lower_limit and lick_wave_limits[1] < lick_wave_lower_limit) or \
           (lick_wave_limits[0] > lick_wave_upper_limit and lick_wave_limits[1] > lick_wave_upper_limit):
            if event == "Process all":
                print("The window band is out of the spectrum range")
            else:
                sg.popup("The window band is out of the spectrum range")
            return None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, params


    # 2) If I want to correct for the emission or doppler correction or sigma correction
        if lick_correct_emission or dop_correction_lick or (correct_ew_sigma and radio_lick_sigma_auto):

            additive_degree_lick = -1 #Not using additive degree to preserve the absorption lines
            multiplicative_degree_lick = 10
            min_wavelength = np.min(lick_wavelength)
            max_wavelength = np.max(lick_wavelength)

            #since ppxf is time consuming, let's limit the wavelength interval to correct the emission
            if max_wavelength > 6500:
                max_wavelength_lick = 6500
            else:
                max_wavelength_lick = max_wavelength
            if min_wavelength < 4000:
                min_wavelength_lick = 4000
            else:
                min_wavelength_lick = min_wavelength

            sigma_guess_lick = 100
            fit_components_lick = ('with_gas')
            lick_ppxf_with_errors = False
            lick_save_plot = False
            regul_err_lick = 0.04 #force some regularization to smooth the fit
            ppxf_pop_noise_lick = 0.0163 #use a random noise value
            ppxf_pop_age_range_lick = np.array([0,16])
            ppxf_pop_met_range_lick = np.array([-2,0.8])

            if lick_correct_emission:
                print ('Removing the emission lines with ppxf...')
            if dop_correction_lick:
                print('Doppler correction with ppxf...')

            emiles_resolution = 2.5

            #setting up the other parameters for pPXF
            stellar_library_lick = 'emiles'
            convolve_temp_lick = True #convolve the templates to galaxy resolution to a better estimation of sigma
            have_user_mask_lick = False
            mask_ranges_lick = None
            custom_emiles_lick = False
            custom_emiles_folder_lick = None
            custom_npz_lick = False
            custom_npz_file_lick = None
            mask_emission_lick = False
            custom_temp_suffix_lick = None
            best_param_lick = False
            best_noise_estimate_lick = True
            frac_chi_lick = 1
            dust_correction_stars_lick = False
            dust_correction_gas_lick = False
            tied_balmer_lick = False
            spec_name_lick = None
            with_plots_lick = False
            ppxf_pop_error_nsim_lick = None
            lick_lg_age = True
            lick_lg_met = True

            #performing ppxf
            if lick_constant_fwhm: # with a fixed delta lambda resolution, no problem
                #issue a warning in case the resolution of EMILES templates is lower than the galaxy
                if spec_lick_res_fwhm < emiles_resolution:
                    print('WARNING: the resolution of the templates is REALLY lower than the galaxy. Consider to first reduce the resolution of your spectrum with the Degrade resolution task')

                kinematics_lick, info_pop_lick, info_pop_mass_lick, mass_light_lick, errors_lick, galaxy_lick, bestfit_flux_lick, bestfit_wave_lick, bestfit_flux_gas_lick, residual_lick, chi_square_lick, age_err_abs_lick, met_err_lick, alpha_err_lick, mass_age_err_abs_lick, mass_met_err_lick, mass_alpha_err_lick, emission_corrected_flux, pop_age, light_weights_age_bin, mass_weights_age_bin, cumulative_mass, light_weights_age_std, mass_weights_age_std, cumulative_light_std, cumulative_mass_std, snr_pop_lick, light_weights_lick, mass_weights_lick, t50_age_lick, t80_age_lick, t50_cosmic_lick, t80_cosmic_lick = span.ppxf_pop(lick_wavelength, lick_flux, min_wavelength_lick, max_wavelength_lick, spec_lick_res_fwhm, z_guess_lick_emission, sigma_guess_lick, fit_components_lick, with_plots_lick, lick_ppxf_with_errors, lick_save_plot, spec_name_lick, regul_err_lick, additive_degree_lick, multiplicative_degree_lick, tied_balmer_lick, stellar_library_lick, dust_correction_stars_lick, dust_correction_gas_lick, ppxf_pop_noise_lick, ppxf_pop_age_range_lick, ppxf_pop_met_range_lick, custom_emiles_lick, custom_emiles_folder_lick, custom_npz_lick, custom_npz_file_lick, mask_emission_lick, custom_temp_suffix_lick, best_param_lick, best_noise_estimate_lick, frac_chi_lick, convolve_temp_lick, have_user_mask_lick, mask_ranges_lick, ppxf_pop_error_nsim_lick, lick_lg_age, lick_lg_met, result_plot_dir)


            #if I don't have a constant FWHM resolution:
            if not lick_constant_fwhm:
                #considering a mean fwhm along the lick band. The errors for the indices used for stellar parameters are about 5%, so it's good for now.
                mean_ref_lick_wavelength = 5080
                spec_lick_res_fwhm = mean_ref_lick_wavelength/spec_lick_res_r

                #issue a warning in case the resolution of EMILES templates is lower than the galaxy
                if spec_lick_res_fwhm < emiles_resolution:
                    print('WARNING: the resolution of the templates is REALLY lower than the galaxy. Consider to first reduce the resolution of your spectrum with the Degrade resolution task')

                kinematics_lick, info_pop_lick, info_pop_mass_lick, mass_light_lick, errors_lick, galaxy_lick, bestfit_flux_lick, bestfit_wave_lick, bestfit_flux_gas_lick, residual_lick, chi_square_lick, age_err_abs_lick, met_err_lick, alpha_err_lick, mass_age_err_abs_lick, mass_met_err_lick, mass_alpha_err_lick, emission_corrected_flux, pop_age, light_weights_age_bin, mass_weights_age_bin, cumulative_mass, light_weights_age_std, mass_weights_age_std, cumulative_light_std, cumulative_mass_std, snr_pop_lick, light_weights_lick, mass_weights_lick, t50_age_lick, t80_age_lick, t50_cosmic_lick, t80_cosmic_lick = span.ppxf_pop(lick_wavelength, lick_flux, min_wavelength_lick, max_wavelength_lick, spec_lick_res_fwhm, z_guess_lick_emission, sigma_guess_lick, fit_components_lick, with_plots_lick, lick_ppxf_with_errors, lick_save_plot, spec_name_lick, regul_err_lick, additive_degree_lick, multiplicative_degree_lick, tied_balmer_lick, stellar_library_lick, dust_correction_stars_lick, dust_correction_gas_lick, ppxf_pop_noise_lick, ppxf_pop_age_range_lick, ppxf_pop_met_range_lick, custom_emiles_lick, custom_emiles_folder_lick, custom_npz_lick, custom_npz_file_lick, mask_emission_lick, custom_temp_suffix_lick, best_param_lick, best_noise_estimate_lick, frac_chi_lick, convolve_temp_lick, have_user_mask_lick, mask_ranges_lick, ppxf_pop_error_nsim_lick, lick_lg_age, lick_lg_met, result_plot_dir)


            if lick_correct_emission:
                lick_wavelength = bestfit_wave_lick
                lick_flux = emission_corrected_flux #using the emission corrected flux from PPXF
                lick_step = lick_wavelength[1] - lick_wavelength[0]

                #rebinning linear
                lick_wavelength, lick_flux, npoint_resampled = spman.resample(lick_wavelength, lick_flux, lick_step)


            if dop_correction_lick: #CAUTION HERE: IF PPXF FINDS NO GAS COMPONENTS?
                #rebinning linear if not done before
                if not lick_correct_emission:
                    #rebinning linear
                    lick_flux = galaxy_lick #using the galaxy flux from pPXF
                    lick_wavelength = bestfit_wave_lick #using the wavelength grid from pPXF
                    lick_step = lick_wavelength[1] - lick_wavelength[0]
                    lick_wavelength, lick_flux, npoint_resampled = spman.resample(lick_wavelength, lick_flux, lick_step) #Rebinning linear
                lick_doppler_vel = (kinematics_lick[0])
                dop_vel = lick_doppler_vel[0]
                lick_wavelength, lick_flux = spman.dopcor(lick_wavelength, lick_flux, dop_vel, True) #doppler correction. The cosmological z correction has been already performed by the span.ppxf_pop function, prior to really run pPXF. Here I correct only for the real velocity component measured by the fit.

            if radio_lick_sigma_auto:
                sigma_lick_ppxf = (kinematics_lick[0])
                sigma_to_correct_lick = sigma_lick_ppxf[1]

        # 3) degrading the resolution, only if smaller than the lick system
        if lick_constant_fwhm and spec_lick_res_fwhm < 8.4:
            lick_degraded_wavelength, lick_degraded_flux = spman.degrade_to_lick(lick_wavelength, lick_flux, spec_lick_res_fwhm, lick_constant_fwhm)
        elif not lick_constant_fwhm and spec_lick_res_r > 6000:
            lick_degraded_wavelength, lick_degraded_flux = spman.degrade_to_lick(lick_wavelength, lick_flux, spec_lick_res_r, lick_constant_fwhm)
        else:
            print('WARNING: The resolution of the spectrum is smaller than the one needed for the Lick/IDS system. I will still calculate the Lick/IDS indices but the results might be inaccurate')
            lick_degraded_wavelength = lick_wavelength
            lick_degraded_flux = lick_flux

        # 4) Measuring the EW and doing plot
        if (event == 'Preview result'):
            lick_single_index = False
            lick_ew_plot = True
            lick_verbose = True
            lick_with_uncertainties = True
            lick_save_plot = False
            lick_normalise_spec = True

            lick_id_array, lick_ew_array, lick_err_array, lick_snr_ew_array, lick_ew_array_mag, lick_err_array_mag = ls.ew_measurement(lick_degraded_wavelength, lick_degraded_flux, lick_index_file, lick_single_index, lick_ew_plot, lick_verbose, lick_with_uncertainties, lick_save_plot, prev_spec, lick_normalise_spec, result_plot_dir)

            print (lick_id_array)
            print (lick_ew_array)
            print ('')
            print ('Raw EW in Mag: ')
            print (np.round(lick_ew_array_mag, decimals = 3))
            print ('SNR: ')
            print (np.round(lick_snr_ew_array), 'per pix')
            print ('')

            # 5) Correcting the EWs
            if correct_ew_sigma and radio_lick_sigma_single:
                corrected_lick_ew_array, corrected_lick_err_array, corrected_lick_ew_mag_array, corrected_lick_err_mag_array = ls.corr_ew_lick(lick_ew_array, lick_err_array, lick_ew_array_mag, sigma_lick_coeff_file, sigma_single_lick)

                print ('Corrected EWs for sigma (A):')
                print (np.round(corrected_lick_ew_array, decimals = 3))
                print ('Corrected Errors for sigma:')
                print (np.round(corrected_lick_err_array, decimals = 3))
                print ('Corrected EWs for sigma (mag):')
                print (np.round(corrected_lick_ew_mag_array, decimals = 3))
                print ('Corrected mag Errors for sigma: ')
                print (np.round(corrected_lick_err_mag_array, decimals = 3))
                print ('')

                #uodating the values
                lick_ew_array = corrected_lick_ew_array
                lick_err_array = corrected_lick_err_array
                lick_ew_array_mag = corrected_lick_ew_mag_array
                lick_err_array_mag = corrected_lick_err_mag_array

            if correct_ew_sigma and radio_lick_sigma_list:
                print ('WARNING: cannot correct for sigma broadening one spectrum with a list of sigmas. Select the single value. Skipping the correction for now...')

            if correct_ew_sigma and radio_lick_sigma_auto:
                corrected_lick_ew_array, corrected_lick_err_array, corrected_lick_ew_mag_array, corrected_lick_err_mag_array = ls.corr_ew_lick(lick_ew_array, lick_err_array, lick_ew_array_mag, sigma_lick_coeff_file, sigma_to_correct_lick)

                print ('Corrected EWs for sigma (A):')
                print (np.round(corrected_lick_ew_array, decimals = 3))
                print ('Corrected Errors for sigma:')
                print (np.round(corrected_lick_err_array, decimals = 3))
                print ('Corrected EWs for sigma (mag):')
                print (np.round(corrected_lick_ew_mag_array, decimals = 3))
                print ('Corrected mag Errors for sigma: ')
                print (np.round(corrected_lick_err_mag_array, decimals = 3))
                print ('')

                #uodating the values
                lick_ew_array = corrected_lick_ew_array
                lick_err_array = corrected_lick_err_array
                lick_ew_array_mag = corrected_lick_ew_mag_array
                lick_err_array_mag = corrected_lick_err_mag_array


        else: #In the case of process selected event, no need to show the plots
            lick_single_index = False
            lick_ew_plot = False
            lick_verbose = True
            lick_with_uncertainties = True
            lick_normalise_spec = True

            lick_id_array, lick_ew_array, lick_err_array, lick_snr_ew_array, lick_ew_array_mag, lick_err_array_mag = ls.ew_measurement(lick_degraded_wavelength, lick_degraded_flux, lick_index_file, lick_single_index, lick_ew_plot, lick_verbose, lick_with_uncertainties, save_plot, prev_spec, lick_normalise_spec, result_plot_dir)

            print (lick_id_array)
            print (lick_ew_array)
            print ('')
            print ('Raw EW in Mag: ')
            print (np.round(lick_ew_array_mag, decimals = 3))
            print ('SNR: ')
            print (np.round(lick_snr_ew_array), 'per pix')
            print ('')

            # 5) Correcting the EWs
            if correct_ew_sigma and radio_lick_sigma_single:
                corrected_lick_ew_array, corrected_lick_err_array, corrected_lick_ew_mag_array, corrected_lick_err_mag_array = ls.corr_ew_lick(lick_ew_array, lick_err_array, lick_ew_array_mag, sigma_lick_coeff_file, sigma_single_lick)

                print ('Corrected EWs for sigma (A):')
                print (np.round(corrected_lick_ew_array, decimals = 3))
                print ('Corrected uncertainties for sigma:')
                print (np.round(corrected_lick_err_array, decimals = 3))
                print ('Corrected EWs for sigma (mag):')
                print (np.round(corrected_lick_ew_mag_array, decimals = 3))
                print ('Corrected mag uncertainties for sigma: ')
                print (np.round(corrected_lick_err_mag_array, decimals = 3))
                print ('')

                #uodating the values
                lick_ew_array = corrected_lick_ew_array
                lick_err_array = corrected_lick_err_array
                lick_ew_array_mag = corrected_lick_ew_mag_array
                lick_err_array_mag = corrected_lick_err_mag_array

            if correct_ew_sigma and radio_lick_sigma_list and event == "Process selected":
                print ('WARNING: cannot correct for sigma broadening one spectrum with a list of sigmas. Select the single value. Skipping the correction for now...')


            # b) Correcting the EWs for a list
            if correct_ew_sigma and radio_lick_sigma_list and event == "Process all":

                # reading the sigma value file
                sigma_values = np.loadtxt(sigma_lick_file, usecols = [1]) #for now it's ok

                #check if the length is the same of the spectra_number to correct
                if len(sigma_values) != spectra_number:
                    sg.popup ('The sigma list file for Lick correction has a length different from the number of spectra you want to correct! I will continue without sigma correction...')

                if len(sigma_values) == spectra_number:
                    corrected_lick_ew_array, corrected_lick_err_array, corrected_lick_ew_mag_array, corrected_lick_err_mag_array = ls.corr_ew_lick(lick_ew_array, lick_err_array, lick_ew_array_mag, sigma_lick_coeff_file, sigma_values[i])

                    #updating the values:
                    lick_ew_array = corrected_lick_ew_array
                    lick_err_array = corrected_lick_err_array
                    lick_ew_array_mag = corrected_lick_ew_mag_array
                    lick_err_array_mag = corrected_lick_err_mag_array


            if correct_ew_sigma and radio_lick_sigma_auto:

                corrected_lick_ew_array, corrected_lick_err_array, corrected_lick_ew_mag_array, corrected_lick_err_mag_array = ls.corr_ew_lick(lick_ew_array, lick_err_array, lick_ew_array_mag, sigma_lick_coeff_file, sigma_to_correct_lick)

                print ('Corrected EWs for sigma (A):')
                print (np.round(corrected_lick_ew_array, decimals = 3))
                print ('Corrected uncertainties for sigma:')
                print (np.round(corrected_lick_err_array, decimals = 3))
                print ('Corrected EWs for sigma (mag):')
                print (np.round(corrected_lick_ew_mag_array, decimals = 3))
                print ('Corrected mag uncertainties for sigma: ')
                print (np.round(corrected_lick_err_mag_array, decimals = 3))
                print('')
                #uodating the values
                lick_ew_array = corrected_lick_ew_array
                lick_err_array = corrected_lick_err_array
                lick_ew_array_mag = corrected_lick_ew_mag_array
                lick_err_array_mag = corrected_lick_err_mag_array


        #6) Constraining the stellar parameters and uncertainties
        if stellar_parameters_lick:

            #assigning meaningful names to the indices used for stellar populations and creating the combined ones
            Hbeta = lick_ew_array[0]
            Hbetae = lick_err_array[0]
            Mg2 = lick_ew_array_mag[1]
            Mg2e = lick_err_array_mag[1]
            Mgb = lick_ew_array[2]
            Mgbe = lick_err_array[2]
            Fe5270 = lick_ew_array[3]
            Fe5270e = lick_err_array[3]
            Fe5335 = lick_ew_array[4]
            Fe5335e = lick_err_array[4]
            Fem = (Fe5270+Fe5335)/2
            Feme = np.sqrt((0.5*Fe5270e)**2+(0.5*Fe5335e)**2)
            MgFe = (np.sqrt(Mgb*(0.72*Fe5270+0.28*Fe5335)))
            MgFe = np.nan_to_num(MgFe, nan=0)
            MgFee = np.sqrt((((Fe5270*18/25+Fe5335*7/25)/(2*np.sqrt(Mgb*(Fe5270*18/25+Fe5335*7/25))))*Mgbe)**2+((Mgb*18/25/(2*np.sqrt(Mgb*(Fe5270*18/25+Fe5335*7/25))))*Fe5270e)**2+((Mgb*7/25/(2*np.sqrt(Mgb*(Fe5270*18/25+Fe5335*7/25))))*Fe5335e)**2)
            MgFee = np.nan_to_num(MgFee, nan=0)

            ssp_lick_indices_list = np.column_stack((Hbeta, MgFe, Fem, Mgb))
            ssp_lick_indices = ssp_lick_indices_list.reshape(-1)
            ssp_lick_indices_err_list = np.column_stack((Hbetae, MgFee, Feme, Mgbe))
            ssp_lick_indices_err = ssp_lick_indices_err_list.reshape(-1)
            lick_for_ssp = np.column_stack((Hbeta, Hbetae, Mg2, Mg2e, Mgb, Mgbe, Fe5270, Fe5270e, Fe5335, Fe5335e, Fem, Feme, MgFe, MgFee))
            lick_for_ssp = lick_for_ssp.reshape(-1,) #collapsing to a 1D array

            age, met, alpha, err_age, err_met, err_alpha = span.lick_pop(ssp_lick_indices, ssp_lick_indices_err, ssp_model, interp_model)

            #doing the plot only for Preview
            if event == 'Preview result':
                span.lick_grids(ssp_model, ssp_lick_indices_list, ssp_lick_indices_err_list, age, True, False, 'none', result_plot_dir)


            print ('Age (Gyr):', round(age,2), '+/-', round(err_age, 2))
            print ('[M/H] (dex):', round(met, 2), '+/-', round(err_met,2))
            print ('[Alpha/Fe]:', round(alpha, 2), '+/-', round(err_alpha,2))
            print('')


        return lick_id_array, lick_ew_array, lick_err_array, lick_snr_ew_array, lick_ew_array_mag, lick_err_array_mag, age, met, alpha, err_age, err_met, err_alpha, lick_for_ssp, ssp_model, ssp_lick_indices_list, ssp_lick_indices_err_list, params

    except Exception as e:
        if event == "Process all":
            print(f"EW measurement failed. Try adjusting the parameters: {e}")
        else:
            sg.popup(f"EW measurement failed. Try adjusting the parameters: {e}")
        return None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, params



def apply_cat_line_fitting(event, save_plot, params):

    """
    Fits the Calcium Triplet lines (CaT) in the NIR portion of the spectrum,
    using the cat_fit function of the spectral_analysis modules. Calculates
    also the equivalent width of the best fit using the ew_measurement task
    of the linestrength module.

    Returns:
    - Central wavelength of the lines
    - wavelength displacement between the central fitted lines and their real wavelength (in air)
    - equivalent widths measured in the best fitting model
    - broadening of the best fitted gaussians (in pixels and km/s). WARNING:
      this is not the velocity dispersion of the spectrum!

    """

    wavelength = params.wavelength
    flux = params.flux
    wave_interval_fit = params.wave_interval_fit
    wave_limits = params.wave_limits
    real_cat1 = params.real_cat1
    real_cat2 = params.real_cat2
    real_cat3 = params.real_cat3
    index_ca1 = params.index_ca1
    index_ca2 = params.index_ca2
    index_ca3 = params.index_ca3
    prev_spec = params.prev_spec
    prev_spec_nopath = params.prev_spec_nopath
    result_plot_dir = params.result_plot_dir
    task_done = params.task_done
    task_done2 = params.task_done2
    task_analysis = params.task_analysis
    wave_limits_cat = params.wave_limits_cat


    # 1) HEADER
    if event == 'Process all':
        task_done2 = 1
    else:
        task_done = 1
        task_analysis = 1

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_done2=task_done2, task_analysis=task_analysis)

    # Check wavelength limits
    #checking limits
    if min(wavelength) > wave_limits_cat[0] or max(wavelength) < wave_limits_cat[1]:
        if event == "Process all":
            print("The window band is out of the spectrum range")
        else:
            sg.popup("The window band is out of the spectrum range")
        return None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, params

    try:
        # Perform the fit
        cat_wave, cat_flux_norm, cat_fit, momentum = span.cat_fitting(wavelength, flux)

        # Identify minima in the fitted spectrum
        min_values_arg = argrelextrema(cat_fit, np.less)
        min_values_arg_extracted = min_values_arg[0]

        # Extract central wavelengths of the fitted CaT lines
        min_wave1 = round(cat_wave[min_values_arg_extracted[0]], 1) if len(min_values_arg_extracted) > 0 else 0.
        min_wave2 = round(cat_wave[min_values_arg_extracted[1]], 1) if len(min_values_arg_extracted) > 1 else 0.
        min_wave3 = round(cat_wave[min_values_arg_extracted[2]], 1) if len(min_values_arg_extracted) > 2 else 0.

        # Compute residuals
        residual_wave1 = min_wave1 - real_cat1
        residual_wave2 = min_wave2 - real_cat2
        residual_wave3 = min_wave3 - real_cat3

        # Convert residuals to velocity space
        c = 299792.458  # Speed of light in km/s
        delta_rv1 = (residual_wave1 / real_cat1) * c
        delta_rv2 = (residual_wave2 / real_cat2) * c
        delta_rv3 = (residual_wave3 / real_cat3) * c
        mean_delta_rv = (delta_rv1 + delta_rv2 + delta_rv3) / 3

        # Extract sigma values in pixels
        sigma_cat1 = momentum[3]
        sigma_cat2 = momentum[9]
        sigma_cat3 = momentum[15]

        # Convert sigma to velocity space
        sigma_cat1_vel = (sigma_cat1 / min_wave1) * c if min_wave1 else 0
        sigma_cat2_vel = (sigma_cat2 / min_wave2) * c if min_wave2 else 0
        sigma_cat3_vel = (sigma_cat3 / min_wave3) * c if min_wave3 else 0

        # Compute residual spectrum
        residual_flux = cat_flux_norm - cat_fit

        # Print results
        print("Central wavelength of the fitted CaT lines (A):", min_wave1, min_wave2, min_wave3)
        print("Real central wavelength of the CaT lines (A):", real_cat1, real_cat2, real_cat3)
        print("Residuals (A):", residual_wave1, residual_wave2, residual_wave3)
        print("Residuals in km/s:", delta_rv1, delta_rv2, delta_rv3)
        print("Mean Delta RV:", mean_delta_rv)
        print("Sigma (A):", sigma_cat1, sigma_cat2, sigma_cat3)
        print("Sigma (km/s):", sigma_cat1_vel, sigma_cat2_vel, sigma_cat3_vel)
        print('')

        #EW measurement of the fitted lines
        cat_single_index = True
        cat_plot = False
        cat_verbose = False
        cat_with_uncertainties = True
        cat_save_plot = False
        cat_normalise_spec = False

        #ca1
        id_array, ew_array_ca1, err_array, snr_ew_array,ew_array_mag, err_array_mag = ls.ew_measurement(cat_wave, cat_fit, index_ca1, cat_single_index, cat_plot, cat_verbose, cat_with_uncertainties, cat_save_plot, prev_spec, cat_normalise_spec, result_plot_dir)

        #ca2
        id_array, ew_array_ca2, err_array, snr_ew_array,ew_array_mag, err_array_mag = ls.ew_measurement(cat_wave, cat_fit, index_ca2, cat_single_index, cat_plot, cat_verbose, cat_with_uncertainties, cat_save_plot, prev_spec, cat_normalise_spec, result_plot_dir)

        #ca3
        id_array, ew_array_ca3, err_array, snr_ew_array,ew_array_mag, err_array_mag = ls.ew_measurement(cat_wave, cat_fit, index_ca3, cat_single_index, cat_plot, cat_verbose, cat_with_uncertainties, cat_save_plot, prev_spec, cat_normalise_spec, result_plot_dir)


        # Print EW measurements
        print("Equivalent Width from the Fit (Angstrom):")
        print("Ca1:", round(ew_array_ca1, 2), "Ca2:", round(ew_array_ca2, 2), "Ca3:", round(ew_array_ca3, 2))
        print('')

        # Plot fitting results
        if event == "Preview result" or (event == 'Process all' and save_plot):
            fig, (ax1, ax2) = plt.subplots(2, figsize=(8.5, 7), gridspec_kw={"height_ratios": [3, 1]})

            fig.suptitle("CaT Lines Fitting")
            ax1.plot(cat_wave, cat_flux_norm, label="Spectrum")
            ax1.plot(cat_wave, cat_fit, label="Fit line spec")
            ax1.set_xlabel("Wavelength (A)")
            ax1.set_ylabel("Norm flux")
            ax1.legend(fontsize=10)

            # Plot residuals
            ax2.plot(cat_wave, residual_flux, linewidth=0.5, label="Residuals")
            ax2.hlines(y=0, xmin=min(cat_wave), xmax=max(cat_wave), linestyles="--", lw=2, linewidth=2, color="r")
            ax2.set_xlabel("Wavelength (A)")
            ax2.set_ylabel("Residuals")
            ax2.legend(fontsize=10)

            if event == "Preview result":
                plt.show()
                plt.close()
            else:
                plt.savefig(result_plot_dir + '/'+ 'cat_fitting_' + prev_spec_nopath + '.png', format='png', dpi=300)
                plt.close()

        return min_wave1, min_wave2, min_wave3, residual_wave1, residual_wave2, residual_wave3, ew_array_ca1, ew_array_ca2, ew_array_ca3, real_cat1, real_cat2, real_cat3, delta_rv1, delta_rv2, delta_rv3, sigma_cat1, sigma_cat2, sigma_cat3, sigma_cat1_vel, sigma_cat2_vel, sigma_cat3_vel, params

    except Exception as e:
        if event == "Process all":
            print(f"CaT fitting failed. Try adjusting the parameters: {str(e)}")
        else:
            sg.popup(f"CaT fitting failed. Try adjusting the parameters: {str(e)}")
        return None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, params



def apply_line_fitting(event, save_plot, params):

    """
    Fits a single line (emission or absorption) with a gaussian function,
    using the cat_fit function of the spectral_analysis modules.

    Returns:
    - Central wavelength of the lines
    - wavelength displacement between the central fitted lines and their real wavelength (in air)
    - equivalent widths measured in the best fitting model
    - broadening of the best fitted gaussians (in pixels and km/s). WARNING:
      this is not the velocity dispersion of the spectrum!

    """

    wavelength = params.wavelength
    flux = params.flux
    wave_interval_fit = params.wave_interval_fit
    guess_param = params.guess_param
    emission_line = params.emission_line
    prev_spec = params.prev_spec
    prev_spec_nopath = params.prev_spec_nopath
    result_plot_dir = params.result_plot_dir
    task_done = params.task_done
    task_done2 = params.task_done2
    task_analysis = params.task_analysis


    # 1) HEADER
    if event == 'Process all':
        task_done2 = 1
    else:
        task_done = 1
        task_analysis = 1

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_done2=task_done2, task_analysis=task_analysis)

    # Check wavelength limits
    wave_limits = np.array([wavelength[0], wavelength[-1]])

    if min(wave_interval_fit) < wave_limits[0] or max(wave_interval_fit) > wave_limits[1]:
        if event == "Process all":
            print("The window band is out of the spectrum range")
        else:
            sg.popup("The window band is out of the spectrum range")
        return None, None, None, params

    try:
        # Perform the fitting
        line_wave, line_flux_norm, line_fit, momentum = span.line_fitting(wavelength, flux, wave_interval_fit, guess_param)

        # Find minimum or maximum value depending on absorption or emission
        if not emission_line:
            min_value_wave_arg = np.argmin(line_fit)
        else:
            min_value_wave_arg = np.argmax(line_fit)

        min_wave = line_wave[min_value_wave_arg]

        # Assign sigma of the line (in pixels)
        sigma_line = abs(momentum[3])

        # Convert sigma to velocity space
        c = 299792.458  # Speed of light in km/s
        sigma_line_vel = (sigma_line / min_wave) * c

        # Compute residual spectrum
        residual_flux = line_flux_norm - line_fit

        # Print results
        print("Central wavelength of the fitted line (A):", min_wave)
        print("Sigma (A):", sigma_line)
        print("Sigma (km/s):", sigma_line_vel)
        print('')
        # Plot fitting results
        if event == "Preview result" or (event == 'Process all' and save_plot):
            fig, (ax1, ax2) = plt.subplots(2, figsize=(8.5, 7), gridspec_kw={"height_ratios": [3, 1]})

            fig.suptitle("Line Fitting")
            ax1.plot(line_wave, line_flux_norm, label="Spectrum")
            ax1.plot(line_wave, line_fit, label="Fit line spec")
            ax1.set_xlabel("Wavelength (A)")
            ax1.set_ylabel("Norm flux")
            ax1.legend(fontsize=10)

            # Plot residuals
            ax2.plot(line_wave, residual_flux, linewidth=0.5, label="Residuals")
            ax2.hlines(y=0, xmin=min(line_wave), xmax=max(line_wave), linestyles="--", lw=2, linewidth=2, color="r")
            ax2.set_xlabel("Wavelength (A)")
            ax2.set_ylabel("Residuals")
            ax2.legend(fontsize=10)

            if event == "Preview result":
                plt.show()
                plt.close()
            else:
                plt.savefig(result_plot_dir + '/'+ 'line_fitting_' + prev_spec_nopath + '.png', format='png', dpi=300)
                plt.close()


        return min_wave, sigma_line, sigma_line_vel, params

    except Exception as e:
        if event == "Process all":
            print(f"Line fitting failed. Try adjusting the parameters: {str(e)}")
        else:
            sg.popup(f"Line fitting failed. Try adjusting the parameters: {str(e)}")
        return None, None, None, params



def apply_ppxf_kinematics(event, save_plot, params):

    """
    Fits the spectrum and retrieves the kinematics moments, using the pPXF algorithm via
    the ppxf_kinematics function of the spectral analysis module.

    Returns:
    - Kinematics moments
    - Kinematics moments formal errors
    - Bestfit template model
    - Number of kinematic components fitted
    - S/N of the fitted spectrum
    - Kinematics moments montecarlo errors (if calculated)

    """

    kin_stars_templates = getattr(params, 'kin_stars_templates', None)
    kin_lam_temp = getattr(params, 'kin_lam_temp', None)
    kin_velscale_templates = getattr(params, 'kin_velscale_templates', None)
    kin_FWHM_gal_cached = getattr(params, 'kin_FWHM_gal_cached', None)

    wavelength = params.wavelength
    flux = params.flux
    wave1_kin = params.wave1_kin
    wave2_kin = params.wave2_kin
    resolution_kin = params.resolution_kin
    constant_resolution_lambda = params.constant_resolution_lambda
    resolution_kin_r = params.resolution_kin_r
    resolution_kin_muse = params.resolution_kin_muse
    redshift_guess_kin = params.redshift_guess_kin
    sigma_guess_kin = params.sigma_guess_kin
    stellar_library_kin = params.stellar_library_kin
    additive_degree_kin = params.additive_degree_kin
    multiplicative_degree_kin = params.multiplicative_degree_kin
    kin_moments = params.kin_moments
    ppxf_kin_noise = params.ppxf_kin_noise
    gas_kin = params.gas_kin
    no_gas_kin = params.no_gas_kin
    kin_best_noise = params.kin_best_noise
    with_errors_kin = params.with_errors_kin
    ppxf_kin_custom_lib = params.ppxf_kin_custom_lib
    ppxf_kin_lib_folder = params.ppxf_kin_lib_folder
    ppxf_kin_custom_temp_suffix = params.ppxf_kin_custom_temp_suffix
    ppxf_kin_generic_lib = params.ppxf_kin_generic_lib
    ppxf_kin_generic_lib_folder = params.ppxf_kin_generic_lib_folder
    ppxf_kin_FWHM_tem_generic = params.ppxf_kin_FWHM_tem_generic
    ppxf_kin_fixed_kin = params.ppxf_kin_fixed_kin
    ppxf_kin_dust_gas = params.ppxf_kin_dust_gas
    ppxf_kin_dust_stars = params.ppxf_kin_dust_stars
    ppxf_kin_tie_balmer = params.ppxf_kin_tie_balmer
    ppxf_kin_two_stellar_components = params.ppxf_kin_two_stellar_components
    ppxf_kin_age_model1 = params.ppxf_kin_age_model1
    ppxf_kin_met_model1 = params.ppxf_kin_met_model1
    ppxf_kin_age_model2 = params.ppxf_kin_age_model2
    ppxf_kin_met_model2 = params.ppxf_kin_met_model2
    ppxf_kin_vel_model1 = params.ppxf_kin_vel_model1
    ppxf_kin_sigma_model1 = params.ppxf_kin_sigma_model1
    ppxf_kin_vel_model2 = params.ppxf_kin_vel_model2
    ppxf_kin_sigma_model2 = params.ppxf_kin_sigma_model2
    ppxf_kin_mask_emission = params.ppxf_kin_mask_emission
    ppxf_kin_have_user_mask = params.ppxf_kin_have_user_mask
    ppxf_kin_mask_ranges = params.ppxf_kin_mask_ranges
    ppxf_kin_mc_sim = params.ppxf_kin_mc_sim
    ppxf_kin_save_spectra = params.ppxf_kin_save_spectra
    prev_spec_nopath = params.prev_spec_nopath
    result_spec = params.result_spec
    result_plot_dir = params.result_plot_dir
    task_done = params.task_done
    task_done2 = params.task_done2
    task_analysis = params.task_analysis


    # 1) HEADER
    if event == 'Process all':
        task_done2 = 1
    else:
        task_done = 1
        task_analysis = 1

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_done2=task_done2, task_analysis=task_analysis)

    # Ensure wave1_kin < wave2_kin
    if wave1_kin > wave2_kin:
        wave1_kin, wave2_kin = wave2_kin, wave1_kin  # Swap values

    wave_limits_kin = (wave1_kin, wave2_kin)
    wave_limits = np.array([wavelength[0], wavelength[-1]])

    # Check wavelength limits
    if wave1_kin < wave_limits[0] or wave2_kin > wave_limits[1]:
        if event == "Process all":
            print("The window band is out of the spectrum range")
        else:
            sg.popup("The window band is out of the spectrum range")
        return None, None, None, None, None, None, None, None, None, None, None, params

    try:
        #fitting with ppxf
        if not ppxf_kin_fixed_kin or not gas_kin:
            kinematics, error_kinematics, bestfit_flux, bestfit_wavelength, bestfit_gas_flux, emission_corrected_flux, gas_without_continuum, kin_component, gas_component, snr_kin, error_kinematics_mc, gas_names, gas_flux, gas_flux_err, updated_templates, kin_lam_temp, kin_velscale_templates, kin_FWHM_gal_cached = span.ppxf_kinematics(wavelength, flux, wave1_kin, wave2_kin, resolution_kin, constant_resolution_lambda, resolution_kin_r, resolution_kin_muse, redshift_guess_kin, sigma_guess_kin, stellar_library_kin, additive_degree_kin, multiplicative_degree_kin, kin_moments, ppxf_kin_noise, gas_kin, no_gas_kin, kin_best_noise, with_errors_kin, ppxf_kin_custom_lib, ppxf_kin_lib_folder, ppxf_kin_custom_temp_suffix, ppxf_kin_generic_lib, ppxf_kin_generic_lib_folder, ppxf_kin_FWHM_tem_generic, ppxf_kin_dust_gas, ppxf_kin_dust_stars, ppxf_kin_tie_balmer, ppxf_kin_two_stellar_components, ppxf_kin_age_model1, ppxf_kin_met_model1, ppxf_kin_age_model2, ppxf_kin_met_model2, ppxf_kin_vel_model1, ppxf_kin_sigma_model1, ppxf_kin_vel_model2, ppxf_kin_sigma_model2, ppxf_kin_mask_emission, ppxf_kin_have_user_mask, ppxf_kin_mask_ranges, ppxf_kin_mc_sim, ppxf_kin_fixed_kin, stars_templates=kin_stars_templates, lam_temp = kin_lam_temp, velscale_cached = kin_velscale_templates, FWHM_gal_cached = kin_FWHM_gal_cached)

        if ppxf_kin_fixed_kin and gas_kin:

            # Fitst fit without gas
            kinematics, error_kinematics, bestfit_flux, bestfit_wavelength, bestfit_gas_flux, emission_corrected_flux, gas_without_continuum, kin_component, gas_component, snr_kin, error_kinematics_mc, gas_names, gas_flux, gas_flux_err, updated_templates, kin_lam_temp, kin_velscale_templates, kin_FWHM_gal_cached = span.ppxf_kinematics(wavelength, flux, wave1_kin, wave2_kin, resolution_kin, constant_resolution_lambda, resolution_kin_r, resolution_kin_muse, redshift_guess_kin, sigma_guess_kin, stellar_library_kin, additive_degree_kin, multiplicative_degree_kin, kin_moments, ppxf_kin_noise, False, True, kin_best_noise, with_errors_kin, ppxf_kin_custom_lib, ppxf_kin_lib_folder, ppxf_kin_custom_temp_suffix, ppxf_kin_generic_lib, ppxf_kin_generic_lib_folder, ppxf_kin_FWHM_tem_generic, ppxf_kin_dust_gas, ppxf_kin_dust_stars, ppxf_kin_tie_balmer, ppxf_kin_two_stellar_components, ppxf_kin_age_model1, ppxf_kin_met_model1, ppxf_kin_age_model2, ppxf_kin_met_model2, ppxf_kin_vel_model1, ppxf_kin_sigma_model1, ppxf_kin_vel_model2, ppxf_kin_sigma_model2, True, ppxf_kin_have_user_mask, ppxf_kin_mask_ranges, ppxf_kin_mc_sim, ppxf_kin_fixed_kin, stars_templates=kin_stars_templates, lam_temp = kin_lam_temp, velscale_cached = kin_velscale_templates, FWHM_gal_cached = kin_FWHM_gal_cached)

            kinematics_fixed = kinematics
            kin_stars_templates_gas = updated_templates #updating the template for the second fit to be like the stellar fit template

            # Second fit for the gas
            kinematics, error_kinematics, bestfit_flux, bestfit_wavelength, bestfit_gas_flux, emission_corrected_flux, gas_without_continuum, kin_component, gas_component, snr_kin, error_kinematics_mc, gas_names, gas_flux, gas_flux_err, updated_templates, kin_lam_temp, kin_velscale_templates, kin_FWHM_gal_cached = span.ppxf_kinematics(wavelength, flux, wave1_kin, wave2_kin, resolution_kin, constant_resolution_lambda, resolution_kin_r, resolution_kin_muse, redshift_guess_kin, sigma_guess_kin, stellar_library_kin, additive_degree_kin, multiplicative_degree_kin, kin_moments, ppxf_kin_noise, gas_kin, no_gas_kin, kin_best_noise, with_errors_kin, ppxf_kin_custom_lib, ppxf_kin_lib_folder, ppxf_kin_custom_temp_suffix, ppxf_kin_generic_lib, ppxf_kin_generic_lib_folder, ppxf_kin_FWHM_tem_generic, ppxf_kin_dust_gas, ppxf_kin_dust_stars, ppxf_kin_tie_balmer, ppxf_kin_two_stellar_components, ppxf_kin_age_model1, ppxf_kin_met_model1, ppxf_kin_age_model2, ppxf_kin_met_model2, ppxf_kin_vel_model1, ppxf_kin_sigma_model1, ppxf_kin_vel_model2, ppxf_kin_sigma_model2, False, ppxf_kin_have_user_mask, ppxf_kin_mask_ranges, ppxf_kin_mc_sim, ppxf_kin_fixed_kin, stars_templates=kin_stars_templates_gas, lam_temp = kin_lam_temp, velscale_cached = kin_velscale_templates, FWHM_gal_cached = kin_FWHM_gal_cached, kinematics_fixed = kinematics_fixed)


        if kin_stars_templates is None:
            params = replace(params, kin_stars_templates=updated_templates, kin_lam_temp=kin_lam_temp, kin_velscale_templates=kin_velscale_templates, kin_FWHM_gal_cached = kin_FWHM_gal_cached)

    #Saving the single stellar component fit results
        if kin_component == 0 and not ppxf_kin_two_stellar_components:
            vel = round(kinematics[0])
            sigma = round(kinematics[1])
            h3 = round(kinematics[2],3)
            h4 = round(kinematics[3],3)
            h5 = round(kinematics[4],3)
            h6 = round(kinematics[5],3)
            err_vel = round(error_kinematics[0])
            err_sigma = round(error_kinematics[1])
            err_h3 = round(error_kinematics[2],3)
            err_h4 = round(error_kinematics[3],3)
            err_h5 = round(error_kinematics[4],3)
            err_h6 = round(error_kinematics[5],3)

            vel_string = str(int(kinematics[0]))
            sigma_string = str(int(kinematics[1]))
            h3_string = str(round(kinematics[2],3))
            h4_string = str(round(kinematics[3],3))

            if with_errors_kin:
                err_rv_kin_mc, err_sigma_kin_mc, err_h3_kin_mc, err_h4_kin_mc, err_h5_kin_mc, err_h6_kin_mc = np.round(error_kinematics_mc[0],3)


        elif kin_component == 0 and ppxf_kin_two_stellar_components:
            vel1 = round(kinematics[0][0])
            sigma1 = round(kinematics[0][1])
            h31 = round(kinematics[0][2],3)
            h41 = round(kinematics[0][3],3)
            h51 = round(kinematics[0][4],3)
            h61 = round(kinematics[0][5],3)
            err_vel1 = round(error_kinematics[0][0])
            err_sigma1 = round(error_kinematics[0][1])
            err_h31 = round(error_kinematics[0][2],3)
            err_h41 = round(error_kinematics[0][3],3)
            err_h51 = round(error_kinematics[0][4],3)
            err_h61 = round(error_kinematics[0][5],3)

            vel2 = round(kinematics[1][0])
            sigma2 = round(kinematics[1][1])
            h32 = round(kinematics[1][2],3)
            h42 = round(kinematics[1][3],3)
            h52 = round(kinematics[1][4],3)
            h62 = round(kinematics[1][5],3)
            err_vel2 = round(error_kinematics[1][0])
            err_sigma2 = round(error_kinematics[1][1])
            err_h32 = round(error_kinematics[1][2],3)
            err_h42 = round(error_kinematics[1][3],3)
            err_h52 = round(error_kinematics[1][4],3)
            err_h62 = round(error_kinematics[1][5],3)

            vel_string1 = str(int(kinematics[0][0]))
            sigma_string1 = str(int(kinematics[0][1]))
            h3_string1 = str(round(kinematics[0][2],3))
            h4_string1 = str(round(kinematics[0][3],3))
            vel_string2 = str(int(kinematics[1][0]))
            sigma_string2 = str(int(kinematics[1][1]))
            h3_string2 = str(round(kinematics[1][2],3))
            h4_string2 = str(round(kinematics[1][3],3))

            if with_errors_kin:
                # extracting the MonteCarlo errors from the error array
                err_rv_kin_mc1, err_sigma_kin_mc1, err_h3_kin_mc1, err_h4_kin_mc1, err_h5_kin_mc1, err_h6_kin_mc1, err_rv_kin_mc2, err_sigma_kin_mc2, err_h3_kin_mc2, err_h4_kin_mc2, err_h5_kin_mc2, err_h6_kin_mc2  = np.round(error_kinematics_mc[0],3)


        else: # with gas
            vel = round(kinematics[0][0])
            sigma = round(kinematics[0][1])
            h3 = round(kinematics[0][2],3)
            h4 = round(kinematics[0][3],3)
            h5 = round(kinematics[0][4],3)
            h6 = round(kinematics[0][5],3)
            err_vel = round(error_kinematics[0][0])
            err_sigma = round(error_kinematics[0][1])
            err_h3 = round(error_kinematics[0][2],3)
            err_h4 = round(error_kinematics[0][3],3)
            err_h5 = round(error_kinematics[0][4],3)
            err_h6 = round(error_kinematics[0][5],3)

            vel_string = str(int(kinematics[0][0]))
            sigma_string = str(int(kinematics[0][1]))
            h3_string = str(round(kinematics[0][2],3))
            h4_string = str(round(kinematics[0][3],3))

            if with_errors_kin:
                err_rv_kin_mc, err_sigma_kin_mc, err_h3_kin_mc, err_h4_kin_mc, err_h5_kin_mc, err_h6_kin_mc = np.round(error_kinematics_mc[0],3)


        #plotting only in the preview and save plots mode
        if event == "Preview result" or (event == 'Process all' and save_plot):
            plt.plot(bestfit_wavelength, bestfit_flux)

            if kin_component == 0 and not ppxf_kin_two_stellar_components:
                plt.title('v = '+ vel_string + ' km/s.  Sigma = '+ sigma_string + ' km/s.  H3 = ' + h3_string + '.  H4 = ' + h4_string)
                if event == 'Preview result':
                    plt.show()
                    plt.close()
                else:
                    plt.savefig(result_plot_dir + '/'+ 'kin_ppxf_' + prev_spec_nopath + '.png', format='png', dpi=300)
                    plt.close()

            elif kin_component == 0 and ppxf_kin_two_stellar_components:
                plt.title('v1 = '+ vel_string1 + ' km/s.  Sigma1 = '+ sigma_string1 + ' km/s.  v2 = ' + vel_string2 + '. Sigma2 = '+ sigma_string2 + ' km/s')
                if event == 'Preview result':
                    plt.show()
                    plt.close()
                else:
                    plt.savefig(result_plot_dir + '/'+ 'kin_ppxf_' + prev_spec_nopath + '.png', format='png', dpi=300)

            else:
                plt.title('v = '+ vel_string + ' km/s.  Sigma = '+ sigma_string + ' km/s.  H3 = ' + h3_string + '.  H4 = ' + h4_string)
                if event == 'Preview result':
                    plt.show()
                    plt.close()
                else:
                    plt.savefig(result_plot_dir + '/'+ 'kin_ppxf_' + prev_spec_nopath + '.png', format='png', dpi=300)

        #saving best fit spec
        if (event == 'Process selected' or event == 'Process all') and ppxf_kin_save_spectra:
            file_fitted_kin = result_spec+'ppxf_kin_bestfit_' + prev_spec_nopath + '.fits'

            uti.save_fits_2d(bestfit_wavelength, bestfit_flux, file_fitted_kin)

            print ('Fit of the spectrum saved to: ', file_fitted_kin)
            print('')

            if gas_kin and ppxf_kin_save_spectra:
                try:
                    file_best_fit_gas = result_spec+'ppxf_kin_bestfit_gas_' + prev_spec_nopath + '.fits'
                    file_emission_corrected_spec = result_spec+'ppxf_kin_emission_corrected_' + prev_spec_nopath + '.fits'
                    file_gas_continuum_subtracted = result_spec+'ppxf_kin_gas_cont_subtracted_' + prev_spec_nopath + '.fits'


                    uti.save_fits_2d(bestfit_wavelength, bestfit_gas_flux, file_best_fit_gas)
                    uti.save_fits_2d(bestfit_wavelength, emission_corrected_flux, file_emission_corrected_spec)
                    uti.save_fits_2d(bestfit_wavelength, gas_without_continuum, file_gas_continuum_subtracted)

                    print ('Best fit gas model, emission corrected spectrum, and gas continuum subtracted spectrum saved to: ', result_spec)
                except Exception:
                    print ('Gas lines not found. I do not create gas spectra')

        plt.close() #closing all the plots

        return kinematics, error_kinematics, bestfit_flux, bestfit_wavelength, kin_component, gas_component, snr_kin, error_kinematics_mc, gas_names, gas_flux, gas_flux_err, params

    except Exception as e:
        if event == "Process all":
            print('Kinematics failed. Common cause: the templates do not cover the wavelength range you want to fit.\nOther possible explanations:\n- The resolution of your spectra is lower than the templates, if you are using the Xshooter templates\n- The templates do not exist. ')
        else:
            sg.popup('Kinematics failed. Common cause: the templates do not cover the wavelength range you want to fit.\nOther possible explanations:\n- The resolution of your spectra is lower than the templates, if you are using the Xshooter templates\n- The templates do not exist. ')
        return None, None, None, None, None, None, None, None, None, None, None, params



def apply_ppxf_stellar_populations(event, save_plot, params):

    """
    Fits the spectrum and retrieves the stellar populations parameters
    and non parametric SFH. Optionally, also the stellar parameters
    via Lick/IDS index measurement can be calculated. Uses the pPXF
    algorithm via the ppxf_pop function of the spectral analysis module
    and the ew_measurement function of the linestrength module.

    Returns:
    - Kinematics moments
    - Stellar popolation parameters (age, metallicity, alpha/Fe in available)
    - Bestfit template model
    - luminosity and mass weights collapsed to the metallicity dimension
    - Mass to light (if available)
    - S/N of the fitted spectrum
    - Errors estimated via bootstrap simulations (if calculated)
    - EW of the Lick/IDS indices used for stellar population studies (if calculated: Hbeta, Mgb, <Fe>, [MgFe]')
    - Stellar population parameters estimated with Lick/IDS analysis (if calculated)

    """

    wavelength = params.wavelength
    flux = params.flux
    wave1_pop = params.wave1_pop
    wave2_pop = params.wave2_pop
    res_pop = params.res_pop
    z_pop = params.z_pop
    sigma_guess_pop = params.sigma_guess_pop
    fit_components = params.fit_components
    with_errors = params.with_errors
    regul_err = params.regul_err
    additive_degree = params.additive_degree
    multiplicative_degree = params.multiplicative_degree
    ppxf_pop_tie_balmer = params.ppxf_pop_tie_balmer
    stellar_library = params.stellar_library
    ppxf_pop_dust_stars = params.ppxf_pop_dust_stars
    ppxf_pop_dust_gas = params.ppxf_pop_dust_gas
    ppxf_pop_noise = params.ppxf_pop_noise
    age_range_array = params.age_range_array
    met_range_array = params.met_range_array
    ppxf_pop_custom_lib = params.ppxf_pop_custom_lib
    ppxf_pop_lib_folder = params.ppxf_pop_lib_folder
    ppxf_pop_custom_npz = params.ppxf_pop_custom_npz
    ppxf_pop_npz_file = params.ppxf_pop_npz_file
    ppxf_pop_mask = params.ppxf_pop_mask
    ppxf_custom_temp_suffix = params.ppxf_custom_temp_suffix
    ppxf_best_param = params.ppxf_best_param
    ppxf_best_noise_estimate = params.ppxf_best_noise_estimate
    ppxf_frac_chi = params.ppxf_frac_chi
    ppxf_pop_convolve = params.ppxf_pop_convolve
    ppxf_pop_want_to_mask = params.ppxf_pop_want_to_mask
    ppxf_pop_mask_ranges = params.ppxf_pop_mask_ranges
    ppxf_pop_error_nsim = params.ppxf_pop_error_nsim
    ppxf_pop_lg_age = params.ppxf_pop_lg_age
    ppxf_pop_lg_met = params.ppxf_pop_lg_met
    stellar_parameters_lick_ppxf = params.stellar_parameters_lick_ppxf
    lick_index_file = params.lick_index_file
    sigma_lick_coeff_file = params.sigma_lick_coeff_file
    ssp_model_ppxf = params.ssp_model_ppxf
    interp_model_ppxf = params.interp_model_ppxf
    prev_spec_nopath = params.prev_spec_nopath
    prev_spec = params.prev_spec
    result_spec = params.result_spec
    result_ppxf_pop_data_dir = params.result_ppxf_pop_data_dir
    # result_ppxf_pop_data_dir_weights = params.result_ppxf_pop_data_dir_weights
    # result_ppxf_pop_data_dir_sfh = params.result_ppxf_pop_data_dir_sfh
    spectra_list_name = params.spectra_list_name
    result_plot_dir = params.result_plot_dir
    task_done = params.task_done
    task_done2 = params.task_done2
    task_analysis = params.task_analysis


    # 1) HEADER
    if event == 'Process all':
        task_done2 = 1
    else:
        task_done = 1
        task_analysis = 1

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_done2=task_done2, task_analysis=task_analysis)

    # Ensure wave1_pop < wave2_pop
    if wave1_pop > wave2_pop:
        wave1_pop, wave2_pop = wave2_pop, wave1_pop

    wave_limits = np.array([np.min(wavelength), np.max(wavelength)])

    # Check wavelength limits
    if wave1_pop < wave_limits[0] or wave2_pop > wave_limits[1]:
        if event == "Process all":
            print("The window band is out of the spectrum range")
        else:
            sg.popup("The window band is out of the spectrum range")
        return None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, params

    # Check if custom templates exist
    if ppxf_pop_custom_lib and not ppxf_pop_custom_npz:
        matching_temp = glob.glob(os.path.join(ppxf_pop_lib_folder, ppxf_custom_temp_suffix))
        if not matching_temp:
            if event == "Process all":
                print("The custom templates do not exist. Stopping")
            else:
                sg.popup("The custom templates do not exist. Stopping")
            return None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, params

    # Show error calculation message if enabled
    if with_errors:
        print("\n********* Calculating the errors for age and metallicity. Please sit and wait... *****\n")


    #setting up the other parameters required:
    if event == "Preview result":
        ppxf_with_plots = True #in the preview mode I show the plot
        ppxf_save_plot = False
        prev_spec_nopath = 'none'
    elif event  == 'Process all' and save_plot:
        ppxf_with_plots = False
        ppxf_save_plot = True

    else:
        ppxf_with_plots = False
        ppxf_save_plot = False

    try:
        kinematics, info_pop, info_pop_mass, mass_light, errors, galaxy, bestfit_flux, bestfit_wave, bestfit_flux_gas, residual_flux, chi_square, age_err_abs, met_err, alpha_err, mass_age_err_abs, mass_met_err, mass_alpha_err, emission_corrected_flux, pop_age, light_weights_age_bin, mass_weights_age_bin, cumulative_mass, light_weights_age_std, mass_weights_age_std, cumulative_light_std, cumulative_mass_std, snr_pop, light_weights, mass_weights, t50_age, t80_age, t50_cosmic, t80_cosmic = span.ppxf_pop(wavelength, flux, wave1_pop, wave2_pop, res_pop, z_pop, sigma_guess_pop, fit_components, ppxf_with_plots, with_errors, ppxf_save_plot, prev_spec_nopath, regul_err, additive_degree, multiplicative_degree, ppxf_pop_tie_balmer, stellar_library, ppxf_pop_dust_stars, ppxf_pop_dust_gas, ppxf_pop_noise, age_range_array, met_range_array, ppxf_pop_custom_lib, ppxf_pop_lib_folder, ppxf_pop_custom_npz, ppxf_pop_npz_file, ppxf_pop_mask, ppxf_custom_temp_suffix, ppxf_best_param, ppxf_best_noise_estimate, ppxf_frac_chi, ppxf_pop_convolve, ppxf_pop_want_to_mask, ppxf_pop_mask_ranges, ppxf_pop_error_nsim, ppxf_pop_lg_age, ppxf_pop_lg_met, result_plot_dir)


        age_ssp =0
        met_ssp = 0
        alpha_ssp = 0
        err_age_ssp = 0
        err_met_ssp = 0
        err_alpha_ssp = 0
        ppxf_lick_params = np.zeros(6)
        ssp_lick_indices_ppxf = np.zeros(5)
        ssp_lick_indices_err_ppxf = np.zeros(5)

        # storing the values in the dataframes and save to disc
        if stellar_library == 'sMILES' and not ppxf_pop_custom_lib:
            alpha = info_pop[2]
            mass_alpha = info_pop_mass[2]

        #If the stellar parameters also with Lick/IDS indices is activated:
        if stellar_parameters_lick_ppxf:

            #Extracting the kinematics from the pPXF fit. THe walues to extract depends whether I have gas or not
            try:
                num_comp_kinematics = len(kinematics)
                kin_stars = np.array(kinematics[0])
                dop_vel_pop_ppxf = kin_stars[0]
                sigma_pop_ppxf = kin_stars[1]

            except (ValueError, IndexError):
                num_comp_kinematics = 0
                kin_stars = kinematics
                dop_vel_pop_ppxf = kinematics[0]
                sigma_pop_ppxf = kinematics[1]

            #Extracting the wave (in A) and emission corrected flux from the ppxf fit
            lick_wavelength_ppxf = bestfit_wave
            lick_flux_ppxf = emission_corrected_flux
            lick_step_ppxf = lick_wavelength_ppxf[1] - lick_wavelength_ppxf[0]

            #rebinning linear
            lick_wavelength_ppxf, lick_flux_ppxf, npoint_resampled = spman.resample(lick_wavelength_ppxf, lick_flux_ppxf, lick_step_ppxf)
            lick_constant_fwhm_ppxf = True

            #Doppler correction from the velociy calculated by pPXF
            lick_wavelength_ppxf, lick_flux_ppxf = spman.dopcor(lick_wavelength_ppxf, lick_flux_ppxf, dop_vel_pop_ppxf, True) #doppler correction. The cosmological z correction has been already performed by the span.ppxf_pop function, prior to really run pPXF. Here I correct only for the real velocity component measured by the fit.

            # 3) degrading the resolution, only if smaller than the lick system
            if res_pop < 8.4:
                lick_degraded_wavelength_ppxf, lick_degraded_flux_ppxf = spman.degrade_to_lick(lick_wavelength_ppxf, lick_flux_ppxf, res_pop, lick_constant_fwhm_ppxf)
            else:
                print('WARNING: The resolution of the spectrum is smaller than the one needed for the Lick/IDS system. I will still calculate the Lick/IDS indices but the results might be inaccurate.')
                lick_degraded_wavelength_ppxf = lick_wavelength_ppxf
                lick_degraded_flux_ppxf = lick_flux_ppxf

            # 4) Measuring the EW and doing plot
            lick_single_index_ppxf = False
            lick_ew_plot_ppxf = False
            lick_verbose_ppxf = False
            lick_with_uncertainties_ppxf = True
            lick_save_plot_ppxf = False
            lick_normalise_spec_ppxf = True

            lick_id_array_ppxf, lick_ew_array_ppxf, lick_err_array_ppxf, lick_snr_ew_array_ppxf, lick_ew_array_mag_ppxf, lick_err_array_mag_ppxf = ls.ew_measurement(lick_degraded_wavelength_ppxf, lick_degraded_flux_ppxf, lick_index_file, lick_single_index_ppxf, lick_ew_plot_ppxf, lick_verbose_ppxf, lick_with_uncertainties_ppxf, lick_save_plot_ppxf, prev_spec, lick_normalise_spec_ppxf, result_plot_dir)

            # 5) Correcting the EWs for sigma by extracting the sigma determined by pPXF. The location changes whether fitting with of tiwthout gas
            sigma_to_correct_lick_ppxf = sigma_pop_ppxf

            #now correcting for real
            corrected_lick_ew_array_ppxf, corrected_lick_err_array_ppxf, corrected_lick_ew_mag_array_ppxf, corrected_lick_err_mag_array_ppxf = ls.corr_ew_lick(lick_ew_array_ppxf, lick_err_array_ppxf, lick_ew_array_mag_ppxf, sigma_lick_coeff_file, sigma_to_correct_lick_ppxf)

            #uodating the values
            lick_ew_array_ppxf = corrected_lick_ew_array_ppxf
            lick_err_array_ppxf = corrected_lick_err_array_ppxf
            lick_ew_array_mag_ppxf = corrected_lick_ew_mag_array_ppxf
            lick_err_array_mag_ppxf = corrected_lick_err_mag_array_ppxf

            #assigning meaningful names to the indices used for stellar populations and creating the combined ones
            Hbeta_ppxf_single = lick_ew_array_ppxf[0]
            Hbetae_ppxf_single = lick_err_array_ppxf[0]
            Mg2_ppxf_single = lick_ew_array_mag_ppxf[1]
            Mg2e_ppxf_single = lick_err_array_mag_ppxf[1]
            Mgb_ppxf_single = lick_ew_array_ppxf[2]
            Mgbe_ppxf_single = lick_err_array_ppxf[2]
            Fe5270_ppxf_single = lick_ew_array_ppxf[3]
            Fe5270e_ppxf_single = lick_err_array_ppxf[3]
            Fe5335_ppxf_single = lick_ew_array_ppxf[4]
            Fe5335e_ppxf_single = lick_err_array_ppxf[4]
            Fem_ppxf_single = (Fe5270_ppxf_single+Fe5335_ppxf_single)/2
            Feme_ppxf_single = np.sqrt((0.5*Fe5270e_ppxf_single)**2+(0.5*Fe5335e_ppxf_single)**2)
            MgFe_ppxf_single = (np.sqrt(Mgb_ppxf_single*(0.72*Fe5270_ppxf_single+0.28*Fe5335_ppxf_single)))
            MgFe_ppxf_single = np.nan_to_num(MgFe_ppxf_single, nan=0)
            MgFee_ppxf_single = np.sqrt((((Fe5270_ppxf_single*18/25+Fe5335_ppxf_single*7/25)/(2*np.sqrt(Mgb_ppxf_single*(Fe5270_ppxf_single*18/25+Fe5335_ppxf_single*7/25))))*Mgbe_ppxf_single)**2+((Mgb_ppxf_single*18/25/(2*np.sqrt(Mgb_ppxf_single*(Fe5270_ppxf_single*18/25+Fe5335_ppxf_single*7/25))))*Fe5270e_ppxf_single)**2+((Mgb_ppxf_single*7/25/(2*np.sqrt(Mgb_ppxf_single*(Fe5270_ppxf_single*18/25+Fe5335_ppxf_single*7/25))))*Fe5335e_ppxf_single)**2)
            MgFee_ppxf_single = np.nan_to_num(MgFee_ppxf_single, nan=0)

            ssp_lick_indices_list_ppxf = np.column_stack((Hbeta_ppxf_single, MgFe_ppxf_single, Fem_ppxf_single, Mgb_ppxf_single))
            ssp_lick_indices_ppxf = ssp_lick_indices_list_ppxf.reshape(-1)
            ssp_lick_indices_err_list_ppxf = np.column_stack((Hbetae_ppxf_single, MgFee_ppxf_single, Feme_ppxf_single, Mgbe_ppxf_single))
            ssp_lick_indices_err_ppxf = ssp_lick_indices_err_list_ppxf.reshape(-1)

            #Determining the stellar parameters
            age_ssp, met_ssp, alpha_ssp, err_age_ssp, err_met_ssp, err_alpha_ssp = span.lick_pop(ssp_lick_indices_ppxf, ssp_lick_indices_err_ppxf, ssp_model_ppxf, interp_model_ppxf)

            #putting in an array:
            ppxf_lick_params = np.array([age_ssp, met_ssp, alpha_ssp, err_age_ssp, err_met_ssp, err_alpha_ssp])

            print ('')
            print ('Age (Gyr):', round(age_ssp,2), '+/-', round(err_age_ssp, 2))
            print ('[M/H] (dex):', round(met_ssp, 2), '+/-', round(err_met_ssp,2))
            print ('[Alpha/Fe]:', round(alpha_ssp, 2), '+/-', round(err_alpha_ssp,2))
            print('')

        #saving the file with the fit
        if event == 'Process selected' or event == 'Process all':
            try:
                #saving the spectra
                file_fit_pop = result_spec+'ppxf_fit_pop_residuals_' + prev_spec_nopath + '.fits'
                file_fit_stellar_template = result_spec+'ppxf_fit_pop_stellar_template_' + prev_spec_nopath + '.fits'
                file_spec_emission_corrected = result_spec+'ppxf_fit_pop_emission_corrected_' + prev_spec_nopath + '.fits'

                #saving the SFH and weights in specific subfolders
                result_ppxf_pop_data_dir_weights = os.path.join(result_ppxf_pop_data_dir, 'weights')
                result_ppxf_pop_data_dir_sfh = os.path.join(result_ppxf_pop_data_dir, 'SFH')

                os.makedirs(result_ppxf_pop_data_dir_weights, exist_ok=True)
                os.makedirs(result_ppxf_pop_data_dir_sfh, exist_ok=True)

                file_sfh = result_ppxf_pop_data_dir_sfh+'/'+spectra_list_name+'_ppxf_fit_pop_SFH_' + prev_spec_nopath + '.dat'
                file_all_light_weights = result_ppxf_pop_data_dir_weights+'/'+spectra_list_name+'_ppxf_fit_pop_light_weights_' + prev_spec_nopath + '.dat'
                file_all_mass_weights = result_ppxf_pop_data_dir_weights+'/'+spectra_list_name+'_ppxf_fit_pop_mass_weights_' + prev_spec_nopath + '.dat'

                bestfit_wave = bestfit_wave

                #in case I don't have gas
                try:
                    #saving the SFH with lg ages or with linear ages
                    lum_cumulative = np.cumsum(light_weights_age_bin)
                    mass_cumulative = np.cumsum(mass_weights_age_bin)
                    if ppxf_pop_lg_age:

                        if with_errors:
                            np.savetxt(file_sfh, np.column_stack([pop_age, light_weights_age_bin, mass_weights_age_bin, lum_cumulative, mass_cumulative, light_weights_age_std, mass_weights_age_std, cumulative_mass_std]), header="lg_age(dex)\tlum_fraction\tmass_fraction\tcumulative_lum\tcumulative_mass\terr_lum\terr_mass\terr_cumul", delimiter='\t')
                        else:
                            np.savetxt(file_sfh, np.column_stack([pop_age, light_weights_age_bin, mass_weights_age_bin, lum_cumulative, mass_cumulative]), header="lg_age(dex)\tlum_fraction\tmass_fraction\tcumulative_lum\tcumulative_mass", delimiter='\t')


                    else:
                        if with_errors:
                            np.savetxt(file_sfh, np.column_stack([pop_age, light_weights_age_bin, mass_weights_age_bin, lum_cumulative, mass_cumulative, light_weights_age_std, mass_weights_age_std, cumulative_mass_std]), header="age(Gyr)\tlum_fraction\tmass_fraction\tcumulative_lum\tcumulative_mass\terr_lum\terr_mass\terr_cumul", delimiter='\t')
                        else:
                            np.savetxt(file_sfh, np.column_stack([pop_age, light_weights_age_bin, mass_weights_age_bin, lum_cumulative, mass_cumulative]), header="age(Gyr)\tlum_fraction\tmass_fraction\tcumulative_lum\tcumulative_mass", delimiter='\t')

                    print ('File containing the luminosity and mass SFH saved: ', file_sfh)

                    #saving the light weights
                    np.savetxt(file_all_light_weights, light_weights.reshape(-1, light_weights.shape[-1]), fmt="%.8e", delimiter=' ', header="Light weights")
                    print ('File containing the light weights saved: ', file_all_light_weights)

                    #saving the mass weights
                    np.savetxt(file_all_mass_weights, mass_weights.reshape(-1, mass_weights.shape[-1]), fmt="%.8e", delimiter=' ', header="Mass weights")
                    print ('File containing the mass weights saved: ', file_all_mass_weights)

                    if bestfit_flux_gas == 0:

                        #saving the residual file
                        uti.save_fits_2d(bestfit_wave, residual_flux, file_fit_pop) #using the save_fits_2d function because the wavelength sampling is not linear

                        #saving the template without gas:
                        uti.save_fits_2d(bestfit_wave, bestfit_flux, file_fit_stellar_template)

                        print ('File containing the residuals of the fit saved: ', file_fit_pop)
                        print ('File containing the stellar fitted template: ', file_fit_stellar_template)
                        print('')

                except ValueError: #considering also the gas template if I receive this error

                    #saving the residual file
                    uti.save_fits_2d(bestfit_wave, residual_flux, file_fit_pop)

                    #saving the best template without gas:
                    stellar_fit_flux = bestfit_flux-bestfit_flux_gas
                    uti.save_fits_2d(bestfit_wave, stellar_fit_flux, file_fit_stellar_template)

                    #saving the emission corrected spectra in linear step
                    uti.save_fits_2d(bestfit_wave, emission_corrected_flux, file_spec_emission_corrected)

                    print ('File containing the residuals of the fit saved: ', file_fit_pop)
                    print ('File containing the stellar fitted template: ', file_fit_stellar_template)
                    print ('File containing the emission corrected spectra (empty if no gas selected): ', file_spec_emission_corrected)
                    print('')
            except TypeError:
                print ('Something went wrong')

        return kinematics, info_pop, info_pop_mass, mass_light, chi_square, met_err, mass_met_err, snr_pop, ppxf_pop_lg_age, ppxf_pop_lg_met, age_err_abs, mass_age_err_abs, alpha_err, mass_alpha_err, t50_age, t80_age, t50_cosmic, t80_cosmic, ssp_lick_indices_ppxf, ssp_lick_indices_err_ppxf, ppxf_lick_params, params

    except Exception:
        if event == "Process all":
            print('Stellar populations and SFH failed. Common cause: The templates do not cover the wavelength range you want to fit.\nOther possible expalations:\n- Check the age-metallicity range\n- If you used custom templates, check if they are on a regular age and metallicity grid.\n Skipping...')
        else:
            sg.popup('Stellar populations and SFH failed. Common cause: The templates do not cover the wavelength range you want to fit.\nOther possible expalations:\n- Check the age-metallicity range\n- If you used custom templates, check if they are on a regular age and metallicity grid.\n Skipping...')
        return None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, params
