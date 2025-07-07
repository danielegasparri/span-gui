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

# Functions to save and restore the parameters, for cleaning the tasks and generating a spectra list

try: #try local import if executed as script
    #GUI import
    from params import SpectraParams
    from FreeSimpleGUI_local import FreeSimpleGUI as sg
    from span_functions import system_span as stm
    from span_modules import misc

except ModuleNotFoundError: #local import if executed as package
    #GUI import
    from .params import SpectraParams
    from span.FreeSimpleGUI_local import FreeSimpleGUI as sg
    from span.span_functions import system_span as stm
    from . import misc

import json
import numpy as np
from dataclasses import replace
import os


def generate_spectra_list(window, params):

    """
    Opens a dialog to select a folder containing spectra and generates a list file.

    Parameters:
    - window: PySimpleGUI window object
    - params: NamedTuple containing the program's operational parameters

    Returns:
    - Updated params with the new spectra list file

    """
    layout, scale_win, fontsize, default_size = misc.get_layout()
    sg.theme('LightBlue')
    list_layout = [
        [sg.Text("Select the folder with the spectra:")],
        [sg.InputText(key='-FOLDER-'), sg.FolderBrowse()],
        [sg.Button('Save')]
    ]

    list_window = sg.Window('Generate spectra list containing 1D spectra', list_layout)

    while True:
        list_event, list_values = list_window.read()

        if list_event == sg.WINDOW_CLOSED:
            break
        elif list_event == 'Save':
            folder_path = list_values['-FOLDER-']
            if folder_path:
                file_list = stm.get_files_in_folder(folder_path)
                output_file = os.path.basename(os.path.normpath(folder_path)) + '_spectra_list.txt'
                stm.save_to_text_file(file_list, output_file)

                sg.Popup('Spectra file list saved in the current working directory',
                         output_file, 'You can now load this list file')

                # Updating the spectra list
                window['spec_list'].update(output_file)
                params = replace(params, spectra_list=output_file)
                params = replace(params, spectra_list_name=os.path.splitext(os.path.basename(params.spectra_list))[0])

    list_window.close()

    return params


def clear_all_tasks(window, params):

    """
    Clears all selected tasks in the GUI and resets the parameters.

    Parameters:
    - window: FreeSimpleGUI window object
    - params: NamedTuple containing the program's operational parameters

    Returns:
    - Updated params with default values

    """

    # Reset checkboxes in the GUI
    window['show_hdr'].Update(value=False)
    window['show_step'].Update(value=False)
    window['show_res'].Update(value=False)
    window['convert_spec'].Update(value=False)
    window['compare_spec'].Update(value=False)
    window['convert_flux'].Update(value=False)
    window['show_snr'].Update(value=False)

    # Reset Spectra Manipulation Panel
    params = replace(params,
                     cropping_spectrum=False,
                     sigma_clipping=False,
                     wavelet_cleaning=False,
                     filter_denoise=False,
                     dop_cor=False,
                     helio_corr=False,
                     rebinning=False,
                     rebinning_log=False,
                     rebinning_linear=True,
                     degrade=False,
                     normalize_wave=False,
                     sigma_broad=False,
                     add_noise=False,
                     continuum_sub=False,
                     average_all=False,
                     norm_and_average=False,
                     do_nothing=True,
                     sum_all=False,
                     normalize_and_sum_all=False,
                     use_for_spec_an=False,
                     subtract_normalized_avg=False,
                     subtract_normalized_spec=False,
                     add_pedestal=False,
                     multiply=False,
                     derivatives=False,
                     reorder_op=False,
                     active_operations=[],
                     reordered_operations=[],
                     current_order=None,
                     save_intermediate_spectra=True,
                     save_final_spectra=False,
                     not_save_spectra=False)

    # Reset Spectral Analysis Frame
    window['bb_fitting'].Update(value=False)
    window['xcorr'].Update(value=False)
    window['sigma_measurement'].Update(value=False)
    window['ew_measurement'].Update(value=False)
    window['line_fitting'].Update(value=False)
    window['ppxf_kin'].Update(value=False)
    window['ppxf_pop'].Update(value=False)
    window['save_plots'].Update(value=False)

    print('All tasks cleared')

    return params


# save all the settings in a JSON file
def save_settings(filename, keys, events, values, params: SpectraParams):
    data = {
        'keys': keys,
        'events': events,
        'values': {
            # Main GUI values
            'spec_list': values.get('spec_list', []),
            'wave_units_nm': values.get('wave_units_nm', False),
            'wave_units_a': values.get('wave_units_a', False),
            'wave_units_mu': values.get('wave_units_mu', False),
            'show_hdr': values.get('show_hdr', False),
            'show_step': values.get('show_step', False),
            'show_res': values.get('show_res', False),
            'lambda_res_left': values.get('lambda_res_left', 0.0),
            'lambda_res_right': values.get('lambda_res_right', 0.0),
            'convert_spec': values.get('convert_spec', False),
            'convert_to_txt': values.get('convert_to_txt', False),
            'convert_to_fits': values.get('convert_to_fits', False),
            'compare_spec': values.get('compare_spec', False),
            'spec_to_compare': values.get('spec_to_compare', ""),
            'convert_flux': values.get('convert_flux', False),
            'convert_to_fnu': values.get('convert_to_fnu', False),
            'convert_to_fl': values.get('convert_to_fl', False),
            'show_snr': values.get('show_snr', False),
            'wave_snr': values.get('wave_snr', 0.0),
            'delta_wave_snr': values.get('delta_wave_snr', 0.0),

            'bb_fitting': values['bb_fitting'],
            'xcorr': values['xcorr'],
            'sigma_measurement': values['sigma_measurement'],
            'ew_measurement': values['ew_measurement'],
            'line_fitting': values['line_fitting'],
            'ppxf_kin': values['ppxf_kin'],
            'ppxf_pop': values['ppxf_pop'],

            'save_plots': values['save_plots'],
            },

        'params': {
            # Params
            'cropping': params.cropping_spectrum,
            'cropping_low_wave': params.cropping_low_wave,
            'cropping_high_wave': params.cropping_high_wave,
            'sigma_clip': params.sigma_clipping,
            'clip_factor': params.clip_factor,
            'res_spec_for_sigma_clip': params.sigma_clip_resolution,
            'single_vel_clip' : params.sigma_clip_single_vel,
            'clip_to_vel': params.sigma_clip_single_value,
            'file_for_clip': params.sigma_clip_have_file,
            'sigma_clip_file': params.sigma_clip_sigma_file,
            'wavelet_cleaning': params.wavelet_cleaning,
            'sigma_wavelets': params.sigma_wavelets,
            'wavelets_layers': params.wavelets_layers,
            'dopcor': params.dop_cor,

            'dopcor_value': params.dop_cor_single_shot_vel,
            'dop_cor_have_vel': params.dop_cor_have_vel,
            'dop_cor_have_z': params.dop_cor_have_z,
            'file_for_dopcor': params.dop_cor_have_file,
            'dopcor_file': params.dop_cor_file,
            'dopcor_single_value': params.dop_cor_single_shot,
            'helio_corr': params.helio_corr,
            'file_for_helio': params.helio_have_file,
            'helio_file': params.helio_file,
            'helio_single_value': params.helio_single_shot,
            'helio_location': params.helio_single_shot_location,
            'helio_date': params.helio_single_shot_date,
            'helio_ra': params.ra_obj,
            'helio_dec': params.dec_obj,

            'rebin': params.rebinning,
            'rebin_pix_lin': params.rebinning_linear,
            'rebin_step_pix': params.rebin_step_pix,
            'rebin_sigma_lin': params.rebinning_log,
            'rebin_step_sigma': params.rebin_step_sigma,

            'degrade_resolution': params.degrade,
            'is_initial_res_r': params.is_initial_res_r,
            'degrade_from_r': params.initial_res_r,
            'res_degrade_to_r': params.res_degrade_to_r,
            'degrade_to_r': params.final_res_r,
            'res_degrade_to_fwhm': params.res_degrade_to_fwhm,
            'final_res_r_to_fwhm': params.final_res_r_to_fwhm,
            'is_initial_res_fwhm': params.is_initial_res_fwhm,
            'degrade_from_l': params.initial_res_fwhm,
            'degrade_to_l': params.final_res_fwhm,
            'res_degrade_muse': params.res_degrade_muse,
            'res_degrade_muse_value': params.res_degrade_muse_value,

            'norm_spec': params.normalize_wave,
            'norm_wave': params.norm_wave,

            'cont_sub': params.continuum_sub,
            'cont_model_filtering' : params.cont_model_filtering,
            'cont_model_poly' : params.cont_model_poly,
            'markers_cont_operations' : params.cont_math_operation,
            'cont_want_to_mask' : params.cont_want_to_mask,
            'cont_mask_ranges' : params.cont_mask_ranges_str,
            'cont_poly_degree' : params.cont_poly_degree,

            'broadening_spec': params.sigma_broad,
            'sigma_to_add': params.sigma_to_add,
            'add_noise': params.add_noise,
            'noise_to_add': params.noise_to_add,

            'filter_denoise' : params.filter_denoise,
            'moving_average' : params.moving_average,
            'box_moving_avg' : params.box_moving_avg,
            'box_moving_avg_size' : params.box_moving_avg_size,
            'gauss_moving_avg' : params.gauss_moving_avg,
            'gauss_moving_avg_kernel' : params.gauss_moving_avg_kernel,
            'low_pass_filter' : params.low_pass_filter,
            'lowpass_cut_off' : params.lowpass_cut_off,
            'lowpass_order' : params.lowpass_order,
            'bandpass_filter' : params.bandpass_filter,
            'bandpass_lower_cut_off' : params.bandpass_lower_cut_off,
            'bandpass_upper_cut_off' : params.bandpass_upper_cut_off,
            'bandpass_order' : params.bandpass_order,
            'avg_all': params.average_all,

            'norm_avg_all': params.norm_and_average,
            'none': params.do_nothing,
            'sum_all': params.sum_all,
            'norm_sum_all': params.normalize_and_sum_all,
            'use_for_spec_an': params.use_for_spec_an,
            'subtract_norm_avg': params.subtract_normalized_avg,
            'subtract_norm_spec': params.subtract_normalized_spec,
            'spec_to_subtract': params.spectra_to_subtract,
            'add_pedestal': params.add_pedestal,
            'pedestal_to_add': params.pedestal_to_add,
            'multiply': params.multiply,
            'multiply_factor': params.multiply_factor,
            'derivatives' : params.derivatives,
            'reorder_op' : params.reorder_op,
            'current_order': params.current_order,
            'active_operations': params.active_operations,
            'reordered_operations': params.reordered_operations,
            'save_intermediate_spectra': params.save_intermediate_spectra,
            'save_final_spectra': params.save_final_spectra,
            'not_save_spectra': params.not_save_spectra,


            'sigma_coeff': params.sigma_coeff,
            'sigma_corr': params.sigma_corr,


            #Blackbody parameters
            'left_wave_bb': params.wave1_bb,
            'right_wave_bb': params.wave2_bb,
            't_guess_bb': params.t_guess,

            #Cross-corr parameters
            'xcorr_template': params.template_crosscorr,
            'xcorr_template_wave_nm': params.lambda_units_template_crosscorr_nm,
            'xcorr_template_wave_a': params.lambda_units_template_crosscorr_a,
            'xcorr_template_wave_mu': params.lambda_units_template_crosscorr_mu,
            'xcorr_smooth_template': params.smooth_template_crosscorr,
            'xcorr_smooth_template_value': params.smooth_value_crosscorr,
            'xcorr_left_lambda': params.low_wave_corr,
            'xcorr_right_lambda': params.high_wave_corr,
            'is_vel_xcorr': params.is_vel_xcorr,
            'is_z_xcorr': params.is_z_xcorr,
            'xcorr_low_vel': params.low_vel_corr,
            'xcorr_high_vel': params.high_vel_corr,
            'low_z_corr': params.low_z_corr,
            'high_z_corr': params.high_z_corr,
            'xcorr_limit_wave_range': params.xcorr_limit_wave_range,
            'xcorr_vel_step': params.xcorr_vel_step,
            'xcorr_z_step': params.xcorr_z_step,

            #Velocity dispersion parameters
            'template_sigma': params.template_sigma,
            'lambda_units_template_sigma_nm': params.lambda_units_template_sigma_nm,
            'lambda_units_template_sigma_a': params.lambda_units_template_sigma_a,
            'lambda_units_template_sigma_mu': params.lambda_units_template_sigma_mu,
            'band_cat': params.band_cat,
            'band_halpha': params.band_halpha,
            'band_nad': params.band_nad,
            'band_h': params.band_h,
            'band_k': params.band_k,
            'resolution_spec': params.resolution_spec,
            'resolution_template': params.resolution_template,
            'band_custom': params.band_custom,
            'low_wave_sigma': params.low_wave_sigma,
            'high_wave_sigma': params.high_wave_sigma,
            'low_wave_cont': params.low_wave_cont,
            'high_wave_cont': params.high_wave_cont,

            #Line-strength parameters
            'ew_idx_file': params.have_index_file,
            'idx_file': params.index_file,
            'single_index': params.single_index,
            'left_wave_blue_cont': params.idx_left_blue,
            'right_wave_blue_cont': params.idx_right_blue,
            'left_wave_red_cont': params.idx_left_red,
            'right_wave_red_cont': params.idx_right_red,
            'left_line': params.idx_left_line,
            'right_line': params.idx_right_line,
            'ew_lick': params.lick_ew,
            'lick_constant_fwhm': params.lick_constant_fwhm,
            'spec_lick_res_fwhm': params.spec_lick_res_fwhm,
            'lick_constant_r': params.lick_constant_r,
            'spec_lick_res_r': params.spec_lick_res_r,
            'lick_correct_emission': params.lick_correct_emission,
            'z_guess_lick_emission': params.z_guess_lick_emission,
            'dop_correction_lick': params.dop_correction_lick,
            'correct_ew_sigma': params.correct_ew_sigma,
            'radio_lick_sigma_auto': params.radio_lick_sigma_auto,
            'radio_lick_sigma_single': params.radio_lick_sigma_single,
            'sigma_single_lick': params.sigma_single_lick,
            'radio_lick_sigma_list': params.radio_lick_sigma_list,
            'sigma_lick_file': params.sigma_lick_file,
            'stellar_parameters_lick': params.stellar_parameters_lick,
            'ssp_model': params.ssp_model,
            'interp_model': params.interp_model,
            'ew_corr_idx_file': params.have_index_file_corr,
            'idx_corr_file': params.index_file_corr,
            'ew_corr_single_idx': params.single_index_corr,
            # 'stellar_spectra_coeff_file': params.stellar_spectra_coeff_file,
            'sigma_coeff_sample_list': params.stellar_spectra_coeff_file,
            'sigma_coeff_sample_list_wave_nm': params.lambda_units_coeff_nm,
            'sigma_coeff_sample_list_wave_a': params.lambda_units_coeff_a,
            'sigma_coeff_sample_list_wave_mu': params.lambda_units_coeff_mu,
            'sigma_coeff_sample_smooth': params.smooth_stellar_sample,
            'sigma_coeff_sample_smooth_sigma': params.smooth_value_sample,
            'sigma_file': params.sigma_vel_file,
            'ew_file_to_correct': params.ew_list_file,
            'coeff_sigma_file': params.sigma_coeff_file,

            #Line(s) fitting parameters
            'cat_fit': params.cat_band_fit,
            'line_fit_single': params.usr_fit_line,
            'emission_line': params.emission_line,
            'left_wave_fitting': params.low_wave_fit,
            'right_wave_fitting': params.high_wave_fit,
            'y0': params.y0,
            'x0': params.x0,
            'a': params.a,
            'sigma': params.sigma,
            'm': params.m,
            'c': params.c,

            #Stars and gas kinematics parameters
            'left_wave_ppxf_kin': params.wave1_kin,
            'right_wave_ppxf_kin': params.wave2_kin,
            'stellar_library_kin': params.stellar_library_kin,
            'constant_resolution_lambda': params.constant_resolution_lambda,
            'resolution_kin_muse': params.resolution_kin_muse,
            'ppxf_resolution': params.resolution_kin,
            'constant_resolution_r': params.constant_resolution_r,
            'ppxf_resolution_r': params.resolution_kin_r,
            'sigma_guess_kin': params.sigma_guess_kin,
            'redshift_guess_kin': params.redshift_guess_kin,
            'additive_degree_kin': params.additive_degree_kin,
            'multiplicative_degree_kin': params.multiplicative_degree_kin,
            'gas_kin': params.gas_kin,
            'no_gas_kin': params.no_gas_kin,
            'kin_best_noise': params.kin_best_noise,
            'with_errors_kin': params.with_errors_kin,
            'kin_moments': params.kin_moments,
            'ppxf_kin_noise': params.ppxf_kin_noise,
            'ppxf_kin_preloaded_lib': params.ppxf_kin_preloaded_lib,
            'ppxf_kin_custom_lib': params.ppxf_kin_custom_lib,
            'ppxf_kin_lib_folder': params.ppxf_kin_lib_folder,
            'ppxf_kin_custom_temp_suffix': params.ppxf_kin_custom_temp_suffix,
            'ppxf_kin_generic_lib': params.ppxf_kin_generic_lib,
            'ppxf_kin_generic_lib_folder': params.ppxf_kin_generic_lib_folder,
            'ppxf_kin_FWHM_tem_generic': params.ppxf_kin_FWHM_tem_generic,
            'ppxf_kin_fixed_kin': params.ppxf_kin_fixed_kin,
            'ppxf_kin_tie_balmer': params.ppxf_kin_tie_balmer,
            'ppxf_kin_dust_stars': params.ppxf_kin_dust_stars,
            'ppxf_kin_dust_gas': params.ppxf_kin_dust_gas,
            'ppxf_kin_two_stellar_components': params.ppxf_kin_two_stellar_components,
            'ppxf_kin_age_model1': params.ppxf_kin_age_model1,
            'ppxf_kin_met_model1': params.ppxf_kin_met_model1,
            'ppxf_kin_age_model2': params.ppxf_kin_age_model2,
            'ppxf_kin_met_model2': params.ppxf_kin_met_model2,
            'ppxf_kin_vel_model1': params.ppxf_kin_vel_model1,
            'ppxf_kin_sigma_model1': params.ppxf_kin_sigma_model1,
            'ppxf_kin_vel_model2': params.ppxf_kin_vel_model2,
            'ppxf_kin_sigma_model2': params.ppxf_kin_sigma_model2,
            'ppxf_kin_mask_emission': params.ppxf_kin_mask_emission,
            'ppxf_kin_have_user_mask': params.ppxf_kin_have_user_mask,
            'ppxf_kin_mask_ranges_str': params.ppxf_kin_mask_ranges_str,

            'ppxf_kin_mc_sim': params.ppxf_kin_mc_sim,
            'ppxf_kin_save_spectra': params.ppxf_kin_save_spectra,

            #Stellar populations and SFH parameters
            'left_wave_ppxf_pop': params.wave1_pop,
            'right_wave_ppxf_pop': params.wave2_pop,
            'resolution_ppxf_pop': params.res_pop,
            'sigma_guess_pop': params.sigma_guess_pop,
            'ppxf_z_pop': params.z_pop,
            'gas_pop': params.pop_with_gas,
            'ppxf_pop_tie_balmer': params.ppxf_pop_tie_balmer,
            'ppxf_pop_dust_stars': params.ppxf_pop_dust_stars,
            'ppxf_pop_dust_gas': params.ppxf_pop_dust_gas,
            'ppxf_pop_noise': params.ppxf_pop_noise,
            'ppxf_min_age': params.ppxf_min_age,
            'ppxf_max_age': params.ppxf_max_age,
            'ppxf_min_met': params.ppxf_min_met,
            'ppxf_max_met': params.ppxf_max_met,
            'no_gas_pop': params.pop_without_gas,
            'regul_err': params.regul_err,
            'additive_degree': params.additive_degree,
            'multiplicative_degree': params.multiplicative_degree,
            'stellar_library': params.stellar_library,
            'ppxf_err_pop': params.with_errors,
            'ppxf_pop_preloaded_lib': params.ppxf_pop_preloaded_lib,
            'ppxf_pop_custom_lib': params.ppxf_pop_custom_lib,
            'ppxf_pop_lib_folder': params.ppxf_pop_lib_folder,
            'ppxf_pop_custom_npz': params.ppxf_pop_custom_npz,
            'ppxf_pop_npz_file': params.ppxf_pop_npz_file,
            'ppxf_pop_mask': params.ppxf_pop_mask,
            'ppxf_custom_temp_suffix': params.ppxf_custom_temp_suffix,
            'ppxf_best_param': params.ppxf_best_param,
            'ppxf_best_noise_estimate': params.ppxf_best_noise_estimate,
            'ppxf_frac_chi': params.ppxf_frac_chi,
            'ppxf_pop_convolve': params.ppxf_pop_convolve,
            'ppxf_pop_want_to_mask': params.ppxf_pop_want_to_mask,
            'ppxf_pop_mask_ranges_str': params.ppxf_pop_mask_ranges_str,
            'ppxf_pop_error_nsim': params.ppxf_pop_error_nsim,
            'ppxf_pop_lg_age': params.ppxf_pop_lg_age,
            'ppxf_pop_lg_met': params.ppxf_pop_lg_met,
            'stellar_parameters_lick_ppxf': params.stellar_parameters_lick_ppxf,
            'ssp_model_ppxf': params.ssp_model_ppxf,
            'interp_model_ppxf': params.interp_model_ppxf,

            # Long-slit (2D) extraction parameters
            'file_path': params.file_path_spec_extr,
            'trace_y_range': params.trace_y_range_str,
            'poly_degree': params.poly_degree_str,
            'extract_y_range': params.extract_y_range_str,
            'snr': params.snr_threshold_str,
            'pix_scale': params.pixel_scale_str,

            #cube extraction parameters
            'ifs_run_id': params.ifs_run_id,
            'ifs_input': params.ifs_input,
            'ifs_redshift': params.ifs_redshift,
            'ifs_routine_read': params.ifs_routine_read_default,
            'ifs_origin': params.ifs_origin,
            'ifs_lmin_tot': params.ifs_lmin_tot,
            'ifs_lmax_tot': params.ifs_lmax_tot,
            'ifs_min_snr_mask': params.ifs_min_snr_mask,
            'ifs_mask': params.ifs_mask,
            'ifs_target_snr': params.ifs_target_snr,
            'ifs_preloaded_routine': params.ifs_preloaded_routine,
            'ifs_user_routine': params.ifs_user_routine,
            'ifs_user_routine_file': params.ifs_user_routine_file,
            'ifs_manual_bin': params.ifs_manual_bin,
            'ifs_voronoi': params.ifs_voronoi,
            'ifs_existing_bin': params.ifs_existing_bin,
            'ifs_existing_bin_folder': params.ifs_existing_bin_folder,
        }
    }
    with open(filename, 'w') as file:
        json.dump(data, file, indent=4)




# function to load the JSON configuration file
def load_settings(filename, params):
    try:
        with open(filename, 'r') as file:
            data = json.load(file)

        keys = data.get('keys', [])
        events = data.get('events', [])
        values = data.get('values', {})
        params_data = data.get('params', {})


        # Now I need to manually update all the parameters of the dataclass
        updated_params = replace(
            params,
                cropping_spectrum=params_data.get('cropping', params.cropping_spectrum),
                cropping_low_wave=params_data.get('cropping_low_wave', params.cropping_low_wave),
                cropping_high_wave=params_data.get('cropping_high_wave', params.cropping_high_wave),
                sigma_clipping=params_data.get('sigma_clip', params.sigma_clipping),
                clip_factor=params_data.get('clip_factor', params.clip_factor),
                sigma_clip_resolution=params_data.get('res_spec_for_sigma_clip', params.sigma_clip_resolution),
                sigma_clip_single_vel=params_data.get('single_vel_clip', params.sigma_clip_single_vel),
                sigma_clip_single_value=params_data.get('clip_to_vel', params.sigma_clip_single_value),
                sigma_clip_have_file=params_data.get('file_for_clip', params.sigma_clip_have_file),
                sigma_clip_sigma_file=params_data.get('sigma_clip_file', params.sigma_clip_sigma_file),
                wavelet_cleaning=params_data.get('wavelet_cleaning', params.wavelet_cleaning),
                sigma_wavelets=params_data.get('sigma_wavelets', params.sigma_wavelets),
                wavelets_layers=params_data.get('wavelets_layers', params.wavelets_layers),
                dop_cor=params_data.get('dopcor', params.dop_cor),

                dop_cor_single_shot_vel=params_data.get('dopcor_value', params.dop_cor_single_shot_vel),
                dop_cor_have_vel=params_data.get('dop_cor_have_vel', params.dop_cor_have_vel),
                dop_cor_have_z=params_data.get('dop_cor_have_z', params.dop_cor_have_z),
                dop_cor_have_file=params_data.get('file_for_dopcor', params.dop_cor_have_file),
                dop_cor_file=params_data.get('dopcor_file', params.dop_cor_file),
                dop_cor_single_shot=params_data.get('dopcor_single_value', params.dop_cor_single_shot),

                helio_corr=params_data.get('helio_corr', params.helio_corr),
                helio_have_file=params_data.get('file_for_helio', params.helio_have_file),
                helio_file=params_data.get('helio_file', params.helio_file),
                helio_single_shot=params_data.get('helio_single_value', params.helio_single_shot),
                helio_single_shot_location=params_data.get('helio_location', params.helio_single_shot_location),
                helio_single_shot_date=params_data.get('helio_date', params.helio_single_shot_date),
                ra_obj=params_data.get('helio_ra', params.ra_obj),
                dec_obj=params_data.get('helio_dec', params.dec_obj),

                rebinning=params_data.get('rebin', params.rebinning),
                rebinning_linear=params_data.get('rebin_pix_lin', params.rebinning_linear),
                rebin_step_pix=params_data.get('rebin_step_pix', params.rebin_step_pix),
                rebinning_log=params_data.get('rebin_sigma_lin', params.rebinning_log),
                rebin_step_sigma=params_data.get('rebin_step_sigma', params.rebin_step_sigma),

                degrade=params_data.get('degrade_resolution', params.degrade),
                is_initial_res_r=params_data.get('is_initial_res_r', params.is_initial_res_r),
                initial_res_r=params_data.get('degrade_from_r', params.initial_res_r),
                res_degrade_to_r=params_data.get('res_degrade_to_r', params.res_degrade_to_r),
                final_res_r=params_data.get('degrade_to_r', params.final_res_r),
                res_degrade_to_fwhm=params_data.get('res_degrade_to_fwhm', params.res_degrade_to_fwhm),
                final_res_r_to_fwhm=params_data.get('final_res_r_to_fwhm', params.final_res_r_to_fwhm),
                is_initial_res_fwhm=params_data.get('is_initial_res_fwhm', params.is_initial_res_fwhm),
                initial_res_fwhm=params_data.get('degrade_from_l', params.initial_res_fwhm),
                final_res_fwhm=params_data.get('degrade_to_l', params.final_res_fwhm),
                res_degrade_muse= params_data.get('res_degrade_muse', params.res_degrade_muse),
                res_degrade_muse_value= params_data.get('res_degrade_muse_value', params.res_degrade_muse_value),

                normalize_wave=params_data.get('norm_spec', params.normalize_wave),
                norm_wave=params_data.get('norm_wave', params.norm_wave),

                continuum_sub=params_data.get('cont_sub', params.continuum_sub),
                cont_model_filtering=params_data.get('cont_model_filtering', params.cont_model_filtering),
                cont_model_poly=params_data.get('cont_model_poly', params.cont_model_poly),
                cont_math_operation=params_data.get('markers_cont_operations', params.cont_math_operation),
                cont_want_to_mask=params_data.get('cont_want_to_mask', params.cont_want_to_mask),
                cont_mask_ranges_str=params_data.get('cont_mask_ranges', params.cont_mask_ranges_str),
                cont_poly_degree=params_data.get('cont_poly_degree', params.cont_poly_degree),

                sigma_broad=params_data.get('broadening_spec', params.sigma_broad),
                sigma_to_add=params_data.get('sigma_to_add', params.sigma_to_add),
                add_noise=params_data.get('add_noise', params.add_noise),
                noise_to_add=params_data.get('noise_to_add', params.noise_to_add),

                filter_denoise=params_data.get('filter_denoise', params.filter_denoise),
                moving_average=params_data.get('moving_average', params.moving_average),
                box_moving_avg=params_data.get('box_moving_avg', params.box_moving_avg),
                box_moving_avg_size=params_data.get('box_moving_avg_size', params.box_moving_avg_size),
                gauss_moving_avg=params_data.get('gauss_moving_avg', params.gauss_moving_avg),
                gauss_moving_avg_kernel=params_data.get('gauss_moving_avg_kernel', params.gauss_moving_avg_kernel),
                low_pass_filter=params_data.get('low_pass_filter', params.low_pass_filter),
                lowpass_cut_off=params_data.get('lowpass_cut_off', params.lowpass_cut_off),
                lowpass_order=params_data.get('lowpass_order', params.lowpass_order),
                bandpass_filter=params_data.get('bandpass_filter', params.bandpass_filter),
                bandpass_lower_cut_off=params_data.get('bandpass_lower_cut_off', params.bandpass_lower_cut_off),
                bandpass_upper_cut_off=params_data.get('bandpass_upper_cut_off', params.bandpass_upper_cut_off),
                bandpass_order=params_data.get('bandpass_order', params.bandpass_order),

                average_all=params_data.get('avg_all', params.average_all),
                norm_and_average=params_data.get('norm_avg_all', params.norm_and_average),
                do_nothing=params_data.get('none', params.do_nothing),
                sum_all=params_data.get('sum_all', params.sum_all),
                normalize_and_sum_all=params_data.get('norm_sum_all', params.normalize_and_sum_all),
                use_for_spec_an=params_data.get('use_for_spec_an', params.use_for_spec_an),
                subtract_normalized_avg=params_data.get('subtract_norm_avg', params.subtract_normalized_avg),
                subtract_normalized_spec=params_data.get('subtract_norm_spec', params.subtract_normalized_spec),
                spectra_to_subtract=params_data.get('spec_to_subtract', params.spectra_to_subtract),
                add_pedestal=params_data.get('add_pedestal', params.add_pedestal),
                pedestal_to_add=params_data.get('pedestal_to_add', params.pedestal_to_add),
                multiply=params_data.get('multiply', params.multiply),
                multiply_factor=params_data.get('multiply_factor', params.multiply_factor),
                derivatives=params_data.get('derivatives', params.derivatives),
                reorder_op=params_data.get('reorder_op', params.reorder_op),
                current_order=params_data.get('current_order', params.current_order),
                active_operations=params_data.get('active_operations', params.active_operations),
                reordered_operations=params_data.get('reordered_operations', params.reordered_operations),
                save_intermediate_spectra=params_data.get('save_intermediate_spectra', params.save_intermediate_spectra),
                save_final_spectra=params_data.get('save_final_spectra', params.save_final_spectra),
                not_save_spectra=params_data.get('not_save_spectra', params.not_save_spectra),

                sigma_coeff=params_data.get('sigma_coeff', params.sigma_coeff),
                sigma_corr=params_data.get('sigma_corr', params.sigma_corr),

                # Blackbody parameters
                wave1_bb=params_data.get('left_wave_bb', params.wave1_bb),
                wave2_bb=params_data.get('right_wave_bb', params.wave2_bb),
                t_guess=params_data.get('t_guess_bb', params.t_guess),

                # Cross-correlation parameters
                template_crosscorr=params_data.get('xcorr_template', params.template_crosscorr),
                lambda_units_template_crosscorr_nm=params_data.get('xcorr_template_wave_nm', params.lambda_units_template_crosscorr_nm),
                lambda_units_template_crosscorr_a=params_data.get('xcorr_template_wave_a', params.lambda_units_template_crosscorr_a),
                lambda_units_template_crosscorr_mu=params_data.get('xcorr_template_wave_mu', params.lambda_units_template_crosscorr_mu),
                smooth_template_crosscorr=params_data.get('xcorr_smooth_template', params.smooth_template_crosscorr),
                smooth_value_crosscorr=params_data.get('xcorr_smooth_template_value', params.smooth_value_crosscorr),
                low_wave_corr=params_data.get('xcorr_left_lambda', params.low_wave_corr),
                high_wave_corr=params_data.get('xcorr_right_lambda', params.high_wave_corr),
                is_vel_xcorr=params_data.get('is_vel_xcorr', params.is_vel_xcorr),
                is_z_xcorr=params_data.get('is_z_xcorr', params.is_z_xcorr),
                low_vel_corr=params_data.get('xcorr_low_vel', params.low_vel_corr),
                high_vel_corr=params_data.get('xcorr_high_vel', params.high_vel_corr),
                low_z_corr=params_data.get('low_z_corr', params.low_z_corr),
                high_z_corr=params_data.get('high_z_corr', params.high_z_corr),
                xcorr_limit_wave_range=params_data.get('xcorr_limit_wave_range', params.xcorr_limit_wave_range),
                xcorr_vel_step=params_data.get('xcorr_vel_step', params.xcorr_vel_step),
                xcorr_z_step=params_data.get('xcorr_z_step', params.xcorr_z_step),

                # Velocity dispersion parameters
                template_sigma=params_data.get('template_sigma', params.template_sigma),
                lambda_units_template_sigma_nm=params_data.get('lambda_units_template_sigma_nm', params.lambda_units_template_sigma_nm),
                lambda_units_template_sigma_a=params_data.get('lambda_units_template_sigma_a', params.lambda_units_template_sigma_a),
                lambda_units_template_sigma_mu=params_data.get('lambda_units_template_sigma_mu', params.lambda_units_template_sigma_mu),
                band_cat=params_data.get('band_cat', params.band_cat),
                band_halpha=params_data.get('band_halpha', params.band_halpha),
                band_nad=params_data.get('band_nad', params.band_nad),
                band_h=params_data.get('band_h', params.band_h),
                band_k=params_data.get('band_k', params.band_k),
                resolution_spec=params_data.get('resolution_spec', params.resolution_spec),
                resolution_template=params_data.get('resolution_template', params.resolution_template),
                band_custom=params_data.get('band_custom', params.band_custom),
                low_wave_sigma=params_data.get('low_wave_sigma', params.low_wave_sigma),
                high_wave_sigma=params_data.get('high_wave_sigma', params.high_wave_sigma),
                low_wave_cont=params_data.get('low_wave_cont', params.low_wave_cont),
                high_wave_cont=params_data.get('high_wave_cont', params.high_wave_cont),

                band_sigma = np.array([params.low_wave_sigma, params.high_wave_sigma]),
                cont_sigma = np.array([params.low_wave_cont, params.high_wave_cont]),

                # Line-strength parameters
                have_index_file=params_data.get('ew_idx_file', params.have_index_file),
                index_file=params_data.get('idx_file', params.index_file),
                single_index=params_data.get('single_index', params.single_index),
                idx_left_blue=params_data.get('left_wave_blue_cont', params.idx_left_blue),
                idx_right_blue=params_data.get('right_wave_blue_cont', params.idx_right_blue),
                idx_left_red=params_data.get('left_wave_red_cont', params.idx_left_red),
                idx_right_red=params_data.get('right_wave_red_cont', params.idx_right_red),
                idx_left_line=params_data.get('left_line', params.idx_left_line),
                idx_right_line=params_data.get('right_line', params.idx_right_line),
                lick_ew=params_data.get('ew_lick', params.lick_ew),
                lick_constant_fwhm=params_data.get('lick_constant_fwhm', params.lick_constant_fwhm),
                spec_lick_res_fwhm=params_data.get('spec_lick_res_fwhm', params.spec_lick_res_fwhm),
                lick_constant_r=params_data.get('lick_constant_r', params.lick_constant_r),
                spec_lick_res_r=params_data.get('spec_lick_res_r', params.spec_lick_res_r),
                lick_correct_emission=params_data.get('lick_correct_emission', params.lick_correct_emission),
                z_guess_lick_emission=params_data.get('z_guess_lick_emission', params.z_guess_lick_emission),
                dop_correction_lick=params_data.get('dop_correction_lick', params.dop_correction_lick),
                correct_ew_sigma=params_data.get('correct_ew_sigma', params.correct_ew_sigma),
                radio_lick_sigma_auto=params_data.get('radio_lick_sigma_auto', params.radio_lick_sigma_auto),
                radio_lick_sigma_single=params_data.get('radio_lick_sigma_single', params.radio_lick_sigma_single),
                sigma_single_lick=params_data.get('sigma_single_lick', params.sigma_single_lick),
                radio_lick_sigma_list=params_data.get('radio_lick_sigma_list', params.radio_lick_sigma_list),
                sigma_lick_file=params_data.get('sigma_lick_file', params.sigma_lick_file),
                stellar_parameters_lick=params_data.get('stellar_parameters_lick', params.stellar_parameters_lick),
                ssp_model=params_data.get('ssp_model', params.ssp_model),
                interp_model=params_data.get('interp_model', params.interp_model),
                have_index_file_corr=params_data.get('ew_corr_idx_file', params.have_index_file_corr),
                index_file_corr=params_data.get('idx_corr_file', params.index_file_corr),
                single_index_corr=params_data.get('ew_corr_single_idx', params.single_index_corr),

                # Stellar spectra coefficients
                stellar_spectra_coeff_file=params_data.get('sigma_coeff_sample_list', params.stellar_spectra_coeff_file),
                lambda_units_coeff_nm=params_data.get('sigma_coeff_sample_list_wave_nm', params.lambda_units_coeff_nm),
                lambda_units_coeff_a=params_data.get('sigma_coeff_sample_list_wave_a', params.lambda_units_coeff_a),
                lambda_units_coeff_mu=params_data.get('sigma_coeff_sample_list_wave_mu', params.lambda_units_coeff_mu),
                smooth_stellar_sample=params_data.get('sigma_coeff_sample_smooth', params.smooth_stellar_sample),
                smooth_value_sample=params_data.get('sigma_coeff_sample_smooth_sigma', params.smooth_value_sample),
                sigma_vel_file=params_data.get('sigma_file', params.sigma_vel_file),
                ew_list_file=params_data.get('ew_file_to_correct', params.ew_list_file),
                sigma_coeff_file=params_data.get('coeff_sigma_file', params.sigma_coeff_file),

                # Line(s) fitting parameters
                cat_band_fit=params_data.get('cat_fit', params.cat_band_fit),
                usr_fit_line=params_data.get('line_fit_single', params.usr_fit_line),
                emission_line=params_data.get('emission_line', params.emission_line),
                low_wave_fit=params_data.get('left_wave_fitting', params.low_wave_fit),
                high_wave_fit=params_data.get('right_wave_fitting', params.high_wave_fit),
                y0=params_data.get('y0', params.y0),
                x0=params_data.get('x0', params.x0),
                a=params_data.get('a', params.a),
                sigma=params_data.get('sigma', params.sigma),
                m=params_data.get('m', params.m),
                c=params_data.get('c', params.c),

                # Stars and gas kinematics parameters
                wave1_kin=params_data.get('left_wave_ppxf_kin', params.wave1_kin),
                wave2_kin=params_data.get('right_wave_ppxf_kin', params.wave2_kin),
                stellar_library_kin=params_data.get('stellar_library_kin', params.stellar_library_kin),
                constant_resolution_lambda=params_data.get('constant_resolution_lambda', params.constant_resolution_lambda),
                resolution_kin=params_data.get('ppxf_resolution', params.resolution_kin),
                resolution_kin_muse = params_data.get('resolution_kin_muse', params.resolution_kin_muse),
                constant_resolution_r=params_data.get('constant_resolution_r', params.constant_resolution_r),
                resolution_kin_r=params_data.get('ppxf_resolution_r', params.resolution_kin_r),
                sigma_guess_kin=params_data.get('sigma_guess_kin', params.sigma_guess_kin),
                redshift_guess_kin=params_data.get('redshift_guess_kin', params.redshift_guess_kin),
                additive_degree_kin=params_data.get('additive_degree_kin', params.additive_degree_kin),
                multiplicative_degree_kin=params_data.get('multiplicative_degree_kin', params.multiplicative_degree_kin),
                gas_kin=params_data.get('gas_kin', params.gas_kin),
                no_gas_kin=params_data.get('no_gas_kin', params.no_gas_kin),
                kin_best_noise=params_data.get('kin_best_noise', params.kin_best_noise),
                with_errors_kin=params_data.get('with_errors_kin', params.with_errors_kin),
                kin_moments=params_data.get('kin_moments', params.kin_moments),
                ppxf_kin_noise=params_data.get('ppxf_kin_noise', params.ppxf_kin_noise),
                ppxf_kin_preloaded_lib=params_data.get('ppxf_kin_preloaded_lib', params.ppxf_kin_preloaded_lib),
                ppxf_kin_custom_lib=params_data.get('ppxf_kin_custom_lib', params.ppxf_kin_custom_lib),
                ppxf_kin_lib_folder=params_data.get('ppxf_kin_lib_folder', params.ppxf_kin_lib_folder),
                ppxf_kin_custom_temp_suffix=params_data.get('ppxf_kin_custom_temp_suffix', params.ppxf_kin_custom_temp_suffix),
                ppxf_kin_generic_lib=params_data.get('ppxf_kin_generic_lib', params.ppxf_kin_generic_lib),
                ppxf_kin_generic_lib_folder=params_data.get('ppxf_kin_generic_lib_folder', params.ppxf_kin_generic_lib_folder),
                ppxf_kin_FWHM_tem_generic=params_data.get('ppxf_kin_FWHM_tem_generic', params.ppxf_kin_FWHM_tem_generic),
                ppxf_kin_fixed_kin=params_data.get('ppxf_kin_fixed_kin', params.ppxf_kin_fixed_kin),

                ppxf_kin_tie_balmer=params_data.get('ppxf_kin_tie_balmer', params.ppxf_kin_tie_balmer),
                ppxf_kin_dust_stars=params_data.get('ppxf_kin_dust_stars', params.ppxf_kin_dust_stars),
                ppxf_kin_dust_gas=params_data.get('ppxf_kin_dust_gas', params.ppxf_kin_dust_gas),
                ppxf_kin_two_stellar_components=params_data.get('ppxf_kin_two_stellar_components', params.ppxf_kin_two_stellar_components),
                ppxf_kin_age_model1=params_data.get('ppxf_kin_age_model1', params.ppxf_kin_age_model1),
                ppxf_kin_met_model1=params_data.get('ppxf_kin_met_model1', params.ppxf_kin_met_model1),
                ppxf_kin_age_model2=params_data.get('ppxf_kin_age_model2', params.ppxf_kin_age_model2),
                ppxf_kin_met_model2=params_data.get('ppxf_kin_met_model2', params.ppxf_kin_met_model2),
                ppxf_kin_vel_model1=params_data.get('ppxf_kin_vel_model1', params.ppxf_kin_vel_model1),
                ppxf_kin_sigma_model1=params_data.get('ppxf_kin_sigma_model1', params.ppxf_kin_sigma_model1),
                ppxf_kin_vel_model2=params_data.get('ppxf_kin_vel_model2', params.ppxf_kin_vel_model2),
                ppxf_kin_sigma_model2=params_data.get('ppxf_kin_sigma_model2', params.ppxf_kin_sigma_model2),
                ppxf_kin_mask_emission=params_data.get('ppxf_kin_mask_emission', params.ppxf_kin_mask_emission),
                ppxf_kin_have_user_mask=params_data.get('ppxf_kin_have_user_mask', params.ppxf_kin_have_user_mask),
                ppxf_kin_mask_ranges_str=params_data.get('ppxf_kin_mask_ranges_str', params.ppxf_kin_mask_ranges_str),
                ppxf_kin_mc_sim=params_data.get('ppxf_kin_mc_sim', params.ppxf_kin_mc_sim),
                ppxf_kin_save_spectra=params_data.get('ppxf_kin_save_spectra', params.ppxf_kin_save_spectra),

                # Stellar populations and SFH parameters
                wave1_pop=params_data.get('left_wave_ppxf_pop', params.wave1_pop),
                wave2_pop=params_data.get('right_wave_ppxf_pop', params.wave2_pop),
                res_pop=params_data.get('resolution_ppxf_pop', params.res_pop),
                sigma_guess_pop=params_data.get('sigma_guess_pop', params.sigma_guess_pop),
                z_pop=params_data.get('ppxf_z_pop', params.z_pop),
                pop_with_gas=params_data.get('gas_pop', params.pop_with_gas),
                ppxf_pop_tie_balmer=params_data.get('ppxf_pop_tie_balmer', params.ppxf_pop_tie_balmer),
                ppxf_pop_dust_stars=params_data.get('ppxf_pop_dust_stars', params.ppxf_pop_dust_stars),
                ppxf_pop_dust_gas=params_data.get('ppxf_pop_dust_gas', params.ppxf_pop_dust_gas),
                ppxf_pop_noise=params_data.get('ppxf_pop_noise', params.ppxf_pop_noise),
                ppxf_min_age=params_data.get('ppxf_min_age', params.ppxf_min_age),
                ppxf_max_age=params_data.get('ppxf_max_age', params.ppxf_max_age),
                ppxf_min_met=params_data.get('ppxf_min_met', params.ppxf_min_met),
                ppxf_max_met=params_data.get('ppxf_max_met', params.ppxf_max_met),
                pop_without_gas=params_data.get('no_gas_pop', params.pop_without_gas),
                regul_err=params_data.get('regul_err', params.regul_err),
                additive_degree=params_data.get('additive_degree', params.additive_degree),
                multiplicative_degree=params_data.get('multiplicative_degree', params.multiplicative_degree),
                stellar_library=params_data.get('stellar_library', params.stellar_library),
                with_errors=params_data.get('ppxf_err_pop', params.with_errors),
                ppxf_pop_preloaded_lib=params_data.get('ppxf_pop_preloaded_lib', params.ppxf_pop_preloaded_lib),
                ppxf_pop_custom_lib=params_data.get('ppxf_pop_custom_lib', params.ppxf_pop_custom_lib),
                ppxf_pop_lib_folder=params_data.get('ppxf_pop_lib_folder', params.ppxf_pop_lib_folder),
                ppxf_pop_custom_npz=params_data.get('ppxf_pop_custom_npz', params.ppxf_pop_custom_npz),
                ppxf_pop_npz_file=params_data.get('ppxf_pop_npz_file', params.ppxf_pop_npz_file),
                ppxf_pop_mask=params_data.get('ppxf_pop_mask', params.ppxf_pop_mask),
                ppxf_custom_temp_suffix=params_data.get('ppxf_custom_temp_suffix', params.ppxf_custom_temp_suffix),
                ppxf_best_param=params_data.get('ppxf_best_param', params.ppxf_best_param),
                ppxf_best_noise_estimate=params_data.get('ppxf_best_noise_estimate', params.ppxf_best_noise_estimate),
                ppxf_frac_chi=params_data.get('ppxf_frac_chi', params.ppxf_frac_chi),
                ppxf_pop_convolve=params_data.get('ppxf_pop_convolve', params.ppxf_pop_convolve),
                ppxf_pop_want_to_mask=params_data.get('ppxf_pop_want_to_mask', params.ppxf_pop_want_to_mask),
                ppxf_pop_mask_ranges_str=params_data.get('ppxf_pop_mask_ranges_str', params.ppxf_pop_mask_ranges_str),
                ppxf_pop_error_nsim=params_data.get('ppxf_pop_error_nsim', params.ppxf_pop_error_nsim),
                ppxf_pop_lg_age=params_data.get('ppxf_pop_lg_age', params.ppxf_pop_lg_age),
                ppxf_pop_lg_met=params_data.get('ppxf_pop_lg_met', params.ppxf_pop_lg_met),
                stellar_parameters_lick_ppxf=params_data.get('stellar_parameters_lick_ppxf', params.stellar_parameters_lick_ppxf),
                ssp_model_ppxf=params_data.get('params.ssp_model_ppxf', params.ssp_model_ppxf),
                interp_model_ppxf=params_data.get('interp_model_ppxf', params.interp_model_ppxf),

                # Long-slit (2D) extraction parameters
                file_path_spec_extr=params_data.get('file_path', params.file_path_spec_extr),
                trace_y_range_str=params_data.get('trace_y_range', params.trace_y_range_str),
                poly_degree_str=params_data.get('poly_degree', params.poly_degree_str),
                extract_y_range_str=params_data.get('extract_y_range', params.extract_y_range_str),
                snr_threshold_str=params_data.get('snr', params.snr_threshold_str),
                pixel_scale_str=params_data.get('pix_scale', params.pixel_scale_str),

                # Cube extraction parameters
                ifs_run_id=params_data.get('ifs_run_id', params.ifs_run_id),
                ifs_input=params_data.get('ifs_input', params.ifs_input),
                ifs_redshift=params_data.get('ifs_redshift', params.ifs_redshift),
                ifs_routine_read_default=params_data.get('ifs_routine_read', params.ifs_routine_read_default),
                ifs_origin=params_data.get('ifs_origin', params.ifs_origin),
                ifs_lmin_tot=params_data.get('ifs_lmin_tot', params.ifs_lmin_tot),
                ifs_lmax_tot=params_data.get('ifs_lmax_tot', params.ifs_lmax_tot),
                ifs_min_snr_mask=params_data.get('ifs_min_snr_mask', params.ifs_min_snr_mask),
                ifs_mask=params_data.get('ifs_mask', params.ifs_mask),
                ifs_target_snr=params_data.get('ifs_target_snr', params.ifs_target_snr),
                ifs_preloaded_routine=params_data.get('ifs_preloaded_routine', params.ifs_preloaded_routine),
                ifs_user_routine=params_data.get('ifs_user_routine', params.ifs_user_routine),
                ifs_user_routine_file=params_data.get('ifs_user_routine_file', params.ifs_user_routine_file),
                ifs_manual_bin=params_data.get('ifs_manual_bin', params.ifs_manual_bin),
                ifs_voronoi=params_data.get('ifs_voronoi', params.ifs_voronoi),
                ifs_existing_bin = params_data.get('ifs_existing_bin', params.ifs_existing_bin),
                ifs_existing_bin_folder = params_data.get('ifs_existing_bin_folder', params.ifs_existing_bin_folder),


        )

        return keys, events, values, updated_params

    except FileNotFoundError:
        print("No settings file found. Using default values.")
        return [], [], {}, params

