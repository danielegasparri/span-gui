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

# Functions to define the parameter windows of the Spectral analysis frame.
# They return the updated params modified by the user

try: #try local import if executed as script
    #GUI import
    from FreeSimpleGUI_local import FreeSimpleGUI as sg

    #SPAN functions import
    from span_functions import system_span as stm
    from span_functions import spec_manipul as spman
    from span_functions import linestrength as ls
    from span_functions import spec_analysis as span
    from span_modules import layouts
    from span_modules import misc
    from params import SpectraParams

except ModuleNotFoundError: #local import if executed as package
    #GUI import
    from span.FreeSimpleGUI_local import FreeSimpleGUI as sg

    #SPAN functions import
    from span.span_functions import system_span as stm
    from span.span_functions import utilities as uti
    from span.span_functions import spec_manipul as spman
    from span.span_functions import linestrength as ls
    from span.span_functions import spec_analysis as span
    from . import layouts
    from . import misc
    from .params import SpectraParams

#python imports
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.backends.backend_tkagg
import time
import os
import glob
from dataclasses import replace


CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(CURRENT_DIR)

fontsize = sg.set_options(font=("Helvetica", 11)) # defaul fontsize



def blackbody_parameters(params: SpectraParams) -> SpectraParams:
    """Handles blackbody fitting parameter input via GUI."""

    sg.theme('LightBlue1')

    # Extract values from params
    wave1_bb = params.wave1_bb
    wave2_bb = params.wave2_bb
    t_guess = params.t_guess

    layout, scale_win, fontsize, default_size = misc.get_layout()

    # Define GUI layout
    bb_layout = [
        [sg.Text('Wave interval (A):'),
         sg.InputText(wave1_bb, size=(6,1), key='left_wave_bb'), sg.Text('-'),
         sg.InputText(wave2_bb, size=(6,1), key='right_wave_bb')],
        [sg.Text('Initial Temperature guess'),
         sg.InputText(t_guess, size=(8,1), key='t_guess_bb')],
        [sg.Push(), sg.Button('Confirm', button_color=('white', 'black'), size=(18,1))]
    ]

    print('*** Blackbody fitting parameters window open. The main panel will be inactive until you close the window ***')

    bb_window = sg.Window('Blackbody fitting parameters', bb_layout)

    while True:
        bb_event, bb_values = bb_window.read()

        if bb_event == sg.WIN_CLOSED:
            break

        try:
            wave1_bb = float(bb_values['left_wave_bb'])
            wave2_bb = float(bb_values['right_wave_bb'])
            t_guess = float(bb_values['t_guess_bb'])

            if t_guess <= 0:
                sg.popup('No blackbody has a negative temperature!')
                continue
            if t_guess > 1e7:
                sg.popup('No stellar blackbody has a temperature greater than 10 million degrees!')
                continue
            if wave1_bb >= wave2_bb:
                sg.popup('The first wavelength cannot be greater than the second!')
                continue
            if wave2_bb - wave1_bb <= 5.:
                sg.popup('The wavelength interval is too small to perform a good fit. Enlarge it!')
                continue

        except ValueError:
            sg.popup('Invalid input parameters!')
            continue

        if bb_event == 'Confirm':
            print('Blackbody parameters confirmed. The main panel is now active again.\n')
            break

    bb_window.close()

    # Update params only with modified values
    return replace(params, wave1_bb=wave1_bb, wave2_bb=wave2_bb, t_guess=t_guess)



def crosscorr_parameters(params: SpectraParams) -> SpectraParams:
    """Handles cross-correlation parameter input via GUI."""

    # Extract parameters from params
    template_crosscorr = params.template_crosscorr
    lambda_units_template_crosscorr = 'nm' if params.lambda_units_template_crosscorr_nm else \
                                      'a' if params.lambda_units_template_crosscorr_a else \
                                      'mu' if params.lambda_units_template_crosscorr_mu else 'nm'
    smooth_template_crosscorr = params.smooth_template_crosscorr
    smooth_value_crosscorr = params.smooth_value_crosscorr
    low_wave_corr = params.low_wave_corr
    high_wave_corr = params.high_wave_corr
    wave_interval_corr = params.wave_interval_corr
    is_vel_xcorr = params.is_vel_xcorr
    low_vel_corr = params.low_vel_corr
    high_vel_corr = params.high_vel_corr
    is_z_xcorr = params.is_z_xcorr
    low_z_corr = params.low_z_corr
    high_z_corr = params.high_z_corr
    xcorr_limit_wave_range = params.xcorr_limit_wave_range
    xcorr_vel_step = params.xcorr_vel_step
    xcorr_z_step = params.xcorr_z_step

    layout, scale_win, fontsize, default_size = misc.get_layout()
    sg.theme('LightBlue1')

    xcorr_layout = [
        [sg.Text('Select a template:', font=('Helvetica', 10, 'bold')),
         sg.InputText(template_crosscorr, size=(32,1), key='xcorr_template'),
         sg.FileBrowse(tooltip='Load a template')],
        [sg.Text('Template wavelength is in:'),
         sg.Radio('nm', "RADIOCORR", default=params.lambda_units_template_crosscorr_nm, key='xcorr_template_wave_nm'),
         sg.Radio('a', "RADIOCORR", default=params.lambda_units_template_crosscorr_a, key='xcorr_template_wave_a'),
         sg.Radio('mu', "RADIOCORR", default=params.lambda_units_template_crosscorr_mu, key='xcorr_template_wave_mu')],
        [sg.Checkbox('Add broadening to template (km/s):', key='xcorr_smooth_template',
                     default=smooth_template_crosscorr),
         sg.InputText(smooth_value_crosscorr, size=(6,1), key='xcorr_smooth_template_value'),
         sg.Button('View template', button_color=('black', 'light blue'))],
        [sg.HorizontalSeparator()],

        [sg.Checkbox('Restrict wavelength range', key='xcorr_limit_wave_range', font = ('', default_size, 'bold'), default = xcorr_limit_wave_range), sg.Text('λ min'), sg.InputText(low_wave_corr, size=(8,1), key='xcorr_left_lambda'), sg.Text('λ max'), sg.InputText(high_wave_corr, size=(8,1), key='xcorr_right_lambda')],

        [sg.Radio('Velocity range (km/s):', "RADIOONLY2", default = is_vel_xcorr, key = 'is_vel_xcorr',tooltip='Use for peculiar Doppler motion or low (<0.01) redshift galaxies', font = ('', default_size, 'bold')), sg.InputText(low_vel_corr, size = (7,1), key = 'xcorr_low_vel'), sg.Text('-'), sg.InputText(high_vel_corr, size = (7,1), key = 'xcorr_high_vel'), sg.Text('Step (km/s)'), sg.InputText(xcorr_vel_step, key = 'xcorr_vel_step', size = (7,1))],
        [sg.Radio('Redshift range (z):', "RADIOONLY2", default = is_z_xcorr, key = 'is_z_xcorr',tooltip='Use for cosmological redshift estimation, where the velocity has no physical meaning', font = ('', default_size, 'bold')), sg.InputText(low_z_corr, size = (7,1), key = 'low_z_corr'), sg.Text('-'), sg.InputText(high_z_corr, size = (7,1), key = 'high_z_corr'), sg.Text('Step (z)'), sg.InputText(xcorr_z_step, key = 'xcorr_z_step', size = (7,1))],
        [sg.Push(), sg.Button('Confirm',button_color= ('white','black'), size = (18,1))]
    ]

    print('*** Cross-corr parameters window open. The main panel will be inactive until you close the window ***')

    xcorr_window = sg.Window('Cross-correlation parameters', xcorr_layout)

    while True:
        xcorr_event, xcorr_values = xcorr_window.read()

        if xcorr_event == sg.WIN_CLOSED:
            break

        # Extract updated values
        is_vel_xcorr = xcorr_values['is_vel_xcorr']
        is_z_xcorr = xcorr_values['is_z_xcorr']
        template_crosscorr = xcorr_values['xcorr_template']
        lambda_units_template_crosscorr = 'nm' if xcorr_values['xcorr_template_wave_nm'] else \
                                          'a' if xcorr_values['xcorr_template_wave_a'] else \
                                          'mu' if xcorr_values['xcorr_template_wave_mu'] else 'nm'
        smooth_template_crosscorr = xcorr_values['xcorr_smooth_template']

        try:
            smooth_value_crosscorr = float(xcorr_values['xcorr_smooth_template_value']) if smooth_template_crosscorr else 0
            if smooth_value_crosscorr < 0:
                sg.popup('Invalid smooth value for template. Must be >= 0!')
                continue
        except ValueError:
            sg.popup('Smooth template value not valid!')
            continue




        # try:
        #     low_wave_corr = float(xcorr_values['xcorr_left_lambda'])
        #     high_wave_corr = float(xcorr_values['xcorr_right_lambda'])
        #     wave_interval_corr = np.array([low_wave_corr,high_wave_corr])
        # except ValueError:
        #     sg.popup('Limit wave values not valid!')
        #     continue


        #cheching the input wavelength range
        xcorr_limit_wave_range = xcorr_values['xcorr_limit_wave_range']
        if xcorr_limit_wave_range:
            try:
                low_wave_corr = float(xcorr_values['xcorr_left_lambda'])
                high_wave_corr = float(xcorr_values['xcorr_right_lambda'])
                wave_interval_corr = np.array([low_wave_corr,high_wave_corr])
                # real_low_wave_corr = np.min(wave_interval_corr)
                # real_high_wave_corr = np.max(wave_interval_corr)
            except Exception:
                sg.popup('Limit wave values not valid!')
                continue





        if is_vel_xcorr:
            try:
                low_vel_corr = float(xcorr_values['xcorr_low_vel'])
                high_vel_corr = float(xcorr_values['xcorr_high_vel'])
                xcorr_vel_step = float(xcorr_values['xcorr_vel_step'])
                if abs(low_vel_corr - high_vel_corr) < 4:
                    sg.popup('Velocity interval too small')
                    continue
            except ValueError:
                sg.popup('Limit velocity values not valid!')
                continue

        if is_z_xcorr:
            try:
                low_z_corr = float(xcorr_values['low_z_corr'])
                high_z_corr = float(xcorr_values['high_z_corr'])
                xcorr_z_step = float(xcorr_values['xcorr_z_step'])

                if low_z_corr < 0 or high_z_corr < 0:
                    sg.popup('Redshift values must be greater than zero!')
                    continue
                if high_z_corr > 10:
                    sg.popup('The maximum redshift available is 10!')
                    continue
                if abs(low_z_corr - high_z_corr) < 0.001:
                    sg.popup('Redshift interval too small!')
                    continue
            except ValueError:
                sg.popup('Redshift values not valid!')
                continue

        if xcorr_event == 'View template':
            if not os.path.isfile(template_crosscorr):
                sg.popup('The template does not exist. I have nothing to show!')
                continue
            try:
                wave_template, flux_template, _, _ = stm.read_spec(template_crosscorr, lambda_units_template_crosscorr)
                if smooth_value_crosscorr > 0:
                    flux_template = spman.sigma_broad(wave_template, flux_template, smooth_value_crosscorr)
            except Exception:
                sg.popup('Cannot read the template!')
                continue

            plt.title(template_crosscorr)
            plt.plot(wave_template, flux_template)
            plt.xlim(low_wave_corr, high_wave_corr)
            plt.xlabel('Wavelength (A)')
            plt.ylabel('Normalised Flux')
            plt.show()
            plt.close()

        if xcorr_event == 'Confirm':
            print('Cross-corr parameters confirmed. The main panel is now active again\n')
            break

    xcorr_window.close()

    # Update params with modified values
    return replace(params,
                   template_crosscorr=template_crosscorr,
                   lambda_units_template_crosscorr_nm=lambda_units_template_crosscorr == 'nm',
                   lambda_units_template_crosscorr_a=lambda_units_template_crosscorr == 'a',
                   lambda_units_template_crosscorr_mu=lambda_units_template_crosscorr == 'mu',
                   lambda_units_template_crosscorr = lambda_units_template_crosscorr,
                   smooth_template_crosscorr=smooth_template_crosscorr,
                   smooth_value_crosscorr=smooth_value_crosscorr,
                   low_wave_corr=low_wave_corr,
                   high_wave_corr=high_wave_corr,
                   is_vel_xcorr=is_vel_xcorr,
                   low_vel_corr=low_vel_corr,
                   high_vel_corr=high_vel_corr,
                   is_z_xcorr=is_z_xcorr,
                   low_z_corr=low_z_corr,
                   high_z_corr=high_z_corr,
                   xcorr_limit_wave_range = xcorr_limit_wave_range,
                   xcorr_vel_step = xcorr_vel_step,
                   xcorr_z_step = xcorr_z_step)



def sigma_parameters(params: SpectraParams) -> SpectraParams:
    """Handles sigma measurement parameter input via GUI."""

    # Extract parameters from params
    template_sigma = params.template_sigma
    lambda_units_template_sigma = 'nm' if params.lambda_units_template_sigma_nm else \
                                  'a' if params.lambda_units_template_sigma_a else \
                                  'mu' if params.lambda_units_template_sigma_mu else 'nm'
    resolution_template = params.resolution_template
    resolution_spec = params.resolution_spec
    band_cat = params.band_cat
    band_halpha = params.band_halpha
    band_nad = params.band_nad
    band_h = params.band_h
    band_k = params.band_k
    band_custom = params.band_custom
    low_wave_sigma = params.low_wave_sigma
    high_wave_sigma = params.high_wave_sigma
    low_wave_cont = params.low_wave_cont
    high_wave_cont = params.high_wave_cont
    band_sigma = params.band_sigma
    cont_sigma = params.cont_sigma

    layout, scale_win, fontsize, default_size = misc.get_layout()
    sg.theme('LightBlue1')

    sigma_layout = [
        [sg.Text('Select a template:', font=('Helvetica', 10, 'bold')),
         sg.InputText(template_sigma, size=(55, 1), key='template_sigma'),
         sg.FileBrowse(tooltip='Load a template')],
        [sg.Text('Template wavelength is in:'),
         sg.Radio('nm', "RADIOSIGMA", default=params.lambda_units_template_sigma_nm, key='sigma_template_wave_nm'),
         sg.Radio('a', "RADIOSIGMA", default=params.lambda_units_template_sigma_a, key='sigma_template_wave_a'),
         sg.Radio('mu', "RADIOSIGMA", default=params.lambda_units_template_sigma_mu, key='sigma_template_wave_mu')],
        [sg.Text('Resolution (R) of the template (0 if (E)MILES)'),
         sg.InputText(resolution_template, size=(5, 1), key='sigma_res_template'),
         sg.Push(), sg.Button('View template', button_color=('black', 'light blue'))],
        [sg.HorizontalSeparator()],
        [sg.Text('Pre-loaded bands to fit for sigma:', font=('Helvetica', 10, 'bold')),
         sg.Radio('CaT', "RADIOBAND", default=band_cat, key='sigma_band_cat'),
         sg.Radio('Ha', "RADIOBAND", default=band_halpha, key='sigma_band_ha'),
         sg.Radio('Nad', "RADIOBAND", default=band_nad, key='sigma_band_nad'),
         sg.Radio('H band', "RADIOBAND", default=band_h, key='sigma_band_h'),
         sg.Radio('K band', "RADIOBAND", default=band_k, key='sigma_band_k')],
        [sg.Radio('Fitting a custom band', "RADIOBAND", default=band_custom, key='sigma_custom_band'),
         sg.Text('Wave interval (A)'),
         sg.InputText(low_wave_sigma, size=(5, 1), key='sigma_left_lambda'),
         sg.Text('-'),
         sg.InputText(high_wave_sigma, size=(5, 1), key='sigma_right_lambda'),
         sg.Text('Cont. interval'),
         sg.InputText(low_wave_cont, size=(5, 1), key='sigma_cont_left_lambda'),
         sg.Text('-'),
         sg.InputText(high_wave_cont, size=(5, 1), key='sigma_cont_right_lambda')],
        [sg.HorizontalSeparator()],
        [sg.Text('Resolution of the spectrum (R)'),
         sg.InputText(resolution_spec, size=(5, 1), key='sigma_spec_res')],
        [sg.Push(), sg.Button('Confirm', button_color=('white', 'black'), size=(18, 1))]
    ]

    print('*** Sigma parameters window open. The main panel will be inactive until you close the window ***')

    sigma_window = sg.Window('Sigma parameters', sigma_layout)

    while True:
        sigma_event, sigma_values = sigma_window.read()

        if sigma_event == sg.WIN_CLOSED:
            break

        # Update parameters from GUI input
        lambda_units_template_sigma = 'nm' if sigma_values['sigma_template_wave_nm'] else \
                                      'a' if sigma_values['sigma_template_wave_a'] else \
                                      'mu' if sigma_values['sigma_template_wave_mu'] else 'nm'
        template_sigma = sigma_values['template_sigma']
        band_cat = sigma_values['sigma_band_cat']
        band_halpha = sigma_values['sigma_band_ha']
        band_nad = sigma_values['sigma_band_nad']
        band_h = sigma_values['sigma_band_h']
        band_k = sigma_values['sigma_band_k']
        band_custom = sigma_values['sigma_custom_band']

        # Assign predefined bands
        predefined_bands = {
            "cat": ([8440., 8702.], [8560., 8604.]),
            "halpha": ([6420., 6610.], [6510., 6540.]),
            "nad": ([5600., 6150.], [5910., 5970.]),
            "h": ([16600., 17200.], [16930., 17040.]),
            "k": ([22700., 23700.], [22700., 22800.])
        }
        for band, (sigma, cont) in predefined_bands.items():
            if locals()[f"band_{band}"]:
                band_sigma, cont_sigma = np.array(sigma), np.array(cont)

        if band_custom:
            try:
                low_wave_sigma = float(sigma_values['sigma_left_lambda'])
                high_wave_sigma = float(sigma_values['sigma_right_lambda'])
                low_wave_cont = float(sigma_values['sigma_cont_left_lambda'])
                high_wave_cont = float(sigma_values['sigma_cont_right_lambda'])

                band_sigma = np.array([low_wave_sigma, high_wave_sigma])
                cont_sigma = np.array([low_wave_cont, high_wave_cont])
            except ValueError:
                sg.popup('Band values for sigma not valid!')
                continue

        try:
            resolution_spec = float(sigma_values['sigma_spec_res'])
            resolution_template = float(sigma_values['sigma_res_template'])
            if resolution_spec <= 0 or resolution_template < 0:
                sg.popup('Invalid resolution values for the spectrum or the template')
                continue
        except ValueError:
            sg.popup('Resolution values for sigma not valid!')
            continue

        if sigma_event == 'View template':
            if not os.path.isfile(template_sigma):
                sg.popup('The template does not exist. I have nothing to show!')
                continue

            # Reading the template
            try:
                wave_template, flux_template, _, _ = stm.read_spec(template_sigma, lambda_units_template_sigma)
            except Exception:
                sg.popup('Cannot read the template!')
                continue

            plt.title(template_sigma)
            plt.plot(wave_template, flux_template)
            plt.xlim(band_sigma[0], band_sigma[1])
            plt.xlabel('Wavelength (A)')
            plt.ylabel('Normalised Flux')
            plt.show()
            plt.close()

        if sigma_event == 'Confirm':
            print('Sigma parameters confirmed. The main panel is now active again\n')
            break

    sigma_window.close()

    # Update params with modified values
    return replace(params,
                   template_sigma=template_sigma,
                   lambda_units_template_sigma_nm=lambda_units_template_sigma == 'nm',
                   lambda_units_template_sigma_a=lambda_units_template_sigma == 'a',
                   lambda_units_template_sigma_mu=lambda_units_template_sigma == 'mu',
                   resolution_template=resolution_template,
                   resolution_spec=resolution_spec,
                   band_cat=band_cat,
                   band_halpha=band_halpha,
                   band_nad=band_nad,
                   band_h=band_h,
                   band_k=band_k,
                   band_custom=band_custom,
                   band_sigma=band_sigma,
                   cont_sigma=cont_sigma,
                   low_wave_sigma=low_wave_sigma,
                   high_wave_sigma=high_wave_sigma,
                   low_wave_cont=low_wave_cont,
                   high_wave_cont=high_wave_cont)



def line_strength_parameters(params: SpectraParams) -> SpectraParams:

    """
    Opens the line strength measurement parameters GUI and updates the values in params.
    """

    # Extract parameters from params
    index_file = params.index_file
    have_index_file = params.have_index_file
    single_index = params.single_index
    idx_left_blue = params.idx_left_blue
    idx_right_blue = params.idx_right_blue
    idx_left_red = params.idx_left_red
    idx_right_red = params.idx_right_red
    idx_left_line = params.idx_left_line
    idx_right_line = params.idx_right_line
    index_usr = params.index_usr
    lick_ew = params.lick_ew
    lick_constant_fwhm = params.lick_constant_fwhm
    lick_constant_r = params.lick_constant_r
    spec_lick_res_fwhm = params.spec_lick_res_fwhm
    spec_lick_res_r = params.spec_lick_res_r
    lick_correct_emission = params.lick_correct_emission
    z_guess_lick_emission = params.z_guess_lick_emission
    correct_ew_sigma = params.correct_ew_sigma
    radio_lick_sigma_auto = params.radio_lick_sigma_auto
    radio_lick_sigma_single = params.radio_lick_sigma_single
    sigma_single_lick = params.sigma_single_lick
    radio_lick_sigma_list = params.radio_lick_sigma_list
    sigma_lick_file = params.sigma_lick_file
    stellar_parameters_lick = params.stellar_parameters_lick
    dop_correction_lick = params.dop_correction_lick
    lick_ssp_models = params.lick_ssp_models
    ssp_model = params.ssp_model
    interp_modes = params.interp_modes
    interp_model = params.interp_model
    sigma_coeff = params.sigma_coeff
    sigma_corr = params.sigma_corr
    stellar_spectra_coeff_file = params.stellar_spectra_coeff_file
    lambda_units_coeff_nm = params.lambda_units_coeff_nm
    lambda_units_coeff_a = params.lambda_units_coeff_a
    lambda_units_coeff_mu = params.lambda_units_coeff_mu
    lambda_units_coeff = params.lambda_units_coeff
    smooth_stellar_sample = params.smooth_stellar_sample
    smooth_value_sample = params.smooth_value_sample
    same_idx_ew_task = params.same_idx_ew_task
    have_index_file_corr = params.have_index_file_corr
    index_file_corr = params.index_file_corr
    single_index_corr = params.single_index_corr
    idx_left_blue_sigma = params.idx_left_blue_sigma
    idx_right_blue_sigma = params.idx_right_blue_sigma
    idx_left_red_sigma = params.idx_left_red_sigma
    idx_right_red_sigma = params.idx_right_red_sigma
    idx_left_line_sigma = params.idx_left_line_sigma
    idx_right_line_sigma = params.idx_right_line_sigma
    sigma_vel_file = params.sigma_vel_file
    ew_list_file = params.ew_list_file
    sigma_coeff_file = params.sigma_coeff_file
    result_sigma_coeff_dir = params.result_sigma_coeff_dir
    spectra_list_name = params.spectra_list_name
    result_plot_dir = params.result_plot_dir
    result_ew_data_dir = params.result_ew_data_dir


    layout, scale_win, fontsize, default_size = misc.get_layout()
    timestamp = time.strftime("%Y%m%d_%H%M%S")

    sg.theme('LightBlue1')
    ew_layout = [
        [sg.Radio('User indices on a list file:', "RADIOEW", default = have_index_file, key = 'ew_idx_file',  font = ('', default_size, 'bold'),tooltip='You need an ASCII file with the index definition. See the readme to how to build this file'), sg.InputText(index_file, size = (48,1), key = 'idx_file'), sg.FileBrowse(size = (10,1))],
        [sg.HorizontalSeparator()],
        [sg.Radio('Single index', "RADIOEW", default = single_index, key = 'ew_single_idx',font = ('', default_size, 'bold'),tooltip='Enter the index definition wavelengths in the windows on the right'), sg.Text('blue cont.:'), sg.InputText(idx_left_blue, size = (5,1), key = 'left_wave_blue_cont'), sg.Text('-'), sg.InputText(idx_right_blue, size = (5,1), key = 'right_wave_blue_cont'), sg.Text('red cont.:'), sg.InputText(idx_left_red, size = (5,1), key = 'left_wave_red_cont'), sg.Text('-'),  sg.InputText(idx_right_red, size = (5,1), key = 'right_wave_red_cont'), sg.Text('line:'), sg.InputText(idx_left_line, size = (5,1), key = 'left_line'), sg.Text('-'), sg.InputText(idx_right_line, size = (5,1), key = 'right_line')],
        [sg.HorizontalSeparator()],

        #Lick/IDS indices
        [sg.Radio('Lick/IDS indices:', "RADIOEW", default = lick_ew, key = 'ew_lick', font = ('', default_size, 'bold'),tooltip='Measure the classical Lick/IDS indices. You do not need to enter their definitions'), sg.Radio('Constant resolution FWHM:', "RADIOLICKRES", default = lick_constant_fwhm, key ='lick_constant_fwhm'), sg.InputText(spec_lick_res_fwhm, key = 'spec_lick_res_fwhm',size = (7,1)), sg.Radio('Constant resolution R:', "RADIOLICKRES", default = lick_constant_r, key ='lick_constant_r'), sg.InputText(spec_lick_res_r, key = 'spec_lick_res_r',size = (8,1))],
        [sg.Text(' '), sg.Checkbox('Emission line(s) correction:', default = lick_correct_emission, key = 'lick_correct_emission'), sg.Text('Redshift guess'), sg.InputText(z_guess_lick_emission, key = 'z_guess_lick_emission', size = (15,1)), sg.Text(' '), sg.Checkbox('Perform Doppler correction', default = dop_correction_lick, key = 'dop_correction_lick')],
        [sg.Text(' '), sg.Checkbox('Correct for sigma:', default = correct_ew_sigma, key = 'correct_ew_sigma'), sg.Radio('Auto', "RADIOLICKSIGMA", default = radio_lick_sigma_auto, key = 'radio_lick_sigma_auto'), sg.Radio('Single (km/s):', "RADIOLICKSIGMA", default = radio_lick_sigma_single, key = 'radio_lick_sigma_single'), sg.InputText(sigma_single_lick, size = (5,1), key = 'sigma_single_lick'), sg.Radio('List:', "RADIOLICKSIGMA", default = radio_lick_sigma_list, key = 'radio_lick_sigma_list'), sg.InputText(sigma_lick_file, key = 'sigma_lick_file', size = (17,1)), sg.FileBrowse(size = (10,1)) ],
        [sg.Text(' '), sg.Checkbox('Estimate stellar parameters with SSP models:', default = stellar_parameters_lick, key = 'stellar_parameters_lick',tooltip='Perform interpolation with SSP model grids to retrieve age, metellicity and alpha enhancement. The Thomas2010 models are the most reliable'), sg.InputCombo(lick_ssp_models, key='ssp_model',default_value=ssp_model, readonly=True, size = (14,1)), sg.Text('Interpolation mode:',tooltip='Interpolate linearly with griddata function or with machine learning Gaussian Process Regression (GPR)'), sg.InputCombo(interp_modes, key='interp_model',default_value=interp_model, readonly=True, size = (14,1))],
        [sg.Push(), sg.Button("Help", size=(12, 1),button_color=('black','orange')), sg.Button('Confirm',button_color= ('white','black'), size = (18,1))],
        [sg.HorizontalSeparator()],

        #7) Determination of the velocity disperion coefficients. Stand alone: need a separate sample of spectra
        [sg.Checkbox('Calculate velocity dispersion coefficients', font = ('', default_size, 'bold'), key = 'sigma_coeff',tooltip='Determination of the velocity dispersion coefficients to correct the EW of the indices of the spectrum', default = sigma_coeff), sg.Push(), sg.Button('Sigma coeff parameters',button_color= ('black','light blue'), size = (22,1)), sg.Button('Compute!',button_color=('black','light grey'), size = (12,1))],

        #8) Apply the velocity dispersion correction coefficientd to the loaded spectra
        [sg.Checkbox('Correct the line-strength for velocity dispersion', font = ('', default_size, 'bold'), key = 'sigma_corr',tooltip='Correct the raw measured EW of the spectra to a zero velocity dispersion frame, using the coefficients calculated in the Sigma coeff determination task', default = sigma_corr), sg.Push(), sg.Button('Sigma corr parameters',button_color= ('black','light blue'), size = (22,1)), sg.Button('Correct!',button_color=('black','light grey'), size = (12,1)) ]

        ]

    print ('*** Line-strength parameters window open. The main panel will be inactive until you close the window ***')
    ew_window = sg.Window('Line-strength parameters', ew_layout)

    while True:
        ew_event, ew_values = ew_window.read()

        if ew_event == sg.WIN_CLOSED:
            break

        # Retrieving the parameters of the GUI
        index_file = ew_values['idx_file']
        have_index_file = ew_values['ew_idx_file']
        single_index = ew_values['ew_single_idx']
        lick_ew = ew_values['ew_lick']
        lick_constant_fwhm = ew_values['lick_constant_fwhm']
        lick_constant_r = ew_values['lick_constant_r']
        lick_correct_emission = ew_values['lick_correct_emission']
        radio_lick_sigma_single = ew_values['radio_lick_sigma_single']
        radio_lick_sigma_list = ew_values['radio_lick_sigma_list']
        radio_lick_sigma_auto = ew_values['radio_lick_sigma_auto']
        sigma_lick_file = ew_values['sigma_lick_file']
        correct_ew_sigma = ew_values['correct_ew_sigma']
        stellar_parameters_lick = ew_values['stellar_parameters_lick']
        dop_correction_lick = ew_values['dop_correction_lick']
        ssp_model = ew_values['ssp_model']
        interp_model = ew_values['interp_model']
        sigma_coeff = ew_values['sigma_coeff']
        sigma_corr = ew_values['sigma_corr']

        try:
            idx_left_blue = float(ew_values['left_wave_blue_cont'])
            idx_right_blue = float(ew_values['right_wave_blue_cont'])
            idx_left_red = float(ew_values['left_wave_red_cont'])
            idx_right_red = float(ew_values['right_wave_red_cont'])
            idx_left_line = float(ew_values['left_line'])
            idx_right_line = float(ew_values['right_line'])
            spec_lick_res_fwhm = float(ew_values['spec_lick_res_fwhm'])
            spec_lick_res_r = float(ew_values['spec_lick_res_r'])
            sigma_single_lick = float(ew_values['sigma_single_lick'])
            z_guess_lick_emission = float(ew_values['z_guess_lick_emission'])

            #building the index
            index_usr = np.array([idx_left_blue, idx_right_blue, idx_left_red, idx_right_red, idx_left_line, idx_right_line]).T
        except ValueError:
            sg.popup('Index values for EW measurement not valid!')
            continue

        if lick_constant_fwhm and (spec_lick_res_fwhm < 0):
            sg.popup ('FWHM resolution must be positive!')
            continue

        if lick_constant_r and (spec_lick_res_r < 0):
            sg.popup ('R resolution must be positive!')
            continue

        if correct_ew_sigma and radio_lick_sigma_single and sigma_single_lick < 0:
            sg.popup ('Sigma value must be positive!')
            continue

        #Check if the sigma file exist
        if correct_ew_sigma and radio_lick_sigma_list:

            #test if file file exist
            cond0000 = (os.path.isfile(sigma_lick_file))
            if not cond0000:
                sg.popup('The sigma list file does not exist. Skipping...')
                continue

        #CLOSING THE WINDOW ONCE I CLICK "CONFIRM"
        if ew_event == 'Confirm':
            print ('Line-strength parameters confirmed. This main panel is now active again')
            print ('')
            break


    #A) SIGMA COEFF DETERMINATION
        if ew_event == 'Sigma coeff parameters':
            sg.theme('LightBlue1')
            sigmacorr_layout = [
                [sg.Radio('Index list file:', "RADIOCOEFF1", default = have_index_file_corr, key = 'ew_corr_idx_file'), sg.InputText(index_file_corr, size = (25,1), key = 'idx_corr_file'), sg.FileBrowse(tooltip='Load a list file containing the names of the indices')],
                [sg.Radio('Usr index', "RADIOCOEFF1", default = single_index_corr, key = 'ew_corr_single_idx'), sg.Text('blue continuum:'), sg.InputText(idx_left_blue_sigma, size = (5,1), key = 'left_wave_blue_cont'), sg.Text('-'), sg.InputText(idx_right_blue_sigma, size = (5,1), key = 'right_wave_blue_cont'), sg.Text('red continuum:'), sg.InputText(idx_left_red_sigma, size = (5,1), key = 'left_wave_red_cont'), sg.Text('-'),  sg.InputText(idx_right_red_sigma, size = (5,1), key = 'right_wave_red_cont'), sg.Text('line:'), sg.InputText(idx_left_line_sigma, size = (5,1), key = 'left_line'), sg.Text('-'), sg.InputText(idx_right_line_sigma, size = (5,1), key = 'right_line')],
                [sg.Text('Sample spectra list:'), sg.InputText(stellar_spectra_coeff_file, size = (17,1), key = 'sigma_coeff_sample_list'), sg.FileBrowse(tooltip='Load a list file containing the names of the sample spectra'), sg.Radio('nm', "RADIOCOEFF", key = 'sigma_coeff_sample_list_wave_nm', default = lambda_units_coeff_nm), sg.Radio('A', "RADIOCOEFF", key = 'sigma_coeff_sample_list_wave_a', default = lambda_units_coeff_a), sg.Radio('mu', "RADIOCOEFF", default = lambda_units_coeff_mu, key = 'sigma_coeff_sample_list_wave_mu'), sg.Checkbox('Add a sigma (km/s): ', key = 'sigma_coeff_sample_smooth', default = smooth_stellar_sample ), sg.InputText(smooth_value_sample, size = (5,1), key = 'sigma_coeff_sample_smooth_sigma')],

                [sg.Push(), sg.Button('Confirm',button_color= ('white','black'), size = (18,1))]
                ]

            sigmacorr_window = sg.Window('Sigma coeff parameters', sigmacorr_layout)

            print (single_index_corr)
            while True:

                sigmacorr_event, sigmacorr_values = sigmacorr_window.read()

                if sigmacorr_event == sg.WIN_CLOSED:
                    break

                # Retrieving the parameters of the GUI
                index_file_corr = sigmacorr_values['idx_corr_file']
                have_index_file_corr = sigmacorr_values['ew_corr_idx_file']
                single_index_corr = sigmacorr_values['ew_corr_single_idx']

                try:
                    idx_left_blue_sigma = float(sigmacorr_values['left_wave_blue_cont'])
                    idx_right_blue_sigma = float(sigmacorr_values['right_wave_blue_cont'])
                    idx_left_red_sigma = float(sigmacorr_values['left_wave_red_cont'])
                    idx_right_red_sigma = float(sigmacorr_values['right_wave_red_cont'])
                    idx_left_line_sigma = float(sigmacorr_values['left_line'])
                    idx_right_line_sigma = float(sigmacorr_values['right_line'])
                    #building the index
                    index_usr_corr = np.array([idx_left_blue_sigma, idx_right_blue_sigma, idx_left_red_sigma, idx_right_red_sigma, idx_left_line_sigma, idx_right_line_sigma]).T #not useful anymore since I have the dataclass now
                except ValueError:
                    sg.popup('Index values not valid!')
                    continue

                stellar_spectra_coeff_file = sigmacorr_values['sigma_coeff_sample_list']
                lambda_units_coeff_nm = sigmacorr_values['sigma_coeff_sample_list_wave_nm']
                lambda_units_coeff_a = sigmacorr_values['sigma_coeff_sample_list_wave_a']
                lambda_units_coeff_mu = sigmacorr_values['sigma_coeff_sample_list_wave_mu']
                #assigning lambda units of the spectra
                if (lambda_units_coeff_nm):
                    lambda_units_coeff = 'nm'
                if (lambda_units_coeff_a):
                    lambda_units_coeff = 'a'
                if (lambda_units_coeff_mu):
                    lambda_units_coeff = 'mu'

                smooth_stellar_sample = sigmacorr_values['sigma_coeff_sample_smooth']
                if smooth_stellar_sample:
                    smooth_value_sample = float(sigmacorr_values['sigma_coeff_sample_smooth_sigma'])
                else:
                    smooth_value_sample = 0.

                if sigmacorr_event == 'Confirm':
                    print ('sigma coeff parameters confirmed. This main panel is now active again')
                    print ('')
                    break

            sigmacorr_window.close()


        #B) CORRECT EWS FOR SIGMA
        if ew_event == 'Sigma corr parameters':
            sg.theme('LightBlue1')
            correw_layout = [
                [sg.Text('Sigma list file'), sg.InputText(sigma_vel_file, size = (14,1), key = 'sigma_file'), sg.FileBrowse(tooltip='Load the ASCII file of the results of the velocity dispersion task, that contains: Spec. name, sigma, err_sigma'), sg.Text('EW file to correct:'), sg.InputText(ew_list_file, size = (14,1), key = 'ew_file_to_correct'), sg.FileBrowse(tooltip='Load the ASCII file of the results of the equivalent width task, that contains: Spec. name, EW, err_ew'), sg.Text('Correction coeff. file'), sg.InputText(sigma_coeff_file, size = (14,1), key = 'coeff_sigma_file'), sg.FileBrowse(tooltip='Load the ASCII file of the results of the Sigma coeff determination task, that contains: coeff. indices, coeff. errors')] ,
                [sg.Push(), sg.Button('Confirm',button_color= ('white','black'), size = (18,1))]
                ]

            print ('*** Sigma corr parameters window open. The main panel will be inactive until you close the window ***')
            correw_window = sg.Window('Sigma correction parameters', correw_layout)

            while True:

                correw_event, correw_values = correw_window.read()

                if correw_event == sg.WIN_CLOSED:
                    break

                # Retrieving the parameters of the GUI
                sigma_vel_file = correw_values['sigma_file']
                ew_list_file = correw_values['ew_file_to_correct']
                sigma_coeff_file = correw_values['coeff_sigma_file']

                if correw_event == 'Confirm':
                    print ('sigma corr parameters confirmed. This main panel is now active again')
                    print ('')
                    break

            correw_window.close()

        # Check conditions
        if (not sigma_coeff and ew_event == 'Compute!'):
            sg.popup ('Cannot compute until you activate the Sigma coeff determination task!')

        if (not sigma_corr and ew_event == 'Correct!'):
            sg.popup ('Cannot correct until you activate the Correct EWs for sigma task!')


        # 8) SIGMA COEFF DETERMINATION
        if (ew_event == 'Compute!' and sigma_coeff):
            #preparing the files to be saved
            if (sigma_coeff and single_index_corr):

                #verify the index definition is in the correct sequence
                if (index_usr_corr[0] > index_usr_corr[1] or index_usr_corr[2] >index_usr_corr[3] or index_usr_corr[4] > index_usr_corr[5]):
                    sg.popup('It seems we have a problem. Did you invert the wavelengths of the indices?')
                    continue

                # Retrieving the parameters of the GUI
                coeff_file = result_sigma_coeff_dir+'/'+spectra_list_name+'_sigma_coeff_' +timestamp + '.dat'
                coeff_id = ['sigma_spline', 'err_spline']
                coeff_number = 4
                coeff_values = np.zeros(coeff_number)
                err_coeff_values = np.zeros(coeff_number)
                coeff_data_array = np.column_stack((coeff_values, err_coeff_values))

                #generating the dataframe and adding the data
                df_coeff = pd.DataFrame(coeff_data_array, columns = coeff_id)

                #writing to a file
                df_coeff.to_csv(coeff_file, index= True, sep=' ')

            elif (sigma_coeff and not single_index_corr):
                cond33 = (os.path.isfile(index_file_corr))
                if not cond33:
                    sg.popup('The index file does not exist. Skipping...')
                    continue

                #exploring the index file for errors
                try:
                    idx_names, indices = ls.read_idx(index_file_corr)
                except ValueError:
                    sg.popup('At least one index in the file is not valid')
                    continue

                if len(indices[:,0]) < 6:
                    sg.popup ('The length of at least one index is not correct')
                    continue

                coeff_file = result_sigma_coeff_dir+'/'+spectra_list_name+'_sigma_coeff_' +timestamp + '.dat'
                coeff_number = 4

                id_array, index = ls.read_idx(index_file_corr)
                num_indices = len(id_array)

                shape = (coeff_number, num_indices)
                coeff_all = np.zeros(shape)
                coeff_err_all = np.zeros(shape)

                err_col_type = np.chararray(num_indices)
                err_col_type = 'e'
                err_col_names = np.char.add(id_array, err_col_type)

                #merge the column array names
                col_names = np.concatenate((id_array, err_col_names))
                coeff_id = col_names
                coeff_data = np.column_stack((coeff_all, coeff_err_all))

                df_coeff = pd.DataFrame(coeff_data, columns = coeff_id)
                df_coeff.to_csv(coeff_file, index= True, sep=' ')


            #1) If I want to measure just one index
            if (sigma_coeff and single_index_corr):

                cond666 = (os.path.isfile(stellar_spectra_coeff_file))
                if not cond666:
                    sg.popup('Stellar spectra file does not exist. Skipping...')
                    continue


                #check to add the absolute path in case the spectra list is given in relative path
                with open(stellar_spectra_coeff_file, 'r') as f:
                    spec_names_sample = []
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith("#"):
                            continue

                        # If I gave the relative path, try to solve the absolute path
                        if not os.path.isabs(line):
                            spec_names_sample.append(os.path.join(BASE_DIR, line))

                        #If I have the absolute path
                        else:
                            # Loading the spectra list
                            spec_names_sample = np.loadtxt(stellar_spectra_coeff_file, dtype = 'str', delimiter = ' ', usecols=[0])
                            break

                # Check if the sample spectra really exist
                sample_spec_not_exist = [spec for spec in spec_names_sample if not os.path.isfile(spec)]

                # If any spectra are missing, show a warning
                if sample_spec_not_exist:
                    sg.popup(
                        "Warning: the following spectra do not exist! Delete them from the list or I will crash!",
                        sample_spec_not_exist
                    )
                    continue

                print ('Running Sigma coefficient determination task...')

                idx_array, ew_coeff_array, err_coeff_array, ew_mean, ew_std, stop_cond = ls.sigma_coeff (stellar_spectra_coeff_file, index_usr_corr, lambda_units_coeff, True, False, smooth_value_sample, True, result_plot_dir)

                df_coeff.iloc[:, 0]= ew_coeff_array
                df_coeff.iloc[:, 1]= err_coeff_array
                df_coeff.columns = ['sigma_spline', 'err_spline']
                df_coeff.to_csv(coeff_file, index= False, sep=' ')

                print ('Sigma correction coefficients for usr index')
                print (ew_coeff_array)
                print ('')
                print ('Error correction coefficients for usr index')
                print (err_coeff_array)
                print('')
                print('I saved the plots for you in the ', result_plot_dir, ' folder' )
                print('')
                print ('File saved: ', coeff_file)


            if (sigma_coeff and not single_index_corr):

                #reading the index file
                try:
                    id_array, index = ls.read_idx(index_file_corr)
                    num_indices = len(id_array)
                except Exception:
                    sg.popup ('The index file does not exist!')
                    continue

                cond666 = (os.path.isfile(stellar_spectra_coeff_file))
                if not cond666:
                    sg.popup('Stellar spectra file does not exist. Skipping...')
                    continue

                #check to add the absolute path in case the spectra list is given in relative path
                with open(stellar_spectra_coeff_file, 'r') as f:
                    spec_names_sample = []
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith("#"):
                            continue

                        # If I gave the relative path, try to solve the absolute path
                        if not os.path.isabs(line):
                            spec_names_sample.append(os.path.join(BASE_DIR, line))

                        #If I have the absolute path
                        else:
                            # Loading the spectra list
                            spec_names_sample = np.loadtxt(stellar_spectra_coeff_file, dtype = 'str', delimiter = ' ', usecols=[0])
                            break

                # Check if the spectra actually exist
                sample_spec_not_exist = [spec for spec in spec_names_sample if not os.path.isfile(spec)]

                # Show a warning if some spectra are missing
                if sample_spec_not_exist:
                    sg.popup(
                        "Warning: the following spectra do not exist! Delete them from the list or I will crash!",
                        sample_spec_not_exist
                    )
                    continue

                print ('Running Sigma coefficient determination task...')

                idx_array, ew_coeff_array, err_coeff_array, ew_mean, ew_std, stop_cond = ls.sigma_coeff (stellar_spectra_coeff_file, index_file_corr, lambda_units_coeff, False, False, smooth_value_sample, True, result_plot_dir)

                for k in range(num_indices):
                    ew_coeff_array_new = ew_coeff_array[:,k]
                    err_coeff_array_new = err_coeff_array[:,k]

                    #Updating the dataframe
                    df_coeff[coeff_id[k]]= ew_coeff_array_new
                    df_coeff[coeff_id[k+num_indices]] = err_coeff_array_new
                    df_coeff.to_csv(coeff_file, index= False, sep=' ')

                print ('Indices')
                print (idx_array)
                print ('')
                print ('Sigma correction coefficients for usr index')
                print (ew_coeff_array)
                print ('')
                print ('Error correction coefficients for usr index')
                print (err_coeff_array)
                print('')
                print('I saved the plots for you in the ', result_plot_dir, ' folder' )
                print('')
                print ('File saved: ', coeff_file)


        #9) CORRECT EWS FOR SIGMA
        if (ew_event == 'Correct!'):
            if sigma_corr:
                task_done2 = 1

                cond77 = (os.path.isfile(sigma_coeff_file))
                if not cond77:
                    sg.popup('Sigma coefficient file does not exist. Skipping...')
                    continue

                cond66 = (os.path.isfile(ew_list_file))
                if not cond66:
                    sg.popup('EW file does not exist. Skipping...')
                    continue

                cond88 = (os.path.isfile(sigma_vel_file))
                if not cond88:
                    sg.popup('Sigma velocities file does not exist. Skipping...')
                    continue

                print ('Running correct EWs task...')

                try:
                    all_idx_np, new_ew_data = ls.corr_ew(ew_list_file, sigma_coeff_file, sigma_vel_file)
                except ValueError:
                    sg.popup ('Cannot compare apples with oranges! The EW file and the correction coefficients file do not have the same indices. ')
                    continue

                #creating and saving here the file since is one shot only
                results_ew = result_ew_data_dir+'/'+spectra_list_name+'_sigma_corr_ew_data_' + timestamp + '.dat'
                df = pd.DataFrame(new_ew_data, columns = all_idx_np)
                df.to_csv(results_ew, index= False, sep=' ')
                print('EWs corrected saved to ', results_ew)
                sg.popup ('Succeed!')
        
        if ew_event == 'Help':
            f = open(os.path.join(BASE_DIR, "help_files", "linestrength.txt"), 'r')
            file_contents = f.read()
            if layout == layouts.layout_android:
                sg.popup_scrolled(file_contents, size=(120, 30))
            else:
                sg.popup_scrolled(file_contents, size=(100, 40))

    ew_window.close()

    # updating the params in the dataclass
    params = replace(params,
        index_file=index_file,
        have_index_file=have_index_file,
        single_index=single_index,
        idx_left_blue=idx_left_blue,
        idx_right_blue=idx_right_blue,
        idx_left_red=idx_left_red,
        idx_right_red=idx_right_red,
        idx_left_line=idx_left_line,
        idx_right_line=idx_right_line,
        lick_ew=lick_ew,
        lick_constant_fwhm=lick_constant_fwhm,
        lick_constant_r=lick_constant_r,
        spec_lick_res_fwhm=spec_lick_res_fwhm,
        spec_lick_res_r=spec_lick_res_r,
        lick_correct_emission=lick_correct_emission,
        z_guess_lick_emission=z_guess_lick_emission,
        correct_ew_sigma=correct_ew_sigma,
        radio_lick_sigma_auto=radio_lick_sigma_auto,
        radio_lick_sigma_single=radio_lick_sigma_single,
        sigma_single_lick=sigma_single_lick,
        radio_lick_sigma_list=radio_lick_sigma_list,
        sigma_lick_file=sigma_lick_file,
        stellar_parameters_lick=stellar_parameters_lick,
        dop_correction_lick=dop_correction_lick,
        lick_ssp_models=lick_ssp_models,
        ssp_model=ssp_model,
        interp_modes=interp_modes,
        interp_model=interp_model,
        sigma_coeff=sigma_coeff,
        sigma_corr=sigma_corr,
        stellar_spectra_coeff_file=stellar_spectra_coeff_file,
        lambda_units_coeff_nm=lambda_units_coeff_nm,
        lambda_units_coeff_a=lambda_units_coeff_a,
        lambda_units_coeff_mu=lambda_units_coeff_mu,
        lambda_units_coeff=lambda_units_coeff,
        smooth_stellar_sample=smooth_stellar_sample,
        smooth_value_sample=smooth_value_sample,
        same_idx_ew_task=same_idx_ew_task,
        have_index_file_corr=have_index_file_corr,
        index_file_corr=index_file_corr,
        single_index_corr=single_index_corr,
        idx_left_blue_sigma=idx_left_blue_sigma,
        idx_right_blue_sigma=idx_right_blue_sigma,
        idx_left_red_sigma=idx_left_red_sigma,
        idx_right_red_sigma=idx_right_red_sigma,
        idx_left_line_sigma=idx_left_line_sigma,
        idx_right_line_sigma=idx_right_line_sigma,
        sigma_vel_file=sigma_vel_file,
        ew_list_file=ew_list_file,
        sigma_coeff_file=sigma_coeff_file
    )

    return params



def line_fitting_parameters(params: SpectraParams) -> SpectraParams:

    """
    Opens the line fitting parameters GUI and updates the values in params.
    """

    # Extract parameters from params
    cat_band_fit = params.cat_band_fit
    usr_fit_line = params.usr_fit_line
    emission_line = params.emission_line
    low_wave_fit = params.low_wave_fit
    high_wave_fit = params.high_wave_fit
    y0 = params.y0
    x0 = params.x0
    a = params.a
    sigma = params.sigma
    m = params.m
    c = params.c

    layout, scale_win, fontsize, default_size = misc.get_layout()
    sg.theme('LightBlue1')
    linefit_layout = [
        [sg.Radio('Automatic fit of the CaT lines', "RADIOFILEFIT", key='cat_fit', default=cat_band_fit,
                  tooltip='Automatic CaT band fitting with three gaussian functions and a line')],
        [sg.Radio('Manual fit of the following line:', "RADIOFILEFIT", default=usr_fit_line, key='line_fit_single',
                  tooltip='Fitting of user defined line with the parameters on the right'),
         sg.Checkbox('Emission', default=emission_line, key='emission_line'),
         sg.Text('Wave:'), sg.InputText(low_wave_fit, size=(4, 1), key='left_wave_fitting'),
         sg.Text('-'), sg.InputText(high_wave_fit, size=(4, 1), key='right_wave_fitting'),
         sg.Text('Guess:', tooltip='Initial guess only for usr line, fitted with a gaussian and a line'),
         sg.Text('Y-off'), sg.InputText(y0, key='y0', size=(3, 1)),
         sg.Text('Line'), sg.InputText(x0, key='x0', size=(4, 1)),
         sg.Text('Height'), sg.InputText(a, key='a', size=(3, 1)),
         sg.Text('Sigma'), sg.InputText(sigma, key='sigma', size=(3, 1)),
         sg.Text('Slope'), sg.InputText(m, key='m', size=(2, 1)),
         sg.Text('Intercept'), sg.InputText(c, key='c', size=(2, 1))],
        [sg.Push(), sg.Button('Confirm', button_color=('white', 'black'), size=(18, 1))]
    ]

    print('*** Line fitting parameters window open. The main panel will be inactive until you close the window ***')
    linefit_window = sg.Window('Line(s) fitting parameters', linefit_layout)

    while True:
        linefit_event, linefit_values = linefit_window.read()

        if linefit_event == sg.WIN_CLOSED:
            break

        # Read values from GUI
        cat_band_fit = linefit_values['cat_fit']
        usr_fit_line = linefit_values['line_fit_single']
        emission_line = linefit_values['emission_line']

        try:
            low_wave_fit = float(linefit_values['left_wave_fitting'])
            high_wave_fit = float(linefit_values['right_wave_fitting'])
            y0 = float(linefit_values['y0'])
            x0 = float(linefit_values['x0'])
            a = float(linefit_values['a'])
            sigma = float(linefit_values['sigma'])
            m = float(linefit_values['m'])
            c = float(linefit_values['c'])

            wave_interval_fit = np.array([low_wave_fit, high_wave_fit])
            guess_param = [y0, x0, a, sigma, m, c]

        except ValueError:
            sg.popup('At least one of the parameters of the Line(s) fitting task is not valid!')
            continue

        if linefit_event == 'Confirm':
            print('Line fitting parameters confirmed. This main panel is now active again')
            print('')
            break

    linefit_window.close()

    # Update only the modified parameters in params
    params = replace(params,
        cat_band_fit=cat_band_fit,
        usr_fit_line=usr_fit_line,
        emission_line=emission_line,
        low_wave_fit=low_wave_fit,
        high_wave_fit=high_wave_fit,
        y0=y0,
        x0=x0,
        a=a,
        sigma=sigma,
        m=m,
        c=c,
    )

    return params



def kinematics_parameters(params: SpectraParams) -> SpectraParams:

    """
    Opens a GUI window to set kinematic fitting parameters
    and updates the params object with the selected values.
    """

    # Extract relevant parameters
    wave1_kin = params.wave1_kin
    wave2_kin = params.wave2_kin
    sigma_guess_kin = params.sigma_guess_kin
    redshift_guess_kin = params.redshift_guess_kin
    constant_resolution_lambda = params.constant_resolution_lambda
    resolution_kin = params.resolution_kin
    constant_resolution_r = params.constant_resolution_r
    resolution_kin_r = params.resolution_kin_r
    resolution_kin_muse = params.resolution_kin_muse
    ppxf_kin_preloaded_lib = params.ppxf_kin_preloaded_lib
    markers_ppxf_kin = params.markers_ppxf_kin
    stellar_library_kin = params.stellar_library_kin
    ppxf_kin_custom_lib = params.ppxf_kin_custom_lib
    ppxf_kin_lib_folder = params.ppxf_kin_lib_folder
    ppxf_kin_custom_temp_suffix = params.ppxf_kin_custom_temp_suffix
    ppxf_kin_generic_lib = params.ppxf_kin_generic_lib
    ppxf_kin_generic_lib_folder = params.ppxf_kin_generic_lib_folder
    ppxf_kin_FWHM_tem_generic = params.ppxf_kin_FWHM_tem_generic
    ppxf_kin_fixed_kin = params.ppxf_kin_fixed_kin
    no_gas_kin = params.no_gas_kin
    gas_kin = params.gas_kin
    ppxf_kin_mask_emission = params.ppxf_kin_mask_emission
    ppxf_kin_two_stellar_components = params.ppxf_kin_two_stellar_components
    ppxf_kin_age_model1 = params.ppxf_kin_age_model1
    ppxf_kin_met_model1 = params.ppxf_kin_met_model1
    ppxf_kin_vel_model1 = params.ppxf_kin_vel_model1
    ppxf_kin_sigma_model1 = params.ppxf_kin_sigma_model1
    ppxf_kin_age_model2 = params.ppxf_kin_age_model2
    ppxf_kin_met_model2 = params.ppxf_kin_met_model2
    ppxf_kin_vel_model2 = params.ppxf_kin_vel_model2
    ppxf_kin_sigma_model2 = params.ppxf_kin_sigma_model2
    kin_moments = params.kin_moments
    additive_degree_kin = params.additive_degree_kin
    multiplicative_degree_kin = params.multiplicative_degree_kin
    ppxf_kin_noise = params.ppxf_kin_noise
    kin_best_noise = params.kin_best_noise
    with_errors_kin = params.with_errors_kin
    ppxf_kin_have_user_mask = params.ppxf_kin_have_user_mask
    ppxf_kin_mask_ranges_str = params.ppxf_kin_mask_ranges_str
    ppxf_kin_mc_sim = params.ppxf_kin_mc_sim
    ppxf_kin_tie_balmer = params.ppxf_kin_tie_balmer
    ppxf_kin_dust_stars = params.ppxf_kin_dust_stars
    ppxf_kin_dust_gas = params.ppxf_kin_dust_gas
    ppxf_kin_save_spectra = params.ppxf_kin_save_spectra
    # prev_spec = params.prev_spec_nopath



    layout, scale_win, fontsize, default_size = misc.get_layout()
    sg.theme('LightBlue1')
    ppxf_kin_layout = [
        [sg.Text('Wavelength interval (A):', font = ('', default_size, 'bold')), sg.InputText(wave1_kin, size = (5,1), key = 'left_wave_ppxf_kin'), sg.Text('-'), sg.InputText(wave2_kin, size = (5,1), key = 'right_wave_ppxf_kin'), sg.Text('Sigma (km/s):',font = ('', default_size, 'bold')), sg.InputText(sigma_guess_kin, size = (5,1), key = 'sigma_guess_kin'), sg.Text('Redshift (z):',font = ('', default_size, 'bold')),sg.InputText(redshift_guess_kin, size = (7,1), key = 'redshift_guess_kin')],
        [sg.HorizontalSeparator()],
        [sg.Radio('Spec. res. FWHM (A):', "RADIORES1", default = constant_resolution_lambda, key = 'constant_resolution_lambda',tooltip='If the spectra have a constant resolution in terms of FWHM (in A)',font = ('', default_size, 'bold')), sg.InputText(resolution_kin , size = (4,1), key = 'ppxf_resolution'), sg.Radio('Spec. res. R:', "RADIORES1", key = 'constant_resolution_r', default = constant_resolution_r,font = ('', default_size, 'bold'), tooltip='If the spectra have a constant resolving power in terms of R'), sg.InputText(resolution_kin_r, size = (6,1), key = 'ppxf_resolution_r'), sg.Radio('Spec. res. MUSE data', "RADIORES1", key = 'resolution_kin_muse', default = resolution_kin_muse,font = ('', default_size, 'bold'), tooltip='If MUSE spectra, using a polynomial LSF')],
        [sg.HorizontalSeparator()],

        [sg.Radio('Preset SPS libraries included with SPAN:', 'RADIOLIBKIN', default = ppxf_kin_preloaded_lib, key = 'ppxf_kin_preloaded_lib', font = ('', default_size, 'bold')), sg.InputCombo(markers_ppxf_kin, key='markers_ppxf_kin',default_value=stellar_library_kin, readonly=True, size = (14,1))],                                                                                                                                                                                                            [sg.Radio('Custom (E)MILES:', 'RADIOLIBKIN', default = ppxf_kin_custom_lib, key = 'ppxf_kin_custom_lib', font = ('', default_size, 'bold'),tooltip='Select a folder containing your set of (E)MILES templates'), sg.InputText(ppxf_kin_lib_folder, size = (27,1), key = 'ppxf_kin_lib_folder'), sg.FolderBrowse(), sg.Text('Prefix:', font = ('', default_size, 'bold'), tooltip='Emiles templates have a suffix, please provide it'), sg.InputText(ppxf_kin_custom_temp_suffix, size = (12,1), key = 'ppxf_kin_custom_temp_suffix') ],
        [sg.Radio('Generic templates:', 'RADIOLIBKIN', default = ppxf_kin_generic_lib, key = 'ppxf_kin_generic_lib', font = ('', default_size, 'bold'),tooltip='Select a folder containing FITS files with linear sample and lamba in A'), sg.InputText(ppxf_kin_generic_lib_folder, size = (27,1), key = 'ppxf_kin_generic_lib_folder'), sg.FolderBrowse(), sg.Text('FWHM tem. (A):', font = ('', default_size, 'bold'), tooltip='Set the resolution of the templates'), sg.InputText(ppxf_kin_FWHM_tem_generic, size = (5,1), key = 'ppxf_kin_FWHM_tem_generic') ],
        [sg.HorizontalSeparator()],

        [sg.Radio('Fitting only stellar kinematics',"RADIOKIN", key = 'no_gas_kin', default = no_gas_kin, font = ('', default_size, 'bold'),tooltip='Considering only the kinematics of stars'), sg.Checkbox('Masking the gas emission lines', default = ppxf_kin_mask_emission, key = 'ppxf_kin_mask_emission', tooltip = ('Activate the masking if gas emission is present!'))],
        [sg.Text(''), sg.Checkbox('Fit two stellar components with the following parameters:', default = ppxf_kin_two_stellar_components, key = 'ppxf_kin_two_stellar_components', tooltip='Enable to fit TWO stellar SSPs with different V and sigma. Only for preset and EMILES SSP', font = ('', default_size, 'bold'))],
        [sg.Text(''), sg.Text('Select SSP model 1:'), sg.Text('Age(Gyr):'), sg.InputText(ppxf_kin_age_model1, size = (4,1), key = 'ppxf_kin_age_model1'), sg.Text('[M/H]:'), sg.InputText(ppxf_kin_met_model1, size = (4,1), key = 'ppxf_kin_met_model1'), sg.Text('Vel.(km/s):'),sg.InputText(ppxf_kin_vel_model1, size = (4,1), key = 'ppxf_kin_vel_model1'), sg.Text('Sigma(km/s):'), sg.InputText(ppxf_kin_sigma_model1, size = (4,1), key = 'ppxf_kin_sigma_model1')],

        [sg.Text(''), sg.Text('Select SSP model 2:'), sg.Text('Age(Gyr):'), sg.InputText(ppxf_kin_age_model2, size = (4,1), key = 'ppxf_kin_age_model2'), sg.Text('[M/H]:'), sg.InputText(ppxf_kin_met_model2, size = (4,1), key = 'ppxf_kin_met_model2'), sg.Text('Vel.(km/s):'),sg.InputText(ppxf_kin_vel_model2, size = (4,1), key = 'ppxf_kin_vel_model2'), sg.Text('Sigma(km/s):'), sg.InputText(ppxf_kin_sigma_model2, size = (4,1), key = 'ppxf_kin_sigma_model2')],
        [sg.Text('', font = ('', 1))],
        [sg.Radio('Fitting gas and stellar kinematics', "RADIOKIN", key = 'gas_kin', default = gas_kin, font = ('', default_size, 'bold'),tooltip='Fitting the kinematics of ONE stellar component and the gas emission'), sg.Checkbox('Fixing stellar kinematics first', default = ppxf_kin_fixed_kin, key = 'ppxf_kin_fixed_kin', tooltip = ('Perform a first fit with stellar kinematics and fix the moments for the gas fitting'))],

        [sg.HorizontalSeparator()],
        [sg.Checkbox('Correct for dust the stars', key = 'ppxf_kin_dust_stars', default = ppxf_kin_dust_stars, tooltip='Applying the default 2-params attenuation curve for stars of Cappellari 2023'), sg.Checkbox('Correct for dust the gas', key = 'ppxf_kin_dust_gas', default = ppxf_kin_dust_gas, tooltip='Applying the Calzetti extinction curve for gas'), sg.Checkbox('Tie Balmer lines', key = 'ppxf_kin_tie_balmer', default = ppxf_kin_tie_balmer)],
        [sg.Checkbox('Mask custom regions (A):', default=ppxf_kin_have_user_mask, key='ppxf_kin_have_user_mask', tooltip='Insert the regions you want to mask.'), sg.InputText(ppxf_kin_mask_ranges_str, size=(30, 1), key='ppxf_kin_mask_ranges'), sg.Button('Graphical masking', key='ppxf_kin_graphical_mask',button_color= ('black','light blue'), tooltip='Opens a plot window to draw the regions to mask')],
        [sg.HorizontalSeparator()],

        [sg.Text('Moments to fit:', font = ('', default_size, 'bold'), tooltip='Moments of the LOSVD. Minimum 2 (V and sigma), maximum 6. Proposed value = 4'), sg.InputText(kin_moments, size = (3,1), key = 'kin_moments'), sg.Text('Add. degree:', font = ('', default_size, 'bold'), tooltip='Additive degree to the fit. Deactivate (-1) if you are interested to gas flux'), sg.InputText(additive_degree_kin, size = (3,1), key = 'additive_degree_kin'), sg.Text('Mult. degree:', font = ('', default_size, 'bold'), tooltip='Multiplicative degree to the fit. Use if you want to use gas flux'), sg.InputText(multiplicative_degree_kin, size = (3,1), key = 'multiplicative_degree_kin'), sg.Text('Noise:', font = ('', default_size, 'bold'), tooltip='Mean noise per pixel of the spectrum'), sg.InputText(ppxf_kin_noise, size = (6,1), key = 'ppxf_kin_noise'), sg.Checkbox('Auto noise', default = kin_best_noise, key = ('kin_best_noise'), tooltip='Auto calculate the noise level of the spectrum for the best formal error estimation')],
        [sg.HorizontalSeparator()],

        [sg.Checkbox('Estimate the uncertainties with MonteCarlo simulations', font = ('', default_size, 'bold'), key = 'with_errors_kin', default = with_errors_kin, tooltip='Calculate the uncertainties in the LOSVD with MonteCarlo simulations'), sg.Text('N. simulations:'), sg.InputText(ppxf_kin_mc_sim, size = (7,1), key = 'ppxf_kin_mc_sim')],
        [sg.Checkbox('Save processed spectra', key = 'ppxf_kin_save_spectra', default = ppxf_kin_save_spectra, tooltip='For each spectrum, save the bestfit model, gas, emission, emission corrected, if available')],
        [sg.Button("Help", size=(12, 1),button_color=('black','orange')), sg.Push(), sg.Button('Confirm',button_color= ('white','black'), size = (18,1))]
        ]

    print ('*** Kinematics parameters window open. The main panel will be inactive until you close the window ***')
    ppxf_kin_window = sg.Window('pPXF Kinematics parameters', ppxf_kin_layout)


    while True:

        ppxf_kin_event, ppxf_kin_values = ppxf_kin_window.read()

        if ppxf_kin_event == sg.WIN_CLOSED:
            break

        constant_resolution_lambda = ppxf_kin_values['constant_resolution_lambda']
        constant_resolution_r = ppxf_kin_values['constant_resolution_r']
        resolution_kin_muse = ppxf_kin_values['resolution_kin_muse']
        stellar_library_kin = ppxf_kin_values['markers_ppxf_kin']
        ppxf_kin_fixed_kin = ppxf_kin_values['ppxf_kin_fixed_kin']

        gas_kin = ppxf_kin_values['gas_kin']
        no_gas_kin = ppxf_kin_values['no_gas_kin']
        kin_best_noise = ppxf_kin_values['kin_best_noise']
        with_errors_kin = ppxf_kin_values['with_errors_kin']
        ppxf_kin_save_spectra = ppxf_kin_values['ppxf_kin_save_spectra']

        #parameters for templates
        ppxf_kin_preloaded_lib = ppxf_kin_values['ppxf_kin_preloaded_lib']
        ppxf_kin_custom_lib = ppxf_kin_values['ppxf_kin_custom_lib']
        ppxf_kin_lib_folder = ppxf_kin_values['ppxf_kin_lib_folder']
        ppxf_kin_custom_temp_suffix = ppxf_kin_values['ppxf_kin_custom_temp_suffix']

        # Generic templates
        ppxf_kin_generic_lib = ppxf_kin_values['ppxf_kin_generic_lib']
        ppxf_kin_generic_lib_folder = ppxf_kin_values['ppxf_kin_generic_lib_folder']

        #checking if the resolution is meaningful
        if ppxf_kin_generic_lib:
            try:
                ppxf_kin_FWHM_tem_generic = float(ppxf_kin_values['ppxf_kin_FWHM_tem_generic'])
                if ppxf_kin_FWHM_tem_generic <=0:
                    sg.popup('Template resolution must be positive!')
                    ppxf_kin_FWHM_tem_generic = params.ppxf_kin_FWHM_tem_generic
                    continue
            except Exception:
                sg.popup('Invalid template resolution!')
                ppxf_kin_FWHM_tem_generic = params.ppxf_kin_FWHM_tem_generic
                continue

        #parameters for dust
        ppxf_kin_tie_balmer = ppxf_kin_values['ppxf_kin_tie_balmer']
        ppxf_kin_dust_stars = ppxf_kin_values['ppxf_kin_dust_stars']
        ppxf_kin_dust_gas = ppxf_kin_values['ppxf_kin_dust_gas']

        # Masking
        ppxf_kin_have_user_mask = ppxf_kin_values['ppxf_kin_have_user_mask']
        if ppxf_kin_have_user_mask:
            try:
                ppxf_kin_mask_ranges_str = ppxf_kin_values['ppxf_kin_mask_ranges']
                ppxf_kin_mask_ranges = eval(ppxf_kin_mask_ranges_str)
            except Exception:
                sg.Popup('Masking values not valid')
                ppxf_kin_mask_ranges_str = params.ppxf_kin_mask_ranges_str
                # ppxf_kin_mask_ranges = ppxf_kin_mask_ranges_default
                continue

        # Graphical masking
        if ppxf_kin_event == 'ppxf_kin_graphical_mask' and ppxf_kin_values['ppxf_kin_have_user_mask']:
            try:
                prev_spec_mask = ''
                prev_spec_mask = params.prev_spec
                touch_mode = (layout == layouts.layout_android)
                wavelength_mask, flux_mask, step_mask, name_mask = stm.read_spec(prev_spec_mask, params.lambda_units)
                current_mask_str = ppxf_kin_window['ppxf_kin_mask_ranges'].get()
                updated_mask = stm.graphical_masking_1D(wavelength_mask, flux_mask, current_mask_str, touch_mode=touch_mode)
                ppxf_kin_window['ppxf_kin_mask_ranges'].update(updated_mask)
            except Exception:
                sg.popup('You need to load at least one spectrum to perform the graphical masking!')
                continue

        #checking the existence of the custom templates in the specified folder
        if ppxf_kin_custom_lib:
            matching_temp_kin = glob.glob(os.path.join(ppxf_kin_lib_folder, ppxf_kin_custom_temp_suffix))
            if not matching_temp_kin:
                sg.Popup('Custom (E)miles templates not found. Check the suffix or the folder')
                continue

        ppxf_kin_mask_emission = ppxf_kin_values['ppxf_kin_mask_emission']
        ppxf_kin_two_stellar_components = ppxf_kin_values['ppxf_kin_two_stellar_components']

        # checking the compatibility with the generic template option
        if ppxf_kin_two_stellar_components and ppxf_kin_generic_lib:
            sg.popup('WARNING: fit with two stellar components is NOT compatible with the generic templates option')
            continue

        #check on the wavelength band
        try:
            wave1_kin = float(ppxf_kin_values['left_wave_ppxf_kin'])
            wave2_kin = float(ppxf_kin_values['right_wave_ppxf_kin'])
            resolution_kin = float(ppxf_kin_values['ppxf_resolution'])
            resolution_kin_r = float(ppxf_kin_values['ppxf_resolution_r'])
            sigma_guess_kin = float(ppxf_kin_values['sigma_guess_kin'])
            redshift_guess_kin = float(ppxf_kin_values['redshift_guess_kin'])
            additive_degree_kin = int(ppxf_kin_values['additive_degree_kin'])
            multiplicative_degree_kin = int(ppxf_kin_values['multiplicative_degree_kin'])

            kin_moments = int(ppxf_kin_values['kin_moments'])
            ppxf_kin_noise = float(ppxf_kin_values['ppxf_kin_noise'])
            ppxf_kin_mc_sim = int(ppxf_kin_values['ppxf_kin_mc_sim'])

            #parameters for two components fitting, checking only if activated
            if ppxf_kin_two_stellar_components:
                ppxf_kin_age_model1 = float(ppxf_kin_values['ppxf_kin_age_model1'])
                ppxf_kin_met_model1 = float(ppxf_kin_values['ppxf_kin_met_model1'])
                ppxf_kin_age_model2 = float(ppxf_kin_values['ppxf_kin_age_model2'])
                ppxf_kin_met_model2 = float(ppxf_kin_values['ppxf_kin_met_model2'])
                ppxf_kin_vel_model1 = float(ppxf_kin_values['ppxf_kin_vel_model1'])
                ppxf_kin_sigma_model1 = float(ppxf_kin_values['ppxf_kin_sigma_model1'])
                ppxf_kin_vel_model2 = float(ppxf_kin_values['ppxf_kin_vel_model2'])
                ppxf_kin_sigma_model2 = float(ppxf_kin_values['ppxf_kin_sigma_model2'])

        except ValueError:
            sg.popup ('Input parameters are not valid numbers!')
            continue

        if kin_moments < 2 or kin_moments > 6:
            sg.popup('Valid kinematics moments are between 2 and 6. Please, input a valid value')
            kin_moments = 4
            continue

        if ppxf_kin_event == 'Confirm':
            print ('Kinematics parameters confirmed. This main panel is now active again')
            print ('')
            break

        if ppxf_kin_event == 'Help':
            f = open(os.path.join(BASE_DIR, "help_files", "help_kinematics.txt"), 'r')
            file_contents = f.read()
            if layout == layouts.layout_android:
                sg.popup_scrolled(file_contents, size=(120, 30))
            else:
                sg.popup_scrolled(file_contents, size=(100, 40))

    ppxf_kin_window.close()

    params = replace(params,
            wave1_kin=wave1_kin,
            wave2_kin=wave2_kin,
            sigma_guess_kin=sigma_guess_kin,
            redshift_guess_kin=redshift_guess_kin,
            resolution_kin=resolution_kin,
            resolution_kin_r=resolution_kin_r,
            resolution_kin_muse=resolution_kin_muse,
            constant_resolution_lambda=constant_resolution_lambda,
            constant_resolution_r=constant_resolution_r,
            ppxf_kin_preloaded_lib=ppxf_kin_preloaded_lib,
            stellar_library_kin=stellar_library_kin,
            ppxf_kin_custom_lib=ppxf_kin_custom_lib,
            ppxf_kin_lib_folder=ppxf_kin_lib_folder,
            ppxf_kin_custom_temp_suffix=ppxf_kin_custom_temp_suffix,
            ppxf_kin_generic_lib = ppxf_kin_generic_lib,
            ppxf_kin_generic_lib_folder = ppxf_kin_generic_lib_folder,
            ppxf_kin_FWHM_tem_generic = ppxf_kin_FWHM_tem_generic,
            ppxf_kin_fixed_kin = ppxf_kin_fixed_kin,
            no_gas_kin=no_gas_kin, gas_kin=gas_kin,
            ppxf_kin_mask_emission=ppxf_kin_mask_emission,
            ppxf_kin_two_stellar_components=ppxf_kin_two_stellar_components,
            ppxf_kin_age_model1 = ppxf_kin_age_model1,
            ppxf_kin_met_model1 = ppxf_kin_met_model1,
            ppxf_kin_vel_model1 =  ppxf_kin_vel_model1,
            ppxf_kin_sigma_model1 = ppxf_kin_sigma_model1,
            ppxf_kin_age_model2 = ppxf_kin_age_model2,
            ppxf_kin_met_model2 = ppxf_kin_met_model2,
            ppxf_kin_vel_model2 = ppxf_kin_vel_model2,
            ppxf_kin_sigma_model2 = ppxf_kin_sigma_model2,
            kin_moments = kin_moments,
            additive_degree_kin = additive_degree_kin,
            multiplicative_degree_kin = multiplicative_degree_kin,
            ppxf_kin_noise = ppxf_kin_noise,
            kin_best_noise = kin_best_noise,
            with_errors_kin = with_errors_kin,
            ppxf_kin_have_user_mask = ppxf_kin_have_user_mask,
            ppxf_kin_mask_ranges_str = ppxf_kin_mask_ranges_str,
            ppxf_kin_mc_sim = ppxf_kin_mc_sim,
            ppxf_kin_tie_balmer = ppxf_kin_tie_balmer,
            ppxf_kin_dust_stars = ppxf_kin_dust_stars,
            ppxf_kin_dust_gas = ppxf_kin_dust_gas,
            ppxf_kin_save_spectra = ppxf_kin_save_spectra
            )

    return params



def population_parameters(params: SpectraParams) -> SpectraParams:

    """
    Opens a GUI window to set stellar population fitting parameters
    and updates the params object with the selected values.
    """

    # Extract relevant parameters
    wave1_pop = params.wave1_pop
    wave2_pop = params.wave2_pop
    res_pop = params.res_pop
    sigma_guess_pop = params.sigma_guess_pop
    z_pop = params.z_pop
    pop_with_gas = params.pop_with_gas
    pop_without_gas = params.pop_without_gas
    fit_components = params.fit_components
    ppxf_pop_dust_stars = params.ppxf_pop_dust_stars
    ppxf_pop_dust_gas = params.ppxf_pop_dust_gas
    ppxf_pop_tie_balmer = params.ppxf_pop_tie_balmer
    ppxf_pop_noise = params.ppxf_pop_noise
    regul_err = params.regul_err
    additive_degree = params.additive_degree
    multiplicative_degree = params.multiplicative_degree
    ppxf_best_noise_estimate = params.ppxf_best_noise_estimate
    ppxf_best_param = params.ppxf_best_param
    ppxf_frac_chi = params.ppxf_frac_chi
    ppxf_pop_convolve = params.ppxf_pop_convolve
    ppxf_pop_mask = params.ppxf_pop_mask
    ppxf_pop_want_to_mask = params.ppxf_pop_want_to_mask
    ppxf_pop_mask_ranges_str = params.ppxf_pop_mask_ranges_str
    ppxf_pop_lg_age = params.ppxf_pop_lg_age
    ppxf_pop_lg_met = params.ppxf_pop_lg_met
    with_errors = params.with_errors
    ppxf_min_age = params.ppxf_min_age
    ppxf_max_age = params.ppxf_max_age
    ppxf_min_met = params.ppxf_min_met
    ppxf_max_met = params.ppxf_max_met
    ppxf_pop_error_nsim = params.ppxf_pop_error_nsim
    ppxf_pop_preloaded_lib = params.ppxf_pop_preloaded_lib
    markers_ppxf = params.markers_ppxf
    stellar_library = params.stellar_library
    ppxf_pop_custom_lib = params.ppxf_pop_custom_lib
    ppxf_pop_lib_folder = params.ppxf_pop_lib_folder
    ppxf_custom_temp_suffix = params.ppxf_custom_temp_suffix
    ppxf_pop_custom_npz = params.ppxf_pop_custom_npz
    ppxf_pop_npz_file = params.ppxf_pop_npz_file
    stellar_parameters_lick_ppxf = params.stellar_parameters_lick_ppxf
    ssp_model_ppxf = params.ssp_model_ppxf
    interp_model_ppxf = params.interp_model_ppxf
    lick_ssp_models_ppxf = params.lick_ssp_models_ppxf
    interp_modes_ppxf = params.interp_modes_ppxf


    layout, scale_win, fontsize, default_size = misc.get_layout()
    sg.theme('LightBlue1')

    ppxf_pop_layout = [
        [sg.Text('Wavelength interval (A):', font = ('', default_size, 'bold')), sg.InputText(wave1_pop, size = (5,1), key = 'left_wave_ppxf_pop'), sg.Text('-'), sg.InputText(wave2_pop, size = (5,1), key = 'right_wave_ppxf_pop'), sg.Text('Spec. resolution FWHM (A):', font = ('', default_size, 'bold'),tooltip='Instrumental resolution (FWHM) in A of the spectral region'), sg.InputText(res_pop , size = (4,1), key = 'resolution_ppxf_pop')],
        [sg.Text('Velocity dispersion guess (km/s):', font = ('', default_size, 'bold')), sg.InputText(sigma_guess_pop, size = (9,1), key = 'sigma_guess_pop'), sg.Text('Redshift guess (z):', font = ('', default_size, 'bold')), sg.InputText(z_pop, size = (8,1), key = 'ppxf_z_pop')],
        [sg.HorizontalSeparator()],

        [sg.Radio('Fitting stars and gas together', "RADIOPOP", key = 'gas_pop', default = pop_with_gas, font = ('', default_size, 'bold')), sg.Radio('Fitting only stars',"RADIOPOP", key = 'no_gas_pop', default = pop_without_gas, font = ('', default_size, 'bold'))],
        [sg.Checkbox('Correct for dust the stars', key = 'ppxf_pop_dust_stars', default = ppxf_pop_dust_stars, tooltip='Applying the default 2-params attenuation curve for stars of Cappellari 2023'), sg.Checkbox('Correct for dust the gas', key = 'ppxf_pop_dust_gas', default = ppxf_pop_dust_gas, tooltip='Applying the Calzetti extinction curve for gas'), sg.Checkbox('Tie Balmer lines', key = 'ppxf_pop_tie_balmer', default = ppxf_pop_tie_balmer)],
        [sg.HorizontalSeparator()],

        [sg.Text('Noise:', font = ('', default_size, 'bold'), tooltip='Mean noise per pixel of the spectrum'), sg.InputText(ppxf_pop_noise, size = (8,1), key = 'ppxf_pop_noise'), sg.Text('Regul. error:', font = ('', default_size, 'bold'), tooltip='Regularization parameter. Higher values = sharper solution. Lower values = smoother fit. Set zero to deactivate.'), sg.InputText(regul_err, size = (5,1), key = 'regul_err'), sg.Text('Add. degree:', font = ('', default_size, 'bold'), tooltip='Additive degree to the fit. Not recommended to use for stellar populations!'), sg.InputText(additive_degree, size = (3,1), key = 'additive_degree'), sg.Text('Mult. degree:', font = ('', default_size, 'bold'), tooltip='Multiplicative degree to the fit. Usually 1 degree every 100 Angstron of spectrum to be fitted'), sg.InputText(multiplicative_degree, size = (3,1), key = 'multiplicative_degree')],
        [sg.Checkbox('Auto noise', default = ppxf_best_noise_estimate, key = 'ppxf_best_noise_estimate', tooltip='Let pPXF to estimate the best noise so that chi2 = 1 for non regularized fit'), sg.Checkbox('Auto noise and Regul. error', default = ppxf_best_param, key = ('ppxf_best_param'), tooltip='DO NOT use with > 200 templates or > 100 nn w. range. The Regul. error is a guess'), sg.Text('Fraction of Dchi2:', tooltip='Fraction of the desired Delta chi2 you want to reach. Max 20-30% if S/N < 40'), sg.Slider(range=(0.1, 1), orientation='h', default_value= ppxf_frac_chi, key='ppxf_frac_chi', resolution = 0.1, size = (13,20))],
        [sg.HorizontalSeparator()],

        [sg.Radio('Preset SPS libraries included with SPAN:', 'RADIOLIBPPXF', default = ppxf_pop_preloaded_lib, key = 'ppxf_pop_preloaded_lib', font = ('', default_size, 'bold'), tooltip='Use the sMILES library to measure also the Alpha/Fe'), sg.InputCombo(markers_ppxf, key='markers_ppxf',default_value=stellar_library, readonly=True, size = (14,1))],                                                                                                                                                                                                            [sg.Radio('Custom (E)MILES:', 'RADIOLIBPPXF', default = ppxf_pop_custom_lib, key = 'ppxf_pop_custom_lib', font = ('', default_size, 'bold'),tooltip='Select a folder containing your set of (E)MILES templates'), sg.InputText(ppxf_pop_lib_folder, size = (21,1), key = 'ppxf_pop_lib_folder'), sg.FolderBrowse(), sg.Text('Prefix:', font = ('', default_size, 'bold'), tooltip='Emiles templates have a suffix, please provide it'), sg.InputText(ppxf_custom_temp_suffix, size = (10,1), key = 'ppxf_custom_temp_suffix') ],
        [sg.Radio('Custom .npz template set:', 'RADIOLIBPPXF', default = ppxf_pop_custom_npz, key = 'ppxf_pop_custom_npz', font = ('', default_size, 'bold'), tooltip='Use your custom .npz template set'), sg.InputText(ppxf_pop_npz_file, size = (21,1), key = 'ppxf_pop_npz_file'), sg.FileBrowse()],

        [sg.Checkbox('Mask the emission lines', default = ppxf_pop_mask, key = 'ppxf_pop_mask', tooltip='If activated, you should perform a fit without gas emission'), sg.Checkbox('Regions to mask (A):', default = ppxf_pop_want_to_mask, key = 'ppxf_pop_want_to_mask'), sg.InputText(ppxf_pop_mask_ranges_str, size = (14,1), key = 'ppxf_pop_mask_ranges'), sg.Button('Graphical masking', key='ppxf_pop_graphical_mask',button_color= ('black','light blue'), tooltip='Opens a plot window to draw the regions to mask')],

        [sg.Checkbox('Convolve templates to galaxy res.', default = ppxf_pop_convolve, key = 'ppxf_pop_convolve', tooltip='If activated, the templates are convolved to the resolution of the galaxy spectrum'), sg.Checkbox('Mean ages in log10', default = ppxf_pop_lg_age, key = 'ppxf_pop_lg_age', tooltip='If de-activated, the mean ages are calculated in the linear grid (in Gyr) instead the default log10 grid of pPXF'), sg.Checkbox('Mean metal in log10', default = ppxf_pop_lg_met, key = 'ppxf_pop_lg_met', tooltip='If de-activated, the mean metallicities (including Alpha/Fe) are calculated in the linear grid instead the default log10 grid of pPXF') ],
        [sg.HorizontalSeparator()],

        [sg.Text('Age range for SPS models (Gyr):', font = ('', default_size, 'bold'))],
        [sg.Text('Min. age:'), sg.Slider(range=(0, 16), orientation='h', default_value= ppxf_min_age, key='ppxf_min_age', resolution = 0.1, size = (20,20)), sg.Text('Max. age:'), sg.Slider(range=(0, 16), orientation='h', default_value= ppxf_max_age, key='ppxf_max_age', resolution = 0.1, size = (20,20))],
        [sg.Text('Metallicity range for SPS models (dex):', font = ('', default_size, 'bold'))],
        [sg.Text('Min. met:'), sg.Slider(range=(-2.5, 0.8), orientation='h', default_value= ppxf_min_met, key='ppxf_min_met', resolution = 0.1, size = (20,20)), sg.Text('Max. met:'), sg.Slider(range=(-2, 0.8), orientation='h', default_value= ppxf_max_met, key='ppxf_max_met', resolution = 0.1, size = (20,20))],
        [sg.Checkbox('Estimate uncertainties for age and met (long process)', font = ('', default_size, 'bold'), key = 'ppxf_err_pop', default = with_errors,tooltip='Calculate the errors for age and metallicity with bootstrap simulations'), sg.Text('N. bootstrap:'), sg.InputText(ppxf_pop_error_nsim, size = (5,1), key = 'ppxf_pop_error_nsim')],
        [sg.HorizontalSeparator()],

        [sg.Checkbox('Lick/IDS analysis with SSP models:', default = stellar_parameters_lick_ppxf, key = 'stellar_parameters_lick_ppxf',tooltip='Use the pPXF results to estimate the stellar parameters also with the Lick/IDS indices', font = ('', default_size, 'bold')), sg.InputCombo(lick_ssp_models_ppxf, key='ssp_model_ppxf',default_value=ssp_model_ppxf, readonly=True, size = (11,1)), sg.Text('Interpolation:',tooltip='Interpolate linearly with griddata function or with machine learning Gaussian Process Regression (GPR)'), sg.InputCombo(interp_modes_ppxf, key='interp_model_ppxf',default_value=interp_model_ppxf, readonly=True, size = (7,1))],

        [sg.Button("Help", size=(12, 1),button_color=('black','orange')), sg.Push(), sg.Button('Confirm',button_color= ('white','black'), size = (18,1))]
        ]

    print ('*** Population parameters window open. The main panel will be inactive until you close the window ***')
    ppxf_pop_window = sg.Window('pPXF Population parameters', ppxf_pop_layout)

    while True:

        ppxf_pop_event, ppxf_pop_values = ppxf_pop_window.read()

        if ppxf_pop_event == sg.WIN_CLOSED:
            break

        ppxf_pop_want_to_mask = ppxf_pop_values['ppxf_pop_want_to_mask']
        ppxf_pop_dust_stars = ppxf_pop_values['ppxf_pop_dust_stars']
        ppxf_pop_dust_gas = ppxf_pop_values['ppxf_pop_dust_gas']
        ppxf_best_noise_estimate = ppxf_pop_values['ppxf_best_noise_estimate']
        ppxf_pop_lg_age = ppxf_pop_values['ppxf_pop_lg_age']
        ppxf_pop_lg_met = ppxf_pop_values['ppxf_pop_lg_met']
        stellar_parameters_lick_ppxf = ppxf_pop_values['stellar_parameters_lick_ppxf']
        ssp_model_ppxf = ppxf_pop_values['ssp_model_ppxf']
        interp_model_ppxf = ppxf_pop_values['interp_model_ppxf']

        if ppxf_pop_want_to_mask:
            try:
                ppxf_pop_mask_ranges_str = ppxf_pop_values['ppxf_pop_mask_ranges']
                ppxf_pop_mask_ranges = eval(ppxf_pop_mask_ranges_str)
            except Exception:
                sg.Popup('Masking values not valid')
                ppxf_pop_mask_ranges_str = params.ppxf_pop_mask_ranges_str
                # ppxf_pop_mask_ranges = ppxf_pop_mask_ranges_default
                continue

        # graphical masking
        if ppxf_pop_event == 'ppxf_pop_graphical_mask' and ppxf_pop_values['ppxf_pop_want_to_mask']:

            try:
                prev_spec_mask_pop = ''
                prev_spec_mask_pop = params.prev_spec
                touch_mode = (layout == layouts.layout_android)
                wavelength_mask_pop, flux_mask_pop, step_mask_pop, name_mask_pop = stm.read_spec(prev_spec_mask_pop, params.lambda_units)
                current_mask_str_pop = ppxf_pop_window['ppxf_pop_mask_ranges'].get()
                updated_mask_pop = stm.graphical_masking_1D(wavelength_mask_pop, flux_mask_pop, current_mask_str_pop, touch_mode=touch_mode)
                ppxf_pop_window['ppxf_pop_mask_ranges'].update(updated_mask_pop)
            except Exception:
                sg.popup('You need to load at least one spectrum to perform the graphical masking!')
                continue


        stellar_library = ppxf_pop_values['markers_ppxf']
        pop_with_gas = ppxf_pop_values['gas_pop']
        if pop_with_gas:
            fit_components = ('with_gas')
            ppxf_pop_tie_balmer = ppxf_pop_values['ppxf_pop_tie_balmer']
        pop_without_gas = ppxf_pop_values['no_gas_pop']
        if pop_without_gas:
            fit_components = ('without_gas')

        ppxf_pop_preloaded_lib = ppxf_pop_values['ppxf_pop_preloaded_lib']
        ppxf_pop_custom_lib = ppxf_pop_values['ppxf_pop_custom_lib']
        ppxf_pop_lib_folder = ppxf_pop_values['ppxf_pop_lib_folder']

        ppxf_pop_custom_npz = ppxf_pop_values['ppxf_pop_custom_npz']
        ppxf_pop_npz_file = ppxf_pop_values['ppxf_pop_npz_file']

        ppxf_pop_mask = ppxf_pop_values['ppxf_pop_mask']
        ppxf_best_param = ppxf_pop_values['ppxf_best_param']
        ppxf_custom_temp_suffix = ppxf_pop_values['ppxf_custom_temp_suffix']
        ppxf_frac_chi = float(ppxf_pop_values['ppxf_frac_chi'])
        ppxf_pop_convolve = ppxf_pop_values['ppxf_pop_convolve']
        with_errors = ppxf_pop_values['ppxf_err_pop']

        #checking the existence of the custom templates in the specified folder
        if ppxf_pop_custom_lib and not ppxf_pop_custom_npz:
            matching_temp = glob.glob(os.path.join(ppxf_pop_lib_folder, ppxf_custom_temp_suffix))
            if not matching_temp:
                sg.Popup('Custom (E)miles templates not found. Check the suffix or the folder')
                continue

        try:
            wave1_pop = float(ppxf_pop_values['left_wave_ppxf_pop'])
            wave2_pop = float(ppxf_pop_values['right_wave_ppxf_pop'])
            res_pop = float(ppxf_pop_values['resolution_ppxf_pop'])
            sigma_guess_pop = float(ppxf_pop_values['sigma_guess_pop'])
            regul_err = float(ppxf_pop_values['regul_err'])
            ppxf_pop_noise = float(ppxf_pop_values['ppxf_pop_noise'])
            additive_degree = int(ppxf_pop_values['additive_degree'])
            multiplicative_degree = int(ppxf_pop_values['multiplicative_degree'])
            z_pop = float(ppxf_pop_values['ppxf_z_pop'])
            ppxf_min_age = float(ppxf_pop_values['ppxf_min_age'])
            ppxf_max_age = float(ppxf_pop_values['ppxf_max_age'])
            ppxf_min_met = float(ppxf_pop_values['ppxf_min_met'])
            ppxf_max_met = float(ppxf_pop_values['ppxf_max_met'])
            age_range_array = np.array([ppxf_min_age, ppxf_max_age])
            met_range_array = np.array([ppxf_min_met, ppxf_max_met])
            ppxf_pop_error_nsim = int(ppxf_pop_values['ppxf_pop_error_nsim'])

        except ValueError:
            sg.popup ('Invalid input parameters!')
            continue

        # consistency check on the input parameters
        if regul_err < 0:
            sg.popup ('Regularization error must be greater than zero. Set zero for non-regularized fit')
            regul_err = 0.02 #If I close the window restore the default value
            continue
        if ppxf_pop_noise <=0:
            sg.popup ('Noise must be positive')
            ppxf_pop_noise = 0.0163
            continue

        if with_errors and ppxf_pop_error_nsim < 3:
            sg.popup('Well, the number of bootstrap simulations cannot be zero, negative or too close to zero! Select at least 3 simulations to perform')
            ppxf_pop_error_nsim = 50
            continue

        if with_errors and ppxf_pop_error_nsim > 100:
            sg.popup('The numnber of bootstrap simulations you entered is quite large. This will require a lot of time (hours). Are you sure?')

        if ppxf_pop_event == 'Confirm':
            print ('Population parameters confirmed. This main panel is now active again')
            print ('')
            break

        if ppxf_pop_event == 'Help':
            f = open(os.path.join(BASE_DIR, "help_files", "help_stellar_pop.txt"), 'r')
            file_contents = f.read()
            if layout == layouts.layout_android:
                sg.popup_scrolled(file_contents, size=(120, 30))
            else:
                sg.popup_scrolled(file_contents, size=(100, 40))

    ppxf_pop_window.close()


    params = replace(params,
                wave1_pop = wave1_pop,
                wave2_pop = wave2_pop,
                res_pop = res_pop,
                sigma_guess_pop = sigma_guess_pop,
                z_pop = z_pop,
                pop_with_gas = pop_with_gas,
                pop_without_gas = pop_without_gas,
                fit_components = fit_components,
                ppxf_pop_dust_stars = ppxf_pop_dust_stars,
                ppxf_pop_dust_gas = ppxf_pop_dust_gas,
                ppxf_pop_tie_balmer = ppxf_pop_tie_balmer,
                ppxf_pop_noise = ppxf_pop_noise,
                regul_err = regul_err,
                additive_degree = additive_degree,
                multiplicative_degree = multiplicative_degree,
                ppxf_best_noise_estimate = ppxf_best_noise_estimate,
                ppxf_best_param = ppxf_best_param,
                ppxf_frac_chi = ppxf_frac_chi,
                ppxf_pop_convolve = ppxf_pop_convolve,
                ppxf_pop_mask = ppxf_pop_mask,
                ppxf_pop_want_to_mask = ppxf_pop_want_to_mask,
                ppxf_pop_mask_ranges_str = ppxf_pop_mask_ranges_str,
                ppxf_pop_lg_age = ppxf_pop_lg_age,
                ppxf_pop_lg_met = ppxf_pop_lg_met,
                with_errors = with_errors,
                ppxf_min_age = ppxf_min_age,
                ppxf_max_age = ppxf_max_age,
                ppxf_min_met = ppxf_min_met,
                ppxf_max_met = ppxf_max_met,
                ppxf_pop_error_nsim = ppxf_pop_error_nsim,
                ppxf_pop_preloaded_lib = ppxf_pop_preloaded_lib,
                markers_ppxf = markers_ppxf,
                stellar_library = stellar_library,
                ppxf_pop_custom_lib = ppxf_pop_custom_lib,
                ppxf_pop_lib_folder = ppxf_pop_lib_folder,
                ppxf_custom_temp_suffix = ppxf_custom_temp_suffix,
                ppxf_pop_custom_npz = ppxf_pop_custom_npz,
                ppxf_pop_npz_file = ppxf_pop_npz_file,
                stellar_parameters_lick_ppxf = stellar_parameters_lick_ppxf,
                ssp_model_ppxf = ssp_model_ppxf,
                interp_model_ppxf = interp_model_ppxf,
                lick_ssp_models_ppxf = lick_ssp_models_ppxf,
                interp_modes_ppxf = interp_modes_ppxf
                )

    return params
