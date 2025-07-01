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

# Function to build up the Spectral manipulation panel and the relative parameter windows.
# Returns all the parameters handled by the Spectral manipulation panel and modified by the user in the GUI.

try: #try local import if executed as script
    #GUI import
    from FreeSimpleGUI_local import FreeSimpleGUI as sg
    from params import SpectraParams
    from span_modules import misc
    from span_modules import layouts

except ModuleNotFoundError: #local import if executed as package
    #GUI import
    from span.FreeSimpleGUI_local import FreeSimpleGUI as sg
    from . import misc
    from . import layouts
    from .params import SpectraParams

#python imports
import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation
import datetime
import os
from dataclasses import replace

# from params import SpectraParams

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(CURRENT_DIR)



def spectra_manipulation(params: SpectraParams) -> SpectraParams:

    # assigning all the params to local variables and work with them.
    current_order = params.current_order
    reorder_op = params.reorder_op
    reordered_operations = params.reordered_operations
    active_operations = params.active_operations

    save_intermediate_spectra = params.save_intermediate_spectra
    save_final_spectra = params.save_final_spectra
    not_save_spectra = params.not_save_spectra

    #DYNAMIC CLEANING PARAMETERS
    clip_factor = params.clip_factor
    sigma_clip_resolution = params.sigma_clip_resolution
    sigma_clip_single_vel = params.sigma_clip_single_vel
    sigma_clip_single_value = params.sigma_clip_single_value
    sigma_clip_have_file = params.sigma_clip_have_file
    sigma_clip_sigma_file = params.sigma_clip_sigma_file

    #DENOISING PARAMETERS
    moving_average = params.moving_average
    box_moving_avg = params.box_moving_avg
    box_moving_avg_size = params.box_moving_avg_size
    gauss_moving_avg = params.gauss_moving_avg
    gauss_moving_avg_kernel = params.gauss_moving_avg_kernel
    low_pass_filter = params.low_pass_filter
    lowpass_cut_off = params.lowpass_cut_off
    lowpass_order = params.lowpass_order
    bandpass_filter = params.bandpass_filter
    bandpass_lower_cut_off = params.bandpass_lower_cut_off
    bandpass_upper_cut_off = params.bandpass_upper_cut_off
    bandpass_order = params.bandpass_order

    #DOPCOR DEFAULT PARAMETERS
    dop_cor_single_shot_vel = params.dop_cor_single_shot_vel
    dop_cor_have_file = params.dop_cor_have_file
    dop_cor_file = params.dop_cor_file
    dop_cor_single_shot = params.dop_cor_single_shot
    dop_cor_have_vel = params.dop_cor_have_vel
    dop_cor_have_z = params.dop_cor_have_z

    #HELIOCOR DEFAULT PARAMETERS
    helio_have_file = params.helio_have_file
    helio_file = params.helio_file
    helio_single_shot = params.helio_single_shot
    helio_single_shot_location = params.helio_single_shot_location
    helio_single_shot_date = params.helio_single_shot_date
    ra_obj = params.ra_obj
    dec_obj = params.dec_obj

    #DEGRADE RESOLUTION PARAMETERS
    is_initial_res_r = params.is_initial_res_r
    initial_res_r = params.initial_res_r
    res_degrade_to_r = params.res_degrade_to_r
    final_res_r = params.final_res_r
    res_degrade_to_fwhm = params.res_degrade_to_fwhm
    final_res_r_to_fwhm = params.final_res_r_to_fwhm
    is_initial_res_fwhm = params.is_initial_res_fwhm
    initial_res_fwhm = params.initial_res_fwhm
    final_res_fwhm = params.final_res_fwhm
    res_degrade_muse = params.res_degrade_muse
    res_degrade_muse_value = params.res_degrade_muse_value

    #CONTINUUM SUBTRACTION PARAMETERS
    markers_cont_operations = params.markers_cont_operations
    cont_math_operation = params.cont_math_operation
    cont_model_filtering = params.cont_model_filtering
    cont_model_poly = params.cont_model_poly
    cont_want_to_mask = params.cont_want_to_mask
    cont_mask_ranges_str = params.cont_mask_ranges_str
    cont_mask_ranges = params.cont_mask_ranges
    cont_poly_degree = params.cont_poly_degree


    cropping_spectrum = params.cropping_spectrum
    cropping_low_wave = params.cropping_low_wave
    cropping_high_wave = params.cropping_high_wave
    sigma_clipping = params.sigma_clipping
    wavelet_cleaning = params.wavelet_cleaning
    sigma_wavelets = params.sigma_wavelets
    wavelets_layers = params.wavelets_layers
    filter_denoise = params.filter_denoise
    dop_cor = params.dop_cor
    helio_corr = params.helio_corr

    # Spectra processing frame default parameters
    rebinning = params.rebinning
    rebinning_log = params.rebinning_log
    rebinning_linear = params.rebinning_linear
    rebin_step_pix = params.rebin_step_pix
    rebin_step_sigma = params.rebin_step_sigma
    degrade = params.degrade
    normalize_wave = params.normalize_wave
    norm_wave = params.norm_wave
    sigma_broad = params.sigma_broad
    sigma_to_add = params.sigma_to_add
    add_noise = params.add_noise
    noise_to_add = params.noise_to_add
    continuum_sub = params.continuum_sub

    # Math frame default parameters
    average_all = params.average_all
    norm_and_average = params.norm_and_average
    do_nothing = params.do_nothing
    sum_all = params.sum_all
    normalize_and_sum_all = params.normalize_and_sum_all
    use_for_spec_an = params.use_for_spec_an
    subtract_normalized_avg = params.subtract_normalized_avg
    subtract_normalized_spec = params.subtract_normalized_spec
    spectra_to_subtract = params.spectra_to_subtract
    add_pedestal = params.add_pedestal
    pedestal_to_add = params.pedestal_to_add
    multiply = params.multiply
    multiply_factor = params.multiply_factor
    derivatives = params.derivatives

    #variables to prevent the program to crash in case some of the list spectra loaded are not valid
    spectra_number = params.spectra_number
    fatal_condition = params.fatal_condition

    layout, scale_win, fontsize, default_size = misc.get_layout()

    if layout == layouts.layout_macos:
        default_font_size = 14
    else:
        default_font_size = 10

    sg.theme('DarkBlue3')
    spec_layout = [



    #Spectra pre-processing
    [sg.Frame('Spectra pre-processing', [
    [sg.Checkbox('Cropping', key ='cropping', font = ('Helvetica', default_font_size, 'bold'), default = cropping_spectrum, tooltip='Crop the spectrum to a user defined wavelength range'), sg.Text('Lower wave'), sg.InputText(cropping_low_wave, key = 'cropping_low_wave', size = (5,1)), sg.Text('Upper wave'), sg.InputText(cropping_high_wave, key = 'cropping_high_wave', size = (5,1))],

    [sg.Checkbox('Dynamic cleaning', font = ('Helvetica', default_font_size, 'bold'), key = 'sigma_clip', default = sigma_clipping,tooltip='Perform sigma clipping to erase spikes'), sg.Push(), sg.Button('Clean parameters',button_color= ('black','light blue'), size = (23,1))],

    [sg.Checkbox('Wavelet cleaning', font = ('Helvetica', default_font_size, 'bold'), key = 'wavelet_cleaning', default = wavelet_cleaning,tooltip='Perform wavelet cleaning of the spectrum'), sg.Text('sigma:'),sg.InputText(sigma_wavelets, key = 'sigma_wavelets', size = (4,1)), sg.Text('Wavelet layers:'), sg.InputText(wavelets_layers, key = 'wavelets_layers', size = (3,1))],

    [sg.Checkbox('Filtering and denoising', font = ('Helvetica', default_font_size, 'bold'), key = 'filter_denoise', default = filter_denoise,tooltip='Filters to smooth the spectrum'), sg.Push(), sg.Button('Denoise parameters',button_color= ('black','light blue'), size = (23,1))],

    [sg.Checkbox('Doppler/z correction', font = ('Helvetica', default_font_size, 'bold'), key = 'dopcor', default = dop_cor,tooltip='Doppler and redshift correction of spectrum, from a list file or from a fixed radial velocity or z value'), sg.Push(), sg.Button('Dopcor parameters',button_color= ('black','light blue'), size = (23,1))],

    [sg.Checkbox('Heliocentric correction', font = ('Helvetica', default_font_size, 'bold'), key = 'helio_corr', default = helio_corr,tooltip='Heliocentric correction, from a formatted file or by inserting the location, time and object coordinates (RA and Dec) manually'), sg.Push(), sg.Button('Heliocor parameters',button_color= ('black','light blue'), size = (23,1))],

    ], font=("Helvetica", 12, 'bold'), title_color = 'lightgreen'),

    #2) spectra processing
    sg.Frame('Spectra processing', [
    [sg.Checkbox('Rebin', font = ('Helvetica', default_font_size, 'bold'), key = 'rebin', default = rebinning,tooltip='Rebinning the spectrum, to a linear wavelength step (A) or to a linear sigma step (km/s)'), sg.Radio('pix lin.', "RADIO1", default=rebinning_linear, key = 'rebin_pix_lin'), sg.InputText(rebin_step_pix, size = (4,1), key = 'rebin_step_pix'), sg.Radio('sigma lin.', "RADIO1", default = rebinning_log, key = 'rebin_sigma_lin'), sg.InputText(rebin_step_sigma, size = (3,1), key = 'rebin_step_sigma')],
    [sg.Checkbox('Degrade resolution', font = ('Helvetica', default_font_size, 'bold'), key = 'degrade_resolution', default = degrade,tooltip='Degrade resolution to a user defined value'), sg.Push(), sg.Button('Degrade parameters',button_color= ('black','light blue'), size = (20,1))],

    [sg.Checkbox('Normalise spectrum to:', font = ('Helvetica', default_font_size, 'bold'), key = 'norm_spec', default = normalize_wave,tooltip='Normalise the flux to a user defined wavelength'), sg.InputText(norm_wave, size = (6,1), key = 'norm_wave'), sg.Text('A')],

    [sg.Checkbox('Sigma broadening', font = ('Helvetica', default_font_size, 'bold'), key = 'broadening_spec', default = sigma_broad,tooltip='Broad the spectrum by adding a user defined sigma (km/s). This will NOT be the total sigma broadening of your spectrum!'), sg.Text('Add sigma (km/s): ', font = ('Helvetica', default_font_size)), sg.InputText(sigma_to_add, size = (4,1), key = 'sigma_to_add')],
    [sg.Checkbox('Add noise', font = ('Helvetica', default_font_size, 'bold'), key = 'add_noise', default = add_noise,tooltip='Adding poissonian noise to the spectrum'), sg.Text('Signal to Noise (S/N) to add:'), sg.InputText(noise_to_add, size = (5,1), key = 'noise_to_add')],

    [sg.Checkbox('Continuum modelling', font = ('Helvetica', default_font_size, 'bold'), key = 'cont_sub', default = continuum_sub,tooltip='Perform the continuum estimation to subtract or divide to the spectrum'), sg.Push(), sg.Button('Continuum parameters',button_color= ('black','light blue'), size = (20,1))],
    [sg.Text('', font = ("Helvetica", 1))],

    ], font=("Helvetica", 12, 'bold'),title_color = 'lightgreen'),

    #3) spectra math
    sg.Frame('Spectra math', [
    [sg.Checkbox('Subtract normalised average', font = ('Helvetica', default_font_size, 'bold'), key = 'subtract_norm_avg', default = subtract_normalized_avg,tooltip='Normalise and subtract to the selected spectrum the normalised average of all the spectra')],
    [sg.Checkbox('Subtract norm. spec.', font = ('Helvetica', default_font_size, 'bold'), key = 'subtract_norm_spec', default = subtract_normalized_spec,tooltip='Normalise and subtract to the selected spectrum a user selected spectrum'), sg.InputText(spectra_to_subtract, size = (17,1), key = 'spec_to_subtract'), sg.FileBrowse(tooltip='Load a spectrum (ASCII or fits) to be normalised and subtracted')],
    [sg.Checkbox('Add constant', font = ('Helvetica', default_font_size, 'bold'), key = 'add_pedestal', default = add_pedestal,tooltip='Simply add a constant value to the spectrum'), sg.InputText(pedestal_to_add, size = (7,1), key = 'pedestal_to_add'), sg.Checkbox('Multiply by:', font = ('Helvetica', default_font_size, 'bold'), key = 'multiply', default = multiply,tooltip='Multiply the spectrum by a constant'), sg.InputText(multiply_factor , size = (7,1), key = 'multiply_factor')],
    [sg.Checkbox('Calculate first and second derivatives', default = derivatives, key = 'derivatives', font = ('Helvetica', default_font_size, 'bold'),tooltip='Calculate the derivative of the spectra')],
    [sg.HorizontalSeparator()],
    [sg.Radio('Average all', "RADIOMATH", key = 'avg_all', default = average_all,tooltip='Average all the loaded spectra'), sg.Radio('Norm. and average all', "RADIOMATH", key = 'norm_avg_all', default = norm_and_average,tooltip='First normalise, then average all the loaded spectra'), sg.Radio('Nothing', "RADIOMATH", key = 'none', default = do_nothing,tooltip='Select this option if you DO NOT want to combine the spectra', font = ('Helvetica', default_font_size, 'bold'))],
    [sg.Radio('Sum all', "RADIOMATH", key = 'sum_all', default = sum_all,tooltip='Sum all the loaded spectra'), sg.Radio('Norm. and sum all', "RADIOMATH", key = 'norm_sum_all', default = normalize_and_sum_all,tooltip='First normalise, then sum all the loaded spectra'), sg.Checkbox('Use for spec. an.', text_color = 'yellow', key = 'use_for_spec_an', default = use_for_spec_an,tooltip='Select this to use the combined spectrum for the spectral analysis', font = ('Helvetica', default_font_size, 'bold'))],
    ],font=("Helvetica", 12, 'bold'))],

    # Bottom parameters
    [sg.Checkbox('Reorder', key = 'reorder_op', default = reorder_op, tooltip='Activate in case you want to perform the spectra manipulation tasks in different order'), sg.Button('Reorder tasks', tooltip='Change the order of the spectra manipulation tasks'), sg.Radio('Save intermediate spectra', "RADIOSAVE", key = 'save_intermediate_spectra', default = save_intermediate_spectra, tooltip='Save a processed spectrum for EACH activated task'), sg.Radio('Save final spectra', "RADIOSAVE", key = 'save_final_spectra', default = save_final_spectra, tooltip='Save only the final processed spectrum after applying the tasks'), sg.Radio('Do not save processed spectra', "RADIOSAVE", key = 'not_save_spectra', default = not_save_spectra, tooltip='Do not save any processed spectrum to disc'), sg.Push(), sg.Button('I need help',button_color=('black','orange'), size = (11,1)), sg.Button('Confirm',button_color= ('white','black'), size = (18,1))]
    ]

    spec_window = sg.Window('Spectra manipulation parameters', spec_layout)


    while True:
        spec_event, spec_values = spec_window.read()

        if spec_event == sg.WIN_CLOSED:
            break

        # Assigning parameters from the GUI to local variables
        reorder_op = spec_values['reorder_op']
        save_intermediate_spectra = spec_values['save_intermediate_spectra']
        save_final_spectra = spec_values['save_final_spectra']
        not_save_spectra = spec_values['not_save_spectra']

        # Spectra pre-processing
        cropping_spectrum = spec_values['cropping']
        sigma_clipping = spec_values['sigma_clip']
        wavelet_cleaning = spec_values['wavelet_cleaning']
        filter_denoise = spec_values['filter_denoise']
        dop_cor = spec_values['dopcor']
        helio_corr = spec_values['helio_corr']

        #spectra processing
        rebinning = spec_values['rebin']
        rebinning_log = spec_values['rebin_sigma_lin']
        rebinning_linear = spec_values['rebin_pix_lin']
        degrade = spec_values['degrade_resolution']
        normalize_wave = spec_values['norm_spec']
        sigma_broad = spec_values['broadening_spec']
        add_noise = spec_values['add_noise']
        continuum_sub = spec_values['cont_sub']

        #math parameters
        average_all = spec_values['avg_all']
        norm_and_average = spec_values['norm_avg_all']
        do_nothing = spec_values['none']
        sum_all = spec_values['sum_all']
        normalize_and_sum_all = spec_values['norm_sum_all']
        use_for_spec_an = spec_values['use_for_spec_an']
        subtract_normalized_avg = spec_values['subtract_norm_avg']
        subtract_normalized_spec = spec_values['subtract_norm_spec']
        spectra_to_subtract = spec_values['spec_to_subtract']
        add_pedestal = spec_values['add_pedestal']
        multiply = spec_values['multiply']
        derivatives = spec_values['derivatives']


    #********** Initializing and checking the variables of the Spectra pre-processing frame **********
        # Creating the dictionary of the spectra manipulation tasks to be user in case of re-ordering for Android devices.
        available_operations = [
            ("Cropping", "cropping_spectrum", cropping_spectrum),
            ("Dynamic cleaning", "sigma_clipping", sigma_clipping),
            ("Wavelet cleaning", "wavelet_cleaning", wavelet_cleaning),
            ("Filtering and denoise", "filter_denoise", filter_denoise),
            ("Doppler/z correction", "dop_cor", dop_cor),
            ("Heliocentric correction", "helio_corr", helio_corr),
            ("Rebinning", "rebining", rebinning),
            ("Degrade resolution", "degrade", degrade),
            ("Normalise spectrum", "normalize_wave", normalize_wave),
            ("Velocity dispersion broadening", "sigma_broad", sigma_broad),
            ("Add noise", "add_noise", add_noise),
            ("Continuun modelling", "continuum_sub", continuum_sub),

            ("Subtract normalised average", "subtract_normalized_avg",subtract_normalized_avg),
            ("Subtract norm. spec","subtract_normalized_spec",subtract_normalized_spec),
            ("Add constant", "add_pedestal", add_pedestal),
            ("Multiply by", "multiply", multiply),
            ("Calculate first and second derivatives", "derivatives",derivatives),
        ]

        # Select only the tasks activated
        active_operations = [op[:2] for op in available_operations if op[2]]

        # Check if a previous current_order exists and match the activated tasks
        if current_order is None or set([op[1] for op in current_order]) != set([op[1] for op in active_operations]):
            current_order = active_operations.copy()

        # Making a copy
        reordered_operations = current_order.copy()

        if spec_event == 'Reorder tasks':

            # reordering windows layout
            layout_reorder = [
                [sg.Text("Please, re-order the tasks:")],
                [sg.Listbox([op[0] for op in reordered_operations], size=(40, 12), key="-OP_LIST-", select_mode=sg.LISTBOX_SELECT_MODE_SINGLE)],
                [sg.Button("Move up"), sg.Button("Move down"), sg.Button("Confirm"), sg.Button("Cancel")],
            ]

            # creating thw window
            window_reorder = sg.Window("Order the tasks", layout_reorder)

            sorting_cond = 0
            while True:
                event_reorder, values_reorder = window_reorder.read()

                if event_reorder == sg.WINDOW_CLOSED or event_reorder == "Cancel":
                    reordered_operations = current_order.copy()
                    break

                elif event_reorder == "Move up":
                    selected = values_reorder["-OP_LIST-"]
                    if selected:
                        idx = [op[0] for op in reordered_operations].index(selected[0])
                        if idx > 0:
                            #change the order of the tasks
                            reordered_operations[idx], reordered_operations[idx - 1] = reordered_operations[idx - 1], reordered_operations[idx]
                            #update the Listbox
                            window_reorder["-OP_LIST-"].update([op[0] for op in reordered_operations])

                elif event_reorder == "Move down":
                    selected = values_reorder["-OP_LIST-"]
                    if selected:
                        idx = [op[0] for op in reordered_operations].index(selected[0])
                        if idx < len(reordered_operations) - 1:
                            #change the order of the tasks
                            reordered_operations[idx], reordered_operations[idx + 1] = reordered_operations[idx + 1], reordered_operations[idx]
                            #update the Listbox
                            window_reorder["-OP_LIST-"].update([op[0] for op in reordered_operations])

                elif event_reorder == "Confirm":
                    sorting_cond = 1
                    #Save the new task order
                    current_order = reordered_operations.copy()
                    active_operations = current_order.copy()

                    #Activate the reorder checkbox automatically
                    spec_window['reorder_op'].update(True)
                    break

            window_reorder.close()

            if sorting_cond == 1 and len(reordered_operations) != 0:
                print("Ordered tasks:")
                for op in reordered_operations:
                    print(op[0])
            if sorting_cond == 1 and len(reordered_operations) == 0:
                print('No active tasks')


        #Initializing the values and parameters of the spectra manipulation panel

        #1) CROPPING PARAMETERS
        if cropping_spectrum:
            try:
                cropping_low_wave = float(spec_values['cropping_low_wave'])
                cropping_high_wave = float(spec_values['cropping_high_wave'])
            except Exception:
                sg.popup ('Cropping parameters not valid')
                continue


        #2) DYNAMIC CLEANING PARAMETERS
        if spec_event  == ('Clean parameters'):

            sg.theme('LightBlue1')
            clean_layout = [
                [sg.Text('Sigma to clip:',tooltip='Clipping factor'), sg.InputText(clip_factor, size = (3,1), key = 'clip_factor'), sg.Text('Res. (R)',tooltip='Spectrum resolution'), sg.InputText(sigma_clip_resolution, size = (5,1), key = 'res_spec_for_sigma_clip'), sg.Radio('Velocity dispersion (km/s)', "RADIOCLIP", default = sigma_clip_single_vel, key = 'single_vel_clip', tooltip='Velocity dispersion'), sg.InputText(sigma_clip_single_value, size = (4,1), key = 'clip_to_vel')],
                [sg.Radio('R and sigma vel file', "RADIOCLIP", default = sigma_clip_have_file, key = 'file_for_clip',tooltip='ASCII file with R and sigma to perform sigma clipping for all the loaded spectra'), sg.InputText(sigma_clip_sigma_file, size=(14, 1), key = 'sigma_clip_file'), sg.FileBrowse(tooltip='Load an ASCII file containing: Name of the spectrum, Resolution (R), sigma (km/s), in the same order of the original spectra list')],
                [sg.Push(), sg.Button('Confirm',button_color= ('white','black'), size = (18,1))]
                ]

            print ('*** Clean spectra window open. The main panel will be inactive until you close the window ***')
            clean_window = sg.Window('Clean spectra parameters', clean_layout)

            while True:
                clean_event, clean_values = clean_window.read()

                if clean_event == sg.WIN_CLOSED:
                    break

                try:
                    clip_factor = float(clean_values['clip_factor'])
                    sigma_clip_single_value = float(clean_values['clip_to_vel'])
                    sigma_clip_resolution = int(clean_values['res_spec_for_sigma_clip'])
                    if clip_factor <=0:
                        sg.popup('Invalid sigma clip factor. Must be > 0!')
                        continue
                except ValueError:
                    sg.popup('Sigma clip factor is not a number!')
                    continue

                sigma_clip_single_vel = clean_values['single_vel_clip']
                sigma_clip_have_file = clean_values['file_for_clip']
                sigma_clip_sigma_file = clean_values['sigma_clip_file']

                if clean_event == 'Confirm':
                    print ('Clean parameters confirmed. This main panel is now active again')
                    print ('')
                    break

            clean_window.close()


        #3) WAVELET PARAMETERS
        if wavelet_cleaning:

            try:
                sigma_wavelets = float(spec_values['sigma_wavelets'])
                wavelets_layers = int(spec_values['wavelets_layers'])
            except Exception:
                sg.popup('Wavelet parameters not valid')
                continue
            if sigma_wavelets <= 0 or wavelets_layers <= 0:
                sg.Popup ('Wavelet parameters must be greater than zero!')
                continue
            if wavelets_layers > 20:
                sg.Popup ('Wavelet layers must be smaller than 20. Try again')
                continue



        #4) FILTERING AND DENOISINS PARAMETERS
        if spec_event == 'Denoise parameters':

            sg.theme('LightBlue1')
            denoise_layout = [
            [sg.Checkbox('Moving average:', font = ('Helvetiva', 11, 'bold'), default = moving_average, key = 'moving_average'), sg.Radio('Simple box:', "MOVAVG", default = box_moving_avg, key = 'box_moving_avg', font = ('Helvetica', 11, 'bold')), sg.Text('Box size (pix):'), sg.InputText(box_moving_avg_size, key = 'box_moving_avg_size', size = (7,1)),sg.Radio('Gaussian kernel:', "MOVAVG", default = gauss_moving_avg, key = 'gauss_moving_avg'), sg.Text('Sigma kernel (pix):'), sg.InputText(gauss_moving_avg_kernel, key = 'gauss_moving_avg_kernel', size = (5,1))],

            [sg.HorizontalSeparator()],
            [sg.Checkbox('Low-pass filter (Butterworth)', font = ('Helvetiva', 11, 'bold'), default = low_pass_filter, key = 'low_pass_filter'), sg.Text('Cut-off:'), sg.InputText(lowpass_cut_off, key = 'lowpass_cut_off', size = (7,1)), sg.Text('Filter order:'), sg.InputText(lowpass_order, key = 'lowpass_order', size = (7,1))],
            [sg.Checkbox('Band-pass filter (Butterworth)', font = ('Helvetiva', 11, 'bold'), default = bandpass_filter, key = 'bandpass_filter'), sg.Text('lower Cut-off:'), sg.InputText(bandpass_lower_cut_off, key = 'bandpass_lower_cut_off', size = (7,1)), sg.Text('upper Cut-off:'), sg.InputText(bandpass_upper_cut_off, key = 'bandpass_upper_cut_off', size = (7,1)), sg.Text('Filter order:'), sg.InputText(bandpass_order, key = 'bandpass_order', size = (7,1))],
            [sg.Push(), sg.Button('Confirm',button_color= ('white','black'), size = (18,1))]
            ]

            print ('*** Denoise window open. The main panel will be inactive until you close the window ***')
            denoise_window = sg.Window('Denoise parameters', denoise_layout)

            while True:
                denoise_event, denoise_values = denoise_window.read()

                if denoise_event == sg.WIN_CLOSED:
                    break

                moving_average = denoise_values['moving_average']
                box_moving_avg = denoise_values['box_moving_avg']
                gauss_moving_avg = denoise_values['gauss_moving_avg']
                low_pass_filter = denoise_values['low_pass_filter']
                bandpass_filter = denoise_values['bandpass_filter']

                try:
                    if moving_average and box_moving_avg:
                        box_moving_avg_size = int(denoise_values['box_moving_avg_size'])
                    if moving_average and gauss_moving_avg:
                        gauss_moving_avg_kernel = float(denoise_values['gauss_moving_avg_kernel'])
                        if gauss_moving_avg_kernel <= 0 or gauss_moving_avg_kernel > 1000:
                            sg.Popup('Gauss kernel must be greater than zero and smaller than 1000')
                            gauss_moving_avg_kernel = 5
                            continue

                    if low_pass_filter:
                        lowpass_cut_off = float(denoise_values['lowpass_cut_off'])
                        if lowpass_cut_off <=0 or lowpass_cut_off >= 1:
                            sg.Popup ('Low-pass cut-off must be greater than 0 and smaller than 1')
                            lowpass_cut_off = 0.1
                            continue
                        lowpass_order = int(denoise_values['lowpass_order'])
                        if lowpass_order <=0 or lowpass_order >30:
                            sg.Popup ('Filter order must be between 1 and 30')
                            lowpass_order = 4
                            continue
                    if bandpass_filter:
                        bandpass_lower_cut_off = float(denoise_values['bandpass_lower_cut_off'])
                        bandpass_upper_cut_off = float(denoise_values['bandpass_upper_cut_off'])
                        if bandpass_lower_cut_off >= bandpass_upper_cut_off:
                            sg.Popup ('Lower cut-off must be greater than lower cut-off')
                            bandpass_lower_cut_off = 0.1
                            bandpass_upper_cut_off = 0.5
                            continue
                        if bandpass_lower_cut_off <=0 or bandpass_lower_cut_off >= 1 or bandpass_upper_cut_off <=0 or bandpass_upper_cut_off >= 1:
                            sg.Popup('Bandpass cut-offs must be greater than 0 and smaller than 1')
                            bandpass_lower_cut_off = 0.1
                            bandpass_upper_cut_off = 0.5
                            continue
                        bandpass_order = int(denoise_values['bandpass_order'])
                        if bandpass_order <=0 or bandpass_order >30:
                            sg.Popup ('Bandpass filter order must be between 1 and 30')
                            bandpass_order = 4
                            continue
                except Exception:
                    sg.Popup ('Parameters not valid!')
                    box_moving_avg_size = 11
                    gauss_moving_avg_kernel = 5
                    lowpass_cut_off = 0.1
                    lowpass_order = 4
                    continue


                if denoise_event == 'Confirm':
                    print ('Denoise parameters confirmed. This main panel is now active again')
                    print ('')
                    break

            denoise_window.close()


        #5) DOPPLER CORRECTION PARAMETERS
        if spec_event  == ('Dopcor parameters'):

            sg.theme('LightBlue1')
            dopcor_layout = [

                [sg.Radio ('I have a velocity value', "RADIODOPVALUE", default = dop_cor_have_vel, key = 'dop_cor_have_vel', font = ('Helvetica', 11, 'bold'), tooltip='Usually a velocity value is used for stars and local (z < 0.01) galaxies'), sg.Radio ('I have a redshift (z) value', "RADIODOPVALUE", default = dop_cor_have_z, key = 'dop_cor_have_z', font = ('Helvetica', 11, 'bold'), tooltip='Usually a redshift (z) value is used for galaxies and has a cosmological meaning')],
                [sg.HorizontalSeparator()],
                [sg.Radio ('I have a list file: ', "RADIODOP", default = dop_cor_have_file, key = 'file_for_dopcor', tooltip='You can imput an ASCII file with the names of your spectra and a velocity or z value to correct. Check the readme_span file for details'), sg.InputText(dop_cor_file, size=(14, 1), key = 'dopcor_file'), sg.FileBrowse(tooltip='Load an ASCII file containing: Name of the spectrum, radial velocity to correct (km/s), in the same order of the original spectra list')],
                [sg.Radio('Single value: ', "RADIODOP", default = dop_cor_single_shot, key ='dopcor_single_value', tooltip='Using the same velocity or z value for all the spectra loaded'), sg.Text('Rec vel (km/s) or z: '), sg.InputText(dop_cor_single_shot_vel, size = (8,1), key = 'dopcor_value')],
                [sg.Push(), sg.Button('Confirm',button_color= ('white','black'), size = (18,1))]
                ]

            print ('*** Dopcor parameters window open. The main panel will be inactive until you close the window ***')
            dopcor_window = sg.Window('Dopcor parameters', dopcor_layout)

            while True:
                dopcor_event, dopcor_values = dopcor_window.read()

                if dopcor_event == sg.WIN_CLOSED:
                    break

                dop_cor_have_file = dopcor_values['file_for_dopcor']
                dop_cor_file = dopcor_values['dopcor_file']
                dop_cor_single_shot = dopcor_values['dopcor_single_value']
                dop_cor_have_vel = dopcor_values['dop_cor_have_vel']
                dop_cor_have_z = dopcor_values['dop_cor_have_z']


                try:
                    dop_cor_single_shot_vel = float(dopcor_values['dopcor_value'])
                    if dop_cor_have_z and dop_cor_single_shot_vel > 100:
                        sg.popup('Warning: the redshift (z) value inserted seems way too large. It is a velocity instead?')
                        continue
                except ValueError:
                    sg.popup('Dopcor value is not a number!')
                    continue

                if dopcor_event == 'Confirm':
                    print ('Dopcor parameters confirmed. This main panel is now active again')
                    print ('')
                    break

            dopcor_window.close()


        #6) HELIOCENTRIC CORRECTION PARAMETERS
        if spec_event  == ('Heliocor parameters'):

            sg.theme('LightBlue1')
            heliocor_layout = [

                [sg.Radio ('I have a file with location, date, RA and Dec. for all the spectra:', "RADIOHEL", default = helio_have_file, key = 'file_for_helio'), sg.InputText(helio_file, size = (37,1), key = 'helio_file'), sg.FileBrowse(tooltip='Load an ASCII file containing: Location, Date (YYYY-MM-DD), RA, Dec.')],
                [sg.Radio('Single correction', "RADIOHEL", default = helio_single_shot, key = 'helio_single_value'), sg.Text('Location:'), sg.InputText(helio_single_shot_location, size = (11,1), key = 'helio_location'), sg.Text('Date:'), sg.InputText(helio_single_shot_date, size = (10,1), key = 'helio_date'), sg.Text('RA:'), sg.InputText(ra_obj, size = (10,1), key = 'helio_ra'), sg.Text('Dec.:'), sg.InputText(dec_obj, size = (10,1), key = 'helio_dec'), sg.Button('loc.list',button_color=('black','light blue'),tooltip='Click to see the pre-loaded location list for heliocentric correction')],
                [sg.Push(), sg.Button('Confirm',button_color= ('white','black'), size = (18,1))]
                ]

            print ('*** Heliocor parameters window open. The main panel will be inactive until you close the window ***')
            heliocor_window = sg.Window('Heliocor parameters', heliocor_layout)

            while True:
                heliocor_event, heliocor_values = heliocor_window.read()

                if heliocor_event == sg.WIN_CLOSED:
                    break

                helio_have_file = heliocor_values['file_for_helio']
                helio_file = heliocor_values['helio_file']
                helio_single_shot = heliocor_values['helio_single_value']
                helio_single_shot_location = heliocor_values['helio_location']
                helio_single_shot_date = heliocor_values['helio_date']
                if helio_corr and helio_single_shot:
                    try:
                        ra_obj = float(heliocor_values['helio_ra'])
                        dec_obj = float(heliocor_values['helio_dec'])
                    except Exception:
                        sg.popup('Coordinates not valid!')
                        continue
                    try:
                        datetime.datetime.strptime(helio_single_shot_date, '%Y-%m-%d')
                    except Exception:
                        sg.popup ('Date format not valid. It must be: YYYY-MM-DD')
                        continue
                    try:
                        location = EarthLocation.of_site(helio_single_shot_location)
                    except Exception:
                        sg.popup ('Location not in the list')
                        continue
        #activating the button location list
                if(heliocor_event == 'loc.list'):
                    try:
                        location_list = EarthLocation.get_site_names()
                        sg.popup_scrolled(location_list, size=(120, 30))
                    except Exception:
                        sg.popup('Location list not available. I need an internet connection')

                if heliocor_event == 'Confirm':
                    print ('Heliocor parameters confirmed. This main panel is now active again')
                    print ('')
                    break

            heliocor_window.close()


    #************ Initializing and checking the variables of the Spectra processing frame ***********

        #1) REBINNING PARAMETERS
        #a) linear rebinning
        if rebinning and rebinning_linear:
            try:
                rebin_step_pix = float(spec_values['rebin_step_pix'])
                if rebin_step_pix <=0:
                    sg.popup('Invalid step. Must be > 0!')
                    continue
            except ValueError:
                sg.popup('Step is not a number!')
                continue


        #b) log rebinning
        if rebinning and rebinning_log:
            try:
                rebin_step_sigma = float(spec_values['rebin_step_sigma'])
                if rebin_step_sigma <1:
                    sg.popup('Invalid step. Must be >= 1!')
                    continue
            except ValueError:
                sg.popup('Step is not a number!')
                continue


        #2) DEGRADE RESOLUTION PARAMETERS
        if spec_event  == ('Degrade parameters'):

            sg.theme('LightBlue1')
            degrade_res_layout = [
                [sg.Radio('From R:', "RADIORESR", default = is_initial_res_r, key = 'is_initial_res_r', font = ('Helvetica', 12, 'bold') ), sg.InputText(initial_res_r, size = (6,1), key = 'degrade_from_r'), sg.Radio('to R:', "RADIORESRTOR", default = res_degrade_to_r, key = 'res_degrade_to_r' ), sg.InputText(final_res_r, size = (6,1), key = 'degrade_to_r'), sg.Radio('to FWHM (A):', "RADIORESRTOR", default = res_degrade_to_fwhm, key = 'res_degrade_to_fwhm'), sg.InputText(final_res_r_to_fwhm, size = (6,1), key = 'final_res_r_to_fwhm')],
                [sg.HorizontalSeparator()],
                [sg.Radio('From FWHM (A):', "RADIORESR", default = is_initial_res_fwhm, key = 'is_initial_res_fwhm', font = ('Helvetica', 12, 'bold')), sg.InputText(initial_res_fwhm, size = (4,1), key = 'degrade_from_l'), sg.Text('to FWHM (A):'), sg.InputText(final_res_fwhm, size = (4,1), key = 'degrade_to_l')],

                [sg.Radio('Degrade MUSE data', "RADIORESR", default = res_degrade_muse, key = 'res_degrade_muse', font = ('Helvetica', 12, 'bold')), sg.Text('to uniform FWHM (A):'), sg.InputText(res_degrade_muse_value, size = (4,1), key = 'res_degrade_muse_value')],
                [sg.Push(), sg.Button('Confirm',button_color= ('white','black'), size = (18,1))]
                ]

            print ('*** Degrade resolution parameters window open. The main panel will be inactive until you close the window ***')
            degrade_res_window = sg.Window('Degrade resolution parameters', degrade_res_layout)

            while True:
                degrade_res_event, degrade_res_values = degrade_res_window.read()

                if degrade_res_event == sg.WIN_CLOSED:
                    break

                try:
                    is_initial_res_r = degrade_res_values['is_initial_res_r']
                    initial_res_r = int(degrade_res_values['degrade_from_r'])
                    res_degrade_to_r = degrade_res_values['res_degrade_to_r']
                    final_res_r = int(degrade_res_values['degrade_to_r'])
                    res_degrade_to_fwhm = degrade_res_values['res_degrade_to_fwhm']
                    final_res_r_to_fwhm = float(degrade_res_values['final_res_r_to_fwhm'])

                    is_initial_res_fwhm = degrade_res_values['is_initial_res_fwhm']
                    initial_res_fwhm = float(degrade_res_values['degrade_from_l'])
                    final_res_fwhm = float(degrade_res_values['degrade_to_l'])
                    res_degrade_muse = degrade_res_values['res_degrade_muse']
                    res_degrade_muse_value = float(degrade_res_values['res_degrade_muse_value'])

                except Exception:
                    sg.Popup('Degrade resolution parameters not valid')
                    continue

                if initial_res_r <=0 or final_res_r <=0 or final_res_r_to_fwhm <=0 or initial_res_fwhm <=0 or final_res_fwhm <=0 or res_degrade_muse_value <= 0:
                    sg.Popup ('Resolution values cannot be negative or zero!')
                    continue

                if final_res_fwhm < initial_res_fwhm or initial_res_r<final_res_r:
                    sg.popup('You want to improve the resolution? That''s impossible!')
                    continue

                if degrade_res_event == 'Confirm':
                    print ('Degrade resolution parameters confirmed. This main panel is now active again')
                    print ('')
                    break

            degrade_res_window.close()


        #3) NORMALISATION PARAMETERS
        if normalize_wave:
            try:
                norm_wave = float(spec_values['norm_wave'])
            except ValueError:
                sg.popup('Normalisation wave not valid!')
                continue


        #4) SIGMA BROADENING PARAMETERS
        if sigma_broad:
            try:
                sigma_to_add = float(spec_values['sigma_to_add'])
                if sigma_to_add < 0:
                    sg.popup('Invalid sigma broadening. Must be >= 0!')
                    continue
            except ValueError:
                sg.popup('Sigma broadening not valid!')
                continue

        #5) ADD NOISE PARAMETERS
        if add_noise:
            try:
                noise_to_add = float(spec_values['noise_to_add'])
                if noise_to_add <=0:
                    sg.popup('Invalid SNR. Must be > 0!')
                    continue
            except ValueError:
                sg.popup('Noise value not valid!')
                continue


        #6) CONTINUUM MODELLING PARAMETERS
        if spec_event == 'Continuum parameters':
            sg.theme('LightBlue1')
            continuum_layout = [
            [sg.Radio('Continuum model: automatic filtering of the spectrum (works good for smooth spectrum and no emission)', "CONTMODE", default = cont_model_filtering, key = 'cont_model_filtering', font = ('Helvetica', 11, 'bold'))],
            [sg.HorizontalSeparator()],
            [sg.Radio('Continuum model: fine-tuning polynomial fitting', "CONTMODE", default = cont_model_poly, key = 'cont_model_poly',font = ('Helvetica', 11, 'bold')), sg.Checkbox('Regions to mask:', default = cont_want_to_mask, key = 'cont_want_to_mask'), sg.InputText(cont_mask_ranges_str, size = (14,1), key = 'cont_mask_ranges'), sg.Text('Polynomial degree:'), sg.InputText(cont_poly_degree, size = (6,1), key = 'cont_poly_degree')],
            [sg.HorizontalSeparator()],
            [sg.Text('Operation on the spectrum:'), sg.InputCombo(markers_cont_operations, key='markers_cont_operations', default_value=cont_math_operation, readonly=True)],
            [sg.Push(), sg.Button('Confirm',button_color= ('white','black'), size = (18,1))]
            ]

            print ('*** Continuum subtraction window open. The main panel will be inactive until you close the window ***')
            continuum_window = sg.Window('Continuum parameters', continuum_layout)

            while True:
                continuum_event, continuum_values = continuum_window.read()

                if continuum_event == sg.WIN_CLOSED:
                    break

                cont_model_filtering = continuum_values['cont_model_filtering']
                cont_model_poly = continuum_values['cont_model_poly']
                cont_math_operation = continuum_values['markers_cont_operations']

                if cont_model_poly:

                    cont_want_to_mask = continuum_values['cont_want_to_mask']
                    if cont_want_to_mask:
                        try:
                            cont_mask_ranges_str = continuum_values['cont_mask_ranges']
                            cont_mask_ranges = eval(cont_mask_ranges_str)
                        except Exception:
                            sg.Popup('Masking values not valid')
                            continue
                    try:
                        cont_poly_degree = int(continuum_values['cont_poly_degree'])
                    except Exception:
                        sg.Popup ('Polynomial degree not valid')
                        continue
                    if cont_poly_degree <0 or cont_poly_degree > 11:
                        sg.Popup('Polynomial degree must be between 0 and 11')
                        cont_poly_degree = 5
                        continue


                if continuum_event == 'Confirm':
                    print ('Continuum parameters confirmed. This main panel is now active again')
                    print ('')
                    break

            continuum_window.close()


        #9) add pedestal and check on the input values
        if add_pedestal:
            try:
                pedestal_to_add = float(spec_values['pedestal_to_add'])
            except ValueError:
                sg.popup('Pedestal value not valid!')
                continue

        #10) multiply by a constant and check on the input values
        if multiply:
            try:
                multiply_factor = float(spec_values['multiply_factor'])
                if multiply_factor <=0:
                    sg.popup('Invalid multiply constant. Must be > 0!')
                    continue
            except ValueError:
                sg.popup('Multiply value not valid!')
                continue

        #Help file
        if spec_event == 'I need help':
            f = open(os.path.join(BASE_DIR, "help_files", "need_help_spec_proc.txt"), 'r')
            file_contents = f.read()
            sg.popup_scrolled(file_contents, size=(120, 30))

        #Confirm the parameters
        if spec_event == 'Confirm':
            print ('Spectra manipulation parameters confirmed. This main panel is now active again')
            print ('')
            break


    #closing the window
    spec_window.close()

    #updating the parameters with the values of the local variables.
    params = replace(params,
        # Spectra pre-processing
        cropping_spectrum=cropping_spectrum,
        cropping_low_wave=cropping_low_wave,
        cropping_high_wave=cropping_high_wave,
        sigma_clipping=sigma_clipping,
        clip_factor=clip_factor,
        sigma_clip_resolution=sigma_clip_resolution,
        sigma_clip_single_vel=sigma_clip_single_vel,
        sigma_clip_single_value=sigma_clip_single_value,
        sigma_clip_have_file=sigma_clip_have_file,
        sigma_clip_sigma_file=sigma_clip_sigma_file,
        wavelet_cleaning=wavelet_cleaning,
        sigma_wavelets=sigma_wavelets,
        wavelets_layers=wavelets_layers,

        # Filtering & denoising
        filter_denoise=filter_denoise,
        moving_average=moving_average,
        box_moving_avg=box_moving_avg,
        box_moving_avg_size=box_moving_avg_size,
        gauss_moving_avg=gauss_moving_avg,
        gauss_moving_avg_kernel=gauss_moving_avg_kernel,
        low_pass_filter=low_pass_filter,
        lowpass_cut_off=lowpass_cut_off,
        lowpass_order=lowpass_order,
        bandpass_filter=bandpass_filter,
        bandpass_lower_cut_off=bandpass_lower_cut_off,
        bandpass_upper_cut_off=bandpass_upper_cut_off,
        bandpass_order=bandpass_order,

        # Doppler & Heliocentric correction
        dop_cor=dop_cor,
        dop_cor_single_shot_vel=dop_cor_single_shot_vel,
        dop_cor_have_file=dop_cor_have_file,
        dop_cor_file=dop_cor_file,
        dop_cor_single_shot=dop_cor_single_shot,
        dop_cor_have_vel=dop_cor_have_vel,
        dop_cor_have_z=dop_cor_have_z,
        helio_corr=helio_corr,
        helio_have_file=helio_have_file,
        helio_file=helio_file,
        helio_single_shot=helio_single_shot,
        helio_single_shot_location=helio_single_shot_location,
        helio_single_shot_date=helio_single_shot_date,
        ra_obj=ra_obj,
        dec_obj=dec_obj,

        # Resolution degradation
        is_initial_res_r=is_initial_res_r,
        initial_res_r=initial_res_r,
        res_degrade_to_r=res_degrade_to_r,
        final_res_r=final_res_r,
        res_degrade_to_fwhm=res_degrade_to_fwhm,
        final_res_r_to_fwhm=final_res_r_to_fwhm,
        is_initial_res_fwhm=is_initial_res_fwhm,
        initial_res_fwhm=initial_res_fwhm,
        final_res_fwhm=final_res_fwhm,
        res_degrade_muse = res_degrade_muse,
        res_degrade_muse_value = res_degrade_muse_value,

        # Continuum operations
        markers_cont_operations=markers_cont_operations,
        cont_math_operation=cont_math_operation,
        cont_model_filtering=cont_model_filtering,
        cont_model_poly=cont_model_poly,
        cont_want_to_mask=cont_want_to_mask,
        cont_mask_ranges_str=cont_mask_ranges_str,
        cont_mask_ranges=cont_mask_ranges,
        cont_poly_degree=cont_poly_degree,

        # Spectra processing
        rebinning=rebinning,
        rebinning_linear=rebinning_linear,
        rebin_step_pix=rebin_step_pix,
        rebinning_log=rebinning_log,
        rebin_step_sigma=rebin_step_sigma,
        degrade=degrade,
        normalize_wave=normalize_wave,
        norm_wave=norm_wave,
        sigma_broad=sigma_broad,
        sigma_to_add=sigma_to_add,
        add_noise=add_noise,
        noise_to_add=noise_to_add,
        continuum_sub=continuum_sub,

        # Spectral math operations
        subtract_normalized_avg=subtract_normalized_avg,
        subtract_normalized_spec=subtract_normalized_spec,
        spectra_to_subtract=spectra_to_subtract,
        add_pedestal=add_pedestal,
        pedestal_to_add=pedestal_to_add,
        multiply=multiply,
        multiply_factor=multiply_factor,
        derivatives=derivatives,
        average_all=average_all,
        norm_and_average=norm_and_average,
        do_nothing=do_nothing,
        sum_all=sum_all,
        normalize_and_sum_all=normalize_and_sum_all,
        use_for_spec_an=use_for_spec_an,

        # Order & active operations
        reorder_op=reorder_op,
        current_order=current_order,
        reordered_operations=reordered_operations,
        active_operations=active_operations,

        # Spectra saving options
        save_intermediate_spectra = save_intermediate_spectra,
        save_final_spectra = save_final_spectra,
        not_save_spectra = not_save_spectra

    )

    return params
