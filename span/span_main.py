#SPectral ANalysis software (SPAN).
#Written by Daniele Gasparri#

"""
    Copyright (C) 2020-2025, Daniele Gasparri

    E-mail: daniele.gasparri@gmail.com

    SPAN is a GUI software that allows to modify and analyze 1D astronomical spectra.

    1. This software is licensed for non-commercial, academic and personal use only.
    2. The source code may be used and modified for research and educational purposes, 
    but any modifications must remain for private use unless explicitly authorized 
    in writing by the original author.
    3. Redistribution of the software in its original, unmodified form is permitted 
    for non-commercial purposes, provided that this license notice is always included.
    4. Redistribution or public release of modified versions of the source code 
    is prohibited without prior written permission from the author.
    5. Any user of this software must properly attribute the original author 
    in any academic work, research, or derivative project.
    6. Commercial use of this software is strictly prohibited without prior 
    written permission from the author.

    DISCLAIMER:
    THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES, OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT, OR OTHERWISE, ARISING FROM, OUT OF, OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""

import numpy as np
import sys
import os
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.backends.backend_tkagg
import time
import json
from dataclasses import replace
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__))) #adding the folder to python path

#SPAN import
from span_imports import *

#Define the base dir of SPAN in your device
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

        
################ LET'S START #######################################################
def main():
    # Loading the GUI theme and the layout for the operating system in use
    layout, scale_win, fontsize, default_size = misc.get_layout()
    
    #Creating the main GUI
    window1 = sg.Window('SPAN - SPectral ANalysis - 7.0 --- Daniele Gasparri ---', layout,finalize=True, resizable=True, scaling = scale_win, modal =False)

    #Allowing elements in the listbox to be deleted
    listbox_widget = window1['-LIST-'].Widget
    listbox_widget.bind("<Delete>", lambda event: window1.write_event_value('-LIST-DELETE-', None))
    listbox_widget.bind("<Double-1>", lambda event: window1.write_event_value('-LIST-DOUBLECLICK-', None))
    last_state = None #placeholder to store the listbox state for undo operations

    #Calling zooming function and use ti
    zm = zoom.ZoomManager.get()
    zm.attach_window(window1)

    # Creating an initial mock preview, to be updated later
    if layout == layouts.layout_windows:
        layout_name = 'windows'
    elif layout == layouts.layout_linux:
        layout_name = 'linux'
    elif layout == layouts.layout_macos:
        layout_name = 'macos'
    elif layout == layouts.layout_android:
        layout_name = 'android'
    else:
        layout_name = 'default'
    fig, ax, _plot_line, hud_text, _preview_canvas = preview_tools.create_preview(layout_name, window1, preview_key='-CANVAS-')
    (_plot_line2,) = ax.plot([], [], lw=0.8, color="red") #creeate second plot for comparison
    
    # wrapping to the provider SNR function which will show the S/N on the preview
    def _snr_provider(lam_x: float):
        return preview_tools.snr_provider(lam_x, ax, preview_interactor)

    # --- Make the preview interactive
    preview_interactor = preview_tools.PreviewInteractor(ax, status_setter=None, get_snr=_snr_provider, hud_text=hud_text, snr_mode="points", snr_halfwin_A=50.0, snr_halfwin_pts=50)
    
    # Redshift estimator overlay
    redshift_shifter = preview_tools.SpectrumShifterInteractor(ax, _plot_line, hud_text=hud_text)
    
    #Loading parameters from the dataclass
    params = SpectraParams()

    # Force centrering the window on the screen
    window1.move(sg.Window.get_screen_size()[0]//2 - window1.size[0]//2, sg.Window.get_screen_size()[1]//2 - window1.size[1]//2)

    # Initializing events, values and keys of the GUI
    keys, events, values = [], [], {}

    # calling the function to check the existence of the SpectralTemplates folder. If not exist, ask the user to download
    misc.check_and_download_spectral_templates()

    #Checking the existence of the default_settings.json file. If not, create it when the GUI opens
    DEFAULT_PARAMS_FILE = os.path.join(BASE_DIR, "system_files", "default_settings.json")
    if not os.path.exists(DEFAULT_PARAMS_FILE):
        window1.write_event_value('-INIT-', None)

    # Prints in the output
    print ('***********************************************')
    print ('********* Welcome to SPAN version 7.0 *********')
    print ('********* Written by Daniele Gasparri *********')
    print ('***********************************************\n')
    print ('SPAN is a software for performing operations and analyses on 1D reduced astronomical spectra.\n')
    print ('If you just click the Load! button, the example files are loaded and you can make some practise.\n')
    print ('NOTE: all the SPAN wavelength units are expressed in Angstrom')
    print ('***********************************************')
    print (f'SPAN will save the results in: {params.result_data}\n')

    # starting the dynamic main GUI window
    while True:

        #taking the actual time to store in the output files
        timestamp = time.strftime("%Y%m%d_%H%M%S")

        #starting all the GUI windows with their values from the SpectraParams dataclass
        window, event, values = sg.read_all_windows()

        # handle Zooming events
        if event == 'Zoom In':
            zm.zoom_in()
            continue
        elif event == 'Zoom Out':
            zm.zoom_out()
            continue
        elif event == 'Reset Zoom':
            zm.set_scale(1.0)
            continue
        
        #---- Listbox related events, only if I load a list of spectra ----
        if event == "-LIST-" and not values['one_spec']:
            params = listbox_events.handle_list_select(event, values, window, params, _plot_line, ax, fig, preview_interactor, redshift_shifter, stm)
            _plot_line2.set_data([], []) #reset second plot, if any
            
        elif event == "-LIST-DELETE-" and not values['one_spec']:
            params, last_state = listbox_events.handle_list_delete(event, values, window, params, _plot_line, ax, fig, preview_interactor, redshift_shifter)

        elif event in ('↑ Move Up', '↓ Move Down', 'Remove') and not values['one_spec']:
            params, last_state = listbox_events.handle_list_menu(event, values, window, params, _plot_line, ax, fig, preview_interactor, redshift_shifter, last_state)

        elif event == "Undo" and not values['one_spec']:
            params = listbox_events.handle_undo(event, values, window, params, last_state)

        elif event == "Save current spectra list...":
            listbox_events.handle_save_list(event, values, window, params)

        elif event == "Compare spectra" and not values['one_spec']:
            params = listbox_events.handle_compare_spectra(event, values, window, params, _plot_line, _plot_line2, ax, fig, preview_interactor, redshift_shifter, stm)

        elif event == "-LIST-DOUBLECLICK-" and not values['one_spec']:
            listbox_events.handle_list_doubleclick(values, window, params, stm)

        #Automatic event to save the default_settings.json file if does not exist.
        if not os.path.exists(DEFAULT_PARAMS_FILE):
            if event == '-INIT-':
                try:
                    settings.save_settings(DEFAULT_PARAMS_FILE, keys, events, values, params)
                    print(f"\nAutomatically stored the default parameters in {DEFAULT_PARAMS_FILE}")
                except json.JSONDecodeError:
                    print('Cannot create the default parameter file')

        #Closing event
        if event in (sg.WIN_CLOSED, 'Exit'):
            break

        #Showing the result folder path
        if event == 'Show result folder':
            sg.popup(f"The SPAN_result folder is in: {params.result_data}")

        # If the user wants to change the SPAN_results directory
        elif event == 'Change result folder...':
            config_file = os.path.join(BASE_DIR, "system_files", "config.json")
            config_folder = misc.load_config(config_file)
            misc.change_result_path(config_folder, config_file)
            params = replace(params, result_path = config_folder["result_path"])
            print ('\nSPAN will now save the results in ', params.result_data)

        # to clean the Output, only if integrated in the GUI
        if event == 'Clean output' and '-OUTPUT-' in window.key_dict:
            window['-OUTPUT-'].update('')
        elif event == 'Clean output' and not '-OUTPUT-' in window.key_dict:
            print("You have external output, cannot clean it!")

        # assign the spectra list value to dataclass parameters
        params = replace(params, spectra_list=values['spec_list'], spectra_list_name=os.path.splitext(os.path.basename(values['spec_list']))[0])
       
       #create a spectra list file
        if event == 'listfile':
            params = settings.generate_spectra_list(window, params)

        #assigning lambda units of the spectra
        if (values['wave_units_nm']):
            params.lambda_units = 'nm'
        elif (values['wave_units_a']):
            params.lambda_units = 'a'
        elif (values['wave_units_mu']):
            params.lambda_units = 'mu'

        # About and Version sections:
        if event == 'About SPAN':
            sg.popup ('SPAN is a Python 3.X software. It can modify the spectra and perform analyses, using both built-in and external (e.g. ppxf) algorithms\n\nSPAN uses FreeSimpleGUI (Copyright (C) 2007 Free Software Foundation, Inc.), which is distributed under the GNU LGPL license. ')
        elif event == 'Version':
            sg.popup ('This is version 7.0 with improved, dynamical, and responsive layout')

        # In the case I want to deselect all the active tasks in the main panel in one click
        elif event == 'Clear all tasks':
            params = settings.clear_all_tasks(window, params)


        # --- Preview frame buttons: Header/Sampling buttons ---
        elif event == '-SHOWHDR-':
            if params.prev_spec == '' or params.spectra_number == 0:
                sg.popup('No spectrum selected. Please, select one in the list.')
            else:
                try:
                    utility_tasks.show_fits_header(params.prev_spec)
                except Exception as e:
                    sg.popup(f'Failed to show header: {e}')
        elif event == '-SHOWSTEP-':
            if params.prev_spec == '' or params.spectra_number == 0:
                sg.popup('No spectrum selected. Please, select one in the list.')
            else:
                try:
                    wl, _, *_ = stm.read_spec(params.prev_spec, params.lambda_units)
                    utility_tasks.show_sampling(wl)
                except Exception as e:
                    sg.popup(f'Failed to show sampling: {e}')
    
    
    #******************* INITIALIZING AND CHECKING THE VARIABLES OF THE UTILITIES FRAME *****************
        # Utilities window
        elif event == 'Utilities':
            if params.spectra_number == 0:
                sg.popup('No spectrum loaded. Please, load one first.')
                continue
            one_spec_flag = values.get('one_spec', False)
            try:
                params = sub_programs.utilities_window(params, one_spec_flag)
            except Exception as e:
                sg.popup(f'Utilities window error: {e}')
            continue

    #********** INITIALIZING AND CHECKING THE SPECTRA MANIPULATION PANEL **********
        elif event == 'Spectra manipulation':
            params = spec_manipulation.spectra_manipulation(params)

        # detecting if at least one task has been activated in the Spectra Manipulation panel and if True change the color of the button
        any_active = any([params.cropping_spectrum, params.sigma_clipping, params.wavelet_cleaning,
                        params.filter_denoise, params.dop_cor, params.helio_corr, params.rebinning,
                        params.degrade, params.normalize_wave, params.sigma_broad, params.add_noise,
                        params.continuum_sub, params.average_all, params.norm_and_average, params.sum_all,
                        params.normalize_and_sum_all, params.subtract_normalized_avg, params.subtract_normalized_spec,
                        params.add_pedestal, params.multiply, params.derivatives,])
        if any_active:
            window['Spectra manipulation'].update(button_color=('white', 'red'))
        else:
            window['Spectra manipulation'].update(button_color= ('black','light blue'))


    #********************************** STAND ALONE SUB-PROGRAMS *************************
        # LONG-SLIT SPECTRA EXTRACTION
        if event == 'Long-slit extraction':
            params = sub_programs.long_slit_extraction(BASE_DIR, layout, params)

        # CUBE EXTRACTION
        elif event == 'DataCube extraction':
            params = sub_programs.datacube_extraction(params)

        # TEXT EDITOR
        elif event == 'Text editor':
            sub_programs.text_editor_window(layout)

        # FITS HEADER EDITOR
        elif event == 'FITS header editor':
            sub_programs.fits_header_window()

        # DATA PLOTTING
        elif event == 'Plot data':
            sub_programs.plot_data_window(BASE_DIR, layout)

        # 2D MAPS PLOTTING
        elif event == 'Plot maps':
            params = sub_programs.plot_maps_window(BASE_DIR, layout, params)


    #********************************** LOADING AND CHECKING THE SPECTRA *************************
        if not values['one_spec']: #I have a list, reading the list.
            #if I press the Load! button, then I print some info and update the list of spectra.
            if event == "Load!":
                params = replace(params, **dict(zip(["spectra_number", "spec_names", "spec_names_nopath", "fatal_condition"],
                    check_spec.load_and_validate_spectra(params.spectra_list, params.lambda_units, window))))
                if params.fatal_condition:
                    continue
            if len(params.spec_names) != params.spectra_number and params.spectra_number > 0:
                sg.popup ('The format of the spectra list file is not correct. Try to adjust the spectra file list')
                continue
            if params.fatal_condition:
                sg.popup ('You did not load any valid spectra. I can do nothing but show you this message until you will load a valid spectra list')
                continue
            #showing a message if no spectra are loaded and button are pressed
            if (params.spec_names[0] == 0 and (event == 'Preview spec.' or event == 'Process selected' or event == 'Show info' or event == 'Preview result' or event == 'Process all' or event == 'Plot' or event == 'One' or event == 'All' or event == 'Compare' or event == 'convert_one' or event == 'convert_all' or event == 'Show snr' or event == 'See plot' or event == 'Save one' or event == 'Save all')):
                sg.popup('Please, load some spectra!')
                continue

            #Define the names to show in the GUI, without the path. Only for visualisation purposes!
            if values['-LIST-']:
                selected_display = values['-LIST-'][0]
                if hasattr(window, "metadata") and window.metadata:
                    sel_full = window.metadata.get(selected_display, "")
                else:
                    sel_full = next((s for s in params.spec_names if isinstance(s, str) and os.path.basename(s) == str(selected_display)),"")
                params = replace(params, prev_spec=sel_full, prev_spec_nopath=os.path.splitext(os.path.basename(sel_full))[0] if sel_full else "")

        # If I load a single spectrum, SPAN needs to check it before loading
        elif values['one_spec']:
            if event == 'Load!':
                # Validate and load spectrum
                params, valid_spec = check_spec.validate_and_load_spectrum(params, window)

                try:
                    single_path = values.get('spec_list', '')
                    if single_path:
                        params = replace(params, prev_spec=single_path)
                    params = replace(params, prev_spec_nopath=os.path.splitext(os.path.basename(params.prev_spec))[0])
                    params = listbox_events.handle_list_select(event, values, window, params, _plot_line, ax, fig, preview_interactor, redshift_shifter, stm, sel_full_override=single_path)
                    
                except Exception as e:
                    print(f"[Preview one_spec] draw failed: {e}")

                if not valid_spec:
                    sg.popup('The format of the spectrum is not correct or you did not load it')
                    continue
            try:
                if not valid_spec:
                    sg.popup("Your spectrum is not valid. Can't do anything")
                    continue
            except Exception:
                sg.popup('You should load your spectrum if you want to use it')
                continue
    
        # Concatenating events to prevent the GUI to crash when no (valid) spectrum is selected or loaded and you want to do something anyway.
        if ( (event == 'Preview spec.' or event == 'Process selected' or event == 'Show info' or event == 'Preview result' or event == 'Plot' or event == 'See plot' or event == 'Save one' or event == 'One' or event == 'All' or event == 'Compare' or event == 'convert_one' or event == 'convert_all' or event == 'Show snr') and params.prev_spec == ''):
            sg.popup('No spectrum selected. Please, select one spectrum in the list. Doing nothing')
            continue


    #****************** SUB WINDOWS DEFINITION AND PARAMETERS OF THE SPECTRAL ANALYSIS FRAME *****************************
        #1) BLACKBODY PARAMETERS
        bb_fit = values['bb_fitting']
        if (event == 'Blackbody parameters'):
            params = param_windows.blackbody_parameters(params)

        #2) CROSS-CORRELATION PARAMETERS
        cross_corr = values['xcorr']
        if (event == 'Cross-corr parameters'):
            params = param_windows.crosscorr_parameters(params)

        # 3) VELOCITY DISPERSION PARAMETERS
        sigma_measurement = values['sigma_measurement']
        if (event == 'Sigma parameters'):
            params = param_windows.sigma_parameters(params)

        #4) EQUIVALENT WIDTH PARAMETERS
        ew_measurement = values['ew_measurement']
        if (event == 'Line-strength parameters'):
            params = param_windows.line_strength_parameters(params)

        #5) LINE(S) FITTING PARAMETERS
        line_fitting = values['line_fitting']
        if event == 'Line fitting parameters':
            params = param_windows.line_fitting_parameters(params)

        #6) KINEMATICS WITH PPXF
        perform_kinematics = values['ppxf_kin']
        if (event == 'Kinematics parameters'):
            params = param_windows.kinematics_parameters(params)

        #7) STELLAR POPULATIONS WITH PPXF
        stellar_pop = values['ppxf_pop']
        if (event == 'Population parameters'):
            params = param_windows.population_parameters(params)

        save_plot = values['save_plots']

     #********************************** PLOT EVENT *************************
        if event == 'Plot':
            wavelength, flux, step, name = stm.read_spec(params.prev_spec, params.lambda_units)
            plt.plot(wavelength, flux)
            plt.xlabel('Wavelength (Å)', fontsize = 9)
            plt.title(params.prev_spec_nopath)
            plt.ylabel('Flux')
            plt.show()
            plt.close()

     #********************************** LOADING AND CHECKING THE SPECTRAHELP FILES *************************
        sg.theme('Lightblue1')
        if event == 'Read me':
            f = open(os.path.join(BASE_DIR, "help_files", "readme_span.txt"), 'r', encoding="utf-8")
            file_contents = f.read()
            if layout == layouts.layout_android:
                sg.popup_scrolled(file_contents, size=(120, 30))
            else:
                sg.popup_scrolled(file_contents, size=(100, 40))

        elif event == 'I need help':
            f = open(os.path.join(BASE_DIR, "help_files", "need_help_spec_proc.txt"), 'r')
            file_contents = f.read()
            if layout == layouts.layout_android:
                sg.popup_scrolled(file_contents, size=(120, 30))
            else:
                sg.popup_scrolled(file_contents, size=(100, 40))

        elif event == 'Help me':
            f = open(os.path.join(BASE_DIR, "help_files", "help_me_spec_analysis.txt"), 'r')
            file_contents = f.read()
            if layout == layouts.layout_android:
                sg.popup_scrolled(file_contents, size=(120, 30))
            else:
                sg.popup_scrolled(file_contents, size=(100, 40))

        elif event == 'Quick start':
            f = open(os.path.join(BASE_DIR, "help_files", "quick_start.txt"), 'r')
            file_contents = f.read()
            if layout == layouts.layout_android:
                sg.popup_scrolled(file_contents, size=(120, 30))
            else:
                sg.popup_scrolled(file_contents, size=(100, 40))

        elif event == 'Tips and tricks':
            f = open(os.path.join(BASE_DIR, "help_files", "tips_and_tricks.txt"), 'r')
            file_contents = f.read()
            if layout == layouts.layout_android:
                sg.popup_scrolled(file_contents, size=(120, 30))
            else:
                sg.popup_scrolled(file_contents, size=(100, 40))
        
        elif event == 'SPAN Manual':
            stm.open_manual()

        # the following lines are just to ensure that if you do not load any spectra but opens the sub-programs, SPAN will not crash
        try:
            original_flux = flux
            original_wavelength = wavelength
        except Exception:
            pass

     
    #************** PREPARING THE ASCII FILES FOR SPECTRAL ANALYSIS RESULTS IN PROCESS ALL MODE *************
        if (event == 'Process all' and not values['one_spec']):

            #1) Blackbody
            if (bb_fit):
                bb_file = files_setup.create_blackbody_file(params.result_bb_dir, params.spectra_list_name, timestamp, params.spectra_number, params.spec_names_nopath)
                df_bb = pd.read_csv(bb_file, sep=' ', index_col=0)

            #2) Cross correlation
            if (cross_corr):
                rv_file = files_setup.create_cross_correlation_file(params.result_xcorr_dir, params.spectra_list_name, timestamp, params.spectra_number, params.spec_names_nopath)
                df_rv = pd.read_csv(rv_file, sep=' ', index_col=0)

            #3) Velocity dispersion measurement
            if (sigma_measurement):
                sigma_file = files_setup.create_velocity_dispersion_file(params.result_vel_disp_dir, params.spectra_list_name, timestamp, params.spectra_number, params.spec_names_nopath)
                df_sigma = pd.read_csv(sigma_file, sep=' ', index_col=0)

            #4) EW measurement
            if (ew_measurement and params.single_index):
                ew_file, ew_file_mag, snr_ew_file = files_setup.create_ew_measurement_files(
                params.result_ew_data_dir, params.spectra_list_name, timestamp, params.spectra_number, params.spec_names_nopath)
                df_ew = pd.read_csv(ew_file, sep=' ', index_col=0)
                df_ew_mag = pd.read_csv(ew_file_mag, sep=' ', index_col=0)
                df_snr_ew = pd.read_csv(snr_ew_file, sep=' ', index_col=0)

            if (ew_measurement and params.have_index_file):
                try:
                    ew_file, ew_file_mag, snr_ew_file, num_indices, ew_id, spectra_id, ew_id_mag, snr_ew_id = files_setup.create_ew_measurement_files_from_index(
                    params.result_ew_data_dir, params.spectra_list_name, timestamp, params.spectra_number, params.spec_names_nopath, params.index_file)
                    df_ew = pd.read_csv(ew_file, sep=' ', index_col=0)
                    df_ew_mag = pd.read_csv(ew_file_mag, sep=' ', index_col=0)
                    df_snr_ew = pd.read_csv(snr_ew_file, sep=' ', index_col=0)
                except Exception:
                    print('At least one index is not valid. Stopping')
                    continue

            if (ew_measurement and params.lick_ew):
                try:
                    ew_lick_file, ew_lick_file_mag, snr_lick_ew_file, num_lick_indices, ew_lick_id, spectra_lick_id, ew_lick_id_mag, snr_lick_ew_id = files_setup.create_lick_ew_measurement_files(
                        params.result_ew_data_dir, params.spectra_list_name, timestamp, params.spectra_number, params.spec_names_nopath, params.lick_index_file)
                    df_ew_lick = pd.read_csv(ew_lick_file, sep=' ', index_col=0)
                    df_ew_lick_mag = pd.read_csv(ew_lick_file_mag, sep=' ', index_col=0)
                    df_snr_lick_ew = pd.read_csv(snr_lick_ew_file, sep=' ', index_col=0)

                    ssp_lick_param_file = files_setup.create_ssp_lick_param_file(params.result_ew_data_dir, params.spectra_list_name, timestamp, params.spectra_number, params.spec_names_nopath)
                    df_lick_param = pd.read_csv(ssp_lick_param_file, sep=' ', index_col=0)
                except Exception:
                    print('The Lick index file in /system_files does not exist. Skipping...')
                    continue

            if (ew_measurement and params.lick_ew and params.stellar_parameters_lick):
                ssp_param_file = files_setup.create_lick_ssp_parameters_file(
                    params.result_ew_data_dir, params.spectra_list_name, timestamp, params.spectra_number, params.spec_names_nopath)
                df_ssp_param = pd.read_csv(ssp_param_file, sep=' ', index_col=0)

            #5) Line(s) fitting
            if (line_fitting and params.cat_band_fit):
                fit_file = files_setup.create_line_fitting_file(
                params.result_line_fitting_dir, params.spectra_list_name, timestamp, params.spectra_number, params.spec_names_nopath)
                df_fit = pd.read_csv(fit_file, sep=' ', index_col=0)

            #5-bis) Line fitting 2
            if (line_fitting and not params.cat_band_fit):
                fit_file = files_setup.create_line_fitting_file_simple(
                params.result_line_fitting_dir, params.spectra_list_name, timestamp, params.spectra_number, params.spec_names_nopath)
                df_fit = pd.read_csv(fit_file, sep=' ', index_col=0)

            #6) kinematics with ppxf
            if (perform_kinematics):
                kinematics_files = files_setup.create_kinematics_files(
                    params.result_ppxf_kin_data_dir, params.spectra_list_name, timestamp,
                    params.spectra_number, params.spec_names_nopath, params.with_errors_kin, params.gas_kin)
                kin_file = kinematics_files.get("stellar")  # Stellar kinematics file
                kin_file_mc = kinematics_files.get("stellar_mc")  # MC errors file (if exists)
                kin_file_gas = kinematics_files.get("gas")  # Gas kinematics file (if exists)
                df_kin = pd.read_csv(kin_file, sep=' ', index_col=0)
                df_kin_mc = None
                df_kin_gas = None
                if kin_file_mc: #kin_mc is not guaranteed to exist, so I must check
                    df_kin_mc = pd.read_csv(kin_file_mc, sep=' ', index_col=0)

            #7) Stellar populations with ppxf
            if (stellar_pop):
                pop_files = files_setup.create_stellar_population_files(
                    params.result_ppxf_pop_data_dir, params.spectra_list_name, timestamp, params.spectra_number,
                    params.spec_names_nopath, params.ppxf_pop_lg_age, params.stellar_parameters_lick_ppxf)
                pop_file = pop_files.get("stellar_population")  # Main stellar population file
                ssp_param_file_ppxf = pop_files.get("ssp_lick_parameters")
                df_pop = pd.read_csv(pop_file, sep=' ', index_col=0)
                df_ssp_param_ppxf = None 
                if ssp_param_file_ppxf: #Lick/IDS stellar pop file is not guaranteed to exist, so I must check
                    df_ssp_param_ppxf = pd.read_csv(ssp_param_file_ppxf, sep=' ', index_col=0)


    #***************************** MEGA EVENT: APPLY TASKS *******************************************
        if (event == 'Preview spec.' or event == 'Preview result' or event == 'Process selected' or event == 'Process all'):
            #resetting the task(s) counter
            params = replace(params, task_done=0, task_done2=0, task_analysis=0, task_spec = 0, task_spec2 =0)

            #If I want to process only the selected spectrum 'prev_spec', I change its name:
            if event == 'Preview spec.' or event == 'Process selected' or event == 'Preview result' or values['one_spec']:
                params = replace(params, spectra_number_to_process=1, spec_names_to_process=[params.prev_spec], spec_names_nopath_to_process=[params.prev_spec_nopath])

            # If I want to process all spectra:
            if event == 'Process all':
                if values['one_spec']:
                    sg.popup('With one spectrum loaded, you need to use ''Process selected''')
                    continue
                else: # renaming the spectra to not overwrite the original list loaded
                    params = replace(params, spec_names_to_process=params.spec_names, spec_names_nopath_to_process=params.spec_names_nopath)
                    #Removing the extension of the file for aesthetic purposes
                    params = replace(params, spec_names_nopath_to_process = [os.path.splitext(name)[0] for name in params.spec_names_nopath_to_process])
                    params = replace(params, spectra_number_to_process = params.spectra_number)

            ###################### Lest's start! --> Cycle for all the spectra #########################
            for i in range(params.spectra_number_to_process):
                print (params.spec_names_nopath_to_process[i])

                if event == "Process all" and not values['one_spec']:
                    params = replace(params, prev_spec = params.spec_names[i])
                    params = replace(params, prev_spec_nopath = os.path.splitext(params.spec_names_nopath[i])[0])

                # READ THE SPECTRA
                params = replace(params, **dict(zip(["wavelength", "flux"], stm.read_spec(params.spec_names_to_process[i], params.lambda_units)[:2])))
                params = replace(params, original_wavelength=params.wavelength, original_flux=params.flux)

                try:
                    params.original_wavelength = params.wavelength
                    params.original_flux = params.flux
                except Exception:
                    if event == 'Process all':
                        print ('Something went wrong')
                    else:
                        sg.popup("You still need to load a valid spectrum. I don't change my mind")
                        continue

                ################################## SPECTRA MANIPULATION TASKS ##########################
                if not params.reorder_op: #without reordering, I execute tasks following the GUI order
                    params = replace( params, reordered_operations = params.active_operations.copy()) #Activates Spectra manipulation tasks in original order

                #performing the spectra manipulation tasks only if I do not perform sum or average of the spectra
                if params.do_nothing:
                    #Cycling to all the active tasks and picks the one activated following the order of the reordered_operations
                    for op_name, op_var in params.reordered_operations:

                        # 1) CROPPING
                        if op_var == "cropping_spectrum":
                            i == 0 and print('\n*** Cropping ***\n')
                            params = apply_spec_tasks.apply_cropping(event, save_plot, params)

                        # 2) DYNAMIC CLEANING
                        elif op_var == "sigma_clipping":
                            i == 0 and print('\n*** Dynamic cleaning ***\n')
                            if not params.sigma_clip_have_file:
                                params = apply_spec_tasks.apply_sigma_clipping(event, save_plot, params)
                            else: #using an external file with R and sigma
                                params = apply_spec_tasks.apply_sigma_clipping_from_file(event, save_plot, params, i)

                        # 3) WAVELET CLEANING
                        elif op_var == "wavelet_cleaning":
                            i == 0 and print('\n*** Wavelet cleaning ***\n')
                            params = apply_spec_tasks.apply_wavelet_cleaning(event, save_plot, params)

                        # 4) DENOISE
                        elif op_var == "filter_denoise":
                            i == 0 and print('\n*** Denoising ***\n')
                            params = apply_spec_tasks.apply_denoising(event, save_plot, params)

                        # 5) DOPPLER CORRECTION
                        elif op_var == "dop_cor":
                            i == 0 and print('\n*** Dopcor/z correazion ***\n')
                            if (params.dop_cor_single_shot):
                                params = apply_spec_tasks.apply_doppler_correction(event, save_plot, params)
                            else: #if I have an external file with dopcor/z values
                                params = apply_spec_tasks.apply_doppler_correction_from_file(event, save_plot, params, i)

                        #6) HELIOCENTRIC CORRECTION
                        elif op_var == "helio_corr":
                            i == 0 and print('\n*** Heliocentric correction ***\n')
                            if (params.helio_single_shot):
                                params = apply_spec_tasks.apply_heliocentric_correction(event, save_plot, params)
                            else: #If I have an external file with heliocentric corrections to apply
                                params = apply_spec_tasks.apply_heliocentric_correction_from_file(event, save_plot, params, i)

                        #7) REBIN
                        elif op_var == "rebinning":
                            i == 0 and print('\n*** Rebinning ***\n')
                            params = apply_spec_tasks.apply_rebinning(event, save_plot, params)

                        # 8) DEGRADE RESOLUTION
                        elif op_var == "degrade":
                            i == 0 and print('\n*** Degrade resolution ***\n')
                            params = apply_spec_tasks.apply_resolution_degradation(event, save_plot, params)

                        # 9) NORMALISE SPECTRUM TO
                        elif op_var == "normalize_wave":
                            i == 0 and print('\n*** Normalize ***\n')
                            params = apply_spec_tasks.apply_normalisation(event, save_plot, params)

                        # 10) SIGMA BROADENING
                        elif op_var == "sigma_broad":
                            i == 0 and print('\n*** Velocity dispersion ***\n')
                            params = apply_spec_tasks.apply_sigma_broadening(event, save_plot, params)

                        # 11) ADD NOISE
                        elif op_var == "add_noise":
                            i == 0 and print('\n*** Add noise ***\n')
                            params = apply_spec_tasks.apply_noise_addition(event, save_plot, params)

                        # 12) CONTINUUM MODELLING
                        elif op_var == "continuum_sub":
                            i == 0 and print('\n*** Continuum modelling ***\n')
                            params = apply_spec_tasks.apply_continuum_subtraction(event, save_plot, params)

                        # 13) SUBTRACT NORMALISED AVERAGE
                        elif op_var == "subtract_normalized_avg":
                            i == 0 and print('\n*** Subtract normalized average ***\n')
                            if not values['one_spec']:
                                params = apply_spec_tasks.apply_subtract_normalised_average(event, save_plot, params)
                            if values['one_spec']:
                                sg.popup('There is no average to subtract!')
                                continue

                        # 14) SUBTRACT NORMALISED SINGLE SPECTRUM
                        elif op_var == "subtract_normalized_spec":
                            i == 0 and print('\n*** Subtract normalized spectrum ***\n')
                            params = apply_spec_tasks.apply_subtract_normalised_spectrum(event, save_plot, params)

                        # 15) ADD CONSTANT (PEDESTAL)
                        elif op_var == "add_pedestal":
                            i == 0 and print('\n*** Add constant ***\n')
                            params = apply_spec_tasks.apply_add_pedestal(event, save_plot, params)

                        # 16) MULTIPLY
                        elif op_var == "multiply":
                            i == 0 and print('\n*** Multiply ***\n')
                            params = apply_spec_tasks.apply_multiplication(event, save_plot, params)

                        # 17) DERIVATIVES
                        elif op_var == "derivatives":
                            i == 0 and print('\n*** Derivatives ***\n')
                            params = apply_spec_tasks.apply_derivatives(event, save_plot, params)

                    #plotting the results
                    if (event == 'Preview spec.'):
                        try:
                            plt.plot(params.original_wavelength, params.original_flux, label = 'Original spec.')
                            plt.plot(params.wavelength, params.flux, label = 'Processed')
                            plt.xlabel('Wavelength (Å)', fontsize = 9)
                            plt.title(params.prev_spec_nopath)
                            plt.ylabel('Flux')
                            plt.legend(fontsize = 10)
                            plt.show()
                            plt.close()
                        except ValueError:
                            print ('Something went wrong, cannot complete the task. Check the spectrum. Tip: it is really a spectrum?')
                            continue

                ################### COMBINE SPECTRA TASKS, not available for 'Process all' #########################
                if (not params.do_nothing or params.use_for_spec_an) and event == 'Process all':
                    i == 0 and sg.popup('Mean and sum of all the spectra require click on process selected')
                    break
                if not params.do_nothing and values['one_spec']:
                    sg.popup ('You just have one spectrum. Cannot do what you want!')
                    continue

                # Apply math combination tasks and discarding the other previous tasks
                if (not params.do_nothing and not values['one_spec']):
                    print ('WARNING: I will discard all the activated tasks to perform this task')
                    params = apply_spec_tasks.combine_spectra(event, save_plot, params)
                    if params.use_for_spec_an: #If I want use the combined spectrum for spectral analysis
                        print ('Using sum or average to spectral analysis')
                        params = replace(params, wavelength=params.proc_wavelength, flux=params.proc_flux)


                #******************************* SPECTRA ANALYSIS TASKS *********************************
                if not event == 'Preview spec.': # Not activated if 'Preview spec.' is pressed

                    #1) BLACKBODY FITTING
                    if (bb_fit):
                        i == 0 and print('\nRunning blackbody fitting task...\n')
                        temperature_bb, residual_bb, T_err, chi2, params = apply_analysis_tasks.apply_blackbody_fitting(event, save_plot, params)
                        #Updating the file
                        if event == 'Process all':
                            try:
                                if temperature_bb is not None:
                                    df_bb.at[i, 'T(K)']= temperature_bb
                                    df_bb.at[i, 'err']= T_err
                                    df_bb.at[i, 'chi2']= chi2
                                    df_bb.to_csv(bb_file, index= False, sep=' ')
                                    i == (params.spectra_number_to_process - 1) and print(f'File saved: {bb_file}\n')
                            except Exception:
                                print ('Cannot write the file')

                    # 2) CROSS-CORRELATION
                    if (cross_corr):
                        i == 0 and print('\nRunning cross-correlation task...\n')
                        value_at_max, error, params = apply_analysis_tasks.apply_cross_correlation(event, save_plot, params)
                        if event == "Process all":
                            # Updating and writing the file
                            file_writer.save_velocity_or_redshift_to_file(i, params, value_at_max, error, df_rv, rv_file)

                    # 3) VELOCITY DISPERSION
                    if (sigma_measurement):
                        i == 0 and print('\nRunning velocity dispersion task...\n')
                        sigma, error, chisqr, params = apply_analysis_tasks.apply_velocity_dispersion(event, save_plot, params)
                        if event == 'Process all':
                            try:
                                ##Updating and writing the file. No need for an external function
                                df_sigma.at[i, 'Sigma(km/s)']= round(sigma,1)
                                df_sigma.at[i, 'err']= round(error,1)
                                df_sigma.to_csv(sigma_file, index= False, sep=' ')
                                i == (params.spectra_number_to_process - 1) and print(f'File saved: {sigma_file}\n')
                            except Exception:
                                print ('Error writing the file')

                    #4a) EQUIVALENT WIDTH, CASE 1: f I want to measure just one index
                    if (ew_measurement and params.single_index):
                        i == 0 and print('\nRunning equivalent width task with a single index...\n')
                        idx, ew, err, snr_ew, ew_mag, err_mag, params = apply_analysis_tasks.apply_ew_measurement_single(event, save_plot, params)
                        if event == 'Process all':
                            # Updating and writing the file
                            file_writer.save_ew_to_file(i, params, ew, err, ew_mag, err_mag, df_ew, ew_file,
                                df_ew_mag, ew_file_mag, df_snr_ew, snr_ew_file, snr_ew)

                    #4b) EQUIVALENT WIDTH, CASE 2: If I have an index list file
                    if (ew_measurement and params.have_index_file):
                        i == 0 and print('\nRunning equivalent width task with an index list...\n')
                        id_array, ew_array, err_array, snr_ew_array, ew_array_mag, err_array_mag, params = apply_analysis_tasks.apply_ew_measurement_list(event, save_plot, params)
                        if event == 'Process all':
                            # Updating and writing the file
                            file_writer.save_ew_indices_to_file(
                                i, params, num_indices, ew_array, err_array, ew_array_mag, err_array_mag,
                                snr_ew_array, df_ew, ew_file, df_ew_mag, ew_file_mag,
                                df_snr_ew, snr_ew_file, ew_id, ew_id_mag, snr_ew_id, spectra_id)

                    #4c) EQUIVALENT WIDTH, CASE 3: If I want to measure the Lick/IDS indices
                    if (ew_measurement and params.lick_ew):
                        i == 0 and print('\nRunning equivalent width task with Lick/IDS indices...\n')
                        lick_id_array, lick_ew_array, lick_err_array, lick_snr_ew_array, lick_ew_array_mag, lick_err_array_mag, age, met, alpha, err_age, err_met, err_alpha, lick_for_ssp, ssp_model, ssp_lick_indices_list, ssp_lick_indices_err_list, params = apply_analysis_tasks.apply_lick_indices_ew_measurement(event, save_plot, i, params)
                        if i == 0:
                            #reading the index file once to retrieve the number of Lick indices in order to fill the file
                            lick_idx_names, lick_indices = ls.read_idx(params.lick_index_file)
                            num_lick_indices = len(lick_idx_names) #19
                        if event == 'Process all':
                            # Updating and writing the file
                            if i == 0:
                                lick_to_plot = [] # define the lists to accomodate the lick indices to be plotted at the end on the index-index grids
                                lick_err_to_plot = []

                            file_writer.save_lick_indices_to_file(
                                i, params, num_lick_indices, lick_ew_array, lick_err_array, lick_ew_array_mag,
                                lick_err_array_mag, lick_snr_ew_array, df_ew_lick, ew_lick_file, df_ew_lick_mag,
                                ew_lick_file_mag, df_snr_lick_ew, snr_lick_ew_file, ew_lick_id, ew_lick_id_mag,
                                snr_lick_ew_id, spectra_lick_id, df_lick_param, ssp_lick_param_file, lick_for_ssp,
                                df_ssp_param, ssp_param_file, age, err_age, met, err_met, alpha, err_alpha, save_plot,
                                ssp_lick_indices_list, ssp_lick_indices_err_list, params.spectra_list_name, params.result_plot_dir,
                                ssp_model, lick_to_plot, lick_err_to_plot)

                    #5) LINE(S) FITTING: CaT
                    if (line_fitting and params.cat_band_fit):
                        i == 0 and print('\nRunning CaT fitting task...\n')
                        min_wave1, min_wave2, min_wave3, residual_wave1, residual_wave2, residual_wave3, ew_array_ca1, ew_array_ca2, ew_array_ca3, real_cat1, real_cat2, real_cat3, delta_rv1, delta_rv2, delta_rv3, sigma_cat1, sigma_cat2, sigma_cat3, sigma_cat1_vel, sigma_cat2_vel, sigma_cat3_vel, params = apply_analysis_tasks.apply_cat_line_fitting(event, save_plot, params)
                        if event == 'Process all':
                        #Updating and writing the file. No need for an external function
                            try:
                                df_fit.at[i, 'ca1_wave']= round(min_wave1,2)
                                df_fit.at[i, 'ca2_wave']= round(min_wave2,2)
                                df_fit.at[i, 'ca3_wave']= round(min_wave3,2)
                                df_fit.at[i, 'dw_ca1']= round(residual_wave1,2)
                                df_fit.at[i, 'dw_ca2']= round(residual_wave2,2)
                                df_fit.at[i, 'dw_ca3']= round(residual_wave3,2)
                                df_fit.at[i, 'ew_ca1']= round(ew_array_ca1,2)
                                df_fit.at[i, 'ew_ca2']= round(ew_array_ca2,2)
                                df_fit.at[i, 'ew_ca3']= round(ew_array_ca3,2)
                                df_fit.to_csv(fit_file, index= False, sep=' ')
                                i == (params.spectra_number_to_process - 1) and print(f'File saved: {fit_file}\n')
                            except Exception:
                                print ('Cannot write the file')

                    # b) LINE(S) FITTING: USER LINE
                    if (line_fitting and not params.cat_band_fit):
                        i == 0 and print('\nRunning line fitting task...\n')
                        min_wave, sigma_line, sigma_line_vel, params = apply_analysis_tasks.apply_line_fitting(event, save_plot, params)
                        if event == 'Process all':
                            #Updating and writing the file
                            try:
                                df_fit.at[i, 'line_wave']= round(min_wave,2)
                                df_fit.to_csv(fit_file, index= False, sep=' ')
                                i == (params.spectra_number_to_process - 1) and print(f'File saved: {fit_file}\n')
                            except Exception:
                                print ('Cannot write the file')

                    # 6) KINEMATICS WITH PPXF
                    if (perform_kinematics):
                        i == 0 and print('\nRunning stars and gas kinematics task...\n')
                        kinematics, error_kinematics, bestfit_flux, bestfit_wavelength, kin_component, kin_gas_component, snr_kin, error_kinematics_mc, kin_gas_names, kin_gas_flux, kin_gas_flux_err, params = apply_analysis_tasks.apply_ppxf_kinematics(event, save_plot, params)
                        if event == 'Process all':
                            #Updating and writing the file(s)
                            df_kin_gas = file_writer.save_kinematics_to_file(i, params, kinematics, error_kinematics, error_kinematics_mc, kin_gas_component, kin_gas_names, kin_gas_flux, kin_gas_flux_err, kin_component, snr_kin, df_kin, kin_file, df_kin_mc, kin_file_mc, df_kin_gas, kin_file_gas)

                    # 7) STELLAR POPULATIONS WITH PPXF
                    if (stellar_pop):
                        i == 0 and print('\nRunning stellar populations and SFH task...\n')
                        kinematics, info_pop, info_pop_mass, mass_light, chi_square, met_err, mass_met_err, snr_pop, ppxf_pop_lg_age, ppxf_pop_lg_met, age_err_abs, mass_age_err_abs, alpha_err, mass_alpha_err, t50_age, t80_age, t50_cosmic, t80_cosmic, ssp_lick_indices_ppxf, ssp_lick_indices_err_ppxf, ppxf_lick_params, params = apply_analysis_tasks.apply_ppxf_stellar_populations(event, save_plot, params)
                        if kinematics is None:
                            print('Kinematics moments are zero, the fit has failed\n')
                        if event == 'Process all':
                            #Updating and writing the file(s)
                            file_writer.save_population_analysis_to_file(
                                i, params, kinematics, info_pop, info_pop_mass, mass_light,
                                chi_square, met_err, mass_met_err, snr_pop, age_err_abs,
                                mass_age_err_abs, alpha_err, mass_alpha_err, t50_age, t80_age, t50_cosmic, t80_cosmic, ssp_lick_indices_ppxf,
                                ssp_lick_indices_err_ppxf, ppxf_lick_params, df_pop, pop_file,
                                df_ssp_param_ppxf, ssp_param_file_ppxf)

                #progress meter and error messages in Process all mode
                if event == "Process all":
                    if not sg.OneLineProgressMeter('Task progress', i+1, params.spectra_number_to_process, 'Processing spectra:', orientation='h',button_color=('white','red')):
                        print ('***CANCELLED***\n')
                        break

                    if (params.save_final_spectra and params.task_spec2 == 1):
                        file_final = params.result_spec+'proc_' + params.prev_spec_nopath + '.fits'
                        uti.save_fits(params.wavelength, params.flux, file_final)
                        #considering also the cont sub task which saves the continuum!
                        if (params.continuum_sub):
                            file_cont = params.result_spec+'cont_' + params.prev_spec_nopath + '.fits'
                            uti.save_fits(params.wavelength, params.continuum_flux, file_cont)
                            print ('File saved: ', file_cont)
                        print(f'File saved: {file_final}\n')

                    elif (params.task_done2 == 0 ):
                        if i == 0:
                            sg.popup('Nothing to process!')
                        if not sg.OneLineProgressMeter('Task progress', i+1, params.spectra_number,  'single', 'Processing spectra:', orientation='h',button_color=('white','red')):
                            break

                # error messages in preview mode
                elif event == 'Process selected' or event == 'Preview result':
                    if (params.task_analysis == 0 and event == 'Preview result'):
                        sg.popup ('No spectral analysis task selected. Nothing to preview!')
                        continue

                    # Save only the final results, without the intermediate files
                    elif (event == 'Process selected' and params.task_done == 0 and params.task_analysis == 0 ):
                        sg.popup ('Nothing to process!')
                        continue
                    if (params.save_final_spectra and event == 'Process selected' and params.task_spec == 1):
                        file_final = params.result_spec+'proc_' + params.prev_spec_nopath + '.fits'
                        uti.save_fits(params.wavelength, params.flux, file_final)

                        #considering also the cont sub task that saves the continuum!
                        if (params.continuum_sub):
                            file_cont = params.result_spec+'cont_' + params.prev_spec_nopath + '.fits'
                            uti.save_fits(params.wavelength, params.continuum_flux, file_cont)
                            print ('File saved: ', file_cont)
                        print(f'File saved: {file_final}\n')
            params = replace(params, kin_stars_templates=None, kin_lam_temp=None, kin_velscale_templates=None)
            
    #************************************** SAVE AND LOAD PARAMETER VALUES *********************************************
        if event == 'Save parameters...':
            # Open a window to select the path to save the file
            filename = sg.popup_get_file('Save file as...', save_as=True, default_extension=".json", file_types=(("JSON Files", "*.json"),))
            if filename:
                try:
                    settings.save_settings(filename, keys, events, values, params)
                    print('User settings saved')
                    sg.popup_ok(f'Configuration file saved:\n{filename}')
                except json.JSONDecodeError:
                    sg.popup_error('Content not valid for JSON.')

        # listing all the parameters to be loaded
        elif event == 'Load parameters...':
            try:
                filename = sg.popup_get_file('Select the file to load...', file_types=(("JSON Files", "*.json"),))
                keys, events, loaded_values, params = settings.load_settings(filename, params)
                values.update(loaded_values)
                #Updating the GUI
                for key, value in loaded_values.items():
                    if key in window.AllKeysDict:
                        window[key].update(value)
                window.refresh()
                sg.Popup('Settings loaded')
                print('Settings loaded')
            except Exception:
                sg.popup('ERROR: Problem loading the parameters')
                print('Settings NOT loaded')

        elif event == 'Restore default parameters':
            try:
                keys, events, loaded_values, params = settings.load_settings(os.path.join(BASE_DIR, "system_files", "default_settings.json"), params)
                values.update(loaded_values)
                #Updating the GUI
                for key, value in loaded_values.items():
                    if key in window.AllKeysDict:
                        window[key].update(value)
                window.refresh()
                sg.Popup('Default parameters restored')
                print('Default parameters restored')
            except Exception:
                sg.Popup('ERROR restoring default parameters')
                print('ERROR restoring default parameters')

    window.close()
    ########################### END OF PROGRAM! ####################################
