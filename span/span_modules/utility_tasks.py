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

# Functions of the Utilities frame

try: #try local import if executed as script
    #GUI import
    from FreeSimpleGUI_local import FreeSimpleGUI as sg
    from span_functions import system_span as stm
    from span_functions import spec_manipul as spman
    from span_functions import utilities as uti

except ModuleNotFoundError: #local import if executed as package
    #GUI import
    from span.FreeSimpleGUI_local import FreeSimpleGUI as sg
    from span.span_functions import system_span as stm
    from span.span_functions import spec_manipul as spman
    from span.span_functions import utilities as uti

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(CURRENT_DIR)

def show_fits_header(prev_spec, layout):

    """
    Displays the FITS header of the selected spectrum in a separate window.

    Returns:
    - None
    """

    print('*** Showing header ***')
    header = uti.show_hdr(prev_spec)

    layout_hdr = [[sg.Multiline(header, size=(120, 30) if layout == "android" else (100, 40),
                                disabled=True, autoscroll=True, key='-MULTILINE-')],
                  [sg.Button('Close')]]

    window_hdr = sg.Window('FITS Header Viewer', layout_hdr, finalize=True)

    while True:
        event_hdr, values_hdr = window_hdr.read()
        if event_hdr in (sg.WIN_CLOSED, 'Close'):
            break

    window_hdr.close()



def show_sampling(wavelength):

    """
    Displays the wavelength sampling step.

    Returns:
    - None

    """

    print('*** Showing sampling ***')
    step_spectrum, is_linear = uti.show_sampling(wavelength)
    sg.popup('Step: ', round(step_spectrum, 4), 'A')



def show_resolution(wavelength, flux, res_wave1, res_wave2):

    """
    Computes and displays the spectral resolution.

    Returns:
    - None
    """

    print('*** Showing the resolution of the spectrum ***')

    if res_wave1 < wavelength[0] or res_wave1 > wavelength[-1] or res_wave2 < wavelength[0] or res_wave2 > wavelength[-1]:
        sg.popup('Wavelength window outside the limits of the spectrum!')
        return

    resolution_R, fwhm, fwhm_err, line_wave, line_flux_spec_norm, line_flux_spec_fit = uti.resolution (wavelength, flux, res_wave1, res_wave2)
    print('Resolution R: ', resolution_R)

    plt.plot(line_wave, line_flux_spec_norm, label='Spectrum')
    plt.plot(line_wave, line_flux_spec_fit, label='Fit line spec')
    plt.xlabel('Wavelength (A)')
    plt.ylabel('Flux')
    plt.title(f'Resolution: {round(resolution_R)}')
    plt.legend()
    plt.show()



def convert_spectrum(wavelength, flux, prev_spec, convert_to_ascii, lambda_units):

    """
    Converts a spectrum to ASCII or FITS format.

    Returns:
    - None
    """
    print('*** Converting spectrum ***')
    type_spec_to_convert = 'ASCII' if convert_to_ascii else 'FITS'

    try:
        uti.convert_spec(wavelength, flux, prev_spec, type_spec_to_convert, lambda_units)
        print(f'Spectrum {os.path.basename(prev_spec)} converted to {type_spec_to_convert} and stored in the same directory.')
    except Exception:
        print('Error')



def compare_spectra(prev_spec, spec_compare_file, lambda_units):

    """
    Reads and compares two spectra, displaying the plots.

    Returns:
    - None
    """

    if not os.path.isfile(spec_compare_file):
        sg.popup('The spectrum to compare does not exist.')
        return

    print('*** Comparing spectra ***')

    try:
        # Read the primary spectrum
        wavelength, flux, _, _ = stm.read_spec(prev_spec, lambda_units)
        wave_compare, flux_compare, _, _ = stm.read_spec(spec_compare_file, lambda_units)

        # Normalisation
        norm_flux = flux/np.median(flux)
        norm_compare_flux = flux_compare/np.median(flux_compare)

        # Plot comparison
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
        fig.suptitle('Compared spectra')
        ax1.plot(wavelength, flux, label='Selected spectrum')
        ax1.plot(wave_compare, flux_compare, label='Comparison spectrum')
        ax1.set_xlabel('Wavelength (A)')
        ax1.set_ylabel('Flux')
        ax1.legend()
        ax2.plot(wavelength, norm_flux, label='Normalised selected spectrum')
        ax2.plot(wave_compare, norm_compare_flux, label='Normalised comparison spectrum')
        ax2.set_xlabel('Wavelength (A)')
        ax2.set_ylabel('Normalised flux')
        ax2.legend()

        plt.show()
    except Exception:
        sg.popup('Cannot compare spectra. Check the comparison spectrum: does it exist and have the same lambda units?')



def convert_flux_task(event, prev_spec, prev_spec_nopath, spec_names, spec_names_nopath,
                      spectra_number, convert_flux, convert_to_flambda, convert_to_fnu,
                      lambda_units, result_spec, result_data, one_spec):
    """
    Convert the flux of one or all spectra and optionally display a plot.

    Returns:
    - None
    """

    try:
        # Determine the type of conversion
        if convert_flux:
            type_to_convert = 'to_flambda' if convert_to_flambda else 'to_fnu' if convert_to_fnu else None

            if event in ['See plot', 'convert_one']:
                # Read the selected spectrum
                wavelength, flux, step, name = stm.read_spec(prev_spec, lambda_units)
                wave_limits = np.array([wavelength[0], wavelength[-1]])

                # Convert flux
                converted_flux = uti.convert_flux(wavelength, flux, prev_spec, type_to_convert, lambda_units)

                # Plot the original and converted spectrum
                if event == 'See plot':
                    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
                    fig.suptitle('Original - Converted Spectrum')
                    ax1.plot(wavelength, flux, label='Original spectrum')
                    ax1.set_xlabel('Wavelength (A)')
                    ax1.set_ylabel('Flux density')
                    ax1.legend(fontsize=10)
                    ax2.plot(wavelength, converted_flux, label='Converted flux')
                    ax2.set_xlabel('Wavelength (A)')
                    ax2.set_ylabel('Flux density')
                    ax2.legend(fontsize=10)
                    plt.show()
                    plt.close()

                # Save converted flux
                if event == 'convert_one':
                    file_converted = result_spec + 'f_converted_' + prev_spec_nopath + '.fits'
                    uti.save_fits(wavelength, converted_flux, file_converted)
                    print('Spectrum converted flux saved in', result_data + '/spec folder')


        # Convert all spectra
        if event == 'convert_all':

            if convert_flux and not one_spec:
                for i in range(spectra_number):
                    # Read the spectrum
                    wavelength, flux, original_step, obj_name = stm.read_spec(spec_names[i], lambda_units)
                    wave_limits = np.array([wavelength[0], wavelength[-1]])

                    # Convert the spectrum flux
                    converted_flux = uti.convert_flux(wavelength, flux, spec_names[i], type_to_convert, lambda_units)

                    # Save the converted spectrum
                    file_converted = result_spec + 'f_converted_' + spec_names_nopath[i] + '.fits'
                    uti.save_fits(wavelength, converted_flux, file_converted)

                print('Spectra converted flux saved in', result_data + '/spec folder')


    except Exception:
        print('Flux conversion failed')



def snr_analysis(event, prev_spec, spec_names, spec_names_nopath, spectra_number,
                 show_snr, snr_wave, epsilon_wave_snr, lambda_units, values,
                 result_snr_dir, spectra_list_name, timestamp):
    """
    Perform SNR analysis on one or multiple spectra, displaying results and saving them if requested.

    Returns:
    - None
    """

    # Show SNR for a single spectrum
    if event == 'Show snr' or event =='Save one':
        # Read the selected spectrum
        wavelength, flux, step, name = stm.read_spec(prev_spec, lambda_units)
        wave_limits = np.array([wavelength[0], wavelength[-1]])

        if show_snr:
            # snr_task = 1
            print('*** Showing the SNR per pixel ***')

            # Validate wavelength range
            if snr_wave < wave_limits[0] or snr_wave > wave_limits[1]:
                sg.popup('Wavelength window outside the limits of the spectrum!')
                return
            if epsilon_wave_snr < 0 or ((snr_wave + epsilon_wave_snr) > wave_limits[1]) or ((snr_wave - epsilon_wave_snr) < wave_limits[0]):
                sg.popup('Wavelength window exceeds the range of the spectrum!')
                return
            if epsilon_wave_snr == 0:
                sg.popup('Wavelength interval cannot be ZERO!')
                return

            # Compute SNR
            try:
                snr_pix, snr_ang = uti.show_snr(wavelength, flux, snr_wave, epsilon_wave_snr)
                if event == 'Show snr':
                    sg.popup(f'SNR per pixel: {int(round(snr_pix))}, SNR per Ångström: {int(round(snr_ang))} at {snr_wave} A')
                if event == 'Save one':
                    print(f'SNR per pixel: {int(round(snr_pix))}, SNR per Ångström: {int(round(snr_ang))} at {snr_wave} A')
            except Exception:
                sg.popup('Failed to calculate the S/N')


    # Save SNR results for all spectra
    if show_snr and event == 'Save all' and not values:
        # snr_task = 1
        snr_pix_array = np.zeros(spectra_number)
        snr_ang_array = np.zeros(spectra_number)

        # Process all spectra
        for i in range(spectra_number):
            wavelength, flux, original_step, obj_name = stm.read_spec(spec_names[i], lambda_units)
            wave_limits = np.array([wavelength[0], wavelength[-1]])

            # Validate wavelength range
            if snr_wave < wave_limits[0] or snr_wave > wave_limits[1]:
                print('*** WARNING *** Wavelength window outside the limits of the spectrum!')
                return
            elif epsilon_wave_snr < 0 or ((snr_wave + epsilon_wave_snr) > wave_limits[1]) or ((snr_wave - epsilon_wave_snr) < wave_limits[0]):
                print('*** WARNING *** Wavelength window exceeds the range of the spectrum!')
                return
            elif epsilon_wave_snr == 0:
                print('Wavelength interval cannot be ZERO!')
                return
            else:
                # Compute SNR
                snr_pix, snr_ang = uti.show_snr(wavelength, flux, snr_wave, epsilon_wave_snr)
                snr_pix_array[i] = int(snr_pix)
                snr_ang_array[i] = int(snr_ang)

        # Save results to files
        file_snr_pix = f"{result_snr_dir}/{spectra_list_name}_SNR_pix_@{snr_wave}A_{timestamp}.dat"
        file_snr_ang = f"{result_snr_dir}/{spectra_list_name}_SNR_ang_@{snr_wave}A_{timestamp}.dat"

        snr_pix_id = ['#Spectrum', f'SNR_per_pix@{snr_wave}A']
        snr_ang_id = ['#Spectrum', f'SNR_per_Ang@{snr_wave}A']

        snr_pix_data_array = np.column_stack((spec_names_nopath, snr_pix_array))
        snr_ang_data_array = np.column_stack((spec_names_nopath, snr_ang_array))

        # Create DataFrame and save to files
        df_snr_pix = pd.DataFrame(snr_pix_data_array, columns=snr_pix_id)
        df_snr_ang = pd.DataFrame(snr_ang_data_array, columns=snr_ang_id)

        df_snr_pix.to_csv(file_snr_pix, index=False, sep=' ')
        df_snr_ang.to_csv(file_snr_ang, index=False, sep=' ')

        print(f'Files {file_snr_pix} and {file_snr_ang} saved')

