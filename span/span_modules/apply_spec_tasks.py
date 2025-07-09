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

# Functions to apply the tasks of the Spectra manipulation panel

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
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.backends.backend_tkagg
from astropy.coordinates import SkyCoord, EarthLocation
import datetime
from dataclasses import replace


CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(CURRENT_DIR)




def apply_cropping(event, save_plot, params):

    """
    Applies the cropping task to a spectrum, limiting its wavelength range

    """

    # Header
    task_done, task_spec, task_done2, task_spec2 = params.task_done, params.task_spec, params.task_done2, params.task_spec2
    if event == 'Process all':
        task_done2, task_spec2 = 1, 1
    else:
        task_done, task_spec = 1, 1

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_spec=task_spec, task_done2=task_done2, task_spec2=task_spec2)

    # Cropping
    try:
        cropping_range = np.array([params.cropping_low_wave, params.cropping_high_wave])
        cropped_wavelength, cropped_flux = spman.crop_spec(params.wavelength, params.flux, cropping_range)

        if len(cropped_wavelength) > 10:
            new_wavelength, new_flux = cropped_wavelength, cropped_flux
        else:
            print("WARNING: The crop window is too small or out of range. Skipping...")
            return params  # Ritorna i params originali senza modifiche

        # Save cropped spectrum if requested
        if params.save_intermediate_spectra and (event in ['Process all', 'Process selected']):
            try:
                file_cropped = os.path.join(params.result_spec, f'crop_{params.prev_spec_nopath}.fits')
                uti.save_fits(new_wavelength, new_flux, file_cropped)
                print("File saved:", file_cropped)
            except ValueError:
                print("Something went wrong, cannot complete the task.")
                return params

        # # Save plot if required # PLOTS  FOR SPECTRAL MANIPULATION PANEL ARE DEACTIVATED BY DEFAULT, BUT IF YOU NEED THEM, JUST UNCOMMENT THESE LINES FOR ALL THE TASKS
        # if event == 'Process all' and save_plot:
        #     plt.plot(params.original_wavelength, params.original_flux, label='Original')
        #     plt.plot(new_wavelength, new_flux, label='Cropped')
        #     plt.xlabel("Wavelength (A)")
        #     plt.ylabel("Flux")
        #     plt.legend(fontsize=10)
        #     plt.title(f'Cropped {params.prev_spec_nopath}')
        #     plt.savefig(os.path.join(params.result_plot_dir, f'cropped_{params.prev_spec_nopath}.png'), dpi=300)
        #     plt.close()

        # Update parameters and return
        return replace(params,
                    wavelength=new_wavelength,
                    flux=new_flux)

    except Exception:
        print('Cropping failed')
        return params



def apply_sigma_clipping(event, save_plot, params):

    """
    Applies sigma clipping to clean the spectrum dynamically

    """

    # Header
    task_done, task_spec, task_done2, task_spec2 = params.task_done, params.task_spec, params.task_done2, params.task_spec2
    if event == 'Process all':
        task_done2, task_spec2 = 1, 1
    else:
        task_done, task_spec = 1, 1

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_spec=task_spec, task_done2=task_done2, task_spec2=task_spec2)

    try:
        # Sigma Clipping
        clip_wavelength, clip_flux = spman.sigma_clip(params.wavelength, params.flux, params.clip_factor,
                                                    params.sigma_clip_resolution, params.sigma_clip_single_value)

        if len(clip_wavelength) < 10: #arbitrary magic number
            print("WARNING: The cleaned spectrum has too few points. Skipping...")
            return params  # Return original params if clipping fails

        # Save clipped spectrum if requested
        if params.save_intermediate_spectra and (event in ['Process all', 'Process selected']):
            file_clipped = os.path.join(params.result_spec, f'clip_{params.prev_spec_nopath}.fits')
            uti.save_fits(clip_wavelength, clip_flux, file_clipped)
            print(f'File saved: {file_clipped}')

        # Save plot if required
        # if event == 'Process all' and save_plot:
        #     plt.plot(params.original_wavelength, params.original_flux, label='Original')
        #     plt.plot(clip_wavelength, clip_flux, label='Cleaned')
        #     plt.xlabel("Wavelength (A)")
        #     plt.ylabel("Flux")
        #     plt.legend(fontsize=10)
        #     plt.title(f'Cleaned {params.prev_spec_nopath}')
        #     plt.savefig(os.path.join(params.result_plot_dir, f'cleaned_{params.prev_spec_nopath}.png'), dpi=300)
        #     plt.close()

        # Update parameters and return
        return replace(params,
                    wavelength=clip_wavelength,
                    flux=clip_flux)

    except Exception:
        print("Sigma clip failed")
        return params



def apply_sigma_clipping_from_file(event, save_plot, params, i):

    """
    Apply sigma clipping using values from an external file

    """

    # Header
    task_done, task_spec, task_done2, task_spec2 = params.task_done, params.task_spec, params.task_done2, params.task_spec2
    if event == 'Process all':
        task_done2, task_spec2 = 1, 1
    else:
        task_done, task_spec = 1, 1

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_spec=task_spec, task_done2=task_done2, task_spec2=task_spec2)

    # Only in Process all mode
    if event != 'Process all':
        sg.popup('This option works only in Process all mode with a list of files')
        return params

    # Check if the file exists
    if not os.path.isfile(params.sigma_clip_sigma_file):
        sg.popup('The sigma clip file does not exist. Skipping...')
        return params

    # Loading the clip file
    try:
        sigma_clip_resolution = np.loadtxt(params.sigma_clip_sigma_file, usecols=[1])
        sigma_clip_vel_value = np.loadtxt(params.sigma_clip_sigma_file, usecols=[2])
    except ValueError:
        sg.popup('Input sigma clip data file is not valid!')
        return params

    # Checking the values of the file
    if len(sigma_clip_vel_value) != params.spectra_number:
        sg.popup('The sigma clip file does not match the spectra list length.')
        return params

    # Clip
    try:
        clip_wavelength, clip_flux = spman.sigma_clip(
            params.wavelength, params.flux, params.clip_factor,
            sigma_clip_resolution[i], sigma_clip_vel_value[i])

        # Saving spectrum
        if params.save_intermediate_spectra and (event == 'Process all' or event == 'Process selected'):
            file_clipped = os.path.join(params.result_spec, f'clip_{params.spec_names_nopath[i]}.fits')
            uti.save_fits(clip_wavelength, clip_flux, file_clipped)
            print(f'File saved: {file_clipped}')

        # Saving plots
        # if event == 'Process all' and save_plot:
        #     plt.title(f'Cleaned {params.spec_names_nopath[i]}')
        #     plt.plot(params.original_wavelength, params.original_flux, label='Original')
        #     plt.plot(clip_wavelength, clip_flux, label='Cleaned')
        #     plt.xlabel("Wavelength (A)")
        #     plt.ylabel("Flux")
        #     plt.legend()
        #     plt.savefig(os.path.join(params.result_plot_dir, f'cleaned_{params.spec_names_nopath[i]}.png'), dpi=300)
        #     plt.close()

        return replace(params,
                    wavelength=clip_wavelength,
                    flux=clip_flux)

    except ValueError:
        print('Sigma clip failed')
        return params



def apply_wavelet_cleaning(event, save_plot, params):

    """
    Applies wavelet-based noise reduction to the spectrum

    """

    # Header
    task_done, task_spec, task_done2, task_spec2 = params.task_done, params.task_spec, params.task_done2, params.task_spec2
    if event == 'Process all':
        task_done2, task_spec2 = 1, 1
    else:
        task_done, task_spec = 1, 1

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_spec=task_spec, task_done2=task_done2, task_spec2=task_spec2)

    # Walelet
    try:
        denoised_flux = spman.wavelet_cleaning(params.wavelength, params.flux, params.sigma_wavelets, params.wavelets_layers)

        # Saving the spectrum
        if params.save_intermediate_spectra and (event == 'Process all' or event == 'Process selected'):
            file_wavelet = os.path.join(params.result_spec, f'wavelet_{params.prev_spec_nopath}.fits')
            uti.save_fits(params.wavelength, denoised_flux, file_wavelet)
            print(f'File saved: {file_wavelet}')

        # Saving plots
        # if event == 'Process all' and save_plot:
        #     plt.title(f'Wavelet Cleaned {params.prev_spec_nopath}')
        #     plt.plot(params.original_wavelength, params.original_flux, label='Original')
        #     plt.plot(params.wavelength, denoised_flux, label='Cleaned')
        #     plt.xlabel("Wavelength (A)")
        #     plt.ylabel("Flux")
        #     plt.legend()
        #     plt.savefig(os.path.join(params.result_plot_dir, f'wavelet_cleaned_{params.prev_spec_nopath}.png'), dpi=300)
        #     plt.close()

        # Updated params
        return replace(params,
                    flux=denoised_flux,
                    task_done=task_done,
                    task_spec=task_spec,
                    task_done2=task_done2,
                    task_spec2=task_spec2)
    except Exception:
        print('Wavelet filtering failed')
        return params



def apply_denoising(event, save_plot, params):

    """
    Applies various denoising techniques to the spectrum

    """

    # Header
    task_done, task_spec, task_done2, task_spec2 = params.task_done, params.task_spec, params.task_done2, params.task_spec2
    if event == 'Process all':
        task_done2, task_spec2 = 1, 1
    else:
        task_done, task_spec = 1, 1

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_spec=task_spec, task_done2=task_done2, task_spec2=task_spec2)

    # Denoising
    try:
        new_flux = params.flux
        if params.moving_average and params.box_moving_avg:
            new_flux = spman.mov_avg(new_flux, params.box_moving_avg_size)
        if params.moving_average and params.gauss_moving_avg:
            new_flux = spman.mov_avg_gauss(params.wavelength, new_flux, params.gauss_moving_avg_kernel)
        if params.low_pass_filter:
            new_flux = spman.lowpass(params.wavelength, new_flux, params.lowpass_cut_off, params.lowpass_order)
        if params.bandpass_filter:
            new_flux = spman.bandpass(params.wavelength, new_flux, params.bandpass_lower_cut_off, params.bandpass_upper_cut_off, params.bandpass_order)

        # Saving the spectrum
        if params.save_intermediate_spectra and (event == 'Process all' or event == 'Process selected'):
            file_denoised = os.path.join(params.result_spec, f'denoised_{params.prev_spec_nopath}.fits')
            uti.save_fits(params.wavelength, new_flux, file_denoised)
            print(f'File saved: {file_denoised}')

        # Saving plots
        # if event == 'Process all' and save_plot:
        #     plt.title(f'Denoised {params.prev_spec_nopath}')
        #     plt.plot(params.original_wavelength, params.original_flux, label='Original')
        #     plt.plot(params.wavelength, new_flux, label='Denoised')
        #     plt.xlabel("Wavelength (A)")
        #     plt.ylabel("Flux")
        #     plt.legend()
        #     plt.savefig(os.path.join(params.result_plot_dir, f'denoised_{params.prev_spec_nopath}.png'), dpi=300)
        #     plt.close()

        return replace(params,
                    flux=new_flux)

    except Exception:
        print("Denoise failed")
        return params



def apply_doppler_correction(event, save_plot, params):

    """
    Applies Doppler correction to the spectrum

    """

    task_done, task_spec, task_done2, task_spec2 = params.task_done, params.task_spec, params.task_done2, params.task_spec2
    if event == 'Process all':
        task_done2, task_spec2 = 1, 1
    else:
        task_done, task_spec = 1, 1

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_spec=task_spec, task_done2=task_done2, task_spec2=task_spec2)

    try:
        # Dopcor function
        new_wavelength, new_flux = spman.dopcor(params.wavelength, params.flux, params.dop_cor_single_shot_vel, params.dop_cor_have_vel)

        # Saving
        if params.save_intermediate_spectra and (event == 'Process all' or event == 'Process selected'):
            file_dopcor = os.path.join(params.result_spec, f'dopcor_{params.prev_spec_nopath}.fits')
            uti.save_fits(new_wavelength, new_flux, file_dopcor)
            print(f'File saved: {file_dopcor}')

        # Salving plot
        # if event == 'Process all' and save_plot:
        #     plt.title(f'Doppler Corrected {params.prev_spec_nopath}')
        #     plt.plot(params.original_wavelength, params.original_flux, label='Original')
        #     plt.plot(new_wavelength, new_flux, label='Dop cor')
        #     plt.xlabel("Wavelength (A)")
        #     plt.ylabel("Flux")
        #     plt.legend()
        #     plt.savefig(os.path.join(params.result_plot_dir, f'dopcor_{params.prev_spec_nopath}.png'), dpi=300)
        #     plt.close()

        # Restituiamo `params` aggiornato
        return replace(params,
                    wavelength=new_wavelength,
                    flux=new_flux)

    except Exception:
        print("Doppler correction failed")
        return params



def apply_doppler_correction_from_file(event, save_plot, params, i):

    """
    Apply Doppler correction using values from an external file

    """

    # Header
    task_done, task_spec, task_done2, task_spec2 = params.task_done, params.task_spec, params.task_done2, params.task_spec2
    if event == 'Process all':
        task_done2, task_spec2 = 1, 1
    else:
        task_done, task_spec = 1, 1

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_spec=task_spec, task_done2=task_done2, task_spec2=task_spec2)

    # This function only works in 'Process all' mode
    if event != 'Process all':
        sg.popup('This option works only in Process all mode with a list of files')
        return params

    # Ensure the Doppler correction file exists
    if not os.path.isfile(params.dop_cor_file):
        sg.popup('The Doppler correction file does not exist. Skipping...')
        return params

    # Attempt to load Doppler correction values
    try:
        dopcor_values = np.loadtxt(params.dop_cor_file, usecols=[1])
    except ValueError:
        sg.popup('Input Doppler correction data file is not valid!')
        return params

    # Check if the file length matches the number of spectra
    if len(dopcor_values) != params.spectra_number:
        sg.popup('The Doppler correction file does not match the spectra list length.')
        return params

    # Apply Doppler correction
    try:
        new_wavelength, new_flux = spman.dopcor(
            params.wavelength, params.flux, dopcor_values[i], params.dop_cor_have_vel)

        # Print for debugging
        print(params.spec_names_nopath[i], dopcor_values[i])

        # Save the corrected spectrum
        if params.save_intermediate_spectra and (event == 'Process all' or event == 'Process selected'):
            file_dopcor = os.path.join(params.result_spec, f'dopcor_{params.spec_names_nopath[i]}.fits')
            uti.save_fits(new_wavelength, new_flux, file_dopcor)
            print(f'File saved: {file_dopcor}')

        # Save plot if required
        # if event == 'Process all' and save_plot:
        #     plt.title(f'Doppler Corrected {params.spec_names_nopath[i]}')
        #     plt.plot(params.original_wavelength, params.original_flux, label='Original')
        #     plt.plot(new_wavelength, new_flux, label='Dop cor')
        #     plt.xlabel("Wavelength (A)")
        #     plt.ylabel("Flux")
        #     plt.legend()
        #     plt.savefig(os.path.join(params.result_plot_dir, f'dopcor_{params.spec_names_nopath[i]}.png'), dpi=300)
        #     plt.close()


        # Return updated params
        return replace(params,
                    wavelength=new_wavelength,
                    flux=new_flux)

    except Exception:
        print("Doppler/z correction failed")
        return params



def apply_heliocentric_correction(event, save_plot, params):

    """
    Applies heliocentric velocity correction

    """

    # Header
    task_done, task_spec, task_done2, task_spec2 = params.task_done, params.task_spec, params.task_done2, params.task_spec2
    if event == 'Process all':
        task_done2, task_spec2 = 1, 1
    else:
        task_done, task_spec = 1, 1

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_spec=task_spec, task_done2=task_done2, task_spec2=task_spec2)

    # Apply heliocentric correction
    try:
        correction, new_wavelength, new_flux = spman.helio_corr(
            params.wavelength, params.flux, params.helio_single_shot_date,
            params.helio_single_shot_location, params.ra_obj, params.dec_obj
        )

        print(f'Heliocentric correction: {correction} km/s')

        # Save the corrected spectrum
        if params.save_intermediate_spectra and (event == 'Process all' or event == 'Process selected'):
            file_helio = os.path.join(params.result_spec, f'helio_{params.prev_spec_nopath}.fits')
            uti.save_fits(new_wavelength, new_flux, file_helio)
            print(f'File saved: {file_helio}')

        # Save plot if required
        # if event == 'Process all' and save_plot:
        #     plt.title(f'Heliocentric Corrected {params.prev_spec_nopath}')
        #     plt.plot(params.original_wavelength, params.original_flux, label='Original')
        #     plt.plot(new_wavelength, new_flux, label='Helio cor')
        #     plt.xlabel("Wavelength (A)")
        #     plt.ylabel("Flux")
        #     plt.legend()
        #     plt.savefig(os.path.join(params.result_plot_dir, f'heliocor_{params.prev_spec_nopath}.png'), dpi=300)
        #     plt.close()

        # Return updated params
        return replace(params,
                    wavelength=new_wavelength,
                    flux=new_flux)

    except Exception:
        print("Heliocentric correction failed")
        return params



def apply_heliocentric_correction_from_file(event, save_plot, params, i):

    """
    Apply heliocentric correction using values from an external file

    """

    # Header
    task_done, task_spec, task_done2, task_spec2 = params.task_done, params.task_spec, params.task_done2, params.task_spec2
    if event == 'Process all':
        task_done2, task_spec2 = 1, 1
    else:
        task_done, task_spec = 1, 1

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_spec=task_spec, task_done2=task_done2, task_spec2=task_spec2)

    # This function only works in 'Process all' mode
    if event != 'Process all':
        sg.popup('This option works only in Process all mode with a list of files')
        return params

    # Ensure the heliocentric correction file exists
    if not os.path.isfile(params.helio_file):
        sg.popup('The heliocentric correction file does not exist. Skipping...')
        return params

    # Attempt to load data from file
    try:
        location, date = np.loadtxt(params.helio_file, dtype=str, usecols=[0, 1]).T
        ra, dec = np.loadtxt(params.helio_file, usecols=[2, 3]).T
    except ValueError:
        sg.popup('Input data in the heliocentric correction file is not valid!')
        return params

    # Validate the file length
    if len(location) != params.spectra_number:
        sg.popup('The heliocentric correction file does not match the spectra list length.')
        return params

    # Validate dates and locations
    date_not_valid = location_not_valid = False
    for s in range(len(location)):
        try:
            datetime.datetime.strptime(date[s], '%Y-%m-%d')
        except ValueError:
            date_not_valid = True
        try:
            EarthLocation.of_site(location[s])
        except Exception:
            location_not_valid = True

    if date_not_valid or location_not_valid:
        sg.popup('Data in the heliocentric correction file is not valid. Please check!')
        return params

    # Apply the heliocentric correction
    try:
        correction, new_wavelength, new_flux = spman.helio_corr(
            params.wavelength, params.flux, date[i],
            location[i], ra[i], dec[i])

        # Output the correction information
        print(params.spec_names_nopath[i], date[i],
              location[i], ra[i], dec[i],
              correction, 'km/s')

        # Save the corrected spectrum
        if params.save_intermediate_spectra and (event == 'Process all' or event == 'Process selected'):
            file_heliocorr = os.path.join(params.result_spec, f'heliocorr_{params.spec_names_nopath[i]}.fits')
            uti.save_fits(new_wavelength, new_flux, file_heliocorr)
            print(f'File saved: {file_heliocorr}')

        # Save plot if required
        # if event == 'Process all' and save_plot:
        #     plt.title(f'Heliocentric Corrected {params.spec_names_nopath[i]}')
        #     plt.plot(params.original_wavelength, params.original_flux, label='Original')
        #     plt.plot(new_wavelength, new_flux, label='Helio cor')
        #     plt.xlabel("Wavelength (A)")
        #     plt.ylabel("Flux")
        #     plt.legend()
        #     plt.savefig(os.path.join(params.result_plot_dir, f'heliocor_{params.spec_names_nopath[i]}.png'), dpi=300)
        #     plt.close()


        # Return updated params
        return replace(params,
                    wavelength=new_wavelength,
                    flux=new_flux)

    except Exception:
        print("Heliocentric correction failed")
        return params



def apply_rebinning(event, save_plot, params):

    """
    Applies rebinning to the spectrum (linear or logarithmic)

    """

    # Header
    task_done, task_spec, task_done2, task_spec2 = params.task_done, params.task_spec, params.task_done2, params.task_spec2
    if event == 'Process all':
        task_done2, task_spec2 = 1, 1
    else:
        task_done, task_spec = 1, 1

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_spec=task_spec, task_done2=task_done2, task_spec2=task_spec2)

    # Apply rebinning
    try:
        if params.rebinning_linear:
            rebinned_wave, rebinned_flux, _ = spman.resample(params.wavelength, params.flux, params.rebin_step_pix)
        elif params.rebinning_log:
            rebinned_wave, rebinned_flux = spman.log_rebin(params.wavelength, params.flux, params.rebin_step_sigma)
        else:
            rebinned_wave, rebinned_flux = params.wavelength, params.flux  # No rebinning applied

        # Save the rebinned spectrum
        if params.save_intermediate_spectra and (event == 'Process all' or event == 'Process selected'):
            file_rebinned = os.path.join(params.result_spec, f'rebinned_{params.prev_spec_nopath}.fits')
            uti.save_fits(rebinned_wave, rebinned_flux, file_rebinned)
            print(f'File saved: {file_rebinned}')

        # Save plot if required
        # if event == 'Process all' and save_plot:
        #     plt.title(f'Rebinned {params.prev_spec_nopath}')
        #     plt.plot(params.original_wavelength, params.original_flux, label='Original')
        #     plt.plot(rebinned_wave, rebinned_flux, label='Binned')
        #     plt.xlabel("Wavelength (A)")
        #     plt.ylabel("Flux")
        #     plt.legend()
        #     plt.savefig(os.path.join(params.result_plot_dir, f'rebin_{params.prev_spec_nopath}.png'), dpi=300)
        #     plt.close()

        # Return updated params
        return replace(params,
                    wavelength=rebinned_wave,
                    flux=rebinned_flux)

    except Exception:
        print("Rebinning failed")
        return params



def apply_resolution_degradation(event, save_plot, params):

    """
    Applies resolution degradation to the spectrum

    """

    # Header
    task_done, task_spec, task_done2, task_spec2 = params.task_done, params.task_spec, params.task_done2, params.task_spec2
    if event == 'Process all':
        task_done2, task_spec2 = 1, 1
    else:
        task_done, task_spec = 1, 1

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_spec=task_spec, task_done2=task_done2, task_spec2=task_spec2)

    try:
        degraded_wave, degraded_flux = params.wavelength, params.flux  # Default, if no degradation is applied

        # Case A: Degrading from R to R
        if params.is_initial_res_r and params.res_degrade_to_r:
            print('*** Degrading resolution (R to R) ***')
            degraded_wave, degraded_flux = spman.degrade(params.wavelength, params.flux, params.initial_res_r, params.final_res_r, True)

            if params.save_intermediate_spectra and (event == 'Process all' or event == 'Process selected'):
                file_degraded = os.path.join(params.result_spec, f'degraded_R{int(round(params.final_res_r))}_{params.prev_spec_nopath}.fits')
                uti.save_fits(degraded_wave, degraded_flux, file_degraded)
                print(f'File saved: {file_degraded}')

        # Case B: Degrading from R to FWHM
        if params.is_initial_res_r and params.res_degrade_to_fwhm:
            print('*** Degrading resolution from R to FWHM ***')
            degraded_wave, degraded_flux = spman.degradeRtoFWHM(params.wavelength, params.flux, params.initial_res_r, params.final_res_r_to_fwhm)

            if params.save_intermediate_spectra and (event == 'Process all' or event == 'Process selected'):
                file_degraded_R_to_FWHM = os.path.join(params.result_spec, f'degraded_FWHM{round(params.final_res_r_to_fwhm, 1)}_{params.prev_spec_nopath}.fits')
                uti.save_fits(degraded_wave, degraded_flux, file_degraded_R_to_FWHM)
                print(f'File saved: {file_degraded_R_to_FWHM}')

        # Case C: Degrading from FWHM to FWHM
        if params.is_initial_res_fwhm:
            print('*** Degrading resolution in FWHM ***')
            degraded_wave, degraded_flux = spman.degrade_lambda(params.wavelength, params.flux, params.initial_res_fwhm, params.final_res_fwhm)

            if params.save_intermediate_spectra and (event == 'Process all' or event == 'Process selected'):
                file_degraded_lambda = os.path.join(params.result_spec, f'degraded_FWHM{round(params.final_res_fwhm, 1)}_{params.prev_spec_nopath}.fits')
                uti.save_fits(degraded_wave, degraded_flux, file_degraded_lambda)
                print(f'File saved: {file_degraded_lambda}')


        # Case D: Degrading MUSE data to a constant FWHM
        if params.res_degrade_muse:
            print ('*** Degrading resolution of MUSE data to a constant FWHM ***')
            degraded_wave, degraded_flux = spman.degrade_muse(params.wavelength, params.flux, params.res_degrade_muse_value)

            if params.save_intermediate_spectra and (event == 'Process all' or event == 'Process selected'):
                file_degraded_res_muse = os.path.join(params.result_spec, f'degraded_FWHM{round(params.res_degrade_muse_value, 1)}_{params.prev_spec_nopath}.fits')
                uti.save_fits(degraded_wave, degraded_flux, file_degraded_res_muse)
                print(f'File saved: {file_degraded_res_muse}')

        # Save plot if required
        # if event == 'Process all' and save_plot:
        #     plt.title(f'Degraded {params.prev_spec_nopath}')
        #     plt.plot(params.original_wavelength, params.original_flux, label='Original')
        #     plt.plot(degraded_wave, degraded_flux, label='Degraded')
        #     plt.xlabel("Wavelength (A)")
        #     plt.ylabel("Flux")
        #     plt.legend()
        #     plt.savefig(os.path.join(params.result_plot_dir, f'degrade_{params.prev_spec_nopath}.png'), dpi=300)
        #     plt.close()

        # Return updated params
        return replace(params,
                    wavelength=degraded_wave,
                    flux=degraded_flux)

    except Exception:
        print("Degrade resolution failed")
        return params



def apply_normalisation(event, save_plot, params):

    """
    Normalises the spectrum at a given wavelength

    """

    # Header
    task_done, task_spec, task_done2, task_spec2 = params.task_done, params.task_spec, params.task_done2, params.task_spec2
    if event == 'Process all':
        task_done2, task_spec2 = 1, 1
    else:
        task_done, task_spec = 1, 1

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_spec=task_spec, task_done2=task_done2, task_spec2=task_spec2)

    if params.norm_wave < params.wave_limits[0] or params.norm_wave > params.wave_limits[1]:
        error_msg = 'ERROR: Normalisation wavelength exceeds the range of the spectrum!'
        print(error_msg) if event == 'Process all' else sg.popup(error_msg)
        return params

    try:
        step = params.wavelength[1] - params.wavelength[0]
        epsilon_norm = step * 10  # Averaging over 10 wavelength steps
        normalised_flux = spman.norm_spec(params.wavelength, params.flux, params.norm_wave, epsilon_norm, params.flux)

        # Save the normalised spectrum
        if params.save_intermediate_spectra and (event == 'Process all' or event == 'Process selected'):
            file_normalised = os.path.join(params.result_spec, f'norm{params.norm_wave}_{params.prev_spec_nopath}.fits')
            uti.save_fits(params.wavelength, normalised_flux, file_normalised)
            print(f'File saved: {file_normalised}')

        # Save plot if required
        # if event == 'Process all' and save_plot:
        #     plt.title(f'Normalised {params.prev_spec_nopath}')
        #     plt.plot(params.original_wavelength, params.original_flux, label='Original')
        #     plt.plot(params.wavelength, normalised_flux, label='Normalised')
        #     plt.xlabel("Wavelength (A)")
        #     plt.ylabel("Flux")
        #     plt.legend()
        #     plt.savefig(os.path.join(params.result_plot_dir, f'norm_{params.prev_spec_nopath}.png'), dpi=300)
        #     plt.close()

        # Return updated params
        return replace(params,
                    flux=normalised_flux)

    except Exception:
        print("Normalisation failed")
        return params



def apply_sigma_broadening(event, save_plot, params):

    """
    Broadens the spectrum by a given sigma value

    """

    # Header
    task_done, task_spec, task_done2, task_spec2 = params.task_done, params.task_spec, params.task_done2, params.task_spec2
    if event == 'Process all':
        task_done2, task_spec2 = 1, 1
    else:
        task_done, task_spec = 1, 1

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_spec=task_spec, task_done2=task_done2, task_spec2=task_spec2)

    try:
        broadened_flux = spman.sigma_broad(params.wavelength, params.flux, params.sigma_to_add)

        # Save the broadened spectrum
        if params.save_intermediate_spectra and (event == 'Process all' or event == 'Process selected'):
            file_broad = os.path.join(params.result_spec, f'broad{int(round(params.sigma_to_add))}_{params.prev_spec_nopath}.fits')
            uti.save_fits(params.wavelength, broadened_flux, file_broad)
            print(f'File saved: {file_broad}')

        # Save plot if required
        # if event == 'Process all' and save_plot:
        #     plt.title(f'Sigma broad {params.prev_spec_nopath}')
        #     plt.plot(params.original_wavelength, params.original_flux, label='Original')
        #     plt.plot(params.wavelength, broadened_flux, label='Broadened')
        #     plt.xlabel("Wavelength (A)")
        #     plt.ylabel("Flux")
        #     plt.legend()
        #     plt.savefig(os.path.join(params.result_plot_dir, f'broad_{params.prev_spec_nopath}.png'), dpi=300)
        #     plt.close()

        # Return updated params
        return replace(params,
                    flux=broadened_flux)

    except Exception:
        print("Velocity broadening failed")
        return params



def apply_noise_addition(event, save_plot, params):

    """
    Adds noise to the spectrum

    """

    # Header
    task_done, task_spec, task_done2, task_spec2 = params.task_done, params.task_spec, params.task_done2, params.task_spec2
    if event == 'Process all':
        task_done2, task_spec2 = 1, 1
    else:
        task_done, task_spec = 1, 1

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_spec=task_spec, task_done2=task_done2, task_spec2=task_spec2)

    try:
        noisy_flux = spman.add_noise(params.wavelength, params.flux, params.noise_to_add)

        # Save the noisy spectrum
        if params.save_intermediate_spectra and (event == 'Process all' or event == 'Process selected'):
            file_noise = os.path.join(params.result_spec, f'SNR{int(round(params.noise_to_add))}_{params.prev_spec_nopath}.fits')
            uti.save_fits(params.wavelength, noisy_flux, file_noise)
            print(f'File saved: {file_noise}')

        # Save plot if required
        # if event == 'Process all' and save_plot:
        #     plt.title(f'Add noise {params.prev_spec_nopath}')
        #     plt.plot(params.original_wavelength, params.original_flux, label='Original')
        #     plt.plot(params.wavelength, noisy_flux, label='Noisy')
        #     plt.xlabel("Wavelength (A)")
        #     plt.ylabel("Flux")
        #     plt.legend()
        #     plt.savefig(os.path.join(params.result_plot_dir, f'SNR_{params.prev_spec_nopath}.png'), dpi=300)
        #     plt.close()

        # Return updated params
        return replace(params,
                    flux=noisy_flux)

    except Exception:
        print("Add noise failed")
        return params



def apply_continuum_subtraction(event, save_plot, params):

    """
    Performs continuum modelling and subtraction

    """

    # Header
    task_done, task_spec, task_done2, task_spec2 = params.task_done, params.task_spec, params.task_done2, params.task_spec2
    if event == 'Process all':
        task_done2, task_spec2 = 1, 1
    else:
        task_done, task_spec = 1, 1

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_spec=task_spec, task_done2=task_done2, task_spec2=task_spec2)

    # Perform continuum subtraction based on selected method
    try:
        if params.cont_model_filtering:
            corrected_flux, continuum_flux = spman.sub_cont(params.wavelength, params.flux, params.cont_math_operation)
        elif params.cont_model_poly:
            preview = event == 'Preview spec.'
            corrected_flux, continuum_flux = spman.continuum(
                params.wavelength, params.flux, params.cont_want_to_mask, params.cont_mask_ranges,
                params.cont_poly_degree, params.cont_math_operation, preview)

        # Save the modified spectrum
        if params.save_intermediate_spectra and (event == 'Process all' or event == 'Process selected'):
            file_cont_sub = os.path.join(params.result_spec, f'cont_sub_{params.prev_spec_nopath}.fits')
            file_cont = os.path.join(params.result_spec, f'cont_{params.prev_spec_nopath}.fits')
            uti.save_fits(params.wavelength, corrected_flux, file_cont_sub)
            uti.save_fits(params.wavelength, continuum_flux, file_cont)
            print(f'Files saved: {file_cont_sub}, {file_cont}')

        # Save plot if required
        # if event == 'Process all' and save_plot:
        #     plt.title(f'Continuum {params.prev_spec_nopath}')
        #     plt.plot(params.original_wavelength, params.original_flux, label='Original')
        #     plt.plot(params.wavelength, corrected_flux, label='Cont. removed')
        #     plt.xlabel("Wavelength (A)")
        #     plt.ylabel("Flux")
        #     plt.legend()
        #     plt.savefig(os.path.join(params.result_plot_dir, f'cont_{params.prev_spec_nopath}.png'), dpi=300)
        #     plt.close()


        # Return updated params
        return replace(params,
                    continuum_flux = continuum_flux,
                    flux=corrected_flux)

    except Exception:
        print("Parameters not valid")
        return params



def apply_subtract_normalised_average(event, save_plot, params):

    """
    Subtracts the normalised average spectrum

    """

    # Header
    task_done, task_spec, task_done2, task_spec2 = params.task_done, params.task_spec, params.task_done2, params.task_spec2
    if event == 'Process all':
        task_done2, task_spec2 = 1, 1
    else:
        task_done, task_spec = 1, 1

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_spec=task_spec, task_done2=task_done2, task_spec2=task_spec2)

    try:
        # Perform subtraction of the normalised average spectrum
        subtracted_flux = spmt.sub_norm_avg(params.wavelength, params.flux, params.lambda_units,
                                            params.spectra_number, params.spec_names)

        # Save the modified spectrum
        if params.save_intermediate_spectra and (event == 'Process all' or event == 'Process selected'):
            file_subtracted_avg = os.path.join(params.result_spec, f'subtracted_average_{params.prev_spec_nopath}.fits')
            uti.save_fits(params.wavelength, subtracted_flux, file_subtracted_avg)
            print(f'File saved: {file_subtracted_avg}')

        # Save plot if required
        # if event == 'Process all' and save_plot:
        #     plt.title(f'Subtract average {params.prev_spec_nopath}')
        #     plt.plot(params.original_wavelength, params.original_flux, label='Original')
        #     plt.plot(params.wavelength, subtracted_flux, label='Subtracted')
        #     plt.xlabel("Wavelength (A)")
        #     plt.ylabel("Flux")
        #     plt.legend()
        #     plt.savefig(os.path.join(params.result_plot_dir, f'sub_avg_{params.prev_spec_nopath}.png'), dpi=300)
        #     plt.close()


        # Return updated params
        return replace(params,
                    flux=subtracted_flux)

    except Exception:
        print("Failed")
        return params



def apply_subtract_normalised_spectrum(event, save_plot, params):

    """
    Subtracts a single normalised spectrum

    """

    # Header
    task_done, task_spec, task_done2, task_spec2 = params.task_done, params.task_spec, params.task_done2, params.task_spec2
    if event == 'Process all':
        task_done2, task_spec2 = 1, 1
    else:
        task_done, task_spec = 1, 1

    if not os.path.isfile(params.spectra_to_subtract):
        print('ERROR: The spectrum to subtract does not exist.')
        return params

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_spec=task_spec, task_done2=task_done2, task_spec2=task_spec2)

    try:
        updated_flux = spmt.sub_norm_single(params.wavelength, params.flux, params.spectra_to_subtract, params.lambda_units)

        # Save the modified spectrum
        if params.save_intermediate_spectra and (event == 'Process all' or event == 'Process selected'):
            file_subtracted_single = os.path.join(params.result_spec, f'subtracted_single_{params.prev_spec_nopath}.fits')
            uti.save_fits(params.wavelength, updated_flux, file_subtracted_single)
            print(f'File saved: {file_subtracted_single}')

        # Save plot if required
        # if event == 'Process all' and save_plot:
        #     plt.title(f'Subtract norm spec {params.prev_spec_nopath}')
        #     plt.plot(params.original_wavelength, params.original_flux, label='Original')
        #     plt.plot(params.wavelength, updated_flux, label='Subtracted')
        #     plt.xlabel("Wavelength (A)")
        #     plt.ylabel("Flux")
        #     plt.legend()
        #     plt.savefig(os.path.join(params.result_plot_dir, f'sub_norm_spec_{params.prev_spec_nopath}.png'), dpi=300)
        #     plt.close()

        # Return updated params
        return replace(params,
                    flux=updated_flux)

    except Exception:
        print("Subtract norm spec failed")
        return params



def apply_add_pedestal(event, save_plot, params):

    """
    Adds a constant pedestal value to the spectrum

    """

    # Header
    task_done, task_spec, task_done2, task_spec2 = params.task_done, params.task_spec, params.task_done2, params.task_spec2
    if event == 'Process all':
        task_done2, task_spec2 = 1, 1
    else:
        task_done, task_spec = 1, 1

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_spec=task_spec, task_done2=task_done2, task_spec2=task_spec2)

    try:
        updated_flux = params.flux + params.pedestal_to_add

        # Save the modified spectrum
        if params.save_intermediate_spectra and (event == 'Process all' or event == 'Process selected'):
            pedestal_suffix = str(int(round(params.pedestal_to_add)))
            file_pedestal = os.path.join(params.result_spec, f'pedestal{pedestal_suffix}_{params.prev_spec_nopath}.fits')
            uti.save_fits(params.wavelength, updated_flux, file_pedestal)
            print(f'File saved: {file_pedestal}')

        # Save plot if required
        # if event == 'Process all' and save_plot:
        #     plt.title(f'Pedestal {params.prev_spec_nopath}')
        #     plt.plot(params.original_wavelength, params.original_flux, label='Original')
        #     plt.plot(params.wavelength, updated_flux, label='Pedestal')
        #     plt.xlabel("Wavelength (A)")
        #     plt.ylabel("Flux")
        #     plt.legend()
        #     plt.savefig(os.path.join(params.result_plot_dir, f'pedestal_{params.prev_spec_nopath}.png'), dpi=300)
        #     plt.close()

        # Return updated params
        return replace(params,
                    flux=updated_flux)

    except Exception:
        print("Add pedestal failed")
        return params



def apply_multiplication(event, save_plot, params):

    """
    Multiplies the spectrum by a constant factor

    """

    # Header
    task_done, task_spec, task_done2, task_spec2 = params.task_done, params.task_spec, params.task_done2, params.task_spec2
    if event == 'Process all':
        task_done2, task_spec2 = 1, 1
    else:
        task_done, task_spec = 1, 1

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_spec=task_spec, task_done2=task_done2, task_spec2=task_spec2)

    try:
        updated_flux = params.flux * params.multiply_factor

        # Save the modified spectrum
        if params.save_intermediate_spectra and (event == 'Process all' or event == 'Process selected'):
            multiplied_suffix = str(int(round(params.multiply_factor)))
            file_multiplied = os.path.join(params.result_spec, f'multiplied{multiplied_suffix}_{params.prev_spec_nopath}.fits')
            uti.save_fits(params.wavelength, updated_flux, file_multiplied)
            print(f'File saved: {file_multiplied}')

        # Save plot if required
        # if event == 'Process all' and save_plot:
        #     plt.title(f'Multiply {params.prev_spec_nopath}')
        #     plt.plot(params.original_wavelength, params.original_flux, label='Original')
        #     plt.plot(params.wavelength, updated_flux, label='Multiplied')
        #     plt.xlabel("Wavelength (A)")
        #     plt.ylabel("Flux")
        #     plt.legend()
        #     plt.savefig(os.path.join(params.result_plot_dir, f'multiply_{params.prev_spec_nopath}.png'), dpi=300)
        #     plt.close()

        # Return updated params
        return replace(params,
                    flux=updated_flux)

    except Exception:
        print("Multiply failed")
        return params



def apply_derivatives(event, save_plot, params):

    """
    Computes the first and second derivatives of the spectrum

    """

    # Header
    task_done, task_spec, task_done2, task_spec2 = params.task_done, params.task_spec, params.task_done2, params.task_spec2
    if event == 'Process all':
        task_done2, task_spec2 = 1, 1
    else:
        task_done, task_spec = 1, 1
        print('WARNING: The derivatives are NOT used for spectral analysis tasks')

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_spec=task_spec, task_done2=task_done2, task_spec2=task_spec2)

    try:
        first_derivative = np.gradient(params.flux, params.wavelength)
        second_derivative = np.gradient(first_derivative, params.wavelength)

        # Plot results if preview is requested
        if event == 'Preview spec.' or (event == 'Process all' and save_plot):
            fig, axes = plt.subplots(3, 1, figsize=(12, 8))
            fig.suptitle('Spectra derivatives')

            axes[0].plot(params.wavelength, params.flux, label='Original spectrum')
            axes[0].set_ylabel('Flux')
            axes[0].legend()

            axes[1].plot(params.wavelength, first_derivative, label='First derivative')
            axes[1].set_ylabel('Flux')
            axes[1].legend()

            axes[2].plot(params.wavelength, second_derivative, label='Second derivative')
            axes[2].set_ylabel('Flux')
            axes[2].legend()

            plt.xlabel('Wavelength ($\AA$)')
            plt.tight_layout()
            if event == 'Preview spec.':
                plt.show()
            else:
                plt.savefig(os.path.join(params.result_plot_dir, f'derivatives_{params.prev_spec_nopath}.png'), dpi=300)
                plt.close()

        # Save derivative spectra if required
        if event in ('Process selected', 'Process all'):
            file_first_derivative = os.path.join(params.result_spec, f'first_deriv_{params.prev_spec_nopath}.fits')
            file_second_derivative = os.path.join(params.result_spec, f'second_deriv_{params.prev_spec_nopath}.fits')
            uti.save_fits(params.wavelength, first_derivative, file_first_derivative)
            uti.save_fits(params.wavelength, second_derivative, file_second_derivative)
            print(f'Derivative spectra saved: {file_first_derivative}, {file_second_derivative}')

        # Return updated params
        return params

    except Exception:
        print(f'Cannot compute the derivatives. Error: {e}')
        return params



def combine_spectra(event, save_plot, params):

    """
    Combine multiple spectra using averaging, summation, and normalisation

    """

    # Header
    task_done, task_spec, task_done2, task_spec2 = params.task_done, params.task_spec, params.task_done2, params.task_spec2
    if event == 'Process all':
        task_done2, task_spec2 = 1, 1
    else:
        task_done, task_spec = 1, 1

    proc_wavelength, proc_flux = None, None

    #updating the params that changed, that is just the check conditions
    params = replace(params, task_done=task_done, task_spec=task_spec, task_done2=task_done2, task_spec2=task_spec2)

    try:
        # Average all spectra
        if params.average_all:
            print('*** Averaging all spectra ***')
            average_spec = spmt.average(params.lambda_units, params.spectra_number, params.spec_names)
            proc_wavelength, proc_flux = average_spec[:, 0], average_spec[:, 1]

            if event == 'Process selected' and params.save_intermediate_spectra:
                file_avg = os.path.join(params.result_spec, 'avg_spectra.fits')
                uti.save_fits(proc_wavelength, proc_flux, file_avg)
                print(f'File saved: {file_avg}')

        # Normalise and average
        if params.norm_and_average:
            print('*** Normalising and averaging all spectra ***')
            average_norm_spec = spmt.average_norm(params.lambda_units, params.wavelength, params.flux, params.spectra_number, params.spec_names)
            proc_wavelength, proc_flux = average_norm_spec[:, 0], average_norm_spec[:, 1]

            if event == 'Process selected' and params.save_intermediate_spectra:
                file_avg_norm = os.path.join(params.result_spec, 'norm_avg_spectra.fits')
                uti.save_fits(proc_wavelength, proc_flux, file_avg_norm)
                print(f'File saved: {file_avg_norm}')

        # Sum all spectra
        if params.sum_all:
            print('*** Summing all spectra ***')
            sum_spec = spmt.sum_spec(params.lambda_units, params.spectra_number, params.spec_names)
            proc_wavelength, proc_flux = sum_spec[:, 0], sum_spec[:, 1]

            if event == 'Process selected' and params.save_intermediate_spectra:
                file_sum = os.path.join(params.result_spec, 'sum_spectra.fits')
                uti.save_fits(proc_wavelength, proc_flux, file_sum)
                print(f'File saved: {file_sum}')

        # Normalise and sum all spectra
        if params.normalize_and_sum_all:
            print('*** Normalising and summing all spectra ***')
            sum_norm_spec = spmt.sum_norm_spec(params.lambda_units, params.spectra_number, params.spec_names)
            proc_wavelength, proc_flux = sum_norm_spec[:, 0], sum_norm_spec[:, 1]

            if event == 'Process selected' and params.save_intermediate_spectra:
                file_sum_norm = os.path.join(params.result_spec, 'norm_sum_spectra.fits')
                uti.save_fits(proc_wavelength, proc_flux, file_sum_norm)
                print(f'File saved: {file_sum_norm}')

        # Plot results
        if event == 'Preview spec.' and proc_wavelength is not None and proc_flux is not None:
            plt.plot(params.original_wavelength, params.original_flux, label='Original spec.')
            plt.plot(proc_wavelength, proc_flux, label='Processed')
            plt.xlabel('Wavelength ($\AA$)', fontsize=9)
            plt.title('Spectral Combination')
            plt.ylabel('Flux')
            plt.legend(fontsize=10)
            plt.show()
            plt.close()

        # Return updated parameters
        return replace(params,
                        proc_wavelength=proc_wavelength,
                        proc_flux=proc_flux)

    except Exception:
        print('Spectra combination failed')
        return params

