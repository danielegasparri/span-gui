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
######################## THE FOLLOWING FUNCTIONS HAVE BEEN INSPIRED BY THE GIST PIPELINE OF BITTNER ET AL., 2019 #########################
############################################# A special thanks to Adrian Bittner ########################################################


#Functions to bin and extract 1D spectra from datacubes, using the GIST pipeline logic.
#The results are fully compatible with the GIST pipeline.

import os
import shutil
import numpy as np
from astropy.io import fits
import sys
import importlib.util
import functools
from vorbin.voronoi_2d_binning import voronoi_2d_binning
import scipy.spatial.distance as dist
import matplotlib.pyplot as plt




#function to create the dictionary (config) following the GIST standard to be passed to the following functions
def buildConfigFromGUI(ifs_run_id, ifs_input, ifs_output, ifs_redshift,
                       ifs_lfs_data, ifs_ow_config, ifs_ow_output,
                       ifs_routine_read, ifs_origin, ifs_lmin_tot, ifs_lmax_tot,
                       ifs_lmin_snr, ifs_lmax_snr, ifs_min_snr_mask,
                       ifs_mask, ifs_bin_method, ifs_target_snr, ifs_covariance,
                       ifs_prepare_method):

    """
    Returns a `configs` dictionary to be read from the following functions of the module

    """

    configs = {
        "INFO": {
            "RUN_NAME": ifs_run_id,
            "INPUT": ifs_input,
            "OUTPUT": ifs_output,
            "REDSHIFT": ifs_redshift,
            "LSF_DATA": ifs_lfs_data,
            "OW_CONFIG": ifs_ow_config,
            "OW_OUTPUT": ifs_ow_output
        },
        "READ": {
            "ROUTINE": ifs_routine_read,
            "ORIGIN": ifs_origin,
            "LMIN_TOT": ifs_lmin_tot,
            "LMAX_TOT": ifs_lmax_tot,
            "LMIN_SNR": ifs_lmin_snr,
            "LMAX_SNR": ifs_lmax_snr
        },
        "MASKING": {
            "MASK_SNR": ifs_min_snr_mask,
            "MASK": ifs_mask
        },
        "BINNING": {
            "VORONOI": ifs_bin_method,
            "TARGET_SNR": ifs_target_snr,
            "COVARIANCE": ifs_covariance
        },
        "EXTRACTING": {
            "MODE": ifs_prepare_method
        }
    }

    return configs


def reading_data(config):

    """
    Reads the datacube using the specified method from the configuration.

    Parameters:
        config (dict): Configuration dictionary containing method and input details.

    Returns:
        cube (object): The loaded datacube, or "SKIP" in case of failure.
    """

    print("Step 1: Reading the datacube")


    method = config.get('READ', {}).get('ROUTINE', '')
    method_nopath = os.path.splitext(os.path.basename(method))[0]

    if not method:
        print("No read-in method specified.")
        return "SKIP"
    try:
        spec = importlib.util.spec_from_file_location("", method)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        print(f"Using the read-in routine for {method_nopath}")
        return module.read_cube(config)
    except Exception as e:
        print(f"Failed to import or execute the read-in routine {method_nopath}: {e}")
        return "SKIP"


def masking(config, cube, preview, manual_bin, existing_bin):

    """
    Applies a spatial mask to the datacube if required.

    Parameters:
        config (dict): Configuration dictionary.
        cube (object): The loaded datacube.
        preview (bool): If True, performs a preview without saving.

    Returns:
        None
    """

    print("\nStep 2: Applying masking, if any")

    output_file = os.path.join(config['INFO']['OUTPUT'], f"{config['INFO']['RUN_NAME']}_mask.fits")

    if (not preview and os.path.isfile(output_file) and not config['INFO'].get('OW_OUTPUT', False) and manual_bin) or existing_bin:
        print("Masking results already exist. Skipping step.")
        return

    generate_and_apply_mask(config, cube)


def binning(config, cube, preview, voronoi, manual_bin, existing_bin):

    """
    Applies spatial binning to the datacube.

    Parameters:
        config (dict): Configuration dictionary.
        cube (object): The loaded datacube.
        preview (bool): If True, performs a preview without saving.
        voronoi (bool): If True, uses Voronoi binning.

    Returns:
        None or "SKIP" in case of failure.
    """

    print("\nStep 3: Applying binning")

    output_file = os.path.join(config['INFO']['OUTPUT'], f"{config['INFO']['RUN_NAME']}_table.fits")

    if (not preview and (os.path.isfile(output_file) and not config['INFO']['OW_OUTPUT'] and manual_bin)) or existing_bin:
        if not existing_bin:
            print("Results of the module are already in the output directory. Module is skipped.")
            return
        if existing_bin:
            print('Using user mask and bin info')
            return

    try:
        generate_bins(config, cube, voronoi)
    except Exception as e:
        print(f"Spatial binning routine {config.get('BINNING', {}).get('VORONOI', 'UNKNOWN')} failed: {e}")
        return "SKIP"


def save_spectra(config, cube, preview, existing_bin):

    """
    Extracts and saves 1D spectra from the datacube.

    Parameters:
        config (dict): Configuration dictionary.
        cube (object): The loaded datacube.
        preview (bool): If True, performs a preview without saving.

    Returns:
        None or "SKIP" in case of failure.
    """

    print("\nStep 4: Saving the extracted 1D spectra")

    output_prefix = os.path.join(config['INFO']['OUTPUT'], config['INFO']['RUN_NAME'])
    output_file = f"{output_prefix}_BinSpectra_linear.fits"

    if (not preview and os.path.isfile(output_file) and not config['INFO'].get('OW_OUTPUT', False)) and not existing_bin:
        print("Spectra extraction results already exist. Skipping step.")
        return

    try:
        prepare_mask_bin(config, cube, preview)
    except Exception as e:
        print(f"Spectra preparation routine {config.get('EXTRACTING', {}).get('MODE', 'UNKNOWN')} failed: {e}")
        return "SKIP"



###############################################################################
################ FUNCTIONS TO PERFORM THE 4 STEPS ABOVE #############

def generate_and_apply_mask(config, cube):

    """
    Creates a combined mask for the datacube, masking defunct spaxels,
    spaxels below a SNR threshold, and spaxels from an external mask file.

    Parameters:
        config (dict): Configuration dictionary.
        cube (dict): Datacube containing spectral and SNR data.

    Returns:
        None (Saves the mask to a FITS file)
    """

    print("Generating spatial mask...")

    # Mask spaxels that contain NaN values or have a non-positive median flux
    spec = cube['spec']
    median_flux = np.nanmedian(spec, axis=0)
    masked_defunct = np.logical_or(np.any(np.isnan(spec), axis=0), median_flux <= 0)
    print(f"Masking defunct spaxels: {np.sum(masked_defunct)} spaxels are rejected.")

    # Mask spaxels based on signal-to-noise ratio (SNR) threshold
    masked_snr = mask_snr(cube['snr'], cube['signal'], config['MASKING']['MASK_SNR'])

    # Mask spaxels based on an external mask file
    mask_filename = config['MASKING'].get('MASK')
    if mask_filename:
        mask_path = os.path.join(os.path.dirname(config['INFO']['INPUT']), mask_filename)

        if os.path.isfile(mask_path):
            mask_data = fits.getdata(mask_path, ext=1).flatten()
            masked_mask = mask_data == 1
            print(f"Masking spaxels according to mask file: {np.sum(masked_mask)} spaxels are rejected.")
        else:
            print(f"Mask file not found: {mask_path}")
            masked_mask = np.zeros(len(cube['snr']), dtype=bool)
    else:
        print("No mask file provided.")
        masked_mask = np.zeros(len(cube['snr']), dtype=bool)

    # Combine all masks
    combined_mask = np.logical_or.reduce((masked_defunct, masked_snr, masked_mask))

    # Save final mask
    save_mask(combined_mask, masked_defunct, masked_snr, masked_mask, config)


def save_mask(combined_mask, masked_defunct, masked_snr, masked_mask, config):

    """
    Saves the final combined mask and its components to a FITS file.

    Parameters:
        combined_mask (np.ndarray): Boolean array of the final combined mask.
        masked_defunct (np.ndarray): Boolean array for defunct spaxels.
        masked_snr (np.ndarray): Boolean array for SNR-masked spaxels.
        masked_mask (np.ndarray): Boolean array for external mask file spaxels.
        config (dict): Configuration dictionary.

    Returns:
        None (Writes the mask to a FITS file)
    """

    output_file = os.path.join(config['INFO']['OUTPUT'], f"{config['INFO']['RUN_NAME']}_mask.fits")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    print(f"Writing mask file: {output_file}")

    # reading and writing the fits
    with fits.HDUList([
        fits.PrimaryHDU(),
        fits.BinTableHDU.from_columns([
            fits.Column(name='MASK', format='I', array=combined_mask.astype(int)),
            fits.Column(name='MASK_DEFUNCT', format='I', array=masked_defunct.astype(int)),
            fits.Column(name='MASK_SNR', format='I', array=masked_snr.astype(int)),
            fits.Column(name='MASK_FILE', format='I', array=masked_mask.astype(int))
        ], name="MASKFILE")
    ]) as hdul:
        # Comments in the header
        hdul[1].header['COMMENT'] = "Value 0 -> Unmasked"
        hdul[1].header['COMMENT'] = "Value 1 -> Masked"

        hdul.writeto(output_file, overwrite=True)

    print(f"Mask file saved successfully: {output_file}")


def mask_snr(snr, signal, min_snr):

    """
    Masks spaxels based on a minimum SNR threshold.

    Parameters:
        snr (np.ndarray): Array of signal-to-noise ratios for each spaxel.
        signal (np.ndarray): Array of signal values for each spaxel.
        min_snr (float): Minimum SNR threshold for masking.

    Returns:
        masked (np.ndarray): Boolean array indicating masked spaxels.
    """

    # Identify spaxels close to the SNR threshold
    idx_snr = np.where(np.abs(snr - min_snr) < 2)[0]

    if len(idx_snr) > 0:
        meanmin_signal = np.mean(signal[idx_snr])
    else:
        meanmin_signal = np.min(signal)  # Fallback if no matching spaxels

    # Mask spaxels below the calculated signal threshold
    masked = signal < meanmin_signal

    if np.all(masked):
        print("No spaxels with S/N above the threshold. Ignoring potential warnings.")

    return masked


def sn_func(index, signal=None, noise=None, covar_vor=0.00):

    """
    Computes the signal-to-noise ratio in a bin for Voronoi binning.

    Parameters:
        index (np.ndarray): Indices of spaxels in the bin.
        signal (np.ndarray): Signal values for each spaxel.
        noise (np.ndarray): Noise values for each spaxel.
        covar_vor (float, optional): Correction factor for spatial correlations.

    Returns:
        sn (float): Estimated signal-to-noise ratio for the bin.
    """

    total_signal = np.sum(signal[index])
    total_noise = np.sqrt(np.sum(noise[index] ** 2))

    if total_noise == 0:
        return 0  # Prevent division by zero

    sn = total_signal / total_noise

    # Apply correction for spatial correlations
    if index.size > 1 and covar_vor > 0:
        sn /= 1 + covar_vor * np.log10(index.size)

    return sn


def generate_bins(config, cube, voronoi):

    """
    Applies Voronoi-binning or treats each spaxel as an individual bin if Voronoi binning is disabled.

    Parameters:
        config (dict): Configuration dictionary.
        cube (dict): Datacube containing spatial and SNR data.
        voronoi (bool): If True, applies Voronoi binning.

    Returns:
        None (Saves the binning results to a FITS file)
    """

    print("Defining the Voronoi bins")

    # Load the mask file safely
    mask_file = os.path.join(config['INFO']['OUTPUT'], f"{config['INFO']['RUN_NAME']}_mask.fits")

    if not os.path.isfile(mask_file):
        print(f"Mask file not found: {mask_file}")
        return "SKIP"

    # Open fits
    with fits.open(mask_file, mode="readonly") as hdul:
        mask_data = hdul[1].data  # reading

    mask = mask_data['MASK']
    idx_unmasked = np.where(mask == 0)[0]
    idx_masked = np.where(mask == 1)[0]

    # Partial function for SNR calculation
    sn_func_covariances = functools.partial(
        sn_func, covar_vor=config['BINNING'].get('COVARIANCE', 0.0))

    if voronoi:
        try:
            # Perform Voronoi binning
            bin_num, x_node, y_node, x_bar, y_bar, sn, n_pixels, _ = voronoi_2d_binning(
                cube['x'][idx_unmasked],
                cube['y'][idx_unmasked],
                cube['signal'][idx_unmasked],
                cube['noise'][idx_unmasked],
                config['BINNING']['TARGET_SNR'],
                plot=False,
                quiet=True,
                pixelsize=cube['pixelsize'],
                sn_func=sn_func_covariances
            )
            print(f"{np.max(bin_num) + 1} Voronoi bins generated!")

        except ValueError as e:
            # Handle case where no binning is needed
            if str(e) == 'All pixels have enough S/N and binning is not needed':
                print("Analysis will continue without Voronoi-binning!")
                bin_num = np.arange(len(idx_unmasked))
                x_node, y_node = cube['x'][idx_unmasked], cube['y'][idx_unmasked]
                sn = cube['snr'][idx_unmasked]
                n_pixels = np.ones(len(idx_unmasked))
            else:
                print(f"Voronoi-binning error: {e}")
                return "SKIP"
    if not voronoi:
        print(f"No Voronoi-binning! {len(idx_unmasked)} spaxels will be treated as individual bins.")
        bin_num = np.arange(len(idx_unmasked))
        x_node, y_node = cube['x'][idx_unmasked], cube['y'][idx_unmasked]
        sn = cube['snr'][idx_unmasked]
        n_pixels = np.ones(len(idx_unmasked))

    # Assign nearest Voronoi bin for masked pixels
    if len(idx_masked) > 0:
        pix_coords = np.column_stack((cube['x'][idx_masked], cube['y'][idx_masked]))
        bin_coords = np.column_stack((x_node, y_node))
        dists = dist.cdist(pix_coords, bin_coords, 'euclidean')
        bin_num_outside = np.argmin(dists, axis=1)
    else:
        bin_num_outside = np.array([])

    # Create extended bin list
    bin_num_long = np.full(len(cube['x']), np.nan)
    bin_num_long[idx_unmasked] = bin_num
    # bin_num_long[idx_masked] = -1 * bin_num_outside this was used in GIST. I changed as followings:
    bin_num_long[idx_masked] = -1  # Assign negative value to unselected spaxels

    # Save binning results
    # if not existing_bin:
    save_bin_info(
        config,
        cube['x'], cube['y'], cube['signal'], cube['snr'],
        bin_num_long, np.unique(bin_num), x_node, y_node, sn, n_pixels, cube['pixelsize'])
    # else:
    #     print('Using the user provided bin and mask info')


def save_bin_info(config, x, y, signal, snr, bin_num_new, ubins, x_node, y_node, sn, n_pixels, pixelsize):

    """
    Saves Voronoi binning results to a GIST-like FITS file.

    Parameters:
        config (dict): Configuration dictionary.
        x, y (np.ndarray): Spaxel coordinates.
        signal, snr (np.ndarray): Signal and SNR values.
        bin_num_new (np.ndarray): Assigned bin number for each spaxel.
        ubins (np.ndarray): Unique bin numbers.
        x_node, y_node (np.ndarray): Coordinates of bin centroids.
        sn (np.ndarray): SNR per bin.
        n_pixels (np.ndarray): Number of spaxels per bin.
        pixelsize (float): Pixel size for FITS metadata.

    Returns:
        None (Writes the FITS file)
    """

    output_file = os.path.join(config['INFO']['OUTPUT'], f"{config['INFO']['RUN_NAME']}_table.fits")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    print(f"Writing: {config['INFO']['RUN_NAME']}_table.fits")

    # Expand data to spaxel level
    x_node_new = np.zeros_like(x)
    y_node_new = np.zeros_like(y)
    sn_new = np.zeros_like(x, dtype=float)
    n_pixels_new = np.zeros_like(x, dtype=int)

    for i, ubin in enumerate(ubins):
        idx = np.where(ubin == np.abs(bin_num_new))[0]
        x_node_new[idx] = x_node[i]
        y_node_new[idx] = y_node[i]
        sn_new[idx] = sn[i]
        n_pixels_new[idx] = n_pixels[i]

    # Create FITS table
    columns = [
        fits.Column(name='ID', format='J', array=np.arange(len(x))),
        fits.Column(name='BIN_ID', format='J', array=bin_num_new),
        fits.Column(name='X', format='D', array=x),
        fits.Column(name='Y', format='D', array=y),
        fits.Column(name='FLUX', format='D', array=signal),
        fits.Column(name='SNR', format='D', array=snr),
        fits.Column(name='XBIN', format='D', array=x_node_new),
        fits.Column(name='YBIN', format='D', array=y_node_new),
        fits.Column(name='SNRBIN', format='D', array=sn_new),
        fits.Column(name='NSPAX', format='J', array=n_pixels_new),
    ]

    # writing fits
    with fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU.from_columns(columns, name="TABLE")]) as hdul:
        hdul.writeto(output_file, overwrite=True)

    fits.setval(output_file, "PIXSIZE", value=pixelsize)

    print(f"Wrote Voronoi table: {output_file}")


def prepare_mask_bin(config, cube, preview):

    """
    Reads spatial bins and mask file, applies binning to spectra, and saves or displays the binned spectra.

    Parameters:
        config (dict): Configuration dictionary.
        cube (dict): Datacube containing spectral data.
        preview (bool): If True, only displays the Voronoi map without saving.

    Returns:
        None
    """

    mask_file = os.path.join(config['INFO']['OUTPUT'], f"{config['INFO']['RUN_NAME']}_mask.fits")
    table_file = os.path.join(config['INFO']['OUTPUT'], f"{config['INFO']['RUN_NAME']}_table.fits")

    if not os.path.isfile(mask_file) or not os.path.isfile(table_file):
        print("Error: Mask or binning table file not found.")
        return

    # Load mask and binning table
    mask = fits.getdata(mask_file, ext=1)['MASK']
    unmasked_spaxels = np.where(mask == 0)[0]
    with fits.open(table_file, mode="readonly") as hdul:
        bin_num = hdul[1].data['BIN_ID'][unmasked_spaxels]

    # Apply Voronoi binning directly (integrated from apply_bin)
    print("Applying spatial bins to linear data...")
    bin_data, bin_error, bin_flux = perform_voronoi(bin_num, cube['spec'][:, unmasked_spaxels], cube['error'][:, unmasked_spaxels])
    print("Applied spatial bins.")

    if not preview:
        save_bin_spec(config, bin_data, bin_error, cube['wave'], "lin")
    else:
        # Display Voronoi map
        try:

            with fits.open(table_file, mode="readonly") as hdul:
                data = hdul[1].data

            x, y, bin_id, signal = data['X'], data['Y'], data['BIN_ID'], data['SNRBIN']

            # Set masked spaxels (bin_id < 0) to zero signal
            signal[bin_id < 0] = 0

            # Create grid
            x_bins, y_bins = np.unique(x), np.unique(y)
            grid_data = np.full((len(y_bins), len(x_bins)), np.nan)

            for i in range(len(x)):
                x_idx = np.searchsorted(x_bins, x[i])
                y_idx = np.searchsorted(y_bins, y[i])
                grid_data[y_idx, x_idx] = signal[i]

            # Plot Voronoi map
            plt.figure(figsize=(8, 6))
            plt.pcolormesh(x_bins, y_bins, grid_data, cmap='inferno', shading='auto')
            plt.colorbar(label="S/N")
            plt.xlabel("R [arcsec]")
            plt.ylabel("R [arcsec]")
            plt.title("Voronoi Map")
            plt.show()

        except Exception as e:
            print(f"Error: Unable to display Voronoi map: {e}")


def perform_voronoi(bin_num, spec, error):

    """
    Aggregates spaxels belonging to the same Voronoi bin.

    Parameters:
        bin_num (np.ndarray): Array of bin numbers for each spaxel.
        spec (np.ndarray): Spectral data array.
        error (np.ndarray): Error array for spectra.

    Returns:
        tuple: Binned spectra, errors, and flux.
    """

    ubins = np.unique(bin_num)
    nbins = len(ubins)
    npix = spec.shape[0]

    bin_data = np.zeros((npix, nbins))
    bin_error = np.zeros((npix, nbins))
    bin_flux = np.zeros(nbins)

    for i, ubin in enumerate(ubins):
        k = np.where(bin_num == ubin)[0]

        if k.size == 1:
            av_spec = spec[:, k].ravel()
            av_err_spec = np.sqrt(error[:, k]).ravel()
        else:
            av_spec = np.nansum(spec[:, k], axis=1)
            av_err_spec = np.sqrt(np.nansum(error[:, k] ** 2, axis=1))

        bin_data[:, i] = av_spec
        bin_error[:, i] = av_err_spec
        bin_flux[i] = np.mean(av_spec)

    return bin_data, bin_error, bin_flux


def save_bin_spec(config, log_spec, log_error, wavelength, flag):

    """
    Saves binned spectra and error spectra to a FITS file.

    Parameters:
        config (dict): Configuration dictionary.
        log_spec (np.ndarray): Array of binned spectra.
        log_error (np.ndarray): Array of error spectra.
        wavelength (np.ndarray): Wavelength array.
        flag (str): Data type flag (e.g., 'lin' for linear binning).

    Returns:
        None (Writes the FITS file)
    """

    output_file = os.path.join(config['INFO']['OUTPUT'], f"{config['INFO']['RUN_NAME']}_BinSpectra_linear.fits")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    print(f"Writing: {output_file}")

    # Opening fits and writing
    with fits.HDUList([
        fits.PrimaryHDU(),
        fits.BinTableHDU.from_columns([
            fits.Column(name='SPEC', format=f"{log_spec.shape[0]}D", array=log_spec.T),
            fits.Column(name='ESPEC', format=f"{log_spec.shape[0]}D", array=log_error.T)
        ], name='BIN_SPECTRA'),
        fits.BinTableHDU.from_columns([
            fits.Column(name='WAVE', format='D', array=wavelength)
        ], name='WAVE')
    ]) as hdul:
        hdul.writeto(output_file, overwrite=True)

    # adding info
    fits.setval(output_file, 'CRPIX1', value=1.0)
    fits.setval(output_file, 'CRVAL1', value=wavelength[0])
    fits.setval(output_file, 'CDELT1', value=wavelength[1] - wavelength[0])

    print(f"Wrote: {output_file}")


def extract(config, preview, voronoi, manual_bin, existing_bin):

    """
    Main function to run the extraction steps in sequence.

    Parameters:
        config (dict): Configuration dictionary.
        preview (bool): If True, runs in preview mode without saving.
        voronoi (bool): If True, applies Voronoi binning.
        manual_bin (bool): If True, namual bin has been selected.

    Returns:
        None
    """

    print("\n--- Starting Extraction Process ---\n")

    # 1) Read the datacube
    cube = reading_data(config)
    if cube == "SKIP":
        print("Extraction aborted: Failed to read the datacube.")
        return

    # 2) Apply spatial mask
    masking(config, cube, preview, manual_bin, existing_bin)

    # 3) Apply Voronoi binning
    binning(config, cube, preview, voronoi, manual_bin, existing_bin)

    # 4) Extract and save spectra
    save_spectra(config, cube, preview, existing_bin)

    print("\n--- Extraction Process Completed ---\n")


























def handle_existing_bin_files(input_folder, output_dir, ifs_run_id):
    """
    Copy and rename existing bin info files from input_folder to output_dir.

    Parameters
    ----------
    input_folder : str
        Directory where the *_table.fits and *_mask.fits files are located.
    output_dir : str
        Destination directory where the renamed files will be copied.
    ifs_run_id : str
        Identifier to replace the original RUN_NAME in the filenames.
    """

    os.makedirs(output_dir, exist_ok=True)

    # Search for *_table.fits and *_mask.fits files
    table_file = None
    mask_file = None
    try:
        for fname in os.listdir(input_folder):
            if fname.endswith('_table.fits'):
                table_file = fname
            elif fname.endswith('_mask.fits'):
                mask_file = fname
    except Exception as e:
        print('The folder does not exist!')

    if not table_file or not mask_file:
        print("Missing required files",
                       "Could not find both *_table.fits and *_mask.fits in the selected folder.")
        return

    # Define full paths
    table_path = os.path.join(input_folder, table_file)
    mask_path = os.path.join(input_folder, mask_file)

    # Define new names
    new_table_name = f"{ifs_run_id}_table.fits"
    new_mask_name = f"{ifs_run_id}_mask.fits"

    # Define new full paths
    new_table_path = os.path.join(output_dir, new_table_name)
    new_mask_path = os.path.join(output_dir, new_mask_name)

    # Copy and rename files
    shutil.copyfile(table_path, new_table_path)
    shutil.copyfile(mask_path, new_mask_path)

    print(f"Files copied and renamed to:\n{new_table_name}\n{new_mask_name}")
