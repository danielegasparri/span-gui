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


#******************************************************************************************
#******************************************************************************************
#************************* FUNCTIONS FOR EXTRACTING DATA CUBES SPECTRA ********************
#******************************************************************************************
#******************************************************************************************

import os
import numpy as np
from astropy.io import fits
import sys
import importlib.util
import functools
from vorbin.voronoi_2d_binning import voronoi_2d_binning
import scipy.spatial.distance as dist
import matplotlib.pyplot as plt


#Functions to bin and extract 1D spectra from datacubes, using the GIST pipeline standard and syntax.
#The results are fully compatible with the GIST pipeline.

#function to create the dictionary (config) to be passed to the following functions
def buildConfigFromGUI(ifs_run_id, ifs_input, ifs_output, ifs_redshift, ifs_parallel,
                       ifs_ncpu, ifs_lfs_data, ifs_template, ifs_ow_config, ifs_ow_output,
                       ifs_routine_read, ifs_debug, ifs_origin, ifs_lmin_tot, ifs_lmax_tot,
                       ifs_lmin_snr, ifs_lmax_snr, ifs_mask_method, ifs_min_snr_mask,
                       ifs_mask, ifs_bin_method, ifs_target_snr, ifs_covariance,
                       ifs_prepare_method):
    """
    Returns a `configs` dictionary to be read from the following functions
    """
    configs = {
        "GENERAL": {
            "RUN_ID": ifs_run_id,
            "INPUT": ifs_input,
            "OUTPUT": ifs_output,
            "REDSHIFT": ifs_redshift,
            "PARALLEL": ifs_parallel,
            "NCPU": ifs_ncpu,
            "LSF_DATA": ifs_lfs_data,
            "LSF_TEMP": ifs_template,
            "OW_CONFIG": ifs_ow_config,
            "OW_OUTPUT": ifs_ow_output
        },
        "READ_DATA": {
            "METHOD": ifs_routine_read,
            "DEBUG": ifs_debug,
            "ORIGIN": ifs_origin,
            "LMIN_TOT": ifs_lmin_tot,
            "LMAX_TOT": ifs_lmax_tot,
            "LMIN_SNR": ifs_lmin_snr,
            "LMAX_SNR": ifs_lmax_snr
        },
        "SPATIAL_MASKING": {
            "METHOD": ifs_mask_method,
            "MIN_SNR": ifs_min_snr_mask,
            "MASK": ifs_mask
        },
        "SPATIAL_BINNING": {
            "METHOD": ifs_bin_method,
            "TARGET_SNR": ifs_target_snr,
            "COVARIANCE": ifs_covariance
        },
        "PREPARE_SPECTRA": {
            "METHOD": ifs_prepare_method
        }
    }

    return configs


#1) READ THE DATACUBE
def reading_data(config):
    """
    This function calls the readData routine specified by the user.
    """
    print("Step 1: Reading the tadacube")

    # Import the chosen readData routine
    try:
        routine_path = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "cube_extract_functions",
            config['READ_DATA']['METHOD'] + '.py'
        )
        spec = importlib.util.spec_from_file_location("", routine_path)
        print(f"Using the read-in routine for {config['READ_DATA']['METHOD']}")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    except Exception as e:
        print(f"Failed to import the read-in routine {config['READ_DATA']['METHOD']}.")

        return "SKIP"

    # Execute the chosen readData routine
    try:
        cube = module.read_cube(config)
    except Exception as e:
        print(f"Read-in routine {config['READ_DATA']['METHOD']} failed to read {config['GENERAL']['INPUT']}.")
        return "SKIP"

    return cube


#2) PERFORM THE SPATIAL MASKING
def masking(config, cube, preview):
    """
    This function calls the spatialMasking routine specified by the user.
    """

    print("Step 2: Applying masking, if any")
    # Check if outputs are already available
    output_file = os.path.join(
        config['GENERAL']['OUTPUT'],
        config['GENERAL']['RUN_ID']
    ) + "_mask.fits"

    if not preview and (os.path.isfile(output_file) and not config['GENERAL']['OW_OUTPUT']):
        print("Results of the spatialMasking module are already in the output directory. Module is skipped.")
        return

    # Execute the chosen spatialMasking routine
    generateSpatialMask(config, cube)



#3) PERFORM THE VORONOI REBINNING
def binning(config, cube, preview, voronoi):
    """
    This function calls the spatialBinning routine specified by the user.
    """
    print("Step 3: Applying binning")

    # Check if outputs are already available
    output_file = os.path.join(
        config['GENERAL']['OUTPUT'],
        config['GENERAL']['RUN_ID']
    ) + "_table.fits"

    if not preview and (os.path.isfile(output_file) and not config['GENERAL']['OW_OUTPUT']):
        print("Results of the module are already in the output directory. Module is skipped.")
        return

    # Execute the chosen spatialBinning routine
    try:
        generateSpatialBins(config, cube, voronoi)
    except Exception as e:
        print(f"Spatial binning routine {config['SPATIAL_BINNING']['METHOD']} failed.")

        return "SKIP"



#4) PREPARE AND SAVE THE SPECTRA
def save_spectra(config, cube, preview):
    """
    This function calls the prepareSpectra routine specified by the user.
    """
    print("Step 4: saving the extracted 1D spectra")

    # Check if outputs are already available
    output_prefix = os.path.join(
        config['GENERAL']['OUTPUT'],
        config['GENERAL']['RUN_ID']
    )
    if not preview and (not config['GENERAL']['OW_OUTPUT'] and
        os.path.isfile(output_prefix + "_BinSpectra_linear.fits")):
        print("Results of the module are already in the output directory. Module is skipped.")
        return

    # Execute the chosen routine
    try:
        prepSpectra(config, cube, preview)
    except Exception as e:
        print(f"Routine {config['PREPARE_SPECTRA']['METHOD']} failed.")

        return "SKIP"







###############################################################################
################ AUXILIARY FUNCTIONS TO PERFORM THE 4 STEPS ABOVE #############

def generateSpatialMask(config, cube):
    """
    Default implementation of the spatialMasking module.

    This function masks defunct spaxels, rejects spaxels with a SNR below
    a given threshold, and masks spaxels according to a provided mask file.
    Finally, all masks are combined and saved.
    """
    # Mask defunct spaxels
    masked_defunct = maskDefunctSpaxels(cube)

    # Mask spaxels with SNR below threshold
    masked_snr = applySNRThreshold(
        cube['snr'],
        cube['signal'],
        config['SPATIAL_MASKING']['MIN_SNR']
    )

    # Mask spaxels according to spatial mask file
    masked_mask = applyMaskFile(config, cube)

    # Create combined mask
    combined_mask_idx = np.where(
        np.logical_or.reduce((masked_defunct, masked_snr, masked_mask))
    )[0]
    combined_mask = np.zeros(len(cube['snr']), dtype=bool)
    combined_mask[combined_mask_idx] = True

    # Save mask to file
    saveMask(combined_mask, masked_defunct, masked_snr, masked_mask, config)


def maskDefunctSpaxels(cube):
    """
    Mask defunct spaxels, in particular those containing np.nan's
    or having a negative median.
    """
    idx_good = np.where(
        np.logical_and(
            np.all(~np.isnan(cube['spec']), axis=0),
            np.nanmedian(cube['spec'], axis=0) > 0.0
        )
    )[0]
    idx_bad = np.where(
        np.logical_or(
            np.any(np.isnan(cube['spec']), axis=0),
            np.nanmedian(cube['spec'], axis=0) <= 0.0
        )
    )[0]

    print(f"Masking defunct spaxels: {len(idx_bad)} spaxels are rejected.")

    masked = np.zeros(len(cube['snr']), dtype=bool)
    masked[idx_bad] = True
    return masked


def applySNRThreshold(snr, signal, min_snr):
    """
    Mask those spaxels that are above the isophote level with a mean
    signal-to-noise ratio of min_snr.
    """

    idx_snr = np.where(np.abs(snr - min_snr) < 2.)[0]
    meanmin_signal = np.mean(signal[idx_snr])
    idx_inside = np.where(signal >= meanmin_signal)[0]
    idx_outside = np.where(signal < meanmin_signal)[0]

    # Fallback in case there's no spaxel inside or outside
    if len(idx_inside) == 0 and len(idx_outside) == 0:
        idx_inside = np.arange(len(snr))
        idx_outside = np.array([], dtype=np.int64)
        print('No spaxels with the S/N trheshold to mask. Please, ignore the possible Python warning')

    masked = np.zeros(len(snr), dtype=bool)
    masked[idx_outside] = True
    return masked


def applyMaskFile(config, cube):
    """
    Select spaxels that are unmasked in the input masking file.
    """
    if (config['SPATIAL_MASKING']['MASK'] in [False, None]):
        print("No mask")
        idx_good = np.arange(len(cube['snr']))
        idx_bad = np.array([], dtype=np.int64)
    else:
        maskfile = os.path.join(
            os.path.dirname(config['GENERAL']['INPUT']),
            config['SPATIAL_MASKING']['MASK']
        )

        if os.path.isfile(maskfile):
            mask_data = fits.open(maskfile)[1].data
            s = np.shape(mask_data)
            mask_data = np.reshape(mask_data, s[0] * s[1])

            idx_good = np.where(mask_data == 0)[0]
            idx_bad = np.where(mask_data == 1)[0]

            print(f"Masking spaxels according to maskfile: {len(idx_bad)} spaxels are rejected.")

        else:
            print(f"No maskfile selected: {maskfile}")
            idx_good = np.arange(len(cube['snr']))
            idx_bad = np.array([], dtype=np.int64)

    masked = np.zeros(len(cube['snr']), dtype=bool)
    masked[idx_bad] = True
    return masked


def saveMask(combined_mask, masked_defunct, masked_snr, masked_mask, config):
    """
    Save the final combined mask to a FITS file.
    """
    outfits = os.path.join(config['GENERAL']['OUTPUT'], config['GENERAL']['RUN_ID'])
    outfits = os.path.abspath(outfits) + "_mask.fits"
    print("Writing: " + config['GENERAL']['RUN_ID'] + "_mask.fits")

    os.makedirs(os.path.dirname(outfits), exist_ok=True)

    # Primary HDU
    pri_hdu = fits.PrimaryHDU()

    # Table HDU with output data (0 → unmasked, 1 → masked)
    cols = [
        fits.Column(name='MASK', format='I', array=combined_mask.astype(int)),
        fits.Column(name='MASK_DEFUNCT', format='I', array=masked_defunct.astype(int)),
        fits.Column(name='MASK_SNR', format='I', array=masked_snr.astype(int)),
        fits.Column(name='MASK_FILE', format='I', array=masked_mask.astype(int))
    ]
    tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
    tbhdu.name = "MASKFILE"
    tbhdu.header['COMMENT'] = "Value 0  -->  unmasked"
    tbhdu.header['COMMENT'] = "Value 1  -->  masked"

    # Create HDU list and write to file
    hdulist = fits.HDUList([pri_hdu, tbhdu])
    hdulist.writeto(outfits, overwrite=True)

    print(f"Wrote mask file: {outfits}")



def sn_func(index, signal=None, noise=None, covar_vor=0.00):
    """
    This function is used by the Voronoi binning routine of Cappellari & Copin
    (2003) to estimate the noise in a bin from the spaxel-level noise.
    """
    # Sum the noise in the spaxels
    sn = np.sum(signal[index]) / np.sqrt(np.sum(noise[index] ** 2))

    # Account for spatial correlations in the noise
    sn /= 1 + covar_vor * np.log10(index.size)
    return sn


def generateSpatialBins(config, cube, voronoi):
    """
    This function applies the Voronoi-binning algorithm of Cappellari & Copin
    2003 to the data.
    """

    # Pass a function for the SNR calculation to the Voronoi-binning algorithm
    sn_func_covariances = functools.partial(
        sn_func, covar_vor=config['SPATIAL_BINNING']['COVARIANCE']
    )

    print("Defining the Voronoi bins")

    # Read maskfile
    maskfile = os.path.join(config['GENERAL']['OUTPUT'], config['GENERAL']['RUN_ID'])
    maskfile = os.path.abspath(maskfile) + "_mask.fits"
    mask = fits.open(maskfile)[1].data.MASK
    idx_unmasked = np.where(mask == 0)[0]
    idx_masked = np.where(mask == 1)[0]

    if voronoi:
        try:
            # Perform Voronoi binning
            bin_num, x_node, y_node, x_bar, y_bar, sn, n_pixels, _ = voronoi_2d_binning(
                cube['x'][idx_unmasked],
                cube['y'][idx_unmasked],
                cube['signal'][idx_unmasked],
                cube['noise'][idx_unmasked],
                config['SPATIAL_BINNING']['TARGET_SNR'],
                plot=False,
                quiet=True,
                pixelsize=cube['pixelsize'],
                sn_func=sn_func_covariances
            )

            print(f"{np.max(bin_num) + 1} Voronoi bins generated!")

        except ValueError as e:
            # If SNR is sufficient and no binning is needed
            if str(e) == 'All pixels have enough S/N and binning is not needed':
                print("Analysis will continue without Voronoi-binning!")

                bin_num, x_node, y_node, sn, n_pixels = noBinning(
                    cube['x'],
                    cube['y'],
                    cube['snr'],
                    idx_unmasked
                )
            else:
                print(f"The Voronoi-binning routine returned the following error:\n{e}")
                return "SKIP"


    else: #when NO voronoi binning is required, i.e. for manual binning
        try:

            print(f"NO Voronoi-binning! {len(idx_unmasked)} spaxels will be treated as Voronoi-bins.")
            bin_num, x_node, y_node, sn, n_pixels = noBinning(
                cube['x'],
                cube['y'],
                cube['snr'],
                idx_unmasked
            )
        except Exception as e:
            print("Error!")

    # Assign nearest Voronoi bin for masked pixels
    bin_num_outside = find_nearest_voronoibin(
        cube['x'], cube['y'], idx_masked, x_node, y_node
    )

    # Generate extended bin list
    ubins = np.unique(bin_num)
    bin_num_long = np.full(len(cube['x']), np.nan)
    bin_num_long[idx_unmasked] = bin_num
    # bin_num_long[idx_masked] = -1 * bin_num_outside
    bin_num_long[idx_masked] = -1 #ASSIGNING A NEGATIVE VALUE FOR ALL THE SPAXELS MASKED! WORKS BETTER LIKE THIS!

    # Save table for ALL spaxels (inside and outside the Voronoi region)
    save_table(
        config,
        cube['x'],
        cube['y'],
        cube['signal'],
        cube['snr'],
        bin_num_long,
        ubins,
        x_node,
        y_node,
        sn,
        n_pixels,
        cube['pixelsize']
    )


def noBinning(x, y, snr, idx_inside):
    """
    In case no Voronoi-binning is required, treat spaxels
    in the input data as Voronoi bins to continue the analysis.
    """
    bin_num = np.arange(len(idx_inside))
    x_node = x[idx_inside]
    y_node = y[idx_inside]
    sn = snr[idx_inside]
    n_pixels = np.ones(len(idx_inside))

    return bin_num, x_node, y_node, sn, n_pixels


def find_nearest_voronoibin(x, y, idx_outside, x_node, y_node):
    """
    This function determines the nearest Voronoi bin for spaxels which do
    not satisfy the minimum SNR threshold (masked spaxels).
    """
    x_out = x[idx_outside]
    y_out = y[idx_outside]
    pix_coords = np.column_stack((x_out, y_out))
    bin_coords = np.column_stack((x_node, y_node))

    dists = dist.cdist(pix_coords, bin_coords, 'euclidean')
    closest = np.argmin(dists, axis=1)
    return closest


def save_table(config, x, y, signal, snr, bin_num_new,
               ubins, x_node, y_node, sn, n_pixels, pixelsize):
    """
    Save all relevant information about the Voronoi binning to a FITS file.
    """
    outfits_table = os.path.join(
        config['GENERAL']['OUTPUT'],
        config['GENERAL']['RUN_ID']
    )
    outfits_table = os.path.abspath(outfits_table) + "_table.fits"

    print(f"Writing: {config['GENERAL']['RUN_ID']}_table.fits")

    # Expand data to spaxel level
    x_node_new = np.zeros(len(x))
    y_node_new = np.zeros(len(x))
    sn_new = np.zeros(len(x))
    n_pixels_new = np.zeros(len(x))

    for i, ubin in enumerate(ubins):
        idx = np.where(ubin == np.abs(bin_num_new))[0]
        x_node_new[idx] = x_node[i]
        y_node_new[idx] = y_node[i]
        sn_new[idx] = sn[i]
        n_pixels_new[idx] = n_pixels[i]

    # Primary HDU
    pri_hdu = fits.PrimaryHDU()

    # Table HDU with output data
    cols = [
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

    tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
    tbhdu.name = "TABLE"

    hdulist = fits.HDUList([pri_hdu, tbhdu])
    hdulist.writeto(outfits_table, overwrite=True)

    hdulist.close()

    fits.setval(outfits_table, "PIXSIZE", value=pixelsize)
    print(f"Wrote Voronoi table: {outfits_table}")



def prepSpectra(config, cube, preview):
    """
    This function performs the following tasks:
     * Read the spatial bins and mask file
     * Apply spatial bins to spectra
     * Save the binned spectra to disk
    """
    maskfile = os.path.join(config['GENERAL']['OUTPUT'], config['GENERAL']['RUN_ID'])
    maskfile = os.path.abspath(maskfile) + "_mask.fits"
    mask = fits.open(maskfile)[1].data.MASK
    idx_unmasked = np.where(mask == 0)[0]

    # Read binning pattern
    tablefile = os.path.join(config['GENERAL']['OUTPUT'], config['GENERAL']['RUN_ID'])
    tablefile = os.path.abspath(tablefile) + "_table.fits"
    bin_num = fits.open(tablefile)[1].data.BIN_ID[idx_unmasked]

    # Apply spatial bins to linear spectra
    bin_data, bin_error, bin_flux = applySpatialBins(
        bin_num,
        cube['spec'][:, idx_unmasked],
        cube['error'][:, idx_unmasked],
        "lin"
    )

    if not preview:
        # Save binned spectra
        saveBinSpectra(config, bin_data, bin_error, cube['wave'], "lin")

    else: # in the preview mode

        # Reading data from the table file generated
        tablefile = os.path.join(config['GENERAL']['OUTPUT'], config['GENERAL']['RUN_ID'])
        tablefile = os.path.abspath(tablefile) + "_table.fits"

        try:
            with fits.open(tablefile) as hdul:
                data = hdul[1].data
                x = data['X']
                y = data['Y']
                bin_id = data['BIN_ID']
                signal = data['SNRBIN']

            #considering the masked spaxels (if any) that have bin_id<0 and set the signal to zero
            signal[bin_id < 0] = 0
            # creating the regular grid
            x_bins = np.unique(x)
            y_bins = np.unique(y)

            # Create empty matrix for the map
            grid_data = np.full((len(y_bins), len(x_bins)), np.nan)

            # Filling the matrix
            for i in range(len(x)):
                x_idx = np.where(x_bins == x[i])[0][0]
                y_idx = np.where(y_bins == y[i])[0][0]
                grid_data[y_idx, x_idx] = signal[i]

            # Creating the coordinates
            X, Y = np.meshgrid(x_bins, y_bins)

            # Plottting the color map
            plt.figure(figsize=(8, 6))
            binplot = plt.pcolormesh(X, Y, grid_data, cmap='inferno', shading='auto')

            # Color bar on the side
            plt.colorbar(binplot, label="S/N")
            plt.xlabel("R [arcsec]")
            plt.ylabel("R [arcsec]")
            plt.title("Voronoi map")
            # plt.gca().invert_xaxis()
            # plt.gca().invert_yaxis()# Inverting the x axis.
            plt.show()

        except Exception as e:
            print("Error. Cannot display the voronoi map", e)



def saveBinSpectra(config, log_spec, log_error, wavelength, flag):
    """
    Save spatially binned spectra and error spectra to disk.
    """
    outfile = os.path.join(config['GENERAL']['OUTPUT'], config['GENERAL']['RUN_ID'])
    outfile = os.path.abspath(outfile)

    outfits_spectra = outfile + '_BinSpectra_linear.fits'
    print(f"Writing: {config['GENERAL']['RUN_ID']}_BinSpectra_linear.fits")

    npix = len(log_spec)

    # Create primary HDU
    pri_hdu = fits.PrimaryHDU()

    # Table HDU for spectra
    cols = [
        fits.Column(name='SPEC', format=str(npix) + 'D', array=log_spec.T),
        fits.Column(name='ESPEC', format=str(npix) + 'D', array=log_error.T)
    ]
    data_hdu = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
    data_hdu.name = 'BIN_SPECTRA'

    # Table HDU for wavelength
    cols = [fits.Column(name='WAVE', format='D', array=wavelength)]
    wavelength_hdu = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
    wavelength_hdu.name = 'WAVE'

    # Create HDU list and save to file
    hdulist = fits.HDUList([pri_hdu, data_hdu, wavelength_hdu])
    hdulist.writeto(outfits_spectra, overwrite=True)

    fits.setval(outfits_spectra, 'CRPIX1', value=1.0)
    fits.setval(outfits_spectra, 'CRVAL1', value=wavelength[0])
    fits.setval(outfits_spectra, 'CDELT1', value=wavelength[1] - wavelength[0])

    print(f"Wrote: {outfits_spectra}")


def applySpatialBins(bin_num, spec, espec, flag):
    """
    Apply the constructed spatial binning scheme to the spectra.
    """
    print(f"Applying the spatial bins to {flag}-data")
    bin_data, bin_error, bin_flux = spatialBinning(bin_num, spec, espec)
    print(f"Applied spatial bins to {flag}-data")
    return bin_data, bin_error, bin_flux


def spatialBinning(bin_num, spec, error):
    """
    Spaxels belonging to the same spatial bin are summed together.
    """
    ubins = np.unique(bin_num)
    nbins = len(ubins)
    npix = spec.shape[0]

    bin_data = np.zeros((npix, nbins))
    bin_error = np.zeros((npix, nbins))
    bin_flux = np.zeros(nbins)

    for i in range(nbins):
        k = np.where(bin_num == ubins[i])[0]
        if len(k) == 1:
            av_spec = spec[:, k]
            av_err_spec = np.sqrt(error[:, k])
        else:
            av_spec = np.nansum(spec[:, k], axis=1)
            av_err_spec = np.sqrt(np.sum(error[:, k], axis=1))

        bin_data[:, i] = av_spec.ravel()
        bin_error[:, i] = av_err_spec.ravel()
        bin_flux[i] = np.mean(av_spec, axis=0)

    return bin_data, bin_error, bin_flux


def extract(config, preview, voronoi):
    """
    Main function to run the extraction steps in sequence.
    """
    # 1) reading the cube
    cube = reading_data(config)
    if cube == "SKIP":
        return

    # 2) apply mask
    _ = masking(config, cube, preview)

    # 3) Voronoi binning
    _ = binning(config, cube, preview, voronoi)

    # 4) Save the spectra
    _ = save_spectra(config, cube, preview)

