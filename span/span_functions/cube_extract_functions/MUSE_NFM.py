from astropy.io import fits
import numpy as np
import os
import logging
import sys

# Adding the parent directory to the Python path for module imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

# Importing necessary functions
from span_functions import der_snr as der_snr

# ======================================
# Function to set DEBUG mode
# ======================================
# def set_debug(cube, xext, yext):
#     """
#     Activate DEBUG mode by restricting the cube to a single row of spaxels.
#
#     Parameters:
#         cube (dict): The cube data structure.
#         xext (int): Number of pixels along the x-axis.
#         yext (int): Number of pixels along the y-axis.
#
#     Returns:
#         dict: The modified cube containing only one row of spaxels.
#     """
#     logging.info("DEBUG mode is activated. Using only one row of spaxels.")
#     mid_row = int(yext / 2)
#     start_idx = mid_row * xext
#     end_idx = (mid_row + 1) * xext
#
#     for key in ['x', 'y', 'snr', 'signal', 'noise']:
#         cube[key] = cube[key][start_idx:end_idx]
#
#     for key in ['spec', 'error']:
#         cube[key] = cube[key][:, start_idx:end_idx]
#
#     return cube

# ======================================
# Function to load MUSE cubes
# ======================================
def read_cube(config):
    """
    Read and process a MUSE datacube, extracting spectral and spatial information.

    Parameters:
        config (dict): Configuration dictionary with input parameters.

    Returns:
        dict: A dictionary containing the processed cube data.
    """
    logging_blanks = (len(os.path.splitext(os.path.basename(__file__))[0]) + 33) * " "

    # Read the MUSE datacube
    logging.info(f"Reading the MUSE-NFM cube: {config['GENERAL']['INPUT']}")
    hdu = fits.open(config['GENERAL']['INPUT'])
    hdr = hdu[1].header
    data = hdu[1].data

    # Reshape the spectral data
    s = np.shape(data)
    spec = np.reshape(data, [s[0], s[1] * s[2]])

    # Handle error spectra
    if len(hdu) == 3:
        logging.info("Reading the error spectra from the cube.")
        stat = hdu[2].data
        espec = np.reshape(stat, [s[0], s[1] * s[2]])
    else:
        logging.info("No error extension found. Estimating error spectra with der_snr algorithm.")
        espec = np.array([der_snr(spec[:, i]) for i in range(spec.shape[1])]).T

    # Extract wavelength information
    wave = hdr['CRVAL3'] + np.arange(s[0]) * hdr['CD3_3']

    # Extract spatial coordinates
    origin = list(map(float, config['READ_DATA']['ORIGIN'].split(',')))
    xaxis = (np.arange(s[2]) - origin[0]) * hdr['CD2_2'] * 3600.0
    yaxis = (np.arange(s[1]) - origin[1]) * hdr['CD2_2'] * 3600.0
    x, y = np.meshgrid(xaxis, yaxis)
    x = x.ravel()
    y = y.ravel()
    pixelsize = hdr['CD2_2'] * 3600.0

    logging.info(
        f"Extracting spatial information:\n"
        f"{logging_blanks}* Spatial coordinates centered at {origin}\n"
        f"{logging_blanks}* Pixel size: {pixelsize} arcsec"
    )

    # Apply redshift correction
    wave /= (1 + config['GENERAL']['REDSHIFT'])
    logging.info(f"Shifting spectra to rest-frame (redshift: {config['GENERAL']['REDSHIFT']}).")

    # Shorten spectra to the specified wavelength range
    lmin = config['READ_DATA']['LMIN_TOT']
    lmax = config['READ_DATA']['LMAX_TOT']
    idx = np.where((wave >= lmin) & (wave <= lmax))[0]
    wave = wave[idx]
    spec = spec[idx, :]
    espec = espec[idx, :]

    logging.info(
        f"Shortening spectra to wavelength range: {lmin} - {lmax} Å."
    )

    # Compute SNR per spaxel
    idx_snr = np.where((wave >= config['READ_DATA']['LMIN_SNR']) & (wave <= config['READ_DATA']['LMAX_SNR']))[0]
    signal = np.nanmedian(spec[idx_snr, :], axis=0)
    noise = np.nanmedian(np.sqrt(espec[idx_snr, :]), axis=0) if len(hdu) == 3 else espec[0, :]
    snr = signal / noise

    logging.info(
        f"Computed SNR in wavelength range: {config['READ_DATA']['LMIN_SNR']} - {config['READ_DATA']['LMAX_SNR']} Å."
    )

    # Store data in a structured dictionary
    cube = {
        'x': x, 'y': y, 'wave': wave, 'spec': spec, 'error': espec,
        'snr': snr, 'signal': signal, 'noise': noise, 'pixelsize': pixelsize
    }

    # # Apply DEBUG mode if specified
    # if config['READ_DATA']['DEBUG']:
    #     cube = set_debug(cube, s[2], s[1])

    logging.info(f"Finished reading MUSE cube: {len(cube['x'])} spectra loaded.")
    return cube
