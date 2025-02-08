from astropy.io import fits
import numpy as np
import os
import logging
import sys

# Add the parent directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from span_functions import der_snr as der_snr



# ======================================
# Function to load MUSE cubes
# ======================================
def read_cube(config):
    """
    Reads a MUSE data cube and extracts relevant spectral and spatial information.

    Parameters:
        config (dict): Configuration dictionary with input file paths and parameters.

    Returns:
        dict: Processed data cube containing spectra, errors, SNR, spatial coordinates, and metadata.
    """
    logging_blanks = (len(os.path.splitext(os.path.basename(__file__))[0]) + 33) * " "

    # Read the MUSE cube
    logging.info(f"Reading the MUSE-WFM cube: {config['GENERAL']['INPUT']}")
    hdu = fits.open(config['GENERAL']['INPUT'])
    hdr = hdu[1].header
    data = hdu[1].data
    shape = data.shape
    spec = data.reshape(shape[0], -1)

    # Handle error spectra or estimate them using DER_SNR
    if len(hdu) == 3:
        logging.info("Reading error spectra from the cube.")
        stat = hdu[2].data
        espec = stat.reshape(shape[0], -1)
    else:
        logging.info("No error extension found. Estimating error spectra with DER_SNR.")
        espec = np.array([der_snr.der_snr(spec[:, i]) for i in range(spec.shape[1])]).T

    # Extract wavelength information
    wave = hdr['CRVAL3'] + np.arange(shape[0]) * hdr['CD3_3']

    # Extract spatial coordinates
    origin = [float(val.strip()) for val in config['READ_DATA']['ORIGIN'].split(',')]
    xaxis = (np.arange(shape[2]) - origin[0]) * hdr['CD2_2'] * 3600.0
    yaxis = (np.arange(shape[1]) - origin[1]) * hdr['CD2_2'] * 3600.0
    x, y = np.meshgrid(xaxis, yaxis)
    x, y = x.ravel(), y.ravel()
    pixelsize = hdr['CD2_2'] * 3600.0

    logging.info(f"Spatial coordinates centered at {origin}, pixel size: {pixelsize:.3f}")

    # De-redshift the spectra
    redshift = config['GENERAL']['REDSHIFT']
    wave /= (1 + redshift)
    logging.info(f"Shifting spectra to rest-frame (redshift: {redshift}).")

    # Limit spectra to the specified wavelength range
    lmin, lmax = config['READ_DATA']['LMIN_TOT'], config['READ_DATA']['LMAX_TOT']
    idx = (wave >= lmin) & (wave <= lmax)
    spec, espec, wave = spec[idx, :], espec[idx, :], wave[idx]
    logging.info(f"Wavelength range limited to {lmin}-{lmax} \u00c5.")

    # Compute SNR per spaxel
    idx_snr = (wave >= config['READ_DATA']['LMIN_SNR']) & (wave <= config['READ_DATA']['LMAX_SNR'])
    signal = np.nanmedian(spec[idx_snr, :], axis=0)
    noise = np.abs(np.nanmedian(np.sqrt(espec[idx_snr, :]), axis=0)) if len(hdu) == 3 else espec[0, :]
    snr = signal / noise
    logging.info(f"Computed SNR in wavelength range {config['READ_DATA']['LMIN_SNR']}-{config['READ_DATA']['LMAX_SNR']} \u00c5.")

    # Store data in a structured dictionary
    cube = {
        'x': x, 'y': y, 'wave': wave, 'spec': spec, 'error': espec,
        'snr': snr, 'signal': signal, 'noise': noise, 'pixelsize': pixelsize
    }

    logging.info(f"Finished reading the MUSE cube. Total spectra: {len(cube['x'])}.")
    print(f"Read {len(cube['x'])} spectra from the MUSE-WFM cube.")

    return cube
