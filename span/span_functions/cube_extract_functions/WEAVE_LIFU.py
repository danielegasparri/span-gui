########################################################################################
# MODIFIED VERSION OF THE SPECTRAL ROUTINE OF THE GIST PIPELINE OF BITTNER ET AL., 2019
######################## A SPECIAL THANKS TO ADRIAN BITTNER ############################
########################################################################################


from astropy.io import fits
import numpy as np
import os
import sys

# Add the parent directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from span_functions import utilities as uti

# ======================================
# Function to load WEAVE cubes. Inspired by the GIST pipeline of Bittner et. al 2019
# ======================================


def read_cube(config):
    """
    Reads a WEAVE data cube and extracts relevant spectral and spatial information.

    Parameters:
        config (dict): Configuration dictionary with input file paths and parameters.

    Returns:
        dict: Processed data cube with spectra, errors, SNR, spatial coordinates, and metadata.
    """

    print(f"Reading the WEAVE cube: {config['INFO']['INPUT']}")
    hdu = fits.open(config['INFO']['INPUT'])
    hdr = hdu[1].header
    data = hdu[1].data
    shape = data.shape  # (nwave, ny, nx)
    spec = data.reshape(shape[0], -1)

    if len(hdu) >=3:
        print("Reading error extension.")
        ivar = hdu[2].data.reshape(shape[0], -1)
        ivar[ivar <= 0] = 1e-5
        error = 1 / np.sqrt(ivar)

    else:
        print("No error extension found. Estimating error spectra from the flux.")
        error = np.array([uti.noise_spec(spec[:, i]) for i in range(spec.shape[1])]).T

    # Extract wavelength information
    wave = hdr['CRVAL3'] + np.arange(shape[0]) * hdr['CD3_3']

    # Spatial grid in arcsec
    origin = [float(val.strip()) for val in config['READ']['ORIGIN'].split(',')]
    try:
        scale = hdr['CD2_2'] * 3600.0  # degrees to arcsec
    except KeyError:
        print('Scale keyword not found. Showing scale in pixels')
        scale = 1.0  # fallback
    xaxis = (np.arange(shape[2]) - origin[0]) * scale
    yaxis = (np.arange(shape[1]) - origin[1]) * scale
    x, y = np.meshgrid(xaxis, yaxis)
    x, y = x.ravel(), y.ravel()

    # De-redshift
    redshift = config['INFO']['REDSHIFT']
    wave /= (1 + redshift)
    print(f"Shifting spectra to rest-frame (z = {redshift}).")

    # Limit to wavelength range
    lmin, lmax = config['READ']['LMIN_TOT'], config['READ']['LMAX_TOT']
    idx = (wave >= lmin) & (wave <= lmax)
    spec, error, wave = spec[idx, :], error[idx, :], wave[idx]
    print(f"Selected wavelength range: {lmin}-{lmax} Å")

    # Compute SNR
    idx_snr = (wave >= config['READ']['LMIN_SNR']) & (wave <= config['READ']['LMAX_SNR'])
    signal = np.nanmedian(spec[idx_snr, :], axis=0)
    noise = np.sqrt(np.nanmedian(error[idx_snr, :]**2, axis=0)) if len(hdu) > 2 else error[0, :]
    # noise = np.abs(np.nanmedian(np.sqrt(error[idx_snr, :]), axis=0)) if len(hdu) > 2 else error[0, :]

    snr = signal / noise
    # snr = np.nan_to_num(signal / noise, nan=0.0, posinf=0.0, neginf=0.0)

    pixelsize = scale

    cube = {
        'x': x, 'y': y, 'wave': wave,
        'spec': spec, 'error': error,
        'snr': snr, 'signal': signal, 'noise': noise,
        'pixelsize': pixelsize
    }

    print(f"Finished reading WEAVE cube. Total spaxels: {len(x)}.")
    return cube
