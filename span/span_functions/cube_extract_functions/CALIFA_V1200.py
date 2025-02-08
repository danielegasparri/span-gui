from astropy.io import fits
import numpy as np
import os
import logging
import sys

# Add the parent directory to the Python path for importing modules
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

# ======================================
# Routine to load CALIFA cubes
# ======================================
def read_cube(config):
    """
    Reads a CALIFA V1200 data cube and processes its spectral and spatial information.

    Parameters:
        config (dict): Configuration dictionary containing the following keys:
            - GENERAL['INPUT']: Path to the FITS file.
            - GENERAL['REDSHIFT']: Redshift to correct the spectra.
            - READ_DATA['ORIGIN']: Origin for spatial coordinates.
            - READ_DATA['LMIN_TOT']: Minimum wavelength for spectra.
            - READ_DATA['LMAX_TOT']: Maximum wavelength for spectra.
            - READ_DATA['LMIN_SNR']: Minimum wavelength for SNR calculation.
            - READ_DATA['LMAX_SNR']: Maximum wavelength for SNR calculation.

    Returns:
        dict: A dictionary containing processed cube data including spatial coordinates,
              wavelengths, spectra, errors, SNR, signal, noise, and pixel size.
    """
    logging_blanks = (len(os.path.splitext(os.path.basename(__file__))[0]) + 33) * " "

    # Log the start of reading the cube
    logging.info(f"Reading the CALIFA V1200 cube: {config['GENERAL']['INPUT']}")
    print("Reading the CALIFA V1200 cube")

    # Open the FITS file
    hdu = fits.open(config['GENERAL']['INPUT'])
    hdr = hdu[0].header
    data = hdu[0].data
    s = data.shape
    spec = data.reshape(s[0], s[1] * s[2])

    # Read error spectra
    logging.info("Reading the error spectra from the cube")
    stat = hdu[1].data
    espec = stat.reshape(s[0], s[1] * s[2])

    # Calculate wavelength array
    wave = hdr['CRVAL3'] + np.arange(s[0]) * hdr['CDELT3']

    # Extract spatial coordinates
    origin = list(map(float, config['READ_DATA']['ORIGIN'].split(',')))
    xaxis = (np.arange(s[2]) - origin[0]) * hdr['CD2_2'] * 3600.0
    yaxis = (np.arange(s[1]) - origin[1]) * hdr['CD2_2'] * 3600.0
    x, y = np.meshgrid(xaxis, yaxis)
    x, y = x.ravel(), y.ravel()
    pixelsize = hdr['CD2_2'] * 3600.0

    logging.info(
        f"Extracting spatial information:\n"
        f"{logging_blanks}* Spatial coordinates centered at {origin}\n"
        f"{logging_blanks}* Spatial pixel size is {pixelsize:.2f} arcseconds"
    )

    # Shift wavelengths to rest-frame
    wave /= (1 + config['GENERAL']['REDSHIFT'])
    logging.info(f"Shifting spectra to rest-frame, assuming a redshift of {config['GENERAL']['REDSHIFT']}")

    # Filter spectra to specified wavelength range
    lmin, lmax = config['READ_DATA']['LMIN_TOT'], config['READ_DATA']['LMAX_TOT']
    idx = (wave >= lmin) & (wave <= lmax)
    wave, spec, espec = wave[idx], spec[idx, :], espec[idx, :]
    logging.info(f"Shortened spectra to the range {lmin} - {lmax} \u00c5.")

    # Convert error spectra to variances
    espec **= 2

    # Compute SNR per spaxel
    idx_snr = (wave >= config['READ_DATA']['LMIN_SNR']) & (wave <= config['READ_DATA']['LMAX_SNR'])
    signal = np.nanmedian(spec[idx_snr, :], axis=0)
    noise = np.abs(np.nanmedian(np.sqrt(espec[idx_snr, :]), axis=0))
    snr = signal / noise
    logging.info(
        f"Computed SNR in the range {config['READ_DATA']['LMIN_SNR']} - {config['READ_DATA']['LMAX_SNR']} \u00c5."
    )

    # Package data into a dictionary
    cube = {
        'x': x,
        'y': y,
        'wave': wave,
        'spec': spec,
        'error': espec,
        'snr': snr,
        'signal': signal,
        'noise': noise,
        'pixelsize': pixelsize,
    }

    # Log completion
    print(f"Finished reading the CALIFA V1200 cube: Read {len(cube['x'])} spectra!")
    logging.info(f"Finished reading the CALIFA V1200 cube: Read {len(cube['x'])} spectra!")

    return cube
