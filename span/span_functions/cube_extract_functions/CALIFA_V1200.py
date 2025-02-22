########################################################################################
# MODIFIED VERSION OF THE SPECTRAL ROUTINE OF THE GIST PIPELINE OF BITTNER ET AL., 2019
######################## A SPECIAL THANKS TO ADRIAN BITTNER ############################
########################################################################################


from astropy.io import fits
import numpy as np
import os
import sys

# Add the parent directory to the Python path for importing modules
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

# ======================================
# Routine to load CALIFA cubes. Inspired by the GIST pipeline of Bittner et. al 2019
# ======================================
def read_cube(config):
    """
    Reads a CALIFA V1200 data cube and processes its spectral and spatial information.

    Parameters:
        config (dict): Configuration dictionary containing the following keys:
            - INFO['INPUT']: Path to the FITS file.
            - INFO['REDSHIFT']: Redshift to correct the spectra.
            - READ['ORIGIN']: Origin for spatial coordinates.
            - READ['LMIN_TOT']: Minimum wavelength for spectra.
            - READ['LMAX_TOT']: Maximum wavelength for spectra.
            - READ['LMIN_SNR']: Minimum wavelength for SNR calculation.
            - READ['LMAX_SNR']: Maximum wavelength for SNR calculation.

    Returns:
        dict: A dictionary containing processed cube data including spatial coordinates,
              wavelengths, spectra, errors, SNR, signal, noise, and pixel size.
    """

    # Log the start of reading the cube
    print(f"Reading the CALIFA V1200 cube: {config['INFO']['INPUT']}")

    # Open the FITS file
    hdu = fits.open(config['INFO']['INPUT'])
    hdr = hdu[0].header
    data = hdu[0].data
    s = data.shape
    spec = data.reshape(s[0], s[1] * s[2])

    # Read error spectra
    print("Reading the error spectra from the cube")
    stat = hdu[1].data
    espec = stat.reshape(s[0], s[1] * s[2])

    # Calculate wavelength array
    wave = hdr['CRVAL3'] + np.arange(s[0]) * hdr['CDELT3']

    # Extract spatial coordinates
    origin = list(map(float, config['READ']['ORIGIN'].split(',')))
    xaxis = (np.arange(s[2]) - origin[0]) * hdr['CD2_2'] * 3600.0
    yaxis = (np.arange(s[1]) - origin[1]) * hdr['CD2_2'] * 3600.0
    x, y = np.meshgrid(xaxis, yaxis)
    x, y = x.ravel(), y.ravel()
    pixelsize = hdr['CD2_2'] * 3600.0


    # De-redshift the spectra
    redshift = config['INFO']['REDSHIFT']
    wave /= (1 + redshift)
    print(f"Shifting spectra to rest-frame (redshift: {redshift}).")

    # Filter spectra to specified wavelength range
    lmin, lmax = config['READ']['LMIN_TOT'], config['READ']['LMAX_TOT']
    idx = (wave >= lmin) & (wave <= lmax)
    wave, spec, espec = wave[idx], spec[idx, :], espec[idx, :]
    print(f"Shortening spectra to wavelength range: {lmin} - {lmax} Å.")

    # Convert error spectra to variances
    espec **= 2

    # Compute SNR per spaxel
    idx_snr = (wave >= config['READ']['LMIN_SNR']) & (wave <= config['READ']['LMAX_SNR'])
    signal = np.nanmedian(spec[idx_snr, :], axis=0)
    noise = np.abs(np.nanmedian(np.sqrt(espec[idx_snr, :]), axis=0))
    snr = signal / noise
    print(f"Computed SNR in wavelength range: {config['READ']['LMIN_SNR']} - {config['READ']['LMAX_SNR']} Å.")

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

    # end
    print(f"Finished reading the CALIFA V1200 cube: Read {len(cube['x'])} spectra!")

    return cube
