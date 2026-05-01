from astropy.io import fits
import numpy as np
import os
import sys

# Add the parent directory to the Python path for importing modules
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

# ======================================
# Routine to load MEGARA cubes
# ======================================
def read_cube(config):
    """
    Reads a MEGARA data cube and processes its spectral and spatial information.

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

    NOTE:
        The MEGARA cubes used here do not have an error/variance extension, so the
        noise is estimated directly from the spectra.
        If the flux unit is Jy, this routine converts it to erg s^-1 cm^-2 A^-1.
    """

    print(f"Reading the MEGARA cube: {config['INFO']['INPUT']}")

    # Open the FITS file
    hdu = fits.open(config['INFO']['INPUT'])
    hdr = hdu[0].header
    data = hdu[0].data
    s = data.shape
    spec = data.reshape(s[0], s[1] * s[2])

    # Calculate wavelength array
    wave = hdr['CRVAL3'] + np.arange(s[0]) * hdr['CDELT3']

    # Extract spatial coordinates
    origin = list(map(float, config['READ']['ORIGIN'].split(',')))
    xaxis = (origin[0] - np.arange(s[2])) * hdr['CDELT1'] * 3600.0
    yaxis = (np.arange(s[1]) - origin[1]) * hdr['CDELT2'] * 3600.0
    
    x, y = np.meshgrid(xaxis, yaxis)
    x, y = x.ravel(), y.ravel()
    pixelsize = np.abs(hdr['CDELT2']) * 3600.0

    # Convert from Jy to Flambda, if needed
    if 'BUNIT' in hdr and hdr['BUNIT'].strip().lower() == 'jy':
        print("Converting flux from Jy to erg/s/cm^2/A")
        conversion = (2.99792458e-5 / (wave**2))[:, np.newaxis]
        spec = spec * conversion

    # De-redshift the spectra
    redshift = config['INFO']['REDSHIFT']
    wave /= (1 + redshift)
    print(f"Shifting spectra to rest-frame (redshift: {redshift}).")

    # Filter spectra to specified wavelength range
    lmin, lmax = config['READ']['LMIN_TOT'], config['READ']['LMAX_TOT']
    idx = (wave >= lmin) & (wave <= lmax)
    wave, spec = wave[idx], spec[idx, :]
    print(f"Shortening spectra to wavelength range: {lmin} - {lmax} Å.")

    # Compute SNR per spaxel
    idx_snr = (wave >= config['READ']['LMIN_SNR']) & (wave <= config['READ']['LMAX_SNR'])

    signal = np.nanmedian(spec[idx_snr, :], axis=0)

    # Robust noise estimate from the spectra
    med_spec = np.nanmedian(spec[idx_snr, :], axis=0)
    noise = 1.4826 * np.nanmedian(np.abs(spec[idx_snr, :] - med_spec), axis=0)

    # Protect against zero / NaN noise
    good_noise = np.isfinite(noise) & (noise > 0)
    if np.any(good_noise):
        fallback_noise = np.nanmedian(noise[good_noise])
    else:
        fallback_noise = 1.0

    noise[~good_noise] = fallback_noise

    snr = signal / noise
    print(f"Computed SNR in wavelength range: {config['READ']['LMIN_SNR']} - {config['READ']['LMAX_SNR']} Å.")

    # Build variance spectra with same shape as spec
    espec = np.tile(noise**2, (len(wave), 1))

    print("Warning: no error/variance extension found in this MEGARA cube. Using a robust noise estimate from the spectra.")

    # Output in SPAN format
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

    print(f"Finished reading the MEGARA cube: Read {len(cube['x'])} spectra!")

    return cube
