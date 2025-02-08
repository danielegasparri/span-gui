from astropy.io import fits
import numpy as np
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from span_functions import der_snr as der_snr

# ======================================
# Routine to load MUSE-cubes. Inspired by the GIST pipeline of Bittner et. al 2019
# ======================================
def read_cube(config):


    # Read MUSE-cube
    print(f"Reading the MUSE-WFM cube: {config['GENERAL']['INPUT']}")

    # Reading the cube
    hdu   = fits.open(config['GENERAL']['INPUT'])
    hdr   = hdu[1].header
    data  = hdu[1].data
    s     = np.shape(data)
    spec  = np.reshape(data,[s[0],s[1]*s[2]])

    # Read the error spectra if available. Otherwise estimate the errors with the der_snr algorithm
    if len(hdu) == 3:
        print("Reading the error spectra from the cube")
        stat  = hdu[2].data
        espec = np.reshape(stat,[s[0],s[1]*s[2]])
    elif len(hdu) == 2:
        print("No error extension found. Estimating the error spectra with the der_snr algorithm")
        espec = np.zeros( spec.shape )
        for i in range( 0, spec.shape[1] ):
            espec[:,i] = der_snr.der_snr( spec[:,i] )

    # Getting the wavelength info
    wave = hdr['CRVAL3']+(np.arange(s[0]))*hdr['CD3_3']
    
    # Getting the spatial coordinates
    origin = [ float(config['READ_DATA']['ORIGIN'].split(',')[0].strip()), float(config['READ_DATA']['ORIGIN'].split(',')[1].strip()) ]
    xaxis = (np.arange(s[2]) - origin[0]) * hdr['CD2_2']*3600.0
    yaxis = (np.arange(s[1]) - origin[1]) * hdr['CD2_2']*3600.0
    x, y  = np.meshgrid(xaxis,yaxis)
    x     = np.reshape(x,[s[1]*s[2]])
    y     = np.reshape(y,[s[1]*s[2]])
    pixelsize = hdr['CD2_2']*3600.0

    # De-redshift the spectra
    redshift = config['GENERAL']['REDSHIFT']
    wave /= (1 + redshift)
    print(f"Shifting spectra to rest-frame (redshift: {redshift}).")

    # Shorten spectra to required wavelength range
    lmin  = config['READ_DATA']['LMIN_TOT']
    lmax  = config['READ_DATA']['LMAX_TOT']
    idx   = np.where( np.logical_and( wave >= lmin, wave <= lmax ) )[0]
    spec  = spec[idx,:]
    espec = espec[idx,:]
    wave  = wave[idx]
    print(f"Shortening spectra to wavelength range: {lmin} - {lmax} Å.")

    # Computing the SNR per spaxel
    idx_snr   = np.where( np.logical_and.reduce([ \
        wave >= config['READ_DATA']['LMIN_SNR'], \
        wave <= config['READ_DATA']['LMAX_SNR'], \
        np.logical_or( wave < 5780/(1+config['GENERAL']['REDSHIFT']), wave > 6050/(1+config['GENERAL']['REDSHIFT'])) ]))[0]
    signal  = np.nanmedian(spec[idx_snr,:],axis=0)
    if len(hdu) == 3:
        noise  = np.abs(np.nanmedian(np.sqrt(espec[idx_snr,:]),axis=0))
    elif len(hdu) == 2:
        noise = espec[0,:]    # DER_SNR returns constant error spectra
    snr    = signal / noise
    print(f"Computed SNR in wavelength range: {config['READ_DATA']['LMIN_SNR']} - {config['READ_DATA']['LMAX_SNR']} Å.")

    # Replacing the np.nan in the laser region by the median of the spectrum
    idx_laser          = np.where( np.logical_and( wave > 5780 / (1+config['GENERAL']['REDSHIFT']), wave < 6050 / (1+config['GENERAL']['REDSHIFT'])) )[0]
    spec[idx_laser,:]  = signal
    espec[idx_laser,:] = noise
    print("Replacing the spectral region affected by the LGS (5780A-6050A) with the median signal of the spectra.")

    # Storing everything into a structure
    cube = {'x':x, 'y':y, 'wave':wave, 'spec':spec, 'error':espec, 'snr':snr, 'signal':signal, 'noise':noise, 'pixelsize':pixelsize}

    print(f"Finished reading MUSE cube: {len(cube['x'])} spectra loaded.")

    return(cube)
