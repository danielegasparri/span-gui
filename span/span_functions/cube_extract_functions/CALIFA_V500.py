from astropy.io import fits
import numpy as np
import os
import logging
import sys

# Aggiunge il percorso per importare moduli da directory superiori
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))


def read_cube(config):
    """
    Legge un datacube CALIFA V500, estrae le informazioni necessarie e calcola i rapporti segnale-rumore (SNR).

    Args:
        config (dict): Configurazione contenente percorsi dei file e parametri per la lettura del datacube.

    Returns:
        dict: Dizionario contenente coordinate spaziali, spettro, errore, SNR e altre informazioni.
    """
    logging_blanks = (len(os.path.splitext(os.path.basename(__file__))[0]) + 33) * " "

    # Legge il datacube CALIFA
    print("Reading the CALIFA V500 cube")
    logging.info("Reading the CALIFA V500 cube: " + config['GENERAL']['INPUT'])

    # Apertura del file FITS
    hdu = fits.open(config['GENERAL']['INPUT'])
    hdr = hdu[0].header
    data = hdu[0].data
    s = np.shape(data)
    spec = np.reshape(data, [s[0], s[1] * s[2]])

    # Lettura degli spettri di errore
    logging.info("Reading the error spectra from the cube")
    stat = hdu[1].data
    espec = np.reshape(stat, [s[0], s[1] * s[2]])

    # Calcolo della lunghezza d'onda
    wave = hdr['CRVAL3'] + (np.arange(s[0])) * hdr['CDELT3']

    # Calcolo delle coordinate spaziali
    origin = [
        float(config['READ_DATA']['ORIGIN'].split(',')[0].strip()),
        float(config['READ_DATA']['ORIGIN'].split(',')[1].strip())
    ]
    xaxis = (np.arange(s[2]) - origin[0]) * hdr['CD2_2'] * 3600.0
    yaxis = (np.arange(s[1]) - origin[1]) * hdr['CD2_2'] * 3600.0
    x, y = np.meshgrid(xaxis, yaxis)
    x = np.reshape(x, [s[1] * s[2]])
    y = np.reshape(y, [s[1] * s[2]])
    pixelsize = hdr['CD2_2'] * 3600.0

    logging.info(
        "Extracting spatial information:\n"
        + logging_blanks + f"* Spatial coordinates are centred to {origin}\n"
        + logging_blanks + f"* Spatial pixelsize is {pixelsize}"
    )

    # Shift degli spettri al rest-frame
    wave = wave / (1 + config['GENERAL']['REDSHIFT'])
    logging.info(f"Shifting spectra to rest-frame, assuming a redshift of {config['GENERAL']['REDSHIFT']}")

    # Riduzione degli spettri al range di lunghezze d'onda richiesto
    lmin = config['READ_DATA']['LMIN_TOT']
    lmax = config['READ_DATA']['LMAX_TOT']
    idx = np.where(np.logical_and(wave >= lmin, wave <= lmax))[0]
    spec = spec[idx, :]
    espec = espec[idx, :]
    wave = wave[idx]
    logging.info(
        f"Shortening spectra to the wavelength range from {lmin}A to {lmax}A."
    )

    # Conversione degli errori in varianze
    espec = espec ** 2

    # Calcolo del rapporto segnale-rumore (SNR) per spaxel
    idx_snr = np.where(
        np.logical_and(
            wave >= config['READ_DATA']['LMIN_SNR'], wave <= config['READ_DATA']['LMAX_SNR']
        )
    )[0]
    signal = np.nanmedian(spec[idx_snr, :], axis=0)
    noise = np.abs(np.nanmedian(np.sqrt(espec[idx_snr, :]), axis=0))
    snr = signal / noise
    logging.info(
        f"Computing the signal-to-noise ratio in the wavelength range from {config['READ_DATA']['LMIN_SNR']}A to {config['READ_DATA']['LMAX_SNR']}A."
    )

    # Creazione del dizionario con i dati del datacube
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

    print(f"Finished reading the CALIFA V500 cube: Read {len(cube['x'])} spectra!")
    logging.info(f"Finished reading the CALIFA V500 cube: Read {len(cube['x'])} spectra!")

    return cube
