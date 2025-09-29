#SPectral ANalysis software (SPAN).
#Written by Daniele Gasparri#

"""
    Copyright (C) 2020-2025, Daniele Gasparri

    E-mail: daniele.gasparri@gmail.com

    SPAN is a GUI software that allows to modify and analyze 1D astronomical spectra.

    1. This software is licensed for non-commercial, academic and personal use only.
    2. The source code may be used and modified for research and educational purposes, 
    but any modifications must remain for private use unless explicitly authorized 
    in writing by the original author.
    3. Redistribution of the software in its original, unmodified form is permitted 
    for non-commercial purposes, provided that this license notice is always included.
    4. Redistribution or public release of modified versions of the source code 
    is prohibited without prior written permission from the author.
    5. Any user of this software must properly attribute the original author 
    in any academic work, research, or derivative project.
    6. Commercial use of this software is strictly prohibited without prior 
    written permission from the author.

    DISCLAIMER:
    THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES, OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT, OR OTHERWISE, ARISING FROM, OUT OF, OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""


from astropy.io import fits
import numpy as np
import os
import sys

# Aggiunge il percorso per importare moduli da directory superiori
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

# ======================================
# Routine to load CALIFA cubes. Inspired by the GIST pipeline of Bittner et. al 2019
# ======================================
def read_cube(config):
    """
    Reads a CALIFA V500 data cube and processes its spectral and spatial information.

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

    # Reading CALIFA datacubes
    print(f"Reading the CALIFA V500 cube: {config['INFO']['INPUT']}")

    # Opening the fits
    hdu = fits.open(config['INFO']['INPUT'])
    hdr = hdu[0].header
    data = hdu[0].data
    s = np.shape(data)
    spec = np.reshape(data, [s[0], s[1] * s[2]])

    # Reading error info
    print("Reading the error spectra from the cube")
    stat = hdu[1].data
    espec = np.reshape(stat, [s[0], s[1] * s[2]])

    # Wavelength
    wave = hdr['CRVAL3'] + (np.arange(s[0])) * hdr['CDELT3']

    # Spatial coordinates
    origin = [
        float(config['READ']['ORIGIN'].split(',')[0].strip()),
        float(config['READ']['ORIGIN'].split(',')[1].strip())
    ]
    xaxis = (np.arange(s[2]) - origin[0]) * hdr['CD2_2'] * 3600.0
    yaxis = (np.arange(s[1]) - origin[1]) * hdr['CD2_2'] * 3600.0
    x, y = np.meshgrid(xaxis, yaxis)
    x = np.reshape(x, [s[1] * s[2]])
    y = np.reshape(y, [s[1] * s[2]])
    pixelsize = hdr['CD2_2'] * 3600.0


    # De-redshift the spectra
    redshift = config['INFO']['REDSHIFT']
    wave /= (1 + redshift)
    print(f"Shifting spectra to rest-frame (redshift: {redshift}).")

    # Cropping
    lmin = config['READ']['LMIN_TOT']
    lmax = config['READ']['LMAX_TOT']
    idx = np.where(np.logical_and(wave >= lmin, wave <= lmax))[0]
    spec = spec[idx, :]
    espec = espec[idx, :]
    wave = wave[idx]
    print(f"Shortening spectra to the wavelength range from {lmin}A to {lmax}A.")

    # Calculating the variance
    espec = espec ** 2

    # S/N
    idx_snr = np.where(
        np.logical_and(
            wave >= config['READ']['LMIN_SNR'], wave <= config['READ']['LMAX_SNR']
        )
    )[0]
    signal = np.nanmedian(spec[idx_snr, :], axis=0)
    noise = np.abs(np.nanmedian(np.sqrt(espec[idx_snr, :]), axis=0))
    snr = signal / noise
    print(f"Computing the signal-to-noise ratio in the wavelength range from {config['READ']['LMIN_SNR']}A to {config['READ']['LMAX_SNR']}A.")

    # DIctionary with datacube info
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

    return cube
