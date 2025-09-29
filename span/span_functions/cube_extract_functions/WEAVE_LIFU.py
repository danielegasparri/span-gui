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
