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

# Adding the parent directory to the Python path for module imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

# Importing necessary functions
from span_functions import utilities as uti

# ======================================
# Function to load MUSE cubes. Inspired by the GIST pipeline of Bittner et. al 2019
# ======================================
def read_cube(config):
    """
    Read and process a MUSE datacube, extracting spectral and spatial information.

    Parameters:
        config (dict): Configuration dictionary with input parameters.

    Returns:
        dict: A dictionary containing the processed cube data.
    """

    # Read the MUSE datacube
    print(f"Reading the MUSE-NFM cube: {config['INFO']['INPUT']}")
    hdu = fits.open(config['INFO']['INPUT'])
    hdr = hdu[1].header
    data = hdu[1].data

    # Reshape the spectral data
    s = np.shape(data)
    spec = np.reshape(data, [s[0], s[1] * s[2]])

    # Handle error spectra
    if len(hdu) == 3:
        print("Reading the error spectra from the cube.")
        stat = hdu[2].data
        espec = np.reshape(stat, [s[0], s[1] * s[2]])
    else:
        print("No error extension found. Estimating error spectra from the flux.")
        espec = np.array([uti.noise_spec(spec[:, i]) for i in range(spec.shape[1])]).T

    # Extract wavelength information
    wave = hdr['CRVAL3'] + np.arange(s[0]) * hdr['CD3_3']

    # Extract spatial coordinates
    origin = list(map(float, config['READ']['ORIGIN'].split(',')))
    xaxis = (np.arange(s[2]) - origin[0]) * hdr['CD2_2'] * 3600.0
    yaxis = (np.arange(s[1]) - origin[1]) * hdr['CD2_2'] * 3600.0
    x, y = np.meshgrid(xaxis, yaxis)
    x = x.ravel()
    y = y.ravel()
    pixelsize = hdr['CD2_2'] * 3600.0

    # De-redshift the spectra
    redshift = config['INFO']['REDSHIFT']
    wave /= (1 + redshift)
    print(f"Shifting spectra to rest-frame (redshift: {redshift}).")

    # Shorten spectra to the specified wavelength range
    lmin = config['READ']['LMIN_TOT']
    lmax = config['READ']['LMAX_TOT']
    idx = np.where((wave >= lmin) & (wave <= lmax))[0]
    wave = wave[idx]
    spec = spec[idx, :]
    espec = espec[idx, :]

    print(f"Shortening spectra to wavelength range: {lmin} - {lmax} Å.")

    # Compute SNR per spaxel
    idx_snr = np.where((wave >= config['READ']['LMIN_SNR']) & (wave <= config['READ']['LMAX_SNR']))[0]
    signal = np.nanmedian(spec[idx_snr, :], axis=0)
    noise = np.nanmedian(np.sqrt(espec[idx_snr, :]), axis=0) if len(hdu) == 3 else espec[0, :]
    snr = signal / noise

    print(f"Computed SNR in wavelength range: {config['READ']['LMIN_SNR']} - {config['READ']['LMAX_SNR']} Å.")

    # Store data in a structured dictionary
    cube = {
        'x': x, 'y': y, 'wave': wave, 'spec': spec, 'error': espec,
        'snr': snr, 'signal': signal, 'noise': noise, 'pixelsize': pixelsize
    }

    print(f"Finished reading MUSE cube: {len(cube['x'])} spectra loaded.")
    return cube
