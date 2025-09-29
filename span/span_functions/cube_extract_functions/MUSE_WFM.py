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
# Function to load MUSE cubes. Inspired by the GIST pipeline of Bittner et. al 2019
# ======================================
def read_cube(config):
    """
    Reads a MUSE data cube and extracts relevant spectral and spatial information.

    Parameters:
        config (dict): Configuration dictionary with input file paths and parameters.

    Returns:
        dict: Processed data cube containing spectra, errors, SNR, spatial coordinates, and metadata.
    """

    # Read the MUSE cube
    print(f"Reading the MUSE-WFM cube: {config['INFO']['INPUT']}")
    hdu = fits.open(config['INFO']['INPUT'])
    hdr = hdu[1].header
    data = hdu[1].data
    shape = data.shape
    spec = data.reshape(shape[0], -1)

    if len(hdu) == 3:
        print("Reading error spectra from the cube.")
        stat = hdu[2].data
        espec = stat.reshape(shape[0], -1)
    else:
        print("No error extension found. Estimating error spectra from the flux.")
        espec = np.array([uti.noise_spec(spec[:, i]) for i in range(spec.shape[1])]).T

    # Extract wavelength information
    wave = hdr['CRVAL3'] + np.arange(shape[0]) * hdr['CD3_3']

    # Extract spatial coordinates
    origin = [float(val.strip()) for val in config['READ']['ORIGIN'].split(',')]
    xaxis = (np.arange(shape[2]) - origin[0]) * hdr['CD2_2'] * 3600.0
    yaxis = (np.arange(shape[1]) - origin[1]) * hdr['CD2_2'] * 3600.0
    x, y = np.meshgrid(xaxis, yaxis)
    x, y = x.ravel(), y.ravel()
    pixelsize = hdr['CD2_2'] * 3600.0

    print(f"Spatial coordinates centered at {origin}, pixel size: {pixelsize:.3f}")

    # De-redshift the spectra
    redshift = config['INFO']['REDSHIFT']
    wave /= (1 + redshift)
    print(f"Shifting spectra to rest-frame (redshift: {redshift}).")

    # Limit spectra to the specified wavelength range
    lmin, lmax = config['READ']['LMIN_TOT'], config['READ']['LMAX_TOT']
    idx = (wave >= lmin) & (wave <= lmax)
    spec, espec, wave = spec[idx, :], espec[idx, :], wave[idx]
    print(f"Wavelength range limited to {lmin}-{lmax} \u00c5.")

    # Compute SNR per spaxel
    idx_snr = (wave >= config['READ']['LMIN_SNR']) & (wave <= config['READ']['LMAX_SNR'])
    signal = np.nanmedian(spec[idx_snr, :], axis=0)
    noise = np.abs(np.nanmedian(np.sqrt(espec[idx_snr, :]), axis=0)) if len(hdu) == 3 else espec[0, :]
    snr = signal / noise
    print(f"Computed SNR in wavelength range {config['READ']['LMIN_SNR']}-{config['READ']['LMAX_SNR']} \u00c5.")

    # Store data in a structured dictionary
    cube = {
        'x': x, 'y': y, 'wave': wave, 'spec': spec, 'error': espec,
        'snr': snr, 'signal': signal, 'noise': noise, 'pixelsize': pixelsize
    }

    print(f"Finished reading the MUSE cube. Total spectra: {len(cube['x'])}.")

    return cube
