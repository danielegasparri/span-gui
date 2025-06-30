#SPectral ANalysis software (SPAN)
#Written by Daniele Gasparri#
#SPectral ANalysis software (SPAN).
#Written by Daniele Gasparri#

"""
    Copyright (C) 2020-2025, Daniele Gasparri

    E-mail: daniele.gasparri@gmail.com

    SPAN is a GUI interface that allows to modify and analyse 1D astronomical spectra.

    1. This software is licensed **for non-commercial use only**.
    2. The source code may be **freely redistributed**, but this license notice must always be included.
    3. Any user who redistributes or uses this software **must properly attribute the original author**.
    4. The source code **may be modified** for non-commercial purposes, but any modifications must be clearly documented.
    5. **Commercial use is strictly prohibited** without prior written permission from the author.

    DISCLAIMER:
    THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES, OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT, OR OTHERWISE, ARISING FROM, OUT OF, OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""

#******************************************************************************************
#******************************************************************************************
#*************************** UTILITIES FUNCTIONS FOR SPAN *********************************
#******************************************************************************************
#******************************************************************************************

try:#Local imports
    from span_functions import spec_manipul as spman

except ModuleNotFoundError: #local import if executed as package
    from . import spec_manipul as spman


#Python imports
import numpy as np
import math as mt


from astropy.io import fits
from astropy.table import Table

from scipy.optimize import curve_fit
import os


#1) Show sampling and identify linear or log spectrum (hopefully!)
def show_sampling(wavelength):

    """
    This function shows the sampling of the selected 1D spectrum.
    Input: wavelength array of the spectrum
    Output: float step value (in A), bool linear (True) or log (False) sampling
    """

    step1 = wavelength[1] - wavelength[0]
    linear_step = True
    return step1, linear_step


#*************************************************************************************************
#2) Show the SNR of a selected window
def show_snr(wavelength, flux, wave_snr, epsilon_wave_snr):

    """
    This function shows the SNR of the selected 1D spectrum.
    Input: wavelength and flux arrays, central wavelength of the
    window to measure the SNR and delta wavelength.
    Output: float SNR per pix and per Angstrom
    """

    step = wavelength[1] - wavelength[0]
    step2 = wavelength[-1] - wavelength[-2]
    epsilon = 1e-4
    if abs(step - step2) > epsilon:
        wavelength, flux, npoint_resampled = spman.resample(wavelength, flux, step)
        print('Spectrum resampled to a linear step')

    mask = (wavelength >= wave_snr - epsilon_wave_snr) & (wavelength <= wave_snr + epsilon_wave_snr)
    flux_snr = flux[mask]
    mean_flux = np.mean(flux_snr)
    snr_pix = mean_flux / np.std(flux_snr)
    snr_ang = snr_pix * np.sqrt(1 / (step))  # Supposing all the units in A!

    return snr_pix, snr_ang


#*************************************************************************************************
#3) Show the header of a fits file
def show_hdr(spec_name):

    """
    This function shows the fits header contained in the primary
    extension of a fits spectrum.
    Input: string name of the spectrum
    Output: string header
    """

    if spec_name.endswith(('.txt', '.dat')):
        return 'ASCII files do not have header!'

    with fits.open(spec_name) as hdu:
        hdr = hdu[0].header
        return repr(hdr)


#*************************************************************************************************
#4) Convert a spectrum to ASCII or binary fits file
def convert_spec(wavelength, flux, spec_name, type_spec_to_convert, lambda_units):

    """
    This function converts the selected spectrum to ASCII or fits file
    Input: wavelength array, flux array, name of the spectrum, type to convert ('ASCII' or 'fits')
    Output: fits or ASCII (.dat) file of the spectrum containing the wavelength and the flux
    """

    if lambda_units == 'nm':
        wavelength = wavelength*10
    if lambda_units == 'mu':
        wavelength = wavelength/10000

    if type_spec_to_convert == 'ASCII':
        new_spec_name = os.path.splitext(spec_name)[0]
        filename = f'{new_spec_name}_SPAN.txt'
        np.savetxt(filename, np.column_stack([wavelength, flux]), header="wavelength\tflux", delimiter='\t')
    elif type_spec_to_convert == 'FITS':
        new_spec_name = os.path.splitext(spec_name)[0]
        filename = f'{new_spec_name}_SPAN.fits'

        CRVAL1 = wavelength[0]
        CDELT1 = wavelength[1]-wavelength[0]
        NAXIS1 = len(wavelength)

        header = fits.Header()
        header['SIMPLE'] = (True, 'conforms to FITS standard')
        header['BITPIX'] = 8
        header['NAXIS'] = 1
        header['NAXIS1'] = NAXIS1
        header['CRVAL1'] = CRVAL1
        header['CDELT1'] = CDELT1
        header['EXTEND'] = True

        hdu = fits.PrimaryHDU(data=flux, header=header)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(filename, overwrite=True)

        # closing fits file
        hdulist.close()


#*************************************************************************************************
#5) Save fits
def save_fits(wavelength, flux, file_name):

    """
    This function is used by SPAN to save the processed spectra
    of the spectra manipulation panel to fits files.
    Input: wavelength array, flux array, name of the spectrum.
    Output: 1D fits file IRAF style.

    """

    CRVAL1 = wavelength[0]
    CDELT1 = wavelength[1]-wavelength[0]
    NAXIS1 = len(wavelength)

    header = fits.Header()
    header['SIMPLE'] = (True, 'conforms to FITS standard')
    header['BITPIX'] = 8
    header['NAXIS'] = 1
    header['NAXIS1'] = NAXIS1
    header['CRVAL1'] = CRVAL1
    header['CDELT1'] = CDELT1
    header['EXTEND'] = True
    header['CREATOR'] = ('SPAN', 'Generated by SPAN')  # Custom keyword for creator
    header.add_comment('Generated by SPAN')

    hdu = fits.PrimaryHDU(data=flux, header=header)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(file_name, overwrite=True)

    # closing fits file
    hdulist.close()


#*************************************************************************************************
#6) Save 2d fits table. Useful if the wavelength sampling is not linear, therefore the wavelength array must be stored entirely.
def save_fits_2d(wavelength, flux, file_name):

    """
    This function saves the processed spectra with non-linear wavelength sampling.
    Input: wavelength array, flux array, name of the spectrum.
    Output: fits file with separate wavelength and flux columns.

    """

    # Creating a table with wavelength and flux values
    col1 = fits.Column(name='WAVELENGTH', array=wavelength, format='D')
    col2 = fits.Column(name='FLUX', array=flux, format='D')
    hdu_table = fits.BinTableHDU.from_columns([col1, col2])

    # Creating the header keywords
    primary_header = fits.Header()
    primary_header['SIMPLE'] = (True, 'conforms to FITS standard')
    primary_header['BITPIX'] = 8
    primary_header['NAXIS'] = 0
    primary_header['EXTEND'] = True
    primary_header['CREATOR'] = ('SPAN', 'Generated by SPAN')
    primary_header.add_comment('Generated by SPAN with non-linear wavelength sampling')
    primary_hdu = fits.PrimaryHDU(header=primary_header)

    # Storing header and data
    hdulist = fits.HDUList([primary_hdu, hdu_table])
    hdulist.writeto(file_name, overwrite=True)

    # closing the file
    hdulist.close()


#*************************************************************************************************
#Gaussian function
"""
The functions below are three gaussians. The first one is
a simple gaussian, the second is a gaussian convolved with
a line, and the third are three gaussians convolved with a
straigth line, used for the line(s) fitting task of the Spectral
Analysis frame.

"""

#7) simple gaussian
def Gauss(x, y0, x0, a, sigma):
    return y0 + a * np.exp(-(x - x0)**2 / (2 * sigma**2))

#8) Gaussian with slope
def Gauss_slope(x, y0, x0, a, sigma, m, c):
    return y0 + a * np.exp(-(x - x0)**2 / (2 * sigma**2)) + m * x + c

#9) Three gaussians with slope
def multiple_gauss(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 6):
        y0, x0, a, sigma, m, c = params[i:i+6]
        y += y0 + a * np.exp(-(x - x0)**2 / (2 * sigma**2)) + m * x + c
    return y


#*************************************************************************************************
# #10) Measure the resolution from an emission (sky) line
# Gaussian function: offset + amplitude * exp[-(x - mu)^2 / (2 * sigma^2)]
def Gauss_res(x, offset, mu, amp, sigma):
    return offset + amp * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))


def resolution(wavelength, flux, wave1, wave2):
    """
    Fit a Gaussian to an emission line and estimate spectral resolution R and FWHM, with uncertainty.

    Parameters:
        wavelength : array-like
            Wavelength array of the 1D spectrum.
        flux : array-like
            Flux array of the 1D spectrum.
        wave1 : float
            Minimum wavelength of the region containing the sky line.
        wave2 : float
            Maximum wavelength of the region containing the sky line.

    Returns:
        resolution_R : int
            Spectral resolution R = lambda / delta_lambda.
        fwhm : float
            Full width at half maximum of the line (in Angstrom).
        fwhm_err : float
            Uncertainty on FWHM.
        line_wave : array
            Wavelength array of the selected region.
        line_flux_spec_norm : array
            Normalised flux of the selected region.
        fitted_gauss : array
            Gaussian model evaluated over line_wave.
    """
    mask = (wavelength >= wave1) & (wavelength <= wave2)
    line_wave = wavelength[mask]
    line_flux_spec = flux[mask]

    # Normalise the flux to median
    line_flux_spec_norm = line_flux_spec / np.median(line_flux_spec)

    # Initial guess for Gaussian parameters
    offset_guess = 1.0
    amp_guess = np.max(line_flux_spec_norm) - offset_guess
    mu_guess = line_wave[np.argmax(line_flux_spec_norm)]
    sigma_guess = (wave2 - wave1) / 6

    try:
        popt, pcov = curve_fit(
            Gauss_res, line_wave, line_flux_spec_norm,
            p0=[offset_guess, mu_guess, amp_guess, sigma_guess],
            bounds=([0, wave1, 0, 0], [np.inf, wave2, np.inf, np.inf])
        )
        offset_fit, mu_fit, amp_fit, sigma_fit = popt
        sigma_err = np.sqrt(np.diag(pcov))[3]  # error on sigma

        fwhm = 2.355 * sigma_fit
        fwhm_err = 2.355 * sigma_err
        resolution_R = int(mu_fit / fwhm)

        print(f"Resolution in A (FWHM): {fwhm:.2f} ± {fwhm_err:.2f}")
        # print(f"Resolution R: {resolution_R}")

        return resolution_R, fwhm, fwhm_err, line_wave, line_flux_spec_norm, Gauss_res(line_wave, *popt)

    except RuntimeError:
        print("Fit did not converge.")
        return None, None, None, line_wave, line_flux_spec_norm, np.zeros_like(line_wave)


#*************************************************************************************************
#11) Convert flux from Jansky to f_lambda or f_nu
def convert_flux(wavelength, flux, spec_name, type_to_convert, lambda_units):

    """
    This function converts the flux from Jansky to density flux (lambda or nu).
    Useful to correctly display in flux density units the IRTF extended spectral
    library of Villaume et al. 2017.
    Input: wavelength array, flux array, path and name of the spectrum, type of
            conversion ('to_flambda' or 'to_fnu'), wavelength scale units ('nm', 'A', 'mu')
    Output: array containing the converted flux values
    """

    flux_points = len(flux)
    converted_flux = np.zeros(flux_points)

    if lambda_units == 'mu':
        wavelength /= 1000
    elif lambda_units == 'A':
        wavelength /= 10

    if type_to_convert == 'to_flambda':
        conversion_factor = 2.999792458e12
        converted_flux = flux * conversion_factor / (wavelength ** 2)
    elif type_to_convert == 'to_fnu':
        conversion_factor = 1e26
        converted_flux = flux * conversion_factor

    return converted_flux


#*************************************************************************************************
#12) Noise estimation of a spectrum
def noise_spec(flux):
    """
    Computes the noise of a spectrum given its flux.

    Parameters:
    -----------
    flux : array-like
        Array containing the flux values of the spectrum.

    Returns:
    --------
    noise : float
        The computed noise of the spectrum.
    """

    # Convert input to a NumPy array with float64 precision
    flux = np.array(flux, dtype=np.float64)

    # Exclude zero values (padded values)
    flux = flux[flux != 0.0]

    n = len(flux)

    # If the spectrum is too short, return 0
    if n < 5:
        return 0.0

    # Compute the signal as the median flux
    signal = np.nanmedian(flux)

    # Compute the noise using the given formula
    noise = 0.6052697 * np.nanmedian(np.abs(2.0 * flux[2:n-2] - flux[0:n-4] - flux[4:n]))

    # Avoid division by zero
    if noise == 0:
        return np.inf if signal > 0 else 0.0

    return float(noise)


#********************** END OF UTILITIES FUNCTIONS ****************************************
#******************************************************************************************
