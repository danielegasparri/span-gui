#SPectral ANalysis software (SPAN)
#Written by Daniele Gasparri#


"""
    Copyright (C) 2018-2024, Daniele Gasparri

    E-mail: daniele.gasparri@gmail.com

    SPAN is a GUI interface that allows to modify and analyse 1D spectra.
    It is provided "as is" without any warranty whatsoever.
    Permission to use for non-commercial purposes is granted.
    Permission to modify for personal or internal use is granted,
    provided that this copyright and disclaimer are included unchanged
    at the beginning of the file. All other rights are reserved.
    In particular, redistribution of the code is not allowed.

"""

#******************************************************************************************
#******************************************************************************************
#*************************** UTILITIES FUNCTIONS FOR SPAN *********************************
#******************************************************************************************
#******************************************************************************************


#Local imports
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
    Output: float step value (in nm), bool linear (True) or log (False) sampling
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
    snr_ang = snr_pix * np.sqrt(1 / (step * 10))  # Supposing all the units in nm!

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

    if lambda_units == 'a':
        wavelength = wavelength*10
    if lambda_units == 'mu':
        wavelength = wavelength/1000

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


#10) Measure the resolution from an emission (sky) line
def resolution(wavelength, flux, wave1, wave2):

    """
    This function fits a gaussian to a sky emission line of the
    1D spectrum and gives the resolution, in R and delta lambda (FWHM)
    Input: wavelength array, flux array, min (wave1) and max (wave2) wavelength that
    identify the range to fit
    Output: int resolution R, arrays of wavelength and flux of the selected wavelength
            range, moments of the fitted gaussian
    """

    step = wavelength[1] - wavelength[0]
    mask = (wavelength >= wave1) & (wavelength <= wave2)
    line_wave = wavelength[mask]
    line_flux_spec = flux[mask]

    wave_norm = line_wave[10]  # guessing for now. fix it later
    epsilon_norm = step * 10
    line_flux_spec_norm = spman.norm_spec(line_wave, line_flux_spec, wave_norm, epsilon_norm, line_flux_spec)

    max_flux_spec = np.argmax(line_flux_spec_norm)
    mean_spec = line_wave[max_flux_spec]
    sigma_spec = 0.1
    offset = 1.

    popt_spec, _ = curve_fit(Gauss, line_wave, line_flux_spec_norm, p0=[offset, mean_spec, max_flux_spec, sigma_spec])
    sigma_fit_spec = popt_spec[3]

    resolution_fwhm = sigma_fit_spec * 2.35
    resolution_R = int(line_wave[0] / resolution_fwhm)

    print('Resolution in A (FWHM): ', round(resolution_fwhm * 10, 2))
    return resolution_R, line_wave, line_flux_spec_norm, Gauss(line_wave, *popt_spec)



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






#********************** END OF UTILITIES FUNCTIONS ****************************************
#******************************************************************************************
