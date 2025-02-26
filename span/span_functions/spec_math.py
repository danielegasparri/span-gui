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
#*************************** SPECTRA MATH FUNCTIONS FOR SPAN ******************************
#******************************************************************************************
#******************************************************************************************

try:#Local imports
    from span_functions import spec_manipul as spman
    from span_functions import system_span as stm

except ModuleNotFoundError: #local import if executed as package
    from . import spec_manipul as spman
    from . import system_span as stm


#Python imports
import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.interpolate import interp1d




#*************************************************************************************************
#1) Average spectra
def average(lambda_units, spectra_number, spec_names):

    """
    This function compiles a simple average of all the spectra loaded
    into SPAN by interpolating to a common wavelength grid before averaging them
    Input: wavelength units ('nm', 'a', 'mu'), int number of spectra,
           path and names of all the spectra.
    Output: 2D array containing the wavelength and the averaged flux

    """

    if spectra_number == 1:
        print('Just one spectra, cannot average anything!')
        return 0

    #Reading the spectra, resample to a common, linear value, then average the fluxes
    for i in range (spectra_number):
        wavelength, flux, step, name = stm.read_spec(spec_names[i], lambda_units)
        if i == 0:
            wavelength_grid = wavelength
            shape = (len(wavelength), spectra_number)
            data_to_df = np.zeros(shape)
            df_to_avg = pd.DataFrame(data_to_df)
            df_to_avg.iloc[:,i] = flux
        if i > 0:

            #Interpolating on the wavelength grid of the first spectrum and saving the fluxes to a file
            interpfunc = interpolate.interp1d(wavelength, flux, kind = 'linear', fill_value='extrapolate')
            interp_flux = (interpfunc(wavelength_grid))
            df_to_avg.iloc[:,i] = interp_flux
    average = df_to_avg.mean(axis=1) #average all over the columns
    average_flux = average.to_numpy(dtype = float)
    average_spec = np.column_stack((wavelength_grid, average_flux))
    return average_spec


#*************************************************************************************************
#2) Normalise and average spectra
def average_norm(lambda_units, wavelength, flux, spectra_number, spec_names):

    """
    This function normalise and average of all the spectra loaded
    into SPAN by interpolating to a common wavelength grid and normalising to the mean
    wavelength value of the selected spectrum before averaging them
    Input: wavelength units ('nm', 'a', 'mu'), wavelength array, flux array of the
           reference spectrum, int number of spectra, path and names of all the spectra.
    Output: 2D array containing the wavelength and the normalised and averaged flux

    """

    if spectra_number == 1:
        print('Just one spectra, cannot average anything!')
        return 0

    wavelength_grid = wavelength
    npoints = len(wavelength_grid)
    step = wavelength_grid[1]-wavelength_grid[0]
    wave_norm = (wavelength_grid[0] + wavelength_grid[npoints-1])/2.
    epsilon_norm = 10*step

    #reference spectrum is the one selected in span:
    norm_flux_reference_spec = spman.norm_spec(wavelength, flux, wave_norm, epsilon_norm, flux)

    shape = (len(wavelength_grid), spectra_number)
    data_to_df = np.zeros(shape)
    df_to_avg = pd.DataFrame(data_to_df)

    for i in range (spectra_number):
        wavelength, flux, step, name = stm.read_spec(spec_names[i], lambda_units)

        #Interpolating on the wavelength grid of the first spectrum and saving the fluxes to a file
        interpfunc = interpolate.interp1d(wavelength, flux, kind = 'linear', fill_value='extrapolate')
        interp_flux = (interpfunc(wavelength_grid))
        norm_interp_flux = spman.norm_spec(wavelength_grid, interp_flux, wave_norm, epsilon_norm, interp_flux)
        df_to_avg.iloc[:,i] = norm_interp_flux

    average = df_to_avg.mean(axis=1) #average all over the columns
    average_flux = average.to_numpy(dtype = float)
    average_norm_spec = np.column_stack((wavelength_grid, average_flux))
    return average_norm_spec


#*************************************************************************************************
#3) Sum the spectra
def sum_spec(lambda_units, spectra_number, spec_names):


    """
    This function compiles the sum of all the spectra loaded
    into SPAN by interpolating to a common wavelength grid before averaging them
    Input: wavelength units ('nm', 'a', 'mu'), int number of spectra,
           path and names of all the spectra.
    Output: 2D array containing the wavelength and the summed flux

    """

    if spectra_number == 1:
        print('Just one spectra, cannot average anything!')
        return 0

    #Reading the spectra, resample to a common, linear value, then average the fluxes
    for i in range (spectra_number):
        wavelength, flux, step, name = stm.read_spec(spec_names[i], lambda_units)
        if i == 0:
            wavelength_grid = wavelength
            shape = (len(wavelength), spectra_number)
            data_to_df = np.zeros(shape)
            df_to_avg = pd.DataFrame(data_to_df)
            df_to_avg.iloc[:,i] = flux

        if i > 0:
            #Interpolating on the wavelength grid of the first spectrum and saving the fluxes to a file
            interpfunc = interpolate.interp1d(wavelength, flux, kind = 'linear', fill_value='extrapolate')
            interp_flux = (interpfunc(wavelength_grid))
            df_to_avg.iloc[:,i] = interp_flux

    sum_all = df_to_avg.sum(axis=1) #average all over the columns
    sum_flux = sum_all.to_numpy(dtype = float)
    sum_spec = np.column_stack((wavelength_grid, sum_flux))
    return sum_spec


#*************************************************************************************************
#4) Normalize and sum the spectra
def sum_norm_spec(lambda_units, spectra_number, spec_names):

    """
    This function normalise and the sum all the spectra loaded
    into SPAN by interpolating to a common wavelength grid before normalising
    to the mean wavelength value of the first spectrum and summing them.
    Input: wavelength units ('nm', 'a', 'mu'), int number of spectra,
           path and names of all the spectra.
    Output: 2D array containing the wavelength and the normalised and summed flux

    """

    if spectra_number == 1:
        print('Just one spectra, cannot average anything!')
        return 0

    #Reading the spectra, resample to a common, linear value, then average the fluxes
    for i in range (spectra_number):
        wavelength, flux, step, name = stm.read_spec(spec_names[i], lambda_units)
        epsilon_norm = 10*step
        npoints = len(wavelength)

        if i == 0:
            wavelength_grid = wavelength
            shape = (len(wavelength), spectra_number)

            #find the normalization wavelength = half way to the beginning and the end of the spectrum
            wave_norm = (wavelength_grid[0] + wavelength_grid[npoints-1])/2.
            norm_flux = spman.norm_spec(wavelength_grid, flux, wave_norm, epsilon_norm, flux)

            #storing data
            data_to_df = np.zeros(shape)
            df_to_avg = pd.DataFrame(data_to_df)
            df_to_avg.iloc[:,i] = norm_flux


        if i > 0:
            #finding the flux at the wavelength_grid location
            interpfunc = interpolate.interp1d(wavelength, flux, kind = 'linear', fill_value='extrapolate')
            interp_flux = (interpfunc(wavelength_grid))

            #normalizing the flux
            norm_flux = spman.norm_spec(wavelength, interp_flux, wave_norm, epsilon_norm, interp_flux)

            #storing data
            df_to_avg.iloc[:,i] = norm_flux

    sum_norm_all = df_to_avg.sum(axis=1) #average all over the columns
    sum_norm_flux = sum_norm_all.to_numpy(dtype = float)
    sum_norm_spec = np.column_stack((wavelength_grid, sum_norm_flux))
    return sum_norm_spec


#*************************************************************************************************
#5) Subtract normalized average

def sub_norm_avg(wavelength, flux, lambda_units, spectra_number, spec_names):

    """
    This function normalises the selected spectrum in SPAN, then normalise
    and average all the spectra and subtract them to the selected one.
    Input: wavelength array, flux array, wavelength units ('nm', 'a', 'mu'),
           int number of spectra, path and names of all the spectra.
    Output: array containing the new flux values

    """

    #averaging all spectra
    average_norm_spec = average_norm(lambda_units, wavelength, flux, spectra_number, spec_names)
    average_wavelength = average_norm_spec[:,0]
    average_flux = average_norm_spec[:,1]

    # interpolating the average flux to the wavelength range of the spectrum
    interpfunc = interpolate.interp1d(average_wavelength, average_flux, kind = 'linear', fill_value='extrapolate')
    interp_average_flux = (interpfunc(wavelength))

    #normalize the spectra
    npoints = len(wavelength)
    wave_normalization = (wavelength[0] + wavelength[npoints-1])/2.
    epsilon_norm = (wavelength[1]-wavelength[0])*10
    norm_flux = spman.norm_spec(wavelength, flux, wave_normalization, epsilon_norm, flux)
    norm_interp_average = spman.norm_spec(wavelength, interp_average_flux, wave_normalization,epsilon_norm, interp_average_flux)
    subtracted_flux = norm_flux - norm_interp_average
    return subtracted_flux


#*************************************************************************************************
#6) Subtract normalised spec
def sub_norm_single(wavelength, flux, spectrum_to_subtract, lambda_units):

    """
    This function norlamise the selected spectrum in SPAN, then subtract
    the normalised spectrum selected by the user.
    IMPORTANT: the normalised spectrum to subtract must have the same wavelength
    units of the considered spectrum.
    Input: wavelength array, flux array, path and name of the spectrum to normalise and subtract,
           wavelength units ('nm', 'a', 'mu').
    Output: array containing the new flux values

    """

    wave_to_sub, flux_to_sub, step, name = stm.read_spec(spectrum_to_subtract, lambda_units)

    # interpolating the average flux to the wavelength range of the spectrum
    interpfunc = interpolate.interp1d(wave_to_sub, flux_to_sub, kind = 'linear', fill_value='extrapolate')
    interp_flux_to_sub = (interpfunc(wavelength))

    #normalise the spectra
    npoints = len(wavelength)
    wave_normalization = (wavelength[0] + wavelength[npoints-1])/2.
    epsilon_norm = (wavelength[1]-wavelength[0])*10

    #normalisation
    norm_flux = spman.norm_spec(wavelength, flux, wave_normalization, epsilon_norm, flux)
    norm_interp_sub = spman.norm_spec(wavelength, interp_flux_to_sub, wave_normalization, epsilon_norm, interp_flux_to_sub)
    subtracted_flux = norm_flux - norm_interp_sub
    return subtracted_flux


#********************** END OF SPECTRA MATH FUNCTIONS *************************************
#******************************************************************************************
