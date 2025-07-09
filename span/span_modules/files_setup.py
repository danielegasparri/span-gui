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

# Functions to set up the ASCII files written for the spectral analysis tasks and the 'Process all' mode

try: #try local import if executed as script
    #GUI import
    from FreeSimpleGUI_local import FreeSimpleGUI as sg
    from span_functions import linestrength as ls

except ModuleNotFoundError: #local import if executed as package
    #GUI import
    from span.FreeSimpleGUI_local import FreeSimpleGUI as sg
    from span.span_functions import linestrength as ls


import numpy as np
import pandas as pd
import os

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(CURRENT_DIR)


def create_blackbody_file(result_bb_dir, spectra_list_name, timestamp, spectra_number, spec_names_nopath):

    """
    Creates and saves a blackbody data file containing temperature values for each spectrum

    """

    # Define the output file path
    bb_file = os.path.join(result_bb_dir, f"{spectra_list_name}_bb_data_{timestamp}.dat")

    # Define column headers
    bb_id = ['#Spectrum', 'T(K)', 'err', 'chi2']

    # Initialise an array with zero values
    bb_values = np.zeros(spectra_number)
    bb_err = np.zeros(spectra_number)
    bb_chi2 = np.zeros(spectra_number)
    bb_data_array = np.column_stack((spec_names_nopath, bb_values, bb_err, bb_chi2))

    # Generate a DataFrame and add the data
    df_bb = pd.DataFrame(bb_data_array, columns=bb_id)

    # Save the DataFrame to a file
    df_bb.to_csv(bb_file, index=True, sep=' ')

    return bb_file



def create_cross_correlation_file(result_xcorr_dir, spectra_list_name, timestamp, spectra_number, spec_names_nopath):

    """
    Creates and saves a cross-correlation results file containing radial velocity (RV) values for each spectrum.

    """

    # Define the output file path
    rv_file = os.path.join(result_xcorr_dir, f"{spectra_list_name}_rv_data_{timestamp}.dat")

    # Define column headers
    rv_id = ['#Spectrum', 'RV(km/s)', 'err']

    # Initialise an array with zero values
    rv_values = np.zeros(spectra_number)
    error_values = np.zeros(spectra_number)
    rv_data_array = np.column_stack((spec_names_nopath, rv_values, error_values))

    # Generate a DataFrame and add the data
    df_rv = pd.DataFrame(rv_data_array, columns=rv_id)

    # Save the DataFrame to a file
    df_rv.to_csv(rv_file, index=True, sep=' ')

    return rv_file



def create_velocity_dispersion_file(result_vel_disp_dir, spectra_list_name, timestamp, spectra_number, spec_names_nopath):

    """
    Creates and saves a velocity dispersion results file containing sigma values and associated errors for each spectrum.

    """

    # Define the output file path
    sigma_file = os.path.join(result_vel_disp_dir, f"{spectra_list_name}_sigma_data_{timestamp}.dat")

    # Define column headers
    sigma_id = ['#Spectrum', 'Sigma(km/s)', 'err']

    # Initialise arrays with zero values for sigma and error
    sigma_values = np.zeros(spectra_number)
    err_values = np.zeros(spectra_number)
    sigma_data_array = np.column_stack((spec_names_nopath, sigma_values, err_values))

    # Generate a DataFrame and add the data
    df_sigma = pd.DataFrame(sigma_data_array, columns=sigma_id)

    # Save the DataFrame to a file
    df_sigma.to_csv(sigma_file, index=True, sep=' ')

    return sigma_file



def create_ew_measurement_files(result_ew_data_dir, spectra_list_name, timestamp, spectra_number, spec_names_nopath):

    """
    Creates and saves equivalent width (EW) measurement files, including EW in Angstroms, EW in magnitudes, and SNR data.

    """

    # Define output file paths
    ew_file = os.path.join(result_ew_data_dir, f"{spectra_list_name}_ew_data_{timestamp}.dat")
    ew_file_mag = os.path.join(result_ew_data_dir, f"{spectra_list_name}_ew_data_mag_{timestamp}.dat")
    snr_ew_file = os.path.join(result_ew_data_dir, f"{spectra_list_name}_snr_ew_data_{timestamp}.dat")

    # Define column headers
    ew_id = ['#Spectrum', 'ew(A)', 'err']
    ew_id_mag = ['#Spectrum', 'ew(Mag)', 'err']
    snr_ew_id = ['#Spectrum', 'SNR']

    # Initialise arrays with zero values for EW, errors, and SNR
    ew_values = np.zeros(spectra_number)
    err_ew_values = np.zeros(spectra_number)
    snr_ew_values = np.zeros(spectra_number)

    # Create data arrays
    ew_data_array = np.column_stack((spec_names_nopath, ew_values, err_ew_values))
    ew_data_array_mag = np.column_stack((spec_names_nopath, ew_values, err_ew_values))
    snr_ew_data_array = np.column_stack((spec_names_nopath, snr_ew_values))

    # Generate DataFrames and add the data
    df_ew = pd.DataFrame(ew_data_array, columns=ew_id)
    df_ew_mag = pd.DataFrame(ew_data_array_mag, columns=ew_id_mag)
    df_snr_ew = pd.DataFrame(snr_ew_data_array, columns=snr_ew_id)

    # Save DataFrames to files
    df_ew.to_csv(ew_file, index=True, sep=' ')
    df_ew_mag.to_csv(ew_file_mag, index=True, sep=' ')
    df_snr_ew.to_csv(snr_ew_file, index=True, sep=' ')

    return ew_file, ew_file_mag, snr_ew_file



def create_ssp_lick_param_file(result_ew_data_dir, spectra_list_name, timestamp, spectra_number, spec_names_nopath):

    """
    Creates and saves the Lick/IDS index file used as stellar population diagnostics

    """

    ssp_lick_param_file = os.path.join(result_ew_data_dir, f"{spectra_list_name}_lick_for_ssp_{timestamp}.dat")

    ssp_lick_param_id = ['#Spectrum', 'Hbeta(A)', 'Hbeta_err(A)', 'Mg2(mag)', 'Mg2_err(mag)', 'Mgb(A)', 'Mgb_err(A)', 'Fe5270(A)', 'Fe5270_err(A)', 'Fe5335(A)', 'Fe5335_err(A)', 'Fem(A)', 'Fem_err(A)', 'MgFe(A)', 'MgFe_err(A)']

    hbeta_values = np.zeros(spectra_number)
    hbetae_values = np.zeros(spectra_number)
    mg2_values = np.zeros(spectra_number)
    mg2e_values = np.zeros(spectra_number)
    mgb_values = np.zeros(spectra_number)
    mgbe_values = np.zeros(spectra_number)
    fe5270_values = np.zeros(spectra_number)
    fe5270e_values = np.zeros(spectra_number)
    fe5335_values = np.zeros(spectra_number)
    fe5335e_values = np.zeros(spectra_number)
    fem_values = np.zeros(spectra_number)
    feme_values = np.zeros(spectra_number)
    mgfe_values = np.zeros(spectra_number)
    mgfee_values = np.zeros(spectra_number)

    ssp_lick_param_array = np.column_stack((spec_names_nopath, hbeta_values, hbetae_values, mg2_values, mg2e_values, mgb_values, mgbe_values, fe5270_values, fe5270e_values, fe5335_values, fe5335e_values, fem_values, feme_values, mgfe_values, mgfee_values))

    df_lick_param = pd.DataFrame(ssp_lick_param_array, columns=ssp_lick_param_id)

    df_lick_param.to_csv(ssp_lick_param_file, index=True, sep=' ')

    return ssp_lick_param_file



def create_ew_measurement_files_from_index(result_ew_data_dir, spectra_list_name, timestamp, spectra_number, spec_names_nopath, index_file):

    """
    Creates and saves equivalent width (EW) measurement files using a predefined index file.

    """

    # Check if the index file exists
    if not os.path.isfile(index_file):
        print('The index file does not exist. Skipping...')
        return None, None, None, None, None, None, None, None

    # Attempt to read the index file
    try:
        idx_names, indices = ls.read_idx(index_file)
    except ValueError:
        print('At least one index in the file is not valid')
        return None, None, None, None, None, None, None, None

    # Ensure the index file has the correct format
    if len(indices[:, 0]) < 6:
        print('The length of at least one index is not correct')
        return None, None, None, None, None, None, None, None

    # Validate index wavelength ranges
    bad_idx = []
    for t in range(len(idx_names)):
        if (indices[0, t] > indices[1, t] or
            indices[2, t] > indices[3, t] or
            indices[4, t] > indices[5, t]):
            bad_idx.append(idx_names[t])

    if bad_idx:
        print('It seems we have a problem. Did you invert the wavelengths of these indices?', bad_idx)
        return None, None, None, None, None, None, None, None

    # Define output file paths
    ew_file = os.path.join(result_ew_data_dir, f"{spectra_list_name}_ew_data_{timestamp}.dat")
    ew_file_mag = os.path.join(result_ew_data_dir, f"{spectra_list_name}_ew_data_mag_{timestamp}.dat")
    snr_ew_file = os.path.join(result_ew_data_dir, f"{spectra_list_name}_snr_ew_data_{timestamp}.dat")

    # Read index data
    id_array, index = ls.read_idx(index_file)
    num_indices = len(id_array)

    # Initialise zero arrays for EW, errors, and SNR
    shape = (spectra_number, num_indices)
    ew_all = np.zeros(shape)
    ew_all_mag = np.zeros(shape)
    snr_ew_all = np.zeros(shape)
    err_all = np.zeros(shape)

    # Define column headers
    spectra_id = ['#Spectrum']

    # Correct way to create a string array for error suffix
    err_col_type = np.full(num_indices, 'e', dtype='<U1')  # Fixed dtype issue

    # Merge index names with error columns
    err_col_names = np.char.add(id_array, err_col_type)  # Now both arrays have the same dtype
    col_names = np.concatenate((id_array, err_col_names))
    ew_id = np.concatenate((spectra_id, col_names))
    ew_data = np.column_stack((spec_names_nopath, ew_all, err_all))
    ew_id_mag = np.concatenate((spectra_id, col_names))
    ew_data_mag = np.column_stack((spec_names_nopath, ew_all_mag, err_all))

    # Create and save DataFrames for EW and EW magnitude
    df_ew = pd.DataFrame(ew_data, columns=ew_id)
    df_ew.to_csv(ew_file, index=True, sep=' ')

    df_ew_mag = pd.DataFrame(ew_data_mag, columns=ew_id_mag)
    df_ew_mag.to_csv(ew_file_mag, index=True, sep=' ')

    # Create and save DataFrame for SNR
    snr_col_names = id_array
    snr_ew_id = np.concatenate((spectra_id, snr_col_names))
    snr_ew_data = np.column_stack((spec_names_nopath, snr_ew_all))
    df_snr_ew = pd.DataFrame(snr_ew_data, columns=snr_ew_id)
    df_snr_ew.to_csv(snr_ew_file, index=True, sep=' ')

    return ew_file, ew_file_mag, snr_ew_file, num_indices, ew_id, spectra_id, ew_id_mag, snr_ew_id



def create_lick_ew_measurement_files(result_ew_data_dir, spectra_list_name, timestamp, spectra_number, spec_names_nopath, lick_index_file):

    """
    Creates and saves equivalent width (EW) measurement files for Lick/IDS indices.

    """

    # Check if the Lick index file exists
    if not os.path.isfile(lick_index_file):
        sg.popup('The index file does not exist. Skipping...')
        return None, None, None

    # Read the index data
    id_lick_array, lick_index = ls.read_idx(lick_index_file)
    num_lick_indices = len(id_lick_array)

    # Define output file paths
    ew_lick_file = os.path.join(result_ew_data_dir, f"{spectra_list_name}_ew_lick_data_{timestamp}.dat")
    ew_lick_file_mag = os.path.join(result_ew_data_dir, f"{spectra_list_name}_ew_lick_data_mag_{timestamp}.dat")
    snr_lick_ew_file = os.path.join(result_ew_data_dir, f"{spectra_list_name}_snr_ew_lick_data_{timestamp}.dat")

    # Initialise zero arrays for EW, errors, and SNR
    shape = (spectra_number, num_lick_indices)
    ew_lick_all = np.zeros(shape)
    ew_lick_all_mag = np.zeros(shape)
    snr_lick_ew_all = np.zeros(shape)
    err_lick_all = np.zeros(shape)

    # Define column headers
    spectra_lick_id = ['#Spectrum']

    # Correct way to create a string array for error suffix
    err_lick_col_type = np.full(num_lick_indices, 'e', dtype='<U1')  # Fixed dtype issue

    # Merge index names with error columns
    err_lick_col_names = np.char.add(id_lick_array, err_lick_col_type)  # Now both arrays have the same dtype
    col_lick_names = np.concatenate((id_lick_array, err_lick_col_names))
    ew_lick_id = np.concatenate((spectra_lick_id, col_lick_names))
    ew_lick_data = np.column_stack((spec_names_nopath, ew_lick_all, err_lick_all))
    ew_lick_id_mag = np.concatenate((spectra_lick_id, col_lick_names))
    ew_lick_data_mag = np.column_stack((spec_names_nopath, ew_lick_all_mag, err_lick_all))

    # Create and save DataFrames for EW and EW magnitude
    df_ew_lick = pd.DataFrame(ew_lick_data, columns=ew_lick_id)
    df_ew_lick.to_csv(ew_lick_file, index=True, sep=' ')

    df_ew_lick_mag = pd.DataFrame(ew_lick_data_mag, columns=ew_lick_id_mag)
    df_ew_lick_mag.to_csv(ew_lick_file_mag, index=True, sep=' ')

    # Create and save DataFrame for SNR
    snr_lick_col_names = id_lick_array
    snr_lick_ew_id = np.concatenate((spectra_lick_id, snr_lick_col_names))
    snr_lick_ew_data = np.column_stack((spec_names_nopath, snr_lick_ew_all))
    df_snr_lick_ew = pd.DataFrame(snr_lick_ew_data, columns=snr_lick_ew_id)
    df_snr_lick_ew.to_csv(snr_lick_ew_file, index=True, sep=' ')

    return ew_lick_file, ew_lick_file_mag, snr_lick_ew_file, num_lick_indices, ew_lick_id, spectra_lick_id, ew_lick_id_mag, snr_lick_ew_id



def create_lick_ssp_parameters_file(result_ew_data_dir, spectra_list_name, timestamp, spectra_number, spec_names_nopath):

    """
    Creates and saves the stellar population parameters file for Lick/IDS index analysis.

    """

    # Define the output file path
    ssp_param_file = os.path.join(result_ew_data_dir, f"{spectra_list_name}_ssp_param_lick_{timestamp}.dat")

    # Define column headers
    ssp_param_id = ['#Spectrum', 'age(Gyr)', 'err_age', 'met', 'err_met', 'alpha', 'err_alpha']

    # Initialise zero arrays for age, metallicity, and alpha enhancement with their errors
    age_ssp = np.zeros(spectra_number)
    age_err_ssp = np.zeros(spectra_number)
    met_ssp = np.zeros(spectra_number)
    met_err_ssp = np.zeros(spectra_number)
    alpha_ssp = np.zeros(spectra_number)
    alpha_err_ssp = np.zeros(spectra_number)

    # Combine the data into a structured array
    ssp_param_array = np.column_stack((spec_names_nopath, age_ssp, age_err_ssp, met_ssp, met_err_ssp, alpha_ssp, alpha_err_ssp))

    # Create a DataFrame and save it to a file
    df_ssp_param = pd.DataFrame(ssp_param_array, columns=ssp_param_id)
    df_ssp_param.to_csv(ssp_param_file, index=True, sep=' ')

    return ssp_param_file



def create_line_fitting_file(result_line_fitting_dir, spectra_list_name, timestamp, spectra_number, spec_names_nopath):

    """
    Creates and saves the line fitting results file.

    """

    # Define the output file path
    fit_file = os.path.join(result_line_fitting_dir, f"{spectra_list_name}_fit_data_{timestamp}.dat")

    # Define column headers
    fit_id = ['#Spectrum', 'ca1_wave', 'ca2_wave', 'ca3_wave', 'dw_ca1', 'dw_ca2', 'dw_ca3', 'ew_ca1', 'ew_ca2', 'ew_ca3']

    # Initialise zero arrays for CaT line properties
    ca1_wave = np.zeros(spectra_number)
    ca2_wave = np.zeros(spectra_number)
    ca3_wave = np.zeros(spectra_number)
    dw_ca1 = np.zeros(spectra_number)
    dw_ca2 = np.zeros(spectra_number)
    dw_ca3 = np.zeros(spectra_number)
    ew_ca1 = np.zeros(spectra_number)
    ew_ca2 = np.zeros(spectra_number)
    ew_ca3 = np.zeros(spectra_number)

    # Combine the data into a structured array
    fit_data_array = np.column_stack((spec_names_nopath, ca1_wave, ca2_wave, ca3_wave, dw_ca1, dw_ca2, dw_ca3, ew_ca1, ew_ca2, ew_ca3))

    # Create a DataFrame and save it to a file
    df_fit = pd.DataFrame(fit_data_array, columns=fit_id)
    df_fit.to_csv(fit_file, index=True, sep=' ')

    return fit_file



def create_line_fitting_file_simple(result_line_fitting_dir, spectra_list_name, timestamp, spectra_number, spec_names_nopath):

    """
    Creates and saves a simplified line fitting results file (single line case).

    """

    # Define the output file path
    fit_file = os.path.join(result_line_fitting_dir, f"{spectra_list_name}_fit_data_{timestamp}.dat")

    # Define column headers
    fit_id = ['#Spectrum', 'line_wave']

    # Initialise an array for the fitted wavelength values
    fit_values = np.zeros(spectra_number)

    # Combine the data into a structured array
    fit_data_array = np.column_stack((spec_names_nopath, fit_values))

    # Create a DataFrame and save it to a file
    df_fit = pd.DataFrame(fit_data_array, columns=fit_id)
    df_fit.to_csv(fit_file, index=True, sep=' ')

    return fit_file



def create_kinematics_files(result_ppxf_kin_data_dir, spectra_list_name, timestamp, spectra_number, spec_names_nopath, with_errors_kin, gas_kin):

    """
    Creates and saves kinematics results files for stellar and gas kinematics.

    """

    kinematics_files = {}

    # Define the output file path for stellar kinematics
    kin_file = os.path.join(result_ppxf_kin_data_dir, f"{spectra_list_name}_kin_data_stars_{timestamp}.dat")

    # Define column headers for stellar kinematics
    kin_id = ['#Spectrum', 'RV(km/s)', 'Sigma(km/s)', 'H3', 'H4', 'H5', 'H6', 'errRV', 'errSigma', 'errH3', 'errH4', 'errH5', 'errH6', 'S/N']

    # Initialise arrays with zero values
    kinematics_data = np.zeros((spectra_number, len(kin_id) - 1))

    # Combine data into a structured array
    kin_data_array = np.column_stack((spec_names_nopath, kinematics_data))

    # Create a DataFrame and save it to a file
    df_kin = pd.DataFrame(kin_data_array, columns=kin_id)
    df_kin.to_csv(kin_file, index=True, sep=' ')
    kinematics_files["stellar"] = kin_file

    # Monte Carlo error estimation file
    if with_errors_kin:
        kin_file_mc = os.path.join(result_ppxf_kin_data_dir, f"{spectra_list_name}_kin_data_stars_mc_errors_{timestamp}.dat")

        # Create the Monte Carlo error file with the same format
        df_kin.to_csv(kin_file_mc, index=True, sep=' ')
        kinematics_files["stellar_mc"] = kin_file_mc

    # # Gas kinematics file
    if gas_kin:
        kin_file_gas = os.path.join(result_ppxf_kin_data_dir, f"{spectra_list_name}_kin_data_gas_{timestamp}.dat")
        kinematics_files["gas"] = kin_file_gas

    return kinematics_files



def create_stellar_population_files(result_ppxf_pop_data_dir, spectra_list_name, timestamp, spectra_number, spec_names_nopath, ppxf_pop_lg_age, stellar_parameters_lick_ppxf):

    """
    Creates and saves stellar population analysis result files.

    """

    pop_files = {}

    # Define the output file path for stellar populations
    pop_file = os.path.join(result_ppxf_pop_data_dir, f"{spectra_list_name}_pop_data_{timestamp}.dat")

    # Define column headers based on whether log age or linear age is used
    if ppxf_pop_lg_age:
        pop_id = ['#Spectrum', 'RV(km/s)', 'Sigma(km/s)', 'H3', 'H4', 'lum_lg_age(dex)', 'lum_met(dex)', 'lum_alpha(dex)',
                  'err_lum_lg_age(dex)', 'err_lum_met(dex)',
                  'err_lum_alpha(dex)', 'M/L', 'mass_lg_age(dex)', 'mass_met(dex)', 'mass_alpha(dex)',
                  'err_mass_lg_age(dex)', 'err_mass_met(dex)',
                  'err_mass_alpha(dex)', 't50_age', 't80_age', 't50_cosmic', 't80_cosmic', 'Chi2', 'S/N']
    else:
        pop_id = ['#Spectrum', 'RV(km/s)', 'Sigma(km/s)', 'H3', 'H4', 'lum_age(Gyr)', 'lum_met(dex)', 'lum_alpha(dex)',
                  'err_lum_age(Gyr)', 'err_lum_met(dex)',
                  'err_lum_alpha(dex)', 'M/L', 'mass_age(Gyr)', 'mass_met(dex)', 'mass_alpha(dex)',
                  'err_mass_age(Gyr)', 'err_mass_met(dex)',
                  'err_mass_alpha(dex)', 't50_age', 't80_age', 't50_cosmic', 't80_cosmic', 'Chi2', 'S/N']

    # Initialise arrays with zero values
    population_data = np.zeros((spectra_number, len(pop_id) - 1))

    # Combine data into a structured array
    pop_data_array = np.column_stack((spec_names_nopath, population_data))

    # Create a DataFrame and save it to a file
    df_pop = pd.DataFrame(pop_data_array, columns=pop_id)
    df_pop.to_csv(pop_file, index=True, sep=' ')
    pop_files["stellar_population"] = pop_file

    # Creating the Lick/IDS stellar parameters file if activated
    if stellar_parameters_lick_ppxf:
        ssp_param_file_ppxf = os.path.join(result_ppxf_pop_data_dir, f"{spectra_list_name}_ssp_param_lick_ppxf_{timestamp}.dat")

        # Define headers for Lick/IDS stellar parameters
        ssp_param_id_ppxf = ['#Spectrum', 'Hbeta(A)', 'Hbeta_err(A)', 'Mg2(mag)', 'Mg2_err(mag)', 'Mgb(A)', 'Mgb_err(A)',
                             'Fem(A)', 'Fem_err(A)', 'MgFe(A)', 'MgFe_err(A)', 'age(Gyr)', 'err_age', 'met', 'err_met',
                             'alpha', 'err_alpha']

        # Initialise arrays with zero values
        ssp_param_data = np.zeros((spectra_number, len(ssp_param_id_ppxf) - 1))

        # Combine data into a structured array
        ssp_param_array_ppxf = np.column_stack((spec_names_nopath, ssp_param_data))

        # Create a DataFrame and save it to a file
        df_ssp_param_ppxf = pd.DataFrame(ssp_param_array_ppxf, columns=ssp_param_id_ppxf)
        df_ssp_param_ppxf.to_csv(ssp_param_file_ppxf, index=True, sep=' ')
        pop_files["ssp_lick_parameters"] = ssp_param_file_ppxf

    return pop_files
