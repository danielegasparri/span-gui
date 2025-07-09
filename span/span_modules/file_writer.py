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

# Functions to save the results of the Spectral analysis panel to ASCII files

import pandas as pd
import numpy as np

try: #try local import if executed as script
    #GUI import
    from params import SpectraParams
    from span_functions import spec_analysis as span

except ModuleNotFoundError: #local import if executed as package
    #GUI import
    from .params import SpectraParams
    from span.span_functions import spec_analysis as span

# from params import SpectraParams



def save_kinematics_to_file(i, params, kinematics, error_kinematics, error_kinematics_mc, gas_component, gas_names, gas_flux, gas_flux_err,
                            kin_component, snr_kin, df_kin, kin_file, df_kin_mc=None, kin_file_mc=None, df_kin_gas=None, kin_file_gas=None):

    """
    Saves kinematics data (stellar/gas) to file(s) using Pandas.

    """

    try:
        number_kin_component = np.max(kin_component)
        if number_kin_component == 0 and not params.ppxf_kin_two_stellar_components:
            vel = round(kinematics[0],3)
            sigma = round(kinematics[1],3)
            h3 = round(kinematics[2],5)
            h4 = round(kinematics[3],5)
            h5 = round(kinematics[4],5)
            h6 = round(kinematics[5],5)
            err_vel = round(error_kinematics[0],3)
            err_sigma = round(error_kinematics[1],3)
            err_h3 = round(error_kinematics[2],5)
            err_h4 = round(error_kinematics[3],5)
            err_h5 = round(error_kinematics[4],5)
            err_h6 = round(error_kinematics[5],5)
            #writing to file
            df_kin.at[i, 'RV(km/s)']= vel
            df_kin.at[i, 'Sigma(km/s)']= sigma
            df_kin.at[i, 'H3']= h3
            df_kin.at[i, 'H4']= h4
            df_kin.at[i, 'H5']= h5
            df_kin.at[i, 'H6']= h6
            df_kin.at[i, 'errRV']= err_vel
            df_kin.at[i, 'errSigma']= err_sigma
            df_kin.at[i, 'errH3']= err_h3
            df_kin.at[i, 'errH4']= err_h4
            df_kin.at[i, 'errH5']= err_h5
            df_kin.at[i, 'errH6']= err_h6
            df_kin.at[i, 'S/N']= round(snr_kin)

            df_kin.to_csv(kin_file, index= False, sep=' ')

            if params.with_errors_kin:
                err_rv_kin_mc, err_sigma_kin_mc, err_h3_kin_mc, err_h4_kin_mc, err_h5_kin_mc, err_h6_kin_mc = np.round(error_kinematics_mc[0],3)
                df_kin_mc.at[i, 'RV(km/s)']= vel
                df_kin_mc.at[i, 'Sigma(km/s)']= sigma
                df_kin_mc.at[i, 'H3']= h3
                df_kin_mc.at[i, 'H4']= h4
                df_kin_mc.at[i, 'H5']= h5
                df_kin_mc.at[i, 'H6']= h6
                df_kin_mc.at[i, 'errRV']= err_rv_kin_mc
                df_kin_mc.at[i, 'errSigma']= err_sigma_kin_mc
                df_kin_mc.at[i, 'errH3']= err_h3_kin_mc
                df_kin_mc.at[i, 'errH4']= err_h4_kin_mc
                df_kin_mc.at[i, 'errH5']= err_h5_kin_mc
                df_kin_mc.at[i, 'errH6']= err_h6_kin_mc
                df_kin_mc.at[i, 'S/N']= round(snr_kin)

                df_kin_mc.to_csv(kin_file_mc, index= False, sep=' ')

        #saving the two component stellar fit results
        elif number_kin_component == 0 and params.ppxf_kin_two_stellar_components:

            # adding the columns to the file
            if i == 0:
                new_component = ['RV_2(km/s)', 'Sigma_2(km/s)', 'H3_2', 'H4_2', 'H5_2', 'H6_2', 'errRV_2','errSigma_2', 'errH3_2','errH4_2', 'errH5_2', 'errH6_2']
                df_kin[new_component] = 0. #filling with zeros

            vel1 = round(kinematics[0][0],3)
            sigma1 = round(kinematics[0][1],3)
            h31 = round(kinematics[0][2],5)
            h41 = round(kinematics[0][3],5)
            h51 = round(kinematics[0][4],5)
            h61 = round(kinematics[0][5],5)
            err_vel1 = round(error_kinematics[0][0],3)
            err_sigma1 = round(error_kinematics[0][1],3)
            err_h31 = round(error_kinematics[0][2],5)
            err_h41 = round(error_kinematics[0][3],5)
            err_h51 = round(error_kinematics[0][4],5)
            err_h61 = round(error_kinematics[0][5],5)

            vel2 = round(kinematics[1][0],3)
            sigma2 = round(kinematics[1][1],3)
            h32 = round(kinematics[1][2],5)
            h42 = round(kinematics[1][3],5)
            h52 = round(kinematics[1][4],5)
            h62 = round(kinematics[1][5],5)
            err_vel2 = round(error_kinematics[1][0],3)
            err_sigma2 = round(error_kinematics[1][1],3)
            err_h32 = round(error_kinematics[1][2],5)
            err_h42 = round(error_kinematics[1][3],5)
            err_h52 = round(error_kinematics[1][4],5)
            err_h62 = round(error_kinematics[1][5],5)

            #filling the dataframe columns for component 1
            df_kin.at[i, 'RV(km/s)']= vel1
            df_kin.at[i, 'Sigma(km/s)']= sigma1
            df_kin.at[i, 'H3']= h31
            df_kin.at[i, 'H4']= h41
            df_kin.at[i, 'H5']= h51
            df_kin.at[i, 'H6']= h61
            df_kin.at[i, 'errRV']= err_vel1
            df_kin.at[i, 'errSigma']= err_sigma1
            df_kin.at[i, 'errH3']= err_h31
            df_kin.at[i, 'errH4']= err_h41
            df_kin.at[i, 'errH5']= err_h51
            df_kin.at[i, 'errH6']= err_h61
            df_kin.at[i, 'S/N']= round(snr_kin)

            #filling the dataframe columns for component 2
            df_kin.at[i, 'RV_2(km/s)']= vel2
            df_kin.at[i, 'Sigma_2(km/s)']= sigma2
            df_kin.at[i, 'H3_2']= h32
            df_kin.at[i, 'H4_2']= h42
            df_kin.at[i, 'H5_2']= h52
            df_kin.at[i, 'H6_2']= h62
            df_kin.at[i, 'errRV_2']= err_vel2
            df_kin.at[i, 'errSigma_2']= err_sigma2
            df_kin.at[i, 'errH3_2']= err_h32
            df_kin.at[i, 'errH4_2']= err_h42
            df_kin.at[i, 'errH5_2']= err_h52
            df_kin.at[i, 'errH6_2']= err_h62

            #writing to file
            df_kin.to_csv(kin_file, index= False, sep=' ')

            # considering also the errorrs with MonteCarlo simulations
            if params.with_errors_kin:
                #updating the dataframe with the second stellar component
                if i == 0:
                    new_component_mc = ['RV_2(km/s)', 'Sigma_2(km/s)', 'H3_2', 'H4_2', 'H5_2', 'H6_2', 'errRV_2','errSigma_2', 'errH3_2','errH4_2', 'errH5_2', 'errH6_2']
                    df_kin_mc[new_component_mc] = 0. #filling with zeros

                # extracting the MonteCarlo errors from the error array
                err_rv_kin_mc1, err_sigma_kin_mc1, err_h3_kin_mc1, err_h4_kin_mc1, err_h5_kin_mc1, err_h6_kin_mc1, err_rv_kin_mc2, err_sigma_kin_mc2, err_h3_kin_mc2, err_h4_kin_mc2, err_h5_kin_mc2, err_h6_kin_mc2  = np.round(error_kinematics_mc[0],5)

                # assigning to the dataframe the first component
                df_kin_mc.at[i, 'RV(km/s)']= vel1
                df_kin_mc.at[i, 'Sigma(km/s)']= sigma1
                df_kin_mc.at[i, 'H3']= h31
                df_kin_mc.at[i, 'H4']= h41
                df_kin_mc.at[i, 'H5']= h51
                df_kin_mc.at[i, 'H6']= h61
                df_kin_mc.at[i, 'errRV']= err_rv_kin_mc1
                df_kin_mc.at[i, 'errSigma']= err_sigma_kin_mc1
                df_kin_mc.at[i, 'errH3']= err_h3_kin_mc1
                df_kin_mc.at[i, 'errH4']= err_h4_kin_mc1
                df_kin_mc.at[i, 'errH5']= err_h5_kin_mc1
                df_kin_mc.at[i, 'errH6']= err_h6_kin_mc1
                df_kin_mc.at[i, 'S/N']= round(snr_kin)

                #assigning to the dataframe the second component
                df_kin_mc.at[i, 'RV_2(km/s)']= vel2
                df_kin_mc.at[i, 'Sigma_2(km/s)']= sigma2
                df_kin_mc.at[i, 'H3_2']= h32
                df_kin_mc.at[i, 'H4_2']= h42
                df_kin_mc.at[i, 'H5_2']= h52
                df_kin_mc.at[i, 'H6_2']= h62
                df_kin_mc.at[i, 'errRV_2']= err_rv_kin_mc2
                df_kin_mc.at[i, 'errSigma_2']= err_sigma_kin_mc2
                df_kin_mc.at[i, 'errH3_2']= err_h3_kin_mc2
                df_kin_mc.at[i, 'errH4_2']= err_h4_kin_mc2
                df_kin_mc.at[i, 'errH5_2']= err_h5_kin_mc2
                df_kin_mc.at[i, 'errH6_2']= err_h6_kin_mc2

                #writing the dataframe to file
                df_kin_mc.to_csv(kin_file_mc, index= False, sep=' ')

        #Saving the stellar and gas fit results
        else:
            vel = round(kinematics[0][0],3)
            sigma = round(kinematics[0][1],3)
            h3 = round(kinematics[0][2],5)
            h4 = round(kinematics[0][3],5)
            h5 = round(kinematics[0][4],5)
            h6 = round(kinematics[0][5],5)
            err_vel = round(error_kinematics[0][0],3)
            err_sigma = round(error_kinematics[0][1],3)
            err_h3 = round(error_kinematics[0][2],5)
            err_h4 = round(error_kinematics[0][3],5)
            err_h5 = round(error_kinematics[0][4],5)
            err_h6 = round(error_kinematics[0][5],5)

            df_kin.at[i, 'RV(km/s)']= vel
            df_kin.at[i, 'Sigma(km/s)']= sigma
            df_kin.at[i, 'H3']= h3
            df_kin.at[i, 'H4']= h4
            df_kin.at[i, 'H5']= h5
            df_kin.at[i, 'H6']= h6
            df_kin.at[i, 'errRV']= err_vel
            df_kin.at[i, 'errSigma']= err_sigma
            df_kin.at[i, 'errH3']= err_h3
            df_kin.at[i, 'errH4']= err_h4
            df_kin.at[i, 'errH5']= err_h5
            df_kin.at[i, 'errH6']= err_h6
            df_kin.at[i, 'S/N']= int(snr_kin)

            df_kin.to_csv(kin_file, index= False, sep=' ')


            #writing also the kin gas file. I need to create it here because I need to know the names of the columns
            if params.gas_kin:

                if df_kin_gas is None or not isinstance(df_kin_gas, pd.DataFrame) and i ==0:

                    kin_id_gas = ['#Spectrum']
                    for name in gas_names:
                        kin_id_gas += [
                            f'RV(km/s)_{name}', f'Sigma(km/s)_{name}',
                            f'H3_{name}', f'H4_{name}', f'H5_{name}', f'H6_{name}',
                            f'Flux_{name}', f'Flux_err_{name}',
                            f'errRV_{name}', f'errSigma_{name}',
                            f'errH3_{name}', f'errH4_{name}', f'errH5_{name}', f'errH6_{name}'
                        ]

                    spectra_number = len(params.spec_names_nopath)
                    gas_data = np.zeros((spectra_number, len(kin_id_gas) - 1))
                    df_kin_gas = pd.DataFrame(np.column_stack((params.spec_names_nopath, gas_data)), columns=kin_id_gas)


                kin_component = np.array(kin_component)

                for t, comp in enumerate(kin_component[gas_component]):
                    name = gas_names[t]
                    df_kin_gas.at[i, f'RV(km/s)_{name}']= round(kinematics[comp][0],3)
                    df_kin_gas.at[i, f'Sigma(km/s)_{name}']= round(kinematics[comp][1],3)
                    df_kin_gas.at[i, f'H3_{name}']= round(kinematics[comp][2],5)
                    df_kin_gas.at[i, f'H4_{name}']= round(kinematics[comp][3],5)
                    df_kin_gas.at[i, f'H5_{name}']= round(kinematics[comp][4],5)
                    df_kin_gas.at[i, f'H6_{name}']= round(kinematics[comp][5],5)
                    #gas flux
                    df_kin_gas.at[i, f'Flux_{name}']= round(gas_flux[t],5)
                    df_kin_gas.at[i, f'Flux_err_{name}']= round(gas_flux_err[t],5)
                    df_kin_gas.at[i, f'errRV_{name}']= round(error_kinematics[comp][0],3)
                    df_kin_gas.at[i, f'errSigma_{name}']= round(error_kinematics[comp][1],3)
                    df_kin_gas.at[i, f'errH3_{name}']= round(error_kinematics[comp][2],5)
                    df_kin_gas.at[i, f'errH4_{name}']= round(error_kinematics[comp][3],5)
                    df_kin_gas.at[i, f'errH5_{name}']= round(error_kinematics[comp][4],5)
                    df_kin_gas.at[i, f'errH6_{name}']= round(error_kinematics[comp][5],5)

                    # df_kin_gas.to_csv(kin_file_gas, index= False, sep=' ')
                if df_kin_gas is not None and kin_file_gas is not None:
                    df_kin_gas.to_csv(kin_file_gas, index=False, sep=' ')

            if params.with_errors_kin:
                err_rv_kin_mc, err_sigma_kin_mc, err_h3_kin_mc, err_h4_kin_mc, err_h5_kin_mc, err_h6_kin_mc = np.round(error_kinematics_mc[0],5)

                df_kin_mc.at[i, 'RV(km/s)']= vel
                df_kin_mc.at[i, 'Sigma(km/s)']= sigma
                df_kin_mc.at[i, 'H3']= h3
                df_kin_mc.at[i, 'H4']= h4
                df_kin_mc.at[i, 'H5']= h5
                df_kin_mc.at[i, 'H6']= h6
                df_kin_mc.at[i, 'errRV']= err_rv_kin_mc
                df_kin_mc.at[i, 'errSigma']= err_sigma_kin_mc
                df_kin_mc.at[i, 'errH3']= err_h3_kin_mc
                df_kin_mc.at[i, 'errH4']= err_h4_kin_mc
                df_kin_mc.at[i, 'errH5']= err_h5_kin_mc
                df_kin_mc.at[i, 'errH6']= err_h6_kin_mc
                df_kin_mc.at[i, 'S/N']= int(snr_kin)

                df_kin_mc.to_csv(kin_file_mc, index= False, sep=' ')
        #print message
        if i == (params.spectra_number_to_process-1):
            print ('File with stellar kinematics saved: ', kin_file)
            if params.gas_kin:
                print ('File with gas kinematics saved: ', kin_file_gas)
            if params.with_errors_kin:
                print ('File with stellar kinematics and MonteCarlo uncertainties saved: ', kin_file_mc)
            print('')

        return df_kin_gas
    except Exception:
        print('Cannot write the files')
        return None



def save_population_analysis_to_file(i, params, kinematics, info_pop, info_pop_mass, mass_light,
                                     chi_square, met_err, mass_met_err, snr_pop, age_err_abs,
                                     mass_age_err_abs, alpha_err, mass_alpha_err, t50_age, t80_age, t50_cosmic, t80_cosmic, ssp_lick_indices_ppxf,
                                     ssp_lick_indices_err_ppxf, ppxf_lick_params, df_pop, pop_file,
                                     df_ssp_param_ppxf, ssp_param_file_ppxf):

    """
    Saves stellar population analysis results to a file

    """

    try:
        # Extracting kinematic values
        try:
            num_comp_kinematics = len(kinematics)
            kin_stars = np.array(kinematics[0])
            rv_pop_ppxf, sigma_pop_ppxf, h3_pop_ppxf, h4_pop_ppxf = kin_stars[:4]
        except (ValueError, IndexError, TypeError):
            num_comp_kinematics = 0
            rv_pop_ppxf, sigma_pop_ppxf, h3_pop_ppxf, h4_pop_ppxf = kinematics[:4]

        # Extracting population parameters
        age, met = info_pop[:2]
        mass_age, mass_met = info_pop_mass[:2]

        # Handling alpha enhancement if using sMILES
        alpha, mass_alpha = (None, None)

        # storing the values in the dataframes and save to disc
        if params.stellar_library == 'sMILES' and not params.ppxf_pop_custom_lib:
            alpha = info_pop[2]
            mass_alpha = info_pop_mass[2]

        df_pop.at[i, 'RV(km/s)']= round(rv_pop_ppxf,2)
        df_pop.at[i, 'Sigma(km/s)']= round(sigma_pop_ppxf,2)
        df_pop.at[i, 'H3']= round(h3_pop_ppxf,3)
        df_pop.at[i, 'H4']= round(h4_pop_ppxf,3)
        df_pop.at[i, 'lum_met(dex)']= round(met,2)
        df_pop.at[i, 'M/L']= round(mass_light,3)
        df_pop.at[i, 'mass_met(dex)']= round(mass_met,3)
        df_pop.at[i, 'Chi2']= round(chi_square,3)
        df_pop.at[i, 'err_lum_met(dex)']= round(met_err,3)
        df_pop.at[i, 'err_mass_met(dex)']= round(mass_met_err,3)
        df_pop.at[i, 't50_age']= round(t50_age,2)
        df_pop.at[i, 't80_age']= round(t80_age,2)
        df_pop.at[i, 't50_cosmic']= round(t50_cosmic,2)
        df_pop.at[i, 't80_cosmic']= round(t80_cosmic,2)
        df_pop.at[i, 'S/N']= round(snr_pop)

        if params.ppxf_pop_lg_age:
            df_pop.at[i, 'lum_lg_age(dex)']= round(age,2)
            df_pop.at[i, 'mass_lg_age(dex)']= round(mass_age,2)
            df_pop.at[i, 'err_lum_lg_age(dex)']= round(age_err_abs,2)
            df_pop.at[i, 'err_mass_lg_age(dex)']= round(mass_age_err_abs,2)
        else:
            df_pop.at[i, 'lum_age(Gyr)']= round(age,1)
            df_pop.at[i, 'mass_age(Gyr)']= round(mass_age,1)
            df_pop.at[i, 'err_lum_age(Gyr)']= round(age_err_abs,2)
            df_pop.at[i, 'err_mass_age(Gyr)']= round(mass_age_err_abs,2)

        #In case I use the sMILES with alpha/Fe
        if params.stellar_library == 'sMILES' and not params.ppxf_pop_custom_lib:
            df_pop.at[i, 'lum_alpha(dex)']= round(alpha,2)
            df_pop.at[i, 'mass_alpha(dex)']= round(mass_alpha,2)
            df_pop.at[i, 'err_lum_alpha(dex)']= round(alpha_err,2)
            df_pop.at[i, 'err_mass_alpha(dex)']= round(mass_alpha_err,2)

        #storing to the file
        df_pop.to_csv(pop_file, index= False, sep=' ')

        # # If I want also to measure stellar parameters with Lick/IDS indices
        if params.stellar_parameters_lick_ppxf:

            #Storing the results to a file
            df_ssp_param_ppxf.at[i, 'Hbeta(A)']= round(ssp_lick_indices_ppxf[0],3)
            df_ssp_param_ppxf.at[i, 'Hbeta_err(A)']= round(ssp_lick_indices_err_ppxf[0],3)
            # df_ssp_param_ppxf.at[i, 'Mg2(mag)']= round(Mg2_ppxf,3)
            # df_ssp_param_ppxf.at[i, 'Mg2_err(mag)']= round(Mg2e_ppxf,3)
            df_ssp_param_ppxf.at[i, 'Mgb(A)']= round(ssp_lick_indices_ppxf[3],3)
            df_ssp_param_ppxf.at[i, 'Mgb_err(A)']= round(ssp_lick_indices_err_ppxf[3],3)
            df_ssp_param_ppxf.at[i, 'Fem(A)']= round(ssp_lick_indices_ppxf[2],3)
            df_ssp_param_ppxf.at[i, 'Fem_err(A)']= round(ssp_lick_indices_err_ppxf[2],3)
            df_ssp_param_ppxf.at[i, 'MgFe(A)']= round(ssp_lick_indices_ppxf[1],3)
            df_ssp_param_ppxf.at[i, 'MgFe_err(A)']= round(ssp_lick_indices_err_ppxf[1],3)
            df_ssp_param_ppxf.at[i, 'age(Gyr)']= round(ppxf_lick_params[0],3)
            df_ssp_param_ppxf.at[i, 'err_age']= round(ppxf_lick_params[3],3)
            df_ssp_param_ppxf.at[i, 'met']= round(ppxf_lick_params[1],4)
            df_ssp_param_ppxf.at[i, 'err_met']= round(ppxf_lick_params[4],4)
            df_ssp_param_ppxf.at[i, 'alpha']= round(ppxf_lick_params[2],4)
            df_ssp_param_ppxf.at[i, 'err_alpha']= round(ppxf_lick_params[5],4)

            #putting nans where needed
            df_ssp_param_ppxf.to_csv(ssp_param_file_ppxf, na_rep='NaN', index= False, sep=' ')

        # at the last spectrum I print some info on the output window
        if i == (params.spectra_number_to_process-1):
            print ('File saved: ', pop_file)
            if params.stellar_parameters_lick_ppxf:
                print ('File with the Lick/IDS stellar parameters saved: ', ssp_param_file_ppxf)
            print('')
    except Exception:
        print('Cannot write the files')



def save_ew_to_file(i, params, ew, err, ew_mag, err_mag, df_ew, ew_file,
                    df_ew_mag, ew_file_mag, df_snr_ew, snr_ew_file, snr_ew):

    """
    Saves Equivalent Width (EW), EW in magnitudes, and Signal-to-Noise Ratio (SNR) to files.

    """

    try:
        #Updating and writing the file
        print ('EW:', ew, '+/-', err)
        print ('EW Mag', ew_mag, '+/-', err_mag)
        print ('SNR: ', snr_ew, 'per pix')
        print ('')

        df_ew.at[i, 'ew(A)']= round(ew,4)
        df_ew.at[i, 'err']= round(err,4)
        df_ew.to_csv(ew_file, index= False, sep=' ')

        df_ew_mag.at[i, 'ew(Mag)']= round(ew_mag,4)
        df_ew_mag.at[i, 'err']= round(err_mag,4)
        df_ew_mag.to_csv(ew_file_mag, index= False, sep=' ')

        df_snr_ew.at[i, 'SNR']= round(snr_ew,4)
        df_snr_ew.to_csv(snr_ew_file, index= False, sep=' ')

        if i == (params.spectra_number_to_process-1):
            print ('File EW saved: ', ew_file)
            print ('File EW in Mag saved: ', ew_file_mag)
            print ('File SNR saved: ', snr_ew_file)
            print('')
    except Exception:
        print('Error writing the file')
        print('')




def save_ew_indices_to_file(i, params, num_indices, ew_array, err_array, ew_array_mag, err_array_mag,
                            snr_ew_array, df_ew, ew_file, df_ew_mag, ew_file_mag,
                            df_snr_ew, snr_ew_file, ew_id, ew_id_mag, snr_ew_id, spectra_id):

    """
    Saves Equivalent Width (EW), EW in magnitudes, and Signal-to-Noise Ratio (SNR) for multiple indices.

    """

    try:
        #Updating and writing the file
        for k in range(num_indices):
            df_ew.at[i,ew_id[k+len(spectra_id)]]= round(ew_array[k], 4)
            df_ew.at[i,ew_id[k+num_indices+ len(spectra_id)]] = round(err_array[k],4)
            df_ew.to_csv(ew_file, index= False, sep=' ')

            df_ew_mag.at[i,ew_id_mag[k+len(spectra_id)]]= round(ew_array_mag[k], 4)
            df_ew_mag.at[i,ew_id_mag[k+num_indices+ len(spectra_id)]] = round(err_array_mag[k],4)
            df_ew_mag.to_csv(ew_file_mag, index= False, sep=' ')

            df_snr_ew.at[i,snr_ew_id[k+len(spectra_id)]]= round(snr_ew_array[k], 4)
            df_snr_ew.to_csv(snr_ew_file, index= False, sep=' ')
        if i == (params.spectra_number_to_process-1):
            print ('File EW saved: ', ew_file)
            print ('File EW in Mag saved: ', ew_file_mag)
            print ('File SNR saved: ', snr_ew_file)
            print('')
    except Exception:
        print('Error writing the file')
        print('')




def save_lick_indices_to_file(i, params, num_lick_indices, lick_ew_array, lick_err_array, lick_ew_array_mag,
                              lick_err_array_mag, lick_snr_ew_array, df_ew_lick, ew_lick_file, df_ew_lick_mag,
                              ew_lick_file_mag, df_snr_lick_ew, snr_lick_ew_file, ew_lick_id, ew_lick_id_mag,
                              snr_lick_ew_id, spectra_lick_id, df_lick_param, ssp_lick_param_file, lick_for_ssp,
                              df_ssp_param, ssp_param_file, age, err_age, met, err_met, alpha, err_alpha, save_plot,
                              ssp_lick_indices_list, ssp_lick_indices_err_list, spectra_list_name, result_plot_dir,
                              ssp_model, lick_to_plot, lick_err_to_plot):

    """
    Saves Lick/IDS indices, equivalent width (EW), EW in magnitudes, and signal-to-noise ratio (SNR) to files.

    """

    try:

        #Updating and writing the file
        for k in range(num_lick_indices):
            df_ew_lick.at[i,ew_lick_id[k+len(spectra_lick_id)]]= round(lick_ew_array[k], 4)
            df_ew_lick.at[i,ew_lick_id[k+num_lick_indices+ len(spectra_lick_id)]] = round(lick_err_array[k],4)
            df_ew_lick.to_csv(ew_lick_file, index= False, sep=' ')
            df_ew_lick_mag.at[i,ew_lick_id_mag[k+len(spectra_lick_id)]]= round(lick_ew_array_mag[k], 4)
            df_ew_lick_mag.at[i,ew_lick_id_mag[k+num_lick_indices+ len(spectra_lick_id)]] = round(lick_err_array_mag[k],4)
            df_ew_lick_mag.to_csv(ew_lick_file_mag, index= False, sep=' ')
            df_snr_lick_ew.at[i,snr_lick_ew_id[k+len(spectra_lick_id)]]= round(lick_snr_ew_array[k], 4)
            df_snr_lick_ew.to_csv(snr_lick_ew_file, index= False, sep=' ')

        #Storing the results to a file
        df_lick_param.at[i, 'Hbeta(A)']= round(lick_for_ssp[0],3)
        df_lick_param.at[i, 'Hbeta_err(A)']= round(lick_for_ssp[1],3)
        df_lick_param.at[i, 'Mg2(mag)']= round(lick_for_ssp[2],3)
        df_lick_param.at[i, 'Mg2_err(mag)']= round(lick_for_ssp[3],3)
        df_lick_param.at[i, 'Mgb(A)']= round(lick_for_ssp[4],3)
        df_lick_param.at[i, 'Mgb_err(A)']= round(lick_for_ssp[5],3)
        df_lick_param.at[i, 'Fe5270(A)']= round(lick_for_ssp[6],3)
        df_lick_param.at[i, 'Fe5270_err(A)']= round(lick_for_ssp[7],3)
        df_lick_param.at[i, 'Fe5335(A)']= round(lick_for_ssp[8],3)
        df_lick_param.at[i, 'Fe5335_err(A)']= round(lick_for_ssp[9],3)
        df_lick_param.at[i, 'Fem(A)']= round(lick_for_ssp[10],3)
        df_lick_param.at[i, 'Fem_err(A)']= round(lick_for_ssp[11],3)
        df_lick_param.at[i, 'MgFe(A)']= round(lick_for_ssp[12],3)
        df_lick_param.at[i, 'MgFe_err(A)']= round(lick_for_ssp[13],3)

        df_lick_param.to_csv(ssp_lick_param_file, na_rep='NaN', index= False, sep=' ')

        if params.stellar_parameters_lick:
            #printing to file
            df_ssp_param.at[i, 'age(Gyr)']= round(age,3)
            df_ssp_param.at[i, 'err_age']= round(err_age,3)
            df_ssp_param.at[i, 'met']= round(met,4)
            df_ssp_param.at[i, 'err_met']= round(err_met,4)
            df_ssp_param.at[i, 'alpha']= round(alpha,4)
            df_ssp_param.at[i, 'err_alpha']= round(err_alpha,4)

            df_ssp_param.to_csv(ssp_param_file, na_rep='NaN', index= False, sep=' ')

            #doing plot pf the index-index-grid
            if save_plot:
                lick_to_plot.append(ssp_lick_indices_list)
                lick_err_to_plot.append(ssp_lick_indices_err_list)

                if i == (params.spectra_number_to_process - 1):
                    lick_to_plot_np = np.vstack(lick_to_plot)
                    lick_err_to_plot_np = np.vstack(lick_err_to_plot)
                    span.lick_grids(ssp_model, lick_to_plot_np, lick_err_to_plot_np, age, False, True, spectra_list_name, result_plot_dir)


        if i == (params.spectra_number_to_process-1):
            print ('File EW saved: ', ew_lick_file)
            print ('File EW in Mag saved: ', ew_lick_file_mag)
            print ('File SNR saved: ', snr_lick_ew_file)
            if params.stellar_parameters_lick:
                print ('File with the stellar parameters saved: ', ssp_param_file)

    except Exception:
        print('Cannot write the files')
        print('')



def save_velocity_or_redshift_to_file(i, params, value_at_max, err, df_rv, rv_file):

    """
    Saves velocity (RV) or redshift (z) to file based on cross-correlation method.

    """

    #Updating and writing the file for velocity
    if params.is_vel_xcorr:
        try:
            df_rv.at[i, 'RV(km/s)']= round(value_at_max,1)
            df_rv.at[i, 'err']= round(err,1)
            df_rv.to_csv(rv_file, index= False, sep=' ')

            if i == (params.spectra_number_to_process-1):
                print ('File saved: ', rv_file)
                print('')
        except Exception:
            print ('Error saving the file')
            print('')

    #Updating and writing the file for z
    if not params.is_vel_xcorr:
        try:
            #Changing che name of the column RV, now is z!
            df_rv.rename(columns={'RV(km/s)': 'z'}, inplace=True)
            #filling the values
            df_rv.at[i, 'z'] = round(value_at_max, 5)
            df_rv.at[i, 'err'] = round(err, 5)
            df_rv.to_csv(rv_file, index= False, sep=' ')

            if i == (params.spectra_number_to_process-1):
                print ('File saved: ', rv_file)
                print('')
        except Exception:
            print ('Error writing the file')
            print('')
