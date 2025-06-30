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
#************** FUNCTIONS TO SELECT SSP TEMPLATES TO BE PASSED TO PPXF ********************
#******************************************************************************************
#******************************************************************************************

import glob
import re
import numpy as np
from scipy import ndimage
from astropy.io import fits
import ppxf.ppxf_util as util
import ppxf.sps_util as lib


###### The following classes and functions have been written to be an interface with the pPXF code.
# They allow to interact with the pre-loaded SSP templates that come with pPXF and to use
# the Xshooter and any (E)MILES based template set with SPAN.


###############################################################################
# Wrapper to use and plot the preloaded models with pPXF
###############################################################################
class SPSLibWrapper:
    def __init__(self, filename, velscale, fwhm_gal=None, age_range=None, lam_range=None,
                 metal_range=None, norm_range=None, norm_type='mean'):
        """
        Wrapper around pPXF's sps_lib, to maintain a consistent interface and
        facilitate quick usage and plotting of preloaded models.
        """
        self.sps_instance = lib.sps_lib(
            filename, velscale, fwhm_gal=fwhm_gal, age_range=age_range,
            lam_range=lam_range, metal_range=metal_range,
            norm_range=norm_range, norm_type=norm_type
        )

    def get_age_grid_2d(self):
        self.age_grid = self.sps_instance.age_grid
        return self.age_grid

    def get_metal_grid_2d(self):
        self.metal_grid = self.sps_instance.metal_grid
        return self.metal_grid

    def mean_age_metal(self, weights, lg_age=True, lg_met = True, quiet=False):
        """
        Computes the weighted mean age and metallicity from pPXF weights.
        :param weights: pPXF weight array
        :param lg_age: If True, compute log10(Age/yr), else linear Gyr
        :param quiet: If True, suppress printout
        :return: (mean_age, mean_metal)
        """
        self.age_grid = self.get_age_grid_2d()
        self.metal_grid = self.get_metal_grid_2d()

        if lg_age:
            lg_age_grid = np.log10(self.age_grid) + 9
            mean_lg_age = np.sum(weights * lg_age_grid) / np.sum(weights)
            mean_age = mean_lg_age
        else:
            lin_age_grid = self.age_grid
            mean_lin_age = np.sum(weights * lin_age_grid) / np.sum(weights)
            mean_age = mean_lin_age

        if lg_met:
            mean_metal = np.sum(weights * self.metal_grid) / np.sum(weights)
        else:
            lin_met_grid = 10**self.metal_grid
            mean_metal = np.sum(weights * lin_met_grid) / np.sum(weights)
            mean_metal = np.log10(mean_metal)




        if not quiet:
            if lg_age:
                print(f'Weighted lg <Age>: {mean_age:.3g}')
                print(f'Weighted <[M/H]>: {mean_metal:.3g}')
            else:
                print(f'Weighted <Age> [Gyr]: {mean_age:.3g}')
                print(f'Weighted <[M/H]>: {mean_metal:.3g}')

        return mean_age, mean_metal

    def plot(self, weights, lg_age=True, nodots=False, colorbar=True, **kwargs):
        """
        Plot the weights in 2D vs age and metallicity.
        :param weights: 2D array of pPXF weights
        :param lg_age: If True, plot vs log10(Age/yr); else linear Age in Gyr
        :param nodots: If True, do not plot the white dots on the map
        :param colorbar: If True, include a colorbar
        """
        ygrid = self.sps_instance.metal_grid
        if lg_age:
            xgrid = np.log10(self.sps_instance.age_grid) + 9
            xlabel = "lg Age (dex)"
        else:
            xgrid = self.sps_instance.age_grid
            xlabel = "Age (Gyr)"

        util.plot_weights_2d(
            xgrid, ygrid, weights,
            xlabel=xlabel, nodots=nodots, colorbar=colorbar, **kwargs
        )

    def get_age_grid(self):
        """
        Return the 1D age array (assuming same ages across metallicity axis).
        """
        return self.sps_instance.age_grid[:, 0]

    def get_metal_grid(self):
        """
        Return the 1D metallicity array (assuming same metals across age axis).
        """
        return self.sps_instance.metal_grid[0, :]


###############################################################################
# For retrieving age and metallicity from the (E)MILES standard
def age_metal_miles(filename):
    """
    Extract age and metallicity from (E)MILES template filename.
    Expected pattern: 'Z[m|p][0-9].[0-9]{2}T[0-9]{2}.[0-9]{4}'
    Example: 'Zp0.06T10.0000', 'Zm0.02T05.0000', etc.
    """
    match = re.search(r'Z[m|p][0-9]\.[0-9]{2}T[0-9]{2}\.[0-9]{4}', filename)
    if not match:
        raise ValueError(f"File format not valid: {filename}")

    s = match.group(0)
    metal = s[:6]  # e.g. 'Zp0.06' or 'Zm0.02'
    age = float(s[7:])  # e.g. '10.0000' -> 10.0

    if "Zm" in metal:
        metal = -float(metal[2:])
    elif "Zp" in metal:
        metal = float(metal[2:])

    return age, metal


###############################################################################
# For retrieving age and metallicity from the filenames of the XSHOOTER SSP library
def age_metal_xshooter(filename):
    """
    Extract age and metallicity from XSHOOTER template filenames.
    Example patterns:
     - 'logT9.20_MH-0.34.fits' -> age ~ 1.58 Gyr, metal ~ -0.34
    """
    age_match = re.search(r'logT(\d+\.\d+)', filename)
    if not age_match:
        raise ValueError(f"Age not found in filename: {filename}")
    # Convert log(Age/yr) -> Age in Gyr
    age = 10 ** float(age_match.group(1)) / 1e9

    metal_match = re.search(r'MH(-?\d+\.\d+)', filename)
    if not metal_match:
        raise ValueError(f"Metallicity not found in filename: {filename}")
    metal = float(metal_match.group(1))

    return age, metal


###############################################################################
# For retrieving age, metallicity, and alpha from the sMILES standard
def age_metal_smiles(filename):
    """
    Extract age, metallicity, and alpha-enhancement from sMILES template filenames.
    Example pattern for age/metal: 'Zp0.06T10.0000'
    Example pattern for alpha: 'aFep10' or 'aFem05'
    """
    # Extract age/metal part
    s = re.findall(r'Z[m|p][0-9]\.[0-9]{2}T[0-9]{2}\.[0-9]{4}', filename)[0]
    metal = s[:6]  # e.g. 'Zp0.06' or 'Zm0.02'
    age = float(s[7:])  # e.g. '10.0000' -> 10.0
    if "Zm" in metal:
        metal = -float(metal[2:])
    elif "Zp" in metal:
        metal = float(metal[2:])

    # Extract alpha part
    aFe_match = re.findall(r'aFe[m|p][0-9]{2}', filename)
    if not aFe_match:
        raise ValueError(f"Alpha enhancement not found in filename: {filename}")

    aFe_find = aFe_match[0]  # e.g. 'aFep10' or 'aFem05'
    if "aFem" in aFe_find:
        afe = -float(aFe_find[4:])
    elif "aFep" in aFe_find:
        afe = float(aFe_find[4:])
    afe = afe / 10.0

    return age, metal, afe


###############################################################################
# Wrapper function for plotting the age-metal weights for all SSP models
def plot_weights(weights, age_grid, metal_grid, lg_age=True, nodots=False, colorbar=True, **kwargs):
    """
    Plot a 2D map of weights vs age and metallicity.
    If weights have dimension=3, they contain alpha as the third axis
    and will be summed over alpha.
    """
    if weights.ndim == 3:
        # Sum over alpha dimension if present
        weights = np.sum(weights, axis=2)
        age_grid = np.mean(age_grid, axis=2)
        metal_grid = np.mean(metal_grid, axis=2)

    ygrid = metal_grid

    if lg_age:
        # Convert age to log10(Age/yr)
        xgrid = np.log10(age_grid) + 9
        xlabel = "lg Age (dex)"
    else:
        xgrid = age_grid
        xlabel = "Age (Gyr)"

    util.plot_weights_2d(
        xgrid, ygrid, weights, xlabel=xlabel,
        nodots=nodots, colorbar=colorbar, **kwargs
    )


###############################################################################
# Wrapper for plotting the alpha weights for the sMILES models
def plot_alpha_weights(weights, alpha_grid, metal_grid, lg_age=True, nodots=False, colorbar=True, **kwargs):
    """
    Plots a 2D map of the summed weights over age dimension, as a function
    of metallicity and alpha-enhancement.
    """
    # Sum over the age dimension
    reduced_weights_alpha = np.sum(weights, axis=0)

    # Compute mean alpha and metallicity grids over age
    plot_alpha_grid = np.mean(alpha_grid, axis=0)
    plot_alpha_grid_metal = np.mean(metal_grid, axis=0)

    ygrid = plot_alpha_grid_metal
    zgrid = plot_alpha_grid

    util.plot_weights_2d(
        ygrid, zgrid, reduced_weights_alpha,
        xlabel="[M/H]", ylabel=r"[$\alpha$/Fe]",
        nodots=nodots, colorbar=colorbar, **kwargs
    )


###############################################################################
# Class for the (E)MILES models (no alpha dimension)
###############################################################################
class miles:
    def __init__(self, pathname, velscale, FWHM_gal=None, FWHM_tem=2.51,
                 age_range=None, metal_range=None, norm_range=None, wave_range=None, R = None):
        """
        Class handling (E)MILES SSP templates without alpha dimension.
        """
        # Find all files matching the provided pattern
        files = glob.glob(pathname)
        if not files:
            raise FileNotFoundError(f"No files found matching pattern: {pathname}")

        # Extract age and metallicity for each file
        age_metal_list = [age_metal_miles(f) for f in files]
        all_ages = np.array([am[0] for am in age_metal_list])
        all_metals = np.array([am[1] for am in age_metal_list])

        # Identify unique ages and metals
        ages = np.unique(all_ages)
        metals = np.unique(all_metals)
        n_ages = len(ages)
        n_metal = len(metals)

        # Create a dictionary for direct access: (age, metal) -> filename
        # This avoids repeated list.index() calls which can be slow
        file_map = {}
        for i, f in enumerate(files):
            key = (all_ages[i], all_metals[i])
            file_map[key] = f

        # Read one file to get initial wavelength info
        with fits.open(files[0]) as hdu:
            ssp = hdu[0].data
            header = hdu[0].header
        lam = header['CRVAL1'] + np.arange(header['NAXIS1']) * header['CDELT1']

        # Cutting the templates to +- 0.02% of the wave_range inserted
        wave_min, wave_max = wave_range
        delta = 0.02 * (wave_max - wave_min)
        cut_min = wave_min - delta
        cut_max = wave_max + delta
        mask = (lam >= cut_min) & (lam <= cut_max)
        lam_cut = lam[mask]
        ssp_cut = ssp[..., mask]

        lam = lam_cut
        ssp = ssp_cut

        lam_range_temp = lam[[0, -1]]

        # Log-rebin to the chosen velocity scale
        ssp_new, ln_lam_temp = util.log_rebin(lam_range_temp, ssp, velscale=velscale)[:2]
        lam_temp = np.exp(ln_lam_temp)

        # Define normalization band if needed
        if norm_range is not None:
            band = (norm_range[0] <= lam_temp) & (lam_temp <= norm_range[1])

        # Prepare arrays
        templates = np.empty((ssp_new.size, n_ages, n_metal))
        age_grid = np.empty((n_ages, n_metal))
        metal_grid = np.empty((n_ages, n_metal))
        flux = np.empty((n_ages, n_metal))


        # Considering three possibilities: Resolution in FWHM, in R and in MUSE LSF
        if FWHM_gal is not None:
            if isinstance(FWHM_gal, (int, float)):
                FWHM_dif = np.sqrt(FWHM_gal**2 - FWHM_tem**2)
                sigma = FWHM_dif/2.355/header['CDELT1']   # Sigma difference in pixels
            elif R is not None:
                FWHM_gal = lam/R #FWHM_gal[mask]
                FWHM_dif = np.sqrt(FWHM_gal**2 - FWHM_tem**2)
                sigma = FWHM_dif/2.355/header['CDELT1']   # Sigma difference in pixels
            else:
                FWHM_gal = 5.866e-8*lam**2-9.187e-4*lam+6.040
                FWHM_dif = np.sqrt(FWHM_gal**2 - FWHM_tem**2)
                sigma = FWHM_dif/2.355/header['CDELT1']   # Sigma difference in pixels

        # Sort and process templates by age and metallicity
        warn_conv_failed = False
        for j, age in enumerate(ages):
            for k, met in enumerate(metals):
                # Retrieve the correct file via our dictionary
                fname = file_map.get((age, met), None)
                if fname is None:
                    raise ValueError(f"Missing template for age={age}, metallicity={met}")
                with fits.open(fname) as hdu:
                    ssp_data = hdu[0].data

                #cutting the templates
                ssp_data = ssp_data[..., mask]

                # Convolving
                if FWHM_gal is not None:
                    if np.isscalar(FWHM_gal):
                        if sigma > 0.01:   # Skip convolution for nearly zero sigma
                            x = np.arange(len(ssp_data))
                            ssp_data = util.varsmooth(x, ssp_data, sigma) # convolution with variable sigma
                        else:
                            if not warn_conv_failed:
                                print('WARNING: The resolution of the templates is lower than galaxy. Skipping convolution')
                                warn_conv_failed = True
                    else:
                        x = np.arange(len(ssp_data))
                        try:
                            ssp_data = util.varsmooth(x, ssp_data, sigma) # convolution with variable sigma
                        except ValueError:
                            if not warn_conv_failed:
                                print('WARNING: The resolution of the templates is lower than galaxy. Skipping convolution')
                                warn_conv_failed = True

                # Log-rebin the template
                ssp_new = util.log_rebin(lam_range_temp, ssp_data, velscale=velscale)[0]

                # Normalize if norm_range is specified
                if norm_range is not None:
                    flux_val = np.mean(ssp_new[band])
                    flux[j, k] = flux_val
                    ssp_new /= flux_val
                else:
                    flux[j, k] = 1.0  # Just for bookkeeping

                # Store the processed template and grids
                templates[:, j, k] = ssp_new
                age_grid[j, k] = age
                metal_grid[j, k] = met

        # Apply optional age range mask
        if age_range is not None:
            age_mask = (age_range[0] <= ages) & (ages <= age_range[1])
            templates = templates[:, age_mask, :]
            age_grid = age_grid[age_mask, :]
            metal_grid = metal_grid[age_mask, :]
            flux = flux[age_mask, :]
            ages = ages[age_mask]

        # Apply optional metal range mask
        if metal_range is not None:
            metal_mask = (metal_range[0] <= metals) & (metals <= metal_range[1])
            templates = templates[:, :, metal_mask]
            age_grid = age_grid[:, metal_mask]
            metal_grid = metal_grid[:, metal_mask]
            flux = flux[:, metal_mask]
            metals = metals[metal_mask]

        # If no normalization range is given, do a global normalization
        if norm_range is None:
            flux_median = np.median(templates[templates > 0])
            templates /= flux_median

        # Save full arrays (before cutting wave range)
        self.templates_full = templates
        self.ln_lam_temp_full = ln_lam_temp
        self.lam_temp_full = lam_temp

        # Apply optional wavelength range
        if wave_range is not None:
            wave_mask = (wave_range[0] <= lam_temp) & (lam_temp <= wave_range[1])
            ln_lam_temp = ln_lam_temp[wave_mask]
            lam_temp = lam_temp[wave_mask]
            templates = templates[wave_mask, :, :]

        # Final attributes
        self.templates = templates
        self.ln_lam_temp = ln_lam_temp
        self.lam_temp = lam_temp
        self.age_grid = age_grid
        self.metal_grid = metal_grid
        self.n_ages, self.n_metal = age_grid.shape
        self.flux = flux



    ###############################################################################
    def mean_age_metal(self, weights, lg_age=True, lg_met = True, quiet=False):
        """
        Calculate mean age and metallicity from pPXF weights.
        :param weights: pPXF output with dimensions [n_ages, n_metal]
        :param lg_age: If True, compute mean log10(Age/yr), else linear mean (Gyr).
        :param quiet: If True, suppress printed output.
        :return: (mean_age, mean_metal)
        """
        if lg_age:
            # Convert age to log10(Age/yr)
            lg_age_grid = np.log10(self.age_grid) + 9
            mean_lg_age = np.sum(weights * lg_age_grid) / np.sum(weights)
            mean_age = mean_lg_age
        else:
            # Linear mean age in Gyr
            lin_age_grid = self.age_grid
            mean_lin_age = np.sum(weights * lin_age_grid) / np.sum(weights)
            mean_age = mean_lin_age

        if lg_met:
            mean_metal = np.sum(weights * self.metal_grid) / np.sum(weights)
        else:
            lin_met_grid = 10**self.metal_grid
            mean_metal = np.sum(weights * lin_met_grid) / np.sum(weights)
            mean_metal = np.log10(mean_metal)


        if not quiet:
            if lg_age:
                print(f'Weighted lg <Age>: {mean_age:.3g}')
                print(f'Weighted <[M/H]>: {mean_metal:.3g}')
            else:
                print(f'Weighted <Age> [Gyr]: {mean_age:.3g}')
                print(f'Weighted <[M/H]>: {mean_metal:.3g}')

        return mean_age, mean_metal

    ###############################################################################
    def get_full_age_grid(self):
        """
        Return the full 2D array of ages.
        """
        return self.age_grid

    ###############################################################################
    def get_full_metal_grid(self):
        """
        Return the full 2D array of metallicities.
        """
        return self.metal_grid


###############################################################################
# Class for sMILES templates WITH alpha dimension
###############################################################################
class smiles:
    def __init__(self, pathname, velscale, FWHM_gal=None, FWHM_tem=2.51,
                 age_range=None, metal_range=None, afe_range=None, norm_range=None, wave_range=None):
        """
        Class handling sMILES SSP templates with alpha dimension.
        """
        files = glob.glob(pathname)
        if not files:
            raise FileNotFoundError(f"No files found matching pattern: {pathname}")

        # Extract age, metallicity, alpha for each file
        ama_list = [age_metal_smiles(f) for f in files]
        all_ages = np.array([x[0] for x in ama_list])
        all_metals = np.array([x[1] for x in ama_list])
        all_alpha = np.array([x[2] for x in ama_list])

        ages = np.unique(all_ages)
        metals = np.unique(all_metals)
        alphas = np.unique(all_alpha)
        n_ages = len(ages)
        n_metal = len(metals)
        n_alpha = len(alphas)

        # Dictionary for file lookup
        file_map = {}
        for i, f in enumerate(files):
            file_map[(all_ages[i], all_metals[i], all_alpha[i])] = f

        # Read the first file for wavelength info
        with fits.open(files[0]) as hdu:
            ssp = hdu[0].data
            header = hdu[0].header
        lam = header['CRVAL1'] + np.arange(header['NAXIS1']) * header['CDELT1']
        lam_range_temp = lam[[0, -1]]

        # Log-rebin
        ssp_new, ln_lam_temp = util.log_rebin(lam_range_temp, ssp, velscale=velscale)[:2]
        lam_temp = np.exp(ln_lam_temp)

        # Normalization band if specified
        if norm_range is not None:
            band = (norm_range[0] <= lam_temp) & (lam_temp <= norm_range[1])

        # Prepare output arrays
        templates = np.empty((ssp_new.size, n_ages, n_metal, n_alpha))
        age_grid = np.empty((n_ages, n_metal, n_alpha))
        metal_grid = np.empty((n_ages, n_metal, n_alpha))
        alpha_grid = np.empty((n_ages, n_metal, n_alpha))
        flux = np.empty((n_ages, n_metal, n_alpha))

        # Compute sigma for convolution if needed
        if FWHM_gal is not None:
            # Avoid sqrt of negative if FWHM_gal < FWHM_tem
            if FWHM_gal**2 > FWHM_tem**2:
                FWHM_dif = np.sqrt(FWHM_gal**2 - FWHM_tem**2)
                sigma = FWHM_dif / 2.355 / header['CDELT1']  # Sigma in pixels
            else:
                print('WARNING: resolution of the templates is lower than the galaxy. Skipping convolution')
                sigma = 0.0


        # Process each template
        for j, age in enumerate(ages):
            for k, met in enumerate(metals):
                for h, afe in enumerate(alphas):
                    fname = file_map.get((age, met, afe), None)
                    if fname is None:
                        raise ValueError(f"Missing template for age={age}, [M/H]={met}, [α/Fe]={afe}")

                    with fits.open(fname) as hdu:
                        ssp_data = hdu[0].data

                    # Convolve if needed
                    if FWHM_gal is not None and sigma > 0.01:
                        ssp_data = ndimage.gaussian_filter1d(ssp_data, sigma)

                    # Log rebin
                    ssp_new = util.log_rebin(lam_range_temp, ssp_data, velscale=velscale)[0]

                    # Normalize if requested
                    if norm_range is not None:
                        flux_val = np.mean(ssp_new[band])
                        flux[j, k, h] = flux_val
                        ssp_new /= flux_val
                    else:
                        flux[j, k, h] = 1.0

                    # Store results
                    templates[:, j, k, h] = ssp_new
                    age_grid[j, k, h] = age
                    metal_grid[j, k, h] = met
                    alpha_grid[j, k, h] = afe

        # Apply optional masks (age, metallicity, alpha)
        if age_range is not None:
            w = (age_range[0] <= ages) & (ages <= age_range[1])
            templates = templates[:, w, :, :]
            age_grid = age_grid[w, :, :]
            metal_grid = metal_grid[w, :, :]
            alpha_grid = alpha_grid[w, :, :]
            flux = flux[w, :, :]
            ages = ages[w]

        if metal_range is not None:
            w = (metal_range[0] <= metals) & (metals <= metal_range[1])
            templates = templates[:, :, w, :]
            age_grid = age_grid[:, w, :]
            metal_grid = metal_grid[:, w, :]
            alpha_grid = alpha_grid[:, w, :]
            flux = flux[:, w, :]
            metals = metals[w]

        if afe_range is not None:
            w = (afe_range[0] <= alphas) & (alphas <= afe_range[1])
            templates = templates[:, :, :, w]
            age_grid = age_grid[:, :, w]
            metal_grid = metal_grid[:, :, w]
            alpha_grid = alpha_grid[:, :, w]
            flux = flux[:, :, w]
            alphas = alphas[w]

        if norm_range is None:
            flux_median = np.median(templates[templates > 0])
            templates /= flux_median

        # Store full arrays
        self.templates_full = templates
        self.ln_lam_temp_full = ln_lam_temp
        self.lam_temp_full = lam_temp

        # Apply optional wave range
        if wave_range is not None:
            wave_mask = (wave_range[0] <= lam_temp) & (lam_temp <= wave_range[1])
            ln_lam_temp = ln_lam_temp[wave_mask]
            lam_temp = lam_temp[wave_mask]
            templates = templates[wave_mask, :, :, :]

        # Final attributes
        self.templates = templates
        self.ln_lam_temp = ln_lam_temp
        self.lam_temp = lam_temp
        self.age_grid = age_grid
        self.metal_grid = metal_grid
        self.alpha_grid = alpha_grid
        self.n_ages, self.n_metal, self.n_alpha = age_grid.shape
        self.flux = flux

    ###############################################################################
    def mean_age_metal(self, weights, lg_age=True, lg_met = True, quiet=False):
        """
        Compute mean age, mean metallicity, and mean alpha from pPXF weights.
        :param weights: pPXF weights array with dimensions [n_ages, n_metal, n_alpha]
        :param lg_age: If True, compute mean log10(Age/yr); else linear Gyr.
        :param quiet: If True, suppress printout.
        :return: (mean_age, mean_metal, mean_afe)
        """
        if lg_age:
            lg_age_grid = np.log10(self.age_grid) + 9
            mean_lg_age = np.sum(weights * lg_age_grid) / np.sum(weights)
            mean_age = mean_lg_age
        else:
            lin_age_grid = self.age_grid
            mean_lin_age = np.sum(weights * lin_age_grid) / np.sum(weights)
            mean_age = mean_lin_age

        if lg_met:
            mean_metal = np.sum(weights * self.metal_grid) / np.sum(weights)
            mean_afe = np.sum(weights * self.alpha_grid) / np.sum(weights)
        else:
            lin_met_grid = 10**self.metal_grid
            mean_metal = np.sum(weights * lin_met_grid) / np.sum(weights)
            mean_metal = np.log10(mean_metal)

            lin_alpha_grid = 10**self.alpha_grid
            mean_afe = np.sum(weights * lin_alpha_grid) / np.sum(weights)
            mean_afe = np.log10(mean_afe)


        if not quiet:
            if lg_age:
                print(f'Weighted lg <Age>: {mean_age:.3g}')
                print(f'Weighted <[M/H]>: {mean_metal:.3g}')
                print(f'Weighted <[α/Fe]>: {mean_afe:.3g}')
            else:
                print(f'Weighted <Age> [Gyr]: {mean_age:.3g}')
                print(f'Weighted <[M/H]>: {mean_metal:.3g}')
                print(f'Weighted <[α/Fe]>: {mean_afe:.3g}')

        return mean_age, mean_metal, mean_afe

    ###############################################################################
    def get_full_age_grid(self):
        return self.age_grid

    def get_full_metal_grid(self):
        return self.metal_grid

    def get_full_alpha_grid(self):
        return self.alpha_grid


###############################################################################
# Class for XSHOOTER templates
###############################################################################
class xshooter:
    """
    Class to handle SSP templates from the XSHOOTER spectral library (XSL).
    """
    def __init__(self, pathname, velscale, FWHM_gal=None, FWHM_tem=0.5,
                 age_range=None, metal_range=None, norm_range=None, wave_range=None, R = None):
        files = glob.glob(pathname)
        if not files:
            raise FileNotFoundError(f"No files found matching pattern: {pathname}")

        # Extract ages and metallicities
        age_metal_list = [age_metal_xshooter(f) for f in files]
        all_ages = np.array([am[0] for am in age_metal_list])
        all_metals = np.array([am[1] for am in age_metal_list])
        ages = np.unique(all_ages)
        metals = np.unique(all_metals)
        n_ages = len(ages)
        n_metal = len(metals)

        # Dictionary for quick file lookup
        file_map = {}
        for i, f in enumerate(files):
            file_map[(all_ages[i], all_metals[i])] = f

        # Open first file to initialize wavelength
        with fits.open(files[0]) as hdu:
            ssp = hdu[0].data
            header = hdu[0].header

        lam = header['CRVAL1'] + np.arange(header['NAXIS1']) * header['CDELT1']

        # Define normalization band
        if norm_range is not None:
            band = (norm_range[0] <= lam) & (lam <= norm_range[1])

        # If wave_range is provided, apply mask to lam and ssp
        if wave_range is not None:
            mask = (lam >= wave_range[0]) & (lam <= wave_range[1])
        else:
            mask = np.ones_like(lam, dtype=bool)

        # Cutting the templates to the galaxy spectra range
        lam = lam[mask]
        ssp = ssp[mask]

        # Recompute FWHM_tem
        lam_range_temp = lam[[0, -1]]

        # Log-rebin
        ssp_new, ln_lam_temp = util.log_rebin(lam_range_temp, ssp, velscale=velscale)[:2]
        lam_temp = np.exp(ln_lam_temp)


        # Initialize arrays
        templates = np.empty((ssp_new.size, n_ages, n_metal))
        age_grid = np.empty((n_ages, n_metal))
        metal_grid = np.empty((n_ages, n_metal))
        flux = np.empty((n_ages, n_metal))


        # Considering three possibilities: Resolution in FWHM, in R and in MUSE LSF
        if FWHM_gal is not None:
            if isinstance(FWHM_gal, (int, float)): # FWHM resolution
                # print('FWHM resolution')
                FWHM_tem = lam/10000.
                FWHM_dif = np.sqrt(FWHM_gal**2 - FWHM_tem**2)
                sigma = FWHM_dif/2.355/header['CDELT1']   # Sigma difference in pixels
            elif R is not None:
                # print('R resolution')
                FWHM_gal = lam/R #FWHM_gal[mask]
                FWHM_tem = lam/10000.
                FWHM_dif = np.sqrt(FWHM_gal**2 - FWHM_tem**2)
                sigma = FWHM_dif/2.355/header['CDELT1']   # Sigma difference in pixels
            else:
                # print('MUSE resolution')
                FWHM_gal = 5.866e-8*lam**2-9.187e-4*lam+6.040
                FWHM_dif = np.sqrt(FWHM_gal**2 - FWHM_tem**2)
                sigma = FWHM_dif/2.355/header['CDELT1']   # Sigma difference in pixels

        # Process each template
        warn_conv_failed = False
        for j, age in enumerate(ages):
            for k, met in enumerate(metals):
                fname = file_map.get((age, met), None)
                if fname is None:
                    raise ValueError(f"Missing template for age={age}, metallicity={met}")

                with fits.open(fname) as hdu:
                    ssp_data_full = hdu[0].data
                    hdr = hdu[0].header

                # Apply the same wavelength masking for each template
                ssp_data = ssp_data_full[mask]

                # Convolving
                if FWHM_gal is not None:
                    if np.isscalar(FWHM_gal):
                        if np.any(sigma) > 0.01:   # Skip convolution for nearly zero sigma
                            try:
                                x = np.arange(len(ssp_data))
                                ssp_data = util.varsmooth(x, ssp_data, sigma)
                            except ValueError:
                                if not warn_conv_failed:
                                    print('WARNING: The resolution of the templates is lower than galaxy. Skipping convolution')
                                    warn_conv_failed = True
                        else:
                            if not warn_conv_failed:
                                print('WARNING: The resolution of the templates is lower than galaxy. Skipping convolution')
                                warn_conv_failed = True
                    else:
                        x = np.arange(len(ssp_data))
                        try:
                            ssp_data = util.varsmooth(x, ssp_data, sigma) # convolution with variable sigma
                        except ValueError:
                            if not warn_conv_failed:
                                print('WARNING: The resolution of the templates is lower than galaxy. Skipping convolution')
                                warn_conv_failed = True


                # Log rebin the masked data
                ssp_new = util.log_rebin(lam_range_temp, ssp_data, velscale=velscale)[0]

                # Normalize
                if norm_range is not None:
                    flux_val = np.mean(ssp_data_full[band]) #Using the NON cutted SSP templates for flux normalization in the V band
                    flux[j, k] = flux_val
                    ssp_new /= flux_val
                else:
                    flux[j, k] = 1.0

                templates[:, j, k] = ssp_new
                age_grid[j, k] = age
                metal_grid[j, k] = met

        # Apply optional age_range filter
        if age_range is not None:
            w = (age_range[0] <= ages) & (ages <= age_range[1])
            templates = templates[:, w, :]
            age_grid = age_grid[w, :]
            metal_grid = metal_grid[w, :]
            flux = flux[w, :]
            ages = ages[w]

        # Apply optional metal_range filter
        if metal_range is not None:
            w = (metal_range[0] <= metals) & (metals <= metal_range[1])
            templates = templates[:, :, w]
            age_grid = age_grid[:, w]
            metal_grid = metal_grid[:, w]
            flux = flux[:, w]
            metals = metals[w]

        # Global normalization if norm_range is None
        if norm_range is None:
            flux_median = np.median(templates[templates > 0])
            templates /= flux_median

        # Store attributes
        self.templates_full = templates
        self.ln_lam_temp_full = ln_lam_temp
        self.lam_temp_full = lam_temp
        self.templates = templates
        self.ln_lam_temp = ln_lam_temp
        self.lam_temp = lam_temp
        self.age_grid = age_grid
        self.metal_grid = metal_grid
        self.n_ages, self.n_metal = age_grid.shape
        self.flux = flux

    ###############################################################################
    def mean_age_metal(self, weights, lg_age=True, lg_met = True, quiet=False):
        """
        Compute mean age and metallicity from pPXF weights.
        :param weights: pPXF weights array with dimensions [n_ages, n_metal]
        :param lg_age: If True, compute log10(Age/yr) mean; else linear (Gyr).
        :param quiet: If True, suppress printout.
        :return: (mean_age, mean_metal)
        """
        if lg_age:
            lg_age_grid = np.log10(self.age_grid) + 9
            mean_lg_age = np.sum(weights * lg_age_grid) / np.sum(weights)
            mean_age = mean_lg_age
        else:
            lin_age_grid = self.age_grid
            mean_lin_age = np.sum(weights * lin_age_grid) / np.sum(weights)
            mean_age = mean_lin_age

        if lg_met:
            mean_metal = np.sum(weights * self.metal_grid) / np.sum(weights)
        else:
            lin_met_grid = 10**self.metal_grid
            mean_metal = np.sum(weights * lin_met_grid) / np.sum(weights)
            mean_metal = np.log10(mean_metal)

        if not quiet:
            if lg_age:
                print(f'Weighted lg <Age>: {mean_age:.3g}')
                print(f'Weighted <[M/H]>: {mean_metal:.3g}')
            else:
                print(f'Weighted <Age> [Gyr]: {mean_age:.3g}')
                print(f'Weighted <[M/H]>: {mean_metal:.3g}')

        return mean_age, mean_metal

    ###############################################################################
    def get_full_age_grid(self):
        """
        Return the 2D array of ages.
        """
        return self.age_grid

    def get_full_metal_grid(self):
        """
        Return the 2D array of metallicities.
        """
        return self.metal_grid


class KinematicTemplates:
    def __init__(self, pathname, velscale, FWHM_gal=None, FWHM_tem=2.51,
                 norm_range=None, wave_range=None, R=None):
        """
        Universal class for loading generic 1D spectral templates for kinematic analysis using pPXF.
        Templates must be 1D FITS spectra covering the same wavelength range as the observed spectrum.
        """

        # Find all FITS files in the given path
        files = sorted(glob.glob(pathname))
        if not files:
            print(f"No template found in: {pathname}")
            return

        # Read the first file to obtain the wavelength axis
        with fits.open(files[0]) as hdu:
            ssp = hdu[0].data
            hdr = hdu[0].header
        lam = hdr['CRVAL1'] + np.arange(hdr['NAXIS1']) * hdr['CDELT1']

        # Apply wavelength trimming (with ±2% margin)
        if wave_range is not None:
            wave_min, wave_max = wave_range
            delta = 0.02 * (wave_max - wave_min)
            cut_min, cut_max = wave_min - delta, wave_max + delta
            mask = (lam >= cut_min) & (lam <= cut_max)
            lam = lam[mask]
        else:
            mask = slice(None)

        # Define the wavelength range of the template
        lam_range_temp = lam[[0, -1]]

        templates = []
        names = []

        warn_conv_failed = False
        for f in files:
            with fits.open(f) as hdu:
                ssp = hdu[0].data[mask]

            # Apply convolution if required
            if FWHM_gal is not None:
                if isinstance(FWHM_gal, (int, float)):
                    FWHM_dif = np.sqrt(FWHM_gal**2 - FWHM_tem**2)
                    sigma = FWHM_dif / 2.355 / hdr['CDELT1']
                    if sigma > 0.01:
                        try:
                            x = np.arange(len(ssp))
                            ssp = util.varsmooth(x, ssp, sigma)
                        except ValueError:
                            if not warn_conv_failed:
                                print('WARNING: The resolution of the templates is lower than galaxy. Skipping convolution')
                                warn_conv_failed = True
                    else:
                        if not warn_conv_failed:
                            print('WARNING: The resolution of the templates is lower than galaxy. Skipping convolution')
                            warn_conv_failed = True

                elif R is not None:
                    lam_local = lam
                    FWHM_gal = lam_local / R
                    FWHM_dif = np.sqrt(FWHM_gal**2 - FWHM_tem**2)
                    sigma = FWHM_dif / 2.355 / hdr['CDELT1']
                    try:
                        x = np.arange(len(ssp))
                        ssp = util.varsmooth(x, ssp, sigma)
                    except ValueError:
                        if not warn_conv_failed:
                            print('WARNING: The resolution of the templates is lower than galaxy. Skipping convolution')
                            warn_conv_failed = True
                else:
                    # Example: parametric MUSE LSF
                    lam_local = lam
                    FWHM_gal = 5.866e-8 * lam_local**2 - 9.187e-4 * lam_local + 6.040
                    FWHM_dif = np.sqrt(FWHM_gal**2 - FWHM_tem**2)
                    sigma = FWHM_dif / 2.355 / hdr['CDELT1']
                    try:
                        x = np.arange(len(ssp))
                        ssp = util.varsmooth(x, ssp, sigma)
                    except ValueError:
                        if not warn_conv_failed:
                            print('WARNING: The resolution of the templates is lower than galaxy. Skipping convolution')
                            warn_conv_failed = True

            # Log-rebin the template
            ssp_rebinned, ln_lam_temp = util.log_rebin(lam_range_temp, ssp, velscale=velscale)[:2]

            # Normalise if a normalisation range is provided
            if norm_range is not None:
                lam_temp = np.exp(ln_lam_temp)
                norm_mask = (lam_temp >= norm_range[0]) & (lam_temp <= norm_range[1])
                mean_flux = np.mean(ssp_rebinned[norm_mask])
                if mean_flux > 0:
                    ssp_rebinned /= mean_flux

            templates.append(ssp_rebinned)
            names.append(f)

        self.templates = np.column_stack(templates)
        self.ln_lam_temp = ln_lam_temp
        self.lam_temp = np.exp(ln_lam_temp)
        self.names = names
        print(f"{len(files)} templates successfully loaded from: {pathname}")
        print('')
