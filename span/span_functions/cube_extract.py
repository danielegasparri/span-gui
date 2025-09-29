#SPectral ANalysis software (SPAN)
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
######################## THE FOLLOWING FUNCTIONS HAVE BEEN INSPIRED BY THE GIST PIPELINE OF BITTNER ET AL., 2019 #########################
############################################# A special thanks to Adrian Bittner ########################################################


#Functions to bin and extract 1D spectra from datacubes, using the GIST pipeline logic.
#The results are fully compatible with the GIST pipeline.

import os
import shutil
import numpy as np
from astropy.io import fits
import sys
import importlib.util
import functools
from vorbin.voronoi_2d_binning import voronoi_2d_binning
import scipy.spatial.distance as dist
import matplotlib.pyplot as plt




#function to create the dictionary (config) following the GIST standard to be passed to the following functions
def buildConfigFromGUI(ifs_run_id, ifs_input, ifs_output, ifs_redshift, ifs_ow_output,
                       ifs_routine_read, ifs_origin, ifs_lmin_tot, ifs_lmax_tot,
                       ifs_lmin_snr, ifs_lmax_snr, ifs_min_snr_mask,
                       ifs_mask, ifs_bin_method, ifs_target_snr, ifs_covariance, ell_pa_astro_deg=None, ell_x0=None, ell_y0=None, ell_q=None, ell_min_dr=0.5, ell_r_max=None):

    """
    Returns a `configs` dictionary to be read from the following functions of the module

    """

    configs = {
        "INFO": {
            "RUN_NAME": ifs_run_id,
            "INPUT": ifs_input,
            "OUTPUT": ifs_output,
            "REDSHIFT": ifs_redshift,
            "OW_OUTPUT": ifs_ow_output
        },
        "READ": {
            "ROUTINE": ifs_routine_read,
            "ORIGIN": ifs_origin,
            "LMIN_TOT": ifs_lmin_tot,
            "LMAX_TOT": ifs_lmax_tot,
            "LMIN_SNR": ifs_lmin_snr,
            "LMAX_SNR": ifs_lmax_snr
        },
        "MASKING": {
            "MASK_SNR": ifs_min_snr_mask,
            "MASK": ifs_mask
        },
        "BINNING": {
            "VORONOI": ifs_bin_method,         
            "BIN_METHOD": ifs_bin_method,           
            "TARGET_SNR": ifs_target_snr,
            "COVARIANCE": ifs_covariance
        },
        "ELLIPTICAL": {
            "PA_ASTR_DEG": ell_pa_astro_deg,       
            "X0": ell_x0, "Y0": ell_y0,            
            "Q": ell_q,                        
            "MIN_DR": float(ell_min_dr),
            "R_MAX": ell_r_max                     
        }
    }

    return configs


def reading_data(config):

    """
    Reads the datacube using the specified method from the configuration.

    Parameters:
        config (dict): Configuration dictionary containing method and input details.

    Returns:
        cube (object): The loaded datacube, or "Failed" in case of failure.
    """

    print("Step 1: Reading the datacube")


    method = config.get('READ', {}).get('ROUTINE', '')
    method_nopath = os.path.splitext(os.path.basename(method))[0]

    if not method:
        print("No read-in method specified.")
        return "Failed"
    try:
        spec = importlib.util.spec_from_file_location("", method)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        print(f"Using the read-in routine for {method_nopath}")
        return module.read_cube(config)
    except Exception as e:
        print(f"Failed to import or execute the read-in routine {method_nopath}: {e}")
        return "Failed"


def masking(config, cube, preview, manual_bin, existing_bin):

    """
    Applies a spatial mask to the datacube if required.

    Parameters:
        config (dict): Configuration dictionary.
        cube (object): The loaded datacube.
        preview (bool): If True, performs a preview without saving.

    Returns:
        None
    """

    print("\nStep 2: Applying masking, if any")

    output_file = os.path.join(config['INFO']['OUTPUT'], f"{config['INFO']['RUN_NAME']}_mask.fits")

    if (not preview and os.path.isfile(output_file) and not config['INFO'].get('OW_OUTPUT', False) and manual_bin) or existing_bin:
        print("Masking results already exist. Skipping step.")
        return

    generate_and_apply_mask(config, cube)


def binning(config, cube, preview, manual_bin, existing_bin):

    """
    Applies spatial binning to the datacube.

    Parameters:
        config (dict): Configuration dictionary.
        cube (object): The loaded datacube.
        preview (bool): If True, performs a preview without saving.

    Returns:
        None or "Failed" in case of failure.
    """

    print("\nStep 3: Applying binning")

    output_file = os.path.join(config['INFO']['OUTPUT'], f"{config['INFO']['RUN_NAME']}_table.fits")

    if (not preview and (os.path.isfile(output_file) and not config['INFO']['OW_OUTPUT'] and manual_bin)) or existing_bin:
        if not existing_bin:
            print("Results of the module are already in the output directory. Module is skipped.")
            return
        if existing_bin:
            print('Using user mask and bin info')
            return

    try:
        generate_bins(config, cube, False)
    except Exception as e:
        print(f"Spatial binning routine {config.get('BINNING', {}).get('VORONOI', 'UNKNOWN')} failed: {e}")
        return "Failed"


def save_spectra(config, cube, preview, existing_bin):

    """
    Extracts and saves 1D spectra from the datacube.

    Parameters:
        config (dict): Configuration dictionary.
        cube (object): The loaded datacube.
        preview (bool): If True, performs a preview without saving.

    Returns:
        None or "Failed" in case of failure.
    """

    print("\nStep 4: Saving the extracted 1D spectra")

    output_prefix = os.path.join(config['INFO']['OUTPUT'], config['INFO']['RUN_NAME'])
    output_file = f"{output_prefix}_BinSpectra_linear.fits"

    if (not preview and os.path.isfile(output_file) and not config['INFO'].get('OW_OUTPUT', False)) and not existing_bin:
        print("Spectra extraction results already exist. Skipping step.")
        return

    try:
        prepare_mask_bin(config, cube, preview)
    except Exception as e:
        # print(f"Spectra preparation routine {config.get('EXTRACTING', {}).get('MODE', 'UNKNOWN')} failed: {e}")
        print(f"Spectra preparation routine failed: {e}")
        return "Failed"

def save_image(config, cube, preview):
    """
    Collapse the cube['signal'] along the spectral axis and save a 2D image (ny, nx) in FITS format.

    Parameters
    ----------
    config : dict
        Configuration dictionary containing 'INFO' > 'OUTPUT' and 'RUN_NAME'.
    cube : dict
        Datacube structure with keys: 'signal', 'x', 'y', 'wave', etc.

    Returns
    -------
    str
        Path to the saved FITS file.
    """
    signal = cube['signal']  # shape: (nz, nspax)
    x = cube['x']  # shape: (nspax,)
    y = cube['y']  # shape: (nspax,)

    # Collapse along wavelength axis (axis 0) → shape: (nspax,)
    collapsed_flux = signal #np.nansum(signal, axis=0)

    # Create 2D image grid
    x_unique = np.sort(np.unique(x))
    y_unique = np.sort(np.unique(y))
    nx = len(x_unique)
    ny = len(y_unique)
    image_2d = np.full((ny, nx), np.nan)

    # Fill 2D image
    for i in range(len(x)):
        xi = np.searchsorted(x_unique, x[i])
        yi = np.searchsorted(y_unique, y[i])
        image_2d[yi, xi] = collapsed_flux[i]

    cube['x_unique'] = x_unique 
    cube['y_unique'] = y_unique 
    cube['white'] = image_2d       # used by geometry estimators
    cube['shape'] = (ny, nx) 
    
    # Prepare output path, only for the extraction and not the preview
    if not preview:
        output_prefix = os.path.join(config['INFO']['OUTPUT'], config['INFO']['RUN_NAME'])
        output_file = f"{output_prefix}_2dimage.fits"
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        # Save to FITS
        hdu = fits.PrimaryHDU(image_2d)
        hdu.header['COMMENT'] = "2D image collapsed along spectral axis"
        hdu.header['HISTORY'] = "Created with SPAN"
        hdu.writeto(output_file, overwrite=True)
        print('')
        print(f"Saved 2D image to: {output_file}")

        return output_file

###############################################################################
################ FUNCTIONS TO PERFORM THE 4 STEPS ABOVE #############

def generate_and_apply_mask(config, cube):

    """
    Creates a combined mask for the datacube, masking defunct spaxels,
    spaxels below a SNR threshold, and spaxels from an external mask file.

    Parameters:
        config (dict): Configuration dictionary.
        cube (dict): Datacube containing spectral and SNR data.

    Returns:
        None (Saves the mask to a FITS file)
    """

    print("Generating spatial mask...")

    # Mask spaxels that contain NaN values or have a non-positive median flux
    spec = cube['spec']
    median_flux = np.nanmedian(spec, axis=0)
    masked_defunct = np.logical_or(np.any(np.isnan(spec), axis=0), median_flux <= 0)
    print(f"Masking defunct spaxels: {np.sum(masked_defunct)} spaxels are rejected.")

    # Mask spaxels based on signal-to-noise ratio (SNR) threshold
    masked_snr = mask_snr(cube['snr'], cube['signal'], config['MASKING']['MASK_SNR'])

    # Mask spaxels based on an external mask file
    mask_filename = config['MASKING'].get('MASK')
    if mask_filename:
        mask_path = os.path.join(os.path.dirname(config['INFO']['INPUT']), mask_filename)

        if os.path.isfile(mask_path):
            mask_data = fits.getdata(mask_path, ext=1).flatten()
            masked_mask = mask_data == 1
            print(f"Masking spaxels according to mask file: {np.sum(masked_mask)} spaxels are rejected.")
        else:
            print(f"Mask file not found: {mask_path}")
            masked_mask = np.zeros(len(cube['snr']), dtype=bool)
    else:
        print("No mask file provided.")
        masked_mask = np.zeros(len(cube['snr']), dtype=bool)

    # Combine all masks
    combined_mask = np.logical_or.reduce((masked_defunct, masked_snr, masked_mask))

    # Save final mask
    save_mask(combined_mask, masked_defunct, masked_snr, masked_mask, config)


def save_mask(combined_mask, masked_defunct, masked_snr, masked_mask, config):

    """
    Saves the final combined mask and its components to a FITS file.

    Parameters:
        combined_mask (np.ndarray): Boolean array of the final combined mask.
        masked_defunct (np.ndarray): Boolean array for defunct spaxels.
        masked_snr (np.ndarray): Boolean array for SNR-masked spaxels.
        masked_mask (np.ndarray): Boolean array for external mask file spaxels.
        config (dict): Configuration dictionary.

    Returns:
        None (Writes the mask to a FITS file)
    """

    output_file = os.path.join(config['INFO']['OUTPUT'], f"{config['INFO']['RUN_NAME']}_mask.fits")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    print(f"Writing mask file: {output_file}")

    # reading and writing the fits
    with fits.HDUList([
        fits.PrimaryHDU(),
        fits.BinTableHDU.from_columns([
            fits.Column(name='MASK', format='I', array=combined_mask.astype(int)),
            fits.Column(name='MASK_DEFUNCT', format='I', array=masked_defunct.astype(int)),
            fits.Column(name='MASK_SNR', format='I', array=masked_snr.astype(int)),
            fits.Column(name='MASK_FILE', format='I', array=masked_mask.astype(int))
        ], name="MASKFILE")
    ]) as hdul:
        # Comments in the header
        hdul[1].header['COMMENT'] = "Value 0 -> Unmasked"
        hdul[1].header['COMMENT'] = "Value 1 -> Masked"

        hdul.writeto(output_file, overwrite=True)

    print(f"Mask file saved successfully: {output_file}")


def mask_snr(snr, signal, min_snr):

    """
    Masks spaxels based on a minimum SNR threshold.

    Parameters:
        snr (np.ndarray): Array of signal-to-noise ratios for each spaxel.
        signal (np.ndarray): Array of signal values for each spaxel.
        min_snr (float): Minimum SNR threshold for masking.

    Returns:
        masked (np.ndarray): Boolean array indicating masked spaxels.
    """

    # Identify spaxels close to the SNR threshold
    idx_snr = np.where(np.abs(snr - min_snr) < 2)[0]

    if len(idx_snr) > 0:
        meanmin_signal = np.mean(signal[idx_snr])
    else:
        meanmin_signal = np.min(signal)  # Fallback if no matching spaxels

    # Mask spaxels below the calculated signal threshold
    masked = signal < meanmin_signal

    if np.all(masked):
        print("No spaxels with S/N above the threshold. Ignoring potential warnings.")

    return masked


def sn_func(index, signal=None, noise=None, covar_vor=0.00):

    """
    Computes the signal-to-noise ratio in a bin for Voronoi binning.

    Parameters:
        index (np.ndarray): Indices of spaxels in the bin.
        signal (np.ndarray): Signal values for each spaxel.
        noise (np.ndarray): Noise values for each spaxel.
        covar_vor (float, optional): Correction factor for spatial correlations.

    Returns:
        sn (float): Estimated signal-to-noise ratio for the bin.
    """

    total_signal = np.sum(signal[index])
    total_noise = np.sqrt(np.sum(noise[index] ** 2))

    if total_noise == 0:
        return 0  # Prevent division by zero

    sn = total_signal / total_noise

    # Apply correction for spatial correlations
    if index.size > 1 and covar_vor > 0:
        sn /= 1 + covar_vor * np.log10(index.size)

    return sn




def generate_bins(config, cube, voronoi):
    """
    Applies spatial binning according to selected method.
    """

    print("Defining spatial bins")

    mask_file = os.path.join(config['INFO']['OUTPUT'], f"{config['INFO']['RUN_NAME']}_mask.fits")
    if not os.path.isfile(mask_file):
        print(f"Mask file not found: {mask_file}")
        return "Failed"

    with fits.open(mask_file, mode="readonly") as hdul:
        mask_data = hdul[1].data

    mask = mask_data['MASK']
    idx_unmasked = np.where(mask == 0)[0]
    idx_masked   = np.where(mask == 1)[0]

    # Partial S/N function (you already have this)
    sn_func_covariances = functools.partial(
        sn_func, covar_vor=config['BINNING'].get('COVARIANCE', 0.0))

    # ---- Choose method ----
    method = config['BINNING'].get('BIN_METHOD', 'VORONOI').upper()
    # Back-compat: if 'voronoi' bool param is used, override
    if voronoi:
        method = 'VORONOI'
    elif method not in ('VORONOI', 'ELLIPTICAL', 'SPAXEL'):
        method = 'SPAXEL'

    if method == 'VORONOI':
        try:
            bin_num, x_node, y_node, x_bar, y_bar, sn, n_pixels, _ = voronoi_2d_binning(
                cube['x'][idx_unmasked],
                cube['y'][idx_unmasked],
                cube['signal'][idx_unmasked],
                cube['noise'][idx_unmasked],
                config['BINNING']['TARGET_SNR'],
                plot=False,
                quiet=True,
                pixelsize=cube['pixelsize'],
                sn_func=sn_func_covariances)
            print(f"{np.max(bin_num) + 1} Voronoi bins generated!")
        except ValueError as e:
            if str(e) == 'All pixels have enough S/N and binning is not needed':
                print("Analysis will continue without Voronoi-binning!")
                bin_num = np.arange(len(idx_unmasked))
                x_node, y_node = cube['x'][idx_unmasked], cube['y'][idx_unmasked]
                sn = cube['snr'][idx_unmasked]
                n_pixels = np.ones(len(idx_unmasked), dtype=int)
            else:
                print(f"Voronoi-binning error: {e}")
                return "Failed"

    elif method == 'ELLIPTICAL':
        print("Using elliptical annuli binning (adaptive S/N).")

        bin_num, x_node, y_node, sn, n_pixels, R_ellipses = generate_bins_elliptical(
            config, cube, idx_unmasked, sn_func_covariances)
        
        
    else:  # 'SPAXEL' = each unmasked spaxel is a bin
        print(f"No binning! {len(idx_unmasked)} spaxels will be treated as individual bins.")
        bin_num = np.arange(len(idx_unmasked))
        x_node, y_node = cube['x'][idx_unmasked], cube['y'][idx_unmasked]
        sn = cube['snr'][idx_unmasked]
        n_pixels = np.ones(len(idx_unmasked), dtype=int)

    # ---- Masked spaxels behaviour (unchanged): keep -1 in BIN_ID ----
    if len(idx_masked) > 0:
        # You compute nearest Voronoi node; we keep same scaffolding for consistency
        pix_coords = np.column_stack((cube['x'][idx_masked], cube['y'][idx_masked]))
        bin_coords = np.column_stack((x_node, y_node))
        dists = dist.cdist(pix_coords, bin_coords, 'euclidean')
        bin_num_outside = np.argmin(dists, axis=1)
    else:
        bin_num_outside = np.array([])

    # ---- Build long BIN_ID over all spaxels ----
    bin_num_long = np.full(len(cube['x']), np.nan)
    bin_num_long[idx_unmasked] = bin_num
    bin_num_long[idx_masked]   = -1

    # ---- Save (same as before) ----
    save_bin_info(
        config,
        cube['x'], cube['y'], cube['signal'], cube['snr'],
        bin_num_long, np.unique(bin_num), x_node, y_node, sn, n_pixels, cube['pixelsize'])



def generate_bins_elliptical(config, cube, idx_unmasked, sn_func_covariances):
    target_snr = float(config['BINNING'].get('TARGET_SNR', 30.0))
    min_dr     = float(config['ELLIPTICAL'].get('MIN_DR', 1.0))

    geo = config.get('ELLIPTICAL', {})
    x0_user  = geo.get('X0', None) 
    y0_user  = geo.get('Y0', None) 
    q_user   = geo.get('Q', None)
    pa_img_user  = geo.get('PA_IMAGE_DEG', None)        
    pa_astro     = geo.get('PA_ASTR_DEG', None)         
    pa_user = None
    if pa_img_user is not None:
        pa_user = float(pa_img_user)                    
    elif pa_astro is not None:
        pa_user = _pa_astro_to_image(float(pa_astro)) 

    white = cube.get('white', None) 
    H, W  = (cube['shape'] if 'shape' in cube else (None, None))

    # --- Grids directly on the stacked image
    X_unique = cube.get('x_unique', None)
    Y_unique = cube.get('y_unique', None)
    if (X_unique is None) or (Y_unique is None):
        X_unique = np.sort(np.unique(cube['x'])) 
        Y_unique = np.sort(np.unique(cube['y'])) 
        # Se vuoi riusarli altrove:
        cube['x_unique'] = X_unique
        cube['y_unique'] = Y_unique

    X_grid, Y_grid = np.meshgrid(X_unique, Y_unique, indexing='xy')

    # mask
    mask_valid_full = np.zeros(H*W, dtype=bool)
    mask_valid_full[idx_unmasked] = True
    mask_valid_2d = mask_valid_full.reshape(H, W)

    # --- Geometry
    if (x0_user is None) or (y0_user is None) or (pa_user is None) or (q_user is None):
        if white is not None:
            x0e, y0e, pae_img, qe = _estimate_centre_pa_q_arcsec(
                white, mask_valid_2d, X_grid, Y_grid,
                x0_user_arcsec=x0_user, y0_user_arcsec=y0_user)
        else:
            # fallback
            x0e = float(np.median(cube['x'][idx_unmasked]))
            y0e = float(np.median(cube['y'][idx_unmasked]))
            pae_img, qe = 0.0, 1.0
        x0 = float(x0_user) if (x0_user is not None) else x0e
        y0 = float(y0_user) if (y0_user is not None) else y0e
        pa = float(pa_user) if (pa_user is not None) else pae_img  
        q  = float(q_user)  if (q_user  is not None) else qe
    else:
        x0 = float(x0_user);  y0 = float(y0_user)
        pa = float(pa_user);  q  = float(q_user)

    # --- elliptical radius
    x_u = cube['x'][idx_unmasked]   # arcsec
    y_u = cube['y'][idx_unmasked]   # arcsec
    r_u = _elliptical_radius(x_u, y_u, x0, y0, pa, q) 

    # --- r_max (arcsec) ---
    r_max_cfg = _as_float_or_none(config.get('ELLIPTICAL', {}).get('R_MAX'))
    if r_max_cfg is not None:
        r_max = float(r_max_cfg)
    else:
        r_max = float(np.nanpercentile(r_u, 99.5))

    # --- Ordering
    order    = np.argsort(r_u)
    r_sorted = r_u[order]

    sig_u = cube['signal'][idx_unmasked]
    noi_u = cube['noise'][idx_unmasked]

    bin_edges = []
    bin_sels  = []

    i0 = 0
    r_in = 0.0
    
    pixelsize =  cube['pixelsize']
    min_n_spaxels = 1
    if min_dr < pixelsize:
        print ('WARNING: minimum dr smaller than the sampling. Adjusting...')
        min_dr = pixelsize
        
    while r_in < r_max:
        # First attempt: an annulus at least min_dr_arcsec thick
        r_out = min(r_in + min_dr, r_max)

        # Expand until at least 1 spaxel is included (handles thin shells vs. sampling)
        while True:
            sel = (r_u >= r_in) & (r_u < r_out)
            if sel.any() or r_out >= r_max:
                break
            r_out = min(r_out + pixelsize, r_max)

        sel = (r_u >= r_in) & (r_u < r_out)
        if not sel.any():
            break  # nothing left to bin

        # Grow bin
        idx_sel = np.where(sel)[0]

        def _sn_for_idx(idxs):
            if sn_func_covariances is not None:
                return float(sn_func_covariances(np.asarray(idxs, dtype=int), sig_u, noi_u))
            else:
                # fallback
                S = sig_u[idxs].sum()
                N = np.sqrt((noi_u[idxs]**2).sum())
                return (S / N) if N > 0 else 0.0

        sn_bin = _sn_for_idx(idx_sel)
        n_spx  = idx_sel.size

        while (sn_bin < target_snr or n_spx < min_n_spaxels) and r_out < r_max:
            r_out = min(r_out + min_dr, r_max)
            sel   = (r_u >= r_in) & (r_u < r_out)
            idx_sel = np.where(sel)[0]
            if idx_sel.size == n_spx and r_out >= r_max:
                break
            n_spx  = idx_sel.size
            sn_bin = _sn_for_idx(idx_sel)

        # Record bin
        if sel.any():
            bin_edges.append((float(r_in), float(r_out)))
            bin_sels.append(sel)

        # Advance
        r_in = r_out
        if r_out >= r_max:
            break

    # --- Outputs ---
    nb = len(bin_sels)
    bin_num  = np.full(x_u.size, -1, dtype=int)
    x_node   = np.zeros(nb, dtype=float)
    y_node   = np.zeros(nb, dtype=float)
    sn_arr   = np.zeros(nb, dtype=float)
    n_sp_arr = np.zeros(nb, dtype=int)

    # Elliptical radius and coordinates along the major axis
    R_flux = np.zeros(nb, dtype=float)
    X_flux = np.zeros(nb, dtype=float)
    Y_flux = np.zeros(nb, dtype=float)
    th = np.deg2rad(pa)

    for b, sel in enumerate(bin_sels):
        idx = np.where(sel)[0]
        bin_num[idx] = b
        n_sp_arr[b]  = idx.size
        x_node[b]    = float(np.mean(x_u[idx]))
        y_node[b]    = float(np.mean(y_u[idx]))
        S = sig_u[idx].sum()
        N = np.sqrt((noi_u[idx]**2).sum())
        sn_arr[b] = (S / N) if N > 0 else 0.0

        # --- R_flux
        r_vals = r_u[idx]
        w      = sig_u[idx] 
        w_ok   = np.isfinite(w) & (w >= 0)

        if np.any(w_ok) and np.sum(w[w_ok]) > 0:
            R_flux[b] = float(np.average(r_vals[w_ok], weights=w[w_ok]))
        else:
            R_flux[b] = float(np.mean(r_vals))  # fallback

        # X and Y coordinates of R_flux on semi-major axis
        X_flux[b] = x0 + R_flux[b] * np.cos(th)
        Y_flux[b] = y0 + R_flux[b] * np.sin(th)

    print("PA:", float(pa+90), "   q:", float(q))
    print(f"{nb} Elliptical bins generated!\n")

    return bin_num, X_flux, Y_flux, sn_arr, n_sp_arr, R_flux

 
def _pa_astro_to_image(pa_astro_deg):
    #Convert astronomical PA (E of N, 0..180) to image-frame PA used by the code
    if pa_astro_deg is None:
        return None
    return float(((pa_astro_deg)-90) % 180.0)


def _estimate_centre_pa_q_arcsec(white, mask_valid_2d, X_grid, Y_grid,
                                 x0_user_arcsec=None, y0_user_arcsec=None):
    """
    Estimate centre (arcsec), position angle (image convention, degrees 0..180),
    and axis ratio q directly in physical coordinates using second-order
    moments of the white-light image.

    Parameters
    ----------
    white : (ny, nx) float
        Collapsed (white-light) image.
    mask_valid_2d : (ny, nx) bool
        True where the pixel is valid.
    X_grid, Y_grid : (ny, nx) float
        2D grids of coordinates in arcsec for each pixel of the white image.
    x0_user_arcsec, y0_user_arcsec : float | None
        User-provided centre in arcsec (same convention as X_grid/Y_grid).

    Returns
    -------
    x0_arc, y0_arc : float
        Adopted centre in arcsec.
    pa_image_deg : float
        Position angle in image convention (0° = +X, counter-clockwise),
        normalised to [0, 180).
    q : float
        Axis ratio b/a, clipped to [0.05, 1.0].
    """
    # Use only valid pixels; clip flux to be non-negative
    flux = np.where(mask_valid_2d, np.clip(white, 0, None), 0.0)
    tot  = flux.sum()
    if tot <= 0:
        # Fallback:
        ydim, xdim = white.shape
        x0_arc = float(np.median(X_grid))
        y0_arc = float(np.median(Y_grid))
        return x0_arc, y0_arc, 0.0, 1.0

    # Centre: use user value if provided (arcsec), otherwise luminosity-weighted centroid (arcsec)
    if (x0_user_arcsec is not None) and (y0_user_arcsec is not None):
        x0_arc = float(x0_user_arcsec)
        y0_arc = float(y0_user_arcsec)
    else:
        x0_arc = float((flux * X_grid).sum() / tot)
        y0_arc = float((flux * Y_grid).sum() / tot)

    # Second moments about the centre (arcsec)
    Xc = X_grid - x0_arc
    Yc = Y_grid - y0_arc
    Ixx = float((flux * Xc * Xc).sum() / tot)
    Iyy = float((flux * Yc * Yc).sum() / tot)
    Ixy = float((flux * Xc * Yc).sum() / tot)

    M = np.array([[Ixx, Ixy], [Ixy, Iyy]], dtype=float)
    evals, evecs = np.linalg.eigh(M)

    # Major axis = eigenvector associated with the largest eigenvalue
    i_max = int(np.argmax(evals))
    vx, vy = evecs[:, i_max]

    # Image PA (0° = +X, CCW), normalised to [0, 180)
    pa_rad = np.arctan2(vy, vx)
    pa_image_deg = float((np.degrees(pa_rad) + 180.0) % 180.0)

    # Axis ratio q = sqrt(min/max) with numerical robustness
    a2 = float(np.max(evals))
    b2 = float(np.min(evals))
    a2 = max(a2, 1e-12)
    b2 = max(b2, 1e-12)
    q = float(np.sqrt(b2 / a2))
    q = float(np.clip(q, 0.05, 1.0))

    print (f'\nOffset from the origin (x,y) [arcsec]: {x0_arc, y0_arc}')
    return x0_arc, y0_arc, pa_image_deg, q


def _elliptical_radius(x, y, x0, y0, pa_deg, q):
    """Elliptical radius for 1D coordinate arrays (x,y)."""
    pa = np.deg2rad(pa_deg)
    cos, sin = np.cos(pa), np.sin(pa)
    xp =  (x - x0) * cos + (y - y0) * sin
    yp = -(x - x0) * sin + (y - y0) * cos
    return np.hypot(xp, yp / max(q, 1e-3))


def _as_float_or_none(v):
    if v is None:
        return None
    if isinstance(v, str) and v.strip() == "":
        return None
    return float(v)







































def save_bin_info(config, x, y, signal, snr, bin_num_new, ubins, x_node, y_node, sn, n_pixels, pixelsize):

    """
    Saves Voronoi binning results to a GIST-like FITS file.

    Parameters:
        config (dict): Configuration dictionary.
        x, y (np.ndarray): Spaxel coordinates.
        signal, snr (np.ndarray): Signal and SNR values.
        bin_num_new (np.ndarray): Assigned bin number for each spaxel.
        ubins (np.ndarray): Unique bin numbers.
        x_node, y_node (np.ndarray): Coordinates of bin centroids.
        sn (np.ndarray): SNR per bin.
        n_pixels (np.ndarray): Number of spaxels per bin.
        pixelsize (float): Pixel size for FITS metadata.

    Returns:
        None (Writes the FITS file)
    """

    output_file = os.path.join(config['INFO']['OUTPUT'], f"{config['INFO']['RUN_NAME']}_table.fits")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    print(f"Writing: {config['INFO']['RUN_NAME']}_table.fits")

    # Expand data to spaxel level
    x_node_new = np.zeros_like(x)
    y_node_new = np.zeros_like(y)
    sn_new = np.zeros_like(x, dtype=float)
    n_pixels_new = np.zeros_like(x, dtype=int)

    # escludi bin = -1 dall’elenco unico
    valid_mask = bin_num_new >= 0
    ubins = np.unique(bin_num_new[valid_mask])

    x_node_new = np.zeros_like(x)
    y_node_new = np.zeros_like(y)
    sn_new = np.zeros_like(x, dtype=float)
    n_pixels_new = np.zeros_like(x, dtype=int)

    for i, ubin in enumerate(ubins):
        idx = np.where(bin_num_new == ubin)[0]
        x_node_new[idx] = x_node[i]
        y_node_new[idx] = y_node[i]
        sn_new[idx] = sn[i]
        n_pixels_new[idx] = n_pixels[i]

    # Create FITS table
    columns = [
        fits.Column(name='ID', format='J', array=np.arange(len(x))),
        fits.Column(name='BIN_ID', format='J', array=bin_num_new),
        fits.Column(name='X', format='D', array=x),
        fits.Column(name='Y', format='D', array=y),
        fits.Column(name='FLUX', format='D', array=signal),
        fits.Column(name='SNR', format='D', array=snr),
        fits.Column(name='XBIN', format='D', array=x_node_new),
        fits.Column(name='YBIN', format='D', array=y_node_new),
        fits.Column(name='SNRBIN', format='D', array=sn_new),
        fits.Column(name='NSPAX', format='J', array=n_pixels_new),
    ]

    # writing fits
    with fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU.from_columns(columns, name="TABLE")]) as hdul:
        hdul.writeto(output_file, overwrite=True)

    fits.setval(output_file, "PIXSIZE", value=pixelsize)

    print(f"Wrote Voronoi table: {output_file}")


def prepare_mask_bin(config, cube, preview):

    """
    Reads spatial bins and mask file, applies binning to spectra, and saves or displays the binned spectra.

    Parameters:
        config (dict): Configuration dictionary.
        cube (dict): Datacube containing spectral data.
        preview (bool): If True, only displays the Voronoi map without saving.

    Returns:
        None
    """

    mask_file = os.path.join(config['INFO']['OUTPUT'], f"{config['INFO']['RUN_NAME']}_mask.fits")
    table_file = os.path.join(config['INFO']['OUTPUT'], f"{config['INFO']['RUN_NAME']}_table.fits")

    if not os.path.isfile(mask_file) or not os.path.isfile(table_file):
        print("Error: Mask or binning table file not found.")
        return
 
    # Load mask and binning table
    mask = fits.getdata(mask_file, ext=1)['MASK']
    unmasked_spaxels = np.where(mask == 0)[0]
    with fits.open(table_file, mode="readonly") as hdul:
        bin_num = hdul[1].data['BIN_ID'][unmasked_spaxels]

    # ---------- FIX: drop -1 bins and compact IDs ----------
    valid = (bin_num >= 0)
    if not np.any(valid):
        print("Warning: no spaxels inside r_max; nothing to stack.")
        return

    # Compact bin ids to 0..nb-1 to be Voronoi-like
    unique_bins = np.sort(np.unique(bin_num[valid]))        # e.g. [0,1,2,...]
    remap = {b: i for i, b in enumerate(unique_bins)}
    bin_num_compact = np.array([remap[b] for b in bin_num[valid]], dtype=int)

    # Subselect spectra/errors to the same valid spaxels
    spec_u = cube['spec'][:, unmasked_spaxels][:, valid]    # shape: (nwave, N_valid)
    err_u  = cube['error'][:, unmasked_spaxels][:, valid]   # shape: (nwave, N_valid)
    # -------------------------------------------------------

    # Apply binning to spectra (your stacker)
    print("Applying spatial bins to linear data...")
    bin_data, bin_error, bin_flux = perform_voronoi(  # <-- unchanged API
        bin_num_compact, 
        spec_u, 
        err_u
    )
    print("Applied spatial bins.")

    if not preview:
        save_bin_spec(config, bin_data, bin_error, cube['wave'])
    else:
        # Display Voronoi map
        try:

            with fits.open(table_file, mode="readonly") as hdul:
                data = hdul[1].data

            x, y, bin_id, signal = data['X'], data['Y'], data['BIN_ID'], data['SNRBIN']

            # Set masked spaxels (bin_id < 0) to zero signal
            signal[bin_id < 0] = 0

            # Create grid
            x_bins, y_bins = np.unique(x), np.unique(y)
            grid_data = np.full((len(y_bins), len(x_bins)), np.nan)

            for i in range(len(x)):
                x_idx = np.searchsorted(x_bins, x[i])
                y_idx = np.searchsorted(y_bins, y[i])
                grid_data[y_idx, x_idx] = signal[i]

            # Plot Voronoi map
            plt.figure(figsize=(8, 6))
            plt.pcolormesh(x_bins, y_bins, grid_data, cmap='inferno', shading='auto')
            plt.colorbar(label="S/N")
            plt.xlabel("R [arcsec]")
            plt.ylabel("R [arcsec]")
            plt.title("Voronoi Map")
            plt.show()

        except Exception as e:
            print(f"Error: Unable to display Voronoi map: {e}")


def perform_voronoi(bin_num, spec, error):

    """
    Aggregates spaxels belonging to the same Voronoi bin.

    Parameters:
        bin_num (np.ndarray): Array of bin numbers for each spaxel.
        spec (np.ndarray): Spectral data array.
        error (np.ndarray): Error array for spectra.

    Returns:
        tuple: Binned spectra, errors, and flux.
    """

    ubins = np.unique(bin_num)
    nbins = len(ubins)
    npix = spec.shape[0]

    bin_data = np.zeros((npix, nbins))
    bin_error = np.zeros((npix, nbins))
    bin_flux = np.zeros(nbins)

    for i, ubin in enumerate(ubins):
        k = np.where(bin_num == ubin)[0]

        if k.size == 1:
            av_spec = spec[:, k].ravel()
            av_err_spec = np.sqrt(error[:, k]).ravel()
        else:
            av_spec = np.nansum(spec[:, k], axis=1)
            av_err_spec = np.sqrt(np.nansum(error[:, k] ** 2, axis=1))

        bin_data[:, i] = av_spec
        bin_error[:, i] = av_err_spec
        bin_flux[i] = np.mean(av_spec)

    return bin_data, bin_error, bin_flux


def save_bin_spec(config, spec, error, wavelength):

    """
    Saves binned spectra and error spectra to a FITS file.

    Parameters:
        config (dict): Configuration dictionary.
        spec (np.ndarray): Array of binned spectra.
        error (np.ndarray): Array of error spectra.
        wavelength (np.ndarray): Wavelength array.

    Returns:
        None (Writes the FITS file)
    """

    output_file = os.path.join(config['INFO']['OUTPUT'], f"{config['INFO']['RUN_NAME']}_BinSpectra_linear.fits")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    print(f"Writing: {output_file}")

    # Opening fits and writing
    with fits.HDUList([
        fits.PrimaryHDU(),
        fits.BinTableHDU.from_columns([
            fits.Column(name='SPEC', format=f"{spec.shape[0]}D", array=spec.T),
            fits.Column(name='ESPEC', format=f"{spec.shape[0]}D", array=error.T)
        ], name='BIN_SPECTRA'),
        fits.BinTableHDU.from_columns([
            fits.Column(name='WAVE', format='D', array=wavelength)
        ], name='WAVE')
    ]) as hdul:
        hdul.writeto(output_file, overwrite=True)

    # adding info
    fits.setval(output_file, 'CRPIX1', value=1.0)
    fits.setval(output_file, 'CRVAL1', value=wavelength[0])
    fits.setval(output_file, 'CDELT1', value=wavelength[1] - wavelength[0])

    print(f"Wrote: {output_file}")


def extract(config, preview, voronoi, elliptical, manual_bin, existing_bin):

    """
    Main function to run the extraction steps in sequence.

    Parameters:
        config (dict): Configuration dictionary.
        preview (bool): If True, runs in preview mode without saving.
        voronoi (bool): If True, applies Voronoi binning.
        manual_bin (bool): If True, namual bin has been selected.

    Returns:
        None
    """

    print("\n--- Starting Extraction Process ---\n")

    # 1) Read the datacube
    cube = reading_data(config)
    if cube == "Failed":
        print("Extraction aborted: Failed to read the datacube.")
        return

    # 2) Creating the collapsed image
    save_image(config, cube, preview)
        
    # 3) Apply spatial mask
    masking(config, cube, preview, manual_bin, existing_bin)

    # 4) Apply Voronoi binning
    binning_result = binning(config, cube, preview, manual_bin, existing_bin)
    if binning_result == "Failed":
        print("Extraction aborted: Failed to perform Voronoi binning.")
        return

    # 5) Extract and save spectra
    save_spectra(config, cube, preview, existing_bin)
    
    print("\n--- Extraction Process Completed ---\n")



def handle_existing_bin_files(input_folder, output_dir, ifs_run_id):
    """
    Copy and rename existing bin info files from input_folder to output_dir.

    Parameters
    ----------
    input_folder : str
        Directory where the *_table.fits and *_mask.fits files are located.
    output_dir : str
        Destination directory where the renamed files will be copied.
    ifs_run_id : str
        Identifier to replace the original RUN_NAME in the filenames.
    """

    os.makedirs(output_dir, exist_ok=True)

    # Search for *_table.fits and *_mask.fits files
    table_file = None
    mask_file = None
    try:
        for fname in os.listdir(input_folder):
            if fname.endswith('_table.fits'):
                table_file = fname
            elif fname.endswith('_mask.fits'):
                mask_file = fname
    except Exception as e:
        print('The folder does not exist!')

    if not table_file or not mask_file:
        print("Missing required files",
                       "Could not find both *_table.fits and *_mask.fits in the selected folder.")
        return

    # Define full paths
    table_path = os.path.join(input_folder, table_file)
    mask_path = os.path.join(input_folder, mask_file)

    # Define new names
    new_table_name = f"{ifs_run_id}_table.fits"
    new_mask_name = f"{ifs_run_id}_mask.fits"

    # Define new full paths
    new_table_path = os.path.join(output_dir, new_table_name)
    new_mask_path = os.path.join(output_dir, new_mask_name)

    # Copy and rename files
    shutil.copyfile(table_path, new_table_path)
    shutil.copyfile(mask_path, new_mask_path)

    print(f"Files copied and renamed to:\n{new_table_name}\n{new_mask_name}")



























