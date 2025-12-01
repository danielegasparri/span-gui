import ppxf.ppxf_util as util
import numpy as np


# simple function to find the nearest FWHM_gal value corresponding to the emission line funcion, in case of kinematics with gas and variable FWHM_gal (i.e. constant resolving power R).
def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx

# 11) Modified emission line list from the ppxf_util.py of the ppxf package.
def emission_lines(ln_lam_temp, lam_range_gal, FWHM_gal, pixel=True,
                   tie_balmer=False, limit_doublets=False, vacuum=False, wave_galaxy = None):
    """
    Generates an array of Gaussian emission lines to be used as gas templates in PPXF.

    Daniele Gasparri:
    Added the 'wave_galaxy' array, which is the galaxy wavelength array, needed to compute the FWHM_gal values in the gas
    emission lines when the FWHM_gal is not constant (i.e. when we are working with spectra with
    a fixed resolving power R. Needed for stars and gas kinematics task when selecting 'Spec. constant R resolution:'.)

    ****************************************************************************
    ADDITIONAL LINES CAN BE ADDED BY EDITING THE CODE OF THIS PROCEDURE, WHICH
    IS MEANT AS A TEMPLATE TO BE COPIED AND MODIFIED BY THE USERS AS NEEDED.
    ****************************************************************************


    Output Parameters
    -----------------

    emission_lines: ndarray
        Array of dimensions ``[ln_lam_temp.size, line_wave.size]`` containing
        the gas templates, one per array column.

    line_names: ndarray
        Array of strings with the name of each line, or group of lines'

    line_wave: ndarray
        Central wavelength of the lines, one for each gas template'

    """
    #        Balmer:     H10       H9         H8        Heps    Hdelta    Hgamma    Hbeta     Halpha
    balmer = np.array([3798.983, 3836.479, 3890.158, 3971.202, 4102.899, 4341.691, 4862.691, 6564.632])  # vacuum wavelengths

    if tie_balmer:

        # Balmer decrement for Case B recombination (T=1e4 K, ne=100 cm^-3)
        # from Storey & Hummer (1995) https://ui.adsabs.harvard.edu/abs/1995MNRAS.272...41S
        # In electronic form https://cdsarc.u-strasbg.fr/viz-bin/Cat?VI/64
        # See Table B.7 of Dopita & Sutherland 2003 https://www.amazon.com/dp/3540433627
        # Also see Table 4.2 of Osterbrock & Ferland 2006 https://www.amazon.co.uk/dp/1891389343/
        wave = balmer
        if not vacuum:
            wave = util.vac_to_air(wave)

        #if FWHM_gal is an array, I need to extract the FWHM values corresponding to the emission lines of the gas template
        if isinstance(FWHM_gal, np.ndarray):
            FWHM_gal_line = np.array([FWHM_gal[find_nearest(wave_galaxy, lw)] for lw in wave])
            gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal_line, pixel)
        else:
            gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel)

        ratios = np.array([0.0530, 0.0731, 0.105, 0.159, 0.259, 0.468, 1, 2.86])
        ratios *= wave[-2]/wave  # Account for varying log-sampled pixel size in Angstrom
        emission_lines = gauss @ ratios
        line_names = ['Balmer']
        w = (lam_range_gal[0] < wave) & (wave < lam_range_gal[1])
        line_wave = np.mean(wave[w]) if np.any(w) else np.mean(wave)

    else:

        line_wave = balmer
        if not vacuum:
            line_wave = util.vac_to_air(line_wave)
        line_names = ['(H10)', '(H9)', '(H8)', '(Heps)', '(Hdelta)', '(Hgamma)', '(Hbeta)', '(Halpha)']

        #if FWHM_gal is an array, I need to extract the FWHM values corresponding to the emission lines of the gas template
        if isinstance(FWHM_gal, np.ndarray):
            FWHM_gal_line = np.array([FWHM_gal[find_nearest(wave_galaxy, lw)] for lw in line_wave])
            emission_lines = util.gaussian(ln_lam_temp, line_wave, FWHM_gal_line, pixel)
        else:
            emission_lines = util.gaussian(ln_lam_temp, line_wave, FWHM_gal, pixel)


    if limit_doublets:

        # The line ratio of this doublet lam3727/lam3729 is constrained by
        # atomic physics to lie in the range 0.28--1.47 (e.g. fig.5.8 of
        # Osterbrock & Ferland 2006 https://www.amazon.co.uk/dp/1891389343/).
        # We model this doublet as a linear combination of two doublets with the
        # maximum and minimum ratios, to limit the ratio to the desired range.
        #       -----[OII]-----
        wave = [3727.092, 3729.875]    # vacuum wavelengths
        if not vacuum:
            wave = util.vac_to_air(wave)
        names = ['[OII]3726_d1', '[OII]3726_d2']

        if isinstance(FWHM_gal, np.ndarray):
            FWHM_gal_line = np.array([FWHM_gal[find_nearest(wave_galaxy, lw)] for lw in wave])
            gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal_line, pixel)
        else:
            gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel)

        doublets = gauss @ [[1, 1], [0.28, 1.47]]  # produces *two* doublets
        emission_lines = np.column_stack([emission_lines, doublets])
        line_names = np.append(line_names, names)
        line_wave = np.append(line_wave, wave)

        # The line ratio of this doublet lam6717/lam6731 is constrained by
        # atomic physics to lie in the range 0.44--1.43 (e.g. fig.5.8 of
        # Osterbrock & Ferland 2006 https://www.amazon.co.uk/dp/1891389343/).
        # We model this doublet as a linear combination of two doublets with the
        # maximum and minimum ratios, to limit the ratio to the desired range.
        #        -----[SII]-----
        wave = [6718.294, 6732.674]    # vacuum wavelengths
        if not vacuum:
            wave = util.vac_to_air(wave)
        names = ['[SII]6731_d1', '[SII]6731_d2']

        #if FWHM_gal is an array, I need to extract the FWHM values corresponding to the emission lines of the gas template
        if isinstance(FWHM_gal, np.ndarray):
            FWHM_gal_line = np.array([FWHM_gal[find_nearest(wave_galaxy, lw)] for lw in wave])
            gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal_line, pixel)
        else:
            gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel)

        doublets = gauss @ [[0.44, 1.43], [1, 1]]  # produces *two* doublets
        emission_lines = np.column_stack([emission_lines, doublets])
        line_names = np.append(line_names, names)
        line_wave = np.append(line_wave, wave)

    else:

        # Here the two doublets are free to have any ratio
        #         -----[OII]-----     -----[SII]-----
        # wave = [3727.092, 3729.875, 6718.294, 6732.674]  # vacuum wavelengths
        wave = [3727.092, 3729.875, 5198.4, 5201.35, 6718.294, 6732.674]  # vacuum wavelengths with NI "empirical"
        # wave = [3727.092, 3729.875, 5196.45, 5198.94, 6718.294, 6732.674] #right NI wavelengths from the emission file of GIST
        if not vacuum:
            wave = util.vac_to_air(wave)
        # names = ['[OII]3726', '[OII]3729', '[SII]6716', '[SII]6731']
        names = ['[OII]3726', '[OII]3729', '[NI]5196', '[NI]5198', '[SII]6716', '[SII]6731']

        #if FWHM_gal is an array, I need to extract the FWHM values corresponding to the emission lines of the gas template
        if isinstance(FWHM_gal, np.ndarray):
            FWHM_gal_line = np.array([FWHM_gal[find_nearest(wave_galaxy, lw)] for lw in wave])
            gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal_line, pixel)
        else:
            gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel)

        emission_lines = np.column_stack([emission_lines, gauss])
        line_names = np.append(line_names, names)
        line_wave = np.append(line_wave, wave)

    # Here the lines are free to have any ratio
    #       -----[NeIII]-----    HeII      HeI
    wave = [3968.59, 3869.86, 4687.015, 5877.243]  # vacuum wavelengths
    if not vacuum:
        wave = util.vac_to_air(wave)
    names = ['[NeIII]3968', '[NeIII]3869', '-HeII4687-', '-HeI5876-']

    if isinstance(FWHM_gal, np.ndarray):
        FWHM_gal_line = np.array([FWHM_gal[find_nearest(wave_galaxy, lw)] for lw in wave])
        gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal_line, pixel)
    else:
        gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel)

    emission_lines = np.column_stack([emission_lines, gauss])
    line_names = np.append(line_names, names)
    line_wave = np.append(line_wave, wave)

    # NIR H lines
    #       paeps      pad      pab
    wave = [10052.1, 10941.1, 12821.6]  # vacuum wavelengths
    if not vacuum:
        wave = util.vac_to_air(wave)
    names = ['-PaEps-', '-Pad-', '-Pab-']

    if isinstance(FWHM_gal, np.ndarray):
        FWHM_gal_line = np.array([FWHM_gal[find_nearest(wave_galaxy, lw)] for lw in wave])
        gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal_line, pixel)
    else:
        gauss = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel)

    emission_lines = np.column_stack([emission_lines, gauss])
    line_names = np.append(line_names, names)
    line_wave = np.append(line_wave, wave)



    ######### Doublets with fixed ratios #########

    # To keep the flux ratio of a doublet fixed, we place the two lines in a single template
    #        -----[OIII]-----
    wave = [4960.295, 5008.240]    # vacuum wavelengths
    if not vacuum:
        wave = util.vac_to_air(wave)


    if isinstance(FWHM_gal, np.ndarray):
        FWHM_gal_line = np.array([FWHM_gal[find_nearest(wave_galaxy, lw)] for lw in wave])
        doublet = util.gaussian(ln_lam_temp, wave, FWHM_gal_line, pixel) @ [0.33, 1]
    else:
        doublet = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel) @ [0.33, 1]



    # doublet = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel) @ [0.33, 1]
    emission_lines = np.column_stack([emission_lines, doublet])
    line_names = np.append(line_names, '[OIII]5007_d')  # single template for this doublet
    line_wave = np.append(line_wave, wave[1])

    # To keep the flux ratio of a doublet fixed, we place the two lines in a single template
    #        -----[OI]-----
    wave = [6302.040, 6365.535]    # vacuum wavelengths
    if not vacuum:
        wave = util.vac_to_air(wave)



    if isinstance(FWHM_gal, np.ndarray):
        FWHM_gal_line = np.array([FWHM_gal[find_nearest(wave_galaxy, lw)] for lw in wave])
        doublet = util.gaussian(ln_lam_temp, wave, FWHM_gal_line, pixel) @ [1, 0.33]
    else:
        doublet = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel) @ [1, 0.33]



    # doublet = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel) @ [1, 0.33]
    emission_lines = np.column_stack([emission_lines, doublet])
    line_names = np.append(line_names, '[OI]6300_d')  # single template for this doublet
    line_wave = np.append(line_wave, wave[0])

    # To keep the flux ratio of a doublet fixed, we place the two lines in a single template
    #       -----[NII]-----
    wave = [6549.860, 6585.271]    # air wavelengths
    if not vacuum:
        wave = util.vac_to_air(wave)



    if isinstance(FWHM_gal, np.ndarray):
        FWHM_gal_line = np.array([FWHM_gal[find_nearest(wave_galaxy, lw)] for lw in wave])
        doublet = util.gaussian(ln_lam_temp, wave, FWHM_gal_line, pixel) @ [0.33, 1]
    else:
        doublet = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel) @ [0.33, 1]


    # doublet = util.gaussian(ln_lam_temp, wave, FWHM_gal, pixel) @ [0.33, 1]
    emission_lines = np.column_stack([emission_lines, doublet])
    line_names = np.append(line_names, '[NII]6583_d')  # single template for this doublet
    line_wave = np.append(line_wave, wave[1])

    # Only include lines falling within the estimated fitted wavelength range.
    #
    w = (lam_range_gal[0] < line_wave) & (line_wave < lam_range_gal[1])
    emission_lines = emission_lines[:, w]
    line_names = line_names[w]
    line_wave = line_wave[w]

    print('Emission lines included in gas templates:')
    print(line_names)

    return emission_lines, line_names, line_wave
