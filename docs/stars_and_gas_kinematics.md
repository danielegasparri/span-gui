SPAN: SPectral ANalysis software V7.1
Daniele Gasparri, October 2025

# Stars and gas kinematics #

![Kinematics](img/kinematics.png)
*Stars and gas kinematics parameters*


SPAN utilizes the capabilities of the pPXF algorithm (Cappellari 2023 and references therein) to perform a comprehensive kinematic analysis of both stellar and gaseous components in galaxy spectra, covering the optical and near-infrared (NIR) regions. This task can measure the line-of-sight velocity moments both of stars and gas, as well as the fluxes of the gas emission lines automatically found in the spectral region considered.

This interface does not display all the parameters and operations available in pPXF, which is a highly versatile tool for various astrophysical applications. However, it provides a user-friendly setup that will be suitable for most science cases.

The default parameters in this window offer a solid starting point for typical galaxy spectra. To perform an initial fit, simply enter the redshift of the spectrum and ensure that the spectral range falls within the selected wavelength range.
Fine-tuning the fit requires adjusting specific parameters. 


## Parameter Overview ##
The parameters are grouped into seven sections, separated by horizontal dividers. Below is a brief description of each:

### First Section: Basic Parameters
- Wavelength Range: Defines the spectral range to be fitted. If your spectra are not rest-frame corrected, the wavelength range to insert is the one displayed by your spectra. A robust estimation of both stellar and gaseous components is achieved by selecting the visible band region 4800-5500 A (or its equivalent range for high-redshift galaxies), which includes key absorption and emission lines. This is a commonly used standard range. In the NIR, the most prominent spectral features for stellar kinematics are the Ca II triplet (CaT) lines in the 8400-8800 A rest-frame range.
- Sigma (km/s): An initial estimate of the expected velocity dispersion of your galaxy spectrum.
- Redshift (z): An approximate redshift value of the spectrum. **Important:** Avoid loading spectra with significantly different redshifts. If necessary, de-redshift your data before performing kinematic analysis. This can be done using the "Doppler/z correction" task in the Spectra manipulation panel.


### Second Section: Spectral Resolution 
This section requires the spectral resolution of the spectra. This is a fundamental parameter in order to perform a correct estimation of kinematic parameters. It is essential to **avoid mixing spectra with different resolutions within the same dataset.**
The resolution can be expressed in terms of:

- FWHM (Full Width at Half Maximum)
- Resolving Power (R = Lambda/DeltaLambda)
- MUSE LSF. SPAN contains the parametrization of MUSE Line Spread Function based on the Eq.8 of Bacon et al. (2017). If you have MUSE data, please select this option. 
If fitting a narrow spectral region (<= 1000 A), the choice between constant FWHM, R or MUSE resolution is not critical. However, for broader wavelength ranges, selecting the appropriate resolution type is crucial. **Important:** For high-redshift galaxies (z > 0.01) and spectra which have not been de-redshifted, SPAN will automatically correct the spectral resolution given in FWHM and MUSE LSF to the corresponding rest-frame wavelength range.


### Third Section: Template Selection
Here you can select the SSP Model Library to use for the fit. The available pre-loaded libraries are:
	1. E-MILES (subsample)
	2. Galaxev (subsample) 
	3. FSPS (subsample)
	4. X-shooter Spectral Library (XSL) (complete sample with Salpeter IMF, better suited for higher resolution spectra, R = 10,000)  
If your spectra have a higher resolution than the templates, you should degrade them to match the template resolution using the "Degrade Resolution" tool in the Spectra manipulation panel.
You can use any kind of EMILES templates or generic ones. The generic templates should have the wavelength scale in linear units and Angstrom. WARNING: Two stellar component fit is not available with a generic template set option (see next section). 


### Fourth Section: Which component to fit

1. "Fitting only stellar kinematics": Considering only the stellar component. To prevent gas emission lines in the optical, the option "Mask emission lines is activated by default". If you are fitting the optical window and are absolutely sure to not have emission (usually in HBeta and HAlpha) you can disable this option. 
2. "Gas and Stars Kinematics": Considers both stellar and gaseous emission lines. The gas components are divided into three families: Balmer lines, forbidden lines, other lines. SPAN will automatically fit all the components found in the spectral range given, using a modified version of the line emission template provided by pPXF and embedded in the "spectral_analysis.py" module in the "span_functions" folder. By default, stars and gas kinematics are calculated together with one fit. However, often is necessary to fix the stellar kinematics to perform a more accurate gas fit. The option "Fixing stellar kinematics first" will perform a first fit only for the stellar components and will use the kinematics moments derived to fix them in the subsequent fit for the gas component.
3. Two stellar component fit option: Here you can also decide to fit two stellar components by activating the "Fit two stellar components with:" checkbox. This option is important only when you have kinematic decoupled stellar components (e.g. counter-rotating discs) with high S/N (> 50). It can be activated both with stellar and gas fit option, but gas will always be fitted with one kinematic component. Dust correction is not available (both for star and gas) in this mode. You need to insert an accurate velocity and velocity dispersion guess for both the stellar components. Beside this, you can select the templates to assign to each component, based on three options: 1) Old (> 5 Gyr) and young (< 5 Gyr) components, 2) Metal rich ([M/H] > 0 ) and metal poor ([M/H] < 0) components, 3) pick two templates with fixed age and metallicity. This option has been used in some works with counter-rotating discs as such Rubino et al., 2021.
**WARNING:** Two component fit is available only by using the pre-set template set of SPAN and with any (E)MILES template. Not available for fully custom templates, since it needs to read the information about age and metallicity of the SSPs.
**IMPORTANT1:** For the two component fit you must always insert a GOOD guess of the velocity and velocity dispersion of the two components you think to see in your spectrum. I stress out that a good guess is necessary (according to Cappellari et al., 2023 and the pPXF documentation), so first take an accurate look at your spectrum. Some tries may be necessary to find the correct parameters and it is better to perform these tests on spectral regions not affected by gas (e.g. the CaT lines) with the option "Fit only stellar kinematics".
**IMPORTANT2:** If "save processed spectra" is activated, in the two component fit mode SPAN will save the bestfit template from pp.bestfit provided directly by pPXF. However, this bestfit template contains the two components summed up. If you want to separate and analyze each stellar component, SPAN will also save the two besftit templates for each stellar component (suffix comp1_ and comp2_) so that comp1_+comp2_ = pp.bestfit. This uses the pp.matrix keyword of pPXF to reconstruct the bestfit models of each stellar component. Be aware that if you intend to use these separated bestfit templates for further analysis besides the kinematics (e.g. line-strength, SFH) you SHOULD perform the two component fit using ONLY multiplicative degree and deactivating the additive (set to -1). If using additive degree, SPAN will reconstruct and save the two separated bestfit templates without considering the additive polynomials you gave, so you should expect comp1_+comp2_ != pp.bestfit. This should preserve the line fluxes anyway, but it's better to simply avoid the additive degree polynomial, where possible. When gas emission is present, is expected that comp1_+comp2_ != pp.bestfit, since here pp.bestfit contains also the emission lines. 
**IMPORTANT3:** pPXF requires spectra to be normalized to unity to avoid instabilities, but the gas flux calculated on normalized spectra is not physically useful. Therefore, during the "Stars and gas kinematics" task with the option "Stellar and gas kinematics", SPAN saves the normalization factor of each spectrum (np.median(galaxy)) and multiplies by this factor the gas flux calculated by pPXF to give a physical flux for the gas. These gas flux values are NOT shown in the pPXF output in the terminal (I cannot modify the pPXF output), but are saved and stored in the ASCII file generated in the "Process all" mode.


### Fifth Section: Dust and masking
Here you can activate the dust/extinction corrections or decide to mask custom regions of your spectra.

The dust correction uses the Cappellari (2023) 2-parameter attenuation model for the stellar component and the Calzetti (2000) 1-parameter attenuation model for the gas. If you are interested only to kinematics (i.e. no gas flux) you can safely neglect the dust parametrization in most situations.
**WARNING:** Dust correction is not available for the two component fit mode.

The masking will act directly on the "goodpixels" keyword of pPXF and can be used also when the emission line masking is activated.
There are two masking options available:

- Manual masking, by inserting the wavelength interval(s) to be masked, separated by commas and parenthesis in the text box. If invalid values are inserted, SPAN will warn you.
- Graphical masking mode, activated by pressing the "Graphical masking" button. In this mode, an interactive Matplotlib window will open displaying the spectrum selected. You can then mask custom portion directly on the spectrum by ctrl+left click and drag. You can deselect the masked region by ctrl+right click and drag. On touchscreen devices (i.e. Android systems), masking and unmasking modes are activated by a tap on the screen and the relative selection is done by tapping and dragging on the spectrum. When graphical masking in done, you can close the Matplotlib window and the text box on the left will update with the new ranges selected.


### Sixth Section: Kinematic Parameters

- Gauss-Hermite Moments: Determines the complexity of the Line-of-Sight Velocity Distribution (LOSVD) model. Minimum moments to be fitted: 2 (radial velocity + velocity dispersion). Maximum moments: 6. Typical values: 2 for regular galaxies and/or for the two stellar component fit, 4 for interacting or asymmetric galaxies. **IMPORTANT:** moments are always referred to the stellar component. Gas moments are set to 2 by default and cannot be changed.
- Polynomial Degree: Specifies the degree of additive and multiplicative polynomials used by pPXF to adjust the continuum level of spectral templates. Additive polynomials give reliable kinematics results. A degree of 4 is a good starting point. **WARNING:** Additive polynomials are fine for pure kinematic analysis but not for measuring the gas line fluxes. If you are interested also to a correct measurement of gas line fluxes, you should deactivate the additive degree polynomials (set the degree to -1) and use ONLY the multiplicative degree polynomials.
- Noise Level: Represents the expected mean noise level in the spectrum (assumed constant across the wavelength range), used to compute the Chi^2 of the fit and derive formal uncertainties.
  As per the pPXF documentation, formal uncertainties are meaningful only if Chi2 = 1.
  If unsure about the noise level, enable "Auto Noise". This feature will:
  1. Perform an initial fit using the user-provided noise level without regularization (bias keyword set to zero).
  2. Compute a refined noise estimate.
  3. Re-run the fit using the updated noise value to obtain more accurate LOSVD parameters.
- Custom bias keyword: To be used only when you fit at least 4 moments. This lets you to change the default bias value that pPXF uses (bias = 0.7*np.sqrt(500./goodpixels.size), from Cappellari & Emsellem 2004) and try to adjust the fit to a better estimate of H3 and H4 moments.

### Seventh Section: Uncertainty Estimation
To estimate uncertainties in stellar kinematics, you can enable Monte Carlo simulations.
This option is recommended if you do not fully trust the formal uncertainties computed by pPXF.
Particularly useful for galaxies with very low velocity dispersion compared to the instrumental resolution (velscale). This operation is very time consuming, so try to not insert too many simulations. 


## Outputs ##
In **"Process selected"** mode, if the option "Save processed spectra" is activated, the task produces:

- bestfit model spectrum
- bestfit gas model spectrum (if kinematics of stars and gas is activated)
- emission corrected spectrum (if kinematics of stars and gas is activated)
- bestfit model spectrum of the stellar component 1 and 2 (if "Fit two stellar components" option is activated)
- continuum subtracted gas spectrum (if kinematics of stars and gas is activated).  
These spectra are saved in the "processed_spectra" subfolder within the "SPAN_results" folder.
If you are not interested in these spectra products, you can just disable the option "Save processed spectra" and save some space on the disk.

In **"Process all"** mode, for each spectrum, if "Save processed spectra" is activated, the task produces all the spectra in the "Process selected" mode, for each spectrum processed.

Moreover, an ASCII file (.dat) is stored in the "stars_and_gas_kinematics" subfolder containing the kinematics moments and the formal errors of the stellar component(s). If the option "Fitting gas and stellar kinematics together" is activated, an ASCII file containing the kinematics, flux and formal errors of the emission lines found is generated.

If "Estimate the uncertainties with MonteCarlo simulations" is activated, another ASCII file is generated containing the stellar kinematics and the MonteCarlo uncertainties. MonteCarlo errors are not calculated for gas.
