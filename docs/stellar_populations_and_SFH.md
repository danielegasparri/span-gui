SPAN: SPectral ANalysis software V6.6
Daniele Gasparri, July 2025

# Stellar populations and SFH #

![Stellar populations](img/stellar_populations.png)
*Stellar populations and SFH parameters*


This module is one of the most critical tools within the spectral analysis framework and extragalactic astronomy in general. It leverages the well-known pPXF algorithm (Cappellari 2023 and referenes therein) to fit a galaxy's spectrum using a set of Single Stellar Population (SSP) models, extracting key stellar population parameters such as age, metallicity, Star Formation History (SFH) and mass growth.

The provided settings allow for an optimal use of pPXF, though not all of its parameters and functionalities are included in this interface.


**Requirements for Reliable Results**  
The full spectral fitting method requires spectra with an acceptable Signal-to-Noise Ratio (S/N):
- S/N > 15 for basic results.
- S/N > 50 for high-accuracy results.


The default settings offer a good starting point for most local galaxy spectra and common spectrographs. To perform an initial fit, simply enter the redshift of the spectrum. For fine-tuning, adjust the parameters based on your specific dataset.


## Parameter Overview ##
The parameters are grouped into six sections, separated by horizontal dividers. Below is a brief description of each:


### First Section: Basic Parameters
- Wavelength Range: Defines the spectral region to be fitted. A good estimate of stellar population parameters in galaxies is obtained by focusing on the 4800-5500 A (rest-frame) which includes age, metallicity, and SFH-sensitive lines such as Hbeta, OIII, Mg, and Fe. Il you have near-infrared (NIR) spectra, the 8400-8800 A (rest-frame) range is good, which includes the Ca II triplet.
- Spectral Resolution FWHM (A): Approximate resolution of your spectrum **in the selected wavelength range**. This is relevant only if including gas emission lines in the fit or if the option "Convolve templates to galaxy resolution" is activated (not mandatory). A precise value is not required, but a good estimate is recommended.
- Velocity Dispersion Guess (km/s): An approximate estimate of the actual stellar velocity dispersion.
- Redshift Guess (z): An initial estimate of the spectrum's redshift. If this value is not within +/- 0.005 the real one, the fit will likely fail.


### Second Section: Gas Emission Lines and Dust Attenuation
- If the spectrum contains emission lines, you should select "Fitting with gas" for more accurate results. Alternatively, you can mask out the emission lines with the "Mask emission lines" option and perform a fit without considering the gas ("Fitting only stars")
- If the spectrum includes multiple Balmer emission lines, you may enable "Tie Balmer" to constrain their ratios. This automatically applies the Calzetti et al. (2000) dust attenuation curve.
Two additional options allow for dust correction. They should used only for spectra with very good flux calibration:
- "Correct for dust the stars": Uses the Cappellari (2023) 2-parameter attenuation model.
- "Correct for dust the gas": Uses the Calzetti (2000) 1-parameter attenuation model.



### Third Section: Noise Estimation, Regularization and polynomials
Noise and Regularization: These are critical for obtaining accurate results. Refer to the pPXF documentation for an in-depth guide.
Recommended Workflow for finding the best parameters from pPXF documentation:

- Perform an unregularized fit (Regul. error = 0).
- Adjust the noise level to achieve chi^2 = 1.
- Set the regularization so that the current delta Chi^2 falls between 20-50% of the desired delta Chi^2.
- If uncertain, use "Auto Noise" option, which performs an initial unregularized fit to determine the optimal noise level.
- For small spectral ranges (<1000 A) and limited templates, "Auto Noise and Regul. error" can be used to optimize both parameters. A good starting point for S/N = 50 is setting "Fraction of Dchi2 to reach" to 0.20.

**Polynomial Adjustments**:  
Additive Degree: Leave disabled (-1) for reliable results, since they may change the absorption line profiles.
Multiplicative Degree: Adjust based on spectral range (about 1 degree per 100 A).
Example: For 4800-5500 A, set Mult. degree = 7.


### Fourth Section: Template Library Selection
Choose the SSP model library for fitting:

- E-MILES, Galaxev, FSPS (pPXF defaults, they are a subsample of the respective libraries).
- X-shooter Spectral Library. A full sample in the optical range with Salpeter IMF and Parsec/Colibri isochrones. They have a resolving power of R = 10,000 and are better suited for medium to high resolutio galaxy spectra.
- sMILES Library (Knowles et al., 2023). Full sample in the optical window with 4 alpha/Fe values and Salpeter IMF. They are useful to extract also the alpha/Fe values besides age and metallicity. 

**Warning** If your spectra have much higher resolution than the templates, apply the "Degrade Resolution" task in the Spectra manipulation panel before fitting may improve the results.


**Template Handling:**  
The sMILES templates are stored in spectralTemplates/sMILES_afeh. This folder can hosts any kind of sMILES templates, not just those provided with SPAN.
If replacing templates, ensure:

- A regular grid of ages and metallicities.
- No mixing of different IMF templates.

You can also insert a custom .npz file following the standard of pPXF and containing any set of templates.

**Important Notes on Stellar Population Models**  
E-MILES, sMILES, and X-shooter work best for quiescent and moderately star forming galaxies but lack very young stellar populations (<30-60 Myr).
For highly star-forming galaxies, consider FSPS (Conroy et al., 2010) or other templates with (also) very young ages. If you use a set of templates that does not properly sample the needed very young regime, you will see a very young, very metal poor population appearing in your fit. That is a red flag!


### Fifth Section: Custom Masking and Stellar Constraints

- You can mask out the emission Lines. In this case, you should select the "Fitting without gas" in Section 3.
- Custom masking is available and can be added also to the automatic masking of the emission lines. There are two custom masking options available: 
    1. Manual masking, by inserting the wavelength interval(s) to be masked, separated by commas and parenthesis in the text box. If invalid values are inserted, SPAN will warn you.
    2. Graphical masking mode, activated by pressing the "Graphical masking" button. In this mode, an interactive Matplotlib window will open displaying the spectrum selected. You can then mask custom portion directly on the spectrum by ctrl+left click and drag. You can deselect the masked region by ctrl+right click and drag. On touchscreen devices (i.e. Android systems), masking and unmasking modes are activated by a tap on the screen and the relative selection is done by tapping and dragging on the spectrum. When graphical masking in done, you can close the Matplotlib window and the text box on the left will update with the new ranges selected.

**Age and Metallicity Ranges:**  
Limit the maximum template age and metallicity based on the galaxy's cosmic age. For high-redshift galaxies, reducing the upper age limit is strongly recommended. Example: we are analyzing a galaxy located at z = 1.20, corresponding to an epoch when the Universe was about 5 Gyr old. You can move the maximum age slider to 5 Gyr to consider only SSP younger than this value for the fit. 


### Sixth Section: Uncertainty Estimation and Lick Indices
Bootstrap Uncertainty Estimation (Kacharov et al., 2018) can be enabled to compute error on age and metallicity.
Suggested simulations: 20-50 (balance between accuracy and speed).

**Lick/IDS Index Analysis:**  
Uses the pPXF emission-corrected spectrum to estimate the stellar parameters also with the Lick/IDS indices.
Select model grids and interpolation mode as in the "Line-Strength Analysis" task.


## Final Notes ##
By default, templates are V-band normalized (5070-5950 A), meaning age, metallicity, and SFH are V luminosity-weighted.
Using pPXF's .flux attribute, SPAN derives mass-weighted results from the same fit.
Best regularization differs for mass-weighted vs. luminosity-weighted results. If mass results are the priority, fine-tune the "Regul. error" accordingly.


## Outputs ##
Every fit generates three plots both in the **"Process selected** and **Process all"** mode:

- Fitted Spectrum
- Age and Metallicity Distribution (luminosity-weighted and mass-weighted).
- Non-parametric SFH (luminosity-weighted and mass-weighted). Carefully analyze these plots to determine whether the fit is physically meaningful. No computer can yet decide if a fit, even if statistically excellent, makes sense in the real universe.
- An ASCII file (.dat) containing the luminosity and mass fraction and cumulative as a function of the age (mass growth). If the "Estimate the uncertainties for age and met (long process)" is activated, also the relative uncertainties are added to the file. This file is stored in the "stellar_population_and_sfh/SFH/" subfolder.
- Two ASCII files containing the luminosity and mass weights respectively. These files are stored in the "stellar_population_and_sfh/weights/" subfolder
- FITS files containing: 1) residual of the fit, 2) bestfit template, 3) emission corrected spectrum (if gas is considered), stored in the "processed_spectra" subfolder.

In **"Process all"** mode, the task stores all the files of the "Process selected" mode for each spectrum. In addition:

- An ASCII file with kinematics, mean luminosity and mass age, metallicity and alpha/Fe (if available) and the relative uncertainties if "Estimate the uncertainties for age and met (long process)" if activated, stored in the "stellar_population_and_sfh" subfolder.
