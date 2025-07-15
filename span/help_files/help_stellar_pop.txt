SPAN: SPectral ANalysis software V6.6
Daniele Gasparri, June 2025

###Stellar populations and SFH ###

This module is one of the most critical tools within the spectral analysis framework and extragalactic astronomy in general. It leverages the well-known pPXF algorithm (Cappellari 2023) to fit a galaxy's spectrum using a set of Single Stellar Population (SSP) models, extracting key stellar population parameters such as age, metallicity, Star Formation History (SFH) and mass growth.

The provided settings allow for an optimal use of pPXF, though not all of its parameters and functionalities are included in this interface.


Requirements for Reliable Results
The full spectral fitting method requires spectra with an acceptable Signal-to-Noise Ratio (S/N):
- S/N > 15 for basic results.
- S/N > 50 for high-accuracy results.


The default settings offer a good starting point for most local galaxy spectra and common spectrographs. To perform an initial fit, simply enter the redshift of the spectrum. For fine-tuning, adjust the parameters based on your specific dataset.


Parameter Overview
The parameters are grouped into six sections, separated by horizontal dividers. Below is a brief description of each:


First Section: Basic Parameters
- Wavelength Range: Defines the spectral region to be fitted. A good estimate of stellar population parameters in galaxies is obtained by focusing on:
4800-5500 A (rest-frame) in the optical, covering age, metallicity, and SFH-sensitive lines such as Hbeta, OIII, Mg, and Fe.
8400-8800 A (rest-frame) in the near-infrared, including the Ca II triplet, which is useful for stellar kinematics.
- Spectral Resolution FWHM (A): Approximate resolution of your spectrum in the selected range. This is relevant only if including gas emission lines in the fit. A precise value is not required, but a good estimate is recommended.
- Velocity Dispersion Guess (km/s): An approximate estimate of the actual stellar velocity dispersion.
- Redshift Guess (z): An initial estimate of the spectrum's redshift. If this value is incorrect, the fit will likely fail.


Second Section: Gas Emission Lines and Dust Attenuation
- If the spectrum contains emission lines, select "Fitting with gas" for improved results.
- If the spectrum includes multiple Balmer emission lines, enable "Tie Balmer" to constrain their ratios. This automatically applies the Calzetti et al. (2000) dust attenuation curve.
Two additional options allow for dust correction:
- "Correct for dust (stars)": Uses the Cappellari (2023) 2-parameter attenuation model.
- "Correct for dust (gas)": Uses the Calzetti (2000) 1-parameter attenuation model.


Third Section: Noise Estimation, Regularization and polynomials
Noise and Regularization: These are critical for obtaining reliable results. Refer to the pPXF documentation for an in-depth guide.
Recommended Workflow:
- Perform an unregularized fit (Regul. error = 0).
- Adjust the noise level to achieve chi^2 = 1.
- Set the regularization so that the current delta Chi^2 falls between 20-50% of the desired delta Chi^2.
- If uncertain, use "Auto Noise", which performs an initial unregularized fit to determine the optimal noise level.
- For small spectral ranges (<1000 A) and limited templates, "Auto Noise and Regul. error" can be used to optimize both parameters. A good starting point for S/N = 50 is setting "Fraction of Dchi2 to reach" to 0.20.

- Polynomial Adjustments:
Additive Degree: Leave disabled (-1) for reliable results.
Multiplicative Degree: Adjust based on spectral range (about 1 degree per 100 A).
Example: For 4800-5500 A, set Mult. degree = 7.


Fourth Section: Template Library Selection
Choose the SSP model library for fitting:
- E-MILES, Galaxev, FSPS (pPXF defaults).
- X-shooter Spectral Library (R = 10,000) (high-resolution spectra).
- sMILES Library (Knowles et al., 2023) (4 alpha/Fe values, Salpeter IMF).
If using sMILES, SPAN will also extract alpha/Fe values.
If your spectra have higher resolution than the templates, apply the "Degrade Resolution" task in the Spectra manipulation panel before fitting.

Template Handling:
The sMILES templates are stored in spectralTemplates/sMILES_afeh. This folder can hosts any kind of sMILES templates, not just those provided with SPAN.
If replacing templates, ensure:
- A regular grid of ages and metallicities.
- No mixing of different IMF templates.

Important Notes on Stellar Population Models
E-MILES, sMILES, and X-shooter work best for quiescent galaxies but lack very young stellar populations (<50-60 Myr).
For star-forming galaxies, consider FSPS (Conroy et al., 2010).


Fifth Section: Custom Masking and Stellar Constraints
You can mask out the emission Lines. In this case, you should select the "Fitting without gas" in Section 3.
Custom masking is available and can be added also to the automatic masking of the emission lines. 
There are two masking options available: the manual ones, by inserting the wavelength interval in the text box or a graphical mode activated by pressing the 'Graphical masking' button. In this mode, an interactive Matplotlib window will open displaying the spectrum selected. You can then mask custom portion directly on the spectrum by ctrl+left click and drag. You can deselect the masked region by ctrl+right click and drag. On touchscreen devices (i.e. Android systems), masking and unmasking modes are activated by a tap on the screen and the relative selection is done by tapping and dragging on the spectrum. When graphical masking in done, you can close the Matplotlib window and the text box will update with the new ranges selected. 

Age and Metallicity Ranges:
Limit the maximum template age based on the galaxy's cosmic age. For high-redshift galaxies, reducing the upper age limit is recommended.


Sixth Section: Uncertainty Estimation and Lick Indices
Bootstrap Uncertainty Estimation (Kacharov et al., 2018) can be enabled to compute error on age and metallicity.
Suggested simulations: 20-50 (balance between accuracy and speed).

Lick/IDS Index Analysis:
Uses the pPXF emission-corrected spectrum to estimate the stellar parameters also with the Lick/IDS indices.
Select model grids and interpolation mode as in the Line-Strength Analysis module.


Final Notes
By default, templates are V-band normalized (5070-5950 A), meaning age, metallicity, and SFH are luminosity-weighted.
Using pPXF's .flux attribute, SPAN derives mass-weighted results from the same fit.
Best regularization differs for mass-weighted vs. luminosity-weighted results. If mass results are the priority, fine-tune the "Regul. error" accordingly.


Task Output
Every fit generates three plots:
- Fitted Spectrum
- Age and Metallicity Distribution (luminosity-weighted and mass-weighted).
- Non-parametric SFH (luminosity-weighted and mass-weighted).

Carefully analyse these plots to determine whether the fit is physically meaningful. No computer can yet decide if a fit, even if statistically excellent, makes sense in the real universe.


Outputs:
In 'Process selected mode' the task produces for the processed spectrum:
	- An ASCII file (.dat) containing the luminosity and mass fraction and cumulative as a function of the age (mass growth). If the 'Estimate the uncertainties for age and met (long process)' is activated, also the relative uncertainties are added to the file. This file is stored in the 'stellar_population_and_sfh/SFH/' subfolder.  
	- Two ASCII files containing the luminosity and mass weights respectively. These files are stored in the 'stellar_population_and_sfh/weights/' subfolder
	- FITS files containing: 1) residual of the fit, 2) bestfit template, 3) emission corrected spectrum (if gas is considered), stored in the 'processed_spectra' subfolder. 

In 'Process all mode', the task stores all the files of the 'Process selected' mode for each spectrum. Moreover:
	- An ASCII file with kinematics, mean luminosity and mass age, metallicity and alpha/Fe (if available) and the relative uncertainties if 'Estimate the uncertainties for age and met (long process)' if activated, stored in the 'stellar_population_and_sfh' subfolder. 