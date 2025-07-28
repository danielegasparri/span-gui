SPAN: SPectral ANalysis software V6.6
Daniele Gasparri, July 2025

# Line-strength analysis #

![Linestrength](img/linestrength.png)
*Line-strength analysis parameters*


SPAN provides a complete and customizable line-strength analysis framework for both stellar and galaxy spectra using built-in routines. These routines have already been used in published works (e.g., Gasparri et al. 2021, 2024) and thoroughly tested. They measure the equivalent width (EW) of absorption and emission lines and their uncertainties with Monte Carlo simulations. 

You can perform two types of line-strength analysis. Let’s take a closer look at both.


## Fully custom analysis using any line-strength index definition
This mode includes the **"User indices on a list file"** and **"Single index"** options. Your spectra **must already be redshift/Doppler corrected** before using them. The wavelength range and the resolution values should be set to the rest-frame reference. 
For **"User indices on a list file"**, you must create and load into SPAN a valid ASCII file containing the names and definitions of the indices you want to measure. An example file is included, and additional index templates can be found in the "example_files" folder. More information can be found in the user manual. For people in a hurry, here some tips:

- The ASCII file can be in .txt or .dat format, with single spaced values.
- The first row of the file must contain the arbitrary names of the indices.
- For any index, starting from the second row, you must provide in the following six rows: blue band limits (min, max), red band limits (min, max), line limits (min, max). Don't worry if the file is not properly formatted: SPAN will notice it and will complain, but it will not crash. 


### From raw EWs to science values
If you measure the EW of your indices on stellar or model spectra, you are set. But if you perform the measurements on galaxy spectra, this is just the first step. You have the raw EW, not necessarily the final, scientific values. The velocity dispersion of the galaxies, indeed, broadends the spectral lines and has effect on the EWs measured. If you are measuring the indices on a single galaxy with constant velocity dispersion, you maybe are ok, but this is a very specific situation that happens few times in the life of an astronomer. In most situations, to assign meaningful values to your measured raw EWs, you should account for the effect of the velocity dispersion broadening and therefore correct the raw values. 

In the custom line-strength analysis you should take care personally of these two steps, but SPAN will help you. Following the approach of Trager et al., 1998 on the Lick/IDS indices, which has been replied (with some modifications) also in the near-infrared (NIR, e.g. Cesetti et al., 2013, Morelli et al., 2020, Gasparri et al., 2021), you can find two dedicated tasks for this:

1. **Calculate velocity dispersion coefficients:** By loading a set of stellar or model templates to SPAN and the same list of indices you measured on your spectra, SPAN will calculate the corrections and the uncertainties to apply for the raw EW od any index as a function of the velocity dispersion broadening. This is achieved by broadening each template spectrum to different velocity dispersion values, from 0 to 400 km/s in 50 km/s steps, and then fitting a spline function to retrieve the relation between the EW and velocity dispersion, for each index. 
2. **Correct line-strength for velocity dispersion:** Once you have the coefficients, this task will help you to correct the raw EWs generated in the first step for the velocity dispersion broadening. You should also know the velocity dispersion of each spectrum. You can measure it with the "Velocity dispersion" or "Stars and gas kinematics" tasks. The "Velocity dispersion" tasks generates ASCII files ready for the "Correct line-strength for velocity dispersion" task. For the "Stars and gas kinematics" output, you just need to delete the RV row. See more details about how the required text files must be formatted in the user manual.


## Automatic Lick/IDS index analysis
In this case, SPAN can do the dirty job for you. Activating all the options available, first a fit with pPXF is performed in order to measure the velocity dispersion and the redshift/Doppler velocity of any single spectrum, considering also the gas emission lines and subtracting them from the spectra. For high redshift galaxies, you do NOT need to de-redshift them. If the corresponding rest-frame wavelength range falls within the Lick/IDS index range, SPAN will notice it, if you inserted the right redshift value and activated the "Perform Doppler correction" option. In this case, the resolution is referred to the wavelength range of your redshifted spectra. If you are working with high redshift spectra (z > 0.01) SPAN will automatically adjust the resolution (if expressed in terms of FWHM) to the equivalent rest-frame wavelength range of the Lick/IDS indices. The definitions of the Lick/IDS indices are stored in the "system_files" folder of SPAN and automatically loaded. Also the resolution adjustement to bring your spectra to the Lick/IDS standard is performed by SPAN. 

Since some of the Lick/IDS indices are used to constrain the properties of the stellar populations in unresolved spectra, you can activate the "Estimate stellar parameters with SSP models" and select which pre-loaded model use and the algorithm to perform the interpolation between data and models. 
The interpolation with the model-based Lick/IDS indices is performed considering the Hbeta–MgFe' indices of Thomas et al., 2003 and the Fe–Mgb index-index grids.  

The interpolation between the measured and model indices is performed with two different approaches:

1. **Linear interpolation**  
   A multi-dimensional linear interpolation using the "griddata" function from the "SciPy" module retrieves age, metallicity, and alpha/Fe (where available) from the observed Lick/IDS indices projected onto the model grids.
2. **Supervised Machine Learning**  
   The **Gaussian Process Regression (GPR)** has proven to be a powerful tool for deriving stellar parameters from stellar spectra (e.g. Bu et al., 2015,2020). To our knowledge, no attempt has yet been made to apply GPR to unresolved galaxy spectra in the context of Lick/IDS line-strength studies. In SPAN, I trained the GPR algorithm with the line-strength predictions of the SSP models. The trained models, stored in the "system_files" folder, are used to constrain the age, metallicity, alpha/Fe and the relative uncertainties of the stellar populations. The GPR method results to be 10 times faster and yields results in excellent agreement with the linear interpolation method for spectra with S/N > 20. For lower S/N spectra, the GPR still converges whereas the simple linear interpolation fails to find a solution or gives very large uncertainties.  


## Outputs
Using the custom line-strength analysis in the "Process all" mode, SPAN will save three ASCII files in the "line-strength_analysis" folder contained in the "SPAN_results" folder:

- A file with the raw EWs and uncertainties expressed in Angstrom.
- A file with the raw EWs and uncertainties expressed in Magnitudes.
- A file with the S/N per pixel calculated in the pseudo-continuum bands for each index.

Using the automatic Lick/IDS analysis in the "Process all" mode, SPAN will save the same ASCII files of the custom analysis, which contains the corrected EWs for the velocity dispersion broadening, if this option ("correct for sigma") has been activated. In addition:

- An ASCII file containing the EWs of the Lick/IDS indices used for stellar population studies (Hbeta, Mg2 (in mag), Mgb, Fe, MgFe') and their uncertainties.
- An ASCII file containing the EWs of these indices, and the stellar parameters estimated (age, metallicity and alpha/Fe) and their uncertainties, if "Estimate stellar parameters with SSP models" option has been activated. 
