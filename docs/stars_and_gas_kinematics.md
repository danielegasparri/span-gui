SPAN: SPectral ANalysis software V6.6
Daniele Gasparri, June 2025

### Stars and gas kinematics ###

SPAN utilizes the capabilities of the pPXF algorithm (Cappellari 2023) to perform a comprehensive kinematic analysis of both stellar and gaseous components in galaxy spectra, covering the optical and near-infrared (NIR) regions.

This interface does not display all the parameters and operations available in pPXF, which is a highly versatile tool for various astrophysical applications. However, it provides a user-friendly setup that will be suitable for most spectral analyses.

The default parameters in this window offer a solid starting point for typical galaxy spectra and state-of-the-art spectrographs. To perform an initial fit, simply enter the redshift of the spectrum and ensure that the spectral range falls within the selected wavelength range.
Fine-tuning the fit requires adjusting specific parameters. 


Parameter Overview
The parameters are grouped into seven sections, separated by horizontal dividers. Below is a brief description of each:

First Section: Basic Parameters
- Wavelength Range: Defines the spectral range to be fitted. A robust estimation of both stellar and gaseous components is achieved by selecting the visible band region 4800-5500 A (or its equivalent range for high-redshift galaxies), which includes key absorption and emission lines. This is a commonly used standard range. In the NIR, the most prominent spectral features for stellar kinematics are the Ca II triplet (CaT) lines in the 8400-8800 A rest-frame range, while for gas emission, the Paschen lines (especially PaB) should be considered.
- Sigma (km/s): An initial estimate of the expected velocity dispersion of your galaxy spectrum.
- Redshift (z): An approximate redshift value of the spectrum. Important: Avoid loading spectra with significantly different redshifts. If necessary, de-redshift your data before performing kinematic analysis. This can be done using the Doppler/z correction task in SPAN's Spectra manipulation panel.


Second Section: Spectral Resolution
This section requires the spectral resolution of the spectra. It is essential to avoid mixing spectra with different resolutions within the same dataset.
Additionally, specify whether the resolution is expressed in terms of:
- FWHM (Full Width at Half Maximum)
- Resolving Power (R = Lambda/DeltaLambda)
- MUSE LSF
If fitting a narrow spectral region (<= 1000 A), the choice between constant FWHM, R or MUSE resolution is not critical. However, for broader wavelength ranges, selecting the appropriate resolution type is crucial. For high-redshift galaxies (z > 0.01), SPAN will automatically correct the spectral resolution to match the rest-frame wavelength range.


Third Section: Template Selection
- SSP Model Library: Choose the template library for the fit. The available pre-loaded libraries are:
	1) E-MILES (subsample)
	2) Galaxev (subsample)
	3) FSPS (subsample)
	4) X-shooter Spectral Library (XSL) (complete sample with Salpeter IMF, better suited for higher resolution spectra, R = 10,000)
  If your spectra have a higher resolution than the templates, you should degrade them to match the template resolution using the "Degrade Resolution" tool in the Spectral Processing module.
  You can use any kind of EMILES templates or even generic templates. These should have the wavelength scale in linear units and Angstrom. 


Fourth Section: Which component to fit
- Fitting Mode:
	1) "Fitting only stellar kinematics": You can mask all potential emission lines to fit only the stellar component. Here you can also decide whether fitting two stellar components by activating the "Fit two stellar components with the following parameters" checkbox. The two component fit is performed by extracting two SSP templates from the library you have chosen above, with defined age and metallicity, following Rubino et al., 2021. You must insert the age and metallicity values of the two SSP to be retrieved, as well as a GOOD guess of the velocity and velocity dispersion of the two components you think to see in your spectrum. I stress out that a good guess is necessary (according to Cappellari et al., 2023 and the pPXF documentation), so first take an accurate look at your spectrum and make some tries.
 
	2) "Gas and Stars Kinematics": Considers both stellar and gaseous emission lines. Here only one stellar component is fitted, along with as many gaseous components SPAN will find on your spectrum (Balmer liner, forbidden lines and other lines). By default, stars and gas are calculated within the same fit. However, often is necessary to fix the stellar kinematics to perform a more accurate gas fit. The option 'Fixing stellar kinematics first' will perform a first fit only for the stellar components and will use the kinematics moments derived to fix them in the subsequent fit for the gas component. 


Fifth Section: Dust and masking
Here you can activate the dust/extinction corrections or decide to mask custom regions of your spectra. This masking will act directly on the 'goodpixels' keyword of pPXF and can be used also when the emission line masking is activated. 
There are two masking options available: the manual ones, by inserting the wavelength interval in the text box or a graphical mode activated by pressing the 'Graphical masking' button. In this mode, an interactive Matplotlib window will open displaying the spectrum selected. You can then mask custom portion directly on the spectrum by ctrl+left click and drag. You can deselect the masked region by ctrl+right click and drag. On touchscreen devices (i.e. Android systems), masking and unmasking modes are activated by a tap on the screen and the relative selection is done by tapping and dragging on the spectrum. When graphical masking in done, you can close the Matplotlib window and the text box will update with the new ranges selected. 


Sixth Section: Kinematic Parameters
- Gauss-Hermite Moments: Determines the complexity of the Line-of-Sight Velocity Distribution (LOSVD) model.
  Minimum: 2 (radial velocity + velocity dispersion).
  Maximum: 6.
  Typical values: 2 for regular galaxies and/or for the two stellar component fit, 4 for interacting or asymmetric galaxies.

- Polynomial Degree: Specifies the degree of additive and multiplicative polynomials used by pPXF to adjust the continuum level of spectral templates. Additive polynomials give reliable kinematics results. A degree of 4 is a good starting point. 
WARNING: Additive polynomials are fine for pure kinematics analysis but not for measuring the gas line fluxes. If you are interested also to a correct measurement of gas line fluxes, you should deactivate the additive degree polynomials (set the degree to -1) and use ONLY the multiplicative degree polynomials. 

- Noise Level: Represents the expected mean noise level in the spectrum (assumed constant across the wavelength range), used to compute the Chi^2 of the fit and derive formal uncertainties.
As per the pPXF documentation (pPXF Documentation), formal uncertainties are meaningful only if Chi2 = 1.
If unsure about the noise level, enable "Auto Noise". This feature will:
	1) Perform an initial fit using the user-provided noise level without regularization (bias keyword set to zero).
	2) Compute a refined noise estimate.
	3) Re-run the fit using the updated noise value to obtain more accurate LOSVD parameters.


Seventh Section: Uncertainty Estimation
To estimate uncertainties in stellar kinematics, you can enable Monte Carlo simulations.
This option is recommended if you do not fully trust the formal uncertainties computed by pPXF.
Particularly useful for galaxies with very low velocity dispersion compared to the instrumental resolution. This operation is very time consuming, so try to not insert too many simulations. 


Outputs:
In 'Process selected' mode, if the option 'Save processed spectra' is activated, the task produces:
	- bestfit model spectrum
	- bestfit gas model spectrum (if kinematics of stars and gas is activated)
	- emission corrected spectrum (if kinematics of stars and gas is activated)
	- continuum subtracted gas spectrum (if kinematics of stars and gas is activated)
These spectra are saved in the 'processed_spectra' subfolder within the 'SPAN_results' folder. 
If you are not interested in these spectra products, you can just disable the option 'Save processed spectra' and save some space on the disc. 

In 'Process all' mode, for each spectrum, if 'Save processed spectra' is activated, the task produces all the spectra in the 'Process selected' mode, for each spectrum processed. 

Moreover, an ASCII file (.dat) is stored in the 'stars_and_gas_kinematics' subfolder containing the kinematics moments and the formal errors of the stellar component(s). If the option 'Fitting gas and stellar kinematics together' is activated, an ASCII file containing the kinematics, flux and formal errors of the emission lines found is generated. 
If 'Estimate the uncertainties with MonteCarlo simulations' is activated, another ASCII file is generated containing the stellar kinematics and the MonteCarlo uncertainties. MonteCarlo errors are not calculated for gas. 