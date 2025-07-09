SPAN: SPectral ANalysis software V6.6
Daniele Gasparri, July 2025

# The Spectral analysis panel#

This panel consists of a single frame containing both basic and advanced spectral analysis tools:

1. Blackbody fitting
2. Cross-correlation
3. Velocity dispersion measurement
4. Line-strength analysis
5. Line fitting
6. Stars and gas kinematics
7. Stellar populations and SFH

Each task operates independently and does not modify the input spectrum.


Here is a brief description of the basic spectral analysis task. For line-strength analysis, Stars and gas kinematics and Stellar populations and SFH, see the dedicated documentation. 


**Planck blackbody fitting**  
Fits the spectral continuum with a Planck blackbody function and returns the best estimation of the effective temperature. The fit is accurate for stellar spectra spanning a wide wavelength range (∼ 3000 Å) and works at best if the peak of the Plank function is included in the spectral region to fit. The results are the best estimation of the effective temperature of the considered stellar spectrum and the best fit Planck function model.

**Cross-Correlation**  
Measures the wavelength shift of the spectra, both in terms of radial velocity and redshift values. 
The result is the radial velocity or z of the spectrum, and the uncertainties estimation via Monte
Carlo simulations.

**Velocity dispersion**  
Performs a brute force least square fit of a selected region of a spectrum with a selected template spectrum and calculates the velocity dispersion. This task produces faster results than the pPXF algorithm and it is suitable in situations where speed is essential. It can be used for an estimation of the velocity dispersion for any situation where an accurate kinematic analysis is not necessary. Spectra needs to de de-redshifter and Doppler corrected.


**Line(s) fitting**  
Fits an absorption or an emission line with a convolution of a straight line for the continuum and a Gaussian function for the line. It can also automatically fit the three Calcium Triplet (CaT) lines in the NIR and gives the EW measured from the best-fitting model, providing de-redshifted and Doppler corrected spectra.


**Working with the Input Spectrum**  
WARNING: The input spectrum may be affected by the tasks activated in the Spectra manipulation panel.
If you wish to analyze the original spectrum, ensure that all the Spectra manipulation tasks are deactivated (they are disabled by default).


**Previewing Results**  
Click "Preview result" to display the analysis output in the output frame and as a Matplotlib plot.
Once satisfied with the preview, you can proceed with processing either a single spectrum or all spectra in the list.


**Processing the Spectra**  

- Process a Single Spectrum:
Click "Process selected" to apply all activated tasks to the currently selected spectrum. The results will be displayed in the output frame.
- Process All Spectra:
Click "Process all" to apply all activated tasks to every spectrum in your list (if using multiple spectra).
Important: This is the only way to save spectral analysis results to a file.


**Saving Plots**  
To save plots generated during the analysis in the "Process all" mode, enable the "Save plots" option.
The program will automatically save a PNG file for each spectrum and each selected analysis task.
Saved plots will be stored in the "SPAN_results/plots" folder.
