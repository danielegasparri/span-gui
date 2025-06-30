SPAN: SPectral ANalysis software V6.6
Daniele Gasparri, June 2025

###The Spectral analysis panel### 

This panel consists of a single frame containing both basic and advanced spectral analysis tools:
1) Blackbody fitting
2) Cross-correlation
3) Velocity dispersion measurement
4) Line-strength analysis
5) Line fitting
6) Stars and gas kinematics
7) Stellar populations and Star Formation History (SFH)

Each task operates independently and does not modify the input spectrum.


Working with the Input Spectrum
WARNING: The input spectrum may be affected by the tasks activated in the Spectra manipulation panel.
If you wish to analyze the original spectrum, ensure that all the Spectra manipulation tasks are deactivated (they are disabled by default).


Previewing Results
Click "Preview result" to display the analysis output in the output frame and as a Matplotlib plot.
Once satisfied with the preview, you can proceed with processing either a single spectrum or all spectra in the list.


Processing the Spectra
- Process a Single Spectrum:
Click "Process selected" to apply all activated tasks to the currently selected spectrum. The results will be displayed in the output frame.

- Process All Spectra:
Click "Process all" to apply all activated tasks to every spectrum in your list (if using multiple spectra).
Important: This is the only way to save spectral analysis results to a file.


Saving Plots
To save plots generated during the analysis in the "Process all" mode:

Enable the "Save plots" option.
The program will automatically save a PNG file for each spectrum and each selected analysis task.
Saved plots will be stored in the "SPAN_results/plots" folder.
