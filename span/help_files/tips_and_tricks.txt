SPAN: SPectral ANalysis software V6.6
Daniele Gasparri, June 2025


Here are some useful tips to help you make the most of SPAN.

1. Setting the Correct Wavelength Units
Due to the variety of 1D spectral formats, SPAN requires you to manually set the wavelength units before loading spectra.
It does not matter whether the spectra use linear or logarithmic wavelength values—you must specify whether the units are in Angstrom, nanometres, or micrometres.
If the wavelength units are incorrect, SPAN will display the spectrum with the wrong scale.
Example: If your optical spectra are in Angstrom but you set nanometres, clicking "Plot" will show a range from 48000 to 55000 A, which is incorrect.
To fix this, select the correct wavelength units and plot one selected spectrum again.


2. Improving Lick/IDS Line-Strength Analysis Accuracy
The "Emission Line(s) Correction" option in the Lick/IDS line-strength analysis uses pPXF to remove gas emissions from the spectra. This process fits the gas component using E-MILES templates and applies a fixed wavelength fitting range (4000-6500 A). If the spectra contain edge artifacts (e.g., from X-shooter UVB and VIS arms), the fit may be imperfect.
Solution: Activate the "Cropping" task in the "Spectra Pre-Processing" panel and crop the spectrum to remove problematic regions. If you want to analyse only the 4800-5500 A range for stellar populations, crop the spectra accordingly.


3. Using "Stellar Populations and SFH" for Kinematics
The "Stellar Populations and SFH" module can be used for kinematics analysis, as an alternative to "Stars and Gas Kinematics". To adapt the settings for kinematics:
- Set Mult. degree = -1 (disable).
- Set Add. degree > 0.
- Enable "Convolve templates to galaxy resolution".
- Enable "Mask emission lines" if focusing on stellar components.
- Set Regul. error = 0 to disable regularization.


4. Saving and Loading Parameters
You can save your session settings, active tasks, and parameters as a .json file via "File --> Save Parameters...".
If you prefer to always load custom default parameters instead of SPAN's defaults, save your parameters via "File --> Save Parameters..." and overwrite "default_settings.json" in the "system_files" subfolder.


5. Managing Spectra Lists and File Organisation
Spectra loaded into SPAN are organised via a spectra list, generated using the "Generate Spectra List containing 1D spectra" button. Spectra files can be stored anywhere on your computer. The spectra list may also contain non-spectral files (e.g., documentation .txt files); SPAN will automatically ignore invalid files when loading the spectra list.
If your spectra folder contains subdirectories, SPAN will scan all subfolders for valid spectra. Use this feature carefully, as it may include unintended files.


6. Managing SPAN Results and Directories
SPAN stores results in the "SPAN_results" directory. The location of this directory must be decided by you the first time SPAN is run. You can move or delete its contents at any time, but it is recommended to do so only when SPAN is not running.
If SPAN detects missing directories upon startup, it will automatically recreate them.


7. Ensuring Correct Spectral Resolution
When using template-based spectral analysis (e.g., SFH or kinematics), ensure that the spectral resolution of your dataset is lower than that of the templates. If your spectra have higher resolution than the templates, use the "Degrade Resolution" task in the "Spectra Manipulation" panel to match the resolution of your spectra to the templates you want to use.


8. Handling High-Redshift Spectra (z > 0.01)
You have two options for handling redshifted spectra:
- Option 1 (Recommended for multiple tasks): De-redshift the spectra using the "Doppler/z Correction" in the "Spectra Manipulation" panel. Then, set redshift = 0 in the analysis tasks.
- Option 2: Enter the real redshift value directly in each spectral analysis task.

Important:
When de-redshifting spectra, the corrected spectral resolution is FWHM_gal/(1 + z). If spectra are not de-redshifted, SPAN will automatically adjust the resolution in kinematics and stellar population analyses.


9. Key Considerations for pPXF and Lick/IDS Analysis
When performing kinematics or line-strength index analysis with pPXF, the resolution of the templates must be higher (i.e., lower FWHM) than the resolution of the spectra.
If this is not the case, degrade the resolution using "Degrade Resolution" in the "Spectra Manipulation" panel.

The Lick/IDS index analysis uses a pre-loaded subset of E-MILES templates (Padova isochrones, FWHM = 2.51 A).
If your spectra have higher resolution (FWHM_spec < 2.51 A), degrade the spectra to at least 2.51 A before fitting.


10. Using sMILES Templates for Stellar Populations and SFH
The sMILES SSP templates included in SPAN are a subset of the full library, using a Salpeter unimodal IMF.
They are the only semi-empirical SSPs that sample a wide range of alpha/Fe values.
You can replace the content of the "sMILES_afeh" folder with a different sMILES SSP dataset (e.g., using a different IMF) and work with ani sMILES template.


11. Processing Spectra with Different Redshifts
If your dataset contains spectra with different redshifts (Delta(z) > 0.002) and you wish to use "Process All", you must first de-redshift them.
Use "Doppler/z Correction" --> "I have a file" In the "Spectra Manipulation" panel to apply individual redshifts.
This ensures that all spectra are in the rest frame before spectral analysis.


12. Handling High-Emission Line Galaxies
If analysing spectra with strong emission lines (i.e. with high star formation), using E-MILES or X-Shooter SSP templates may produce an artificially young and metal-poor population.
If this occurs, use the FSPS library, which includes SSPs as young as 1 Myr (compared to 30-60 Myr in E-MILES and X-Shooter).


13. Speeding Up Spectral Analysis
Cropping spectra to include only the wavelength regions of interest can significantly speed up spectral analysis.


14. Stellar populations with Lick/IDS indices
The determination of the stellar parameters with the Lick/IDS indices and the GPR method uses predictive machine-learning models and the scikit-learn module. The default models stored in the system_files folder have been trained with the scikit-learn 1.4.2 version. If you are using different version for this library, the predictions may be inaccurate and a warning is issued in the output window. In this case you should delete the trained models stored in the system_files folder (the files have the .pkl extension) and run again SPAN. If the Lick/IDS task with GPR stellar parameters determination does not find the trained models, it will create again the first time is executed, using the scikit-learn version you have installed. This will require a little time, but then you can use these new trained models until you change the version of the scikit-learn module. You should perform this operation when compiling the source code for Android systems using the Pydroid3 app, since the latest compatible version of the scikit-learn module is 1.3.0.
