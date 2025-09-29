SPAN: SPectral ANalysis software V7.0
Daniele Gasparri, September 2025

# SPAN user manual


## Purpose
SPAN is a Python 3.X multi-platform graphical interface program designed to perform operations and analyses on astronomical wavelength calibrated 1D spectra.
SPAN has been developed and optimized to analyze galaxy and stellar spectra in the optical and near infrared (NIR) atmospheric windows.
SPAN accepts as input ASCII and FITS spectra files.
The 1D spectra files required can be generated also with SPAN, both from long-slit 2D FITS images and 3D data cube (e.g. MUSE data) fully reduced and wavelength calibrated.
   
SPAN deals with linear sampled spectra, with wavelength in physical units (A, nm and mu). If you don't have linear sampled spectra, SPAN will try to read the spectra, will convert them automatically to linear sampling and will assign a physical wavelength scale, if needed. If these operations fails, your spectra will show a strange wavelength scale in the Preview window on the right or when clicking "Plot". If that is the case, you will need to adjust them with other software before load to SPAN.

The program has been tested with IRAF-reduced spectra, SDSS spectra, IRTF (also extended version) spectra, SAURON spectra, X-Shooter library spectra, JWST spectra, MUSE, CALIFA, JWST Nirspec, and WEAVE LIFU data cubes, (E)MILES, GALAXEV and FSPS stellar libraries, and complies with the ESO standard for 1D spectra. 
SPAN DOES NOT accept ASCII spectra files with Fortran scientific notation, like the PHOENIX synthetic stellar spectra. In this case, you will need to open the file and substitute the scientific notation of flux and wavelength "D" with "E" (you can do this operation even with the embedded text editor of SPAN).

Currently, SPAN considers only the wavelength and the flux, discarding the (potential) column with uncertainties.


## What do you need to run SPAN

- In order to run the source code, you need Python >=3.10 and the following modules installed (pip3 install <library>):
    1) Numpy
    2) Astropy
    3) Pandas
    4) Matplotlib
    5) Scipy
    6) scikit-image
    7) PyWavelets
    8) joblib
    9) scikit-learn
    10) ppxf
    11) vorbin
    12) certifi
    13) emcee
    
SPAN is optimized and can run also on most Android devices using the Pydroid3 app. The list and versions of packages needed is stored in the "README_ANDROID.txt" file.


 - A screen resolution of at least 1600X900 is required, otherwise the panel will be truncated. Optimal resolution: 1920X1080. 
 
    
## How SPAN works
SPAN can work with just one 1D spectrum, either in FITS or ASCII format, with the first column to be wavelength and the second flux. The wavelength and the flux of the FITS files must be in the primary HDU. 

SPAN accepts and processes a list of n 1D spectra, where n must be greater than 1. In order to do this, you need to create and load a text file containing the relative path of the spectra (with respect to the location of the main SPAN program), or the absolute path, and the complete spectra names. The first row of this list file must be commented with # and usually contains something like that: #Spectrum. You can put any type of 1D spectra in this file list, but I strongly suggest to insert spectra with at least the SAME wavelength unit scale.
It seems difficult, but don't worry: the button "Generate spectra list containing 1D spectra" will help you to create a spectra list file by selecting a folder containing the spectra you want to process.

You can find example file lists in the example_files directory. They are:

1. xshooter_vis_sample_list_spectra.dat, already preloaded in the main application (you just need to click "Load!"), contains 5 spectra of the central regions of nearby galaxies observed with the VIS arm of ESO XShooter spectrograph at resolution of R = 5000. Sampling is linear and the wavelength units to set ("Wavelength of the spectra is in:") are "nm";
2. ngc5806_bins.dat contains the spatial bins of a spiral galaxy observed with the TNG telescope at resolution FWHM = 3.5 A from 4700 to 6700 A. Sampling is logarithmic and wavelengths are in log(A). SPAN will take care of everything; you just need to set "A" in the "Wavelength of the spectra is in:" option of the "Prepare and load" frame before clicking "Load!";


## Quick start
If you installed SPAN as a Python package (pip3 install span-gui), just type in the terminal "span-gui".

At the first run, SPAN will ask you to download the auxiliary SSP spectral templates, which do not come with the Pypi or GIThub distribution for size issues. You can skip the download and SPAN will work, but the spectral analysis tasks devoted to full spectral fitting will use only the SSP sample provided by pPXF (EMILES, FSPS, GALAXEV, and, of course, any of the template that you will provide!).

At the first run, SPAN will also ask you to select the location of the SPAN_results folder, that is the folder where ALL the results will be saved. 


Once everything is set, press the "Load!" button to load the example files.

The spectra loaded will appear in the upper central frame (the white window). Just select one spectrum with the mouse, then look at the Preview on the right and interact with it. If you want more control on the plot, you can also click on the button "Plot". Close the plot to activate again the main panel.
You can analyze the selected spectrum by activating any of the spectral analysis tasks and/or you can modify the spectrum by opening the "Spectra manipulation" panel in the "Preview" frame on the right.

Let's open the "Spectra manipulation" panel and activate one of the many tasks, for example the "Add noise", then we confirm the choice by pressing the "Confirm" button. Now, we are back to the main panel and we press the "Preview spec." button to see the result. The preview will not change because it will always show you the original selected spectrum. If you like it, you can click the "Process selected" button to save this new noisy spectrum (but first you need to close the plot window!). If you press the "Process all" button, you will apply the task selected to all the loaded spectra. The results will be stored in the folder "SPAN_results/processed_spectra", located in the folder you have selected the first time you opened SPAN. The output window is your friend: it will tell you all the things that the program is doing.

Now, let's try something in the "Spectral analysis" frame. We activate the "Line-strength analysis" task and take a look at the parameters to set by clicking the button "Line-strength parameters". We select the "Single index option" and confirm the selection by clicking the button "Confirm". The "Line-strength analysis" window will close automatically and we are back to the main panel. Now, we click the "Preview result" button to see the result of this task.
The spectrum does look strange to you? Did you deactivate the "Add noise" task above? If not, the "Line-strength analysis" task is analyzing this noisy spectrum and not the original one.
The "Spectral Analysis" frame will consider the spectrum (or the spectra) processed by the tasks activated in the "Spectra manipulation" panel. If you activate 10 tasks, the spectrum processed in this frame will be the sum of all the activated tasks in the "Spectra manipulation" panel. If you don't activate any task, the spectrum processed will be the original loaded.

If you activated so many tasks that the entropy of the program tends to infinite, don't panic. Just click on the menu: "Edit --> Clear all tasks" to start from fresh. If you want to restore also the default parameters, you can do it with: "File --> Restore default parameters". If you want to save your personal parameters, you can do with: "File --> Save parameters" and load them again whenever you want.

If you want to change the location of the SPAN_results folder, you can do it with: "Edit --> Change result folder..."

You can play with other sample 1D spectra by loading the ready to use spectra list files provided with SPAN, for example the "ngc5806_bins.dat" located in the "example_files" folder. This spectra list contains 39 1D spectra in the optical window of the galaxy NGC 5806. Just browse this spectra list in the "3. Browse the spectra list or one spectrum" section of the "Prepare and load spectra" frame.


## General description and usage

SPAN is composed by a main graphical window that shows most of the spectral analysis tasks that the user can perform on 1D spectra.
In this window you will find two main panels: the top and the bottom one. 

### The upper panel
The top panel is divided in three frames. Here is a description of each one.

#### Prepare and Load Spectra frame
Any operation begins within the upper-left frame, called "Prepare and load spectra". There are four basic steps to load the spectra to SPAN and start the analysis.

1. **(Optional) Extract spectra from 2D or 3D fits:**
    This step is mandatory if you do not still have the 1D spectra needed by SPAN. It allows you to extract 1D spectra either from 2D fully reduced fits images or 3D fully reduced fits images, i.e. datacubes.
    - If you have 2D fits of long-slit spectra with the dispersion axis along the X axis of the image, press the "Long-slit extraction" button. There you can fit the trace, correct for distortion and/or slope and extract a 1D spectrum or a series of 1D spectra binned in order to reach a Signal to Noise (S/N) threshold;
    - If you have MUSE, CALIFA, WEAVE LIFU or JWST NIRSpec IFU datacubes, press the "DataCube extraction" button. To achieve the extraction, SPAN uses routines inspired to the famous GIST Pipeline (Bittner et al. 2019).
    Both these extraction routines will save the 1D spectra in the "SPAN_results" folder and a spectra file list in the directory SPAN_results/spectra_lists, ready to be loaded.

2. **Generate a spectra list with 1D spectra.** If you already have 1D spectra stored in a folder (and the relative subfolders, if any), you should click on the button "Generate a spectra list containing 1D spectra". You should then browse the folder where you stored your spectra. SPAN will read all the spectra contained in the selected folder and in any eventual subfolder and will create an ASCII file with their names and paths. The spectra list generated will be automatically loaded in the "Browse the spectra list or one spectrum". In case you want to load just a single 1D spectrum, you can skip this step.

3. **Browse the spectra list or one spectrum.** 
If you generated a spectra list in the previous step, this has been automatically loaded here. In this case, you should only select the wavelength units of the spectra contained in the spectra list. It is therefore important that all your spectra share the same wavelength units. It doesn't matter whether they are linearly or logarithmically rebinned, SPAN will read them correctly as far as you select the correct wavelength units.
In case your spectra list is already in your device (i.e. you skipped the step 2.) you should browse it, then select the right wavelength units of the spectra.
In case you just want to load a single 1D spectrum, just browse the spectrum and activate the option "I browsed a single spectrum".

4. **Finally load the spectra to SPAN.** 
This step is self explicative. Once you browsed the spectra list or the single spectrum and set the right wavelength units in step 3., here you need to press the "Load!" button to effectively load your spectra (or a spectrum) in the listbox on the right. Once done, select one spectrum in the listbox and check if everything is ok by pressing the "Plot" button. Since the official wavelength units of SPAN are Angstrom, you should check if the wavelength scale reproduced in the plot is actually correct. If not, you probably made a mistake in step 3., by setting the wrong wavelength units of your spectra. Try again with a different unit and press the "Plot" button again. Now the spectrum should be in the correct wavelength range.


#### Loaded Spectra
In this central frame you will find the listbox populated by your loaded spectra. You can select with the mouse any spectrum to activate it and perform operations and analysis. Once selected, SPAN will automatically load the spectrum and it will be ready for manipulation and analysis. You can use also the keyboard arrows to fast scroll the spectra in the list, which are activated by pressing the space bar of your keyboard. This listbox is static: you cannot modify or re-order the content unless you load a different spectra list file or spectrum.  


#### Preview
Since version 7, SPAN embeds a Preview window in the main interface to provide a fast way to visualize and inspect the currently loaded spectrum. This frame is based on an embedded Matplotlib figure that updates dynamically when you select a spectrum from the list.
When you hover the mouse over the previewed spectrum, some basic info are showed: wavelength, flux (expressed in the original units of your spectra) and an estimation of the S/N if performed, by analyzing in real time a small window around the mouse cursor (+- 50 points). You can also zoom using the touch capabilities (if available) of the trackpad, or the scroll button of your mouse. You can pan by left clicking and dragging. You can restore the original view by double click on the spectrum. 

**Visual redshift estimation in the Preview window**  

SPAN introduces a unique manual redshift estimation tool directly integrated into the preview window. By right-clicking and dragging the spectrum against fixed rest-frame line markers, you can quickly align prominent spectral features and obtain an approximate redshift value.

This functionality is not intended to replace quantitative fitting, but rather to provide a first-guess estimate in cases where the redshift is completely unknown. Such an estimate is essential for spectral analysis: for example, full spectral fitting with pPXF requires the input redshift to be known within approximately Delta(z) < 0.03 to ensure convergence and reliable results.

In the era of AI and automated black-box pipelines, there is still no substitute for visually inspecting spectra to understand where key features are located and to establish a sensible starting point. A “by eye” redshift estimation may still be the more accurate way to start with a solid base.  


**How it works**  

The basic principle is simple: resample the spectrum in a log spaced grid, therefore the redshift will be simply the linear (on the log space) shift of the whole spectrum from the starting point to the end point. The starting point is obvious, but what about the end point? We need some rest-frame reference lines.  
**Reference lines:** Fixed vertical markers are shown at the rest-frame wavelengths of OII, Hbeta, OIII, Halpha and CaT lines. These lines remain at fixed wavelengths in the preview and serve as anchors.  
**Right-click drag:** By pressing and holding the right mouse button, you can drag the entire spectrum horizontally, aligning its observed lines with the fixed rest-frame markers.  
**Cumulative shifting:** Multiple drags accumulate. Therefore, if you are out of dragging space you can pause, pan and zoom and start shifting the spectrum again until you find a match with the anchors.  
**HUD feedback:** While dragging, the HUD at the bottom-right of the preview shows the current estimated redshift, computed as the average shift relative to the reference lines.
**Reset:** A double right-click resets the spectrum to its original position and clears the estimated redshift. If you change spectrum, the redshift (and the shift) will automatically reset.


#### The Spectra manipulation panel
Since version 6.3, all the tasks devoted to the manipulation of the spectra have been grouped in the "Spectra manipulation" panel.
Here you will find some useful tasks that can be performed on spectra, grouped into the "Spectra pre-processing," "Spectra processing," and "Spectra math" frames. Any task executed within these frames modifies the selected spectrum and will have effects on the "Spectral analysis" frame.
You can choose multiple tasks (e.g., rebinning, dopcor, adding noise...) without limitations. The "Preview spec." button allows you to observe the effect of the task(s) performed.
The spectrum displayed and used in the "Spectral analysis" frame will be the one resulting from the selected tasks.

By default, the tasks are performed in series, following their order in the panel. No intermediate graphical information is available: if you activate three tasks, you will see the combined effect of all when you click the "Preview spec." button in the main panel. If you don't perform any task, don't worry: the original spectrum will be visible and ready to be used for spectral analysis.
You can change the order of the tasks performed. Activate the tasks you want to apply, then click the button "Reorder tasks" and change their order as you wish, then confirm the selection.

The four math tasks in the "Spectra math" frame that involve all the spectra ("Average all," "Norm. and average all," "Sum all," "Norm. and sum all") act on all the original spectra loaded (and don't work if you have loaded just one spectrum), and remain insensitive to other tasks performed. By activating the "Use for spec. an." option, you force the program to use the result of these operations for the spectral analysis, disregarding any other task performed on individual spectra. Be cautious in managing this option. In any case, a message in the terminal window will appear, indicating that you are using the combined original spectra for spectral analysis.


### The bottom panel
This panel is composed by two frames. The left one is the "Spectral analysis" frame, which contains the following tasks: 1) Blackbody fitting, 2) Cross-correlation, 3) Velocity dispersion, 4) Line-strength analysis, 5) Line(s) fitting, 6) Stars and gas kinematics, 7) Stellar populations and SFH.
Each task is independent from the others and does not modify the spectra.

The "Preview result" button will display the task(s) result on the selected spectrum in a graphic Matplotlib window and in the output frame on the right. If no task is selected, a warning message will pop-up when clicking the button.

The right frame displays the text output of the software. This is how SPAN communicates with you. This panel reproduces the computer terminal and shows the output of the operations performed, including errors and warnings.


### Apply the tasks
Once you are satisfied with your work, you can process the spectra or the single selected spectrum. The "Process selected" button will perform all the tasks activated in the "Spectra manipulation" panel and in the "Spectral analysis" frame, saving the new processed spectrum to a fits file. By default, the program will save intermediate spectra if more than one Spectra manipulation task is activated, i.e. one version for each activated task in the Spectra manipulation panel. For example, if you have selected rebinning, sigma broadening, and add noise, the program will save a spectrum with rebinning done, a second spectrum with rebinning + sigma broadening applied, and a third with rebinning + sigma broadening + add noise applied.
If you are not interested in saving all the intermediate spectra files modified in the "Spectra manipulation" panel, you can select the "Save final spectra" option at the very bottom of the "Spectra manipulation" panel, and only the spectrum at the end of the selected tasks (if any) will be saved. This is strongly recommended to do if you are applying more than one task with the reorder option activated. 
If you are planning to use the tasks of the "Spectra manipulation" panel just as preparatory phases to the spectral analysis, maybe you do not want to save the processed spectra every time you perform a spectral analysis task. In this case, in the "Spectra manipulation" panel you can select the option "Do not save processed spectra".
**IMPORTANT:** In the "Process selected mode", the results of the spectral analysis frame will be written only in the output frame.

By clicking "Process all", you will apply all the tasks to all the spectra in your list. This is the only way to save the results of the "Spectral analysis" frame in an ACII file. You can also store the plots generated during the spectral analysis by activating the option "Save spectral analysis plots" at the very bottom of the SPAN panel. The plots will be saved in high resolution PNG format and stored in the "plots" subdirectory of the "SPAN_results" folder.

### Zooming
Since version 7, SPAN includes a cross-platform zooming system designed to improve the readability of the GUI and the embedded Matplotlib preview. The zoom affects the entire application, including, GUI elements (buttons, text labels, checkboxes, frames, etc.), subwindows opened during the workflow, Matplotlib figures opened in external windows.

The zoom is controlled through the View → Zoom In / Zoom Out / Reset Zoom menu entries. Zoom In increases the size of all fonts and widgets. Zoom Out decreases them. Reset Zoom restores everything to the original size as defined by the layout.

Because Tkinter behaves differently across platforms, the implementation includes a few platform-specific adjustments. Windows (tested on Windows 10/11): Zooming and reset work as expected with no additional adjustments. Linux (tested on Ubuntu 22.04, X11 and Wayland): By default, Tk tends to stretch canvases after a resize. SPAN explicitly disables this behavior at reset, so the embedded preview canvas returns to its original dimensions. macOS: The zooming behavior is expected to follow Windows closely, since Tk for macOS manages widget geometry in a similar way. This has not yet been extensively tested, so user feedback is welcome.


## The input files

In order to work properly, SPAN sometimes needs input text files containing information about your data. To see how they must be formatted, please take a look at those coming with SPAN and already set by default in the graphic interface.

**IMPORTANT:** The text files MUST always have the first line as header, identified by # (e.g. #spectrum)
          
1. **Spectra list file:**
    It is essential. If you don't believe it, try to perform any task without upload the spectra and you will see the effects! It is just an ASCII file containing the path (relative if they are in a subfolder of SPAN, absolute if they are elsewhere) and the complete names (with file extension) of the spectra you want to process. You can use any spectra you want, with different format (fits, ASCII...) and resolutions, but it is mandatory to use spectra with the same wavelength units. If you just want to play with one spectrum, then load the ASCII or fits 1D spectrum and activate the option "I browsed a single spectrum" before clicking the button "Load!".

                                            example_list.dat
                                            
                                            #filename ---> header: always necessary!
                                            [path/]spectrum1.fits
                                            [path/]spectrum2.fits
                                            [path/]spectrum3.fits
                                            [path/]spectrum4.fits


Other ASCII files may be needed in the "Spectra manipulation" panel for some specific tasks. They are:
                                            
2. **Doppler correction file** for the "Doppler/z correction" task and the "I have a list file" option selected:
    It is an ASCII file containing two columns: 1) Name of the spectrum and 2) Radial velocity to correct to the spectrum. This file has the same format of the output text file generated by the Cross-correlation task, so you can directly use it. 
                    
                                            example_dopcor.dat
                                        
                                        #spectrum       RV(km/s) ---> header: always necessary!
                                        [path/]spectrum1.fits  1000
                                        [path/]spectrum2.fits  1001
                                        [path/]spectrum3.fits  1002
                                        [path/]spectrum4.fits  1003
                                            
3. **Heliocentric correction file** for the "Heliocentric correction" task and the "I have a file with location..." option selected:
    It is an ASCII file containing three columns, separated by a space, following the SAME order of your spectra in the loaded spectra list file: 1) Name of the location, 2) Date of the observation (just year, month, day, not the hour), 3) RA of the object (format: degree.decimal), 4) Dec. of the object (format: degree.decimal).
    
                                            example_heliocorr.dat
                                        
                                #where  date        RA          Dec
                                paranal 2016-6-4    4.88375     35.0436389
                                paranal 2016-6-30   10.555      1.11121
                                aao     2011-12-24  -50.034     55.3232
                                aao     2018-2-13   -11.443     11.2323
                                SRT     2020-7-31   70.234      55.32432


Some external files may be needed for specific options of the "Spectral analysis" tasks. They are:
                            
4. **Cross-correlation and velocity dispersion tasks:**
    These task require a single template, in fits or ASCII format (i.e.. just a spectrum!).
    
    
5. **Line-strength analysis task** and the option "User indices on a list file" selected:
    It is an ASCII text file containing the index definitions. One index per column. Don't mess it up with the index file, otherwise you will obtain inconsistent results! Luckily, you can always test a single index and see the graphical preview before running the wrong indices on 240913352 spectra and waste one year of your life.
    
                                            example_idx_list_file.dat

                                    #Idx1    Idx2  ---> header: always necessary!
                                    8474   8474 ---> row2: left blue continuum band, in A
                                    8484   8484 ---> row3: right blue continuum band, in A
                                    8563   8563 ---> row4: left red continuum band, in A
                                    8577   8577 ---> row5: right red continuum band, in A
                                    8461   8484 ---> row6: left line limits, in A
                                    8474   8513 ---> row7: right line limits, in A


6. **Calculate velocity dispersion coefficients**, located in the "Line-strength parameters" sub-window :
    It determines 4 spline correction coefficients in order to correct the equivalent width of galaxy spectra broadened by the velocity dispersion. It needs a sample of unbroadened spectra that are a good match of the expected stellar populations of the galaxy spectra you want to correct to the zero velocity dispersion frame. The input file is just an ASCII file containing the list of the spectra used as sample, i.e., a normal spectra list!

                                            example_coeff_determ.dat
                                            
                                            #filename ---> header: always necessary!
                                            [path/]stellar_spectrum1.fits
                                            [path/]stellar_spectrum2.fits
                                            [path/]stellar_spectrum3.fits
                                            [path/]stellar_spectrum4.fits

                                    
7. **Correct the line-strength for velocity dispersion task:** 
    To apply the velocity dispersion coefficients and correct the raw equivalent widths to the zero velocity dispersion frame, you need this task and three files: 
    - Sigma list file: a file containing the name of the spectra, the velocity dispersion and the relative uncertainties. It has the same format of the output file generated by the Velocity dispersion task. 
        
                                            example_sigma_vel.dat
                                            
                                            #Spectrum       Sigma(km/s) err ---> Header: always necessary
                                            spectrum_name1  166.2       3.0
                                            spectrum_name2  241.5       3.1
                                            spectrum_name3  335.1       6.2
                                            spectrum_name4  241.5       3.2
        
    - EW file list to correct: the text file containing the raw equivalent widths you want to correct. It has the same format of the output file generated by the Line-strength measurement task. BE CAREFULL to check that the indices are in the EXACT same order of those you used in the "Calculate velocity dispersion coefficients" task for the correction coefficient determination.
        
                                            
                                            example_uncorrected_ew.dat
                                            
                            #Spectrum       idx1    idx2    idx3   idx1err idx2err idx3err
                            spectrum_name1  0.27    1.38    3.56    0.01     0.01    0.02
                            spectrum_name2  0.15    1.32    3.43    0.01     0.02    0.02
                            spectrum_name3  0.08    0.75    2.81    0.01     0.02    0.02
                            spectrum_name4  0.14    1.25    3.18    0.01     0.01    0.01

        
    - Correction coefficients file: it is the output file generated by the "Calculate velocity dispersion coefficients task".
                                    
                                            example_correction_coeff.dat
                        #Pa1          Ca1           Ca2          Pa1e         Ca1e         Ca2e
                        4.3282e-08   1.06712e-08  -2.7344e-09  -5.7463e-09   2.2911e-09   2.8072e-10
                       -2.9602e-05  -1.2012e-05   -3.5782e-07   3.9353e-06  -1.9246e-06  -2.9293e-07
                        0.0017       0.0021        8.5793e-05  -0.0001       0.0004       9.9212e-05
                       -0.0029      -0.0085       -0.0016       0.0053      -0.0003      -0.0002

                                    

## File organization
SPAN generates different types of **results**, which are all stored in the "SPAN_results" folder:

- **Extracted spectra** from the "long-slit extraction" and the "Datacube extraction" sub-programs. The spectra extracted from long-slit data are stored in the "longslit_extracted" folder. The spectra extracted from datacube data are stored in "RUN_NAME" folder, where "RUN_NAME" is the arbitrary name you had to set in the "Datacube extraction" sub-program.
- **Processed spectra** in FITS format, both in the "Process selected" and "Process all" mode. These are processed spectra from the "Spectra Manipulation" panel or auxiliary spectra generated from Spectral analysis tasks (e.g. best fit and residuals from the "Stars and gas kinematics" and "Stellar populations and SFH" tasks). These spectra are stored in the "processed_spectra" folder. 
- **ASCII files** in plain text .dat format, containing the results of the Spectral analysis tasks, which are generated only in the "Process all" mode. These products are saved in specific folders with the same name of the spectral analysis task applied.
- **Plots** in high resolution (300 dpi) PNG images. They are generated only for the Spectral analysis tasks and are the plots displayed also in the "Preview result" mode. These plots are saved only in "Process all" mode and if the option "Save spectral analysis plots" is activated. If you just need one specific plot for one spectrum in the list, you can save it directly from the Matplotlib window that opens in the "Preview result" mode. These plots are stored in the "plots" folder. 

The **spectra list** files generated by the "Generate spectra list containing 1D spectra" are saved in the "spectra_lists" folder within the "SPAN_results" main folder.



## List of operations you can perform

SPAN can perform many operations on the spectra.

**WARNING:** All the wavelengths of SPAN are given in A, in air, and all the velocities are in km/s.

Here is a description of any task:

1. **Utilities (Window --> Utilities):**  
    - Show the header of the selected spectrum = shows the header of the selected spectrum, only for FITS files;
    - Show the wavelength step of the spectrum = shows the step of the selected spectrum;
    - Estimate the resolution = calculates the resolution of the selected spectrum by trying to fit an emission sky line. In The W1 and W2 you should put a small wavelength interval containing a sky line: it's up to you!
    - Convert the spectrum to = converts the selected spectrum to ASCII or Fits;
    - Compare spectrum with = compares the selected spectrum with another one selected by you. This comparison spectrum should have the same wavelength units;
    - Convert Flux = converts the flux from frequency to lambda and vice-versa. The buttons "see plot", "save one" and "save all" are active to see and save the results for one or all the spectra;
    - S/N = measures the Signal to Noise in the selected spectrum, in the central wavelength selected (W.). The buttons "save one" and "save all" are active to save one or all the S/N computed for the spectra.

2. **Spectra manipulation panel:**  
    - Spectra pre-processing frame
        - Cropping = performs a simple cropping of the spectra. If the wavelength window to crop is outside the spectrum, SPAN will ignore the task and will not perform the crop;
        - Dynamic cleaning = performs a sigma clipping on the spectra. The sigma clip factor, the resolving power of the spectrum and the velocity dispersion (instrumental and/or intrinsic) of the selected spectrum is required in order to perform a better cleaning. For the "Process all" mode, the option "R and sigma vel file" is available in order to have R (resolution) and sigma values for all the spectra to be processed. Be VERY careful to use this task with strong emission line spectra;
        - Wavelet cleaning = performs a wavelet denoise of the spectra. The mean standard deviation of the spectra continuum (sigma) and the number of wavelet layers to consider are required. You don't need to measure it, just try different values. Be careful to not delete the signal;
        - Filtering and denoising = smooths the spectra by performing some denoising filters: box window moving average, gaussian kernel moving average, low-pass Butterworth filter and band-pass Butterworth filter;
        - Dopcor/z correction = performs the doppler or z correction of the spectra. Single shot option with user input value of radial velocity (in km/s) or z is available both for one or all the spectra. "I have a file" option only works with the "Process all" mode: you need a text file with the spectra name and the recession velocities or z values. This file can be generated by the "Cross-correlation" task in "Process all" mode;
        - Heliocentric correction = performs the heliocentric correction on the spectra. The "Single" option requires a location which can be selected from the "loc.list" button (it requires an internet connection the first time!). The other fields are the date in the format YYYY-MM-DD and the RA and Dec. of the observed object (in decimals). In the "I have a file" option, available only for the "Process all" mode, a list file with location, date, RA and Dec. coordinates for each object is required.
    - Spectra processing frame
        - Rebin = performs a rebin/resample of the spectra in linear wavelength step ("pix.lin" option, with the step in A) and in sigma linear step ("sigma lin." option, with the sigma step in km/s);
        - Degrade resolution = degrades the resolution of the spectra from R to R, from R to FWHM and from FWHM to FWHM;
        - Normalize spectrum to = normalizes the spectra to the wavelength provided (in A);
        - Sigma broadening = broads the spectra by convolving with a gaussian function with the standard deviation provided by you, in km/s. Remember that the real broadening of the spectra will be the quadratic sum between the broadening and the instrumental sigma;
        - Add noise = adds a random Poisson noise to the spectra with a S/N defined by you. Remember that the final S/N of the spectra will the sum in quadrature between the added noise and the intrinsic S/N of the spectra;
        - Continuum modelling = models the continuum shape with two options: 1) Simple filtering of the continuum by reducing the spectrum to a very small resolution (R = 50), and 2) polynomial fitting, with the possibility to mask emission/contaminated regions. Both the continuum models can be divided or subtracted to the original spectrum;
    - Spectra math
        - Subtract normalized average = subtracts to the spectra the normalized average made from all the spectra loaded;
        - Subtract norm. spec. = subtracts to the spectra a normalized spectrum selected by you, which shares the same wavelength units;
        - Add constant = adds a constant to the spectra;
        - Multiply by a constant = multiplies the spectra by constant value.
        - Calculate first and second derivatives = automatic calculation of the derivatives of the spectra. This task does not modify the original spectra and the derivative spectra cannot be directly used for spectral analysis.
        - Average all = averages all the spectra (only available in "Process selected" mode);
        - Norm. and average all = normalizes to a common wavelength and average all the spectra (only available in "Process selected" mode);
        - Sum all = sums all the spectra (only available in "Process selected");
        - Norm. and sum all = Normalizes and sum all the spectra (only available in "Process selected"). The option "Use for spec. an." forces the program to use the result of one of these 4 operations for the following spectral analysis.


3. **Spectral analysis:**  
    - Blackbody fitting = performs a fit of the spectrum with Planck's blackbody equation and gives the temperature estimation. It works with any type of spectra but it performs better for stellar spectra, with wide (at least 5000 A) wavelength range;
    b)  Cross-correlation = performs a cross-correlation of the spectra with any template. You can smooth the template to a velocity dispersion value in order to improve the cross-correlation and should identify a narrow region of the spectrum to be cross-correlated (tip: the Calcium triplet lines are the best features in the NIR);
    - Velocity dispersion = performs the measurement of the velocity dispersion of the spectra with a rough (but fats) fit it with any template. Some pre-loaded bands in the visible and NIR are shown but you can input any band. The routine succeeds with strong features (the CaT is the best). It is a little rough but very fast and gives reasonably accurate results;
    - Line-strength analysis = performs the equivalent width measurement of the spectra, with a single index, with a list of indices or with the Lick/IDS system. The results are provided in Angstrom. MonteCarlo simulations are run for the uncertainties estimation. The calculation of the Lick/IDS indices can be personalized in many ways: you can correct for the emission, for the velocity dispersion and the recession velocity. You can also perform a linear interpolation with the SSP models of Thomas et al. 2010, Xshooter, MILES and sMILES to retrieve the age, metallicity and alpha-enhancement (not available for the Xshooter models) of the stellar populations via linear interpolation or with machine-learning pre-trained models (Gaussian Process Regression). From the "Line-strength parameters" window, it is possible also to perform the "Calculate velocity dispersion coefficients" task. This task broadens a sample of SSP spectra up to 400 km/s and calculates the deviation of the equivalent width of the indices contained in the index file provided. It works only by pressing the "Compute!" button and creates a text file with a third order polynomial curve that fits the behavior of the broadened index (or indices). The "Correct the line-strength for velocity dispersion" task performs the correction of the equivalent widths based on the coefficients estimated with the "Calculate velocity dispersion coefficients" task. It works only by pressing the "Correct!" button and requires the raw equivalent width measurements stored in the ASCII file generated previously, with the same indices in the same order to that considered in the "Calculate velocity dispersion coefficients". The output files of the "Line-strength analysis", "Calculate velocity dispersion coefficients" and "Velocity dispersion" are ready to be used for this task, if we are considering the same spectra and indices;
    - Line(s) fitting = performs the fitting of a line in the inserted wavelength range using a combination of a Gaussian function to model the spectral line and straight line for the continuum. If "CaT lines" is selected, the task will perform an automatic fitting of the Calcium Triplet lines, assuming they have been previously corrected for redshift and/or Doppler velocity;
    - Stars and gas kinematics = uses the pPXF algorithm of Cappellari et al. 2023 to fit a wavelength region of the spectra with a combination of templates. You can select the template library you prefer among the pre-loaded EMILES, GALAXEV, FSPS and XSHOOTER, or you can use your custom set of templates. You can decide how many moments to fit, whether fit only the stellar component or also the gas, whether estimate or not the uncertainties with MonteCarlo simulations and much more. It returns the radial velocity, the velocity dispersion and the higher moments up to H6 (if needed, and a nice plot courtesy of Cappellari), as well as the spectra product (best fit, residuals, emission corrected spectra, if any, gas spectra, if any, best fit gas spectra, if any);
    - Stellar populations and SFH = uses pPXF to fit a wavelength region of the spectra with a combination of templates. You can select the template library you prefer among the EMILES, GALAXEV, FSPS, XSHOOTER and sMILES, add any EMILES custom library, or any .npz file following the pPXF standard. You can decide whether include the gas emission or not, the reddening, the order of multiplicative and additive polynomials of the fit, the age and metallicity range of the templates, and much more. It returns a beautiful plot, the kinematics, the weighted age (in luminosity and mass), metallicity (in luminosity and mass), the M/L, the SFH and saves the best fit template and the emission corrected spectra (if any). Works great in the visible and in the NIR, but this depends on the quality of your spectra.
    


# The sub-programs

The two light-blue buttons in the upper left corner of SPAN (in the "Prepare and load spectra" frame) are sub-programs that might help you to generate the 1D spectra needed. Here is how they works:

1. **Long-slit extraction:** Allows the extraction of a single 1D spectrum or a series of 1D spectra from a reduced and wavelength calibrated 2D fits image containing the long-slit spectrum of a source, with dispersion axis along the X-axis and the spatial axis along the Y-axis.
Before proceed to the extraction, you need to load a valid 2D fits image, then you need to:
    - Open the spectrum and see if everything is ok;
    - Fit the photometric trace in order to find the maximum along the dispersion axis. You need to set the degree of polynomial curve that will be used to fit the trace and correct the distortion and slope of the spectrum;
    - Correct the spectrum for distortion and slope using the model trace obtained in the previous step. Then, you can:
    - Extract and save only one 1D spectrum within the selected Y range (useful for point sources);
    - Extract and save a series of n 1D spectra covering all the spatial axis and obtained by binning contiguous rows in order to reach the desired S/N. A spectra list file ready to be loaded to SPAN is also generated, as well as a text file containing the position of the bins relative to the central region of the galaxy and the S/N. 
    The S/N threshold that you must insert is just a very rough estimation of the real S/N. A good starting value to produce 1D spectra with bins with realistic S/N > 30 is 20. Adjust the SNR Threshold to your preference by looking at the real S/N of the bins.
    The pixel scale parameter is optional. If you set to zero it will not be considered. This option is useful if you have the spectrum of an extended source (e.g. a galaxy) and want to sample different regions.

2. **DataCube extraction:** Following the GIST pipeline standard (Bittner at al., 2019, 2021), this sub-program allows you to extract 1D spectra from MUSE and CALIFA DataCubes using the Voronoi binning (Cappellari et al., 2003) or manual binning. It also allows to visualize the DataCube loaded and dynamically create a mask (if needed).


In the menu bar you can find more sub-programs that might help you in the difficult work of analyzing and processing astronomical spectra. They work independently from the main program, so you can also not load any spectra if you don't need to perform tasks on them. Here is how they works:

1. **Text editor:** A simple ASCII file editor where you can create, read or modify ASCII files, included those generated by the SPAN tasks. Some basics operations are available, such find, replace and merge rows;

2. **FITS header editor:** An header editor to add, remove and save the keywords of fits header files. You can select between: "Single fits header editor" to work with the keywords of one fits file, "List of fits header editor" to modify the keywords of a list of fits files, "Extract keyword from list" to extract and save in an ASCII file one or more keywords from the headers of a list of fits files;

3. **Plot data:** A sub-program to plot the data generated by the "Spectral analysis" frame and, in general, all the data stored in ASCII space-separated data. Once you browse for the text file and click the "Load" button, the program will automatically recognize the column names. Select a name for the x and y axis and plot the data to see them in an IDL style plot.
You can personalize the plot by adding the error bars, set the log scale, add a linear fit (simple fit without considering the uncertainties), set the labels, the range, the font size, size and colors of the markers and decide if visualize the legend or not. You may also save the plot in high resolution PNG image format, in the directory where you run SPAN.
If any error occur, the program will warn you.

4. **Plot maps:** A sub-program to plot 2D maps from extracted quantities from datacubes. If you extracted datacube spectra and performed some spectral analysis with SPAN, you can load the RUNNAME_table.fits file stored in the extracted spectra folder within the "SPAN_results" folder, and any of the ASCII files generated by the spectral analysis tasks in the "Process all" mode. SPAN will show you the quantities available from this file and generate beautiful and customizable 2D maps.

5. **Utilities:** A standalone frame that allows you to find out information about the selected spectrum, such as the header, the sampling (in A), the S/N, or simply convert the spectrum to ASCII or binary fits. 



# Tricks in the menu bar

The menu bar was introduced in version 4.5 of SPAN, offering several helpful options to enhance your experience with spectral analysis. Here is a detailed overview of some options that you won't find in the main panel (unless you are using the Android version):

1. File --> Save Parameters...: Allows to save all parameters and values from the main panel and the various parameter windows of the tasks in a .json file.
This feature is very useful as it enables you to preserve any modifications made to parameters, facilitating the restoration of your session each time you reopen SPAN;
2. File --> Load Parameters...: Allows to load the parameters saved in the .json file;
3. File --> Restore Default Parameters: Resets all the parameters to their default values. Useful if numerous parameter modifications during a lengthy session have resulted in issues, allowing you to start from fresh;
4. Edit --> Clear All Tasks: Immediately deactivates all tasks activated during the session, enabling a clean restart;
5. Edit --> Clean Output: Deletes the content of the output window. Particularly useful during extended sessions where the generated output may become quite large.
6. Edit --> Show result folder: shows the location of the "SPAN_results" folder in case you forgot;
7. Edit --> Change result folder...: create a new "SPAN_results" folder wherever you want.

Please, report any bug or comment to daniele.gasparri@gmail.com
Have fun!

Daniele Gasparri  
2025-09-20  
Greetings from the Atacama desert!
