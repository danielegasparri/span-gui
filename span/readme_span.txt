SPAN: SPectral ANalysis software V6.6
Daniele Gasparri, June 2025

****Purpose****
SPAN is a Python 3.X multi-platform graphical interface program designed to perform operations and analyses on astronomical wavelength calibrated 1D spectra.
SPAN has been developed and optimised to analyse galaxy and stellar spectra in the optical and near infrared (NIR) bands.
SPAN accepts as input ASCII and fits spectra files.
The 1D spectra files required can be generated also with SPAN, both from long-slit 2D fits images and 3D data cube (e.g. MUSE data) fully reduced and wavelength calibrated.
   
SPAN deals with linear sampled spectra, with wavelength in physical units (A, nm and mu). If you don't have linear sampled spectra, SPAN will try to read the spectra, will convert them automatically to linear sampling and will assign a physical wavelength scale, if needed. If these operations fails, your spectra will show a strange wavelength scale when clicking 'Plot'. If that is the case, you will need to adjust them with other software before load to SPAN.

The program has been tested with IRAF-reduced spectra, SDSS spectra, IRTF (also extended version) spectra, SAURON spectra, X-Shooter library spectra, JWST spectra, MUSE and CALIFA data cubes, (E)MILES, GALAXEV and FSPS stellar libraries, and complies with the ESO standard for 1D spectra. 
SPAN DOES NOT accept ASCII spectra files with Fortran scientific notation, like the PHOENIX synthetic stellar spectra. In this case, you will need to open the file and substitute the scientific notation of flux and wavelength 'D' with 'E' (you can do this operation even with the embedded text editor of SPAN).

Currently, SPAN considers only the wavelength and the flux, discarding the (potential) column with uncertainties.


****What do you need to run SPAN****
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
    
    SPAN is optimised and can run also on most Android devices using the Pydroid3 app. The list and versions of packages needed is stored in the 'README_ANDROID.txt' file.


 - A screen resolution of at least 1600X900 is required, otherwise the panel will be truncated.
 
    
****How SPAN works****
SPAN can work with just one 1D spectrum, either in Fits or ASCII format, with the first column to be wavelength and the second flux. The wavelength and the flux of the fits files must be in the primary HDU. 

SPAN can load and process a list of n 1D spectra, where n must be greater than 1. In order to do this, you need to create and load a text file containing the relative path of the spectra (with respect to the location of the main SPAN program), or the absolute path, and the complete spectra names. The first row of this list file must be commented with # and usually contains something like that: #Spectrum. You can put any type of 1D spectra in this file list, but I strongly suggest to insert spectra with at least the SAME wavelength unit scale.
It seems difficult, but don't worry: the button 'Generate spectra list containing 1D spectra' will help you to create a spectra list file by selecting a folder containing the spectra you want to process.

You can find example file lists in the example_files directory. They are:
1) xshooter_vis_sample_list_spectra.dat, already preloaded in the main application (you just need to click 'Load!'), contains 5 spectra of the central regions of nearby galaxies observed with the VIS arm of ESO XShooter spectrograph at resolution of R = 5000. Wavelength units are in nm. Sampling is linear and the wavelength units to set ('Wavelength of the spectra is in:') are 'nm';
2) ngc5806_bins.dat contains the spatial bins of a spiral galaxy observed with the TNG telescope at resolution FWHM of 3.5 A from 4700 to 6700 A. Sampling is logarithmic and wavelengths are in log(A). SPAN will take care of everything; you just need to set 'A' in the 'Wavelength of the spectra is in:' option of the 'Prepare and load' frame before clicking 'Load!';


****Quick start****
If you installed SPAN as a Python package (pip3 install span-gui), just type in the terminal 'span-gui' in your working directory. This can be any folder on your device but don't forget its location, because SPAN will save the spectra list files here. 
If you want to compile the source code go to the root folder of SPAN and type in the terminal: python3 __main__.py. In this case, the spectral list files will be saved here. 

At the first run, SPAN will ask you to download the auxiliary SSP spectral templates, which do not come with the Pypi of GIThub distribution for size issues. You can skip the download and SPAN will work, but the spectral analysis tasks devoted to full spectral fitting will use only the SSP sample provided by pPXF (EMILES, FSPS, GALAXEV, and, of course, any of the template that you will provide!).

At the first run, SPAN will also ask you to select the location of the SPAN_results folder, that is the folder where ALL the results will be saved. 


Once everything is set, press the 'Load!' button to load the example files.

The spectra loaded will appear in the upper central frame (the white window). Just select one spectrum with the mouse, then click 'Plot' to see it. Close the plot to activate again the main panel.
You can analyse the selected spectrum by activating any of the spectral analysis tasks and/or you can modify the spectrum by opening the 'Spectra manipulation' panel in the 'Utilities' frame on the right.
Let's open the 'Spectra manipulation' panel and activate one of the many tasks, for example the 'Add noise', then we confirm the choice by pressing the 'Confirm' button. Now, we are back to the main panel and we press the 'Preview spec.' button to see the result. If you like it, you can click the 'Process selected' button to save this new noisy spectrum (but first you need to close the plot window!). If you press the 'Process all' button, you will apply the task selected to all the loaded spectra. The results will be stored in the folder 'SPAN_results/spec', located in the folder you have selected the first time you opened SPAN. The output window is your friend: it will tell you all the things that the program is doing.

Now, let's try something in the 'Spectral analysis' frame. We activate the 'Line-strength analysis' task and take a look at the parameters to set by clicking the button 'Line-strength parameters'. We select the 'Single index option' and confirm the selection by clicking the button 'Confirm'. The 'Line-strength analysis' window will close automatically and we are back to the main panel. Now, we click the 'Preview result' button to see the result of this task.
The spectrum does look strange to you? Did you deactivate the 'Add noise' task above? If not, the 'Line-strength analysis' task is analysing this noisy spectrum and not the original one.
The 'Spectral Analysis' frame will consider the spectrum (or the spectra) processed by the tasks activated in the 'Spectra manipulation' panel. If you activate 10 tasks, the spectrum processed in this frame will be the sum of all the activated tasks in the 'Spectra manipulation' panel. If you don't activate any task, the spectrum processed will be the original loaded.

If you activated so many tasks that the entropy of the program tends to infinite, don't panic. Just click on the menu: 'Edit --> Clear all tasks' to start from fresh. If you want to restore also the default parameters, you can do it with: 'File --> Restore default parameters'. If you want to save your personal parameters, you can do with: 'File --> Save parameters' and load them again whenever you want.

If you want to change the location of the SPAN_results folder, you can do it with: 'Edit --> Change result folder...'

You can play with other sample 1D spectra by loading the ready to use spectra list files provided with SPAN, for example the 'ngc5806_bins.dat' located in the 'example_files' folder. This spectra list contains 39 1D spectra in the optical window of the galaxy NGC 5806.Just browse this spectra list in the '3. Browse the spectra list or one spectrum' section of the 'Prepare and load spectra' frame.


#################################################################################################
################################### General description and usage ###############################

SPAN is composed by a main graphical window that shows most of the spectral analysis tasks that the user can perform on 1D spectra.
In this window you will find two panels separated by a horizontal line.

###The upper panel###
This top panel is divided in three frames.
Any operation begins within the upper-left frame, called 'Prepare and load spectra'. There are four basic steps to load the spectra to SPAN and start the analysis.
- 1. Extract spectra from 2D or 3D fits:
    This step is mandatory if you do not still have the 1D spectra needed by SPAN. It allows you to extract 1D spectra either from 2D fully reduced fits images or 3D fully reduced fits images, i.e. datacube.
    - If you have 2D fits of long-slit spectra with the dispersion axis along the X axis of the image, press the 'Long-slit extraction' button. There you can fit the trace, correct for distortion and/or slope and extract a 1D spectrum or a series of 1D spectra binned in order to reach a Signal to Noise (S/N) threshold;
    - If you have MUSE or CALIFA data cubes, press the 'DataCube extraction' button. To achieve the extraction, SPAN uses routines inspired to the famous GIST Pipeline (Bittner et al. 2019).
Both these extraction routines will save the 1D spectra in the 'SPAN_results' folder and a spectra file list in the directory where you are running SPAN, ready to be loaded.

- 2. Generate a spectra list with 1D spectra. If you already have 1D spectra stored in a folder (and the relative subfolders, if any), you should click now on the button 'Generate a spectra list containing 1D spectra'. You should then browse the folder where you stored your spectra. SPAN will read all the spectra contained in the selected folder and in any eventual subfolder and will create an ASCII file with their names and paths. The spectra list generated will be automatically loaded in the 'Browse the spectra list or one spectrum'. In case you want to load just a single 1D spectrum, you can skip this step.

- 3. Browse the spectra list or one spectrum.
If you generated a spectra list in the previous step, this has been automatically loaded here. In this case you should only select the wavelength units of the spectra contained in the spectra list. It is therefore important that all your spectra share the same wavelength units. It doesn't matter whether they are linearly or logarithmically rebinned, SPAN will read them correctly as far as you select the correct wavelength units.
In case your spectra list is already in your device (and skipped the step 2.) you should browse it, then select the right wavelength units of the spectra.
In case you just want to load a single 1D spectrum, just browse the spectrum and activate the option 'I browsed a single spectrum'.

- 4. Finally load the spectra to SPAN.
This step is self explicative. Once you browsed the spectra list or the single spectrum and set the right wavelength units in step 3., here you need to press the 'Load!' button to effectively load your spectra (or a spectrum) in the listbox on the right. Once done, select one spectrum in the listbox and check if everything is ok by pressing the 'Plot' button. Since the official wavelength units of SPAN are Angtrom, you should check if the wavelength scale reproduced in the plot is actually correct. If not, you probably made a mistake in step 3., by setting the wrong wavelength units of your spectra. Try again with a different unit and press the 'Plot' button again. Now the spectrum should be in the correct wavelength range.


The right frame (Utilities) is a standalone frame that allows you to find out information about the selected spectrum, such as the header, the sampling (in A), the S/N, or simply convert the spectrum to ASCII or binary fits. Finally, here you will find the 'Spectra manipulation panel' and some other sub-programs that you may find useful. We'll talk about them later.


### The Spectra manipulation panel ###
Since version 6.3, all the tasks devoted to the manipulation of the spectra have been grouped in the 'Spectra manipulation' panel.
Here you will find some useful tasks that can be performed on spectra, grouped into the 'Spectra pre-processing,' 'Spectra processing,' and 'Spectra math' frames. Any task executed within these frames modifies the selected spectrum and will have effects on the 'Spectral analysis' frame.
You can choose multiple tasks (e.g., rebinning, dopcor, adding noise...) without limitations. The 'Preview spec.' button allows you to observe the effect of the task(s) performed. 
The spectrum displayed and used in the 'Spectral analysis' frame will be the one resulting from the selected tasks.

By default, the tasks are performed in series, following their order in the panel. No intermediate graphical information is available: if you activate three tasks, you will see the combined effect of all when you click the 'Preview spec.' button in the main panel. If you don't perform any task, don't worry: the original spectrum will be visible and ready to use for spectral analysis. 
You can change the order of the tasks performed. Activate the tasks you want to use, then click the button 'Reorder tasks' and change their order as you wish, then confirm the selection.

The four math tasks in the 'Spectra math' frame that involve all the spectra ('Average all,' 'Normalise and average all,' 'Sum all,' 'Norm. and sum all') act on all the original spectra loaded (and don't work if you have loaded just one spectrum), and remain insensitive to other tasks performed. By activating the 'Use for spec. an.' option, you force the program to utilise the result of these operations for the spectral analysis, disregarding any other task performed on individual spectra. Be cautious in managing this option. In any case, a message in the terminal window will appear, indicating that you are using the combined original spectra for spectral analysis.


###The bottom panel###
This panel is composed by two frames. The left one contains basic and more advanced spectral analysis tools.
The 'Spectral analysis frame' contains the following tasks: 1) Blackbody fitting, 2) Cross-correlation, 3) Velocity dispersion, 4) Line-strength analysis, 5) Line(s) fitting, 6) Stars and gas kinematics, 7) Stellar populations and SFH. 
Each task is independent from the others and does not modify the spectra.

The 'Preview result' button will display the task(s) result on the selected spectrum in a graphic Matplotlib window and in the output frame on the right. If no task is selected, a warning message will pop-up when clicking the button.

The right frame displays the text output of the software. This is how SPAN communicates with you. This panel reproduces the computer terminal and shows the output of the operations performed, including errors and warnings.


###Apply the tasks###
Once you are satisfied with your work, you can process the spectra or the single selected spectrum. The 'Process selected' button will perform all the tasks activated in the 'Spectra manipulation' panel and in the 'Spectral analysis' frame, saving the new processed spectrum to a fits file. By default, the program will save intermediate spectra if more than one Spectra manipulation task is activated, i.e. one version for each activated task in the Spectra manipulation panel. For example, if you have selected rebinning, sigma broadening, and add noise, the program will save a spectrum with rebinning done, a second spectrum with rebinning + sigma broadening applied, and a third with rebinning + sigma broadening + add noise applied.
If you are not interested in saving all the intermediate spectra files modified in the 'Spectra manipulation' panel, you can select the 'Save final spectra' option at the very bottom of the 'Spectra manipulation' panel, and only the spectrum at the end of the selected tasks (if any) will be saved. This is strongly recommended to do if you are applying more than one task with the reorder option activated. If you are planning to use the tasks of the 'Spectra manipulation' panel just as preparatory phases to the spectral analysis, maybe you do not want to save the processed spectra every time you perform a spectral analysis task. In this case, in the 'Spectra manipulation' panel you can select the option 'Do not save processed spectra'.
WARNING: In the 'Process selected mode', the results of the spectral analysis frame will be written only in the output frame.

By clicking 'Process all', you will apply all the tasks to all the spectra in your list. This is the only way to save the results of the 'Spectral analysis' frame in an ACII file. You can also store the plots generated during the spectral analysis by activating the option 'Save plots' at the very bottom of the SPAN panel. The plots will be saved in high resolution PNG format and stored in the 'plots' subdirectory of the 'SPAN_results' folder.


### The sub programs ###
At the bottom of the 'Utilities' frame, you will find 3 light blue buttons. These are sub-programs that will allow you to: open, create and modify an ASCII file (Text editor), add and modify the keywords in the header of fits files (FITS header editor), and plot the data generated by the Spectral Analysis in the 'Process all' mode.


********************************************************************************************************
*************************************** The input files ************************************************

In order to work properly, the program needs some text files containing information about your data. To see how they must be formatted, please take a look at those coming with SPAN and already set by default in the graphic interface.

IMPORTANT: The text files MUST always have the first line as header, identified by # (e.g. #spectrum)
          
1) Spectra file list task:
    It is essential. If you don't believe it, try to perform any task without upload the spectra and you will see the effects! It is just an ASCII file containing the path (relative if they are in a subfolder of SPAN, absolute if they are elsewhere) and the complete names (with file extension) of the spectra you want to process. You can use any spectra you want, with different format (fits, ASCII...) and resolutions, but it is mandatory to use spectra with the same wavelength units. If you just want to play with one spectrum, then load the ASCII or fits 1D spectrum and activate the option 'I browsed a single spectrum' before clicking the button 'Load!'.

                                            example_list.dat
                                            
                                            #filename ---> header: always necessary!
                                            [path/]spectrum1.fits
                                            [path/]spectrum2.fits
                                            [path/]spectrum3.fits
                                            [path/]spectrum4.fits


Other ASCII files may be needed in the 'Spectra manipulation' panel for some specific tasks. They are:
                                            
2) Doppler correction file for the 'Doppler/z correction' task and the 'I have a list file' option selected:
    It is an ASCII file containing two columns: 1) Name of the spectrum and 2) Radial velocity to correct to the spectrum. This file has the same format of the output text file generated by the Cross-correlation task, so you can directly use it. 
                    
                                            example_dopcor.dat
                                        
                                        #spectrum       RV(km/s) ---> header: always necessary!
                                        [path/]spectrum1.fits  1000
                                        [path/]spectrum2.fits  1001
                                        [path/]spectrum3.fits  1002
                                        [path/]spectrum4.fits  1003
                                            
3) Heliocentric correction file for the 'Heliocentric correction' task and the 'I have a file with location...' option selected: 
    It is an ASCII file containing three columns, separated by a space: 1) Name of the location, 2) Date of the observation (just year, month, day, not the hour), 3) RA of the object (format: degree.decimal), 4) Dec. of the object (format: degree.decimal).
    
                                            example_heliocorr.dat
                                        
                                #where  date        RA          Dec
                                paranal 2016-6-4    4.88375     35.0436389
                                paranal 2016-6-30   10.555      1.11121
                                aao     2011-12-24  -50.034     55.3232
                                aao     2018-2-13   -11.443     11.2323
                                SRT     2020-7-31   70.234      55.32432


Some external files may be needed for specific options of the 'Spectral analysis' tasks. They are:
                            
4) Cross-correlation and velocity dispersion tasks:
    These task require a single template, in fits or ASCII format (ie. just a spectrum!)
    
    
5) Line-strength analysis task and the option 'User indices on a list file' selected:
    It is an ASCII text file containing the index definitions. One index per column. Don't mess it up with the index file, otherwise you will obtain inconsistent results! Luckily, you can always test a single index and see the graphical preview before running the wrong indices on 240913352 spectra and waste one year of your life.
    
                                            example_idx_list_file.dat

                                    #Idx1    Idx2  ---> header: always necessary!
                                    8474   8474 ---> row2: left blue continuum band, in A
                                    8484   8484 ---> row3: right blue continuum band, in A
                                    8563   8563 ---> row4: left red continuum band, in A
                                    8577   8577 ---> row5: right red continuum band, in A
                                    8461   8484 ---> row6: left line limits, in A
                                    8474   8513 ---> row7: right line limits, in A


6) Calculate velocity dispersion coefficients, located in the 'Line-strength parameters' sub-window : 
    It determines 4 spline correction coefficients in order to correct the equivalent width of galactic spectra broadened by the velocity dispersion. It needs a sample of unbroadened spectra that are a good match of the expected stellar populations of the galaxy spectra you want to correct to the zero velocity dispersion frame. The input file is just an ASCII file containing the list of the spectra used as sample. By default, the program has stored the spectra of 31 giant K and early M (<5) stars of the IRTF catalog. This means that this sample is suitable only for the NIR band (8500-24000 A).
Why this wavelength-limited sample? Because the Lick/IDS indices in the optical already have their own correction coefficients, safely stored in SPAN and ready to be used. There is no need to calculate them again! 

                                            example_coeff_determ.dat
                                            
                                            #filename ---> header: always necessary!
                                            [path/]stellar_spectrum1.fits
                                            [path/]stellar_spectrum2.fits
                                            [path/]stellar_spectrum3.fits
                                            [path/]stellar_spectrum4.fits

                                    
7) Correct the line-strength for velocity dispersion task: 
    To apply the velocity dispersion coefficients and correct the raw equivalent widths to the zero velocity dispersion frame, you need this task and three files: 
        1) Sigma list file: a file containing the name of the spectra, the velocity dispersion and the relative uncertainties. It has the same format of the output file generated by the Velocity dispersion task. 
        
                                            example_sigma_vel.dat
                                            
                                            #Spectrum       Sigma(km/s) err ---> Header: always necessary
                                            spectrum_name1  166.2       3.0
                                            spectrum_name2  241.5       3.1
                                            spectrum_name3  335.1       6.2
                                            spectrum_name4  241.5       3.2
        
        2) EW file list to correct: the text file containing the raw equivalent widths you want to correct. It has the same format of the output file generated by the Line-strength measurement task. BE CAREFULL to check that the indices are in the EXACT same order of those you used in the 'Calculate velocity dispersion coefficients' task for the correction coefficient determination.
        
                                            
                                            example_uncorrected_ew.dat
                                            
                            #Spectrum       idx1    idx2    idx3   idx1err idx2err idx3err
                            spectrum_name1  0.27    1.38    3.56    0.01     0.01    0.02
                            spectrum_name2  0.15    1.32    3.43    0.01     0.02    0.02
                            spectrum_name3  0.08    0.75    2.81    0.01     0.02    0.02
                            spectrum_name4  0.14    1.25    3.18    0.01     0.01    0.01

        
        3) Correction coefficients file: it is the output file generated by the 'Calculate velocity dispersion coefficients task'. 
                                    
                                            example_correction_coeff.dat
                        #Pa1          Ca1           Ca2          Pa1e         Ca1e         Ca2e
                        4.3282e-08   1.06712e-08  -2.7344e-09  -5.7463e-09   2.2911e-09   2.8072e-10
                       -2.9602e-05  -1.2012e-05   -3.5782e-07   3.9353e-06  -1.9246e-06  -2.9293e-07
                        0.0017       0.0021        8.5793e-05  -0.0001       0.0004       9.9212e-05
                       -0.0029      -0.0085       -0.0016       0.0053      -0.0003      -0.0002

                                    




##################################################################################################
################################## List of operations you can perform ############################
SPAN can perform many operations on the spectra.

WARNING: All the wavelengths of SPAN are given in A, in air, and all the velocities are in km/s.

Here is a description of the functions:

1) Utilities frame: stand alone frame with action buttons on the right.
    a) Show the header of the selected spectrum = shows the header of the selected spectrum, both fits and ASCII;
    b) Show the wavelength step of the spectrum = shows the step of the selected spectrum;
    c) Estimate the resolution = shows the resolution of the selected spectrum by trying to fit an emission sky line. In The W1 and W2 you should put a small wavelength interval containing a sky line: it's up to you!
    d) Convert the spectrum to = converts the selected spectrum to ASCII of Fits;
    e) Compare spectrum with = Compares the selected spectrum with another one selected by the user;
    f) Convert Flux = converts the flux from frequency to lambda and vice-versa. The buttons 'see plot', 'save one' and 'save all' are active to see and save the results for one or all the spectra;
    g) S/N = measures the Signal to Noise in the selected spectrum, in the W. central wavelength selected by the user. The buttons 'save one' and 'save all' are active to save one or all the SNR calculated for the spectra.

2) Spectra manipulation panel
    a) Spectra pre-processing frame
        a) Cropping = performs a simple cropping of the spectra. If the wavelength window to crop is outside the spectrum, SPAN will ignore the task and will not perform the crop;
        b) Dynamic cleaning = performs a sigma clipping on the spectra. The sigma clip factor, the resolving power of the spectrum and the velocity dispersion (instrumental and/or intrinsic) of the selected spectrum is required in order to perform a better cleaning. For the 'Process all' mode, the option 'R and sigma vel file' is available in order to have R (resolution) and sigma values for all the spectra to be processed. Be VERY careful to use this task with strong emission line spectra;
        c) Wavelet cleaning = performs a wavelet denoise of the spectra. The mean standard deviation of the spectra continuum (sigma) and the number of wavelet layers to consider are required. You don't need to measure it, just try different values. Be careful to not delete the signal;
        d) Filtering and denoising = smooths the spectra by performing some denoising filters: box window moving average, gaussian kernel moving average, low-pass Butterworth filter and band-pass Butterworth filter;
        e) Dopcor/z correction = performs the doppler or z correction of the spectra. Single shot option with user input value of radial velocity (in km/s) or z is available both for one or all the spectra. 'I have a file' option only works with the 'Process all' mode: you need a text file with the spectra name and the recession velocities or z values. This file can be generated by the 'Cross-correlation' task in 'Process all' mode;
        f) Heliocentric correction = performs the heliocentric correction on the spectra. The 'Single' option require a location, that can be selected from the 'loc.list' button (it requires an internet connection the first time!). The other fields are the date in the format YYYY-MM-DD and the RA and Dec. of the observed object (in decimals). In the 'I have a file' option, available only for the 'Process all' mode, a list file with location, date, RA and Dec. coordinates for each object is required.

    b) Spectra processing frame
        a) Rebin = performs a rebin/resample of the spectra in linear wavelength step ('pix.lin' option, with the step in A) and in sigma linear step ('sigma lin.' option, with the sigma step in km/s);
        b) Degrade resolution = degrades the resolution of the spectra from R to R, from R to FWHM and from FWHM to FWHM;
        c) Normalise spectrum to = normalises the spectra to the wavelength provided by the user (in A);
        d) Sigma broadening = broads the spectra by convolving with a gaussian function with the standard deviation provided by the user, in km/s. Remember that the real broadening of the spectra will be the quadratic sum between the broadening and the instrumental sigma of the spectra;
        e) Add noise = adds a random Poisson noise to the spectra with a SNR defined by the user. Remember that the final SNR of the spectra will the sum in quadrature between the added noise and the intrinsic SNR of the spectra;
        f) Continuum modelling = models the continuum shape with two options: 1) Simple filtering of the continuum by reducing the spectrum to a very small resolution (R = 50), and 2) polynomial fitting, with the possibility to masks emission/contaminated regions. Both the continuum models can be divided or subtracted to the original spectrum;


    c) Spectra math
        a) Subtract normalised average = subtracts to the spectra the normalised average made from all the spectra loaded;
        b) Subtract norm. spec. = subtracts to the spectra a normalised spectrum selected by the user;
        c) Add constant = add a constant to the spectra;
        d) Multiply by a constant = multiplies the spectra by a user defined constant value.
        e) Calculate first and second derivatives = automatic calculation of the derivatives of the spectra. This task does not modify the original spectra and the derivative spectra cannot be directly used for spectral analysis.
        f) Average all = averages all the spectra (only available in 'Process selected' mode);
        g) Norm. and average all = normalizes to a common wavelength and average all the spectra (only available in 'Process selected' mode);
        h) Sum all = sums all the spectra (only available in 'Process selected');
        i) Norm. and sum all = Normalizes and sum all the spectra (only available in 'Process selected'). The option 'Use for spec. an.' forces the program to use the result of one of these 4 operations for the following spectral analysis;


2) Spectral analysis. This is the core of SPAN. Each task can be fine-tuned by clicking on the relative parameters button on the right:
    a) Blackbody fitting = performs a fit of the spectrum with Planck's blackbody equation and gives the temperature estimation. It works with any type of spectra but it performs better for stellar spectra, with wide (at least 5000 A) wavelength range;
    b)  Cross-correlation = performs a cross-correlation of the spectra with a user selected template. This uses the function crosscorrRV of pyastronomy. The user can smooth the template to a velocity dispersion value in order to improve the cross-correlation and should identify a narrow region of the spectrum to be cross-correlated (tip: the Calcium triplet lines are the best features);
    c) Velocity dispersion = performs the measurement of the velocity dispersion of the spectra by fitting with a user provided template. Some pre-loaded bands in the visible and NIR are shown but the user can select an independent band. The routine succeeds with strong features (the CaT is the best). It is a little rough but very fast and gives reasonably accurate results;
    d) Line-strength analysis = performs the equivalent width measurement of the spectra, with a single user provided index, with a list of indices or the Lick/IDS system. The results are provided in Angstrom. MonteCarlo simulations are run for the uncertainties estimation. The calculation of the Lick/IDS indices can be personalised in many ways: you can correct for the emission, for the velocity dispersion and the recession velocity. You can also perform a linear interpolation with the SSP models of Thomas et al. 2010, xshooter, MILES and sMILES to retrieve the age, metallicity and alpha-enhancement (not available for the xshooter models) of the stellar populations via linear interpolation or with machinne-learning pre-trained models (Gaussian Process Regression). From the 'Line-strength parameters' window, it is possible also to perform the 'Calculate velocity dispersion coefficients' task. This task broadens a sample of K and early M stars of the IRTF library up to 400 km/s and calculates the deviation of the equivalent width of the index/index file provided in the EW measurement task. It works only by pressing the 'Compute!' button and creates a text file with a third order polynomial curve that fits the behaviour of the broadened index (or indices). The 'Correct the line-strength for velocity dispersion' task performs the correction of the equivalent widths based on the coefficients estimated with the 'Calculate velocity dispersion coefficients' task. It works only by pressing the 'Correct!' button and require an EW measurement files with the same indices in the same order to that considered in the 'Calculate velocity dispersion coefficients'. The output files of the 'Line-strength analysis', 'Calculate velocity dispersion coefficients' and 'Velocity dispersion' are ready to be used for this task, if we are considering the same spectra and indices;
    e) Line(s) fitting = performs the fitting of an user defined line with user defined parameters and a combination of gaussian model the spectral lines and straight line for the continuum. If 'CaT lines' is selected, the task will perform an automatic fitting of the Calcium Triplet lines, assuming they have been previously corrected to the rest frame velocity;
    f) Stars and gas kinematics = uses the known ppxf algorithm of Cappellari et al. 2004 to fit a user defined wavelength region of the spectra with a combination of templates. You can select the template library you prefer among the EMILES, GALAXEV, FSPS and XSHOOTER, the moments to fit, whether fit only the stellar component or also the gas, and whether estimate or not the uncertainties with MonteCarlo simulations. It returns the radial velocity, the velocity dispersion and the higher moments up to H6 (if needed, and a nice plot courtesy of Cappellari);
    g) Stellar populations and SFH = uses ppxf to fit a user defined wavelength region of the spectra with a combination of templates. You can select the template library you prefer among the EMILES, GALAXEV, FSPS, XSHOOTER and sMILES, or add any EMILES custom library. The user can decide whether include the gas emission or not, the reddening and the order of multiplicative and additive polynomials of the fit. The age and metallicity range of the templates can be set. It returns a beautiful plot, the kinematics, the weighted age (in luminosity and mass), metallicity (weighted in luminosity and mass), the M/L, the SFH and saves the best fit template and the emission corrected spectra (if any). Works great in the visible and in the NIR, but this depends on the quality of your spectra.
    

##################################################################################################
################################## The sub-programs ############################

The two light-blue buttons in the upper left corner of SPAN (in the 'Prepare and load spectra' frame) are sub-programs that might help you to generate the 1D spectra needed. Here is how they works:

1) Long-slit extraction: allows the extraction of a single 1D spectrum or a series of 1D spectra from a reduced and wavelength calibrated 2D fits image containing the long-slit spectrum of a source, with dispersion axis along the X-axis and the spatial axis along the Y-axis.
Before proceed to the extraction, you need to load a valid 2D fits image, then you need to:
    1) Open the spectrum and see if everything is ok;
    2) Fit the photometric trace in order to find the maximum along the dispersion axis. You need to set the degree of polynomial curve that will be used to fit the trace and correct the distortion and slope of the spectrum;
    3) Correct the spectrum for distortion and slope using the model trace obtained in the previous step.
Then, you can:
    a) extract and save only one 1D spectrum within the selected Y range (useful for point sources);
    b) extract and save a series of n 1D spectra covering all the spatial axis and obtained by binning contiguous rows in order to reach the desired S/N. A spectra list file ready to be loaded to SPAN is also generated, as well as a text file containing the position of the bins relative to the central region of the galaxy and the S/N. 
    The S/N threshold that you must insert is just a very rough estimation of the real S/N. A good starting value to produce 1D spectra with bins with realistic S/N > 30 is 20. Adjust the SNR Threshold to your preference by looking at the real S/N of the bins.
    The pixel scale parameter is optional. If you set to zero it will not be considered. This option is useful if you have the spectrum of an extended source (e.g. a galaxy) and want to sample different regions.

2) DataCube extraction: following the GIST pipeline standard (Bittner at al., 2019), this sub-program allows you to extract 1D spectra from MUSE and CALIFA DataCubes using the Voronoi binning (Cappellari et al., 2003) or manual binning. It also allows to visualise the DataCube loaded and dynamically create a mask (if needed).

The three light-blue buttons at the bottom of the 'Utilities' frame are sub-programs that might help you in the difficult task of analysing and processing astronomical spectra. They work independently from the main program, so you can also not load spectra if you don't need to perform tasks on them. Here is how they works:

1) Text editor: a simple ASCII file editor where you can create, read or modify ASCII files, included those generated by the SPAN tasks. Some basics operations are available, such find, replace and merge rows;

2) FITS header editor: an header editor to add, remove and save the keywords of fits header files. You can select between: 'Single fits header editor' to work with the keywords of one fits file, 'List of fits header editor' to modify the keywords of a list of fits files, 'Extract keyword from list' to extract and save in an ASCII file one or more keywords from the headers of a list of fits files;

3) Plot data: a sub-program to plot the data generated by the 'Spectral analysis' frame and, in general, all the data stored in ASCII space-separated data. Once you browse for the text file and click the 'Load' button, the program will automatically recognise the column names. Select a name for the x and y axis and plot the data to see them in an IDL style plot.
You can personalise the plot by adding the error bars, set the log scale, add a linear fit (simple fit without considering the uncertainties), set the labels, the range, the font size, size and colours of the markers and decide if visualise the legend or not. You may also save the plot in high resolution PNG image format, in the directory where you run SPAN.
If any error occur, the program will warn you.



##################################################################################################
################################## Tricks in the menu bar ############################

The menu bar was introduced in version 4.5 of SPAN, offering several helpful options to enhance your experience with spectral analysis. Here is a detailed overview of some options that you won't find in the main panel (unless you are using the Android version):

1) File --> Save Parameters...: Allows you to save all parameters and values from the main panel and the various parameter windows of the tasks in a .json file.
This feature is very useful as it enables you to preserve any modifications made to parameters, facilitating the restoration of your session each time you reopen SPAN;
2) File --> Load Parameters...: Allows to load the parameters saved in the .json file. This functionality allows you to resume your work with personalised parameters instead of modifying the default ones every time;
3) File --> Restore Default Parameters: Resets all the parameters to their default values. Useful if numerous parameter modifications during a lengthy session have resulted in issues, allowing you to start from fresh;
4) Edit --> Clear All Tasks: Immediately deactivates all tasks activated during the session, enabling a clean restart;
5) Edit --> Clean Output: Deletes the content of the output window. Particularly useful during extended sessions where the generated output may become quite large.
6) Edit --> Show result folder: shows the location of the 'SPAN_results' folder in case you forgot;
7) Edit --> Change result folder...: create a new 'SPAN_results' folder wherever you want to store the results.

Please, report any bug or comment to daniele.gasparri@gmail.com
Have fun!

Daniele Gasparri
2025-06-05
Greetings from the Canary islands!
