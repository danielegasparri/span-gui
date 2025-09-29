SPAN: SPectral ANalysis software V7.0
Daniele Gasparri, September 2025

# DataCube extraction #

![Datacube extraction](img/datacube_extraction.png)
*Datacube extraction parameters*


This module allows you to extract a series of n 1D spectra from a 3D FITS image (i.e., a datacube) following the GIST pipeline standard (Bittner et al., 2019, 2021).

The required parameters are similar to those used in the GIST pipeline (and for any datacube extraction in general). The data product generated (i.e. spectra and bin table) are fully compatible with GIST.


## Overview
The essential files and parameters required for this module to function correctly are highlighted in bold. The first step is to load a fully reduced, valid FITS datacube. Given the lack of a standardized convention for FITS keywords and data storage extensions, it is strongly recommended to use the "View datacube" button to verify that the datacube is correctly read.

This module has been tested with MUSE, CALIFA, JWST NIRSpec IFU, and the new WEAVE LIFU data products. Other formats may only be partially supported. If at least the flux values are properly read, the "View datacube" button should display the image, though the wavelength slider in the Matplotlib window may show only a generic "Wavelength index" instead of actual wavelength values.


## Basic parameters for extraction
Once the datacube is correctly loaded, configure the following basic parameters:

- A run name (arbitrary).
- The source redshift (use zero if no de-redshifting is required).
- The wavelength range for extraction. **WARNING**: if de-redshift is needed (i.e. inserted redshift value different than 0), the wavelength range to extract is referred to the equivalent rest-frame range. Example: we have observations in the near-infrared from 9200 to 12500 A from a galaxy located at z = 1.17. If redshift value is set to zero (no de-redshift) we can insert this wavelength range (or a subrange) for extraction, but if we want to perform the de-redshift, we must insert the equivalent rest-frame wavelength range for extraction. In this case, it will be roughly 4200-5800 A. 
- Selecting the extraction routine. In the middle panel, choose the routine for reading and extracting the data. Pre-loaded routines are available for MUSE, CALIFA, WEAVE LIFU, and JWST NIRSpec IFU datacubes. You can write your own extraction routines as .py files following the structure of the available ones and load using the option "using a user defined routine for extraction".
- Configure the zero point for spatial coordinates. This is typically the center of the datacube spatial axes but can also be set to the spaxel coordinates of the galaxy's center. The "View datacube" option is useful for retrieving this information.
- Wavelength range to consider for S/N estimation. If left empy, the S/N for the Voronoi rebinning will be estimated using the provided extraction wavelength range. Otherwise, enter the wavelength range where estimate the S/N. This field has effect only for Voronoi binning mode.


## Applying a spatial mask
Locate the required "Select a FITS mask" field to browse and load a valid mask FITS file. This mask must have the same spatial dimensions as the datacube. If mask is not needed, simply leave this field empty.
If no mask is available but you need it, generate one by clicking "Generate mask". A Matplotlib window will open, displaying the datacube with a slider to navigate across wavelength indices (note: these are indices, not actual wavelengths!). 
To mask specific spaxels (e.g., to exclude contaminated regions and/or bad spaxels), Ctrl+left-click to mask and Ctrl+right-click to unmask. You can also mask/unmask larger areas by clicking and dragging with the Ctrl+left or right mouse button. On touchscreen devices (i.e. Android), masking is performed by tapping and dragging. Unmasking is not yet available. Ensure you remain within the plot boundaries while dragging, or the selection may not be applied.
Once satisfied, close the Matplotlib window. The mask will be automatically saved and loaded in the relative field.
If you open the "Generate mask" window but no masking is needed, simply close the window. An blank mask file will be created but it will not have effect on your data.
Optionally, you can also mask spaxels with low S/N values by specifying a threshold.
**Important**: bad and/or invalid spaxels are automatically masked out by SPAN. 


## Binning
Binning spaxels together to extract spectra with higher signal is a very common practice when dealing with datacubes. SPAN offers three options for binning:
1. **Voronoi binning:** The standard for datacube data. This option bins contiguous spaxels using the Voronoi method (Cappellari et al., 2003) to achieve the desired S/N; 
2. **Elliptical radial binning:** Creates concentric annuli around the selected center with a minimum thickness of dR (in arcsec) and a minimum S/N threshold. The center, orientation and eccentricity of the ellipses can be set manually. In this case, the center is taken as the zero point for spatial coordinates given above ("Origin (in pixel) of the coordinates"), while the position angle (PA, assuming standard astronomical orientation of the cube, i.e. north up, east left) and the eccentricity (q = b/a) are inserted here. If automatic placing is selected ("Auto center"), SPAN will find the photometric center for you from the 2D collapsed image automatically generated and will center the ellipses around this point. If automatic geometry is selected ("Auto ellipses") SPAN will fit the ellipses following the photometric profiles of the 2D collapsed image and will find the best values of PA and q for you. If you have reliable data at hand, please use the manual method, at least for the estimation of the photometric center (the option "View datacube" may be very useful). The coordinates of the generated bins written in the _table.fits (XBIN, YBIN) will be the coordinates of the elliptical bins on the major axis corresponding to the elliptical radius R_flux. If the center of your coordinate system corresponds to the photometric center of the ellipses, these coordinates just give you the distance between the center and the corresponding bin point (signal weighted mean) on the major axis. If the origin of your coordinate system is different than the center of the ellipses, then you must account for the offset to treat XBIN and YBIN as distances. **TIP:** Enable the "Plot radial profile" option in the plot maps sub-program if you want to plot directly the radial profiles of your elliptical binned and analyzed spectra.
3. **Manual binning:** Defines manual bins with mouse multiple selection. The manual bin selection lets you to have the full control on the spaxels to bin. By clicking to the "Manual binning" button, an iterative Matplotlib will appear showing your DataCube. You can select any region to be binned by Ctrl+left-click or dragging on the image. Deselecting can be performed by Ctrl+right-click or dragging. On touchscreen devices (i.e. Android), manual binning regions are selected by tapping and dragging. Deselecting is not yet available. Any region fully separated by the others by at least one spaxel, is considered by SPAN as a single bin. You can draw as many regions as you want. Once done, close the matplotlib window and the binning info will be saved.
4. **Use pre-existing bin info:** This is useful, for example, if you have a datacube split between blue and red arm (e.g. WEAVE LIFU) and you want to sample exactly the same regions, either with the Voronoi or manual binning. You just need to select the folder where the bin info and masking files are stored from the previous run that you want to take as reference for the binning info. SPAN will find the required files in the folder, copy them to the folder where you are saving the new extracted spectra and rename them according to the folder name. Therefore, the extraction routine will use these bin and mask info instead of generating new ones.

You can always preview the results by clicking on the "Preview bins" button. If Voronoi binning is activated, you will see the Voronoi bins color coded by their S/N. If Manual binning is activated, you will see the regions you selected that will be binned. In this case you will see still the single spaxels inside these regions, color coded by their S/N.


## Starting the Extraction
Once all parameters are set, click "Extract!" and grab a coffee! The extraction process may take several minutes for large datacubes (i.e. MUSE).

Upon completion, both GIST-standard and SPAN-standard spectra will be saved.

**IMPORTANT:** Once the extraction has been performed, to prevent accidental overwriting of the files any operation of the panel is locked until you will change the "Name of the run", including the preview. If you want to change any parameters and see the result, you MUST first change the "Name of the run". If you try to perform another extraction without changing the "Name of the run", a terminal message will warn that the results in this folder are already available and no extraction is performed.
