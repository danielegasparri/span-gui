SPAN for Android has been tested using the latest (as of September 2025) Pydroid3 app (the free version is fine) with Android 11-14 OS on Samsung and Xiaomi devices. 
Please, download the latest Pydroid3 app and the Pydroid3 repository (they are two separate apps) from Google Play.
DO NOT install span-gui with pip. Copy the content of the folder to your device (usually in Downloads or Documents).

Before compiling, you need to manually install with the embedded pip of Pydroid3 the following modules, in this order:
  numpy
  astropy
  pandas
  matplotlib
  scipy
  scikit-image
  PyWavelets
  joblib
  scikit-learn
  ppxf
  vorbin
  tk
  certifi
  emcee

Once done, you need to compile the __main__.py.
Put your mobile device in landscape mode (horizontal), otherwise the GUI panel will be truncated. 
Enjoy!

If you experiment an error during the installation of the required modules, try to open the emulated terminal of Pydroid3 and install cython as follow:
pip3 install cython
Then proceed with the installation of the required modules via pip.

WARNING: SPAN for Android has been adapted to meet the touch capabilities of Android based devices, but some interactive functions in the Preview frame of the main GUI are limited. This does not prevent you to fully exploit the analysis capabilities of SPAN, but for a full experience you may find useful to connect a trackpad or a mouse to your device. 

IMPORTANT: due to different screen resolutions of Android devices, you may need to adjust the scaling of the GUI layout if it doesn't fit you screen (i.e. too big). In this case you should open the misc.py module in the span_modules folder, locate the def get_layout function (line 118) and modify the scaling factor for Android at line 153: scale_win = 2.25. Change this factor until you get a pleasant view on your device.  