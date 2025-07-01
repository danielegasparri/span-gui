SPAN for Android has been tested using the Pydroid3 app (the free version is fine).
Please, download the latest Pydroid3 app and the Pydroid3 repository (they are two separate apps) from Google Play.
DO NOT install span-gui with pip. Copy the content of the folder to your device (usually in Download or Documents), open the __main__.py in Pydroid3 and add the tk dependency in __main__.py: 
 import tkinter as tk

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

Once done, you need to compile the __main__.py.
Put your mobile device in landscape mode (horizontal), otherwise the GUI panel will be truncated. 
Enjoy!

If you experiment an error during the installation of the required modules, try to open the emulated terminal of Pydroid3 and install cython as follow:
pip3 install cython
Then proceed with the installation of the required modules via pip