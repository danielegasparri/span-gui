"""
  Copyright (C) 2014-2021, Michele Cappellari
  E-mail: michele.cappellari_at_physics.ox.ac.uk
  http://purl.org/cappellari

  V1.0.0: Created to emulate my IDL procedure with the same name.
        Michele Cappellari, Oxford, 28 March 2014
  V1.1.0: Included reversed colormap. MC, Oxford, 9 August 2015
  V1.1.1: Register colormaps in Matplotlib. MC, Oxford, 29 March 2017
  V1.2.0: Do not re-register colormap. Start x coordinate from zero.
        Reduced numbers duplication to make it easier to modify.
        MC, Oxford, 14 September 2021

"""
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

def register_sauron_colormap():
    if 'sauron' in mpl.colormaps:
        return

    x = np.array([0, 42.5, 85, 105, 117.5, 127.5, 137.5, 150, 170, 212.5, 255]) / 255
    r = [0.0, 0.0, 0.4,  0.5, 0.3, 0.0, 0.7, 1.0, 1.0, 1.0, 0.9]
    g = [0.0, 0.0, 0.85, 1.0, 1.0, 0.9, 1.0, 1.0, 0.85, 0.0, 0.9]
    b = [0.0, 1.0, 1.0,  1.0, 0.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9]

    cdict = {
        'red':   list(zip(x, r, r)),
        'green': list(zip(x, g, g)),
        'blue':  list(zip(x, b, b))
    }

    sauron = LinearSegmentedColormap('sauron', segmentdata=cdict)
    sauron_r = LinearSegmentedColormap('sauron_r', segmentdata={
        'red':   list(zip(x, r[::-1], r[::-1])),
        'green': list(zip(x, g[::-1], g[::-1])),
        'blue':  list(zip(x, b[::-1], b[::-1]))
    })

    mpl.colormaps.register(name='sauron', cmap=sauron)
    mpl.colormaps.register(name='sauron_r', cmap=sauron_r)

##############################################################################

# Usage example for the SAURON colormap.

if __name__ == '__main__':

    n = 41 
    x, y = np.ogrid[-n:n, -n:n]
    img = x**2 - 2*y**2

    register_sauron_colormap()
    
    plt.clf()

    plt.subplot(121)
    plt.imshow(img, cmap='sauron')
    plt.title("SAURON colormap")

    plt.subplot(122)
    plt.imshow(img, cmap='sauron_r')
    plt.title("reversed colormap")
