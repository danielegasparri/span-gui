import ppxf.sps_util as lib
import ppxf.ppxf_util as util
import numpy as np

class SPSLibWrapper:
    def __init__(self, filename, velscale, fwhm_gal=None, age_range=None, lam_range=None,
                 metal_range=None, norm_range=None, norm_type='mean'):
        """
        Inizializza l'istanza della classe sps_lib del modulo sps_util.
        """
        self.sps_instance = lib.sps_lib(
            filename, velscale, fwhm_gal=fwhm_gal, age_range=age_range,
            lam_range=lam_range, metal_range=metal_range, norm_range=norm_range, norm_type=norm_type
        )

    def get_age_grid_2d(self):
        """
        Restituisce la griglia delle età (`age_grid`) in Gyr.
        """
        self.age_grid = self.sps_instance.age_grid
        return self.age_grid #[:,0]

    def get_metal_grid_2d(self):
        """
        Restituisce la griglia delle metallicità (`metal_grid`).
        """
        self.metal_grid = self.sps_instance.metal_grid
        return self.metal_grid #[:,0]




    def mean_age_metal(self, weights, lg_age= True, quiet=False):
        """
        Compute the weighted ages and metallicities, given the weights returned
        by pPXF. The output population will be light or mass-weighted,
        depending on whether the input is light or mass weights.
        The normalization of the weights is irrelevant as it cancels out.
        """
        # assert weights.ndim == 2, "`weights` must be 2-dim"
        # assert self.age_grid.shape == self.metal_grid.shape == weights.shape, \
        #     "Input weight dimensions do not match"
        self.age_grid = self.get_age_grid_2d()
        self.metal_grid= self.get_metal_grid_2d()
        if lg_age == True:
            lg_age_grid = np.log10(self.age_grid) + 9
            mean_lg_age = np.sum(weights*lg_age_grid)/np.sum(weights)
            mean_age = mean_lg_age

        if lg_age == False:
            lin_age_grid = self.age_grid
            mean_lin_age = np.sum(weights*lin_age_grid)/np.sum(weights)
            mean_age = mean_lin_age

        metal_grid = self.metal_grid
        mean_metal = np.sum(weights*metal_grid)/np.sum(weights)

        if not quiet:
            if lg_age == True:
                print('Weighted lg <Age>: %#.3g' % (mean_age))
                print('Weighted <[M/H]>: %#.3g' % mean_metal)
            if lg_age == False:
                print('Weighted <Age> [Gyr]: %#.3g' % ((mean_age)))
                print('Weighted <[M/H]>: %#.3g' % mean_metal)
        return mean_age, mean_metal


        #Modified by Daniele Gasparri in order to plot either the mean lg ages or the linear mean ages.
    def plot(self, weights, lg_age= True, nodots=False, colorbar=True, **kwargs):

        assert weights.ndim == 2, "`weights` must be 2-dim"
        assert self.age_grid.shape == self.metal_grid.shape == weights.shape, \
            "Input weight dimensions do not match"

        ygrid = self.metal_grid

        if lg_age == True:
            xgrid = np.log10(self.age_grid) + 9
            util.plot_weights_2d(xgrid, ygrid, weights, xlabel="lg Age (dex)",
                             nodots=nodots, colorbar=colorbar, **kwargs)
        if lg_age == False:
            xgrid = self.age_grid
            util.plot_weights_2d(xgrid, ygrid, weights, xlabel="Age (Gyr)",
                             nodots=nodots, colorbar=colorbar, **kwargs)


    def get_age_grid(self):
        # lin_age_grid = 10**(self.age_grid)/ 1e9
        return self.age_grid[:,0]#lin_age_grid[:,0]

    def get_metal_grid(self):
        return self.metal_grid[:,0]
