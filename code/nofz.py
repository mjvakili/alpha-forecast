import fitsio
import numpy as np
from fitsio import FITS
import os
import pandas as pd
class nz_extractor():
   
    def __init__(self, zmax, obj_type, incomp):
        '''
	zmax (float): maximum redshift of the lens tomographic bin
        obj_type(string): gal or source
	incomp: True or False
	'''

        self.zmax = zmax
        self.obj_type = obj_type
        self.incomp = incomp

        return None

    def fname(self):

        if self.obj_type == "gal":
            
            fname = "gal_zmax_{0:.1f}.fits".format(self.zmax)
            if self.incomp == True:
                fname = "incomplete_"+fname

        elif self.obj_type == "shear":

            fname = "shear_zmax_{0:.1f}.fits".format(self.zmax)
            if self.incomp == True:
                fname = "incomplete_"+fname

        else:
            raise ValueError
        
        fname = os.path.join('data', fname)

        return fname

    def redshift_reader(self):

        redshift = fitsio.read(self.fname(), columns = ["redshift"])["redshift"]
        
        bins = np.linspace(0, 2.5, 500)
        hist, bin_edges = np.histogram(redshift, bins, density = True)
        bin_cens = 0.5*(bin_edges[1:] + bin_edges[:-1])
        
        return bin_cens, hist

class nz_compiler(nz_extractor):

    def __init__(self, zmin_list):

        self.zmin_list = zmin_list

    def nz_saver(self, incomp):

        for z_index, zmin in enumerate(self.zmin_list):

            nz_extractor_gal = nz_extractor(zmax = zmin + 0.1 , obj_type = "gal", incomp = incomp)
            nz_extractor_shear = nz_extractor(zmax = zmin + 0.1 , obj_type = "shear", incomp = incomp)
            
            gal_z, gal_nz = nz_extractor_gal.redshift_reader()
            shear_z, shear_nz = nz_extractor_shear.redshift_reader()
            
            nz_dict = {"z": shear_z, "gal_nz": gal_nz, "shear_nz": shear_nz}
            nz_pd = pd.DataFrame(nz_dict)
            nz_fname = "nz_zmin_{0:.1f}_zmax_{1:.1f}_incomp_{2:s}".format(zmin, zmin+0.1, str(incomp))
            nz_fname = os.path.join('data', nz_fname)
            nz_pd.to_csv(nz_fname)

    def nz_all(self):

        for incomp in [False, True]:

            self.nz_saver(incomp) 

compiler = nz_compiler([0.9, 1.0])
compiler.nz_all()
