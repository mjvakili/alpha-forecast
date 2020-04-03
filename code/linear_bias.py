import numpy as np
import os
import pyccl as ccl
import h5py
import pandas as pd
from scipy.interpolate import CubicSpline
from scipy.interpolate import interp1d
from astropy.cosmology import LambdaCDM

def shear_extractor(zmin = 1.0, incomp = True, shape_noise = 0.3):
    """
    extract the shear signal and its covariances 
    per tomographic bin
    """
    zmax = zmin + 0.1
    fname = "gtheta_zmin_{0:.1f}_zmax_{1:.1f}.hdf5".format(zmin, zmax)
    if incomp == True:
        fname = "incomplete_"+fname
    shear_data = h5py.File(fname, "r")
    
    theta = shear_data["theta"][:]
    xi = -1.*shear_data["xi"][:]
    # jacknife covariance without shot noise
    jk_cov = shear_data["cov"][:]
    # Npairs per angular bin
    npairs = shear_data["npairs"][:]
    #Poisson noise
    shot_cov = np.diag(shape_noise**2./npairs)
    #total covariance
    total_cov = jk_cov + shot_cov
    #total signal-to-noise ratio (snr)
    total_snr = np.dot(xi.T, np.linalg.solve(total_cov, xi))**0.5
    # jacknife snr
    jk_snr = np.dot(xi.T, np.linalg.solve(jk_cov, xi))**0.5
    # Poisson snr
    shot_snr = np.dot(xi.T, np.linalg.solve(shot_cov, xi))**0.5
    # every relevant info in a dictionary placeholder 
    shear_dict = {"theta": theta, 
                  "xi" : xi, 
		  "shot_cov": shot_cov,
		  "jk_cov": jk_cov,
		  "total_cov": total_cov,
		  "total_snr": total_snr,
		  "jk_snr": jk_snr,
		  "shot_snr": shot_snr}
    
    return shear_dict		  

def z_angular_mpc(zmean):
    '''returns the minimum angular cut to be applied 
       in order to discard nonlinearities
    '''
    cosmo = LambdaCDM(H0 = 67, Om0 = 0.319, Ode0 = 0.681, Ob0 = 0.049)
    angle = cosmo.arcsec_per_kpc_comoving(zmean).value * 12000./60
    
    return angle


def load_nz(zmin, incomp):
     
    nz_fname = "nz_zmin_{0:.1f}_zmax_{1:.1f}_incomp_{2:s}.csv".format(zmin, zmin+0.1, str(incomp))
    nz_fname = os.path.join('data', nz_fname)
    nz_df = pd.read_csv(nz_fname)
    z, gal_nz, shear_nz = nz_df["z"], nz_df["gal_nz"], nz_df["shear_nz"]
    zmean = np.sum(z * gal_nz)/np.sum(gal_nz)    
    theta_min = z_angular_mpc(zmean)    
 
    return z, gal_nz, shear_nz, theta_min

class cosmo_model:
    def __init__(self, Omega_m = .319, Omega_b = .049, h = .67, sigma8 = .83, n_s = .96):
        '''
        initializing a default LCDM model with the parameters taken from the flagship database
        '''
        self.Omega_m = Omega_m
        self.Omega_b = Omega_b
        self.sigma8 = sigma8
        self.n_s = n_s
        
    def model(self):

        return ccl.Cosmology(self.Omega_m, self.Omega_b, self.h, self.sigma8, self.n_s)


class linear_model(cosmo_model):

    def __init__(self, zmin, lmin, lmax, nl, incomp, cosmo_dict):
        
        self.zmin = zmin
        self.ell = np.linspace(lmin, lmax, nl)
        self.incomp = incomp
        self.theta = shear_extractor(zmin = zmin, 
                                     incomp = incomp, 
                                     shape_noise = 0.3)["theta"]
        self.z, self.gal_nz, self.shear_nz, self.theta_min = load_nz(zmin, incomp)
        print("minimum possible angular scale is ", self.theta_min)
        self.theta = self.theta[self.theta > self.theta_min]
        self.Omega_m = cosmo_dict["Omega_m"]
        self.Omega_b = cosmo_dict["Omega_b"]
        self.h = cosmo_dict["h"]
        self.sigma8 = cosmo_dict["sigma8"]
        self.n_s = cosmo_dict["n_s"]
        super().__init__(self.Omega_m, self.Omega_b, self.h, self.sigma8, self.n_s)        
        self.model = super().model()
    
    def __call__(self, theta):

        return self._sum_stat(theta)

    def _sum_stat(self, theta):
        
        lens = ccl.NumberCountsTracer(self.model, 
                                      dndz = (self.z, self.gal_nz), 
                                      has_rsd = False, 
                                      bias = (self.z, theta*np.ones(len(self.z))))
        source = ccl.WeakLensingTracer(self.model, 
                                       dndz=(self.z, self.shear_nz), 
                                       has_shear=True, 
                                       ia_bias=None)
        
        cl_gm = ccl.angular_cl(self.model, lens, source, self.ell)
        xi = ccl.correlation(self.model, self.ell, cl_gm, 
                             self.theta/60, corr_type = 'gl', 
                             method = 'Bessel')
        
        return xi

if __name__ == '__main__':
    cosmo_dict = {"Omega_m" : 0.319, "Omega_b" : 0.04, "sigma8" : 0.83, "h" : 0.67, "n_s" : 0.96}
    lm = linear_model(0.9, 10, 1000, 1000, False, cosmo_dict)
    print("predicted GGL signal is = ", lm(1.0))
