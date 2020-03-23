import numpy as np
import os.path as path
import pyccl as ccl
import h5py
from scipy.interpolate import CubicSpline
from scipy.interpolate import interp1d

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


def cosmo_model():
    '''
    initializing a default LCDM model with the parameters taken from 
    the flagship database
    '''
    model = ccl.Cosmology(Omega_m = 0.319, Omega_b = 0.049, h = 0.67, sigma8 = 0.83, n_s = 0.96)

    return model

def load_nz(zmin, incomp):
    
    nz_fname = "nz_zmin_{0:.1f}_zmax_{1:.1f}_incomp_{2:s}.csv".format(zmin, zmin+0.1, str(incomp))
    nz_fname = os.path.join('data', nz_fname)
    nz_df = pd.read_csv(nz_fname)
    z, gal_nz, shear_nz = nz_df["z"], nz_df["gal_nz"], nz_df["shear_nz"]
    
    return z, gal_nz, shear_nz

class linear_model():

    def __init__(self, zmin, lmin, lmax, nl, incomp):
        
        self.zmin = zmin
        self.model = cosmo_model()
        self.thetas, self.z, self.nz = angular_range(self.zmin)
        self.ell = np.linspace(lmin, lmax, nl)
        self.incomp = incomp
        self.theta = shear_extractor(zmin = zmini, 
                                     incomp = incomp, 
                                     shape_noise = 0.3)["theta"]
        self.z, self.gal_nz, self.shear_nz = load_nz(zmin, incomp)
    
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
        xi = ccl.correlation(self.model, self.ell, cl_gg, 
                             self.theta/60, corr_type = 'GS', 
                             method = 'Bessel')
        
        return xi
