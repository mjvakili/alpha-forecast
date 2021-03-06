import numpy as np
from single_cross import cross_corr
import fitsio
import sys
import h5py
import kmeans_radec
import fitsio
import numpy as np
from fitsio import FITS

def gtheta(config, zmax, out_file_name, NJK):

    gal = fitsio.read("data/gal_zmax_"+str(zmax)+".fits",
    		      columns = ["RA", "DEC"])
    shear = fitsio.read("data/shear_zmax_"+str(zmax)+".fits",
    		      columns = ["RA", "DEC", "gamma1", "gamma2"])
    random = fitsio.read('flagship_randoms_v2.fits')

    print("Done with reading the full-sky catalogs")
    theta, xi_t, xi_x, npairs, xi_tr, xi_xr, npairs_r = cross_corr(gal, shear, random, config)
    print("Done with the correlation function")
    print("Done with closing the random, gal, and shear files")
    
    xi_jk_holder = []
    print("Starting the jacknife resampling")
    for jk in range(NJK):
    	
	gal_jk = fitsio.read("data/gal_zmax_"+str(zmax)+"_jk_"+str(jk)+".fits",
    			  columns = ["RA", "DEC"])
    	shear_jk = fitsio.read("data/shear_zmax_"+str(zmax)+"_jk_"+str(jk)+".fits",
    		            columns = ["RA", "DEC", "gamma1", "gamma2"])
    	random_jk = fitsio.read("data/random_jk_"+str(jk)+".fits",
    		            columns = ["RA", "DEC"])
	
	#compute the cross-correlation for the jackknife region
	theta_jk, xi_t_jk, xi_x_jk, npairs_jk, xi_tr_jk, xi_xr_jk, npairs_r_jk = cross_corr(gal_jk, shear_jk, random_jk, config)
           
	xi_jk_holder.append(xi_t_jk - xi_tr_jk)
	print("done with jk resampling = ", jk)
    	   
    xi_jk_holder = np.array(xi_jk_holder)
    if xi_jk_holder.shape[0] < NJK:
      xi_jk_holder = np.zeros((NJK, config["nbins"]))
    cov  =  ((NJK - 1)**2./NJK)*np.cov(xi_jk_holder.T)
    corr_file = h5py.File(out_file_name, "w")
    corr_file["theta"] = theta
    corr_file["xi"] = xi_t - xi_tr
    corr_file["npairs"] = npairs
    corr_file["cov"] = cov
    corr_file["xi_samples"] = xi_jk_holder
    corr_file.close()
    print("done with saving the results in the output file") 
    
    return None

if __name__ == '__main__':
    
    zmin = np.float(sys.argv[1])
    print('minimum redshift = ', zmin)
    zmax = np.float(sys.argv[2])
    print('maximum redshift  = ', zmax)
    nbins = int(sys.argv[3])
    print('number of angular bins', nbins)
    tmin = np.float(sys.argv[4])
    print('minimum angle = ', tmin)
    tmax = np.float(sys.argv[5])
    print('maximum angle  = ', tmax)
    NJK = int(sys.argv[6])
    print('number of jacknife resampling  = ', NJK)
    
    print("done with loading the cats")
    config = {'min_sep': tmin, 
    	      'max_sep': tmax, 
	      'units': 'arcmin', 
	      'nbins': nbins, 
	      'bin_slop': 0.001}

    out_file_name = 'gtheta_zmin_'+str(zmin)+'_zmax_'+str(zmax)+'.hdf5'
    sample_file = h5py.File(out_file_name , 'w')
    sample_file.create_dataset("theta", (nbins, ), data = np.zeros((nbins)))
    sample_file.create_dataset("xi", (nbins, ), data = np.zeros((nbins)))
    sample_file.create_dataset("npairs", (nbins, ), data = np.zeros((nbins)))
    sample_file.create_dataset("cov", (nbins, nbins), data = np.zeros((nbins, nbins)))
    sample_file.create_dataset("xi_samples", (nbins, 100), data = np.zeros((nbins, NJK)))
    sample_file.close()
    print("done with output initialization")
    gtheta(config, zmax, out_file_name, NJK)
