import numpy as np
from single_cross import cross_corr
from kmeans import gal_loader, shear_loader
import fitsio
import sys
import h5py

from matplotlib import gridspec
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib
matplotlib.rcParams['xtick.major.size'] = 7
matplotlib.rcParams['xtick.labelsize'] = 'x-large'
matplotlib.rcParams['ytick.major.size'] = 7
matplotlib.rcParams['ytick.labelsize'] = 'x-large'
matplotlib.rcParams['xtick.top'] = False
matplotlib.rcParams['ytick.right'] = False
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['font.size'] = 15
matplotlib.rcParams['figure.figsize'] = [7,7]
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.switch_backend('Agg')
graph_dir = 'graphs/'

def gtheta(gal, shear, random, config, out_file_name, NJK):
  
    theta, xi_t, xi_x, npairs, xi_tr, xi_xr, npairs_r = cross_corr(gal, shear, random, config, weight = "False")
    print("Done with the correlation function")
    xi_jk_holder = []
    print("Starting the jacknife resampling")
    for jk in range(NJK):
	   
        mask_random = random["JK_LABEL"] != jk
 	random_jk = {"RA": random["RA"][mask_random], "DEC": random["DEC"][mask_random]}
	mask_gal = gal["JK_LABEL"] != jk
	gal_jk = {"RA": gal["RA"][mask_gal], "DEC": gal["DEC"][mask_gal]}
	mask_shear = shear["JK_LABEL"] != jk
	shear_jk = {"RA": shear["RA"][mask_shear], "DEC": shear["DEC"][mask_shear]}
        theta_jk, xi_t_jk, xi_x_jk, npairs_jk, xi_tr_jk, xi_xr_jk, npairs_r_jk = cross_corr(gal_jk, shear_jk, random_jk, config, weight = "False")
           
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
    vizualise_results(theta, xi_t, xi_tr, cov)
    print("done with vizualisation")
    return None

def vizualise_results(theta, xi, cov):

    plt.figure().set_size_inches(8, 6)
    plt.errorbar(theta, xi_t, np.diag(cov)**.5, fmt = "o", capsize = 5, linewidth = 0)
    plt.errorbar(theta, xi_tr, np.diag(cov)**.5, fmt = "o", capsize = 5, linewidth = 0)
    plt.errorbar(theta, xi_t - xi_tr, np.diag(cov)**.5, fmt = "o", capsize = 5, linewidth = 0)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'$w(\theta)$')
    plt.tight_layout()
    plt.savefig(graph_dir+'gtheta_'+str(zmin)+'.png')
    plt.close()
    
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
    gal, random  = gal_loader(zmin, zmax) 
    shear  = shear_loader(zmin, zmax) 
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
    gtheta(gal, shear, random, config, out_file_name, NJK)
