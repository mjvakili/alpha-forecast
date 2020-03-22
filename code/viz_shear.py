import h5py
import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
plt.switch_backend("Agg")


class shear_summarizer():
    
    def __init__(self, z_list, incomp, shape_noise):
        '''
	z_list: List of redshifts; the lower bounds of tomographic bins
	incomp: True(False); whether the Halpha sampe is incomplete or not
	shape_noise: shape_noise of background sources
	'''
        self.z_list = z_list
	self.incomp = incomp
	self.shape_noise = shape_noise
        shear_clctn = []
	for z in z_list:
	    
	    shear_dict = shear_extractor(z, self.incomp, self.shape_noise)
	    shear_clctn.append(shear_dict)
        
	self.shear_clctn = shear_clctn

	return None

    def plot_shear(self):

        fig, ax = plt.subplots()
        for z_indx, z in enumerate(self.z_list):

	    shear_dict = self.shear_clctn[z_indx]
            theta, xi, xi_cov = shear_dict["theta"], shear_dict["xi"], shear_dict["total_cov"]
	    label = label_decoder(z)
            ax.errorbar(theta, xi, np.diag(xi_cov)**.5,
	          fmt = "o", capsize = 4, elinewidth = 2, label = label)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim([0.9, 120])
        ax.set_ylim([1e-6, .003])
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        plt.tick_params(
          axis='x',          # changes apply to the x-axis
          which='minor',      # both major and minor ticks are affected
          bottom=False,      # ticks along the bottom edge are off
          top=False,         # ticks along the top edge are off
          labelbottom=False)    
        plt.setp(ax, xticks=[1, 10, 100], xticklabels=['1', '10', '100'])
        ax.set_xlabel(r'$\theta$ [arcmin]', fontsize = 20)
        ax.set_ylabel(r'$\gamma_{t}(\theta)$', fontsize = 20)
        if self.incomp == True:
	    title = "Pessimistic lensing measurements"
        elif self.incomp == False:
	    title = "Idealistic lensing measurements"
        plt.title(title)	
	plt.legend(fontsize = 10, frameon = False, shadow = False)
	plt.tight_layout()
        plt.savefig('/home/vakili/public_html/lrg_paper2/shear_completeness_'+str(self.incomp)+'.png', dpi = 600)
        
	return None

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
    total_snr = np.dot(xi, np.linalg.solve(total_cov, xi))
    # jacknife snr
    jk_snr = np.dot(xi, np.linalg.solve(jk_cov, xi))
    # Poisson snr
    shot_snr = np.dot(xi, np.linalg.solve(shot_cov, xi))
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

def matplotlib_settings():
    '''
    configuration of matplotlib settings
    '''
    matplotlib.rcParams['xtick.major.size'] = 7
    matplotlib.rcParams['xtick.labelsize'] = 'x-large'
    matplotlib.rcParams['ytick.major.size'] = 7
    matplotlib.rcParams['ytick.labelsize'] = 'x-large'
    matplotlib.rcParams['xtick.top'] = False
    matplotlib.rcParams['ytick.right'] = False
    matplotlib.rcParams['ytick.direction'] = 'in'
    matplotlib.rcParams['xtick.direction'] = 'in'
    matplotlib.rcParams['font.size'] = 15
    matplotlib.rcParams['figure.figsize'] = [7, 7]
    
    return None

def label_decoder(zmin):
    """figure friendly labels 
    for the redshift bins
    """
    zmax, zs = zmin + 0.1, zmin + 0.1
    label = r"${0:.1f}<z<{1:.1f}, \; z>{2:.1f}$".format(zmin, zmax, zs)

    return label

if __name__ == '__main__':
    
    matplotlib_settings()
    z_list = [0.9, 1.1, 1.2, 1.3, 1.4, 1.5]
    incomp = False
    shape_noise = 0.3
    summarizer = shear_summarizer(z_list, incomp, shape_noise)
    summarizer.plot_shear()
