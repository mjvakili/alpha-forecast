import h5py
import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib


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
          axis='x',         
          which='minor',               
	  bottom=False,     
          top=False,                   
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
	fname = os.path.join(graph_dir, 'shear_completeness_'+str(self.incomp)+'.png')
        plt.savefig(fname, dpi = 600)
        plt.close()

	return None
    
    def plot_noise(self, z_pivot):
        
        if z_pivot in self.z_list:
	    print("dissecting the noise contributions at pivot redshift.{0:.1f}".format(z_pivot))
            fig, ax = plt.subplots()
            pivot_mask = self.z_list == z_pivot
            shear_clctn = self.shear_clctn[pivot_mask]
	    theta = shear_clctn["theta"]
	    jk_noise = np.diag(shear_clctn["jk_cov"])**.5
	    shot_noise = np.diag(shear_clctn["shot_cov"])**.5
	    total_noise = np.diag(shear_clctn["total_cov"])**.5
	    if self.incomp == True:
	        title = "Pessimistic error bars at {0:.1f} < {1:s} < {2:.1f}".format(z_pivot, r"$z_{\rm H\alpha}$", z_pivot+.1)
                color = "C3"
            elif self.incomp == False:
	        title = "Idealistic error bars at {0:.1f} < {1:s} < {2:.1f}".format(z_pivot, r"$z_{\rm H\alpha}$", z_pivot+.1)
                color = "C0" 	
	    ax.plot(theta, total_noise, linewidth = 3, 
		        linestyle = "solid", color = color, 
			label = "sqrt of the diagonal elements of the total covariance")
	    ax.plot(theta, jk_noise, linewidth = 3, 
		        linestyle = "dashed", color = color, 
			label = "sqrt of the diagonal elements of the jacknife covariance")
	    ax.plot(theta, shot_noise, linewidth = 3, 
		        linestyle = "dashdot", color = color, 
			label = "shot-noise")
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim([0.9, 120])
            ax.set_ylim([5e-7, 1e-4])
            ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            plt.tick_params(
                    axis='x',          
                    which='minor',      
		    bottom=False,     
                    top=False,         
		    labelbottom=False)    
            plt.setp(ax, xticks=[1, 10, 100], xticklabels=['1', '10', '100'])
            ax.set_xlabel(r'$\theta$ [arcmin]', fontsize = 20)
            ax.set_ylabel(r'$\sigma [\gamma_{t}](\theta)$', fontsize = 20)
            plt.title(title)	
	    plt.legend(fontsize = 10, frameon = False, shadow = False)
	    plt.tight_layout()
	    fname = os.path.join(graph_dir, 'noise_completeness_'+str(self.incomp)+'.png')
            plt.savefig(fname, dpi = 600)
            plt.close()
	else:
	    print("the pivot redshift needs to be in the provided list")
	    raise ValueError


    def plot_SNR(self):
        '''
	shows the redshift evolution of 
	the signal-to-noise and different 
	contributions to it
	'''
	total_snr, jk_snr, shot_snr = [], [], []
        fig, ax = plt.subplots()
        for z_indx, z in enumerate(self.z_list):

	    shear_dict = self.shear_clctn[z_indx]
            total_snr.append(shear_dict["total_snr"])
            jk_snr.append(shear_dict["jk_snr"])
            shot_snr.append(shear_dict["shot_snr"])
        
	if self.incomp == True:
	    title = "Pessimistic lensing signal-to-noise ratio"
            color = "C3"
        elif self.incomp == False:
	    title = "Idealistic lensing signal-to-noise ratio"
            color = "C0" 	
	ax.plot(np.array(self.z_list) +.05, total_snr, linewidth = 3, linestyle = "solid", color = color, label = "signal-to-noise ratio")
	ax.plot(np.array(self.z_list) +.05, jk_snr, linewidth = 3, linestyle = "dashed", color = color, label = "signal-to-noise ratio jacknife only")
	ax.plot(np.array(self.z_list) +.05, shot_snr, linewidth = 3, linestyle = "dashdot", color = color, label = "signal-to-noise ratio shot-noise only")
        ax.set_xlim([0.85, 1.65])
        ax.set_ylim([8, 170])
        ax.set_yscale('log')
	ax.set_xlabel(r'$z_{\rm H\alpha}$', fontsize = 20)
        ax.set_ylabel(r'S/N $[\gamma_{t}]$', fontsize = 20)
        plt.title(title)	
	plt.legend(fontsize = 10, frameon = False, shadow = False)
	plt.tight_layout()
        fname = os.path.join(graph_dir, 'snr_completeness_'+str(self.incomp)+'.png')
	plt.savefig(fname, dpi = 600)
        plt.close()

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

def matplotlib_settings():
    '''
    configuration of matplotlib settings
    '''
    plt.switch_backend("Agg")
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
    zmax, zs = zmin + 0.1, zmin + 0.2
    label = "{0:.1f}< {3:s} <{1:.1f},  {4:s}>{2:.1f}".format(zmin, zmax, zs, r"$z_{\rm H\alpha}$", r"$z_{\rm s}$")

    return label

if __name__ == '__main__':
    

    matplotlib_settings()
    z_list = [0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
    shape_noise = 0.3
    graph_dir = "/home/vakili/public_html/lrg_paper2"
    for incomp in [False, True]:
        
	summarizer = shear_summarizer(z_list, incomp, shape_noise)
        summarizer.plot_shear()
        summarizer.plot_SNR()
	summarizer.plot_noise(1.0)
