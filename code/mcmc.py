from linear_bias import linear_model, shear_extractor
import os 
import sys
import numpy as np
import emcee
from numpy.linalg import solve
import h5py

def lnPost(theta, **kwargs):

    def lnprior(theta, **kwargs):
        '''log prior 
        '''
        obs = kwargs['data']
        obs_cov = kwargs['data_cov']
        kwargs.pop('data', None)
        kwargs.pop('data_cov', None)
        prior_min = 0.2
        prior_max = 4.0
        
        if (prior_min < theta[0] < prior_max): 
            return 0.0
        else:
            return -np.inf
     
    def lnlike(theta, **kwargs):
        '''log likelihood
        '''
        obs = kwargs['data']
        obs_cov = kwargs['data_cov']
        kwargs.pop('data', None)
        kwargs.pop('data_cov', None)
        # Likelihood
        model_obs = generator(theta)
        res = model_obs - obs
        neg_chisq = -0.5 * np.sum(np.dot(res , np.linalg.solve(obs_cov , res)))
        #print("neg_chi_tot" , neg_chisq*2./len(res))

        return neg_chisq
    
    lp = lnprior(theta , **kwargs)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, **kwargs)

def mcmc_mpi(Nwalkers, Niters, zmin, incomp, shape_noise, chain_file_name): 
    '''
    Parameters
    -----------
    - Nwalker : 
        Number of walkers
    - Nchains : 
        Number of MCMC chains   
    '''
    shear_data = shear_extractor(zmin, incomp, shape_noise)
    Ndim = 1
    random_guess = np.array([1.2])
    #initializing the positions of the walkers
    pos0 = np.repeat(random_guess, Nwalkers).reshape(Ndim, Nwalkers).T + \
                         1.e-4 * np.random.randn(Ndim * Nwalkers).reshape(Nwalkers, Ndim)
    print("initial positions of the walkers = ", pos0)
    #setting up the data kwargs to be passed to log-likelihood
    data_kwargs = { 'data': shear_data["xi"], 
                    'data_cov': shear_data["total_cov"]}
    #setting up the MCMC sampler
    sampler = emcee.EnsembleSampler(Nwalkers, Ndim, lnPost, kwargs = data_kwargs)  
    sampler.sample(pos0, iterations = Niters)
    cnt = 0

    for result in sampler.sample(pos0, iterations = Niters):
        position = list(result)[0]
        sample_file = h5py.File(chain_file_name)
        sample_file["mcmc"][cnt] = position
        sample_file.close()
        print(cnt)
        cnt += 1
        pass
    
    return None

def chain_fname(zmin, incomp):

    chain_file_name = 'mcmc_bias_wtheta_zmin_'+str(zmin)+'_incompleteness_'+str(incomp)+'.hdf5'
    sample_file = h5py.File(chain_file_name , 'w')
    sample_file.create_dataset("mcmc",(Niters, Nwalkers, 1), data = np.zeros((Niters, Nwalkers , 1)))
    sample_file.close()

    return chain_file_name 

if __name__=="__main__": 

    Nwalkers = int(sys.argv[1])
    print('N walkers = ', Nwalkers)
    Niters = int(sys.argv[2])
    print('N iterations = ', Niters)
    zmin = np.float(sys.argv[3])
    print('zmin = ', np.float(zmin))
    #setting the l-range for integration
    lmin, lmax, nl = 10, 10000, 1000
    #incompleteness
    incomp = True
    #shape noise
    shape_noise = 0.3
    # assumed cosmo model
    cosmo_dict = {"Omega_m" : 0.319, "Omega_b" : 0.04, 
                  "sigma8" : 0.83, "h" : 0.67, "n_s" : 0.96}
    #linear bias model on top of the nonlinear cosmo model
    generator = linear_model(zmin, lmin, lmax, nl, incomp, cosmo_dict)
    #MCMC output initialization 
    chain_file_name = chain_fname(zmin, incomp)
    #run the MCMC
    mcmc_mpi(Nwalkers, Niters, zmin, incomp, shape_noise, chain_file_name)
