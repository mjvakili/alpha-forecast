import treecorr

def cross_corr(gal, shear, random, config):
    '''
    measure the cross-correlation between a gal catalog and a shear catalog
    '''
    min_sep = config['min_sep']
    max_sep = config['max_sep']
    units = config['units']
    nbins = config['nbins']
    bin_slop = config['bin_slop']

    gal_cat = treecorr.Catalog(ra = gal["RA"] , dec = gal["DEC"], ra_units='deg', dec_units='deg')
    shear_cat = treecorr.Catalog(ra = shear["RA"] , dec = shear["DEC"], g1 = shear["gamma1"], g2 = -1.*shear["gamma2"], ra_units='deg', dec_units='deg')
    ran_cat = treecorr.Catalog(ra = random["RA"] , dec = random["DEC"], ra_units='deg', dec_units='deg')
    ng = treecorr.NGCorrelation(min_sep = min_sep, max_sep = max_sep, nbins = nbins, sep_units = units, bin_slop = bin_slop)
    rg = treecorr.NGCorrelation(min_sep = min_sep, max_sep = max_sep, nbins = nbins, sep_units = units, bin_slop = bin_slop)
    ng.process(gal_cat, shear_cat, num_threads = 60)
    rg.process(ran_cat, shear_cat, num_threads = 60)
    theta, xi_t , xi_x , w , npairs = ng.meanr, ng.xi, ng.xi_im, ng.weight, ng.npairs
    theta_r, xi_tr , xi_xr , wr , npairs_r = rg.meanr, rg.xi, rg.xi_im, rg.weight, rg.npairs
    
    return theta, xi_t, xi_x, npairs, xi_tr, xi_xr, npairs_r
