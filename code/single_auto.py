import treecorr

def auto_corr(gal, random, config, weight = "False"):
    '''
    measure the auto-correlation with LS estimator using a 
    gal catalog and a random catalog
    '''
    min_sep = config['min_sep']
    max_sep = config['max_sep']
    units = config['units']
    nbins = config['nbins']
    bin_slop = config['bin_slop']

    if weight == "True":
        gal_cat = treecorr.Catalog(ra = gal["RA"] , dec = gal["DEC"], w = gal["W"], ra_units='deg', dec_units='deg')
    if weight == "False":
        gal_cat = treecorr.Catalog(ra = gal["RA"] , dec = gal["DEC"], ra_units='deg', dec_units='deg')

    ran_cat = treecorr.Catalog(ra = random["RA"] , dec = random["DEC"], ra_units='deg', dec_units='deg')
    
    nn = treecorr.NNCorrelation(min_sep = min_sep, max_sep = max_sep, nbins = nbins, sep_units = units, bin_slop = bin_slop)
    rr = treecorr.NNCorrelation(min_sep = min_sep, max_sep = max_sep, nbins = nbins, sep_units = units, bin_slop = bin_slop)
    dr = treecorr.NNCorrelation(min_sep = min_sep, max_sep = max_sep, nbins = nbins, sep_units = units, bin_slop = bin_slop)
    nn.process(gal_cat, gal_cat, num_threads = 60)
    rr.process(ran_cat, ran_cat, num_threads = 60)
    dr.process(gal_cat, ran_cat, num_threads = 60)
    xi, varxi = nn.calculateXi(rr, dr)
    theta = nn.meanr
    
    return theta, xi, varxi
