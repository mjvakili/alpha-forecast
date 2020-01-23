'''
a simple script for generating uniform randoms and 
jackknife regions assuming the KiDS-1000 mask
'''
import fitsio
import numpy as np
import healpy as hp
import kmeans_radec
from kmeans_radec import KMeans, kmeans_sample
from astropy.table import Table

lens = fitsio.read("lens.fits", columns = ["ra_gal", "dec_gal", "observed_redshift_gal"])

ra, dec = lens["ra_gal"], lens["dec_gal"]

ra_min, dec_min, ra_max, dec_max = 0, 0, 90, 90

Nr= 50000000
ran_ra = np.random.uniform(0,360,Nr)
ran_dec = np.degrees(np.arcsin(np.random.uniform(-1,1,Nr)))

ran_mask = (ran_ra > ra_min)&(ran_ra < ra_max)&(ran_dec > dec_min)&(ran_dec < dec_max)
ran_ra, ran_dec = ran_ra[ran_mask], ran_dec[ran_mask]

randoms = {'ra': ran_ra,
           'dec': ran_dec}

coord = np.vstack([randoms['ra'], randoms['dec']]).T
ncen = 100
km = kmeans_sample(coord, ncen, maxiter=30, tol=1.0e-4)

labels = km.find_nearest(coord)

table = Table([coord[:,0], coord[:,1], labels], names=('RA', 'DEC', 'JK_LABEL'))
table.write('flagship_randoms_v2.fits', format='fits')
np.savetxt("flagship_jk_centers_v2.txt", km.centers)
