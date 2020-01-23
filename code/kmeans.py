import kmeans_radec
import fitsio
import numpy as np

def gal_loader(zmin, zmax):

    lens = fitsio.read("lens.fits", columns = ["ra_gal", "dec_gal", "observed_redshift_gal"])
    ra, dec, z = lens["ra_gal"], lens["dec_gal"], lens["observed_redshift_gal"]
    mask = (z>zmin)&(z<zmax)
    ra, dec = ra[mask], dec[mask]
    print("length of the catalog after applying the cut", len(ra))
    coord = np.vstack([ra, dec]).T
    centers = np.loadtxt("flagship_jk_centers_v2.txt")
    NJK = centers.shape[0]
    print("Segmentation begins!")
    gal_labels_jk = kmeans_radec.find_nearest(coord, centers)
    print("Done with assigning jacknife labels to galaxies")
    
    gals = {"RA": ra, 
    	    "DEC": dec, 
	    "JK_LABEL": gal_labels_jk}
    
    randoms = fitsio.read('flagship_randoms_v2.fits') 
    
    return gals, randoms

def shear_loader(zmin, zmax):

    shear = fitsio.read("source.fits", columns = ["ra_gal", "dec_gal", "observed_redshift_gal"])
    ra, dec, z = shear["ra_gal"], shear["dec_gal"], shear["observed_redshift_gal"]
    mask = (z>zmax)
    ra, dec = ra[mask], dec[mask]
    print("length of the catalog after applying the cut", len(ra))
    coord = np.vstack([ra, dec]).T
    centers = np.loadtxt("flagship_jk_centers_v2.txt")
    NJK = centers.shape[0]
    print("Segmentation begins!")
    gal_labels_jk = kmeans_radec.find_nearest(coord, centers)
    print("Done with assigning jacknife labels to galaxies")
    
    gals = {"RA": ra, 
    	    "DEC": dec, 
	    "JK_LABEL": gal_labels_jk}
    
    return gals
