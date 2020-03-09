import kmeans_radec
import fitsio
import numpy as np
from fitsio import FITS

def gal_preprocess(zmin, zmax):

    lens = fitsio.read("lens.fits", columns = ["ra_gal", "dec_gal", "observed_redshift_gal"])
    ra, dec, z = lens["ra_gal"], lens["dec_gal"], lens["observed_redshift_gal"]
    mask = (z>zmin)&(z<zmax)
    ra, dec, z = ra[mask], dec[mask], z[mask]
    print("length of the catalog after applying the cut", len(ra))
    coord = np.vstack([ra, dec]).T
    centers = np.loadtxt("flagship_jk_centers_v2.txt")
    NJK = centers.shape[0]
    print("Segmentation begins!")
    gal_labels_jk = kmeans_radec.find_nearest(coord, centers)
    print("Done with assigning jacknife labels to galaxies")
    
    gals = {"RA": ra, 
    	    "DEC": dec, 
	    "redshift": z,
	    "JK_LABEL": gal_labels_jk}
    
    fits = FITS('data/gal_zmax_'+str(zmax)+'.fits', 'rw')
    fits.write(gals, names = ["RA", "DEC", "redshift", "JK_LABEL"])
    fits.close()
    
    for jk in range(len(centers)):
    	
        gal_jk_mask = gals["JK_LABEL"] != jk
	gal_jk = {"RA": ra[gal_jk_mask], 
	    "DEC": dec[gal_jk_mask], 
	    "redshift": z[gal_jk_mask]}
        
	fits = FITS('data/gal_zmax_'+str(zmax)+'_jk_'+str(jk)+'.fits','rw')
	fits.write(gal_jk, names = ["RA", "DEC", "redshift"])
        fits.close()
    
    return None

def random_preprocess():

    randoms = fitsio.read('flagship_randoms_v2.fits')
    ra, dec, jk_label = randoms["RA"], randoms["DEC"], randoms["JK_LABEL"]
    centers = np.loadtxt("flagship_jk_centers_v2.txt")
    
    for jk in range(len(centers)):
    	
        jk_mask =  jk_label!= jk
	random_jk = {"RA": ra[jk_mask], 
	          "DEC": dec[jk_mask]} 
 
	fits = FITS('data/random_jk_'+str(jk)+'.fits','rw')
	fits.write(random_jk, names = ["RA", "DEC"])
        fits.close()

    return None	

def shear_preprocess(zmin, zmax):

    shear = fitsio.read("source_v2.fits", columns = ["ra_gal", "dec_gal", "observed_redshift_gal", "gamma1", "gamma2"])
    ra, dec, z = shear["ra_gal"], shear["dec_gal"], shear["observed_redshift_gal"]
    g1, g2 = shear["gamma1"], shear["gamma2"]
    mask = (z>zmax+0.1)
    z = z[mask]
    ra, dec = ra[mask], dec[mask]
    g1, g2 = g1[mask], g2[mask]
    print("length of the catalog after applying the cut", len(ra))
    coord = np.vstack([ra, dec]).T
    centers = np.loadtxt("flagship_jk_centers_v2.txt")
    NJK = centers.shape[0]
    print("Segmentation begins!")
    gal_labels_jk = kmeans_radec.find_nearest(coord, centers)
    print("Done with assigning jacknife labels to galaxies")
    
    gals = {"RA": ra, 
    	"DEC": dec, 
	"gamma1": g1, 
	"gamma2": g2,
	"redshift": z,
	"JK_LABEL": gal_labels_jk} 
    
    fits = FITS('data/shear_zmax_'+str(zmax)+'.fits', 'rw')
    fits.write(gals, names = ["RA", "DEC", "gamma1", "gamma2", "redshift", "JK_LABEL"])
    fits.close()

    for jk in range(len(centers)):
    	
        gal_jk_mask = gals["JK_LABEL"] != jk
	gal_jk = {"RA": ra[gal_jk_mask], 
	    "DEC": dec[gal_jk_mask], 
	    "gamma1": g1[gal_jk_mask], 
	    "gamma2": g2[gal_jk_mask],
	    "redshift": z[gal_jk_mask]}
        
	fits = FITS('data/shear_zmax_'+str(zmax)+'_jk_'+str(jk)+'.fits','rw')
	fits.write(gal_jk, names = ["RA", "DEC", "gamma1", "gamma2", "redshift"])
        fits.close()
	
    return None

if __name__  == '__main__':
    
    base_z = 0.9
    delta_z = 0.1
    random_preprocess()
    for zmin in base_z + delta_z*np.arange(9):

        zmax = zmin + delta_z
        gal_preprocess(zmin, zmax)
	print("done with gal segmentation")
        shear_preprocess(zmin, zmax)
	print("done with shear segmentation")
