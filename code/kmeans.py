import kmeans_radec
import fitsio

def gal_loader(zmin, zmax):

    gal = fitsio.read("lens.fits", columns = ["ra_gal", "dec_gal", "observed_redshift_gal"])
    ra, dec, z = lens["ra_gal"], lens["dec_gal"], lens["observed_redshift_gal"]
    coord = np.vstack([ra, dec]).T
    centers = np.loadtxt("flagship_jk_centers_v2.txt")
    NJK = centers.shape[0]
    gal_labels_jk = kmeans_radec.find_nearest(coord, centers)
    print("Done with assigning jacknife labels to galaxies")
    
    gals = {"RA": ra, 
    	    "DEC": dec, 
	    "JK_LABEL": gal_labels_jk}
    
    randoms = fitsio.read('flagship_randoms_v2.fits') 
    
    return gals, randoms
