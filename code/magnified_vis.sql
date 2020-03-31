SELECT `ra_gal`, `dec_gal`,
`ra_gal_mag`, `dec_gal_mag`,
`kappa`,`gamma1`,`gamma2`
`bulge_angle`, `bulge_axis_ratio`, `disk_angle`, `disk_axis_ratio`
, 1./(1. - `kappa`*`kappa` - (`gamma1`*`gamma1` + `gamma2`*`gamma2`)) as `mu`
, -2.5 * log10(`euclid_vis`)- 48.6 as `mag_vis` 
FROM flagship_mock_1_8_4_s
WHERE (x_gal>0. AND y_gal>0. AND z_gal>0.)
AND (observed_redshift_gal>=1.0 AND observed_redshift_gal<1.1)
AND (-2.5 * log10(euclid_vis) -2.5* log10(1./(1. - `kappa`*`kappa` - (`gamma1`*`gamma1` + `gamma2`*`gamma2`))) - 48.6 < 24.5)
