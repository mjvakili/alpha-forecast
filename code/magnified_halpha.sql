SELECT `ra_gal`, `dec_gal`, `ra_gal_mag`, `dec_gal_mag`, `observed_redshift_gal`, 
       `x_gal`, `y_gal`, `z_gal`, `vx_gal`, `vy_gal`, `vz_gal`, 
       1./(1. - `kappa`*`kappa` - (`gamma1`*`gamma1` + `gamma2`*`gamma2`)) as `mu` 

FROM flagship_mock_1_8_4_s 

WHERE (x_gal>0. AND y_gal>0. AND z_gal>0.) AND (observed_redshift_gal>=0.8) AND 
      (observed_redshift_gal<=2.8) AND (-2.5 * log10(euclid_nisp_h) -48.6 < 24) AND 
      (logf_halpha_model3_ext + log10(1./(1. - `kappa`*`kappa` - (`gamma1`*`gamma1` + `gamma2`*`gamma2`))) >-15.699) AND 
      (kind = 0)
