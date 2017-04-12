## This a schema browser for the fields in the catalogue (output_results.csv).

| Column                 | Type                  | Description                                     |
|------------------------|-----------------------|-------------------------------------------------|
| Objid                  | 1590 non-null int64   | SDSS DR7 Object ID                              |
| plate                  | 1590 non-null int64   | Spectroscopic Plate SDSS DR7                    |
| mjd                    | 1590 non-null int64   | Modified Jullian Date SDSS DR7                  |
| fiberid                | 1590 non-null int64   | Fiber ID SDSS DR7                               |
| ra                     | 1590 non-null float64 | Right ascenscion J2000                          |
| dec                    | 1590 non-null float64 | Declination J2000                               |
| fiberMag_z             | 1590 non-null float64 | Fiber correction calculated using the *z* band  |
| modelMag_u             | 1590 non-null float64 | Model mag from SDSS DR7 - u band                |
| modelMag_g             | 1590 non-null float64 | Model mag from SDSS DR7 - g band                |
| modelMag_r             | 1590 non-null float64 | Model mag from SDSS DR7 - r band                |
| modelMag_i             | 1590 non-null float64 | Model mag from SDSS DR7 - i band                |
| modelMag_z             | 1590 non-null float64 | Model mag from SDSS DR7 - z band                |
| extinction_u           | 1590 non-null float64 | u band extinction from SDSS DR7                 | 
| extinction_g           | 1590 non-null float64 | g band extinction from SDSS DR7                 |
| extinction_r           | 1590 non-null float64 | r band extinction from SDSS DR7                 |
| extinction_i           | 1590 non-null float64 | i band extinction from SDSS DR7                 |
| extinction_z           | 1590 non-null float64 | z band extinction from SDSS DR7                 |
| fuv_mag                | 1590 non-null float64 | FUV magnitude from GALEX DR6 Plus 7             |
| nuv_mag                | 1590 non-null float64 | NUV magnitude from GALEX DR6 Plus 7             |
| dered_u                | 1590 non-null float64 | SDSS magnitude corrected by extinction - u band |
| dered_g                | 1590 non-null float64 |             |
| dered_r                | 1590 non-null float64 |             |
| dered_i                | 1590 non-null float64 |             |
| dered_z                | 1590 non-null float64 |             |
| petroR90_r             | 1590 non-null float64 |             |
| modelMagErr_u          | 1590 non-null float64 |             |
| modelMagErr_g          | 1590 non-null float64 |             |
| modelMagErr_r          | 1590 non-null float64 |             |
| modelMagErr_i          | 1590 non-null float64 |             |
| modelMagErr_z          | 1590 non-null float64 |             |
| fuv_magerr             | 1590 non-null float64 |             |
| nuv_magerr             | 1590 non-null float64 |             |
| e_bv                   | 1590 non-null float64 |             |
| redshift               | 1590 non-null float64 |             |
| s2n_r                  | 1590 non-null float64 |             |
| sn_fuv_auto            | 1590 non-null float64 |             |
| survey                 | 1590 non-null object  |             |
| morph_type             | 1590 non-null int64   |             |
| mag_ab_u               | 1590 non-null float64 |             |
| mag_ab_g               | 1590 non-null float64 |             |
| mag_ab_r               | 1590 non-null float64 |             |
| mag_ab_i               | 1590 non-null float64 |             |
| mag_ab_z               | 1590 non-null float64 |             |
| av_mw                  | 1590 non-null float64 |             |
| a_lambda_u             | 1590 non-null float64 |             |
| a_lambda_g             | 1590 non-null float64 |             |
| a_lambda_r             | 1590 non-null float64 |             |
| a_lambda_i             | 1590 non-null float64 |             |
| a_lambda_z             | 1590 non-null float64 |             |
| a_lambda_fuv           | 1590 non-null float64 |             |
| a_lambda_nuv           | 1590 non-null float64 |             |
| fiber_corr             | 1590 non-null float64 |             |
| flux_u(E-17)           | 1590 non-null float64 |             |
| flux_g(E-17)           | 1590 non-null float64 |             |
| flux_r(E-17)           | 1590 non-null float64 |             |
| flux_i(E-17)           | 1590 non-null float64 |             |
| flux_z(E-17)           | 1590 non-null float64 |             |
| flux_fuv(E-17)         | 1590 non-null float64 |             |
| flux_nuv(E-17)         | 1590 non-null float64 |             |
| flux_u_err(E-17)       | 1590 non-null float64 |             |
| flux_g_err(E-17)       | 1590 non-null float64 |             |
| flux_r_err(E-17)       | 1590 non-null float64 |             |
| flux_i_err(E-17)       | 1590 non-null float64 |             |
| flux_z_err(E-17)       | 1590 non-null float64 |             |
| flux_fuv_err(E-17)     | 1590 non-null float64 |             |
| flux_nuv_err(E-17)     | 1590 non-null float64 |             |
| flux_u_esc(E-17)       | 1590 non-null float64 |             |
| flux_g_esc(E-17)       | 1590 non-null float64 |             |
| flux_r_esc(E-17)       | 1590 non-null float64 |             |
| flux_i_esc(E-17)       | 1590 non-null float64 |             |
| flux_z_esc(E-17)       | 1590 non-null float64 |             |
| flux_fuv_esc(E-17)     | 1590 non-null float64 |             |
| flux_nuv_esc(E-17)     | 1590 non-null float64 |             |
| scaling_factor         | 1590 non-null float64 |             |
| chi^2(goodness)        | 1590 non-null float64 |             |
| absolute_goodness      | 1590 non-null float64 |             |
| Mag_Abs_u              | 1590 non-null float64 |             |
| Mag_Abs_g              | 1590 non-null float64 |             |
| Mag_Abs_r              | 1590 non-null float64 |             |
| Mag_Abs_i              | 1590 non-null float64 |             |
| Mag_Abs_z              | 1590 non-null float64 |             |
| NUV_closest_synth_flux | 1590 non-null float64 |             |
| FUV_closest_synth_flux | 1590 non-null float64 |             |
