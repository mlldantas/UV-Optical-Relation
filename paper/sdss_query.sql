SELECT res.Objid, res.plate, res.mjd, res.fiberid, res.ra, res.dec, glx.modelMag_u, glx.modelMag_g, glx.modelMag_r, glx.modelMag_i, glx.modelMag_z,res.fuv_mag, res.nuv_mag, glx.dered_u, glx.dered_g, glx.dered_r, glx.dered_i, glx.dered_z, glx.petroR90_r, res.fuv_magerr, res.nuv_magerr, res.e_bv, s.sn_1 as s2n_r, s.z AS redshift, res.sn_fuv_auto, res.survey, glx.type as morph_type
INTO MyDB.MyData_GALEX_SDSS
FROM MyDB.MyData_GALEX AS res, Galaxy AS glx, SpecObj AS s
WHERE glx.nchild = 0 AND (glx.Flags & (0x0000000000040000 + 0x0000000000000002 + 0x0000080000000000) ) = 0 AND glx.dered_r<=22.0 AND (glx.dered_u between 5 AND 30) AND (glx.dered_g between 5 AND 30) AND glx.dered_r > 5 AND (glx.dered_i between 5 AND 30 AND (glx.dered_z between 5 AND 30) AND s.specClass = 2 AND s.zStatus > 1 AND s.zConf > 0.9 AND s.zWarning = 0 AND res.Objid = glx.objid AND glx.objid = s.bestObjID AND s.z BETWEEN 0.05 AND 0.075 AND (glx.petroR90_r BETWEEN 1.5 AND 3.5))
ORDER BY res.plate, res.mjd, res.fiberid
