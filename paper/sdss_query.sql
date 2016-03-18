SELECT res.Objid, res.plate, res.mjd, res.fiberid, res.ra, res.dec, glx.fiberMag_z,
  glx.modelMag_u, glx.modelMag_g, glx.modelMag_r, glx.modelMag_i, glx.modelMag_z, glx.extinction_u,
  glx.extinction_g, glx.extinction_r, glx.extinction_i, glx.extinction_z, res.fuv_mag, res.nuv_mag, glx.dered_u, glx.dered_g, glx.dered_r, glx.dered_i,
  glx.dered_z, glx.petroR90_r, glx.modelMagErr_u, glx.modelMagErr_g, glx.modelMagErr_r, 
  glx.modelMagErr_i, glx.modelMagErr_z, res.fuv_magerr, res.nuv_magerr, res.e_bv, s.z AS redshift,
  s.sn_1 AS s2n_r, res.sn_fuv_auto, res.survey, glx.type AS morph_type  
INTO MyDB.MyData_GALEX_SDSS
FROM MyDB.MyData_GALEX AS res, Galaxy AS glx, SpecObjAll AS s
WHERE glx.nchild = 0 AND (glx.Flags & (0x0000000000040000 + 0x0000000000000002 + 0x0000080000000000)) = 0 
  AND (glx.dered_u BETWEEN 16 AND 25) AND (glx.dered_g BETWEEN 16 AND 25) AND (glx.dered_r BETWEEN 16 AND 22) AND (glx.dered_i BETWEEN 16 AND 25) 
    AND (glx.dered_z BETWEEN 16 AND 25) AND (s.specClass = 2) AND (s.zStatus > 1) AND (s.zConf > 0.9) AND (s.zWarning = 0) AND (res.Objid = glx.objid) 
    AND (glx.objid = s.bestObjID) AND (s.z BETWEEN 0.05 AND 0.075) AND (glx.petroR90_r BETWEEN 1.5 AND 3.5)
ORDER BY res.plate, res.mjd, res.fiberid
