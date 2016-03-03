SELECT r.Objid, r.plate, r.mjd, r.fiberid, p.ra, p.dec, p.nuv_mag, p.fuv_mag, p.nuv_magerr, p .fuv_magerr, p.e_bv, p.fuv_s2n as sn_fuv_auto, pe.mpstype as survey 
INTO mydb.Galex_results_crosscheck
FROM PhotoObjAll as p, XSDSSDR7 as s, MyDB.starlight_params as r, photoextract as pe
WHERE (p.objid = s.objid AND s.SDSSObjid = r.ObjID AND pe.photoExtractID = p.photoExtractID AND p.nuv_mag > -99 AND (p.fuv_mag > -99) AND p.band=3 AND s.multipleMatchCount = 1 )
ORDER BY r.plate, r.mjd, r.fiberid
