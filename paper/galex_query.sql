SELECT r.Objid, r.plate, r.mjd, r.fiberid, p.ra, p.dec, p.nuv_mag, p.fuv_mag, p.nuv_magerr, p.fuv_magerr, p.e_bv, p.fuv_s2n AS sn_fuv_auto, pe.mpstype AS survey 
INTO MyDB.MyData_GALEX
FROM PhotoObjAll AS p, XSDSSDR7 AS s, MyDB.starlight_params AS r, photoextract AS pe
WHERE (p.objid=s.objid) AND (s.SDSSObjid=r.ObjID) AND (pe.photoExtractID=p.photoExtractID) AND (p.nuv_mag > -99) AND (p.fuv_mag > -99) AND (p.band=3) AND (s.multipleMatchCount=1) AND (pe.mpstype='AIS') AND (p.nuv_fwhm_world < 0.01)
ORDER BY r.plate, r.mjd, r.fiberid
