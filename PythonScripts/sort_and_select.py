#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
    Sort program
    @author:  Maria Luiza L. Dantas
    @date:    2016.03.03
    @version: 0.0.1

"""

import numpy as np
import csv

# ======================================================================================================================
# Main thread ----------------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    data = np.loadtxt('MyData_GALEX_SDSS.csv', delimiter=',', dtype=object)

    dictionary = {}
    for i in range(len(data[0, :])):                             # Converting numpy array into dictionary
        dictionary[data[0, i]] = np.array(data[0 + 1:, i], dtype=object)

    myheader = dictionary.keys()

    objectid = dictionary['Objid'].astype(str)
    s2n_fuv  = dictionary['sn_fuv_auto'].astype(float)
    
    # first_column = data[1:, 0].astype(long)
    # lines        = data[1:, :].astype(str)
    # columns      = data[1:, :].astype(str)
    
    table = data[1:, :]
    
    # Criando lista de ids não repetidas -------------------------------------------------------------------------------

    ids = set(list(objectid))

    new_list = []
    for i in ids:
        mask = (table[:, 0] == i)
        s2n_fuv_max = np.argmax(s2n_fuv*mask)
        new_list.append(table[s2n_fuv_max, :])
    new_table = np.array(new_list)
    

np.savetxt('MyData_GALEX_SDSS_clean.csv', new_table, fmt='%s', delimiter=',', newline='\n', header='Objid,plate,mjd,fiberid,ra,dec,fiberMag_z,modelMag_u,modelMag_g,modelMag_r,modelMag_i,modelMag_z,fuv_mag,nuv_mag,dered_u,dered_g,dered_r,dered_i,dered_z,petroR90_r,modelMagErr_u,modelMagErr_g,modelMagErr_r,modelMagErr_i,modelMagErr_z,fuv_magerr,nuv_magerr,e_bv,s2n_r,sn_fuv_auto,survey,morph_type')
    
__author__ = 'mldantas'
