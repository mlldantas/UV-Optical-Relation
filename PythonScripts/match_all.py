#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
    Tables Match (my table with STARLIGHT SDSS DR7 tables)
    @author:  Maria Luiza Linhares Dantas
    @date:    2016.21.07
    @version: 0.0.1

"""

from __future__ import division
import numpy as np


# Main thread
if __name__ == '__main__':

    # Configuring the inputs -------------------------------------------------------------------------------------------
    my_data       = '/home/mldantas/Dropbox/MSc_paper/Programs/Results/output_results.txt'
    lines         = '/home/mldantas/Dropbox/STARLIGHT/lines.txt'
    syn01       = '/home/mldantas/Dropbox/STARLIGHT/SYN01_MALU.txt'
    # syn02       = '/home/mldantas/Dropbox/STARLIGHT/SYN02_MALU.txt'
    # syn03       = '/home/mldantas/Dropbox/STARLIGHT/SYN03_MALU.txt'
    # syn04       = '/home/mldantas/Dropbox/STARLIGHT/SYN04_MALU.txt'
    dn4000_txt    = '/home/mldantas/Dropbox/STARLIGHT/dn4000_MALU.txt'

    # The program itself -----------------------------------------------------------------------------------------------
    # Creating a dictionary --------------------------------------------------------------------------------------------
    my_data      = np.loadtxt(my_data, dtype=object)
    lines_table  = np.loadtxt(lines, dtype=object)
    extinction   = np.loadtxt(syn01, dtype=object)
    dn4000_table = np.loadtxt(dn4000_txt, dtype=object)

    my_dictionary = {}
    for i in range(len(my_data[0, :])):                                         # Converting numpy array into dictionary
        my_dictionary[my_data[0, i]] = np.array(my_data[0 + 1:, i], dtype=str)

    extinction_dictionary = {}
    for j in range(len(extinction[0,:])):
        extinction_dictionary[extinction[0, j]] = np.array(extinction[0 + 1:, j], dtype=str)

    dn4000_dictionary = {}
    for k in range(len(dn4000_table[0, :])):
        dn4000_dictionary[dn4000_table[0, k]] = np.array(dn4000_table[0 + 1:, k], dtype=str)

    lines_dictionary = {}
    for j in range((lines_table[0, :]).size):
        lines_dictionary[lines_table[0, j]] = np.array(lines_table[0 + 1:, j], dtype=str)

    # Reading the data and performing the cross-match ------------------------------------------------------------------
    my_plate     = my_dictionary['plate'].astype(int)
    my_mjd       = my_dictionary['mjd'].astype(int)
    my_fiberid   = my_dictionary['fiberid'].astype(int)
    my_fuv_obs   = my_dictionary['flux_fuv_esc(E-17)'].astype(float)
    my_fuv_synth = my_dictionary['FUV_closest_synth_flux'].astype(float)

    ids_extinction      = extinction_dictionary['SC5-output'].astype(object)
    internal_extinction = extinction_dictionary['A_V'].astype(float)

    dn4000_ids       = dn4000_dictionary['SC5-output_file'].astype(str)
    dn4000_obs_break = dn4000_dictionary['Dn4000(obs)'].astype(float)
    dn4000_syn_break = dn4000_dictionary['Dn4000(syn)'].astype(float)

    lines_plate   = lines_dictionary['plate'].astype(int)
    lines_mjd     = lines_dictionary['mjd'].astype(int)
    lines_fiberid = lines_dictionary['fiberID'].astype(int)

    ## Extinction cross-match ------------------------------------------------------------------------------------------
    extinction_plate   = []
    extinction_mjd     = []
    extinction_fiberid = []
    for a in range(internal_extinction.size):
        extinction_plate_i   = ids_extinction[a][0:4]
        extinction_mjd_i     = ids_extinction[a][5:10]
        extinction_fiberid_i = ids_extinction[a][11:14]
        extinction_plate.append(int(extinction_plate_i))
        extinction_mjd.append(int(extinction_mjd_i))
        extinction_fiberid.append(int(extinction_fiberid_i))
    extinction_plate   = np.array(extinction_plate)
    extinction_mjd     = np.array(extinction_mjd)
    extinction_fiberid = np.array(extinction_fiberid)

    indexes = np.arange(my_plate.size)
    extinction_index = []
    for b in range(extinction_plate.size):
        extinction_index_m = indexes[(my_plate == extinction_plate[b]) * (my_mjd == extinction_mjd[b]) *
                                     (my_fiberid == extinction_fiberid[b])]
        if extinction_index_m.size is 0:
            continue
        extinction_index.append(b)
    my_extinction_plate   = extinction_plate[extinction_index]
    my_extinction_mjd     = extinction_mjd[extinction_index]
    my_extinction_fiberid = extinction_fiberid[extinction_index]
    my_internal_extinction = internal_extinction[extinction_index]

    ## Dn4000 crossmatch -----------------------------------------------------------------------------------------------
    dn4000_plate   = []
    dn4000_mjd     = []
    dn4000_fiberid = []
    for l in range(dn4000_ids.size):
        dn4000_plate_i   = dn4000_ids[l].split('.')[0]
        dn4000_mjd_i     = dn4000_ids[l].split('.')[1]
        dn4000_fiberid_i = dn4000_ids[l].split('.')[2]
        dn4000_plate.append(int(dn4000_plate_i))
        dn4000_mjd.append(int(dn4000_mjd_i))
        dn4000_fiberid.append(int(dn4000_fiberid_i))
    dn4000_plate   = np.array(dn4000_plate)
    dn4000_mjd     = np.array(dn4000_mjd)
    dn4000_fiberid = np.array(dn4000_fiberid)


    dn4000_indexes = np.arange(my_plate.size)
    dn4000_data_index = []
    for m in range(dn4000_plate.size):
        dn4000_data_index_m = dn4000_indexes[(my_plate == dn4000_plate[m]) * (my_mjd == dn4000_mjd[m]) *
                                             (my_fiberid == dn4000_fiberid[m])]
        if dn4000_data_index_m.size is 0:
            continue
        dn4000_data_index.append(m)
    my_dn4000_synth   = dn4000_syn_break[dn4000_data_index]
    my_dn4000_plate   = dn4000_plate[dn4000_data_index]
    my_dn4000_mjd     = dn4000_mjd[dn4000_data_index]
    my_dn4000_fiberid = dn4000_fiberid[dn4000_data_index]


    ## Lines crossmatch ------------------------------------------------------------------------------------------------
    indexes = np.arange(my_plate.size)
    new_index = []
    for i in range(lines_plate.size):
        index = indexes[(my_plate == lines_plate[i]) * (my_mjd == lines_mjd[i]) * (my_fiberid == lines_fiberid[i])]
        if index.size is 0:
            continue
        new_index.append(i)
    my_h_alpha    = lines_dictionary['F_Halpha'].astype(float)[new_index]
    my_ew_h_alpha = lines_dictionary['EW_Halpha'].astype(float)[new_index]
    my_h_beta     = lines_dictionary['F_Hbeta'].astype(float)[new_index]
    my_oiii       = lines_dictionary['F_oiii'].astype(float)[new_index]
    my_nii        = lines_dictionary['F_nii'].astype(float)[new_index]

    # Saving the resulting dataset -------------------------------------------------------------------------------------
    np.savetxt('/home/mldantas/Dropbox/MSc_paper/Programs/Results/my_data_match.csv',
               np.column_stack((my_plate, my_mjd, my_fiberid, my_fuv_obs, my_fuv_synth, my_dn4000_synth,
                                my_internal_extinction, my_h_alpha, my_ew_h_alpha, my_h_beta, my_oiii, my_nii)),
               fmt="%d,%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f", delimiter=',', newline='\n',
               header='plate,mjd,fiber_id,fuv_flux_obs,fuv_flux_syn,dn4000_synth,extinction,H_alpha,EW_H_alpha,H_beta,'
                      'OIII,NII')
