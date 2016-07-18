#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
    Tables Match (my table with STARLIGHT SDSS DR7 tables) and Plots
    @author:  Maria Luiza Linhares Dantas
    @date:    2016.23.06
    @version: 0.0.1

"""

from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import os
import csv
import scipy.stats as s

# ======================================================================================================================

# Main thread
if __name__ == '__main__':

    # Configuring the inputs -------------------------------------------------------------------------------------------
    my_data       = '/home/mldantas/Dropbox/MSc_paper/Programs/Results/output_results.txt'
    lines         = '/home/mldantas/Dropbox/STARLIGHT/lines.txt'
    # syn01       = '/home/mldantas/Dropbox/STARLIGHT/SYN01_MALU.txt'
    # syn02       = '/home/mldantas/Dropbox/STARLIGHT/SYN02_MALU.txt'
    # syn03       = '/home/mldantas/Dropbox/STARLIGHT/SYN03_MALU.txt'
    # syn04       = '/home/mldantas/Dropbox/STARLIGHT/SYN04_MALU.txt'
    dn4000_txt    = '/home/mldantas/Dropbox/STARLIGHT/dn4000_MALU.txt'

    # The program itself -----------------------------------------------------------------------------------------------
    # Creating a dictionary --------------------------------------------------------------------------------------------

    my_data      = np.loadtxt(my_data, dtype=object)
    lines_table  = np.loadtxt(lines, dtype=object)
    dn4000_table = np.loadtxt(dn4000_txt, dtype=object)

    my_dictionary = {}
    for i in range(len(my_data[0, :])):                                         # Converting numpy array into dictionary
        my_dictionary[my_data[0, i]] = np.array(my_data[0 + 1:, i], dtype=str)

    lines_dictionary = {}
    for j in range((lines_table[0, :]).size):
        lines_dictionary[lines_table[0, j]] = np.array(lines_table[0 + 1:, j], dtype=str)

    dn4000_dictionary = {}
    for k in range(len(dn4000_table[0, :])):
        dn4000_dictionary[dn4000_table[0, k]] = np.array(dn4000_table[0 + 1:, k], dtype=str)

    my_plate   = my_dictionary['plate'].astype(int)
    my_mjd     = my_dictionary['mjd'].astype(int)
    my_fiberid = my_dictionary['fiberid'].astype(int)

    lines_plate   = lines_dictionary['plate'].astype(int)
    lines_mjd     = lines_dictionary['mjd'].astype(int)
    lines_fiberid = lines_dictionary['fiberID'].astype(int)

    dn4000_ids     = dn4000_dictionary['SC5-output_file'].astype(str)
    dn4000_obs_break = dn4000_dictionary['Dn4000(obs)'].astype(float)
    dn4000_syn_break = dn4000_dictionary['Dn4000(syn)'].astype(float)

    output_file = open('/home/mldantas/Dropbox/MSc_paper/Programs/Results/matched_data.txt', 'w')

    # Crossmatching my data with the Dn4000 table of STARLIGHT ---------------------------------------------------------
    dn4000_plate = []
    dn4000_mjd = []
    dn4000_fiberid = []
    dn4000_index = []
    dn4000_indexes = np.arange(dn4000_ids.size)
    for l in range(dn4000_ids.size):
        dn4000_plate_i   = dn4000_ids[l].split('.')[0]
        dn4000_mjd_i     = dn4000_ids[l].split('.')[1]
        dn4000_fiberid_i = dn4000_ids[l].split('.')[2]
        dn4000_plate.append(dn4000_plate_i)
        dn4000_mjd.append(dn4000_mjd_i)
        dn4000_fiberid.append(dn4000_fiberid_i)
    dn4000_plate   = np.array(dn4000_plate)
    dn4000_mjd     = np.array(dn4000_mjd)
    dn4000_fiberid = np.array(dn4000_fiberid)

    input_file  = open('/home/mldantas/Dropbox/MSc_paper/Programs/Results/dn4000_matched.txt')
    reader_file = csv.reader(input_file)

    print "{:<5} {:<7} {:<4}".format('plate', 'mjd', 'fiberid')

    indexes = np.arange(my_plate.size)
    dn4000_data_index = []
    for m in range(dn4000_ids.size):
        dn4000_data_index_m = \
            indexes[(my_plate == dn4000_plate[m]) * (my_mjd == dn4000_mjd[m]) * (my_fiberid == dn4000_fiberid[m])]
        if dn4000_data_index is 0:
            continue
        print "{:<5} {:<7} {:<4}".format(dn4000_plate[m], dn4000_mjd[m], dn4000_fiberid[m])
        dn4000_data_index.append(m)
    dn4000_data_index = np.array(dn4000_data_index)
    exit()


    # Crossmatching my data with the lines' table of STARLIGHT ---------------------------------------------------------
    indexes = np.arange(my_plate.size)
    new_index = []
    for i in range(lines_plate.size):
        index = indexes[(my_plate == lines_plate[i]) * (my_mjd == lines_mjd[i]) * (my_fiberid == lines_fiberid[i])]
        if index.size is 0:
            continue
        new_index.append(i)



    output_file.close()


    # Observational Plots ----------------------------------------------------------------------------------------------
    # Fiber Correction histogram ---------------------------------------------------------------------------------------
    plt.hist(my_dictionary['fiber_corr'].astype(float), bins=10, color='#00CED1')
    plt.xlabel(r"$10^{mag(z)_{fiber}-mag(z)_{model}}$", fontsize=25)
    plt.ylabel(r"Counts", fontsize=25)
    plt.tick_params('both', labelsize='20')
    plt.grid(alpha=0.40)
    plt.show()

    # g-r histogram ----------------------------------------------------------------------------------------------------
    mag_u_dered = my_dictionary['dered_u'].astype(float)
    mag_r_dered = my_dictionary['dered_r'].astype(float)
    plt.hist(mag_u_dered-mag_r_dered, bins=40, color='#00CED1')
    plt.xlabel(r"u-r", fontsize=25)
    plt.ylabel(r"Counts", fontsize=25)
    plt.tick_params('both', labelsize='20')
    plt.grid(alpha=0.40)
    plt.show()

    # FUV-NUV histogram ------------------------------------------------------------------------------------------------
    mag_fuv_dered = my_dictionary['fuv_mag'].astype(float)
    mag_nuv_dered = my_dictionary['nuv_mag'].astype(float)
    plt.hist(mag_fuv_dered-mag_nuv_dered, bins=40, color='#00CED1')
    plt.xlabel(r"FUV-NUV", fontsize=25)
    plt.ylabel(r"Counts", fontsize=25)
    plt.tick_params('both', labelsize='20')
    plt.grid(alpha=0.40)
    plt.show()

    # g-r x FUV-r ------------------------------------------------------------------------------------------------------
    mag_g_dered = my_dictionary['dered_g'].astype(float)
    plt.plot(mag_g_dered-mag_r_dered, mag_fuv_dered-mag_r_dered, 'o', color='#00CED1')
    plt.xlabel(r"g-r", fontsize=25)
    plt.ylabel(r"FUV-r", fontsize=25)
    plt.tick_params('both', labelsize='20')
    plt.grid(alpha=0.40)
    plt.show()

    # g-r x FUV-r ------------------------------------------------------------------------------------------------------
    mag_g_dered = my_dictionary['dered_g'].astype(float)
    plt.plot(mag_g_dered-mag_r_dered, mag_fuv_dered, 'o', color='#00CED1')
    plt.xlabel(r"g-r", fontsize=25)
    plt.ylabel(r"FUV", fontsize=25)
    plt.tick_params('both', labelsize='20')
    plt.grid(alpha=0.40)
    plt.show()

    # STARLIGHT output plots -------------------------------------------------------------------------------------------
    # BPT and WHAN diagrams --------------------------------------------------------------------------------------------
    my_h_alpha = lines_dictionary['F_Halpha'].astype(float)[new_index]
    my_h_beta  = lines_dictionary['F_Hbeta'].astype(float)[new_index]
    my_oiii    = lines_dictionary['F_oiii'].astype(float)[new_index]
    my_nii     = lines_dictionary['F_nii'].astype(float)[new_index]

    ## BPT -------------------------------------------------------------------------------------------------------------

    xbpt_01 = np.linspace(-1.8, -0.1, 1000)
    xbpt_02 = np.linspace(-1.8, 0.4, 1000)
    xbpt_03 = np.linspace(-1.8, 0.3, 1000)
    schawinski_x = np.linspace(-0.182, 2.0, 1000)

    ybpt_01 = []
    for j in range(len(xbpt_01)):
        ybpt_01j = 0.61 / (xbpt_01[j] - 0.05) + 1.3      #Kauffman03
        ybpt_01.append(ybpt_01j)
    ybpt_01 = np.array(ybpt_01)

    ybpt_02 = []
    for k in range(len(xbpt_02)):
        ybpt_02k = 0.61 / (xbpt_02[k] - 0.47) + 1.19     #Kewley01
        ybpt_02.append(ybpt_02k)
    ybpt_02 = np.array(ybpt_02)

    ybpt_03 = []
    for n in range(len(xbpt_03)):
        ybpt_03n = (-30.787 + 1.1358 * xbpt_03[n] +0.27297) * np.tanh(5.7409 * xbpt_03) - 31.093   # Stasinska
        ybpt_03.append(ybpt_03n)
    ybpt_03 = np.array(ybpt_03)


    schawinski_y = []
    for l in range(len(schawinski_x)):
        schawinski_yl = 1.05 * schawinski_x[l] + 0.45
        schawinski_y.append(schawinski_yl)
    schawinski_y = np.array(schawinski_y)

    # markerplot = 500 * np.arcsinh((my_dictionary['FUV_closest_synth_flux'].astype(float) -
    #                          my_dictionary['flux_fuv_esc(E-17)'].astype(float))/
    #                         my_dictionary['flux_fuv_esc(E-17)'].astype(float))

    markerplot = 100 * my_dictionary['flux_fuv_esc(E-17)'].astype(float)**2

    # markerplot = 200 * dn4000_syn_break[dn4000_data_index]
    # markerplot = markerplot**2

    plot01 = plt.scatter(np.log10(my_nii/my_h_alpha), np.log10(my_oiii/my_h_beta), s=markerplot, c='#00CED1', alpha=0.3)
    plot02, = plt.plot(xbpt_01, ybpt_01, '--', color='black')
    plot03, = plt.plot(xbpt_02, ybpt_02, '-', color='black')
    plot04, = plt.plot(schawinski_x, schawinski_y, '', color='black')
    plot05, = plt.plot(xbpt_03, ybpt_03,  '-.', color='black')
    plt.legend([plot01, plot02, plot03, plot04, plot05], [r"$\propto \, Dn4000_{synth}$ break", r"Kauffman+03",
                                                          r"Kewley+01", r"Schawinski+07", r"Stasinska+06"], numpoints=1,
               loc='lower left', fontsize=20)
    plt.xlabel(r"$\log \left[NIII\right]/H \alpha $", fontsize=20)
    plt.ylabel(r"$\log \left[OIII\right]/H \beta $", fontsize=20)
    plt.text(-0.8, -0.5, r"Star Forming", fontsize=20)
    plt.text(-0.6, 1.0, r"AGN", fontsize=20)
    plt.text(-0.1, -1.0, r"Composite", fontsize=20)
    plt.text(0.05, 0.35, r"Schawinski+07", fontsize=20)
    plt.xlim([-1.8, 0.4])
    plt.ylim([-1.5, 1.5])
    plt.minorticks_on()
    plt.tick_params('both', labelsize='20')
    plt.grid(alpha=0.5)
    plt.savefig(os.path.join('/home/mldantas/Dropbox/MSc_paper/PaperFigs', 'bpt.pdf'), format='pdf', dpi = 100)
    plt.show()

    plot01 = plt.hexbin(np.log10(my_nii/my_h_alpha), np.log10(my_oiii/my_h_beta))
    plt.xlabel(r"$\log \left[NIII\right]/H \alpha $", fontsize=25)
    plt.ylabel(r"$\log \left[OIII\right]/H \beta $", fontsize=25)
    plt.text(-0.8, -0.5, r"SF", fontsize=25)
    plt.text(-0.6, 1.0, r"AGN", fontsize=25)
    plt.xlim([-1.8, 0.4])
    plt.ylim([-1.5, 1.5])
    plt.minorticks_on()
    plt.tick_params('both', labelsize='30')
    plt.grid(alpha=0.5)
    plt.show()

