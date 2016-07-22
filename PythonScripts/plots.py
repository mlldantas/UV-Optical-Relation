#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
    Plots from the table generated at match_all.py
    @author:  Maria Luiza Linhares Dantas
    @date:    2016.23.06
    @version: 0.0.1

"""

from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import pylab

# Main thread
if __name__ == '__main__':

    # Configuring the inputs -------------------------------------------------------------------------------------------
    my_data         = np.loadtxt('/home/mldantas/Dropbox/MSc_paper/Programs/Results/output_results.txt', dtype=str)
    my_matched_data = np.loadtxt('/home/mldantas/Dropbox/MSc_paper/Programs/Results/my_data_match.csv', delimiter=',',
                                 dtype=str)

    my_dictionary = {}
    for i in range(len(my_data[0, :])):                                         # Converting numpy array into dictionary
        my_dictionary[my_data[0, i]] = np.array(my_data[0 + 1:, i], dtype=str)

    my_matched_dictionary = {}
    for j in range(len(my_matched_data[0, :])):                                 # Converting numpy array into dictionary
        my_matched_dictionary[my_matched_data[0, j]] = np.array(my_matched_data[0 + 1:, j], dtype=str)

    # # Observational Plots ----------------------------------------------------------------------------------------------
    # # Fiber Correction histogram ---------------------------------------------------------------------------------------
    # plt.hist(my_dictionary['fiber_corr'].astype(float), bins=10, color='#00CED1')
    # plt.xlabel(r"$10^{mag(z)_{fiber}-mag(z)_{model}}$", fontsize=25)
    # plt.ylabel(r"Counts", fontsize=25)
    # plt.tick_params('both', labelsize='20')
    # plt.grid(alpha=0.40)
    # plt.show()
    #
    # # g-r histogram ----------------------------------------------------------------------------------------------------
    # mag_u_dered = my_dictionary['dered_u'].astype(float)
    # mag_r_dered = my_dictionary['dered_r'].astype(float)
    # plt.hist(mag_u_dered-mag_r_dered, bins=40, color='#00CED1')
    # plt.xlabel(r"u-r", fontsize=25)
    # plt.ylabel(r"Counts", fontsize=25)
    # plt.tick_params('both', labelsize='20')
    # plt.grid(alpha=0.40)
    # plt.show()
    #
    # # FUV-NUV histogram ------------------------------------------------------------------------------------------------
    # mag_fuv_dered = my_dictionary['fuv_mag'].astype(float)
    # mag_nuv_dered = my_dictionary['nuv_mag'].astype(float)
    # plt.hist(mag_fuv_dered-mag_nuv_dered, bins=40, color='#00CED1')
    # plt.xlabel(r"FUV-NUV", fontsize=25)
    # plt.ylabel(r"Counts", fontsize=25)
    # plt.tick_params('both', labelsize='20')
    # plt.grid(alpha=0.40)
    # plt.show()
    #
    # # FUV-r histogram --------------------------------------------------------------------------------------------------
    # plt.hist(mag_fuv_dered-mag_r_dered, bins=60, color='#00CED1')
    # plt.xlabel(r"FUV-r", fontsize=25)
    # plt.ylabel(r"Counts", fontsize=25)
    # plt.tick_params('both', labelsize='20')
    # plt.grid(alpha=0.40)
    # plt.show()
    #
    # # g-r x FUV-r ------------------------------------------------------------------------------------------------------
    # mag_g_dered = my_dictionary['dered_g'].astype(float)
    # plt.plot(mag_g_dered-mag_r_dered, mag_fuv_dered-mag_r_dered, 'o', color='#00CED1')
    # plt.xlabel(r"g-r", fontsize=25)
    # plt.ylabel(r"FUV-r", fontsize=25)
    # plt.tick_params('both', labelsize='20')
    # plt.grid(alpha=0.40)
    # plt.show()
    #
    # # g-r x FUV-r ------------------------------------------------------------------------------------------------------
    # mag_g_dered = my_dictionary['dered_g'].astype(float)
    # plt.plot(mag_g_dered-mag_r_dered, mag_fuv_dered, 'o', color='#00CED1')
    # plt.xlabel(r"g-r", fontsize=25)
    # plt.ylabel(r"FUV", fontsize=25)
    # plt.tick_params('both', labelsize='20')
    # plt.grid(alpha=0.40)
    # plt.show()

    # STARLIGHT output plots -------------------------------------------------------------------------------------------
    # BPT and WHAN diagrams --------------------------------------------------------------------------------------------
    my_h_alpha    = my_matched_dictionary['H_alpha'].astype(float)
    my_h_beta     = my_matched_dictionary['H_beta'].astype(float)
    my_oiii       = my_matched_dictionary['OIII'].astype(float)
    my_nii        = my_matched_dictionary['NII'].astype(float)
    my_ew_h_alpha = my_matched_dictionary['EW_H_alpha'].astype(float)

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
        ybpt_03n = (-30.787 + (1.1358 * xbpt_03[n]) + 0.27297) * np.tanh(5.7409 * xbpt_03[n]) - 31.093   # Stasinska
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



    # markerplot = 100 * my_dictionary['flux_fuv_esc(E-17)'].astype(float)**2

    index = np.where(np.log10(my_nii/my_h_alpha)!=0)

    fuv_max = my_dictionary['flux_fuv_esc(E-17)'].astype(float)[index].max()
    markerplot_fuv = []
    for d in range(my_dictionary['flux_fuv_esc(E-17)'].astype(float)[index].size):
        markerplot_d = np.log(1./(my_dictionary['flux_fuv_esc(E-17)'].astype(float)[index][d]/fuv_max))
        markerplot_fuv.append(markerplot_d)
    markerplot_fuv = np.array(markerplot_fuv)

    dn4000_syn_break = my_matched_dictionary['dn4000_synth'].astype(float)[index]
    markerplot = 20 * dn4000_syn_break ** 5

    plt.subplot(1,2,1)
    pylab.colors = markerplot_fuv
    plot01 = plt.scatter(np.log10(my_nii/my_h_alpha)[index], np.log10(my_oiii/my_h_beta)[index], c=markerplot_fuv,
                         s=markerplot, alpha=0.8)
    plot02, = plt.plot(xbpt_01, ybpt_01, ':', color='black')
    plot03, = plt.plot(xbpt_02, ybpt_02, '-', color='black')
    plot04, = plt.plot(schawinski_x, schawinski_y, '--', color='black')
    plot05, = plt.plot(xbpt_03, ybpt_03,  '-.', color='black')
    plot01.set_alpha(0.45)
    plt.legend([plot01, plot02, plot03, plot04, plot05], [r"$\propto \, D_{n}4000_{synth}$ break", r"Kauffman+03",
                                                          r"Kewley+01", r"Schawinski+07", r"Stasinska+06"], numpoints=1,
               loc='lower left', fontsize=20)
    plt.xlabel(r"$\log ([NII]/H{\alpha})$", fontweight='bold', fontsize = 30)
    plt.ylabel(r"$\log (\left[OIII\right]/H \beta) $", fontweight='bold', fontsize=30)
    plt.text(-1.5, -0.5, r"Star Forming", fontsize=25)
    plt.text(-0.8, 1.0, r"AGN", fontsize=25)
    plt.text(-0.8, 0.9, r"Seyfert", fontsize=20)
    plt.text(0.05, 0.0, r"AGN(?)", fontsize=25)
    plt.text(0.05, -0.1, r"LINER", fontsize=20)
    plt.text(-0.25, -1.3, r"Composite", fontsize=25)
    plt.xlim([-1.8, 0.4])
    plt.ylim([-1.5, 1.5])
    plt.minorticks_on()
    plt.tick_params('both', labelsize='30')
    plt.grid(alpha=0.5)

    plt.subplot(1,2,2)
    pylab.colors = markerplot_fuv
    plot01 = plt.scatter(np.log10(my_nii/my_h_alpha)[index], np.log10(my_ew_h_alpha)[index], s = markerplot,
                         c=markerplot_fuv, alpha = 0.7)
    plot02, = plt.plot([-0.4, -0.4], [3.5, 0.5], '-', color='black',  linewidth=1.5)
    plot03, = plt.plot([-0.4, 3.5], [0.78, 0.78], '-', color='black',  linewidth=1.5)
    plt.legend([plot01], [r"$\propto \, D_{n}4000_{synth}$ break"], numpoints=1, loc='lower left', fontsize = 20)
    #plt.axvline(x=-0.4, color='black', linewidth=1.5)
    plt.axhline(y=+0.5, color='black', linewidth=0.5)
    #plt.axhline(y=0.78, xmin=0.635, xmax=1, color='black', linewidth=1.5)
    plt.xlabel(r"$\log ([NII]/H{\alpha})$", fontweight='bold', fontsize = 30)
    plt.ylabel(r"$\log EW(H{\alpha})$", fontweight='bold', fontsize = 30)
    plt.text(-1.5, 0., r"Retired/Passive", fontsize=25)
    plt.text(0.1, 1.0, r"sAGN", fontsize=25)
    plt.text(0.1, 0.6, r"wAGN", fontsize=25)
    plt.text(-1.5, 0.6, r"Star Forming", fontsize=25)
    plt.xlim([-1.8, 0.4])
    plt.ylim([-1.0, 3.5])
    plt.minorticks_on()
    plt.tick_params('both', labelsize='30')
    plt.grid(alpha=0.5)
    cbar_r = plt.colorbar(label=r"Flux $(FUV_{max})$")
    cbar_r.ax.set_ylabel(r"$\log_{e} \left( \frac{1}{\rm Flux (FUV)/ Flux (FUV)_{MAX}} \right)$", fontsize=35)
    cbar_r.ax.tick_params(labelsize=30)
    plt.show()

    # plot01 = plt.hexbin(np.log10(my_nii/my_h_alpha), np.log10(my_oiii/my_h_beta))
    # plt.xlabel(r"$\log \left[NIII\right]/H \alpha $", fontsize=25)
    # plt.ylabel(r"$\log \left[OIII\right]/H \beta $", fontsize=25)
    # plt.text(-0.8, -0.5, r"SF", fontsize=25)
    # plt.text(-0.6, 1.0, r"AGN", fontsize=25)
    # plt.xlim([-1.8, 0.4])
    # plt.ylim([-1.5, 1.5])
    # plt.minorticks_on()
    # plt.tick_params('both', labelsize='30')
    # plt.grid(alpha=0.5)
    # plt.show()

