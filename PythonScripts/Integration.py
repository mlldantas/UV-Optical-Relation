#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
    Data Reduction and Adjustments
    @author:  Maria Luiza Linhares Dantas
    @date:    2016.03.10
    @version: 0.1.9
    This program uses synthetic spectra produced by STARLIGHT and adjusts observational data using chi^2 minimization
    method.

    This version takes into account the fiber correction in the flux calculation, and also takes into account SDSS'
    deredenned
    magnitudes and Fitzpatrick's Law to correct the UV magnitudes, and filter responses.

    In this version we completed the loop in order to analyze all the data and implement a single output.
    For more information about the input data, please read the schema text file.
"""

from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import os
import cosmocalc as cosmo
import scipy.stats as s
import csv
from fitzpatrick import fm_unred
from cardelli import cardelli_redlaw

# ======================================================================================================================

# Main thread
if __name__ == '__main__':

    # Configuring the inputs -------------------------------------------------------------------------------------------

    results_path      = '/home/mldantas/Dropbox/MSc_paper/Programs/Fig'
    spec_files_path   = '/home/mldantas/Dropbox/MSc_paper/Programs/EspectrosRebinadosSelfZ'
    data_obs          = '/home/mldantas/Dropbox/MSc_paper/Dataset/MyData_GALEX_SDSS_clean.csv'
    #data_obs          = '/home/mldantas/Dropbox/MSc_paper/Dataset/new_teste.txt'
    filelist          = '/home/mldantas/Dropbox/MSc_paper/Programs/filelistnorm.txt'
    pivot_wavelengths = '/home/mldantas/Dropbox/MSc_paper/Programs/pivot_lambdas.txt'

    # Constants --------------------------------------------------------------------------------------------------------
    c = 2.99792458E18	# Light Speed in Angstrom/s
    r_v = 3.1		    # Milky Way

    # The program itself -----------------------------------------------------------------------------------------------
    # Creating a dictionary --------------------------------------------------------------------------------------------

    my_data = np.loadtxt(data_obs, delimiter=',', dtype=object)

    dictionary = {}
    for i in range(len(my_data[0, :])):                             # Converting numpy array into dictionary
        dictionary[my_data[0, i]] = np.array(my_data[0 + 1:, i], dtype=str)

    object_id       = dictionary['Objid'].astype(long)
    plate           = dictionary['plate'].astype(int)
    mjd             = dictionary['mjd'].astype(int)
    fiberid         = dictionary['fiberid'].astype(int)
    fiber_mag_z     = dictionary['fiberMag_z'].astype(float)
    model_mag_u     = dictionary['modelMag_u'].astype(float)
    model_mag_g     = dictionary['modelMag_g'].astype(float)
    model_mag_r     = dictionary['modelMag_r'].astype(float)
    model_mag_i     = dictionary['modelMag_i'].astype(float)
    model_mag_z     = dictionary['modelMag_z'].astype(float)
    mag_fuv         = dictionary['fuv_mag'].astype(float)
    mag_nuv         = dictionary['nuv_mag'].astype(float)
    extinction_u    = dictionary['extinction_u'].astype(float)
    extinction_g    = dictionary['extinction_g'].astype(float)
    extinction_r    = dictionary['extinction_r'].astype(float)
    extinction_i    = dictionary['extinction_i'].astype(float)
    extinction_z    = dictionary['extinction_z'].astype(float)
    mag_dered_sdss  = my_data[1:, 14:19].astype(float)
    e_bv            = dictionary['e_bv'].astype(float)
    redshift        = dictionary['redshift'].astype(float)
    model_mag_u_err = dictionary['modelMagErr_u'].astype(float)
    model_mag_g_err = dictionary['modelMagErr_g'].astype(float)
    model_mag_r_err = dictionary['modelMagErr_r'].astype(float)
    model_mag_i_err = dictionary['modelMagErr_i'].astype(float)
    model_mag_z_err = dictionary['modelMagErr_z'].astype(float)
    mag_fuv_err     = dictionary['fuv_magerr'].astype(float)
    mag_nuv_err     = dictionary['nuv_magerr'].astype(float)


    # Loading filelist -------------------------------------------------------------------------------------------------
    filelist = np.loadtxt(filelist, dtype=str)
    indexes  = np.arange(plate.size)

    # SDSS Magnitude Offset correction (AB Magnitudes) ------------------------------------------------------------------
    mag_ab_u = model_mag_u - extinction_u - 0.04
    mag_ab_g = model_mag_g - extinction_g + 0.01
    mag_ab_r = model_mag_r - extinction_r + 0.01
    mag_ab_i = model_mag_i - extinction_i + 0.01
    mag_ab_z = model_mag_z - extinction_z + 0.02

    mag_ab_sdss        = np.array([mag_ab_u, mag_ab_g, mag_ab_r, mag_ab_i, mag_ab_z]).T
    mag_ab_sdss_griz   = np.array([mag_ab_g, mag_ab_r, mag_ab_i, mag_ab_z]).T
    model_mag_err_sdss = np.array([model_mag_u_err, model_mag_g_err, model_mag_r_err, model_mag_i_err,
                                   model_mag_z_err]).T


    # Matching observational data with synthetic data ------------------------------------------------------------------
    for spec_file in filelist:
        basename = os.path.split(spec_file)[-1]     # Get filename
        basename = os.path.splitext(basename)[0]    # Remove file extension

        plates    = basename.split('.')[0]
        mjds      = basename.split('.')[1]
        fiberids  = basename.split('.')[2]

        plates    = int(plates)
        mjds      = int(mjds)
        fiberids  = int(fiberids)

        index = indexes[(plate == plates) * (mjd == mjds) * (fiberid == fiberids)]
        if index.size is 0:
            continue

        # Observed data reduction and adjustment -----------------------------------------------------------------------
        # Loading effective wavelengths --------------------------------------------------------------------------------
        wavelengths_all        = np.loadtxt(pivot_wavelengths, usecols=[0, 1, 2, 3, 4, 5, 6], dtype=float)
        wavelengths_sdss       = np.loadtxt(pivot_wavelengths, usecols=[0, 1, 2, 3, 4], dtype=float)
        wavelengths_sdss_griz  = np.loadtxt(pivot_wavelengths, usecols=[1, 2, 3, 4], dtype=float)
        wavelengths_galex      = np.loadtxt(pivot_wavelengths, usecols=[5, 6], dtype=float)
        fuv_obs_wavelength     = np.loadtxt(pivot_wavelengths, usecols=[5], dtype=float)
        nuv_obs_wavelength     = np.loadtxt(pivot_wavelengths, usecols=[6], dtype=float)

        # Calculating Milky Way's extinction parameters using Cardelli's Extinction Law --------------------------------
        a_v = ['A_V']
        a_v = e_bv[index] * r_v
        q = cardelli_redlaw(wavelengths_all, r_v)
        a_lambda = np.zeros(wavelengths_all.size)
        for j in range(wavelengths_all.size):
            a_lambda[j] = q[j] * a_v

        # Calculating all SDSS fluxes ----------------------------------------------------------------------------------
        ## The following fluxes already account for extinction and K correction ----------------------------------------
        flux_sdss = ((c/(wavelengths_sdss ** 2.)) * 10.**(0.4 * (fiber_mag_z[index]-model_mag_u[index])) *
                     10. ** (-0.4 * (mag_ab_sdss[index][0] + 48.60))/1E-17)
                                                                                # AB mag into flux + 1e-17 normalization
        # Calculating g,r,i,z SDSS fluxes ------------------------------------------------------------------------------
        ## The following fluxes already account for extinction and K correction ----------------------------------------
        flux_sdss_griz = ((c/(wavelengths_sdss_griz ** 2.)) * 10.**(0.4 * (fiber_mag_z[index]-model_mag_u[index])) *
                     10. ** (-0.4 * (mag_ab_sdss_griz[index][0] + 48.60))/1E-17)
                                                                                # AB mag into flux + 1e-17 normalization

        # Calculating GALEX fluxes and correcting it for numerous effects ----------------------------------------------
        ## NUV band correction -----------------------------------------------------------------------------------------
        flux_nuv_fitzpatrick = ['flux_nuv_fitzpatrick']
        flux_nuv = (10. ** (-0.4 * (mag_nuv[index] - 20.08))) * 2.06 * 10 ** (-16)                   # AB mag into flux
        flux_nuv_fitzpatrick = fm_unred(nuv_obs_wavelength, flux_nuv, e_bv[index])    # NUV Flux corrected by extinction
        flux_nuv_fitzpatrick = flux_nuv_fitzpatrick / (1E-17)                                      # 1e-17 normalization

        ## FUV band correction -----------------------------------------------------------------------------------------
        flux_fuv_fitzpatrick = ['flux_fuv_fitzpatrick']
        flux_fuv = (10. ** (-0.4 * (mag_fuv[index] - 18.82))) * 1.40 * 10 ** (-15)                   # AB mag into flux
        flux_fuv_fitzpatrick = fm_unred(fuv_obs_wavelength, flux_fuv, e_bv[index])    # FUV Flux corrected by extinction
        flux_fuv_fitzpatrick = flux_fuv_fitzpatrick / (1E-17)                                      # 1e-17 normalization

        flux_galex = np.array([flux_fuv_fitzpatrick, flux_nuv_fitzpatrick])

        # Reading the synthetic spectra files --------------------------------------------------------------------------
        synthetic_wavelengths = np.loadtxt(os.path.join(spec_files_path, spec_file), usecols=[0])
        synthetic_fluxes      = np.loadtxt(os.path.join(spec_files_path, spec_file), usecols=[1])

        # Setting NUV and FUV mags dered using the original flux-magnitude transformation ------------------------------
        mag_nuv_dered_fm = -2.5 * np.log10(flux_nuv_fitzpatrick * (1E-17) / (2.06 * (10**(-16)))) + 20.08
        mag_fuv_dered_fm = -2.5 * np.log10(flux_fuv_fitzpatrick * (1E-17) / (1.40 * (10**(-15)))) + 18.82

        # Flux correction using chi^2 minimization ---------------------------------------------------------------------
        numerator_sum   = 0
        denominator_sum = 0
        for p in range(wavelengths_sdss_griz.size):
            index2 = np.abs(wavelengths_sdss_griz[p] - synthetic_wavelengths).argmin()
            numerator_sum   = numerator_sum + (flux_sdss_griz[p] * synthetic_fluxes[index2])
            denominator_sum = denominator_sum + (flux_sdss_griz[p] ** 2.)
        correction_factor = (numerator_sum / denominator_sum)

        flux_sdss_corr  = flux_sdss * correction_factor
        flux_nuv_corr   = flux_nuv_fitzpatrick * correction_factor
        flux_fuv_corr   = flux_fuv_fitzpatrick * correction_factor
        flux_galex_corr = np.array([float(flux_fuv_corr), float(flux_nuv_corr)])
        flux_all_corr   = np.array([float(flux_sdss_corr[0]), float(flux_sdss_corr[1]), float(flux_sdss_corr[2]),
                                  float(flux_sdss_corr[3]), float(flux_sdss_corr[4]), float(flux_fuv_corr),
                                  float(flux_nuv_corr)])

        # Chi^2 selection  ---------------------------------------------------------------------------------------------
        flux_sdss_synth = []
        for n in range(wavelengths_sdss_griz.size):
            index_flux_synth = np.abs(wavelengths_sdss_griz[n] - synthetic_wavelengths).argmin()
            flux_synth = synthetic_fluxes[index_flux_synth]
            flux_sdss_synth.append(flux_synth)
        flux_sdss_synth = np.array(flux_sdss_synth)

        goodness = ['goodness']
        goodness = s.chisquare(flux_sdss_corr[1:], flux_sdss_synth)
        absolute_goodness = goodness[0]/3


        # UV Synthetic equivalent fluxes -------------------------------------------------------------------------------
        k = np.abs(wavelengths_galex[0] - synthetic_wavelengths).argmin()
        fuv_synth_flux = synthetic_fluxes[k]

        m = np.abs(wavelengths_galex[1] - synthetic_wavelengths).argmin()
        nuv_synth_flux = synthetic_fluxes[m]

        # Absolute Magnitudes using Luminosity Distance ----------------------------------------------------------------
        luminosity_distance = cosmo.cosmocalc(redshift[index])['DL_Mpc']
        mag_fuv_abs         = mag_fuv_dered_fm - 5 * np.log10(luminosity_distance) - 25
        mag_nuv_abs         = mag_nuv_dered_fm - 5 * np.log10(luminosity_distance) - 25
        mag_u_abs           = mag_ab_u[index] - 5 * np.log10(luminosity_distance) - 25
        mag_g_abs           = mag_ab_g[index] - 5 * np.log10(luminosity_distance) - 25
        mag_r_abs           = mag_ab_r[index] - 5 * np.log10(luminosity_distance) - 25
        mag_i_abs           = mag_ab_i[index] - 5 * np.log10(luminosity_distance) - 25
        mag_z_abs           = mag_ab_z[index] - 5 * np.log10(luminosity_distance) - 25

        # SDSS Fluxes Uncertainty --------------------------------------------------------------------------------------
        model_mag_err_sdss_index = model_mag_err_sdss[index][0]
        sigma_flux=[]
        for i in range(mag_dered_sdss[index].size):
            sigma_magnitude_square = model_mag_err_sdss_index[i] ** 2
            sigma_flux_i = \
                np.sqrt(((((c/(wavelengths_sdss[i]**2))) * (0.4**2) * (mag_ab_sdss[i][0] + 48.60) *
                          10 ** (-0.4 * (mag_ab_sdss[i][0] + 48.60) - 1.)) ** 2) * sigma_magnitude_square) \
                * wavelengths_sdss[i]
            sigma_flux.append(float(sigma_flux_i))
        sigma_flux = np.array(sigma_flux)

        # GALEX Fluxes Uncertainty -------------------------------------------------------------------------------------
        ## FUV ---------------------------------------------------------------------------------------------------------
        sigma_fuv_flux = ['error_flux_fuv']
        sigma_fuv_flux = ((1.4 * 16 * (mag_fuv_dered_fm - 18.82) *
                          (10 ** (-0.4 * (mag_fuv_dered_fm - 18.82) - 1.))) * mag_fuv_err[index]) * fuv_obs_wavelength


        ## NUV ---------------------------------------------------------------------------------------------------------
        sigma_nuv_flux = ['error_flux_nuv']
        sigma_nuv_flux = ((20.6 * (0.4**2) * (mag_nuv_dered_fm - 20.08) *
                          (10 ** (-0.4 * (mag_nuv_dered_fm - 20.08) - 1.))) * mag_nuv_err[index]) * nuv_obs_wavelength


        sigma_uv = np.array([float(sigma_fuv_flux), float(sigma_nuv_flux)])


        # Plotting the results -----------------------------------------------------------------------------------------
        plt.semilogy(synthetic_wavelengths, synthetic_fluxes * synthetic_wavelengths, '-', linewidth=1., color='black',
                     alpha=1., label = r" Starlight Spectrum Object {0:s}".format(basename))
        plt.semilogy(wavelengths_sdss, flux_sdss_corr * wavelengths_sdss, 'o', markersize=15., color="#FF4500",
                     alpha=0.5, label = 'SDSS fluxes')
        plt.errorbar(wavelengths_sdss, flux_sdss_corr * wavelengths_sdss, yerr=sigma_flux, fmt='.', color='#FF4500',
                     elinewidth=3)
        plt.semilogy(wavelengths_galex, flux_galex_corr * wavelengths_galex, 'o', markersize=15., color="#0000FF",
                     alpha=0.5, label = 'GALEX fluxes')
        plt.errorbar(wavelengths_galex, flux_galex_corr * wavelengths_galex, yerr=sigma_uv, fmt='.', color='#0000FF',
                     elinewidth=3)
        plt.xlabel('$\lambda$ ($\AA$)', fontsize=25)
        plt.ylabel(r"$F_{\lambda} \lambda$ ($10^{-17} erg \, s^{-1} cm^{-2}$)", fontsize=25)
        plt.legend(loc='best', numpoints=1, fontsize=20, frameon=False)
        plt.grid(alpha=0.40)
        plt.tick_params('both', labelsize='20')
        plt.minorticks_on()
        plt.xlim([1000, 10000])
        figure = plt.gcf() # get current figure
        figure.set_size_inches(12, 8)
        # plt.arrow(wavelengths_galex[0], flux_fuv_corr[0] * wavelengths_galex[0], 0,
        #           fuv_synth_flux * wavelengths_galex[0]**(9/10), linewidth=1.5,head_width=100, head_length=(3 * 10**4),
        #           color='red', fc='red', ec='red', alpha=0.6)
        # plt.arrow(wavelengths_galex[1], flux_nuv_corr[0] * wavelengths_galex[1], 0,
        #           nuv_synth_flux * wavelengths_galex[1]**(8.2/10), linewidth=1.5,head_width=100, head_length=(1 * 10**4),
        #           color='red', fc='red', ec='red', alpha=0.6)
        # plt.tick_params('both', labelsize='28')
        plt.savefig(os.path.join(results_path, basename+'_final.png'), dpi = 100)
        # plt.show()
        plt.clf()
        # exit()

    #parameters_out.close()


__author__ = 'mldantas'
