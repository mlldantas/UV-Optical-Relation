#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
    Data Reduction and Adjustments

    @author:  Maria Luiza Linhares Dantas
    @date:    2015.01.08
    @version: 0.1.8

    This program uses synthetic spectra produced by STARLIGHT and adjusts observational data using chi^2 minimization
    method.

    In this version we completed the loop in order to analyze all the data and implement a single output.

    For more information about the input data, please read the schema text file.

    Main modification in this version: we are using SDSS' dereddenned magnitudes and Fitzpatrick's Law to correct the
    UV magnitudes.
    
    This version accounts the fiber correction in the flux calculation.

"""

from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import os
from fitzpatrick import fm_unred
from cardelli import cardelli_redlaw
import cosmocalc as cosmo
import scipy.stats as s

# ======================================================================================================================

# Main thread
if __name__ == '__main__':

    # Configuring the inputs -------------------------------------------------------------------------------------------
    spec_files_path   = '/home/mldantas/Documentos/Programas/Integracao/EspectrosRebinadosSelfZ'
    results_path      = '/home/mldantas/Documentos/Programas/Integracao/Results/IntegracaoV10'
    data_obs_whan     = '/home/mldantas/Documentos/Programas/Integracao/Results/Car_Fitz/data_match_002.txt'
    filelist          = '/home/mldantas/Documentos/Programas/Integracao/filelistnorm.txt'
    pivot_wavelengths = '/home/mldantas/Documentos/Programas/Integracao/pivot_lambdas.txt'
    fmagz             = '/home/mldantas/Documentos/Programas/Integracao/fibermagz.csv'

    # Constants --------------------------------------------------------------------------------------------------------
    c = 2.99792458E18	# Light Speed in Angstrom/s
    r_v = 3.1		    # Milky Way

    # The program itself -----------------------------------------------------------------------------------------------
    # Loading observational data ---------------------------------------------------------------------------------------
    objid                = np.loadtxt(data_obs_whan, usecols=[0], dtype = long)
    plate                = np.loadtxt(data_obs_whan, usecols=[1])
    mjd                  = np.loadtxt(data_obs_whan, usecols=[2])
    fiberid              = np.loadtxt(data_obs_whan, usecols=[3])
    alpha                = np.loadtxt(data_obs_whan, usecols=[4])
    delta                = np.loadtxt(data_obs_whan, usecols=[5])
    mag_obs_all          = np.loadtxt(data_obs_whan, usecols=[6, 7, 8, 9, 10, 11, 12])
    mag_obs_sdss         = np.loadtxt(data_obs_whan, usecols=[6, 7, 8, 9, 10])
    mag_obs_sdss_griz    = np.loadtxt(data_obs_whan, usecols=[7, 8, 9, 10])
    mag_obs_galex        = np.loadtxt(data_obs_whan, usecols=[11, 12])
    u_band               = np.loadtxt(data_obs_whan, usecols=[6])
    g_band               = np.loadtxt(data_obs_whan, usecols=[7])
    r_band               = np.loadtxt(data_obs_whan, usecols=[8])
    i_band               = np.loadtxt(data_obs_whan, usecols=[9])
    z_band               = np.loadtxt(data_obs_whan, usecols=[10])
    fuv_band             = np.loadtxt(data_obs_whan, usecols=[11])
    nuv_band             = np.loadtxt(data_obs_whan, usecols=[12])
    mag_dered_sdss       = np.loadtxt(data_obs_whan, usecols=[13, 14, 15, 16, 17])
    mag_dered_u          = np.loadtxt(data_obs_whan, usecols=[13])
    mag_dered_g          = np.loadtxt(data_obs_whan, usecols=[14])
    mag_dered_r          = np.loadtxt(data_obs_whan, usecols=[15])
    mag_dered_i          = np.loadtxt(data_obs_whan, usecols=[16])
    mag_dered_z          = np.loadtxt(data_obs_whan, usecols=[17])
    mag_ab_sdss_corr_str = np.loadtxt(data_obs_whan, usecols=[18, 19, 20, 21, 22])
    mag_ab_sdss_corr_griz= np.loadtxt(data_obs_whan, usecols=[19, 20, 21, 22])
    mag_ab_u_corr_str    = np.loadtxt(data_obs_whan, usecols=[18]) #Cardelli extinction, offset and k corrections made
    mag_ab_g_corr_str    = np.loadtxt(data_obs_whan, usecols=[19]) #Cardelli extinction, offset and k corrections made
    mag_ab_r_corr_str    = np.loadtxt(data_obs_whan, usecols=[20]) #Cardelli extinction, offset and k corrections made
    mag_ab_i_corr_str    = np.loadtxt(data_obs_whan, usecols=[21]) #Cardelli extinction, offset and k corrections made
    mag_ab_z_corr_str    = np.loadtxt(data_obs_whan, usecols=[22]) #Cardelli extinction, offset and k corrections made
    mag_ab_sdss_corr_m   = np.loadtxt(data_obs_whan, usecols=[23, 24, 25, 26, 27])
    mag_ab_u_corr_m      = np.loadtxt(data_obs_whan, usecols=[23]) #Cardelli extinction, Offset Corrected
    mag_ab_g_corr_m      = np.loadtxt(data_obs_whan, usecols=[24]) #Cardelli extinction, Offset Corrected
    mag_ab_r_corr_m      = np.loadtxt(data_obs_whan, usecols=[25]) #Cardelli extinction, Offset Corrected
    mag_ab_i_corr_m      = np.loadtxt(data_obs_whan, usecols=[26]) #Cardell/ (1.1**3)i extinction, Offset Corrected
    mag_ab_z_corr_m      = np.loadtxt(data_obs_whan, usecols=[27]) #Cardelli extinction, Offset Corrected
    mmag_sdss_err        = np.loadtxt(data_obs_whan, usecols=[28, 29, 30, 31, 32])
    mmag_u_err           = np.loadtxt(data_obs_whan, usecols=[28]) #Model Mag u measurement error
    mmag_g_err           = np.loadtxt(data_obs_whan, usecols=[29]) #Model Mag u measurement error
    mmag_r_err           = np.loadtxt(data_obs_whan, usecols=[30]) #Model Mag u measurement error
    mmag_i_err           = np.loadtxt(data_obs_whan, usecols=[31]) #Model Mag u measurement error
    mmag_z_err           = np.loadtxt(data_obs_whan, usecols=[32]) #Model Mag u measurement error
    mag_fuv_err          = np.loadtxt(data_obs_whan, usecols=[33]) #Model Mag u measurement error
    mag_nuv_err          = np.loadtxt(data_obs_whan, usecols=[34])
    extinction_u         = np.loadtxt(data_obs_whan, usecols=[35])
    extinction_g         = np.loadtxt(data_obs_whan, usecols=[36])
    extinction_r         = np.loadtxt(data_obs_whan, usecols=[37])
    extinction_i         = np.loadtxt(data_obs_whan, usecols=[38])
    extinction_z         = np.loadtxt(data_obs_whan, usecols=[39])
    petro90_r            = np.loadtxt(data_obs_whan, usecols=[40])
    e_bv                 = np.loadtxt(data_obs_whan, usecols=[41])
    s2n_r                = np.loadtxt(data_obs_whan, usecols=[42])
    s2n_fuv              = np.loadtxt(data_obs_whan, usecols=[43])
    z                    = np.loadtxt(data_obs_whan, usecols=[44])
    z_err                = np.loadtxt(data_obs_whan, usecols=[45])
    av_starlight         = np.loadtxt(data_obs_whan, usecols=[46])
    dn4000_synth         = np.loadtxt(data_obs_whan, usecols=[47])
    dn4000_obs           = np.loadtxt(data_obs_whan, usecols=[48])
    morphological_type   = np.loadtxt(data_obs_whan, usecols=[49])
    survey               = np.loadtxt(data_obs_whan, usecols=[50], dtype = str)
    x_whan               = np.loadtxt(data_obs_whan, usecols=[51])
    y_whan               = np.loadtxt(data_obs_whan, usecols=[52])
    class_whan           = np.loadtxt(data_obs_whan, usecols=[53], dtype = int)
    young_pop_ratio      = np.loadtxt(data_obs_whan, usecols=[54])
    interm_pop_ratio     = np.loadtxt(data_obs_whan, usecols=[55])
    old_pop_ratio        = np.loadtxt(data_obs_whan, usecols=[56])

    fibermag_z           = np.loadtxt(fmagz, delimiter=',', usecols=[5])

    # Loading filelist -------------------------------------------------------------------------------------------------
    filelist = np.loadtxt(filelist, dtype=str)
    indexes  = np.arange(plate.size)

    # Creating output file> > parameters_out ---------------------------------------------------------------------------
    parameters_out \
        = open('/home/mldantas/Documentos/Programas/Integracao/Results/IntegracaoV10/outputdata.txt', 'w')

    print >> parameters_out, '%s' % '#', '%11s' % 'ObjectID', '%13s' % 'plate', '%6s' % 'mjd', '%10s' % 'fiberID', \
        '%8s' % 'RA', '%12s' % 'Dec', '%17s' % 'model_mag_u', '%15s' % 'model_mag_g', '%14s' % 'model_mag_r', \
        '%14s' % 'model_mag_i', '%14s' % 'model_mag_z', '%12s' % 'mag_fuv', '%12s' % 'mag_nuv', \
        '%17s' % 'mmag_extinc_u', '%17s' % 'mmag_extinc_g', '%17s' % 'mmag_extinc_r', '%17s' % 'mmag_extinc_i',\
        '%17s' % 'mmag_extinc_z', '%15s' % 'mmag_err_u', '%15s' % 'mmag_err_g','%15s' % 'mmag_err_r',\
        '%15s' % 'mmag_err_i','%15s' % 'mmag_err_z', '%16s' % 'mag_err_fuv', '%16s' % 'mag_err_nuv', \
        '%15s' % 'mag_u_dered', '%14s' % 'mag_g_dered', '%14s' % 'mag_r_dered', '%14s' % 'mag_i_dered', \
        '%14s' % 'mag_z_dered', '%20s' % 'mag_fuv_dered_fm', '%20s' % 'mag_nuv_dered_fm',\
        '%15s' % 'mag_ab_u', '%15s' % 'mag_ab_g', '%15s' % 'mag_ab_r', '%15s' % 'mag_ab_i', '%15s' % 'mag_ab_z',\
        '%16s' % 'mag_u_abs', '%18s' % 'mag_g_abs', '%18s' % 'mag_r_abs', '%18s' % 'mag_i_abs', '%18s' % 'mag_z_abs', \
        '%18s' % 'mag_fuv_abs', '%18s' % 'mag_nuv_abs', '%14s' % 'flux_u', '%16s' % 'flux_g', \
        '%16s' % 'flux_r', '%16s' % 'flux_i', '%16s' % 'flux_z', '%14s' % 'flux_fuv', '%14s' % 'flux_nuv',\
        '%18s' % 'flux_sca_u', '%18s' % 'flux_sca_g', '%18s' % 'flux_sca_r', '%18s' % 'flux_sca_i',\
        '%18s' % 'flux_sca_z', '%20s' % 'flux_sca_fuv', '%20s' % 'flux_sca_nuv', '%20s' % 'flux_fuv_synth', \
        '%17s' % 'flux_nuv_synth', '%9s' % 'R90_r', '%13s' % 'E(B-V)', '%12s' % 'S2N_r', '%17s' % 'S2N_FUV_auto', \
        '%12s' % 'redshift', '%17s' % 'redshift_err', '%10s' % 'AV_MK', '%16s' % 'AV_GAL_synth', '%12s' % 'A_lamb_u', \
        '%12s' % 'A_lamb_g', '%12s' % 'A_lamb_r', '%12s' % 'A_lamb_i', '%12s' % 'A_lamb_z', '%12s' % 'A_lamb_fuv', \
        '%12s' % 'A_lamb_nuv', '%16s' % 'Dn4000_Synth', '%12s' % 'Dn4000_Obs', '%14s' % 'Morph_Type', '%9s' % 'Survey',\
        '%20s' % 'Scaling_Factor', '%17s' % 'log10(NII/Ha)', '%15s' % 'log10(EW(Ha))', '%13s' % 'WHAN_Class', \
        '%16s' % 'young_pop_ratio', '%18s' % 'interm_pop_ratio', '%15s' % 'old_pop_ratio', \
        '%26s' % 'Luminosity_Distance(Mpc)', '%17s' % 'Chi2_Goodness', '%27s' % 'Chi2_Relative_Goodness'


    # SDSS Magnitude Offset correction ---------------------------------------------------------------------------------
    mag_ab_u = u_band - extinction_u - 0.04
    mag_ab_g = g_band - extinction_g + 0.01
    mag_ab_r = r_band - extinction_r + 0.01
    mag_ab_i = i_band - extinction_i + 0.01
    mag_ab_z = z_band - extinction_z + 0.02

    mag_ab_sdss      = np.array([mag_ab_u, mag_ab_g, mag_ab_r, mag_ab_i, mag_ab_z]).T
    mag_ab_sdss_griz = np.array([mag_ab_g, mag_ab_r, mag_ab_i, mag_ab_z]).T


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
        a_v = e_bv[index] * r_v
        q = cardelli_redlaw(wavelengths_all, r_v)
        a_lambda = np.zeros(wavelengths_all.size)
        for j in range(wavelengths_all.size):
            a_lambda[j] = q[j] * a_v

        # Calculating all SDSS fluxes ----------------------------------------------------------------------------------
        ## The following fluxes already account for extinction and K correction ----------------------------------------
        flux_sdss = ((c/(wavelengths_sdss ** 2.)) * 10.**(0.4 * (fibermag_z[index]-u_band[index])) *
                     10. ** (-0.4 * (mag_ab_sdss[index][0] + 48.60))/1E-17)
                                                                                # AB mag into flux + 1e-17 normalization
        # Calculating g,r,i,z SDSS fluxes ------------------------------------------------------------------------------
        ## The following fluxes already account for extinction and K correction ----------------------------------------
        flux_sdss_griz = \
            ((c/(wavelengths_sdss_griz ** 2.)) * 10.**(0.4 * (fibermag_z[index]-u_band[index])) *
             10. ** (-0.4 * (mag_ab_sdss_griz[index][0] + 48.60))/1E-17)
                                                                                # AB mag into flux + 1e-17 normalization

        # Calculating GALEX fluxes and correcting it for numerous effects ----------------------------------------------
        ## NUV band correction -----------------------------------------------------------------------------------------
        flux_nuv = (10. ** (-0.4 * (nuv_band[index] - 20.08))) * 2.06 * 10 ** (-16)                   # AB mag into flux
        flux_nuv_fitzpatrick = fm_unred(nuv_obs_wavelength, flux_nuv, e_bv[index])    # NUV Flux corrected by extinction
        flux_nuv_fitzpatrick = flux_nuv_fitzpatrick / (1E-17)                                      # 1e-17 normalization

        ## FUV band correction -----------------------------------------------------------------------------------------
        flux_fuv = (10. ** (-0.4 * (fuv_band[index] - 18.82))) * 1.40 * 10 ** (-15)                   # AB mag into flux
        flux_fuv_fitzpatrick = fm_unred(fuv_obs_wavelength, flux_fuv, e_bv[index])    # FUV Flux corrected by extinction
        flux_fuv_fitzpatrick = flux_fuv_fitzpatrick / (1E-17)                                      # 1e-17 normalization

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
            numerator_sum   = numerator_sum + (flux_sdss[p] * synthetic_fluxes[index2])
            denominator_sum = denominator_sum + (flux_sdss[p] ** 2.)
        correction_factor = (numerator_sum / denominator_sum)

        flux_sdss_corr = flux_sdss * correction_factor
        flux_nuv_corr  = flux_nuv_fitzpatrick * correction_factor
        flux_fuv_corr  = flux_fuv_fitzpatrick * correction_factor
        flux_galex_corr = np.array([float(flux_fuv_corr), float(flux_nuv_corr)])
        flux_all_corr = np.array([float(flux_sdss_corr[0]), float(flux_sdss_corr[1]), float(flux_sdss_corr[2]),
                                  float(flux_sdss_corr[3]), float(flux_sdss_corr[4]), float(flux_fuv_corr),
                                  float(flux_nuv_corr)])

        # UV Synthetic equivalent fluxes -------------------------------------------------------------------------------
        k = np.abs(fuv_obs_wavelength - synthetic_wavelengths).argmin()
        fuv_synth_flux = synthetic_fluxes[k]

        m = np.abs(nuv_obs_wavelength - synthetic_wavelengths).argmin()
        nuv_synth_flux = synthetic_fluxes[m]

        # Absolute Magnitudes using Luminosity Distance ----------------------------------------------------------------
        luminosity_distance = cosmo.cosmocalc(z[index])['DL_Mpc']
        mag_fuv_abs         = mag_fuv_dered_fm - 5 * np.log10(luminosity_distance) - 25
        mag_nuv_abs         = mag_nuv_dered_fm - 5 * np.log10(luminosity_distance) - 25
        mag_u_abs           = mag_ab_u[index] - 5 * np.log10(luminosity_distance) - 25
        mag_g_abs           = mag_ab_g[index] - 5 * np.log10(luminosity_distance) - 25
        mag_r_abs           = mag_ab_r[index] - 5 * np.log10(luminosity_distance) - 25
        mag_i_abs           = mag_ab_i[index] - 5 * np.log10(luminosity_distance) - 25
        mag_z_abs           = mag_ab_z[index] - 5 * np.log10(luminosity_distance) - 25

        # Chi^2 selection  ---------------------------------------------------------------------------------------------
        flux_sdss_synth = []
        for n in range(wavelengths_sdss_griz.size):
            index_flux_synth = np.abs(wavelengths_sdss_griz[n] - synthetic_wavelengths).argmin()
            flux_synth = synthetic_fluxes[index_flux_synth]
            flux_sdss_synth.append(flux_synth)
        flux_sdss_synth = np.array(flux_sdss_synth)

        goodness = s.chisquare(flux_sdss_corr[1:], flux_sdss_synth)
        absolute_goodness = goodness[0]/3

        # SDSS Fluxes Uncertainty --------------------------------------------------------------------------------------
        mmag_sdss_err_index = mmag_sdss_err[index][0]
        sigma_flux = []
        for i in range(mag_dered_sdss[index].size):
            sigma_magnitude_square = mmag_sdss_err_index[i] ** 2
            sigma_flux_i = \
                np.sqrt(((((c/(wavelengths_sdss[i]**2)) * 1E17) * (0.4**2) * (mag_ab_sdss[i][0] + 48.60) *
                          10 ** (-0.4 * (mag_ab_sdss[i][0] + 48.60) - 1.)) ** 2) * sigma_magnitude_square) \
                * wavelengths_sdss[i]
            sigma_flux.append(float(sigma_flux_i))
        sigma_flux = np.array(sigma_flux)

        # GALEX Fluxes Uncertainty -------------------------------------------------------------------------------------
        ## FUV ---------------------------------------------------------------------------------------------------------

        sigma_fuv_flux = ((1.4 * 16 * (mag_fuv_dered_fm - 18.82) *
                          (10 ** (-0.4 * (mag_fuv_dered_fm - 18.82) - 1.))) * mag_fuv_err[index]) \
                        * fuv_obs_wavelength


        ## NUV ---------------------------------------------------------------------------------------------------------
        sigma_nuv_flux = ((20.6 * (0.4**2) * (mag_nuv_dered_fm - 20.08) *
                          (10 ** (-0.4 * (mag_nuv_dered_fm - 20.08) - 1.))) * mag_nuv_err[index]) \
                         * nuv_obs_wavelength


        sigma_uv = np.array([float(sigma_fuv_flux), float(sigma_nuv_flux)])


        #Saving the results in a new file -----------------------------------------------------------------------------
        print >> parameters_out, '%18d' % objid[index], '%7d' % plate[index], '%8d' % mjd[index], '%7d' % fiberid[index],\
            '%13.4f' % alpha[index], '%13.7f' % delta[index], '%12.5f' % u_band[index], '%14.5f' % g_band[index], \
            '%14.5f' % r_band[index], '%14.5f' % i_band[index], '%14.5f' % z_band[index], '%14.5f' % fuv_band[index], \
            '%12.5f' % nuv_band[index], '%14.5f' % extinction_u[index], '%17.5f' % extinction_g[index], \
            '%17.5f' % extinction_r[index], '%17.5f' % extinction_i[index], '%17.5f' % extinction_z[index], \
            '%16.5f' % mmag_u_err[index], '%15.5f' % mmag_g_err[index], '%15.5f' % mmag_r_err[index], \
            '%15.5f' % mmag_i_err[index], '%15.5f' % mmag_z_err[index], '%16.5f' % mag_fuv_err[index], \
            '%16.5f' % mag_nuv_err[index], '%15.5f' % mag_dered_u[index], '%14.5f' % mag_dered_g[index], \
            '%14.5f' % mag_dered_r[index], '%14.5f' % mag_dered_i[index], '%14.5f' % mag_dered_z[index], \
            '%17.5f' % mag_fuv_dered_fm, '%20.5f' % mag_nuv_dered_fm, '%20.5f' % mag_ab_u[index], \
            '%15.5f' % mag_ab_g[index], '%15.5f' % mag_ab_r[index], '%15.5f' % mag_ab_i[index], \
            '%15.5f' % mag_ab_z[index], '%15.5f' % mag_u_abs, '%18.5f' % mag_g_abs, '%18.5f' % mag_r_abs, \
            '%18.5f' % mag_i_abs, '%18.5f' % mag_z_abs, '%18.5f' % mag_fuv_abs, '%18.5f' % mag_nuv_abs, \
            '%16.5f' % flux_sdss[0], '%16.5f' % flux_sdss[1], '%16.5f' % flux_sdss[2], '%16.5f' % flux_sdss[3], \
            '%16.5f' % flux_sdss[4], '%13.5f' % flux_fuv_fitzpatrick, '%14.5f' % flux_nuv_fitzpatrick, \
            '%17.5f' % flux_sdss_corr[0], '%18.5f' % flux_sdss_corr[1], '%18.5f' % flux_sdss_corr[2], \
            '%18.5f' % flux_sdss_corr[3], '%18.5f' % flux_sdss_corr[4], '%20.5f' % flux_fuv_corr, \
            '%20.5f' % flux_nuv_corr, '%18.5f' % fuv_synth_flux, '%16.5f' % nuv_synth_flux, '%14.5f' % petro90_r[index],\
            '%14.8f' % e_bv[index], '%11.4f' % s2n_r[index], '%14.6f' % s2n_fuv[index], '%15.7f' % z[index], \
            '%14.2e' % z_err[index], '%12.4f' % a_v, '%13.4f' % av_starlight[index], '%13.3f' % a_lambda[0], \
            '%12.3f' % a_lambda[1],  '%12.3f' % a_lambda[2], '%12.3f' % a_lambda[3], '%12.3f' % a_lambda[4],\
            '%12.3f' % a_lambda[5], '%12.3f' % a_lambda[6], '%14.3f' % dn4000_synth[index], \
            '%13.3f' % dn4000_obs[index], '%12d' % morphological_type[index], '%15s' % survey[index], \
            '%15.5f' % correction_factor, '%18.5f' % x_whan[index], '%15.5f' % y_whan[index],\
            '%12d' % class_whan[index], '%16.5f' % young_pop_ratio[index], '%18.5f' % interm_pop_ratio[index], \
            '%15.5f' % old_pop_ratio[index], '%17.2f' % luminosity_distance, '%27.5f' % float(goodness[0]), \
            '%25.5f' % absolute_goodness


        # Plotting the results -----------------------------------------------------------------------------------------
        plt.semilogy(synthetic_wavelengths, synthetic_fluxes * synthetic_wavelengths, '-', linewidth=1., color='black',
                  label = r" Starlight Spectrum Object {0:s}".format(basename))
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
        #plt.show()
        plt.clf()
        #exit()


    parameters_out.close()


__author__ = 'mldantas'
