#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
    Data Reduction and Adjustments
    @author:  Maria Luiza Linhares Dantas
    @date:    2016.21.06
    @version: 0.1.10
    This program uses synthetic spectra produced by STARLIGHT and adjusts observational data using chi^2 minimization
    method.
    In this version we completed the loop in order to analyze all the data and implement a single output.
    For more information about the input data, please read the schema text file.
    Main modification in this version: fiber correction implemented correctly as well as errors accounting for these
    mesurements.
"""

from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import os
from fitzpatrick import fm_unred
from cardelli import cardelli_redlaw
import cosmocalc as cosmo
import scipy.stats as s
import csv

# ======================================================================================================================

# Main thread
if __name__ == '__main__':

    # Configuring the inputs -------------------------------------------------------------------------------------------
    results_path      = '/home/mldantas/Dropbox/MSc_paper/Programs/Results'
    figs_path         = '/home/mldantas/Dropbox/MSc_paper/Programs/Results/Fig'
    spec_files_path   = '/home/mldantas/Dropbox/MSc_paper/Dataset/Rebinned_Spec_SelfZ'
    my_data           = '/home/mldantas/Dropbox/MSc_paper/Dataset/MyData_GALEX_SDSS_clean_sort.csv'
    filelist          = '/home/mldantas/Dropbox/MSc_paper/Dataset/filelist_individual_spectra.txt'
    pivot_wavelengths = '/home/mldantas/Dropbox/MSc_paper/Programs/pivot_lambdas.txt'
    fuv_response      = '/home/mldantas/Dropbox/MSc_paper/UV_filter_response/fuv_filter_response_1angstrom.txt'
    nuv_response      = '/home/mldantas/Dropbox/MSc_paper/UV_filter_response/nuv_filter_response_1angstrom.txt'

    # Constants --------------------------------------------------------------------------------------------------------
    c = 2.99792458E18	# Light Speed in Angstrom/s
    r_v = 3.1		    # Milky Way

    # The program itself -----------------------------------------------------------------------------------------------
    # Creating a dictionary --------------------------------------------------------------------------------------------

    my_data = np.loadtxt(my_data, delimiter=',', dtype=object)

    dictionary = {}
    for i in range(len(my_data[0, :])):                                         # Converting numpy array into dictionary
        dictionary[my_data[0, i]] = np.array(my_data[0 + 1:, i], dtype=str)

    objid         = dictionary['Objid'].astype(long)
    plate         = dictionary['plate'].astype(int)
    mjd           = dictionary['mjd'].astype(int)
    fiberid       = dictionary['fiberid'].astype(int)
    ra            = dictionary['ra'].astype(float)
    dec           = dictionary['dec'].astype(float)
    fibermag_z    = dictionary['fiberMag_z'].astype(float)
    modelmag_u    = dictionary['modelMag_u'].astype(float)
    modelmag_g    = dictionary['modelMag_g'].astype(float)
    modelmag_r    = dictionary['modelMag_r'].astype(float)
    modelmag_i    = dictionary['modelMag_i'].astype(float)
    modelmag_z    = dictionary['modelMag_z'].astype(float)
    extinction_u  = dictionary['extinction_u'].astype(float)
    extinction_g  = dictionary['extinction_g'].astype(float)
    extinction_r  = dictionary['extinction_r'].astype(float)
    extinction_i  = dictionary['extinction_i'].astype(float)
    extinction_z  = dictionary['extinction_z'].astype(float)
    fuv_mag       = dictionary['fuv_mag'].astype(float)
    nuv_mag       = dictionary['nuv_mag'].astype(float)
    dered_u       = dictionary['dered_u'].astype(float)
    dered_g       = dictionary['dered_g'].astype(float)
    dered_r       = dictionary['dered_r'].astype(float)
    dered_i       = dictionary['dered_i'].astype(float)
    dered_z       = dictionary['dered_z'].astype(float)
    petro_r90_r   = dictionary['petroR90_r'].astype(float)
    mod_mag_err_u = dictionary['modelMagErr_u'].astype(float)
    mod_mag_err_g = dictionary['modelMagErr_g'].astype(float)
    mod_mag_err_r = dictionary['modelMagErr_r'].astype(float)
    mod_mag_err_i = dictionary['modelMagErr_i'].astype(float)
    mod_mag_err_z = dictionary['modelMagErr_z'].astype(float)
    fuv_magerr    = dictionary['fuv_magerr'].astype(float)
    nuv_magerr    = dictionary['nuv_magerr'].astype(float)
    e_bv          = dictionary['e_bv'].astype(float)
    redshift      = dictionary['redshift'].astype(float)
    s2n_r         = dictionary['s2n_r'].astype(float)
    sn_fuv        = dictionary['sn_fuv_auto'].astype(float)
    survey        = dictionary['survey'].astype(str)
    morph_type    = dictionary['morph_type'].astype(float)

    fuv_filter_wavelength = np.loadtxt(fuv_response, usecols=[0])
    fuv_filter_response   = np.loadtxt(fuv_response, usecols=[1])

    nuv_filter_wavelength = np.loadtxt(nuv_response, usecols=[0])
    nuv_filter_response   = np.loadtxt(nuv_response, usecols=[1])

    # Creating output file> > parameters_out ---------------------------------------------------------------------------
    output_file = open('/home/mldantas/Dropbox/MSc_paper/Programs/Results/output_results.csv', 'w')

    input_file  = open('/home/mldantas/Dropbox/MSc_paper/Dataset/MyData_GALEX_SDSS_clean_sort.csv')
    reader_file = csv.reader(input_file)

    lines   = my_data[1:, :].astype(object)
    columns = my_data[1:, :].astype(object)

    lines_no_objid = my_data[1:, 1:].astype(object)

    print >> output_file, '#', '%10s' % 'Objid', '%14s' % 'plate', '%6s' % 'mjd', '%10s' % 'fiberid', '%8s' % 'ra', \
        '%13s' % 'dec', '%16s' % 'fiberMag_z', '%12s' % 'modelMag_u', '%12s' % 'modelMag_g', \
        '%12s' % 'modelMag_r', '%12s' % 'modelMag_i', '%12s' % 'modelMag_z', '%14s' % 'extinction_u', \
        '%14s' % 'extinction_g', '%14s' % 'extinction_r', '%14s' % 'extinction_i', '%14s' % 'extinction_z', \
        '%10s' % 'fuv_mag', '%12s' % 'nuv_mag', '%12s' % 'dered_u', '%12s' % 'dered_g',  '%12s' % 'dered_r', \
        '%12s' % 'dered_i', '%12s' % 'dered_z', '%13s' %'petroR90_r', '%14s' % 'modelMagErr_u', '%14s' % 'modelMagErr_g', \
        '%14s' % 'modelMagErr_r', '%14s' % 'modelMagErr_i', '%14s' % 'modelMagErr_z', '%12s' % 'fuv_magerr', \
        '%12s' % 'nuv_magerr', '%9s' % 'e_bv', '%13s' % 'redshift', '%9s' % 's2n_r', '%15s' % 'sn_fuv_auto', \
        '%8s' % 'survey', '%12s' % 'morph_type', '%10s' % 'mag_ab_u', '%10s' % 'mag_ab_g', '%10s' % 'mag_ab_r', \
        '%10s' % 'mag_ab_i', '%10s' % 'mag_ab_z', '%8s' % 'av_mw', '%13s' % 'a_lambda_u', '%11s' % 'a_lambda_g', \
        '%11s' % 'a_lambda_r', '%11s' % 'a_lambda_i', '%11s' % 'a_lambda_z', '%13s' % 'a_lambda_fuv', \
        '%13s' % 'a_lambda_nuv', '%11s' % 'fiber_corr', '%13s' % 'flux_u(E-17)', '%13s' % 'flux_g(E-17)', \
        '%13s' % 'flux_r(E-17)', '%13s' % 'flux_i(E-17)', '%13s' % 'flux_z(E-17)', '%15s' % 'flux_fuv(E-17)', \
        '%15s' % 'flux_nuv(E-17)', '%17s' % 'flux_u_err(E-17)', '%17s' % 'flux_g_err(E-17)', \
        '%17s' % 'flux_r_err(E-17)', '%17s' % 'flux_i_err(E-17)', '%17s' % 'flux_z_err(E-17)', \
        '%19s' % 'flux_fuv_err(E-17)', '%19s' % 'flux_nuv_err(E-17)', '%17s' % 'flux_u_esc(E-17)', \
        '%17s' % 'flux_g_esc(E-17)', '%17s' % 'flux_r_esc(E-17)', '%17s' % 'flux_i_esc(E-17)', \
        '%17s' % 'flux_z_esc(E-17)', '%19s' % 'flux_fuv_esc(E-17)', '%19s' % 'flux_nuv_esc(E-17)', \
        '%15s' % 'scaling_factor', '%16s' % 'chi^2(goodness)', '%18s' % 'absolute_goodness', '%10s' % 'Mag_Abs_u', \
        '%11s' % 'Mag_Abs_g', '%11s' % 'Mag_Abs_r', '%11s' % 'Mag_Abs_i', '%11s' % 'Mag_Abs_z'

    # Loading filelist -------------------------------------------------------------------------------------------------
    filelist = np.loadtxt(filelist, dtype=str)
    indexes  = np.arange(plate.size)

    # SDSS Magnitude Offset correction (AB Magnitudes) -----------------------------------------------------------------
    mag_ab_u = modelmag_u - extinction_u - 0.04
    mag_ab_g = modelmag_g - extinction_g + 0.01
    mag_ab_r = modelmag_r - extinction_r + 0.01
    mag_ab_i = modelmag_i - extinction_i + 0.01
    mag_ab_z = modelmag_z - extinction_z + 0.02

    mag_ab_sdss        = np.array([mag_ab_u, mag_ab_g, mag_ab_r, mag_ab_i, mag_ab_z]).T
    mag_ab_sdss_griz   = np.array([mag_ab_g, mag_ab_r, mag_ab_i, mag_ab_z]).T
    model_mag_err_sdss = np.array([mod_mag_err_u, mod_mag_err_g, mod_mag_err_r, mod_mag_err_i, mod_mag_err_z]).T
    sdss_mag_dered     = np.array([dered_u, dered_g, dered_r, dered_i, dered_z]).T
    mod_mag_all        = np.array([modelmag_u, modelmag_g, modelmag_r, modelmag_i, modelmag_z]).T

    # Matching observational data with synthetic data ------------------------------------------------------------------
    for spec_file in filelist:
        basename = os.path.split(spec_file)[-1]     # Get filename
        basename = os.path.splitext(basename)[0]    # Remove file extension

        plates   = basename.split('.')[0]
        mjds     = basename.split('.')[1]
        fiberids = basename.split('.')[2]

        plates   = int(plates)
        mjds     = int(mjds)
        fiberids = int(fiberids)

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

        # Calculating the local fiber correction -----------------------------------------------------------------------
        if ((fibermag_z[index]-modelmag_z[index]) < 1):
            fiber_correction = 1
        else:
            fiber_correction = 10**(0.4*(fibermag_z[index]-modelmag_z[index]))

        # Calculating all SDSS fluxes ----------------------------------------------------------------------------------
        ## The following fluxes already account for extinction and K correction ----------------------------------------
        flux_sdss = fiber_correction * \
                    ((c/(wavelengths_sdss ** 2.)) * 10. ** (-0.4 * (mag_ab_sdss[index][0] + 48.60))/1E-17)
                                                                                # AB mag into flux + 1e-17 normalization
        # Calculating g,r,i,z SDSS fluxes ------------------------------------------------------------------------------
        ## The following fluxes already account for extinction and K correction ----------------------------------------
        flux_sdss_griz = fiber_correction * ((c/(wavelengths_sdss_griz ** 2.))
                                             * 10. ** (-0.4 * (mag_ab_sdss_griz[index][0] + 48.60))/1E-17)
                                                                                # AB mag into flux + 1e-17 normalization

        # Calculating GALEX fluxes and correcting it for numerous effects ----------------------------------------------
        ## NUV band correction -----------------------------------------------------------------------------------------
        flux_nuv = (10. ** (-0.4 * (nuv_mag[index] - 20.08))) * 2.06 * 10 ** (-16)                   # AB mag into flux
        flux_nuv_fitzpatrick = fm_unred(nuv_obs_wavelength, flux_nuv, e_bv[index])    # NUV Flux corrected by extinction
        flux_nuv_fitzpatrick = flux_nuv_fitzpatrick / (1E-17)                                      # 1e-17 normalization

        ## FUV band correction -----------------------------------------------------------------------------------------
        flux_fuv = (10. ** (-0.4 * (fuv_mag[index] - 18.82))) * 1.40 * 10 ** (-15)                   # AB mag into flux
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

        goodness = s.chisquare(flux_sdss_corr[1:], flux_sdss_synth)
        absolute_goodness = goodness[0]/3

        # UV Synthetic equivalent fluxes -------------------------------------------------------------------------------
        ## Only necessary for the example plot (with arrows)
        # k = np.abs(wavelengths_galex[0] - synthetic_wavelengths).argmin()
        # fuv_synth_flux = synthetic_fluxes[k]
        #
        # m = np.abs(wavelengths_galex[1] - synthetic_wavelengths).argmin()
        # nuv_synth_flux = synthetic_fluxes[m]

        # Absolute Magnitudes using Luminosity Distance ----------------------------------------------------------------
        # This is optional to calculate,so for now it will be kept commented -------------------------------------------
        luminosity_distance = cosmo.cosmocalc(redshift[index])['DL_Mpc']
        # mag_fuv_abs         = mag_fuv_dered_fm - 5 * np.log10(luminosity_distance) - 25
        # mag_nuv_abs         = mag_nuv_dered_fm - 5 * np.log10(luminosity_distance) - 25
        mag_u_abs           = mag_ab_u[index] - 5 * np.log10(luminosity_distance) - 25
        mag_g_abs           = mag_ab_g[index] - 5 * np.log10(luminosity_distance) - 25
        mag_r_abs           = mag_ab_r[index] - 5 * np.log10(luminosity_distance) - 25
        mag_i_abs           = mag_ab_i[index] - 5 * np.log10(luminosity_distance) - 25
        mag_z_abs           = mag_ab_z[index] - 5 * np.log10(luminosity_distance) - 25

        # SDSS Fluxes Uncertainty --------------------------------------------------------------------------------------
        model_mag_err_sdss_index = model_mag_err_sdss[index][0]
        sigma_flux= []
        for i in range(sdss_mag_dered[index].size):
            if (fiber_correction == 1):
                sigma_fiber_corr_square = 0.
            else:
                sigma_fiber_corr_square = (-0.4 * (modelmag_z[index] - fibermag_z[index])
                                           * 10 ** (0.4*(fibermag_z[index] - modelmag_z[index]) - 1)) ** 2 \
                                          + mod_mag_err_z[index] ** 2
            sigma_magnitude_square = model_mag_err_sdss_index[i] ** 2
            sigma_flux_i = np.sqrt(((((c/(wavelengths_sdss[i]**2))) * 2 * (0.4**2) * (mag_ab_sdss[i][0] + 48.60)
                                     * 10 ** (-0.4 * (mag_ab_sdss[i][0] + 48.60) - 1.)) ** 2)
                                   * sigma_magnitude_square + sigma_fiber_corr_square)
            sigma_flux.append(float(sigma_flux_i))
        sigma_flux = np.array(sigma_flux)


        # GALEX Fluxes Uncertainty -------------------------------------------------------------------------------------
        ## FUV ---------------------------------------------------------------------------------------------------------
        sigma_fuv_flux = ((1.4 * (0.4**2) * 10**(-15) * (mag_fuv_dered_fm - 18.82)
                           * (10 ** (-0.4 * (mag_fuv_dered_fm - 18.82) - 1.))) * fuv_magerr[index])/(1E-17)


        ## NUV ---------------------------------------------------------------------------------------------------------
        sigma_nuv_flux = ((2.06 * (0.4**2) *  10**(-16) *(mag_nuv_dered_fm - 20.08)
                           * (10 ** (-0.4 * (mag_nuv_dered_fm - 20.08) - 1.))) * nuv_magerr[index])/(1E-17)

        ## UV sigmas
        sigma_uv = np.array([float(sigma_fuv_flux), float(sigma_nuv_flux)])


        # Saving the results in a new file -----------------------------------------------------------------------------
        print >> output_file, '%18d' % objid[index], '%7d' % plate[index], '%8d' % mjd[index], '%7d' % fiberid[index],\
            '%13.4f' % ra[index], '%13.7f' % dec[index], '%12.5f' % fibermag_z[index], '%12.5f' % modelmag_u[index], \
            '%12.5f' % modelmag_g[index], '%12.5f' % modelmag_r[index], '%12.5f' % modelmag_i[index], \
            '%12.5f' % modelmag_z[index], '%12.5f' % extinction_u[index], '%14.5f' % extinction_g[index], \
            '%14.5f' % extinction_r[index], '%14.5f' % extinction_i[index], '%14.5f' % extinction_z[index], \
            '%14.5f' % fuv_mag[index], '%12.5f' % nuv_mag[index], '%12.5f' % dered_u[index], '%12.5f' % dered_g[index], \
            '%12.5f' % dered_r[index], '%12.5f' % dered_i[index], '%12.5f' % dered_z[index], \
            '%10.5f' % petro_r90_r[index], '%13.5f' % mod_mag_err_u[index], '%14.5f' % mod_mag_err_g[index], \
            '%14.5f' % mod_mag_err_r[index], '%14.5f' % mod_mag_err_i[index], '%14.5f' % mod_mag_err_z[index], \
            '%13.5f' % fuv_magerr[index], '%12.5f' % nuv_magerr[index], '%14.8f' % e_bv[index], \
            '%11.7f' % redshift[index], '%9.4f' % s2n_r[index], '%12.6f' % sn_fuv[index], '%11s' % survey[index], \
            '%7d' % morph_type[index], '%14.5f' % mag_ab_u[index], '%10.5f' % mag_ab_g[index], \
            '%10.5f' % mag_ab_r[index], '%10.5f' % mag_ab_i [index], '%10.5f' % mag_ab_z[index], '%10.5f' % a_v, \
            '%10.5f' % a_lambda[0], '%11.5f' % a_lambda[1], '%11.5f' % a_lambda[2], '%11.5f' % a_lambda[3], \
            '%11.5f' % a_lambda[4], '%11.5f' % a_lambda[5], '%13.5f' % a_lambda[6], '%13.5f' % fiber_correction, \
            '%12.5f' % flux_sdss[0], '%13.5f' % flux_sdss[1], '%13.5f' % flux_sdss[2], '%13.5f' % flux_sdss[3], \
            '%13.5f' % flux_sdss[4], '%13.5f' % flux_fuv_fitzpatrick, '%15.5f' % flux_nuv_fitzpatrick, \
            '%15.5f' % sigma_flux[0], '%17.5f' % sigma_flux[1], '%17.5f' % sigma_flux[2], '%17.5f' % sigma_flux[3], \
            '%17.5f' % sigma_flux[4], '%19.5f' % sigma_fuv_flux, '%19.5f' % sigma_nuv_flux, '%19.5f' % flux_sdss_corr[0], \
            '%17.5f' % flux_sdss_corr[1], '%17.5f' % flux_sdss_corr[2], '%17.5f' % flux_sdss_corr[3], \
            '%17.5f' % flux_sdss_corr[4], '%17.5f' % flux_fuv_corr, '%19.5f' % flux_nuv_corr, \
            '%17.5f' % correction_factor, '%16.5f' % goodness[0], '%16.5f' % absolute_goodness, '%16.5f' % mag_u_abs, \
            '%11.5f' % mag_g_abs, '%11.5f' % mag_r_abs, '%11.5f' % mag_i_abs, '%11.5f' % mag_z_abs

        # Plotting the results -----------------------------------------------------------------------------------------
        plt.semilogy(synthetic_wavelengths, synthetic_fluxes * synthetic_wavelengths, '-', linewidth=1., color='black',
                     alpha=1., label=r" Starlight Spectrum Object {0:s}".format(basename))
        plt.semilogy(wavelengths_sdss, flux_sdss_corr * wavelengths_sdss, 'o', markersize=15., color="#FF4500",
                     alpha=0.5, label = 'SDSS fluxes')
        plt.errorbar(wavelengths_sdss, flux_sdss_corr * wavelengths_sdss, yerr=(sigma_flux * wavelengths_sdss), fmt='.',
                     color='#FF4500', elinewidth=3)
        plt.semilogy(wavelengths_galex, flux_galex_corr * wavelengths_galex, 'o', markersize=15., color="#0000FF",
                     alpha=0.5, label = 'GALEX fluxes')
        plt.errorbar(wavelengths_galex, flux_galex_corr * wavelengths_galex, yerr=(sigma_uv * wavelengths_galex),
                     fmt='.', color='#0000FF', elinewidth=3)
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
        plt.savefig(os.path.join(figs_path, basename+'.png'), dpi = 100)
        plt.clf()


    output_file.close()


__author__ = 'mldantas'
