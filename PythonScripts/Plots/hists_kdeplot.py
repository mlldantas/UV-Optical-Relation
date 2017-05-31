#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
    This code aims on running GMM on different UV-optical colours and check galaxy distributions.
    @author:  Maria Luiza Linhares Dantas
    @date:    2017.04.11
    @version: 0.0.1
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Main thread
if __name__ == '__main__':

    # Configuring the inputs -------------------------------------------------------------------------------------------
    my_data = pd.read_csv('/home/mldantas/Dropbox/MSc_paper/Programs/Results/output_results.txt', delim_whitespace=True)

    optical_colour_gr = my_data['mag_ab_g']-my_data['mag_ab_r']
    optical_colour_ur = my_data['mag_ab_u']-my_data['mag_ab_r']
    uv_colour_fuv_nuv = my_data['fuv_mag']-my_data['nuv_mag']
    optical_uv_colour = my_data['nuv_mag']-my_data['mag_ab_r']


    # Plotting the results ---------------------------------------------------------------------------------------------
    with sns.axes_style('whitegrid', {'axes.grid': True}):

        plot01 = plt.subplot(2, 2, 1)
        bins = np.arange(optical_colour_gr.min(), optical_colour_gr.max(), 0.03)
        plt.hist(optical_colour_gr, bins=bins, normed=True)
        sns.kdeplot(optical_colour_gr, shade=True, color='r')
        plt.xlabel(r"g-r", fontsize=15)
        plt.ylabel(r"Frequency", fontsize=15)
        plt.tick_params('both', labelsize='15')

        plt.subplot(2, 2, 2)
        bins = np.arange(optical_colour_ur.min(), optical_colour_ur.max(), 0.03)
        plt.hist(optical_colour_ur, bins=bins, normed=True)
        sns.kdeplot(optical_colour_ur, shade=True, color='r')
        plt.xlabel(r"u-r", fontsize=15)
        plt.ylabel(r"Frequency", fontsize=15)
        plt.xlim(optical_colour_ur.min(), 4)
        plt.tick_params('both', labelsize='15')

        plt.subplot(2, 2, 3)
        bins = np.arange(optical_uv_colour.min(), optical_uv_colour.max(), 0.03)
        plt.hist(optical_uv_colour, bins=bins, normed=True)
        sns.kdeplot(optical_uv_colour, shade=True, color='r')
        plt.xlabel(r"NUV-r", fontsize=15)
        plt.ylabel(r"Frequency", fontsize=15)
        plt.tick_params('both', labelsize='15')

        plt.subplot(2, 2, 4)
        bins = np.arange(uv_colour_fuv_nuv.min(), uv_colour_fuv_nuv.max(), 0.03)
        plt.hist(uv_colour_fuv_nuv, bins=bins, normed=True)
        sns.kdeplot(uv_colour_fuv_nuv, shade=True, color='r')
        plt.xlabel(r"FUV-NUV", fontsize=15)
        plt.ylabel(r"Frequency", fontsize=15)
        plt.tick_params('both', labelsize='15')
        
    plt.show()


    # temp_gr = np.array(optical_colour_gr).reshape(1, -1)
    # kernel_density = KernelDensity(kernel='gaussian', bandwidth=0.75).fit(temp_gr)
    # gr_density = kernel_density.score_samples(temp_gr)
    # x_axis = np.linspace(optical_colour_gr.min(), optical_colour_gr.max(), step)
    # print x_axis.shape
    # x_axis = np.array(x_axis).reshape(step,)
    # y_axis = np.exp(gr_density)
    # print x_axis.shape
    # print y_axis
    # print y_axis.shape
    # print optical_colour_gr.shape
    # plt.plot(x_axis, gr_density, '-')
    # plt.fill(x_axis[:, 0], np.exp(gr_density), fc='black')

    #
    # # Plot the progression of histograms to kernels
    # np.random.seed(1)
    # N = 20
    # X = np.concatenate((np.random.normal(0, 1, 0.3 * N),
    #                     np.random.normal(5, 1, 0.7 * N)))[:, np.newaxis]
    # print X.shape
    # print X
    # exit()
    #
    # X_plot = np.linspace(-5, 10, 1000)[:, np.newaxis]
    # bins = np.linspace(-5, 10, 10)
    #
    # fig, ax = plt.subplots(1, 1, sharex=True, sharey=True)
    #
    #
    # # Gaussian KDE
    # kde = KernelDensity(kernel='gaussian', bandwidth=0.75).fit(X)
    # log_dens = kde.score_samples(X_plot)
    # ax.fill(X_plot[:, 0], np.exp(log_dens), fc='#AAAAFF')
    #
    # print X_plot[:, 0].shape
    # print np.exp(log_dens).shape



    plt.show()

