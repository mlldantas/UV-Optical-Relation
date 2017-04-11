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
from scipy.stats import norm
from sklearn.neighbors.kde import KernelDensity

# Main thread
if __name__ == '__main__':

    # Configuring the inputs -------------------------------------------------------------------------------------------
    my_data = pd.read_csv('/home/mldantas/Dropbox/MSc_paper/Programs/Results/output_results.txt', delim_whitespace=True)

    # Using kernel density from scikit-learn to check for double gaussians ---------------------------------------------
    optical_colour_gr = my_data['dered_g']-my_data['dered_r']
    temp_gr = np.array(optical_colour_gr).reshape(1, -1)
    kernel_density = KernelDensity(kernel='gaussian', bandwidth=0.75).fit(temp_gr)
    gr_density = kernel_density.score_samples(temp_gr)

    # Plotting the results ---------------------------------------------------------------------------------------------
    bins = np.arange(optical_colour_gr.min(), optical_colour_gr.max(), 0.03)
    step = 10000
    x_axis = np.linspace(optical_colour_gr.min(), optical_colour_gr.max(), step)[:, np.newaxis]
    x_axis = np.array(x_axis).reshape(step,)
    y_axis = np.exp(gr_density)
    print x_axis.shape
    print y_axis
    print y_axis.shape
    plt.hist(optical_colour_gr, bins=bins)
    print optical_colour_gr.shape
    # sns.kdeplot(optical_colour_gr, shade=True, color='r')
    # plt.plot(x_axis, gr_density, '-')
    # plt.fill(x_axis[:, 0], np.exp(gr_density), fc='black')
    plt.xlabel(r"g-r")
    plt.ylabel(r"Frequency")
    plt.show()



    # Plot the progression of histograms to kernels
    np.random.seed(1)
    N = 20
    X = np.concatenate((np.random.normal(0, 1, 0.3 * N),
                        np.random.normal(5, 1, 0.7 * N)))[:, np.newaxis]
  
    X_plot = np.linspace(-5, 10, 1000)[:, np.newaxis]
    bins = np.linspace(-5, 10, 10)

    fig, ax = plt.subplots(1, 1, sharex=True, sharey=True)


    # Gaussian KDE
    kde = KernelDensity(kernel='gaussian', bandwidth=0.75).fit(X)
    log_dens = kde.score_samples(X_plot)
    ax.fill(X_plot[:, 0], np.exp(log_dens), fc='#AAAAFF')

    print X_plot[:, 0].shape
    print np.exp(log_dens).shape



    plt.show()

