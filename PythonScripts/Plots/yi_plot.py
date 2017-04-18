#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
    This program plots the Yi et al. (2011) diagram.
    @author:  Maria Luiza Linhares Dantas
    @date:    2017.17.04
    @version: 0.0.1
"""

import matplotlib.pyplot as plt
import numpy as np

# ======================================================================================================================

# Main thread
if __name__ == '__main__':

    my_data = np.loadtxt('/home/mldantas/Dropbox/MSc_paper/Programs/Results/output_results_with_lines_d4000.csv',
                         delimiter=',', dtype=str)
    my_dictionary = {}
    for i in range(len(my_data[0, :])):  # Converting numpy array into dictionary
        my_dictionary[my_data[0, i]] = np.array(my_data[0 + 1:, i], dtype=str)

    nii = my_dictionary['F_nii'].astype(float)
    halpha = my_dictionary['F_Halpha'].astype(float)
    xbpt = (nii / halpha)*100

    # uv_upturn_index = np.where((((my_dictionary['nuv_mag'].astype(float)) -
    #                             (my_dictionary['dered_r'].astype(float)))> 5.4))
    # print uv_upturn_index

    plt.scatter(((my_dictionary['nuv_mag'].astype(float)) - (my_dictionary['dered_r'].astype(float))),
                ((my_dictionary['fuv_mag'].astype(float)) - (my_dictionary['nuv_mag'].astype(float))), s=xbpt, alpha=0.3)
    # plt.scatter(((my_dictionary['nuv_mag'].astype(float)[uv_upturn_index]) - (my_dictionary['dered_r'].astype(float)[uv_upturn_index])),
    #             ((my_dictionary['fuv_mag'].astype(float)[uv_upturn_index]) - (my_dictionary['nuv_mag'].astype(float)[uv_upturn_index])), s=xbpt, c='red',
    #             alpha=0.7)
    plt.axvline(x=5.4, color='black', linewidth=2.)
    plt.axhline(y=0.9, color='black', linewidth=2.)
    plt.xlabel("NUV-r")
    plt.ylabel("FUV-NUV")
    plt.show()
