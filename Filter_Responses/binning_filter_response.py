#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
    Data Reduction and Adjustments
    @author:  Maria Luiza Linhares Dantas
    @date:    2016.03.08
    @version: 0.0.1
    This program makes the new binning for one all the FUV and NUV filter responses.
"""
# ======================================================================================================================

from __future__ import division
import numpy as np
import scipy.interpolate as s
import os

# ======================================================================================================================

# Main thread ----------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    
    fuv_response = np.loadtxt('FUV_filter_response.txt')
    nuv_response = np.loadtxt('NUV_filter_response.txt')

    wavelength_fuv      = fuv_response[:, 0]
    filter_response_fuv = fuv_response[:, 1]
    
    wavelength_nuv      = nuv_response[:, 0]
    filter_response_nuv = nuv_response[:, 1]
    
    # For the FUV band:
    
    interp_function_fuv = s.interp1d(wavelength_fuv, filter_response_fuv)
    new_wavelength_fuv  = np.linspace(np.ceil(wavelength_fuv.min()), np.ceil(wavelength_fuv.max()), np.ceil(wavelength_fuv.max()-wavelength_fuv.min()+1))
    new_wavelength_fuv  = np.arange(np.ceil(wavelength_fuv.min()), np.ceil(wavelength_fuv.max()), 1)
    new_filter_response_fuv = np.arange(filter_response_fuv.min(), filter_response_fuv.max(), np.ceil(wavelength_fuv.max()-wavelength_fuv.min()+1))
    new_filter_response_fuv = interp_function_fuv(new_wavelength_fuv)

    np.savetxt('fuv_filter_response_1angstrom.txt', np.column_stack((new_wavelength_fuv, new_filter_response_fuv)),
               fmt='%s', delimiter='       ', header='wavelength   response', comments='#')
    
    # For the NUV band:
    interp_function_nuv = s.interp1d(wavelength_nuv, filter_response_nuv)
    new_wavelength_nuv  = np.linspace(np.ceil(wavelength_nuv.min()), np.ceil(wavelength_nuv.max()), np.ceil(wavelength_nuv.max()-wavelength_nuv.min()+1))
    new_wavelength_nuv  = np.arange(np.ceil(wavelength_nuv.min()), np.ceil(wavelength_nuv.max()), 1)
    new_filter_response_nuv = np.arange(filter_response_nuv.min(), filter_response_nuv.max(), np.ceil(wavelength_nuv.max()-wavelength_nuv.min()+1))
    new_filter_response_nuv = interp_function_nuv(new_wavelength_nuv)

    np.savetxt('nuv_filter_response_1angstrom.txt', np.column_stack((new_wavelength_nuv, new_filter_response_nuv)),
              fmt='%s', delimiter='       ', header='wavelength   response', comments='#')

__author__ = 'mldantas'
