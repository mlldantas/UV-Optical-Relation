#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
    Data Reduction and Adjustments

    @author:  Maria Luiza Linhares Dantas
    @date:    2015.09.18
    @version: 0.0.2


"""
# ======================================================================================================================
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import scipy.stats as s

# Main thread ----------------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    data_path = '/home/mldantas/Documentos/Programas/Synthesis_Results/Output'
    data_av   = '/home/mldantas/Documentos/Programas/Sintese'
    results   = '/home/mldantas/Documentos/Programas/Synthesis_Results/'
    filelist  = '/home/mldantas/Documentos/Programas/Synthesis_Results/filelist.txt'
    f1517_base= '/home/mldantas/Documentos/Programas/Synthesis_Results/base150.txt'

    cardelli1517 = 2.642

    filelist = np.loadtxt(filelist, dtype=str)
    f1517    = np.loadtxt(f1517_base, usecols=[3])
    base     = np.loadtxt(f1517_base, usecols=[0], dtype=str)
    age      = np.loadtxt(f1517_base, usecols=[1], dtype=float)

    m_number = []
    base_number = []
    for i in range(base.size):
        m_number_k = str.split(base[i], '_')[2]
        m_number_k = int(m_number_k[1:])
        basename_base = str.split(base[i], '_')[-1]                      # Getting final information about the base
        basename_base = str.split(basename_base, '.')[0]                 # Getting base number
        m_number.append(int(m_number_k))
        base_number.append(int(basename_base))
    base_number = np.array(base_number)                                  # Array of all templates' numbers
    m_number = np.array(m_number)

    av_out = open('/home/mldantas/Documentos/Programas/Synthesis_Results/avs.txt', 'w')
    print >> av_out, '#plate', '%s' % 'mjd', '%s' % 'fiberid','%s' % 'av_out'
    plates  =[]
    mjds    =[]
    fiberids=[]
    av = []
    for infile in glob.glob(os.path.join(data_av, '*')):
        basename_av = os.path.split(infile)[-1]
        basename_av = os.path.splitext(basename_av)[0]
        basename_av = os.path.splitext(basename_av)[0]    # Remove file extension
        basename_av = os.path.splitext(basename_av)[0]    # Remove file extension
        basename_av = os.path.splitext(basename_av)[0]    # Remove file extension
        basename_av = os.path.splitext(basename_av)[0]    # Remove file extension
        platesi, mjdsi, fiberidsi = basename_av.split('.')
        plates.append(platesi)
        mjds.append(mjdsi)
        fiberids.append(fiberidsi)
        review_file = open(infile, 'r').read(3011)
        avi = review_file[3003:3011]
        av.append(float(avi))
        print >> av_out, avi, platesi, mjdsi, fiberidsi
    av       = np.array(av)
    plates   = np.array(plates)
    mjds     = np.array(mjds)
    fiberids = np.array(fiberids)


    for infile in filelist:
        file = np.loadtxt(os.path.join(data_path, infile), dtype=str)
        basename = os.path.split(infile)[-1]
        basename = os.path.splitext(basename)[0]    # Remove file extension
        plate, mjd, fiberid = basename.split('.')
        av_i = av[np.where((plate==plates) * (mjd==mjds) * (fiberid==fiberids))]
        print av_i, plate, mjd, fiberid
        parameters_out = open('/home/mldantas/Documentos/Programas/Synthesis_Results/'+basename+'_age.txt', 'w')
        print >> parameters_out, '#component_j', '%9s' % 'x_j(%)', '%7s' % 'y_i',  '%13s' % 'Age'

        x      = file[0:, 1].astype(float)
        base_j = file[0:, 9].astype(str)
        m_number_ii =[]
        basename_results = []
        for i in range(base_j.size):
            basename_results_i = str.split(base_j[i], '_')[0][3:]
            m_number_i = str.split(base_j[i], '_')[1]
            m_number_i = int(m_number_i[1:])
            basename_results.append(int(basename_results_i))
            m_number_ii.append(m_number_i)
        basename_results = np.array(basename_results)
        m_number_ii = np.array(m_number_ii)
        m_number_j = m_number_ii
        # index = np.where(x!=0)
        # xi = x[index]
        componentj = basename_results
        base_j = base_j

        indexes  = np.arange(base_number.size)
        indexx = indexes[(componentj == base_number) * (m_number_j == m_number)]
        if indexx is 0:
            continue
        weights = []
        fluxes = []
        components = []
        for indexx in range(base_number.size):
            mi = m_number_j[indexx]
            ci = componentj[indexx]
            wi = x[indexx]
            f1517_i = f1517[indexx]
            fi = wi * f1517_i *(10 ** (-0.4 * av_i * cardelli1517))
            if (fi == 0):
                continue
            print >> parameters_out, '%5s' % base_j[indexx], '%10.4f' % wi, '%10.4f' % fi, '%15.2e' % age[indexx]
            weights.append(wi)
            fluxes.append(fi)
            components.append(ci)
        weights    = np.array(weights)
        fluxes     = np.array(fluxes)
        components = np.array(components)

        sum = 0.
        for i in range(weights.size):
            sum = sum + (fluxes[i] * weights[i])
        print >> parameters_out, '\n#sum = %.4f' % sum
        print >> parameters_out, '#av = %.4f' % av_i
        i = i+1

        parameters_out.close()


__author__ = 'mldantas'
