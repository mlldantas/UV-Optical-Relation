#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
    Sort program
    @author:  Maria Luiza L. Dantas
    @date:    2015.30.11
    @version: 0.0.1

"""

import numpy as np

# ======================================================================================================================
# Main thread
if __name__ == '__main__':

    data = 'tabelaexemplo.txt'

    parametro1 = np.loadtxt(data, usecols=[0]) #e.g. identifiers
    parametro2 = np.loadtxt(data, usecols=[1])
    parametro3 = np.loadtxt(data, usecols=[2])
    table      = np.loadtxt(data)

    # Criando lista de ids não repetidas -------------------------------------------------------------------------------

    parametro1 = set(list(parametro1))

    print table

    new_list = []
    for i in parametro1:
        mask = (table[:, 0] == i)
        print i,mask
        #subtable = table[mask,:]
        #par_max = np.argmax(subtable[:, 2])
        par_max = np.argmax(parametro3*mask)
        #print subtable
        print par_max
        #new_list.append(subtable[par_max,:])
        new_list.append(table[par_max, :])
    print new_list
    new_table = np.array(new_list)

    print new_table


__author__ = 'mldantas'
