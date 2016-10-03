#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
    Sort program
    @author:  Maria Luiza L. Dantas
    @date:    2016.10.03
    @version: 0.0.2
    This version we select the best NUV magnitude error from GALEX in the GAMA catalog.

"""

import numpy as np

# ======================================================================================================================
# Main thread ----------------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    data = np.loadtxt('/home/mldantas/Dropbox/DoutoradoIAG/GAMAZOO/GAMA_zoo_GALEX.csv', delimiter=',', dtype=object)

    dictionary = {}
    for i in range(len(data[0, :])):                             # Converting numpy array into dictionary
        dictionary[data[0, i]] = np.array(data[0 + 1:, i], dtype=object)

    objectid = dictionary['OBJID'].astype(str)
    nuv_err  = dictionary['BEST_FLUXERR_NUV'].astype(float)

    table = data[1:, :]

    # Criando lista de ids não repetidas -------------------------------------------------------------------------------

    ids = set(list(objectid))

    new_list = []
    for i in ids:
        mask = (table[:, 0] == i)
        nuv_err_min = np.argmin(nuv_err*mask)
        new_list.append(table[nuv_err_min, :])
    new_table = np.array(new_list)


    np.savetxt('/home/mldantas/Dropbox/DoutoradoIAG/GAMAZOO/GAMA_zoo_GALEX_clean.csv', new_table, fmt='%s',
               delimiter=',', newline='\n',
               header='CATAID,OBJID,RA,DEC,THETA_IMAGE,THETA_J2000,THRESHOLD,MU_THRESHOLD,FLUX_PETRO_U,FLUXERR_PETRO_U,'
                      'MAG_PETRO_U,MAGERR_PETRO_U,FLUX_PETRO_G,FLUXERR_PETRO_G,MAG_PETRO_G,MAGERR_PETRO_G,FLUX_PETRO_R,'
                      'FLUXERR_PETRO_R,MAG_PETRO_R,MAGERR_PETRO_R,FLUX_PETRO_I,FLUXERR_PETRO_I,MAG_PETRO_I,'
                      'MAGERR_PETRO_I,FLUX_PETRO_Z,FLUXERR_PETRO_Z,MAG_PETRO_Z,MAGERR_PETRO_Z,FLUX_PETRO_Y,'
                      'FLUXERR_PETRO_Y,MAG_PETRO_Y,MAGERR_PETRO_Y,FLUX_PETRO_J,FLUXERR_PETRO_J,MAG_PETRO_J,'
                      'MAGERR_PETRO_J,FLUX_PETRO_H,FLUXERR_PETRO_H,MAG_PETRO_H,MAGERR_PETRO_H,FLUX_PETRO_K,'
                      'FLUXERR_PETRO_K,MAG_PETRO_K,MAGERR_PETRO_K,Column1,Column2,Column3,Column4,EFF_EXPTIME_NUV,'
                      'EFF_EXPTIME_FUV,BEST_FLUX_NUV,BEST_FLUXERR_NUV,BEST_MAG_NUV,BEST_MAGERR_NUV,BEST_FLUX_FUV,'
                      'BEST_FLUXERR_FUV,BEST_MAG_FUV,BEST_MAGERR_FUV,BEST_METHOD,NMATCHUV,NMATCHOPT,NUVFLAG,FUVFLAG,'
                      'TOT_FLUX_NUV,TOT_FLUXERR_NUV,TOT_MAG_NUV,TOT_MAGERR_NUV,TOT_FLUX_FUV,TOT_FLUXERR_FUV,TOT_MAG_FUV,'
                      'TOT_MAGERR_FUV,SPLIT_FLUX_NUV,SPLIT_FLUXERR_NUV,SPLIT_MAG_NUV,SPLIT_MAGERR_NUV,SPLIT_FLUX_FUV,'
                      'SPLIT_FLUXERR_FUV,SPLIT_MAG_FUV,SPLIT_MAGERR_FUV,COG_FLUX_NUV,COG_FLUXERR_NUV,COG_MAG_NUV,'
                      'COG_MAGERR_NUV,COG_FLUX_FUV,COG_FLUXERR_FUV,COG_MAG_FUV,COG_MAGERR_FUV,COG_RE_CRIT_NUV,'
                      'COG_HP_RAD_NUV,COG_RE_CRIT_FUV,COG_HP_RAD_FUV,COG_CONFSOURCE_NUV,COG_CONFSOURCE_FUV,NN_GGOID,'
                      'NN_RA_GALEX,NN_DEC_GALEX,NN_DIST,NN_NMATCH4,NN_MANY2ONE,NN_FLUX_NUV,NN_FLUXERR_NUV,NN_MAG_NUV,'
                      'NN_MAGERR_NUV,NN_FLUX_FUV,NN_FLUXERR_FUV,NN_MAG_FUV,NN_MAGERR_FUV,NN_ARTIFACT_NUV,'
                      'NN_ARTIFACT_FUV,NN_SFLAGS_NUV,NN_SFLAGS_FUV,NN_FLUX50_RADIUS_NUV,NN_FLUX50_RADIUS_FUV,'
                      'NN_SEMIMAJOR_NUV,NN_SEMIMINOR_NUV,NN_SEMIMAJORERR_NUV,NN_SEMIMINORERR_NUV,NN_POSANG_NUV,'
                      'NN_POSANGERR_NUV')

__author__ = 'mldantas'
