#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 21:33:30 2017

@author: filipmichalsky
"""
import sys
sys.path.append('../')

import sqlite3
import numpy as np

class Thermo():
    """ Class of thermodynanmics
    """

    def __init__ (self, species, T, ki, vi_p, vi_dp, db_name='NASA_coef.sqlite'):
        self.species = species
        self.T = T
        self.ki = ki
        self.b_ki = None
        self.vi_p = vi_p
        self.vi_dp = vi_dp
        self.dao = ThermoDAO(db_name)
        
    def get_backward_coefs (self):
        vi_p = np.matrix(self.vi_p).T
        vi_dp = np.matrix(self.vi_dp).T
        vi = vi_dp - vi_p  # calculate overall stoicheometric coefficients

        H_species = []
        S_species = []
        for s in self.species:
            if self.T >= 1000: # high temperature range
                NASA_coefs = self.dao.get_coeffs(s, 'high')
            else: # low temperature range
                NASA_coefs = self.dao.get_coeffs(s, 'low')

            # Compute enthalpy/(RT) of a specie
            H_T_arr = np.array([1, 1/2*self.T, 1/3*(self.T**2), 1/4*(self.T**3), 1/5*(self.T**4), 1/self.T])
            H_coef_arr = np.array([NASA_coefs[0], NASA_coefs[1], NASA_coefs[2], NASA_coefs[3], NASA_coefs[4], NASA_coefs[5]])
            H_species.append(np.dot(H_T_arr, H_coef_arr))

            # Compute entropy/R of a specie
            S_T_arr = np.array([np.log(self.T), self.T, 1/2*(self.T**2), 1/3*(self.T**3), 1/4*(self.T**4), 1])
            S_coef_arr = np.array([NASA_coefs[0], NASA_coefs[1], NASA_coefs[2], NASA_coefs[3], NASA_coefs[4], NASA_coefs[6]])
            S_species.append(np.dot(S_T_arr, S_coef_arr))

        delta_H = np.squeeze(np.asarray(np.dot(H_species, vi))) # delta enthalpy of a reaction
        delta_S = np.squeeze(np.asarray(np.dot(S_species, vi))) # delta entropy of a reaction

        gamma = np.sum(vi, axis=0)
        
        p0 = 10e5
        R = 8.314
        ke = np.squeeze(np.asarray(np.multiply(np.power(p0/(R*self.T), gamma), np.exp(delta_S-delta_H))))

        self.b_ki = self.ki/ke

        return self.b_ki


class ThermoDAO():
    """ Database Access Object for thermodynanmics
    """

    def __init__ (self, db_name='NASA_coef.sqlite'):
        self.db = sqlite3.connect(db_name)
        self.cursor = self.db.cursor()


    def get_coeffs (self, species_name, temp_range):
        query = '''SELECT COEFF_1, COEFF_2, COEFF_3, COEFF_4, COEFF_5, COEFF_6, COEFF_7
                    FROM {} 
                    WHERE SPECIES_NAME = "{}"'''.format(temp_range.upper(), species_name)
        coeffs = list(self.cursor.execute(query).fetchall()[0])
        return coeffs