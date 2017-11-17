import os.path
import sqlite3
import numpy as np
from chemkin.chemkin_errors import ChemKinError

class Thermo():
    """ Class of thermodynanmics
    """

    def __init__ (self, species, T, ki, vi_p, vi_dp, db_name='NASA_coef.sqlite'):
        self.p0 = 10e5
        self.R = 8.314
        self.species = species
        self.T = T
        self.ki = ki
        self.b_ki = None
        self.vi_p = vi_p
        self.vi_dp = vi_dp
        self.vi = np.matrix(self.vi_dp).T - np.matrix(self.vi_p).T  # calculate overall stoicheometric coefficients
        self.gamma = np.sum(self.vi, axis=0)
        self.dao = ThermoDAO(db_name)


    def get_backward_coefs (self):
        species_high = self.dao.get_species(self.T, 'high')
        species_low = self.dao.get_species(self.T, 'low')

        H_over_RT_species = []
        S_over_R_species = []
        for s in self.species:
            if self.T >= 1000: # high temperature range
                if s not in species_high:
                    raise ChemKinError('Thermo.get_backward_coefs()',
                        'The specie {}\'s high temperature bound is not higher than the current T={}.'.format(s, self.T))
                NASA_coefs = self.dao.get_coeffs(s, 'high')
            else: # low temperature range
                if s not in species_low:
                    raise ChemKinError('Thermo.get_backward_coefs()',
                        'The specie {}\'s low temperature bound is not lower than the current T={}.'.format(s, self.T))
                NASA_coefs = self.dao.get_coeffs(s, 'low')

            # Compute enthalpy/(RT) of a specie
            H_T_arr = np.array([1, 1/2*self.T, 1/3*(self.T**2), 1/4*(self.T**3), 1/5*(self.T**4), 1/self.T])
            H_coef_arr = np.array([NASA_coefs[0], NASA_coefs[1], NASA_coefs[2], NASA_coefs[3], NASA_coefs[4], NASA_coefs[5]])
            H_over_RT_species.append(np.dot(H_T_arr, H_coef_arr))

            # Compute entropy/R of a specie
            S_T_arr = np.array([np.log(self.T), self.T, 1/2*(self.T**2), 1/3*(self.T**3), 1/4*(self.T**4), 1])
            S_coef_arr = np.array([NASA_coefs[0], NASA_coefs[1], NASA_coefs[2], NASA_coefs[3], NASA_coefs[4], NASA_coefs[6]])
            S_over_R_species.append(np.dot(S_T_arr, S_coef_arr))

        delta_H_over_RT = np.squeeze(np.asarray(np.dot(H_over_RT_species, self.vi))) # delta enthalpy of a system of reactions
        delta_S_over_R = np.squeeze(np.asarray(np.dot(S_over_R_species, self.vi))) # delta entropy of a system of reactions
        
        ke = np.squeeze(np.asarray(np.multiply(np.power(self.p0/(self.R*self.T), self.gamma), np.exp(delta_S_over_R-delta_H_over_RT))))

        self.b_ki = self.ki/ke

        return self.b_ki

class ThermoDAO():
    """ Database Access Object for thermodynanmics
    """

    def __init__ (self, db_name):
        BASE_DIR = os.path.dirname(os.path.abspath(__file__))
        self.db_path = os.path.join(BASE_DIR, db_name)

    def get_coeffs (self, species_name, temp_range):
        db = sqlite3.connect(self.db_path)
        cursor = db.cursor()
        query = '''SELECT COEFF_1, COEFF_2, COEFF_3, COEFF_4, COEFF_5, COEFF_6, COEFF_7
                    FROM {} 
                    WHERE SPECIES_NAME = "{}"'''.format(temp_range.upper(), species_name)
        coeffs = list(cursor.execute(query).fetchall()[0])
        db.commit()
        db.close()
        return coeffs

    def get_species(self, temp, temp_range):
        db = sqlite3.connect(self.db_path)
        cursor = db.cursor()
        if temp_range == 'low': # temp_range == 'low'
            query = '''SELECT SPECIES_NAME FROM {} WHERE TLOW < {}'''.format(temp_range.upper(), temp)
        else: # temp_range == 'high'
            query = '''SELECT SPECIES_NAME FROM {} WHERE THIGH > {}'''.format(temp_range.upper(), temp)
        species = []
        for s in cursor.execute(query).fetchall():
            species.append(s[0])
        db.commit()
        db.close()
        return species

