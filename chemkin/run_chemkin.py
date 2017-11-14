"""
A demo of using chemkin.py to calculate the reaction rates of a system of 
irreversible elementary reactions by reading reaction data stored in XML
files.

"""
import sys
sys.path.append('../')

from reaction.reaction_coefficients import *
from reaction.elementary_rxn import *

import numpy as np

from preprocessing.parse_xml import *

from thermodynamics.thermo import *
#import preprocessing.tests.test_parse_xml



xml_file = './xml-files/rxns_rev.xml'
xml_parser = XmlParser(xml_file)


species, rxn_data_list = xml_parser.load()
n_species = len(species)

xi = [2.0, 1.0, 0.5, 1.0, 1.0] # specii concentration

species_idx_dict = {} # build the dictionary of key = species_name, value = species_index
for i, s in enumerate(species):
	species_idx_dict[s] = i


Ti = [750, 1500, 2500]

for T in Ti:
	sys_vi_p = [] # list of reactant Stoichiometric coefficients in each rxn
	sys_vi_dp = [] # list of product Stoichiometric coefficients in each rxn
	ki = [] # list of reation rate coefficients in each rxn
	is_reversible = None # indicator of the system of reactions being irreversible/reversible

	for rxn_data in rxn_data_list: # 1 rxn per rxn_data
		if rxn_data.type != RxnType.Elementary: 
			raise TypeError('Non-elementary reactions not implemented yet.')

		if is_reversible == None:
			is_reversible = rxn_data.reversible # set the indicator of the system of reactions to be irreversible/reversible

		if rxn_data.reversible != is_reversible: # the system of reactions in the XML file must be all irreversible/reversible
			raise TypeError('The system of reactions are inconsistent in reversibility.')
		
		rxn_id = rxn_data.rxn_id # save id

		rxn_vi_p = np.zeros((n_species,)) # save the Stoichiometric coefficients of the reactants in this rxn
		for s, vi in rxn_data.reactants.items():
			idx = species_idx_dict[s] # get index of the specii
			rxn_vi_p[idx] = vi
		sys_vi_p.append(list(rxn_vi_p))

		rxn_vi_dp = np.zeros((n_species,)) # save the Stoichiometric coefficients of the products in this rxn
		for s, vi in rxn_data.products.items():
			idx = species_idx_dict[s] # get index of the specii
			rxn_vi_dp[idx] = vi
		sys_vi_dp.append(list(rxn_vi_dp))
		
		coef_params = rxn_data.rate_coeff
		if isinstance(coef_params, list):
			if len(coef_params) == 3: # modified arrhenius coef
				A = coef_params[0]
				b = coef_params[1]
				E = coef_params[2]
				ki.append(ModArrheniusCoef(A, b, E, T).get_coef())
			else: # arrhenius coef
				A = coef_params[0]
				E = coef_params[1]
				ki.append(ArrheniusCoef(A, E, T).get_coef())
		else: # const coef
			ki.append(ConstCoef(coef_params).get_coef())

	# print(sys_vi_p)
	# print(sys_vi_dp)
	# print(ki)
	# print(IrrevElemRxn(ki, xi, sys_vi_p, sys_vi_dp))
	if is_reversible == False:
		rxn_rates = IrrevElemRxn(ki, xi, sys_vi_p, sys_vi_dp).reaction_rate()
	else:
		b_ki = Thermo(species, T, ki, sys_vi_p, sys_vi_dp).get_backward_coefs()
		rxn_rates = RevElemRxn(ki, b_ki, xi, sys_vi_p, sys_vi_dp).reaction_rate()

	

	print('------At Temperature', T, 'K------')
	for s, rate in zip(species, rxn_rates):
		print('    ', s, ':', rate)
	print('--------------------------------')
