import numpy as np
from chemkin import pckg_xml_path
from chemkin.preprocessing.parse_xml import XmlParser
from chemkin.solver.ODEint_solver import ODE_int_solver
from chemkin.viz import summary

#########
#FOR TESTING THE SOLVER - testing for ONE TEMP and ONE SPECIE
Ti = [1500] #temp held constant
xi = [2., 1., .5, 1., 1., 1., .5, 1.] # specie concentrations 'rxns_reversible.xml'
xml_parser = XmlParser(pckg_xml_path('rxns_reversible'))

parsed_data_list = xml_parser.parsed_data_list(Ti)
summary.plot_species_concentration(parsed_data_list, xi)

# species = parsed_data_list[0]['species']
# ki = parsed_data_list[0]['ki']
# sys_vi_p = parsed_data_list[0]['sys_vi_p']
# sys_vi_dp = parsed_data_list[0]['sys_vi_dp']
# is_reversible = parsed_data_list[0]['is_reversible']
# T = parsed_data_list[0]['T']
# b_ki = parsed_data_list[0]['b_ki']

# n_steps = 101
# time_steps = np.linspace(0, 1, n_steps)
# my_solver = ODE_int_solver(T, xi, ki, b_ki, sys_vi_p, sys_vi_dp)
# sol = my_solver.solve(time_steps)
# print(sol.shape)

#########