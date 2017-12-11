"""
Tests for the ODE Solver ODEint_solver.py module
"""

import sys
import numpy as np
from chemkin import pckg_xml_path
from chemkin.preprocessing.parse_xml import XmlParser
from chemkin.solver.ODEint_solver import ODE_int_solver
from chemkin.reaction.elementary_rxn import ElementaryRxn

def test_ODE_solver_functionality():
    
    """
    Test basic ODE solver functionality.
    """
    
    Ti = [1500]
    xi = [2., 1., .5, 1., 1., 1., .5, 1.] # specie concentrations 'rxns_reversible.xml'
    xml_parser = XmlParser(pckg_xml_path('rxns_reversible'))
    parsed_data_list = xml_parser.parsed_data_list(Ti)
    #print(parsed_data_list)
    for parsed_data in parsed_data_list:
        
        species = parsed_data['species']
        #print(species)
        ki = parsed_data['ki']
        sys_vi_p = parsed_data['sys_vi_p']
        sys_vi_dp = parsed_data['sys_vi_dp']
        T= parsed_data['T']
        b_ki=parsed_data['b_ki']

        rxn = ElementaryRxn(ki, b_ki, xi, sys_vi_p, sys_vi_dp)
        my_solver = ODE_int_solver(T, rxn)
        
        error_msg = "Solver is not initiating properly"
        assert my_solver != None and isinstance(my_solver,ODE_int_solver),error_msg
 
def test_ODE_solver_solve():
    
    """
    Tests ODE solver solve() method.
    """
    
    Ti = [1500]
    xi = [2., 1., .5, 1., 1., 1., .5, 1.] # specie concentrations 'rxns_reversible.xml'
    xml_parser = XmlParser(pckg_xml_path('rxns_reversible'))
    parsed_data_list = xml_parser.parsed_data_list(Ti)

    for parsed_data in parsed_data_list:
        
        species = parsed_data['species']
        ki = parsed_data['ki']
        sys_vi_p = parsed_data['sys_vi_p']
        sys_vi_dp = parsed_data['sys_vi_dp']
        T= parsed_data['T']
        b_ki=parsed_data['b_ki']       

        rxn = ElementaryRxn(ki, b_ki, xi, sys_vi_p, sys_vi_dp)
        my_solver = ODE_int_solver(T, rxn)

        time_int = time_steps = np.linspace(0, 100, 101)
        
        s , t_c , t_o =my_solver.solve(time_int)
        
        error_msg1 = "The dimensions of input do not match output"
        
        assert s.shape == (len(s),len(xi)),error_msg1
        
        assert len(t_c)==len(ki) and len(t_c)==len(b_ki),error_msg1
        
        #error_msg2 = "Reaction has not reached an equilibrium"
        #assert all(i>0 for i in t_c) and t_o >0, error_msg2
