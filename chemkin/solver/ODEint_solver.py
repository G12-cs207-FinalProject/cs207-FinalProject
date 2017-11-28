"""
The ODEint_solver.py contains a class ODE_int_solver which will return time evolution
 of a specie i over specified time range t.
"""

from chemkin.reaction.elementary_rxn import ReversibleElementaryRxn
import numpy as np
#########
#FOR TESTING THE SOLVER - testing for ONE TEMP and ONE SPECIE
from chemkin import pckg_xml_path
from chemkin.preprocessing.parse_xml import XmlParser

temp = [750] #temp held constant
xi = [2., 1., .5, 1., 1., 1., .5, 1.] # specie concentrations 'rxns_reversible.xml'

xml_parser = XmlParser(pckg_xml_path('rxns_reversible'))

parsed_data_list = xml_parser.parsed_data_list(temp)

species = parsed_data_list[0]['species']
ki = parsed_data_list[0]['ki']
sys_vi_p = parsed_data_list[0]['sys_vi_p']
sys_vi_dp = parsed_data_list[0]['sys_vi_dp']
is_reversible = parsed_data_list[0]['is_reversible']
T = parsed_data_list[0]['T']
b_ki = parsed_data_list[0]['b_ki']
"""
if is_reversible == False:
	rxn_rates = IrreversibleElementaryRxn(ki, xi, sys_vi_p, sys_vi_dp).reaction_rate()
		
else:
	b_ki = parsed_data['b_ki']
        if str(b_ki) == 'Not Defined':
		test_flag = 1 # reaction rates cannot be printed because T is not in some specie's temperature range 
		print('------At Temperature', T, 'K------')
		print('Backward reaction coefficients not defined: T={} is not in some specie\'s temperature range.'.format(T))
		print('--------------------------------')
		continue
"""
rxn_rates = ReversibleElementaryRxn(ki, b_ki, xi, sys_vi_p, sys_vi_dp).reaction_rate()[0]
print(species)
print(rxn_rates)
print(xi)
#########

class ODE_int_solver():
    #Integrates the time evolution of a concentration of a specie i over specified time range
    
    def __init__(self,temp,xi,ki, b_ki, sys_vi_p, sys_vi_dp):
        
        self.temp = temp #temperature will be held constant in the calculation
        
        self.xi = xi #need initial concentration of ALL species
        self.ki = ki
        self.b_ki = b_ki
        self.sys_vi_p = sys_vi_p
        self.sys_vi_dp = sys_vi_dp
        
    def solve(self,time_int):
        
        #for starters, focus on time evolution of ALL species
        
        from scipy.integrate import odeint
        
        rnx = ReversibleElementaryRxn(self.ki, self.b_ki, self.xi, self.sys_vi_p, self.sys_vi_dp)
        
        i = 0
        def rxn_rate(x,t):
            nonlocal i
            rnx.xi=x
            i+=1
            print(i,x)
            return rnx.reaction_rate()
            
        #func = ReversibleElementaryRxn(ki, b_ki, xi, sys_vi_p, sys_vi_dp).reaction_rate()

        #need to add t to reaction_rate method!!!!
        
        sol = odeint(func=rxn_rate,y0=self.xi,t=time_int)
        
        
        return sol
    
my_solver = ODE_int_solver(T,xi,ki, b_ki, sys_vi_p, sys_vi_dp)
my_solver.solve(np.linspace(0.00000001,0.0000001 ,10))