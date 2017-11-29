import numpy as np
import matplotlib.pyplot as plt
from chemkin.reaction.elementary_rxn import ElementaryRxn
from chemkin.solver.ODEint_solver import ODE_int_solver

def print_reaction_rate(parsed_data_list, xi):
	test_flag = 0 # reation rates can be printed
	
	for parsed_data in parsed_data_list:

		species = parsed_data['species']
		ki = parsed_data['ki']
		sys_vi_p = parsed_data['sys_vi_p']
		sys_vi_dp = parsed_data['sys_vi_dp']
		T = parsed_data['T']

		b_ki = parsed_data['b_ki']
		if str(b_ki) == 'Not Defined':
			test_flag = 1 # reaction rates cannot be printed because T is not in some specie's temperature range

			print('------At Temperature', T, 'K------')
			print('Backward reaction coefficients not defined: T={} is not in some specie\'s temperature range.'.format(T))
			print('--------------------------------')
			continue

		rxn_rates = ElementaryRxn(ki, b_ki, xi, sys_vi_p, sys_vi_dp).reaction_rate()

		
		print('------At Temperature', T, 'K------')
		for s, rate in zip(species, rxn_rates):
			print('    ', s, ':', rate)
		print('--------------------------------')

	return test_flag

def plot_species_concentration(parsed_data_list, xi):
	test_flag = 0 # species_concentrations can be plotted
	for parsed_data in parsed_data_list:

		species = parsed_data['species']
		ki = parsed_data['ki']
		sys_vi_p = parsed_data['sys_vi_p']
		sys_vi_dp = parsed_data['sys_vi_dp']
		T = parsed_data['T']
		
		b_ki = parsed_data['b_ki']
		if str(b_ki) == 'Not Defined':
			test_flag = 1 # rspecies_concentrations cannot be plotted because T is not in some specie's temperature range

			print('------At Temperature', T, 'K------')
			print('Backward reaction coefficients not defined: T={} is not in some specie\'s temperature range.'.format(T))
			print('--------------------------------')
			continue

		n_steps = 101
		time_steps = np.linspace(0, 10, n_steps)
		print('hhhhhh', xi)
		my_solver = ODE_int_solver(T, xi, ki, b_ki, sys_vi_p, sys_vi_dp)
		sol, critical_t = my_solver.solve(time_steps)
		print(np.min(sol), np.max(sol))

		print('Critical t: {}'.format(critical_t))

		for i, s in enumerate(species):
			plt.plot(time_steps, sol[:, i], label='{}'.format(s))

		plt.legend()
		plt.ylim(0, 10)
		plt.savefig('evolution_{}.png'.format(T))
		plt.show()

	return test_flag

