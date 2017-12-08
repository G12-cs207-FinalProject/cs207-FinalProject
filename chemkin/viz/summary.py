import numpy as np
import matplotlib
matplotlib.use('agg') # must be called before importing pyplot
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

def plot_species_concentration(parsed_data_list, xi, n_steps=101, end_t=1e-12):
	test_flag = 0 # species_concentrations can be plotted
	for parsed_data in parsed_data_list:

		species = parsed_data['species']
		ki = parsed_data['ki']
		sys_vi_p = parsed_data['sys_vi_p']
		sys_vi_dp = parsed_data['sys_vi_dp']
		T = parsed_data['T']
		
		b_ki = parsed_data['b_ki']
		if str(b_ki) == 'Not Defined':
			test_flag = 1 # species_concentrations cannot be plotted because T is not in some specie's temperature range

			print('------At Temperature', T, 'K------')
			print('Backward reaction coefficients not defined: T={} is not in some specie\'s temperature range.'.format(T))
			print('--------------------------------\n')
			continue

		time_steps = np.linspace(0, end_t, n_steps)
		my_solver = ODE_int_solver(T, xi, ki, b_ki, sys_vi_p, sys_vi_dp)
		sol, _, _ = my_solver.solve(time_steps)
		
		# print(np.min(sol), np.max(sol))

		f1, ax1 = plt.subplots(1, 1)
		print('------At Temperature', T, 'K------')
		print('Specie Concentration at the Start and the End')
		for i, s in enumerate(species):
			# Plot the evolution of all species' concentration
			ax1.plot(time_steps, sol[:, i], label='{}'.format(s))
			print('  {}: start = {}, end = {}'.format(s, sol[0, i], sol[-1, i]))
		print('--------------------------------\n')

		ax1.legend()
		ax1.set_ylim(0, np.max(sol)+1)
		ax1.set_xlabel('Time (sec)')
		ax1.set_ylabel('Species Concentration')
		ax1.set_title('Evolution of Species Concentration over Time')
		f1.savefig('evolution_{}K.png'.format(T))

	return test_flag


def plot_time_to_equilibrium(parsed_data_list, xi, n_steps=101, end_t=100):
	test_flag = 0 # time_to_equilibrium can be plotted
	for parsed_data in parsed_data_list:

		ki = parsed_data['ki']
		sys_vi_p = parsed_data['sys_vi_p']
		sys_vi_dp = parsed_data['sys_vi_dp']
		T = parsed_data['T']
		
		b_ki = parsed_data['b_ki']
		if str(b_ki) == 'Not Defined':
			test_flag = 1 # time_to_equilibrium cannot be plotted because T is not in some specie's temperature range

			print('------At Temperature', T, 'K------')
			print('Backward reaction coefficients not defined: T={} is not in some specie\'s temperature range.'.format(T))
			print('--------------------------------\n')
			continue

		time_steps = np.linspace(0, end_t, n_steps)
		my_solver = ODE_int_solver(T, xi, ki, b_ki, sys_vi_p, sys_vi_dp)
		_, critical_t, overall_critical_t = my_solver.solve(time_steps)

		print('------At Temperature', T, 'K------')
		print('Time to Equilibrium')
		for i, t in enumerate(critical_t):
			if t == -2:
				t = 'Not yet reached equilibrium.'
			print('  Reaction #{}: {}'.format(i, t))

		print()
		if overall_critical_t == -2:
			overall_critical_t = 'The system has not yet reached equilibrium.'
		print('Overall Time to Equilibrium: {}'.format(overall_critical_t))
		print('--------------------------------\n')

		# Plot Log-scale Time to Equilibrium
		f2, ax2 = plt.subplots(1, 1)
		ax2.bar(np.arange(len(critical_t)), np.log(critical_t+2))
		ax2.set_xlabel('Reaction equation')
		ax2.set_ylabel('log(Time)')
		#print(parsed_data['equations'])
		ax2.xaxis.set_ticks(np.arange(0,11,1))
		ax2.set_xticklabels(zip([i for i in range(1,len(parsed_data['equations'])+1)],parsed_data['equations']),minor=False)
		for tick in ax2.get_xticklabels():
			tick.set_rotation(90)
		ax2.locator_params(tight=False,nbins=11)
		ax2.set_title('Log-scale Time to Equilibrium')
		f2.savefig('Time to Equilibrium_{}K.png'.format(T))

	return test_flag


# def plot_species_concentration(parsed_data_list, xi):
#     test_flag = 0 # species_concentrations can be plotted
#     for parsed_data in parsed_data_list:
        
#         species = parsed_data['species']
#         ki = parsed_data['ki']
#         sys_vi_p = parsed_data['sys_vi_p']
#         sys_vi_dp = parsed_data['sys_vi_dp']
#         T = parsed_data['T']
		
#         b_ki = parsed_data['b_ki']
#         if str(b_ki) == 'Not Defined':
#             test_flag = 1 # rspecies_concentrations cannot be plotted because T is not in some specie's temperature range

#             print('------At Temperature', T, 'K------')
#             print('Backward reaction coefficients not defined: T={} is not in some specie\'s temperature range.'.format(T))
#             print('--------------------------------')
#             continue

#         n_steps = 101
#         time_steps = np.linspace(0, 0.000000000001, n_steps)
#         my_solver = ODE_int_solver(T, xi, ki, b_ki, sys_vi_p, sys_vi_dp)
#         sol, critical_t = my_solver.solve(time_steps)
#         print(np.min(sol), np.max(sol))

#         print('Critical t: {}'.format(critical_t))
        
#         for i, s in enumerate(species):
#             plt.plot(time_steps, sol[:, i], label='{}'.format(s))
            
#         plt.legend()
#         plt.ylim(0, 3)
#         plt.savefig('evolution_{}.png'.format(T))

#         plt.show()

#     return test_flag




