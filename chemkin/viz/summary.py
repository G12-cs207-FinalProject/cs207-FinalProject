import os
import os.path
import numpy as np
import matplotlib
matplotlib.use('agg') # must be called before importing pyplot
import matplotlib.pyplot as plt
from chemkin.reaction.elementary_rxn import ElementaryRxn

def print_reaction_rate(parsed_data_list, xi):
	''' Function to print the reaction rates
	'''
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


def print_species_concentration(parsed_data_list, xi, n_steps=101, end_t=1e-12):
	''' Function to print the species concentration at a given time: end_t
	'''
	test_flag = 0 # species_concentrations can be printed
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

		species_concentration = ElementaryRxn(ki, b_ki, xi, sys_vi_p, sys_vi_dp).species_concentration(T, end_t, n_steps)
		
		# print(np.min(sol), np.max(sol))

		print('------At Temperature', T, 'K------')
		print('Specie Concentration at time = {}'.format(end_t))
		for i, s in enumerate(species):
			print('  {}: start = {}, end = {}'.format(s, xi[i], species_concentration[i]))
		print('--------------------------------\n')

	return test_flag


def plot_species_concentration(parsed_data_list, xi, n_steps=101, end_t=1e-12):
	''' Function to plot the evolution of species concentration from start to an end time: end_t
	'''
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
		species_concentration_evolution = ElementaryRxn(ki, b_ki, xi, sys_vi_p, sys_vi_dp).species_concentration_evolution(T, end_t, n_steps)
		

		f1, ax1 = plt.subplots(1, 1)
		for i, s in enumerate(species):
			# Plot the evolution of all species' concentration
			ax1.plot(time_steps, species_concentration_evolution[:, i], label='{}'.format(s))
		
		ax1.legend()
		ax1.set_ylim(0, np.max(species_concentration_evolution)+1)
		ax1.set_xlabel('Time (sec)')
		ax1.set_ylabel('Species Concentration')
		ax1.set_title('Evolution of Species Concentration over Time')

		BASE_DIR = os.path.dirname(os.path.abspath(__file__))
		image_dir = os.path.join(BASE_DIR, 'examples')
		if not os.path.exists(image_dir):
			os.makedirs(image_dir)
		img_path = os.path.join(image_dir, 'evolution_{}K.png'.format(T))
		f1.savefig(img_path)

	return test_flag


def print_time_to_equilibrium(parsed_data_list, xi, n_steps=101):
	''' Function to print the time to equilibrium of all reactions
	'''
	test_flag = 0 # time_to_equilibrium can be printed
	for parsed_data in parsed_data_list:

		ki = parsed_data['ki']
		sys_vi_p = parsed_data['sys_vi_p']
		sys_vi_dp = parsed_data['sys_vi_dp']
		T = parsed_data['T']
		
		b_ki = parsed_data['b_ki']
		if str(b_ki) == 'Not Defined':
			test_flag = 1 # time_to_equilibrium cannot be printed because T is not in some specie's temperature range

			print('------At Temperature', T, 'K------')
			print('Backward reaction coefficients not defined: T={} is not in some specie\'s temperature range.'.format(T))
			print('--------------------------------\n')
			continue

		end_t, critical_t, overall_critical_t = ElementaryRxn(ki, b_ki, xi, sys_vi_p, sys_vi_dp).time_to_equilibrium(T, n_steps)
		time_steps = np.linspace(0, end_t, n_steps)

		print('------At Temperature', T, 'K------')
		print('Time to Equilibrium (end_t = {})'.format(end_t))
		for i, t in enumerate(critical_t):
			if t == -100:
				t = 'Not yet reached equilibrium.'
			print('  Reaction #{}: {}'.format(i, t))

		print()
		if overall_critical_t == -100:
			overall_critical_t = 'The system has not yet reached equilibrium.'
		print('Overall Time to Equilibrium: {}'.format(overall_critical_t))
		print('--------------------------------\n')

	return test_flag



def plot_time_to_equilibrium(parsed_data_list, xi, n_steps=101):
	''' Function to plot the time to equilibrium of all reactions
	'''
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

		end_t, critical_t, overall_critical_t = ElementaryRxn(ki, b_ki, xi, sys_vi_p, sys_vi_dp).time_to_equilibrium(T, n_steps)
		time_steps = np.linspace(0, end_t, n_steps)
		
		# Plot Log-scale Time to Equilibrium
		f2, ax2 = plt.subplots(1, 1)
		xpos = np.arange(0, len(parsed_data['equations']))
		
		# Take Log-transform of each rxn's time-to-equilibrium
		yval = np.copy(critical_t)
		for j, v in enumerate(yval):
			if v == -100:
				yval[j] = 0
			else:
				yval[j] = np.log10(v+1)
		labels = parsed_data['equations']

		for j, (x, y) in enumerate(zip(xpos, yval)):
			ax2.bar(x, y, width=0.5, bottom=0.0, align='center', alpha=0.6, label=labels[j]+':'+str('{:0.3e}'.format(y)))
		
		ax2.legend(fontsize=9)
		ax2.set_xlabel('Reaction equation')
		ax2.set_ylabel('log(Time + 1)')
		ax2.set_xticks(np.arange(0, len(parsed_data['equations'])))
		ax2.set_title('Log-scale Time to Equilibrium')

		BASE_DIR = os.path.dirname(os.path.abspath(__file__))
		image_dir = os.path.join(BASE_DIR, 'examples')
		if not os.path.exists(image_dir):
			os.makedirs(image_dir)
		img_path = os.path.join(image_dir, 'Time_to_Equilibrium_{}K.png'.format(T))
		f2.savefig(img_path)

	return test_flag



