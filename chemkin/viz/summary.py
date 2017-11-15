from chemkin.reaction.elementary_rxn import IrreversibleElementaryRxn
from chemkin.reaction.elementary_rxn import ReversibleElementaryRxn

def print_reaction_rate(parsed_data_list, xi):
	for parsed_data in parsed_data_list:

		species = parsed_data['species']
		ki = parsed_data['ki']
		sys_vi_p = parsed_data['sys_vi_p']
		sys_vi_dp = parsed_data['sys_vi_dp']
		is_reversible = parsed_data['is_reversible']
		T = parsed_data['T']
	
		if is_reversible == False:
			rxn_rates = IrreversibleElementaryRxn(ki, xi, sys_vi_p, sys_vi_dp).reaction_rate()
		
		else:
			b_ki = parsed_data['b_ki']
			if str(b_ki) == 'Not Defined':
				print('------At Temperature', T, 'K------')
				print('Backward reaction coefficients not defined: T={} is not in some specie\'s temperature range.'.format(T))
				print('--------------------------------')
				continue
			
			rxn_rates = ReversibleElementaryRxn(ki, b_ki, xi, sys_vi_p, sys_vi_dp).reaction_rate()

		
		print('------At Temperature', T, 'K------')
		for s, rate in zip(species, rxn_rates):
			print('    ', s, ':', rate)
		print('--------------------------------')