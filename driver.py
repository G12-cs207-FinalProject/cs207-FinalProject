##################################
# High level use of the library #
##################################

from chemkin import pckg_xml_path
from chemkin.preprocessing.parse_xml import XmlParser
from chemkin.viz import summary
from chemkin.reaction.elementary_rxn import IrreversibleElementaryRxn

# Reversible Reactions
Ti = [100, 750, 1500, 2500, 5000]
xi = [2., 1., .5, 1., 1., 1., .5, 1.] # specie concentrations 'rxns_reversible.xml'
xml_parser = XmlParser(pckg_xml_path('rxns_reversible'))

parsed_data_list = xml_parser.parsed_data_list(Ti)
summary.print_reaction_rate(parsed_data_list, xi)



#########################################
# More user-oriented use of the library #
#########################################
# Irreversible Reactions
Ti2 = [750, 1500, 2500]
xi2 = [2.0, 1.0, 0.5, 1.0, 1.0]  # specie concentrations for 'rxns_hw5.xml'
xml_parser = XmlParser(pckg_xml_path('rxns_hw5'))

parsed_data_list = xml_parser.parsed_data_list(Ti2)

for i, parsed_data in enumerate(parsed_data_list):

    species = parsed_data['species']
    ki = parsed_data['ki']
    sys_vi_p = parsed_data['sys_vi_p']
    sys_vi_dp = parsed_data['sys_vi_dp']
    is_reversible = parsed_data['is_reversible']
    T = parsed_data['T']

    rxn = IrreversibleElementaryRxn(ki, xi2, sys_vi_p, sys_vi_dp)
    progress_rates = rxn.progress_rate()
    rxn_rates = rxn.reaction_rate()

    print('At T = {}: '.format(Ti2[i]))
    print('The progress rates: ')
    print(progress_rates)
    print('The reaction rates: ')
    print(rxn_rates)
    print('==================')