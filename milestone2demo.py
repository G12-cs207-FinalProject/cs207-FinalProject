from chemkin import pckg_xml_path
from chemkin.preprocessing.parse_xml import XmlParser
from chemkin.viz import summary
from chemkin.reaction.elementary_rxn import ElementaryRxn

#########################################
# Demo use of the library #
#########################################

# Reversible Reactions
#Ti = [100, 750, 1500, 2500, 5000]

Ti2 = input("Please enter a list of the temperatures for your species, separated by commas: ")
try:
    Ti2 = [float(i.strip()) for i in Ti2.split(",")]
except ValueError:
    import sys
    sys.exit("The entered temperatures are not valid")

#xi = [2., 1., .5, 1., 1., 1., .5, 1.] # specie concentrations 'rxns_reversible.xml'
xi2 = input("Please enter a list of the concentrations for your species (same order as temperatures and separated by commas): ")

try:
    xi2 = [float(i.strip()) for i in xi2.split(",")]
except ValueError:
    import sys
    sys.exit("The entered concentrations are not valid")
    
curr_path = input("Please input path for the xml with your reaction data: ")
print("\n")

xml_parser = XmlParser(curr_path)

parsed_data_list = xml_parser.parsed_data_list(Ti2)

for i, parsed_data in enumerate(parsed_data_list):

    species = parsed_data['species']
    ki = parsed_data['ki']
    b_ki = parsed_data['b_ki']
    sys_vi_p = parsed_data['sys_vi_p']
    sys_vi_dp = parsed_data['sys_vi_dp']
    is_reversible = parsed_data['is_reversible']
    T = parsed_data['T']

    rxn = ElementaryRxn(ki, b_ki, xi2, sys_vi_p, sys_vi_dp)
    progress_rates = rxn.progress_rate()
    rxn_rates = rxn.reaction_rate()

    print('At T = {}: '.format(Ti2[i]))
    print('The progress rates: ')
    print(progress_rates)
    print('The reaction rates: ')
    print(rxn_rates)
    print('==================')
