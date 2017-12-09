#####################################################
# Additional Feature of the chemkin library 
# 1. Solving the evolution of species concentration
# 2. Visualizing the solver's results
#####################################################

from chemkin import pckg_xml_path
from chemkin.preprocessing.parse_xml import XmlParser
from chemkin.viz import summary


Ti = [900, 2500]
xi = [2., 1., .5, 1., 1., 1., .5, 1.] # specie concentrations
xml_parser = XmlParser(pckg_xml_path('rxns_reversible'))
parsed_data_list = xml_parser.parsed_data_list(Ti)

# Print the species concentration
summary.print_species_concentration(parsed_data_list, xi, end_t=1e-12)

# Plot the species concentration
summary.plot_species_concentration(parsed_data_list, xi, end_t=1e-12)

# Print the time to equilibrium of all reactions
summary.print_time_to_equilibrium(parsed_data_list, xi)

# Plot the time to equilibrium of all reactions
summary.plot_time_to_equilibrium(parsed_data_list, xi)



