from chemkin import pckg_xml_path
from chemkin.preprocessing.parse_xml import XmlParser
from chemkin.viz import summary


Ti = [100, 1500] # temp held constant
xi = [2., 1., .5, 1., 1., 1., .5, 1.] # specie concentrations 'rxns_reversible.xml'
xml_parser = XmlParser(pckg_xml_path('rxns_reversible'))

parsed_data_list = xml_parser.parsed_data_list(Ti)

# summary.print_species_concentration(parsed_data_list, xi)

# summary.plot_species_concentration(parsed_data_list, xi)

# summary.print_time_to_equilibrium(parsed_data_list, xi)

summary.plot_time_to_equilibrium(parsed_data_list, xi)



