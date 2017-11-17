from chemkin import pckg_xml_path
from chemkin.preprocessing.parse_xml import XmlParser
from chemkin.viz import summary

# Ti = [750, 1500, 2500]
# xi = [2.0, 1.0, 0.5, 1.0, 1.0]  # specie concentrations for 'rxns_hw5.xml'
# xml_parser = XmlParser(pckg_xml_path('rxns_hw5'))

Ti = [100, 750, 1500, 2500, 5000]
xi = [2., 1., .5, 1., 1., 1., .5, 1.] # specie concentrations 'rxns_reversible.xml'
xml_parser = XmlParser(pckg_xml_path('rxns_reversible'))

parsed_data_list = xml_parser.parsed_data_list(Ti)
summary.print_reaction_rate(parsed_data_list, xi)