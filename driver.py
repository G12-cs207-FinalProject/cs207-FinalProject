from chemkin import XmlParser, pckg_xml_path
from chemkin.viz import summary

Ti = [750, 1500, 2500]
# xi = [2.0, 1.0, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0] # specii concentration
xi = [2., 1., .5, 1., 1., 1., .5, 1.] 


xml_parser = XmlParser(pckg_xml_path('rxns_hw5'))
parsed_data_list = xml_parser.parsed_data_list(Ti)
summary.print_reaction_rate(parsed_data_list, xi)
