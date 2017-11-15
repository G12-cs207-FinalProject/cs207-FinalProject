from chemkin.preprocessing.parse_xml import XmlParser
from chemkin.viz import summary

Ti = [100, 750, 1500, 2500, 5000] 
xi = [2.0, 1.0, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0] # specii concentration

xml_file = './chemkin/xml-files/rxns_reversible.xml'
xml_parser = XmlParser(xml_file)
parsed_data_list = xml_parser.parsed_data_list(Ti)
summary.print_reaction_rate(parsed_data_list, xi)