from chemkin import pckg_xml_path
from chemkin.preprocessing.parse_xml import XmlParser
from chemkin.viz import summary

# test xml files for irreversible reactions
irrev_xmls = ['rxns_hw5', 'rxns_ideal', 'rxns_neg_A', 'rxns_neg_A_2', 'rxns_neg_k, rxns_non_elementary']
# Ti = [750, 1500, 2500]
# xi = [2.0, 1.0, 0.5, 1.0, 1.0]
# xml_parser = XmlParser(pckg_xml_path(irrev_xmls[0]))

# test xml files for reversible reactions
rev_xmls = ['rxns_reversible, rxns_rev_neg_A, rxns_rev_neg_k, rxns_rev_non_elementary, rxns_rev_inconsistent']
Ti = [100, 750, 1500, 2500, 5000]
xi = [2., 1., .5, 1., 1., 1., .5, 1.] # specie concentrations 'rxns_reversible.xml'
xml_parser = XmlParser(pckg_xml_path(rev_xmls[0]))


parsed_data_list = xml_parser.parsed_data_list(Ti)
summary.print_reaction_rate(parsed_data_list, xi)