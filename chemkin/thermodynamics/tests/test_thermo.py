###############################################################################
# Tests for chemkin.thermodynamics.thermo module
###############################################################################

from chemkin import pckg_xml_path
from chemkin.chemkin_errors import ChemKinError
from chemkin.preprocessing.parse_xml import XmlParser

def test_get_backward_coefs_normal():
	Ti = [750, 1500]
	xml_parser = XmlParser(pckg_xml_path('rxns_reversible'))
	parsed_data_list = xml_parser.parsed_data_list(Ti) # calling Thermo().get_backward_coefs()
	assert parsed_data_list[0]['T'] == 750
	assert parsed_data_list[1]['T'] == 1500

def test_get_backward_coefs_high_range_err():
	Ti = [5000]
	xml_parser = XmlParser(pckg_xml_path('rxns_reversible'))
	parsed_data_list = xml_parser.parsed_data_list(Ti) # calling Thermo().get_backward_coefs()
	assert parsed_data_list[0]['b_ki'] == 'Not Defined'

def test_get_backward_coefs_low_range_err():
	Ti = [100]
	xml_parser = XmlParser(pckg_xml_path('rxns_reversible'))
	parsed_data_list = xml_parser.parsed_data_list(Ti) # calling Thermo().get_backward_coefs()
	assert parsed_data_list[0]['b_ki'] == 'Not Defined'