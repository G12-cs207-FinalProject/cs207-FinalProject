###############################################################################
# Tests for chemkin.viz.summary module
###############################################################################

import matplotlib.pyplot as plt
from chemkin import pckg_xml_path
from chemkin.preprocessing.parse_xml import XmlParser
from chemkin.viz import summary

def test_print_normal_irreversible():
	Ti = [2500]
	xi = [2.0, 1.0, 0.5, 1.0, 1.0]  # specie concentrations for 'rxns_hw5.xml'
	xml_parser = XmlParser(pckg_xml_path('rxns_hw5'))
	parsed_data_list = xml_parser.parsed_data_list(Ti)
	test_flag = summary.print_reaction_rate(parsed_data_list, xi)
	assert test_flag == 0

def test_print_normal_reversible():
	Ti = [2500]
	xi = [2., 1., .5, 1., 1., 1., .5, 1.] # specie concentrations 'rxns_reversible.xml'
	xml_parser = XmlParser(pckg_xml_path('rxns_reversible'))
	parsed_data_list = xml_parser.parsed_data_list(Ti)
	test_flag = summary.print_reaction_rate(parsed_data_list, xi)
	assert test_flag == 0

def test_print_abnormal_reversible(): 
	Ti = [10000]
	xi = [2., 1., .5, 1., 1., 1., .5, 1.] # specie concentrations 'rxns_reversible.xml'
	xml_parser = XmlParser(pckg_xml_path('rxns_reversible'))
	parsed_data_list = xml_parser.parsed_data_list(Ti)
	test_flag = summary.print_reaction_rate(parsed_data_list, xi)
	assert test_flag == 1

def test_print_species_concentration_normal():
	Ti = [2500]
	xi = [2., 1., .5, 1., 1., 1., .5, 1.] # specie concentrations 'rxns_reversible.xml'
	xml_parser = XmlParser(pckg_xml_path('rxns_reversible'))
	parsed_data_list = xml_parser.parsed_data_list(Ti)
	test_flag = summary.print_species_concentration(parsed_data_list, xi)
	assert test_flag == 0

def test_print_species_concentration_abnormal():
	Ti = [10]
	xi = [2., 1., .5, 1., 1., 1., .5, 1.] # specie concentrations 'rxns_reversible.xml'
	xml_parser = XmlParser(pckg_xml_path('rxns_reversible'))
	parsed_data_list = xml_parser.parsed_data_list(Ti)
	test_flag = summary.print_species_concentration(parsed_data_list, xi)
	assert test_flag == 1


def test_plot_species_concentration_normal():
	Ti = [2500]
	xi = [2., 1., .5, 1., 1., 1., .5, 1.] # specie concentrations 'rxns_reversible.xml'
	xml_parser = XmlParser(pckg_xml_path('rxns_reversible'))
	parsed_data_list = xml_parser.parsed_data_list(Ti)
	test_flag = summary.plot_species_concentration(parsed_data_list, xi)
	assert test_flag == 0

def test_plot_species_concentration_abnormal():
	Ti = [10]
	xi = [2., 1., .5, 1., 1., 1., .5, 1.] # specie concentrations 'rxns_reversible.xml'
	xml_parser = XmlParser(pckg_xml_path('rxns_reversible'))
	parsed_data_list = xml_parser.parsed_data_list(Ti)
	test_flag = summary.plot_species_concentration(parsed_data_list, xi)
	assert test_flag == 1

def test_print_time_to_equilibrium_normal():
	Ti = [2500]
	xi = [2., 1., .5, 1., 1., 1., .5, 1.] # specie concentrations 'rxns_reversible.xml'
	xml_parser = XmlParser(pckg_xml_path('rxns_reversible'))
	parsed_data_list = xml_parser.parsed_data_list(Ti)
	test_flag = summary.print_time_to_equilibrium(parsed_data_list, xi)
	assert test_flag == 0

def test_print_time_to_equilibrium_abnormal():
	Ti = [10]
	xi = [2., 1., .5, 1., 1., 1., .5, 1.] # specie concentrations 'rxns_reversible.xml'
	xml_parser = XmlParser(pckg_xml_path('rxns_reversible'))
	parsed_data_list = xml_parser.parsed_data_list(Ti)
	test_flag = summary.print_time_to_equilibrium(parsed_data_list, xi)
	assert test_flag == 1

def test_plot_time_to_equilibrium_normal():
	Ti = [2500]
	xi = [2., 1., .5, 1., 1., 1., .5, 1.] # specie concentrations 'rxns_reversible.xml'
	xml_parser = XmlParser(pckg_xml_path('rxns_reversible'))
	parsed_data_list = xml_parser.parsed_data_list(Ti)
	test_flag = summary.plot_time_to_equilibrium(parsed_data_list, xi)
	assert test_flag == 0

def test_plot_time_to_equilibrium_abnormal():
	Ti = [10]
	xi = [2., 1., .5, 1., 1., 1., .5, 1.] # specie concentrations 'rxns_reversible.xml'
	xml_parser = XmlParser(pckg_xml_path('rxns_reversible'))
	parsed_data_list = xml_parser.parsed_data_list(Ti)
	test_flag = summary.plot_time_to_equilibrium(parsed_data_list, xi)
	assert test_flag == 1
