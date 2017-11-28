###############################################################################
# Tests for RxnData and XmlParser classes.

# Class methods are used on different XML files located in xml-files directory.
# Each version of the reactions XML file tests different features / cases:
#   rxns_ideal: Contains data in correct manner (the "ideal" example file).
#   rxns_neg_A: Same as rxns_ideal except for negative A rateCoeff in first
#       reaction.
###############################################################################

from pytest import approx

from chemkin import pckg_xml_path
from chemkin.chemkin_errors import ChemKinError
from chemkin.preprocessing.parse_xml import RxnType, XmlParser


def test_parse_basic_functionality ():
    """Ensures number of reactions returned is correct and that attributes
    of the <reaction> element in the XML file are parsed correctly.
    """
    xml = XmlParser(pckg_xml_path('rxns_ideal.xml'))
    species, rxns = xml.load()

    # Correct number of reactions returned.
    err_msg = 'Expected 2 reactions but received {}.'.format(len(rxns))
    assert len(rxns) == 2, err_msg

    # Attributes of <reaction> element parsed properly.
    err_msg = 'reversible attribute not parsed properly.'
    assert rxns[0].reversible == False, err_msg
    assert rxns[1].reversible == False, err_msg

    err_msg = 'rxn_id attribute not parsed properly.'
    assert rxns[0].rxn_id == 'reaction01', err_msg
    assert rxns[1].rxn_id == 'reaction02', err_msg

    err_msg = 'type attribute not parsed properly.'
    assert rxns[0].type == RxnType.Elementary, err_msg
    assert rxns[1].type == RxnType.Elementary, err_msg


def test_parse_reactants_products ():
    """Ensures reactants and products parsed correctly."""
    xml = XmlParser(pckg_xml_path('rxns_ideal.xml'))
    species, rxns = xml.load()

    err_msg = 'reactants not parsed correctly.'
    assert rxns[0].reactants == {'H':1, 'O2':1}, err_msg
    assert rxns[1].reactants == {'H2':1, 'O':1}, err_msg

    err_msg = 'products not parsed correctly.'
    assert rxns[0].products == {'OH':1, 'O':1}, err_msg
    assert rxns[1].products == {'OH':1, 'H':1}, err_msg


def test_parse_rxn_coeff ():
    """Ensures reaction coefficients are what we expect them to be."""
    xml = XmlParser(pckg_xml_path('rxns.xml'))
    species, rxns = xml.load()

    rxn1_coeff = rxns[0].rate_coeff
    rxn2_coeff = rxns[1].rate_coeff
    rxn3_coeff = rxns[2].rate_coeff

    # reaction01
    assert rxn1_coeff[0] == approx(3.52e+10)
    assert rxn1_coeff[1] == approx(7.14e+04)

    # reaction02
    assert rxn2_coeff[0] == approx(5.06e-2)
    assert rxn2_coeff[1] == approx(2.7)
    assert rxn2_coeff[2] == approx(2.63e+04)

    # reaction03
    assert rxn3_coeff == approx(1.0e+03)


def test_rxndata_equation ():
    """Ensures equation representation of RxnData is correct."""
    xml = XmlParser(pckg_xml_path('rxns_ideal.xml'))
    species, rxns = xml.load()

    err_msg = 'equation() method result different than expected.'
    expected_0 = 'H + O2 [=] O + OH'
    print(rxns[0].equation())
    assert rxns[0].equation() == expected_0, err_msg

    err_msg = 'equation() method result different than expected.'
    expected_1 = 'H2 + O [=] H + OH'
    assert rxns[1].equation() == expected_1, err_msg


def test_badparse_negative_A_arr ():
    """Ensures ChemKinError raised when A coefficient for one of the
    reactions is negative.
    """
    xml = XmlParser(pckg_xml_path('rxns_neg_A_2'))
    try:
        rxns = xml.load()
    except ChemKinError as err:
        assert type(err) == ChemKinError
        assert str(err).find(
              'A coeff < 0 in reaction with id = reaction02') != -1


def test_badparse_negative_A_modarr ():
    """Ensures ChemKinError raised when A coefficient for one of the
    reactions is negative.
    """
    xml = XmlParser(pckg_xml_path('rxns_neg_A'))
    try:
        rxns = xml.load()
    except ChemKinError as err:
        assert type(err) == ChemKinError
        assert str(err).find(
              'A coeff < 0 in reaction with id = reaction01') != -1
