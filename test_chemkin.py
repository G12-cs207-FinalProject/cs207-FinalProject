"""
Test suite for the chemkin.py module

"""
import numpy as np
import pytest
from pytest import approx
import warnings
from reaction_coefs import *
import chemkin as ck

# tests for RxnCoef() base class
def test_base_get_coef():
    try:
        RxnCoef().get_coef()
    except NotImplementedError as err:
        assert(type(err) == NotImplementedError)


# tests for ConstCoef() class
def test_const_get_coef_result():
    assert ConstCoef(100).get_coef() == 100

def test_const_get_coef_neg_k():
    try:
        ConstCoef(-10).get_coef()
    except ValueError as err:
        assert(type(err) == ValueError)


# tests for ArrheniusCoef() class
def test_arr_get_coef_result():
    assert ArrheniusCoef(A=np.power(10, 4), E=100, T=300, R=8.314).get_coef() == 9607.0007471344579

def test_arr_get_coef_neg_T():
    try:
        ArrheniusCoef(A=np.power(10, 4), E=100, T=-300, R=8.314).get_coef()
    except ValueError as err:
        assert (type(err) == ValueError)

def test_arr_get_coef_neg_A():
    try:
        ArrheniusCoef(A=-np.power(10, 4), E=100, T=300, R=8.314).get_coef()
    except ValueError as err:
        assert (type(err) == ValueError)

def test_arr_get_coef_neg_R():
    try:
        ArrheniusCoef(A=np.power(10, 4), E=100, T=300, R=-8.314).get_coef()
    except ValueError as err:
        assert (type(err) == ValueError)

def test_arr_get_coef_overflow():
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        try:
            ArrheniusCoef(A=np.power(1, 256), E=-10000000., T=300, R=8.314).get_coef()
        except OverflowError as err:
            assert(str(err) == "The result is too large/small.")


# tests for ModArrheniusCoef() class
def test_modarr_get_coef_result():
    assert ModArrheniusCoef(A=np.power(10, 4), b=0.5, E=100, T=300, R=8.314).get_coef() == 166398.13402389048

def test_modarr_get_coef_neg_T():
    try:
        ModArrheniusCoef(A=np.power(10, 4), b=0.5, E=100, T=-300, R=8.314).get_coef()
    except ValueError as err:
        assert (type(err) == ValueError)

def test_modarr_get_coef_neg_A():
    try:
        ModArrheniusCoef(A=-np.power(10, 4), b=0.5, E=100, T=300, R=8.314).get_coef()
    except ValueError as err:
        assert (type(err) == ValueError)

def test_modarr_get_coef_neg_R():
    try:
        ModArrheniusCoef(A=np.power(10, 4), b=0.5, E=100, T=300, R=-8.314).get_coef()
    except ValueError as err:
        assert (type(err) == ValueError)

def test_modarr_get_coef_complex_b():
    try:
        ModArrheniusCoef(A=np.power(10, 4), b=(3+1j), E=100, T=300, R=8.314).get_coef()
    except ValueError as err:
        assert (type(err) == ValueError)

def test_modarr_get_coef_overflow():
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        try:
            ModArrheniusCoef(A=np.power(10, 4), b=256.0, E=100, T=300, R=8.314).get_coef()
        except OverflowError as err:
            assert(str(err) == "The result is too large/small.")
            
import rxns as rx
import io
import sys

# tests for Rxn() base class
def test_Rxn_progress_rate_not_implemented():
    try:
        reac1 = rx.Rxn([10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]], [[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
        reac1.progress_rate()
    except NotImplementedError as err:
        assert(type(err) == NotImplementedError)


def test_Rxn_reaction_rate_not_implemented():
    try:
        reac1 = rx.Rxn([10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]], [[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
        reac1.progress_rate()
    except NotImplementedError as err:
        assert(type(err) == NotImplementedError)


def test_Rxn_len_result():
    reac1 = rx.Rxn([10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]], [[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    assert len(reac1) == 2


def test_Rxn_repr_result():
    old_stdout = sys.stdout
    capturedOutput = io.StringIO()
    sys.stdout = capturedOutput
    reac1 = rx.Rxn([10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]], [[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    print(reac1)
    sys.stdout = old_stdout
    assert capturedOutput.getvalue() == "Rxn(ki=[10, 10], xi=[1.0, 2.0, 1.0], vi_p=[[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]], vi_dp=[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])\n"


# tests for IrrevElemRxn() class
def test_IrrevElem_len_result():
    reac1 = rx.IrrevElemRxn([10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]], [[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    assert len(reac1) == 2


def test_IrrevElem_repr_result():
    old_stdout = sys.stdout
    capturedOutput = io.StringIO()
    sys.stdout = capturedOutput
    reac1 = rx.IrrevElemRxn(ki=[10, 10], xi=[1.0, 2.0, 1.0], vi_p=[[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]], vi_dp=[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])

    print(reac1)
    sys.stdout = old_stdout
    assert capturedOutput.getvalue() == "IrrevElemRxn(ki=[10, 10], xi=[1.0, 2.0, 1.0], vi_p=[[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]], vi_dp=[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])\n"


def test_IrrevElem_progress_rate_result():
    reac1 = rx.IrrevElemRxn([10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]],
                            [[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    assert (reac1.progress_rate() == [40.0, 10.0]).all()


def test_IrrevElem_progress_rate_neg_ki():
    reac1 = rx.IrrevElemRxn([-10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    try:
        reac1.progress_rate()
    except ValueError as err:
        assert(type(err) == ValueError)

        
def test_IrrevElem_progress_rate_neg_xi():
    reac1 = rx.IrrevElemRxn([10, 10], [-1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    try: 
        reac1.progress_rate()
    except ValueError as err:
        assert(type(err) == ValueError)


def test_IrrevElem_progress_rate_dim_compat():
    reac1 = rx.IrrevElemRxn([10, 10], [1.0, 2.0, 1.0], [[ 2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    try:
        reac1.progress_rate()
    except ValueError as err:
        assert(type(err) == ValueError)
    

def test_IrrevElem_reaction_rate_result():
    reac1 = rx.IrrevElemRxn([10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])

    assert (reac1.reaction_rate() == [-60., -70., 70.]).all()


def test_IrrevElem_reaction_rate_neg_ki():
    reac1 = rx.IrrevElemRxn([-10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    try:
        reac1.reaction_rate()
    except ValueError as err:
        assert(type(err) == ValueError)
        

def test_IrrevElem_reaction_rate_neg_xi():
    reac1 = rx.IrrevElemRxn([10, 10], [-1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    try:
         reac1.reaction_rate()
    except ValueError as err:
        assert(type(err) == ValueError)


def test_IrrevElem_reaction_rate_dim_compat():
    reac1 = rx.IrrevElemRxn([10, 10], [1.0, 2.0, 1.0], [[2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    try:
        reac1.reaction_rate()
    except ValueError as err:
        assert(type(err) == ValueError)


###############################################################################
# Tests for RxnData and XmlParser classes.

# Class methods are used on different XML files located in xml-files directory.
# Each version of the reactions XML file tests different features / cases:
#   rxns_ideal: Contains data in correct manner (the "ideal" example file).
#   rxns_neg_A: Same as rxns_ideal except for negative A rateCoeff in first
#       reaction.
###############################################################################

def test_parse_basic_functionality ():
    """ Ensures number of reactions returned is correct and that attributes
    of the <reaction> element in the XML file are parsed correctly.
    """
    xml = ck.XmlParser('xml-files/rxns_ideal.xml')
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
    assert rxns[0].type == ck.RxnType.Elementary, err_msg
    assert rxns[1].type == ck.RxnType.Elementary, err_msg


def test_parse_reactants_products ():
    """ Ensures reactants and products parsed correctly. """
    xml = ck.XmlParser('xml-files/rxns_ideal.xml')
    species, rxns = xml.load()

    err_msg = 'reactants not parsed correctly.'
    assert rxns[0].reactants == {'H':1, 'O2':1}, err_msg
    assert rxns[1].reactants == {'H2':1, 'O':1}, err_msg

    err_msg = 'products not parsed correctly.'
    assert rxns[0].products == {'OH':1, 'O':1}, err_msg
    assert rxns[1].products == {'OH':1, 'H':1}, err_msg


def test_parse_rxn_coeff ():
    """ Ensures reaction coefficients are what we expect them to be. """
    xml = ck.XmlParser('xml-files/rxns_ideal.xml')
    species,rxns = xml.load()

    rxn1_coeff = rxns[0].rate_coeff
    rxn2_coeff = rxns[1].rate_coeff

    # A
    assert rxn1_coeff[0] == approx(3.52e+10)
    assert rxn2_coeff[0] == approx(5.06e-2)

    # b
    assert rxn1_coeff[1] == approx(-0.7)
    assert rxn2_coeff[1] == approx(2.7)

    # E
    assert rxn1_coeff[2] == approx(7.14e+04)
    assert rxn2_coeff[2] == approx(2.63e+04)


def test_rxndata_equation ():
    """ Ensures equation representation of RxnData is correct. """
    xml = ck.XmlParser('xml-files/rxns_ideal.xml')
    species,rxns = xml.load()

    err_msg = 'equation() method result different than expected.'
    expected_0 = 'H + O2 [=] O + OH'
    print(rxns[0].equation())
    assert rxns[0].equation() == expected_0, err_msg

    err_msg = 'equation() method result different than expected.'
    expected_1 = 'H2 + O [=] H + OH'
    assert rxns[1].equation() == expected_1, err_msg


def test_badparse_negative_A ():
    """ Ensures ChemKinError raised when A coefficient for one of the
    reactions is negative.
    """
    xml = ck.XmlParser('xml-files/rxns_neg_A.xml')

    with pytest.raises(ck.ChemKinError, message='Expecting ChemKinError '
                                                'because of negative A coeff.'):
        rxns = xml.load()