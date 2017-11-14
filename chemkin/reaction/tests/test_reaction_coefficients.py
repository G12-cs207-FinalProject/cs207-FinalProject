"""
Test suite for the reaction_coefficients.py module

"""

import warnings

import numpy as np

from chemkin.reaction.reaction_coefficients import ArrheniusCoefficient, \
    ConstantCoefficient, ModifiedArrheniusCoefficient, RxnCoefficientBase


# tests for RxnCoefficientBase() base class
def test_base_get_coef ():
    try:
        RxnCoefficientBase().get_coef()
    except NotImplementedError as err:
        assert (type(err) == NotImplementedError)


def test_base_repr_ ():
    assert repr(RxnCoefficientBase()) == 'RxnCoefficientBase()'


# tests for ConstantCoefficient() class
def test_const_get_coef_result ():
    assert ConstantCoefficient(100).get_coef() == 100


def test_const_get_coef_neg_k ():
    try:
        ConstantCoefficient(-10).get_coef()
    except ValueError as err:
        assert (type(err) == ValueError)


def test_const_repr_ ():
    assert repr(ConstantCoefficient(5.0)) == 'ConstantCoefficient(k = 5.0)'


# tests for ArrheniusCoefficient() class
def test_arr_get_coef_result ():
    assert ArrheniusCoefficient(A=np.power(10, 4), E=100, T=300,
                                R=8.314).get_coef() == 9607.0007471344579


def test_arr_get_coef_neg_T ():
    try:
        ArrheniusCoefficient(A=np.power(10, 4), E=100, T=-300,
                             R=8.314).get_coef()
    except ValueError as err:
        assert (type(err) == ValueError)


def test_arr_get_coef_neg_A ():
    try:
        ArrheniusCoefficient(A=-np.power(10, 4), E=100, T=300,
                             R=8.314).get_coef()
    except ValueError as err:
        assert (type(err) == ValueError)


def test_arr_get_coef_neg_R ():
    try:
        ArrheniusCoefficient(A=np.power(10, 4), E=100, T=300,
                             R=-8.314).get_coef()
    except ValueError as err:
        assert (type(err) == ValueError)


def test_arr_get_coef_overflow ():
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        try:
            ArrheniusCoefficient(A=np.power(1, 256), E=-10000000., T=300,
                                 R=8.314).get_coef()
        except OverflowError as err:
            assert (str(err) == "The result is too large/small.")


def test_arr_repr_ ():
    assert repr(ArrheniusCoefficient(10, 100,
                                     150)) == 'ArrheniusCoefficient(A=10, ' \
                                              'E=100, T=150, R=8.314)'


# tests for ModifiedArrheniusCoefficient() class
def test_modarr_get_coef_result ():
    assert ModifiedArrheniusCoefficient(A=np.power(10, 4), b=0.5, E=100, T=300,
                                        R=8.314).get_coef() == \
           166398.13402389048


def test_modarr_get_coef_neg_T ():
    try:
        ModifiedArrheniusCoefficient(A=np.power(10, 4), b=0.5, E=100, T=-300,
                                     R=8.314).get_coef()
    except ValueError as err:
        assert (type(err) == ValueError)


def test_modarr_get_coef_neg_A ():
    try:
        ModifiedArrheniusCoefficient(A=-np.power(10, 4), b=0.5, E=100, T=300,
                                     R=8.314).get_coef()
    except ValueError as err:
        assert (type(err) == ValueError)


def test_modarr_get_coef_neg_R ():
    try:
        ModifiedArrheniusCoefficient(A=np.power(10, 4), b=0.5, E=100, T=300,
                                     R=-8.314).get_coef()
    except ValueError as err:
        assert (type(err) == ValueError)


def test_modarr_get_coef_complex_b ():
    try:
        ModifiedArrheniusCoefficient(A=np.power(10, 4), b=(3 + 1j), E=100,
                                     T=300, R=8.314).get_coef()
    except ValueError as err:
        assert (type(err) == ValueError)


def test_modarr_get_coef_overflow ():
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        try:
            ModifiedArrheniusCoefficient(A=np.power(10, 4), b=256.0, E=100,
                                         T=300, R=8.314).get_coef()
        except OverflowError as err:
            assert (str(err) == "The result is too large/small.")


def test_modarr_repr_ ():
    assert repr(ModifiedArrheniusCoefficient(A=10, b=0.5, E=100,
                                             T=300)) == 'ModifiedArrheniusCoefficient(A=10, b=0.5, E=100, T=300, R=8.314)'
