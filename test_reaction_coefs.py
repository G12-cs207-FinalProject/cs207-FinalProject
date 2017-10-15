import numpy as np
import warnings
from reaction_coefs import *

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