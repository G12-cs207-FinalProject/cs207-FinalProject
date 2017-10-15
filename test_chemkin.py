"""
Test suite for the chemkin.py module

"""
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