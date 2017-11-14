"""
Test suite for the elementary_rxn.py module

"""

import chemkin.reaction.elementary_rxn as er
import io
import sys


# tests for IrrevElemRxn() class
def test_IrrevElem_len_result():
    reac1 = er.IrreversibleElementaryRxn([10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]], [[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    assert len(reac1) == 2


def test_IrrevElem_repr_result():
    old_stdout = sys.stdout
    capturedOutput = io.StringIO()
    sys.stdout = capturedOutput
    reac1 = er.IrreversibleElementaryRxn(ki=[10, 10], xi=[1.0, 2.0, 1.0], vi_p=[[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]], vi_dp=[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])

    print(reac1)
    sys.stdout = old_stdout
    assert capturedOutput.getvalue() == "IrreversibleElementaryRxn(ki=[10, 10], xi=[1.0, 2.0, 1.0], vi_p=[[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]], vi_dp=[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])\n"


def test_IrrevElem_progress_rate_result():
    reac1 = er.IrreversibleElementaryRxn([10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]],
                         [[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    assert (reac1.progress_rate() == [40.0, 10.0]).all()


def test_IrrevElem_progress_rate_neg_ki():
    reac1 = er.IrreversibleElementaryRxn([-10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]]
                         ,[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    try:
        reac1.progress_rate()
    except ValueError as err: assert(type(err) == ValueError)


def test_IrrevElem_progress_rate_neg_xi():
    reac1 = er.IrreversibleElementaryRxn([10, 10], [-1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    try:
        reac1.progress_rate()
    except ValueError as err:
        assert(type(err ) == ValueError)


def test_IrrevElem_progress_rate_dim_compat():
    reac1 = er.IrreversibleElementaryRxn([10, 10], [1.0, 2.0, 1.0], [[ 2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    try:
        reac1.progress_rate()
    except ValueError as err:
        assert(type(err ) == ValueError)


def test_IrrevElem_reaction_rate_result():
    reac1 = er.IrreversibleElementaryRxn([10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0,
                         0.0, 2.0], [0.0, 1.0, 1.0]])

    assert (reac1.reaction_rate() == [-60., -70., 70.]).all()


def test_IrrevElem_reaction_rate_neg_ki():
    reac1 = er.IrreversibleElementaryRxn([-10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0
                         , 2.0], [0.0, 1.0, 1.0]])
    try:
        reac1.reaction_rate()
    except ValueError as err:
        assert(type(err) == ValueError)


def test_IrrevElem_reaction_rate_neg_xi():
    reac1 = er.IrreversibleElementaryRxn([10, 10], [-1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0
                         , 2.0], [0.0, 1.0, 1.0]])
    try:
        reac1.reaction_rate()
    except ValueError as err:
        assert(type(err) == ValueError)


def test_IrrevElem_reaction_rate_dim_compat():
    reac1 = er.IrreversibleElementaryRxn([10, 10], [1.0, 2.0, 1.0], [[2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0 , 2.0], [0.0, 1.0, 1.0]])
    try:
        reac1.reaction_rate()
    except ValueError as err:
        assert(type(err) == ValueError)