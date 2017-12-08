"""
Test suite for the base_rxn.py module

"""
import chemkin.reaction.base_rxn as rxn
import sys
import io

# tests for RxnBase() base class
def test_RxnBase_progress_rate_not_implemented():
    try:
        reac1 = rxn.RxnBase([10, 10], [20, 20], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]], [[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
        reac1.progress_rate()
    except NotImplementedError as err:
        assert(type(err) == NotImplementedError)


def test_RxnBase_reaction_rate_not_implemented():
    try:
        reac1 = rxn.RxnBase([10, 10], [20, 20], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]], [[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
        reac1.progress_rate()
    except NotImplementedError as err:
        assert(type(err) == NotImplementedError)


def test_RxnBase_len_result():
    reac1 = rxn.RxnBase([10, 10], [20, 20], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]], [[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    assert len(reac1) == 2


def test_RxnBase_repr_result():
    old_stdout = sys.stdout
    capturedOutput = io.StringIO()
    sys.stdout = capturedOutput
    reac1 = rxn.RxnBase([10, 10], [20, 20], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]], [[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    print(reac1)
    sys.stdout = old_stdout
    assert capturedOutput.getvalue() == "RxnBase(ki=[10, 10], b_ki=[20, 20], xi=[1.0, 2.0, 1.0], vi_p=[[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]], vi_dp=[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])\n"

def test_RxnBase_reaction_rate():
    try:
        reac1 = rxn.RxnBase([10, 10], [20, 20], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]], [[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
        reac1.reaction_rate()
    except NotImplementedError as err:
        assert(type(err) == NotImplementedError)