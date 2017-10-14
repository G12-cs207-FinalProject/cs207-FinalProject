import numpy as np
from reaction_coefs import *

def test_base_compute_coef():
    try:
        RxnCoef().compute_coef()
    except NotImplementedError as err:
        assert(type(err) == NotImplementedError)

def test_const_compute_coef():
    assert ConstCoef(100).compute_coef() == 100
