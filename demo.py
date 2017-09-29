"""
Demo
----

Temporary file containing funcs used to ensure Travis CI and Coveralls are set up
properly with repo.
"""


def add_num (x1, x2):
    return x1 + x2


def multiply (x1, x2):
    return x1 * x2


def test_add_num ():
    assert add_num(2, 8) == 10


def test_multiply ():
    assert multiply(2, 5) == 10
