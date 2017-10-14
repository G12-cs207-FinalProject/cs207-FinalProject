#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 14 16:17:28 2017

@author: filipmichalsky
"""
import rxns as rx

def test_reaction_object():
    reac1 = rx.IrrevElemRxn([10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    return reac1
    
def test_progress_rate_result():
    assert (reac1.progress_rate) == [40.0, 10.0]).all()

def test_progress_rate_multi_neg_ki():
    try:
        progress_rate_multi(ki=[-10, 20], xi=[3.0, 2.0, 1.0], vi_p=[[1.0, 2.0, 0.0], [1.0, 0.0, 2.0]], vi_dp=[[0.0, 0.0, 2.0], [0.0, 1.0, 2.0]])
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_progress_rate_multi_neg_xi():
    try: 
        progress_rate_multi(ki=[10, 20], xi=[-3.0, 2.0, 1.0], vi_p=[[1.0, 2.0, 0.0], [1.0, 0.0, 2.0]], vi_dp=[[0.0, 0.0, 2.0], [0.0, 1.0, 2.0]])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_rate_multi_dim_compat():
    try:
        progress_rate_multi(ki=[10, 20], xi=[1.0, 2.0, 1.0], vi_p=[[1.0, 0.0], [1.0, 0.0, 2.0]], vi_dp=[[0.0, 0.0, 2.0], [0.0, 1.0, 2.0]])
    except ValueError as err:
        assert(type(err) == ValueError)
        