#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 14 16:17:28 2017

@author: filipmichalsky
"""
import rxns as rx


# tests for Rxn() base class
def test_base_progress_rate():
    try:
        rx.Rxn().progress_rate()
    except NotImplementedError as err:
        assert(type(err) == NotImplementedError)

def test_base_reaction_rate():
    try:
        rx.Rxn().reaction_rate()
    except NotImplementedError as err:
        assert(type(err) == NotImplementedError)




# tests for IrrevElemRxn() class    
def test_progress_rate_result_IrrevElem():
    reac1 = rx.IrrevElemRxn([10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    assert (reac1.progress_rate() == [40.0, 10.0]).all()


def test_progress_rate_neg_ki_IrrevElem():
    reac1 = rx.IrrevElemRxn([-10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    try:
        reac1.progress_rate()
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_progress_rate_neg_xi_IrrevElem():
    reac1 = rx.IrrevElemRxn([10, 10], [-1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    try: 
        reac1.progress_rate()
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_rate_dim_compat_IrrevElem():
    reac1 = rx.IrrevElemRxn([10, 10], [1.0, 2.0, 1.0], [[ 2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    try:
        reac1.progress_rate()
    except ValueError as err:
        assert(type(err) == ValueError)

    

def test_reaction_rate_result_IrrevElem():
    reac1 = rx.IrrevElemRxn([10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])

    assert (reac1.reaction_rate() == [-60., -70., 70.]).all()

def test_reaction_rate_neg_ki_IrrevElem():
    reac1 = rx.IrrevElemRxn([-10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    try:
        reac1.reaction_rate()
    except ValueError as err:
        assert(type(err) == ValueError)
        

def test_reaction_rate_neg_xi_IrrevElem():
    reac1 = rx.IrrevElemRxn([10, 10], [-1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    try:
         reac1.reaction_rate()
    except ValueError as err:
        assert(type(err) == ValueError)


def test_reaction_rate_multi_dim_compat_IrrevElem():
    reac1 = rx.IrrevElemRxn([10, 10], [1.0, 2.0, 1.0], [[ 2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
    try:
        reac1.reaction_rate()
    except ValueError as err:
        assert(type(err) == ValueError)

