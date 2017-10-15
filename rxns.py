#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 14 15:23:22 2017

@author: filipmichalsky
"""
import numpy as np

class Rxn():
    """Base class for reaction"""
    
    def __init__(self):
        self.ki = None
        self.xi = None
        self.vi_p = None
        self.vi_dp = None
        self.wi = None
        self.rates = None
        self.name = None


    def progress_rate(self):
        """Not-implemented method to return the rate coefficient at base class"""
        raise NotImplementedError('Subclass must implement this method')

    def reaction_rate(self):
        """Not-implemented method to return the rate coefficient at base class"""
        raise NotImplementedError('Subclass must implement this method')
        
        
class ElemRxn(Rxn):
    """Subclass of Rxn for Elementary Reaction """
    pass

class NonElemRxn(Rxn):
    """Subclass of Rxn for NonElementary Reaction """
    pass

class RevElemRxn(ElemRxn):
    """Subclass of Rxn for Reversible Elementary Reaction """
    pass

class IrrevElemRxn(ElemRxn):
    """Subclass of Rxn for Irriversible Elementary Reaction """
    
    def __init__(self,ki,xi,vi_p,vi_dp):
        self.ki = ki
        self.xi = xi
        self.vi_p = vi_p
        self.vi_dp = vi_dp
        self.wi = None
        self.rates = None
        self.name = None

    
    def progress_rate(self):
        """
        Returns the progress rate w for a system of irreversible elementary reaction of the form: 
                v_11' A + v_21' B -> v_31" C
                v_12' A + v_32' C -> v_22" B + v_32"C
        
        INPUTS
        =======
        ki: a list of floats, required
            Reaction rate coefficients, units = (1/time) · mol^(1-q) · volume^(q-1) where q is the order of reaction
        xi: a list of floats, required,
            Concentrations of moleuclar species, units = mol/volume
        vi_p: a list of floats, optional, default value is [[1.0, 1.0, 0.0], [1.0, 0.0, 1.0]]
            Stoichiometric coefficients of the reactants
        vi_dp: a numpy matrix of floats, optional, default value is [[0.0, 0.0, 1.0], [0.0, 1.0, 1.0]]
            Stoichiometric coefficients of the products
        
        RETURNS
        ========
        w: a numpy array of floats,
           Has the form w unless k <= 0 or xi < 0
           in which cases a ValueError exception is raised
        
        NOTES
        =====
        PRE: 
             - ki, xi, vi_p, and vi_dp have list or np.array type
             - xi is ordered in the form [[A], [B], [C]]
             - vi_p is ordered in the form [[v_11', v_21', v_31'], [v_12', v_22', v_32']]
             - vi_dp is ordered in the form [[v_11", v_21", v_31"], [v_12", v_22", v_32"]]
             - four or fewer inputs
        POST:
             - ki, xi, vi_p, and vi_dp are not changed by this function
             - raises a ValueError exception if any(ki <= 0)
             - raises a ValueError exception if any(xi < 0)
             - returns a numpy array of floats of progress rate w
        
        EXAMPLES
        =========
        >>> IrrevElemRxn([10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]]).progress_rate()
        array([ 40.,  10.])
        """
        
        
        # check value conditions
        if (any(i <= 0 for i in self.ki)): # check reaction coefficients 
            raise ValueError("reaction rate coefficients ki must be positive.")
        elif (any(i < 0 for i in self.xi)): # check concentration array
            raise ValueError("concentrations xi cannot be negative.")
        else:
            # Convert inputs into numpy column vectors/matrices
            xi = np.array(self.xi).reshape(-1,1)
            vi_p = np.matrix(self.vi_p).T
            vi_dp = np.matrix(self.vi_dp).T
            ki = np.array(self.ki).reshape(-1,1)
        
            prodi = np.prod(np.power(xi, vi_p), axis=0) # calculate the product of xi^(vi_p) for each reaction
            self.wi = np.squeeze(ki* np.array(prodi).reshape(-1,1)) # multiply ki by prodi


            return self.wi
        
    def reaction_rate(self):
        """Returns the progress rate w for a system of irreversible elementary reaction of the form: 
                v_11' A + v_21' B -> v_31" C
                v_32' C -> v_12' A + v_22' B
        
        INPUTS
        =======
        ki: a list of floats, required
            Reaction rate coefficients, units = (1/time) · mol^(1-q) · volume^(q-1) where q is the order of reaction
        xi: a list of floats, required,
            Concentrations of moleuclar species, units = mol/volume
        vi_p: a list of floats, optional, default value is [[1.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
            Stoichiometric coefficients of the reactants
        vi_dp: a numpy matrix of floats, optional, default value is [[0.0, 0.0, 1.0], [1.0, 1.0, 0.0]]
            Stoichiometric coefficients of the products
        
        RETURNS
        ========
        rates: a numpy array of floats,
           Has the form rates unless k <= 0 or xi < 0
           in which cases a ValueError exception is raised
    
        NOTES
        =====
        PRE: 
             - ki, xi, vi_p, and vi_dp have list or np.array type
             - xi is ordered in the form [[A], [B], [C]]
             - vi_p is ordered in the form [[v_11', v_21', v_31'], [v_12', v_22', v_32']]
             - vi_dp is ordered in the form [[v_11", v_21", v_31"], [v_12", v_22", v_32"]]
             - four or fewer inputs
        POST:
             - ki, xi, vi_p, and vi_dp are not changed by this function
             - raises a ValueError exception if any(ki <= 0)
             - raises a ValueError exception if any(xi < 0)
             - returns a numpy array of floats of reaction rates, rates
        
        EXAMPLES
        =========
        >>> IrrevElemRxn([10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]]).reaction_rate()
        array([-60., -70.,  70.])
        """
        import numpy as np
        # check value conditions
        if (any(i <= 0 for i in self.ki)): # check reaction coefficients 
            raise ValueError("reaction rates ki must be positive.")
        elif (any(i < 0 for i in self.xi)): # check concentration array
            raise ValueError("concentrations xi cannot be negative.")
        else:
            # Convert inputs into numpy column vectors/matrices
            xi = np.array(self.xi).reshape(-1,1)
            vi_p = np.matrix(self.vi_p).T
            vi_dp = np.matrix(self.vi_dp).T
            ki = np.array(self.ki).reshape(-1,1)
        
            prodi = np.prod(np.power(xi, vi_p), axis=0) # calculate the product of xi^(vi_p) for each reaction
            w=self.progress_rate()
            #w = np.squeeze(ki* np.array(prodi).reshape(-1,1)).reshape(-1,1) # multiply ki by prodi
            vi = vi_dp - vi_p 
            self.rates = np.squeeze(np.array(np.dot(vi,w)))
            return self.rates

#reac1 = IrrevElemRxn([10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])
#print(reac1.progress_rate())
#print(reac1.reaction_rate())
#print(repr(reac1))
#progress_rate_multi(ki=[10, 10], xi=[1.0, 2.0, 1.0], vi_p=[[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]], vi_dp=[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]])

