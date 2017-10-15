import numpy as np

class Rxn():
    """Base class of reactions

    ATTRIBUTES:
    ========
    self.ki: float or a list of floats, required
        Reaction rate coefficient
        Initialized with constructor argument ki

    self.xi: a list of floats, required,
        Concentrations of molecular species
        Initialized with constructor argument xi

    self.vi_p: a list of floats, required,
        Stoichiometric coefficients of the reactants
        Initialized with constructor argument vi_p

    self.vi_dp: a list of floats, required,
        Stoichiometric coefficients of the products
        Initialized with constructor argument vi_dp

    METHODS:
    ========
    __len__(): Returns the number of reactions
    __repr__(): Prints the class name with its attributes

    progress_rate(): Calculates and returns the progress rate
        This should be implemented in its subclass.

    reaction_rate(): Calculates and returns the reaction rate
        This should be implemented in its subclass.
    """
    
    def __init__(self, ki, xi, vi_p, vi_dp):
        self.ki = ki
        self.xi = xi
        self.vi_p = vi_p
        self.vi_dp = vi_dp
        self.wi = None
        self.rates = None

    # def __len__(self):
    #     """Returns the number of species in the reaction"""
    #     return len(self.xi)

    def __len__(self):
        """Returns the number of reactions"""
        return len(self.vi_p)

    def __repr__(self):
        return 'Rxn(ki={}, xi={}, vi_p={}, vi_dp={})'.format(self.ki, self.xi, self.vi_p,  self.vi_dp)

    def progress_rate(self):
        """ Calculates and returns the progress rate

        Method is not implemented in this base class, but should be implemented in its subclasses.
        """
        raise NotImplementedError('Subclass must implement this method')

    def reaction_rate(self):
        """ Calculates and returns the reaction rate

        Method is not implemented in this base class, but should be implemented in its subclasses.
        """
        raise NotImplementedError('Subclass must implement this method')
        
        
class ElemRxn(Rxn):
    """Class of elementary reactions
    Subclass of Rxn
    """
    pass

class NonElemRxn(Rxn):
    """Class of non-elementary reactions
    Subclass of Rxn
    """
    pass

class RevElemRxn(ElemRxn):
    """Class of reversible elementary reactions
    Subclass of ElemRxn
    """
    pass

class IrrevElemRxn(ElemRxn):
    """Class of irreversible elementary reactions
    Subclass of ElemRxn
    Default form of
            v_11' A + v_21' B -> v_31" C
            v_32' C -> v_12' A + v_22' B

    ATTRIBUTES:
    ========
    self.ki: float or a list of floats, optional, default value = [10.0, 10.0]
        Reaction rate coefficient
        Initialized with constructor argument ki

    self.xi: a list of floats, optional, default value = [1.0, 1.0, 1.0]
        Concentrations of molecular species
        Initialized with constructor argument xi

    self.vi_p: a list of floats, optional, default value = [[1.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        Stoichiometric coefficients of the reactants
        Initialized with constructor argument vi_p

    self.vi_dp: a list of floats, optional, default value = [[0.0, 0.0, 1.0], [1.0, 1.0, 0.0]]
        Stoichiometric coefficients of the products
        Initialized with constructor argument vi_dp

    NOTES
    =====
    PRE:
        - ki, xi, vi_p, and vi_dp have list or np.array type
        - xi is ordered in the form [[A], [B], [C]]
        - vi_p is ordered in the form [[v_11', v_21', v_31'], [v_12', v_22', v_32']]
        - vi_dp is ordered in the form [[v_11", v_21", v_31"], [v_12", v_22", v_32"]]
        - four or fewer inputs
    METHODS:
    ========
    __len__(): Returns the number of reactions
    __repr__(): Prints the class name with its attributes
    progress_rate(): Calculates and returns the progress rate
    reaction_rate(): Calculates and returns the reaction rate
    """
    
    def __init__(self, ki=[10.0, 10.0], xi=[1.0, 1.0, 1.0], vi_p= [[1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], vi_dp = [[0.0, 0.0, 1.0], [1.0, 1.0, 0.0]]):
        self.ki = ki
        self.xi = xi
        self.vi_p = vi_p
        self.vi_dp = vi_dp
        self.wi = None
        self.rates = None


    # def __len__(self):
    #     """Returns the number of species in the reaction"""
    #     return len(self.xi)

    def __len__(self):
        """Returns the number of reactions"""
        return len(self.vi_p)

    def __repr__(self):
        return 'IrrevElemRxn(ki={}, xi={}, vi_p={}, vi_dp={})'.format(self.ki, self.xi, self.vi_p,  self.vi_dp)
    
    def progress_rate(self):
        """
        Returns the progress rate w for a system of irreversible elementary reaction
        
        RETURNS
        ========
        self.wi: a numpy array of floats,
           Has the form self.w unless self.ki <= 0 or self.xi < 0
           in which cases a ValueError exception is raised
        
        NOTES
        =====
        POST:
             - self.ki, self.xi, self.vi_p, and self.vi_dp are not changed by this function
             - raises a ValueError exception if any(self.ki <= 0)
             - raises a ValueError exception if any(self.xi < 0)
             - returns a numpy array of floats of progress rate self.wi
        
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
            ki = np.array(self.ki).reshape(-1,1)
        
            prodi = np.prod(np.power(xi, vi_p), axis=0) # calculate the product of xi^(vi_p) for each reaction
            self.wi = np.squeeze(ki* np.array(prodi).reshape(-1,1)) # multiply ki by prodi


            return self.wi
        
    def reaction_rate(self):
        """Returns the progress rate w for a system of irreversible elementary reaction

        RETURNS
        ========
        self.rates: a numpy array of floats,
           Has the form self.rates unless self.ki <= 0 or self.xi < 0
           in which cases a ValueError exception is raised
    
        NOTES
        =====
        POST:
             - self.ki, self.xi, self.vi_p, and self.vi_dp are not changed by this function
             - raises a ValueError exception if any(self.ki <= 0)
             - raises a ValueError exception if any(self.xi < 0)
             - returns a numpy array of floats of reaction rates, self.rates
        
        EXAMPLES
        =========
        >>> IrrevElemRxn([10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]]).reaction_rate()
        array([-60., -70.,  70.])
        """
        import numpy as np
        # check value conditions
        if any(i <= 0 for i in self.ki): # check reaction coefficients
            raise ValueError("reaction rates ki must be positive.")
        elif any(i < 0 for i in self.xi): # check concentration array
            raise ValueError("concentrations xi cannot be negative.")
        else:
            # Convert inputs into numpy column vectors/matrices
            vi_p = np.matrix(self.vi_p).T
            vi_dp = np.matrix(self.vi_dp).T

            w = self.progress_rate()  # calculate progress rate
            vi = vi_dp - vi_p  # calculate overall stoicheometric coefficients

            self.rates = np.squeeze(np.array(np.dot(vi, w)))  # calculate reaction rate

            return self.rates
