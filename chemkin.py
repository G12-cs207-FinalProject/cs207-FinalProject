"""
Chemical Kinetics
-----------------

Module chemkin was created for calculating progress rates and reaction rates for a single reaction/set of reactions.
The module handles both reversible and irreversible reactions, as well as elementary/nonelementary reactions.

The module contains classes to calculate reaction coefficients and reaction base classes (RxnCoef and
Reaction) as well as subclasses of RxnCoeff class (ConstCoef,ArrheniusCoef,ModArrheniusCoef)
and subclasses of Rxn class (ElemRxn and NonElemRxn). Each of the latter also has subclasses 
(IrrElemRxn,RevElemRxn,IrrNonElemRxn, and RevNonElemRxn)


demo.py file serves as a demo of the code.

Created by Michelle Ho, Jasmine Tong, Filip Michalsky and Nathaniel Stein

#to be added - RevNonElemRxn and IrrNonElemRxn classes
"""


import numpy as np

class RxnCoef():
	"""Base class of reaction rate coefficients

	ATTRIBUTES:
	========
	self.k: float
		Reaction rate coefficient
		Initialized to None in base class and should be calculated and returned in its subclass.

	METHODS:
	========
	get_coef(): Calculates and returns the reaction rate coefficient
		This should be implemented in its subclass.
	"""

	def __init__(self):
		self.k = None


	def __repr__(self):
		return 'RxnCoef()'

	def get_coef(self):
		""" Calculates and returns the reaction rate coefficient

		Method is not implemented in this base class, but should be implemented in its subclasses.
		"""
		raise NotImplementedError('Subclass must implement this method')


class ConstCoef(RxnCoef):
	"""Class of constant reaction rate coefficients
	Subclass of RxnCoef

	ATTRIBUTES:
	========
	self.k: float, required
		Reaction rate coefficient
		Initialized with constructor argument k

	METHODS:
	========
	get_coef(): Returns the constant rate coefficient
	"""
	def __init__(self, k):

		"""
		NOTES
		=====
		PRE:
			- self.k have numeric type
			- one input
		"""
		self.k = k

	def __repr__(self):
		return 'ConstCoef(k = {})'.format(self.k)

	def get_coef(self):
		"""Returns the constant rate coefficient

		RETURNS:
		========
		self.k: float
			Constant reaction rate coefficient
		NOTES
		=====
		POST:
			 - self.k is not changed by this function
			 - raises a ValueError exception if k <= 0
			 - returns a float of self.k
		
		EXAMPLES:
		=========
		>>> ConstCoef(10.0).get_coef()
		10.0
		"""
		if self.k < 0:
			raise ValueError("Negative reaction rate coefficients are prohibited.")
		return self.k

class ArrheniusCoef(RxnCoef):
	"""Class of Arrhenius reaction rate coefficient
	Subclass of RxnCoef

	ATTRIBUTES:
	========
	self.k: float
		Reaction rate coefficient
		Calculated by self.get_coef()
	self.A: float, required
		Arrhenius prefactor
		Initialized with constructor argument A
	self.E: float, required
		Activation energy
		Initialized with constructor argument E
	self.T: float, required
		Temperature (in Kelvin)
		Initialized with constructor argument T
	self.R: float, optional, default value = 8.314
		Ideal gas constant
		Initialized with constructor argument R
		Default value of R should not be changed except for unit conversion

	METHODS:
	========
	get_coef(): Calculates and returns the Arrhenius rate coefficient
	"""
	def __init__(self, A, E, T, R=8.314):

		"""
		NOTES
		=====
		PRE:
			- self.A, self.E, self.T and self.R have numeric type
			- four or fewer inputs
		"""
		super().__init__()
		self.A = A
		self.E = E
		self.T = T
		self.R = R

	def __repr__(self):
		return 'ArrheniusCoef(A={}, E={}, T={}, R={})'.format(self.A, self.E, self.T, self.R)

	def get_coef(self):
		""" Method to calculate the Arrhenius reaction rate coefficient

		RETURNS:
		========
		self.k: float
			Arrhenius reaction rate coefficient

		NOTES
		=====
		POST:
			 - self.A, self.E, self.T and self.R are not changed by this function
			 - raises a ValueError exception if A <= 0 or T <= 0
			 - returns a float of self.k

		EXAMPLES:
		=========
		>>> ArrheniusCoef(2.0, 3.0, 100.0).get_coef()
		1.9927962618542914
		"""
		if self.A < 0.0:
			raise ValueError("A = {0:18.16e}:  Negative Arrhenius prefactor is prohibited!".format(self.A))

		if self.T < 0.0:
			raise ValueError("T = {0:18.16e}:  Negative temperatures are prohibited!".format(self.T))

		if self.R < 0.0:
			raise ValueError("R = {0:18.16e}:  Negative ideal gas constant is prohibited!".format(self.R))

		else:
			try:
				self.k = self.A * np.exp(-self.E / (self.R * self.T))
			except Warning:
				raise OverflowError("The result is too large/small.")
			return self.k


class ModArrheniusCoef(ArrheniusCoef):
	"""Class of Modified Arrhenius reaction rate coefficient
		Subclass of ArrheniusCoef

		ATTRIBUTES:
		========
		self.k: float
			Reaction rate coefficient
			Calculated by self.get_coef()
		self.A: float, required
			Arrhenius prefactor
			Initialized with constructor argument A
		self.b: float, required
			Modified Arrhenius parameter
			Initialized with constructor argument b
		self.E: float, required
			Activation energy
			Initialized with constructor argument E
		self.T: float, required
			Temperature (in Kelvin)
			Initialized with constructor argument T
		self.R: float, optional, default value = 8.314
			Ideal gas constant
			Initialized with constructor argument R
			Default value of R should not be changed except for unit conversion

		METHODS:
		========
		get_coef(): Calculates and returns the Modified Arrhenius rate coefficient
	"""
	def __init__(self, A, b, E, T, R=8.314):

		"""
		NOTES
		=====
		PRE:
			- self.A, self.b, self.E, self.T and self.R have numeric type
			- five or fewer inputs
		"""
		super().__init__(A, E, T, R)
		self.b = b

	def __repr__(self):
		return 'ModArrheniusCoef(A={}, b={}, E={}, T={}, R={})'.format(self.A, self.b, self.E, self.T, self.R)

	def get_coef(self):
		""" Method to calculate the Modified Arrhenius reaction rate coefficient

		RETURNS:
		========
		self.k: float
		Modified Arrhenius reaction rate coefficient

		NOTES
		=====
		POST:
			 - self.A, self.E, self.T and self.R are not changed by this function
			 - raises a ValueError exception if A <= 0, T <= 0 or b is not a real number
			 - returns a float of self.k

		EXAMPLES:
		=========
		>>> ModArrheniusCoef(2.0, -0.5, 3.0, 100.0).get_coef()
		0.19927962618542916
		"""
		if self.A < 0.0:
			raise ValueError("A = {0:18.16e}:  Negative Arrhenius prefactor is prohibited!".format(self.A))

		if self.T < 0.0:
			raise ValueError("T = {0:18.16e}:  Negative temperatures are prohibited!".format(self.T))

		if self.R < 0.0:
			raise ValueError("R = {0:18.16e}:  Negative ideal gas constant is prohibited!".format(self.R))

		if not np.isreal(self.b):
			raise ValueError('Modified Arrhenius parameter b must be real!')

		else:
			try:
				self.k = self.A * np.power(self.T, self.b) * np.exp(-self.E / (self.R * self.T))
			except Warning:
				raise OverflowError("The result is too large/small.")

			return self.k

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
