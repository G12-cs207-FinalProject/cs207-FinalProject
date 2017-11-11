import numpy as np

class RxnCoefficientBase():
    """Base class of reaction rate coefficients

	ATTRIBUTES:
	========
	k: float
		Reaction rate coefficient
		Initialized to None in base class and should be calculated and
		returned in its subclass.

	METHODS:
	========
	get_coef(): Calculates and returns the reaction rate coefficient
		This should be implemented in its subclass.
	"""

    def __init__ (self):
        self.k = None

    def __repr__ (self):
        return 'RxnCoefficientBase()'

    def __eq__ (self, other):
        """ Returns the true if 2 coefficients have the same value"""
        if self.get_coef() == other.get_coef():
            return True
        else:
            return False

    def get_coef (self):
        """ Calculates and returns the reaction rate coefficient

        Method is not implemented in this base class, but should be
        implemented in its subclasses.
        """
        raise NotImplementedError('Subclass must implement this method')


class ConstantCoefficient(RxnCoefficientBase):
    """Class of constant reaction rate coefficients
	Subclass of RxnCoefficientBase

	ATTRIBUTES:
	========
	k: float, required
		Reaction rate coefficient
		Initialized with constructor argument k

	METHODS:
	========
	get_coef(): Returns the constant rate coefficient
	"""

    def __init__ (self, k):
        """
        NOTES
        =====
        PRE:
            - k have numeric type
            - one input
        """
        self.k = k

    def __repr__ (self):
        return 'ConstantCoefficient(k = {})'.format(self.k)

    def get_coef (self):
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
        >>> ConstantCoefficient(10.0).get_coef()
        10.0
        """
        if self.k < 0:
            raise ValueError(
                  "Negative reaction rate coefficients are prohibited.")
        return self.k


class ArrheniusCoefficient(RxnCoefficientBase):
    """Class of Arrhenius reaction rate coefficient
	Subclass of RxnCoefficientBase

	ATTRIBUTES:
	========
	k: float
		Reaction rate coefficient
		Calculated by self.get_coef()
	A: float, required
		Arrhenius prefactor
		Initialized with constructor argument A
	E: float, required
		Activation energy
		Initialized with constructor argument E
	T: float, required
		Temperature (in Kelvin)
		Initialized with constructor argument T
	R: float, optional, default value = 8.314
		Ideal gas constant
		Initialized with constructor argument R
		Default value of R should not be changed except for unit conversion

	METHODS:
	========
	get_coef(): Calculates and returns the Arrhenius rate coefficient
	"""

    def __init__ (self, A, E, T, R=8.314):

        """
        NOTES
        =====
        PRE:
            - A, E, T and R have numeric type
            - four or fewer inputs
        """
        super().__init__()
        self.A = A
        self.E = E
        self.T = T
        self.R = R

    def __repr__ (self):
        return 'ArrheniusCoefficient(A={}, E={}, T={}, R={})'.format(self.A, self.E,
                                                              self.T, self.R)

    def get_coef (self):
        """ Method to calculate the Arrhenius reaction rate coefficient

        RETURNS:
        ========
        self.k: float
            Arrhenius reaction rate coefficient

        NOTES
        =====
        POST:
             - self.A, self.E, self.T and self.R are not changed by this
             function
             - raises a ValueError exception if A <= 0 or T <= 0
             - returns a float of self.k

        EXAMPLES:
        =========
        >>> ArrheniusCoefficient(2.0, 3.0, 100.0).get_coef()
        1.9927962618542914
        """
        if self.A < 0.0:
            raise ValueError(
                  "A = {0:18.16e}:  Negative Arrhenius prefactor is "
                  "prohibited!".format(
                        self.A))

        if self.T < 0.0:
            raise ValueError(
                  "T = {0:18.16e}:  Negative temperatures are "
                  "prohibited!".format(
                        self.T))

        if self.R < 0.0:
            raise ValueError(
                  "R = {0:18.16e}:  Negative ideal gas constant is "
                  "prohibited!".format(
                        self.R))

        else:
            try:
                self.k = self.A * np.exp(-self.E / (self.R * self.T))
            except Warning:
                raise OverflowError("The result is too large/small.")
            return self.k


class ModifiedArrheniusCoefficient(ArrheniusCoefficient):
    """Class of Modified Arrhenius reaction rate coefficient
		Subclass of ArrheniusCoefficient

		ATTRIBUTES:
		========
		k: float
			Reaction rate coefficient
			Calculated by self.get_coef()
		A: float, required
			Arrhenius prefactor
			Initialized with constructor argument A
		b: float, required
			Modified Arrhenius parameter
			Initialized with constructor argument b
		E: float, required
			Activation energy
			Initialized with constructor argument E
		T: float, required
			Temperature (in Kelvin)
			Initialized with constructor argument T
		R: float, optional, default value = 8.314
			Ideal gas constant
			Initialized with constructor argument R
			Default value of R should not be changed except for unit conversion

		METHODS:
		========
		get_coef(): Calculates and returns the Modified Arrhenius rate
		coefficient
	"""

    def __init__ (self, A, b, E, T, R=8.314):

        """
		NOTES
		=====
		PRE:
			- A, b, E, T and R have numeric type
			- five or fewer inputs
		"""
        super().__init__(A, E, T, R)
        self.b = b

    def __repr__ (self):
        return 'ModifiedArrheniusCoefficient(A={}, b={}, E={}, T={}, R={})'.format(self.A,
                                                                       self.b,
                                                                       self.E,
                                                                       self.T,
                                                                       self.R)

    def get_coef (self):
        """ Method to calculate the Modified Arrhenius reaction rate coefficient

        RETURNS:
        ========
        self.k: float
        Modified Arrhenius reaction rate coefficient

        NOTES
        =====
        POST:
             - self.A, self.E, self.T and self.R are not changed by this
             function
             - raises a ValueError exception if A <= 0, T <= 0 or b is not a
             real number
             - returns a float of self.k

        EXAMPLES:
        =========
        >>> ModifiedArrheniusCoefficient(2.0, -0.5, 3.0, 100.0).get_coef()
        0.19927962618542916
        """
        if self.A < 0.0:
            raise ValueError(
                  "A = {0:18.16e}:  Negative Arrhenius prefactor is "
                  "prohibited!".format(
                        self.A))

        if self.T < 0.0:
            raise ValueError(
                  "T = {0:18.16e}:  Negative temperatures are "
                  "prohibited!".format(
                        self.T))

        if self.R < 0.0:
            raise ValueError(
                  "R = {0:18.16e}:  Negative ideal gas constant is "
                  "prohibited!".format(
                        self.R))

        if not np.isreal(self.b):
            raise ValueError('Modified Arrhenius parameter b must be real!')

        else:
            try:
                self.k = self.A * np.power(self.T, self.b) * np.exp(
                      -self.E / (self.R * self.T))
            except Warning:
                raise OverflowError("The result is too large/small.")

            return self.k