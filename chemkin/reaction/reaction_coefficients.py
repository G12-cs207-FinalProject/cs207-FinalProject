import numpy as np
from chemkin.thermodynamics.thermo import ThermoDAO
from chemkin.chemkin_errors import ChemKinError

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
        """Calculates and returns the reaction rate coefficient

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

class BackwardCoefficient():
    """ Class of BackwardCoefficient
    """

    def __init__ (self, species, T, ki, is_reversible, vi_p, vi_dp, db_name='NASA_coef.sqlite'):
        self.p0 = 10e5
        self.R = 8.314
        self.species = species
        self.T = T
        self.ki = ki
        self.b_ki = np.zeros((len(is_reversible,)))
        self.is_reversible = is_reversible
        self.vi_p = vi_p
        self.vi_dp = vi_dp
        self.vi = np.matrix(self.vi_dp).T - np.matrix(self.vi_p).T  # calculate overall stoicheometric coefficients
        self.gamma = np.sum(self.vi, axis=0)
        self.dao = ThermoDAO(db_name)


    def get_backward_coefs (self):
        species_high = self.dao.get_species(self.T, 'high')
        species_low = self.dao.get_species(self.T, 'low')

        H_over_RT_species = []
        S_over_R_species = []
        for s in self.species:
            if self.T >= 1000: # high temperature range
                if s not in species_high:
                    raise ChemKinError('Thermo.get_backward_coefs()',
                        'The specie {}\'s high temperature bound is not higher than the current T={}.'.format(s, self.T))
                NASA_coefs = self.dao.get_coeffs(s, 'high')
            else: # low temperature range
                if s not in species_low:
                    raise ChemKinError('Thermo.get_backward_coefs()',
                        'The specie {}\'s low temperature bound is not lower than the current T={}.'.format(s, self.T))
                NASA_coefs = self.dao.get_coeffs(s, 'low')

            # Compute enthalpy/(RT) of a specie
            H_T_arr = np.array([1, 1/2*self.T, 1/3*(self.T**2), 1/4*(self.T**3), 1/5*(self.T**4), 1/self.T])
            H_coef_arr = np.array([NASA_coefs[0], NASA_coefs[1], NASA_coefs[2], NASA_coefs[3], NASA_coefs[4], NASA_coefs[5]])
            H_over_RT_species.append(np.dot(H_T_arr, H_coef_arr))

            # Compute entropy/R of a specie
            S_T_arr = np.array([np.log(self.T), self.T, 1/2*(self.T**2), 1/3*(self.T**3), 1/4*(self.T**4), 1])
            S_coef_arr = np.array([NASA_coefs[0], NASA_coefs[1], NASA_coefs[2], NASA_coefs[3], NASA_coefs[4], NASA_coefs[6]])
            S_over_R_species.append(np.dot(S_T_arr, S_coef_arr))

        delta_H_over_RT = np.squeeze(np.asarray(np.dot(H_over_RT_species, self.vi))) # delta enthalpy of a system of reactions
        delta_S_over_R = np.squeeze(np.asarray(np.dot(S_over_R_species, self.vi))) # delta entropy of a system of reactions

        factor = np.asarray(np.power(self.p0 / (self.R * self.T), self.gamma))[0]

        for index, rxn_is_rev in enumerate(self.is_reversible):
            if rxn_is_rev:
                ke = factor[index] * np.exp(delta_S_over_R[index] - delta_H_over_RT[index])

                self.b_ki[index] = self.ki[index]/ke

        return self.b_ki
        
if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)