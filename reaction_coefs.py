import numpy as np
import doctest

class RxnCoef():
	"""Base class for reaction coefficients"""

	def __init__(self):
		self.k = None


	def get_coef(self):
		"""Not-implemented method to return the rate coefficient at base class"""
		raise NotImplementedError('Subclass must implement this method')


class ConstCoef(RxnCoef):
	"""Subclass of RxnCoef for constant reaction coefficients"""
	def __init__(self, k):
		self.k = k

	def get_coef(self):
		"""Simply returns the constant rate coefficient

		RETURNS:
		========
		k: float
			Constant reaction rate coefficient
		
		EXAMPLES:
		=========
		>>> ConstCoef(10.0).get_coef()
		10.0
		"""
		if self.k < 0:
			raise ValueError("Negative reaction rate coefficients are prohibited.")
		return self.k

class ArrheniusCoef(RxnCoef):
	"""Subclass of RxnCeof for Arrhenius reaction coefficients """
	def __init__(self, A, E, T, R=8.314):
		self.A = A
		self.E = E
		self.T = T
		self.R = R
		super().__init__()


	def get_coef(self):
		""" Method to calculate the Arrhenius reaction rate coefficient

		RETURNS:
		========
		k: float
			Arrhenius reaction rate coefficient

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
	def __init__(self, A, b, E, T, R=8.314):
		super().__init__(A, E, T, R)
		self.b = b

	def get_coef(self):
		""" Method to calculate the Modified Arrhenius reaction rate coefficient

		RETURNS:
		========
		k: float
		Modified Arrhenius reaction rate coefficient

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

