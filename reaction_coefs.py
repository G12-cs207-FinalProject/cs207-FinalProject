import numpy as np
import doctest

class RxnCoef():
	"""Base class for reaction coefficients"""
	def __init__(self):
		self.coef = None

	def compute_coef(self):
		"""Not-implemented method to return the rate coefficient at base class"""
		raise NotImplementedError('Subclass must implement this method')


class ConstCoef(RxnCoef):
	"""Subclass for constant reaction coefficients"""
	def __init__(self, k):
		self.coef = k

	def compute_coef(self):
		"""Simply returns the constant rate coefficient

		RETURNS:
		========
		coef: float
			Constant reaction rate coefficient
		
		EXAMPLES:
		=========
		>>> test_coef_const.compute_coef()
		10.0
		"""
		return self.coef

#doctest.testmod(extraglobs={'test_coef_const': ConstCoef(10.0)}, verbose=True)