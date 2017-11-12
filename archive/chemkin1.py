"""
Chemical Kinetics
-----------------

Module chemkin was created for calculating progress rates and reaction rates
for a single reaction/set of reactions.
The module handles both reversible and irreversible reactions, as well as
elementary/nonelementary reactions.

The module contains classes to calculate reaction coefficients and reaction
base classes (RxnCoef and Rxn) as well as subclasses of RxnCoeff class (
ConstCoef,ArrheniusCoef, ModArrheniusCoef) and subclasses of Rxn class (
ElemRxn and NonElemRxn). Each of the latter also has subclasses (IrrevElemRxn,
RevElemRxn, IrrevNonElemRxn, and RevNonElemRxn)

The module also contains 2 helper classes (XMLParser and RxnData) that allow
you to work with reaction data stored in XML files.


run_chemkin.py file serves as a demo of the code.

Created by Michelle Ho, Jasmine Tong, Filip Michalsky and Nathaniel Stein

#to be added - RevNonElemRxn, IrrNonElemRxn and RevNonElemRxn classes
"""

from enum import Enum
import xml.etree.ElementTree as ET

import numpy as np


class RxnType(Enum):
    Elementary = 1


class ChemKinError(Exception):
    """ Encapsulates errors encountered in this module. """

    def __init__ (self, method, info=None):
        msg = 'Error encountered in chemkin.py method: {0}.'.format(method)
        if info is not None:
            msg = msg + ' ' + info
        Exception.__init__(self, msg)

        self.method = method
        self.info = info


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
        return 'RxnCoef()'

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


class RxnBase():
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

    def __init__ (self, ki, xi, vi_p, vi_dp):
        self.ki = ki
        self.xi = xi
        self.vi_p = vi_p
        self.vi_dp = vi_dp
        self.wi = None
        self.rates = None

    # def __len__(self):
    #     """Returns the number of species in the reaction"""
    #     return len(self.xi)

    def __len__ (self):
        """Returns the number of reactions"""
        return len(self.vi_p)

    def __repr__ (self):
        return 'Rxn(ki={}, xi={}, vi_p={}, vi_dp={})'.format(self.ki, self.xi,
                                                             self.vi_p,
                                                             self.vi_dp)

    def progress_rate (self):
        """ Calculates and returns the progress rate

        Method is not implemented in this base class, but should be
        implemented in its subclasses.
        """
        raise NotImplementedError('Subclass must implement this method')

    def reaction_rate (self):
        """ Calculates and returns the reaction rate

        Method is not implemented in this base class, but should be
        implemented in its subclasses.
        """
        raise NotImplementedError('Subclass must implement this method')


class ElemRxn(RxnBase):
    """Class of elementary reactions
    Subclass of Rxn
    """
    pass


class NonElemRxn(RxnBase):
    """Class of non-elementary reactions
    Subclass of Rxn
    """
    pass


class RevElemRxn(ElemRxn):
    """Class of reversible elementary reactions
    Subclass of ElemRxn
    """
    pass


class IrreversibleElementaryRxn(ElemRxn):
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

    self.vi_p: a list of floats, optional, default value = [[1.0, 1.0, 0.0],
    [0.0, 0.0, 1.0]]
        Stoichiometric coefficients of the reactants
        Initialized with constructor argument vi_p

    self.vi_dp: a list of floats, optional, default value = [[0.0, 0.0, 1.0],
    [1.0, 1.0, 0.0]]
        Stoichiometric coefficients of the products
        Initialized with constructor argument vi_dp

    NOTES
    =====
    PRE:
        - ki, xi, vi_p, and vi_dp have list or np.array type
        - xi is ordered in the form [[A], [B], [C]]
        - vi_p is ordered in the form [[v_11', v_21', v_31'], [v_12', v_22',
        v_32']]
        - vi_dp is ordered in the form [[v_11", v_21", v_31"], [v_12", v_22",
        v_32"]]
        - four or fewer inputs
    METHODS:
    ========
    __len__(): Returns the number of reactions
    __repr__(): Prints the class name with its attributes
    progress_rate(): Calculates and returns the progress rate
    reaction_rate(): Calculates and returns the reaction rate
    """

    def __init__ (self, ki=[10.0, 10.0], xi=[1.0, 1.0, 1.0],
                  vi_p=[[1.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                  vi_dp=[[0.0, 0.0, 1.0], [1.0, 1.0, 0.0]]):
        self.ki = ki
        self.xi = xi
        self.vi_p = vi_p
        self.vi_dp = vi_dp
        self.wi = None
        self.rates = None

    # def __len__(self):
    #     """Returns the number of species in the reaction"""
    #     return len(self.xi)

    def __len__ (self):
        """Returns the number of reactions"""
        return len(self.vi_p)

    def __repr__ (self):
        return 'IrrevElemRxn(ki={}, xi={}, vi_p={}, vi_dp={})'.format(self.ki,
                                                                      self.xi,
                                                                      self.vi_p,
                                                                      self.vi_dp)

    def progress_rate (self):
        """
        Returns the progress rate w for a system of irreversible elementary
        reaction
        
        RETURNS
        ========
        self.wi: a numpy array of floats,
           Has the form self.w unless self.ki <= 0 or self.xi < 0
           in which cases a ValueError exception is raised
        
        NOTES
        =====
        POST:
             - self.ki, self.xi, self.vi_p, and self.vi_dp are not changed by
             this function
             - raises a ValueError exception if any(self.ki <= 0)
             - raises a ValueError exception if any(self.xi < 0)
             - returns a numpy array of floats of progress rate self.wi
        
        EXAMPLES
        =========
        >>> IrrevElemRxn([10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0,0.0, 2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]]).progress_rate()
        array([ 40.,  10.])
        """

        # check value conditions
        if (any(i <= 0 for i in self.ki)):  # check reaction coefficients
            raise ValueError("reaction rate coefficients ki must be positive.")
        elif (any(i < 0 for i in self.xi)):  # check concentration array
            raise ValueError("concentrations xi cannot be negative.")
        else:
            # Convert inputs into numpy column vectors/matrices
            xi = np.array(self.xi).reshape(-1, 1)
            vi_p = np.matrix(self.vi_p).T
            ki = np.array(self.ki).reshape(-1, 1)

            prodi = np.prod(np.power(xi, vi_p),
                            axis=0)  # calculate the product of xi^(vi_p) for
            #  each reaction
            self.wi = np.squeeze(
                  ki * np.array(prodi).reshape(-1, 1))  # multiply ki by prodi

            return self.wi

    def reaction_rate (self):
        """Returns the progress rate w for a system of irreversible
        elementary reaction

        RETURNS
        ========
        self.rates: a numpy array of floats,
           Has the form self.rates unless self.ki <= 0 or self.xi < 0
           in which cases a ValueError exception is raised
    
        NOTES
        =====
        POST:
             - self.ki, self.xi, self.vi_p, and self.vi_dp are not changed by
             this function
             - raises a ValueError exception if any(self.ki <= 0)
             - raises a ValueError exception if any(self.xi < 0)
             - returns a numpy array of floats of reaction rates, self.rates
        
        EXAMPLES
        =========
        >>> IrrevElemRxn([10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0,0.0, 2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]]).reaction_rate()
        array([-60., -70.,  70.])
        """
        # check value conditions
        if any(i <= 0 for i in self.ki):  # check reaction coefficients
            raise ValueError("reaction rates ki must be positive.")
        elif any(i < 0 for i in self.xi):  # check concentration array
            raise ValueError("concentrations xi cannot be negative.")
        else:
            # Convert inputs into numpy column vectors/matrices
            vi_p = np.matrix(self.vi_p).T
            vi_dp = np.matrix(self.vi_dp).T

            w = self.progress_rate()  # calculate progress rate
            vi = vi_dp - vi_p  # calculate overall stoicheometric coefficients

            self.rates = np.squeeze(
                  np.array(np.dot(vi, w)))  # calculate reaction rate

            return self.rates


class RxnData():
    """ Container for individual reaction data.

    Attributes
    ----------
    rxn_id : str
        id attribute of <reaction> element.
    reversible : bool
        True if reaction is reversible; False if irreversible.
    reactants : Dict[str, int]
        Mapping of species to stoichiometric coefficients for reactants.
        Example: {'H2':1, 'O':1}
    products : Dict[str, int]
        Mapping of species to stoichiometric coefficients for products.
    rate_coeff : List[float] or float
        Reaction rate coefficients contained depend on the type of rate
        coefficients in the XML file, dictated by the child of the <rateCoeff>
        element:
            Arrhenius: [A, E]
            modifiedArrhenius: [A, b, E]
            Constant: k
    type : RxnType
        Enum value for reaction type.
    """

    def __init__ (self, rxn_id=None, reversible=None, reactants=None,
                  products=None, rate_coeff=None, type=None):
        self.rxn_id = rxn_id
        self.reversible = reversible
        self.reactants = reactants
        self.products = products
        self.rate_coeff = rate_coeff
        self.type = type

    def equation (self):
        """ Returns equation representation of reactants and products.

        Species are listed in alphabetical order on both reactant and product
        side.
        """
        reactant_side = self.__build_equation_side(self.reactants)
        product_side = self.__build_equation_side(self.products)
        return '{0} [=] {1}'.format(reactant_side, product_side)

    @staticmethod
    def __build_equation_side (species_dict):
        result = ''
        species_count = -1
        for species, conc in sorted(species_dict.items()):
            species_count += 1

            # Handle first species differently than rest, then proceed to
            # next species.
            if species_count == 0:
                if abs(conc) == 1:
                    conc_part = ''
                else:
                    conc_part = conc
                result += '{0}{1}'.format(conc_part, species)
                continue

            # Non-first species.
            if conc < 0:
                operator = ' - '
            else:
                operator = ' + '

            if abs(conc) == 1:
                conc_part = ''
            else:
                conc_part = abs(conc)

            result += '{0}{1}{2}'.format(operator, conc_part, species)

        return result


class XmlParser():
    """ Core method load() produces list of RxnData from XML file contents.

    Notes
    -----
    ChemKinError raised when invalid values are encountered in the XML file.
    """

    def __init__ (self, path):
        """ Ensures path contains .xml file extension. """
        if path[-4:] != '.xml':
            path += '.xml'
        self.path = path

    def load (self):
        """ Parses XML file contents to create list of RxnData objects
        representing the reactions in the file.
        """
        tree = ET.parse(self.path)
        root = tree.getroot()
        species = []
        for species_i in root.find('phase'):
            new_species = species_i.text.strip().split()
            species.extend(new_species)

        results = []
        for rxn in root.find('reactionData').findall('reaction'):
            rxn_data = self.__extract_data_from_reaction_element(rxn)
            results.append(rxn_data)
        return species, results

    def __extract_data_from_reaction_element (self, rxn):
        """ Returns RxnData object containing data from <reaction> XML element.
        Raises ChemKinError for invalid attribute/element values.
        """
        result = RxnData()
        result.rxn_id = rxn.get('id')

        # reversible
        reversible = rxn.get('reversible').lower().strip()
        if reversible in ['no', 'n', 'false', 'f']:
            result.reversible = False
        elif reversible in ['yes', 'y', 'true', 't']:
            result.reversible = True
        else:
            raise ChemKinError(
                  'XmlParser.load()',
                  'Invalid reversible attribute in reaction {}'.format(
                      result.rxn_id))

        # type
        rnx_type = rxn.get('type').lower().strip()
        if rnx_type == 'elementary':
            result.type = RxnType.Elementary
        else:
            raise ChemKinError(
                  'XmlParser.load()',
                  'Invalid type attribute in reaction {}'.format(
                      result.rxn_id))

        # rate_coeff
        rate_coeff = rxn.find('rateCoeff')
        if rate_coeff is None:
            raise ChemKinError('XmlParser.load()',
                               'No <rateCoeff> element found in one of the '
                               'reactions.')

        if rate_coeff.find('Arrhenius') is not None:
            arrhenius = rate_coeff.find('Arrhenius')
            A = float(arrhenius.find('A').text.strip())
            if A < 0:
                raise ChemKinError('XmlParser.load()',
                      'A coeff < 0 in reaction with id = {}'.format(
                            result.rxn_id))
            E = float(arrhenius.find('E').text.strip())
            result.rate_coeff = [A, E]
        elif rate_coeff.find('modifiedArrhenius') is not None:
            mod_arrhenius = rate_coeff.find('modifiedArrhenius')
            A = float(mod_arrhenius.find('A').text.strip())
            if A < 0:
                raise ChemKinError(
                      'A coeff < 0 in reaction with id = {}'.format(
                            result.rxn_id))
            b = float(mod_arrhenius.find('b').text.strip())
            E = float(mod_arrhenius.find('E').text.strip())
            result.rate_coeff = [A, b, E]
        elif rate_coeff.find('Constant') is not None:
            const = rate_coeff.find('Constant')
            result.rate_coeff = float(const.find('k').text.strip())
        else:
            raise ChemKinError('XmlParser.load()',
                               'No recognized child of <rateCoeff> found '
                               'from which to parse coefficients.')

        # reactants / products
        result.reactants = self.__map_conc_to_species(rxn.find('reactants'))
        result.products = self.__map_conc_to_species(rxn.find('products'))

        return result

    @staticmethod
    def __map_conc_to_species (tag):
        """ Creates dict mapping the species to concentration where tag is
        either a <reactants> or <products> XML tag.
        """
        result = {}
        for item in tag.text.strip().split():
            species, conc = item.strip().split(':')
            conc = int(conc)
            species = species.upper()
            result[species] = conc
        return result
