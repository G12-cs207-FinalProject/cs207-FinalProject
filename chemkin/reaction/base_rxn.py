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
        return 'RxnBase(ki={}, xi={}, vi_p={}, vi_dp={})'.format(self.ki, self.xi,
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