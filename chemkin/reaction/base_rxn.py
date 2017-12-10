import numpy as np
from chemkin.solver.ODEint_solver import ODE_int_solver

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

    def __init__(self, ki, b_ki, xi, vi_p, vi_dp):
        self.ki = ki
        self.b_ki = b_ki
        self.xi = xi
        self.vi_p = vi_p
        self.vi_dp = vi_dp
        self.wi = None
        self.f_wi = None
        self.b_wi = None
        self.rates = None
        self.is_reversible = []

        # Determine if reaction i is reversible
        for b_k in b_ki:
            if b_k != 0:
                self.is_reversible.append(True)
            else:
                self.is_reversible.append(False)

    # def __init__ (self, ki, xi, vi_p, vi_dp):
    #     self.ki = ki
    #     self.xi = xi
    #     self.vi_p = vi_p
    #     self.vi_dp = vi_dp
    #     self.wi = None
    #     self.rates = None

    # def __len__(self):
    #     """Returns the number of species in the reaction"""
    #     return len(self.xi)

    def __len__ (self):
        """Returns the number of reactions"""
        return len(self.vi_p)

    def __repr__ (self):
        return 'RxnBase(ki={}, b_ki={}, xi={}, vi_p={}, vi_dp={})'.format(self.ki, self.b_ki, self.xi, self.vi_p, self.vi_dp)

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


    def species_concentration(self, T, end_t, n_steps=101):
        """ Return the list of the species concentration at Temperatrue = T and time = end_t
        """
        time_steps = np.linspace(0, end_t, n_steps)
        # solver = ODE_int_solver(T, self.xi, self.ki, self.b_ki, self.vi_p, self.vi_dp)
        solver = ODE_int_solver(T, self)
        sol, _, _ = solver.solve(time_steps)
        return sol[-1, :]

    def species_concentration_evolution(self, T, end_t, n_steps=101):
        """ Return the list of the species concentration evolution at Temperatrue = T and from start to end_t
        """
        time_steps = np.linspace(0, end_t, n_steps)
        # solver = ODE_int_solver(T, self.xi, self.ki, self.b_ki, self.vi_p, self.vi_dp)
        solver = ODE_int_solver(T, self)
        sol, _, _ = solver.solve(time_steps)
        return sol

    def time_to_equilibrium(self, T, n_steps=101):
        """ Return the list of time to equilibrium of all the reactions and the time to equilibrium of the overall system
        """
        end_t = 1e10
        time_steps = np.linspace(0, end_t, n_steps)
        solver = ODE_int_solver(T, self)
        _, critical_t, overall_critical_t = solver.solve(time_steps)
        return end_t, critical_t, overall_critical_t


