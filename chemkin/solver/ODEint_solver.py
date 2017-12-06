"""
Contains class ODE_int_solver to compute the evolution of reaction species'
concentrations over specified time range.
"""
import numpy as np
from scipy.integrate import odeint

from chemkin.reaction.elementary_rxn import ElementaryRxn


class ODE_int_solver():
    """Integrates the time evolution of a concentration of a specie i over
    specified time range.

    Attributes:
        temp (float): Temperature for reaction (assumed to be held constant).
        species_equil_thresh (float, default 1e-5): Species concentration
            evolution  is defined to reach equilibrium once np.abs(bw - fw) <
            species_equil_thresh, where (bw - fw) is the difference in
            backward and forward progress rates for that species.
        overall_equil_thresh (float, default 1e-2): Overall reaction is
            defined to reach equilibrium once np.linalg.norm(prograte_diff) <
            overall_equil_thresh, where prograte_diff is the difference in
            the backward and forward progress rate vectors of the reaction.
        critical_t (list of floats of len(ki)): Stores time(s) at which
            reaction's component species' concentrations reach equilibrium.
        overall_critical_t (float): Stores time at which overall reaction
            reaches equilibrium.
        xi, ki, b_ki, sys_vi_p, sys_vi_dp: Constructor args for ElementaryRxn
            class.
    """

    def __init__ (self, temp, xi, ki, b_ki, sys_vi_p, sys_vi_dp, equil_thresh
    =1e-5, overall_equil_thresh=1e-2):
        self.temp = temp
        self.xi = xi
        self.ki = ki
        self.b_ki = b_ki
        self.sys_vi_p = sys_vi_p
        self.sys_vi_dp = sys_vi_dp
        self.species_equil_thresh = equil_thresh
        self.overall_equil_thresh = overall_equil_thresh
        # Init critical time attrs with dummy values.
        self.critical_t = -2*np.ones((len(ki),))
        self.overall_critical_t = -2

    def solve (self, time_int):
        """Solves evolution of specie concentration over specified time range.

        Args:
            time_int (list of floats): Time steps over which to solve evolution
                of species concentrations.

        Returns:
            sol (numpy array, shape (len(time_int), len(self.xi)): Results of
                calling scipy.integrate.odeint().
            self.critical_t
            self.overall_critical_t
        """
        rxn = ElementaryRxn(self.ki, self.b_ki, self.xi, self.sys_vi_p,
                            self.sys_vi_dp)

        i = 0  # iteration count

        def rxn_rate (x, t):
            nonlocal  i
            rxn.xi = x

            if i != 0:
                for j, (bw, fw) in enumerate(zip(rxn.b_wi, rxn.f_wi)):
                    if np.abs(bw - fw) < self.species_equil_thresh:
                        if self.critical_t[j] == -2:
                            self.critical_t[j] = t

                if np.linalg.norm(
                            rxn.b_wi - rxn.f_wi) < self.overall_equil_thresh:
                    if self.overall_critical_t == -2:
                        self.overall_critical_t = t
            i += 1
            return rxn.reaction_rate()

        sol = odeint(func=rxn_rate, y0=self.xi, t=time_int, mxstep=5000000)
        return sol, self.critical_t, self.overall_critical_t
