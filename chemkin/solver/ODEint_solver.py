"""
The ODEint_solver.py contains a class ODE_int_solver which will return time
evolution of a specie i over specified time range t.
"""
import numpy as np
from scipy.integrate import odeint
from chemkin.reaction.elementary_rxn import ElementaryRxn


class ODE_int_solver():
    """Integrates the time evolution of a concentration of a specie i over
    specified time range.

    Attributes:
        critical_t (float): Stores time at which reaction reaches equilibrium
            after integrating.
        equil_thresh (float, default 0.01): The reaction is defined to reach
            equilibrium once np.linalg.norm(prograte_diff) < equil_thresh,
            where prograte_diff is the difference in forward and backward
            progress rates for reaction.
        critical_i (int): Number of time steps used by ODE solver before
            solving critical_t.
    """

    def __init__ (self, temp, xi, ki, b_ki, sys_vi_p, sys_vi_dp, equil_thresh
    =0.01):

        self.temp = temp  # temperature will be held constant in the calculation
        self.critical_t = None
        self.xi = xi  # need initial concentration of ALL species
        self.ki = ki
        self.b_ki = b_ki
        self.sys_vi_p = sys_vi_p
        self.sys_vi_dp = sys_vi_dp
        self.equil_thresh = equil_thresh
        self.critical_i = None

    def solve (self, time_int):
        """Solves evolution of specie concentration over specified time range.

        Args:
            time_int (List[float]): Time steps over which to solve evolution
                of concentration.

        Returns:
            sol (array, shape (len(t), len(y0)): Results of calling
                scipy.integrate.odeint().
            critical_t (float): Time at which reaction reaches equilibrium
                after integrating.
        """

        # for starters, focus on time evolution of ALL species
        rxn = ElementaryRxn(self.ki, self.b_ki, self.xi, self.sys_vi_p,
                            self.sys_vi_dp)

        i = 0  # iteration count

        def rxn_rate (x, t):
            nonlocal  i
            rxn.xi = x
            if i != 0:
                if np.linalg.norm(rxn.b_wi - rxn.f_wi) < self.equil_thresh:
                    if self.critical_t is None:
                        self.critical_t = t
                        self.critical_i = i
            i += 1
            return rxn.reaction_rate()

        sol = odeint(func=rxn_rate, y0=self.xi, t=time_int, mxstep=5000000)
        return sol, self.critical_t
