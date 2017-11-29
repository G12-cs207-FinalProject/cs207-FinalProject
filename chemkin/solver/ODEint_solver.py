"""
The ODEint_solver.py contains a class ODE_int_solver which will return time
evolution
 of a specie i over specified time range t.
"""
from scipy.integrate import odeint

from chemkin.reaction.elementary_rxn import ElementaryRxn


class ODE_int_solver():
    # Integrates the time evolution of a concentration of a specie i over
    # specified time range

    def __init__ (self, temp, xi, ki, b_ki, sys_vi_p, sys_vi_dp):

        self.temp = temp  # temperature will be held constant in the calculation
        self.critical_t = None
        self.xi = xi  # need initial concentration of ALL species
        self.ki = ki
        self.b_ki = b_ki
        self.sys_vi_p = sys_vi_p
        self.sys_vi_dp = sys_vi_dp

    def equil_met(self):
        return False

    def solve (self, time_int):

        # for starters, focus on time evolution of ALL species
        rxn = ElementaryRxn(self.ki, self.b_ki, self.xi, self.sys_vi_p,
                            self.sys_vi_dp)


        def rxn_rate (x, t):
            rxn.xi = x
            if self.equil_met():
                if self.critical_t is None:
                    self.critical_t = t

            return rxn.reaction_rate()

        sol = odeint(func=rxn_rate, y0=self.xi, t=time_int, mxstep=5000000)

        return sol, self.critical_t
