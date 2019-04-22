from DefaultEfficiency import DefaultEfficiency
from PhysicalParameters import PhysicalParameters
import numpy as np


class TimeDependentEfficiency(DefaultEfficiency):
    """docstring for TimeDependentEfficiency"""

    def __init__(self, s12, s13, p, z_min=0, z_max = 2200):
        self.base_eff = np.ones((len(s12), len(s13)))

        default_param = PhysicalParameters()
        self.param = default_param

        c = 300 # mm/ns

        self.t_min = z_min * default_param.mK / (c * p)
        self.t_min[self.t_min == None] = 0

        self.t_max = z_max * default_param.mK / (c * p)
        self.t_max[self.t_max == None] = 0

    def get_eff(self, t):
        eff = np.copy(self.base_eff)
        eff[t<self.t_min] = 0
        eff[t>self.t_max] = 0
        return eff

    def default_time_range(self):
        # Returns the 
        # return np.array([0])
        return np.linspace(np.nanmin(self.t_min), np.nanmax(self.t_max), 150)
