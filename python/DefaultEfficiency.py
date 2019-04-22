import autograd.numpy as np
from PhysicalParameters import PhysicalParameters

class DefaultEfficiency():
    """
    DefaultEfficiency


    A flat efficiency with a time-acceptance between 0 and 4 KS lifetimes
    """
    def __init__(self, s12, s13):
        self.base_eff = np.ones((len(s12), len(s13)))
        self.param = PhysicalParameters()
        self.t_min = 0
        self.t_max = 4*self.param.KL_tau

    def get_eff(self, t):
        if t<self.t_min or t > self.t_max:
            return 0*self.base_eff
        return self.base_eff

    def default_time_range(self):
        # Returns the 
        # return np.array([0])

        time_range = np.linspace(0, 5*self.param.KS_tau, 100)
        time_range = np.append(time_range,  np.linspace(5*self.param.KS_tau, 4*self.param.KL_tau, 50))
        return time_range

    def get_time_averaged_eff(self):
        # return the efficiency average over the default time range
        # weighted with exp(t/tau_KS)
        ts = self.default_time_range()
        dts = ts[1:]-ts[:-1]
        dts = np.append(dts[0], dts)
        effs_at_ts = [ self.get_eff(t)*dts[i]*np.exp(-t/self.param.KS_tau) for i, t in enumerate(ts) ]
        effs = np.array(effs_at_ts).sum(0)
        return effs/effs.sum()*effs.size # return efficiency relative to average (which 1/effs.size)



    def __call__(self, t):
        return self.get_eff(t)

if __name__ == "__main__":
    s12 = [1, 2]
    d = DefaultEfficiency(s12, s12)
    print np.shape(d.get_time_averaged_eff())
    print (d.get_time_averaged_eff())