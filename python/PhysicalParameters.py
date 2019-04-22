import math, cmath, sys
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

class PhysicalParameters(object):
    """PhysicalParameters is just a container that holds the physical parameters which 
    governs the kaon regeneration

    ALL BASIC parameters are just object properties and can be set to non-default values with
    pp = PhysicalParameters() # init with default values
    pp.A = 197 # calculate mixing in gold-197 instead of default aluminium-27

    OF SPECIAL INTEREST is 
    pp.include_mixing = False # completely ignore mixing

    IMPORTANT!: When updating a parameter call
    pp.update()

    to make sure that all calculated parameters are updated accordingly.
    (sadly there is no way to enforce this since everything is public in python)

    """
    def __init__(self, 
        include_mixing=True, 
        use_low_p_dchi=False,
        sum_chi=None):


        self._A = 27. # Aluminium for a start
        self._p = 30. # GeV/c
        self._arg_f = -124.7/180.*math.pi 

        #rho = 20.70e-3*1e6 # 10 * density of aluminium, in kg/m^3
        self._rho = 0.02*2.70 # 2 % aluminium, in g/cm^3
        self.N_A = 6.02e+23 # in 
        self.mK = 0.497 # neutral kaon mass in GeV/c^2
        self.c = 30 # speed of light in cm/ns

        # not including the minus from Eq (61) because expression for df gives (f_bar - f)
        self._g_S = 1/0.0896
        self._g_L = 1/51.2
        self._m_S = 497.6 # absolute size does not matter, only delta m
        self._dm  = 5.290 # 0.474*self.g_S, unit is in ns

        self.include_mixing = include_mixing

        self.default_eps = cmath.rect(2.2e-3, math.pi/4)

        self._get_dchi = self.get_dchi_simple
        if use_low_p_dchi:
            self.load_dsigma_proton()
            self._get_dchi = self.get_dchi_low_p

        self._sum_sig = 0
        if sum_chi=="low_p":
            self._sum_sig = 2*373
        elif sum_chi=="high_p":
            self._sum_sig = 713+428
        
        self.update()

    # Inputs defined as properties to make sure
    # update() is called whenever they are changed, so that
    # derived quantities are recalculated

    @property
    def A(self):
        return self._A  # Aluminium for a start
    @A.setter
    def A(self, A):
        self._A = A
        self.update()

    @property
    def p(self):
        return self._p  # GeV/c
    @p.setter
    def p(self, p):
        self._p = p
        self.update()

    @property
    def arg_f(self):
        return self._arg_f  
    @arg_f.setter
    def arg_f(self, arg_f):
        self._arg_f = arg_f
        self.update()
                #rho = 20.70e-3*1e6 # 10 * density of aluminium, in kg/m^3
    @property
    def rho(self):
        return self._rho # 2 % aluminium, in g/cm^3
    @rho.setter
    def rho(self, rho):
        self._rho = rho
        self.update()
                # not including the minus from Eq (61) because expression for df gives (f_bar - f)
    @property
    def g_S(self):
        return self._g_S 
    @g_S.setter
    def g_S(self, g_S):
        self._g_S = g_S
        self.update()

    @property
    def g_L(self):
        return self._g_L 
    @g_L.setter
    def g_L(self, g_L):
        self._g_L = g_L
        self.update()

    @property
    def m_S(self):
        return self._m_S  # absolute size does not matter, only delta m
    @m_S.setter
    def m_S(self, m_S):
        self._m_S = m_S
        self.update()

    @property
    def dm(self):
        return self._dm   # 0.474*self.g_S, unit is in ns
    @dm.setter
    def dm(self, dm):
        self._dm = dm
        self.update()


    def get_dchi_simple(self):
        self.dsig = 23.2*self.A**(0.758)/self.p**0.614*1e-27 # in GeV/c * cm^-2
        self.Im_f = (self.p)/(4*math.pi)*self.dsig
        self.abs_f = self.Im_f/math.sin(self.arg_f)
        self.delta_f = self.abs_f*np.cos(self.arg_f) + 1j * self.abs_f*np.sin(self.arg_f)
        self.dchi = 2*math.pi*self.N/self.mK*self.delta_f*self.c # in 1/ns
        return self.dchi

    def get_dchi_low_p(self):
        try:
            dsig_m = self.sig_f['kmp'](self.p)
            dsig_p = self.sig_f['kpp'](self.p)
        except ValueError as e:
            print e
            print np.min(self.p)
            sys.exit()
        self.dsig = np.power(self.A, 0.758)/(1+1.25*np.exp(-1.814*self.p))*(dsig_m - dsig_p)*1e-27
        self.Im_f = (self.p)/(4*math.pi)*self.dsig
        self.abs_f = self.Im_f/math.sin(self.arg_f)
        self.delta_f = self.abs_f*np.cos(self.arg_f) + 1j * self.abs_f*np.sin(self.arg_f)
        self.dchi = 2*math.pi*self.N/self.mK*self.delta_f*self.c # in 1/ns
        return self.dchi

    def update(self):
        self.N = self.rho * self.N_A/self.A # scattering density in cm^-3

        self.dchi = self._get_dchi()
        if not self.include_mixing:
            self.dchi = 0

        self.KS_tau = 1/self.g_S
        self.KL_tau = 1/self.g_L

        self.dg  = self.g_L - self.g_S
        self.m_L = self.m_S + self.dm

        self.l_S = self.m_S - 1j/2*self.g_S
        self.l_L = self.m_L - 1j/2*self.g_L
        self.dl  = self.l_L - self.l_S

        # Use the total absorbtion set in __init__
        self.sum_f_im = -self.p/(4*math.pi)*self._sum_sig*1e-27
        self.sum_chi = 1j*self.sum_f_im*2*math.pi*self.N/self.mK*self.c 

        self.Sigma = 0.5*(self.l_S + self.l_L + self.sum_chi)# + self.chi + self.chib)
        self.Omega = 0.5*np.sqrt(self.dl**2 + self.dchi**2)

    def load_dsigma_proton(self):
        self.df = {}
        args = {'sep':"\s+", 
                'skiprows':11, 
                'usecols':[1,4,5,6,7,8,9], 
                'names':['P', 'SIG', 'STA+', 'STA-', 'SYS+', 'SYS-', 'FLAG']}
        # load PDG data 
        self.df["kpp"] = pd.read_csv("../pdg_data/rpp2016-kpp_total.dat", **args)
        self.df["kmp"] = pd.read_csv("../pdg_data/rpp2016-kmp_total.dat", **args)

        self.new_df = {}

        def wm(x, w):
            return np.average(x, weights=w)
        def wme(x, w):
            if len(x)==1:
                return np.sqrt(1/(w))
            wm_x = wm(x, w)
            var_x = np.sum(np.power(x-wm_x,2)*w)
            sum_w = np.sum(w)
            sum_w2 = np.sum(np.power(w,2))
            
            return np.sqrt(var_x*(sum_w2)/(sum_w*sum_w) * sum_w / (sum_w*sum_w-sum_w2))

        for k in self.df.keys():
            self.df[k]['TOTAL ERR'] = (self.df[k]['STA+'] + self.df[k]['STA-'])/2 + self.df[k]['SIG']*(self.df[k]['SYS+'] + self.df[k]['SYS-'])/(200.)
            self.df[k]['WEIGHT'] = 1/(self.df[k]['TOTAL ERR']*self.df[k]['TOTAL ERR'])
            lambda_wm = lambda x: wm(x, self.df[k].loc[x.index, 'WEIGHT'])
            lambda_wme = lambda x: wme(x, self.df[k].loc[x.index, 'WEIGHT'])
            
            f = {'SIG': {'MEAN' : lambda_wm, 'ERR': lambda_wme} }
            # Groupby and aggregate with your dictionary:
            self.new_df[k] = self.df[k].groupby(['P'], group_keys=False).agg(f)
            self.new_df[k].columns = list(map('_'.join, self.new_df[k].columns.values))
            self.new_df[k]['P'] = self.new_df[k].index

        self.sig_f = {}
        for k in self.df.keys():
            self.sig_f[k] = interp1d(self.new_df[k]['P'], self.new_df[k]['SIG_MEAN'],
                bounds_error=False, fill_value='extrapolate')

if __name__ == "__main__":
    p = PhysicalParameters(use_low_p_dchi=True, sum_chi=None)
    rho_Si = 2.32
    
    print "Default:"
    p.p = 23.4
    p.rho = 0.024*rho_Si
    print '  |r| =', 0.5*abs(p.dchi/p.dl) 
    print '  ph  =', np.angle((p.dchi/p.dl) )/np.pi*180.

    print "LHCb LL:"
    p.p = 23.4
    p.rho = 0.024*rho_Si
    print '  |r| =', 0.5*abs(p.dchi/p.dl) 
    print '  ph  =', np.angle((p.dchi/p.dl) )/np.pi*180.
    
    print "\nLHCb DD:"
    p.p = 36.6
    p.rho = 0.016*rho_Si
    print '  |r|=', 0.5*abs(p.dchi/p.dl) 
    print '  ph  =', np.angle((p.dchi/p.dl) )/np.pi*180.

    print "\nBelle II:"
    p.p = 1.1
    p.rho = 0.040*rho_Si
    print '  |r|=', 0.5*abs(p.dchi/p.dl) 
    print '  ph  =', np.angle((p.dchi/p.dl) )/np.pi*180.
    

