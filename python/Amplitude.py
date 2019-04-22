import cmath, math, os
import numpy as np
from UtilityFunctions import load_amplitude
from PhysicalParameters import PhysicalParameters

class Amplitude:
    """
    Calculates the KS and KL amplitudes at a given time, given the values at t = 0

    Takes the arguments
        KS amplitude at t = 0, 2D numpy array over phase space called A_L/S (same size ofc)
        s12 1D array with the s12 values for which amplitudes are given in KS
        s13 1D array with the s13 values for which amplitudes are given in KS
        optional: KL amplitude at t = 0, 2D numpy array over phase space called A_L/S (same size ofc)
                    if not passed, all KL contributions are ignored
        optional: eps: complex number for eps K0 CPV parameter, default is 2e-3 * exp(i*pi/4) which is close to actual value

    No assumptions are made on whether the amplitudes are for D0 or D0b or whether they are for
    KS = K1 or KS = K1 + eps*K2. That all goes into defining the starting amplitudes, and should be handled 
    elsewhere. ie. THIS CLASS KNOWS NOTHING ABOUT GAMMA OR KS CPV - it can only do flavour tagged D0 decays.
    """

    def __init__(self, A_KS, s12, s13, A_KL=np.array([0]), eps=PhysicalParameters().default_eps, parameters=PhysicalParameters()):
        self.eps = eps
        self.parameters = parameters
        if not A_KL.any(): # .any() returns true if a single element returns is non-zero
            A_KL = A_KS*0 # ensure same dimensions
            self.eps = 0. # if no KL is supplied, set eps = 0 for consistency

        # These maps hold the cashed amplitude values as a function of time
        # each amplitude is a 2D numpy array holding python complex numbers
        self.A_KS = LimitedSizeDict(size_limit=3) # save last 2 calculated times in cache
        self.A_KL = LimitedSizeDict(size_limit=3)
        self.A_KS_0 = A_KS
        self.A_KL_0 = A_KL

        #  save s coordinates corresponding to matrix
        self.s12 = s12 # s(Ks, pi+-)
        self.s13 = s13 # s(Ks, pi+-)

    def get_A(self, t):
        ''' returns the amplitude for a decay to pipi at time t
            that is A_KS + eps*A_KL
        '''
        if not self.A_KS.has_key(t):
            self.calculate_amplitudes(t)

        if self.eps != 0:
            return self.A_KS[t] + self.eps*self.A_KL[t]
        return self.A_KS[t]

    def get_A_KS(self, t):
        ''' returns the KS amplitude component for time t'''
        if not self.A_KS.has_key(t):
            self.calculate_amplitudes(t)
        return self.A_KS[t]

    def get_A_KL(self, t):
        ''' returns the KL amplitude component for time t'''
        if not self.A_KL.has_key(t):
            self.calculate_amplitudes(t)
        return self.A_KL[t]

    def calculate_amplitudes(self, t):
        '''Calculates amplitude for both KL and KS as they have a lot of computations in common'''

        # These factors are the same all over phsp and may as well be calculated once
        # print "AA:", self.parameters.Omega*t
        cosOt = np.cos(self.parameters.Omega*t)
        sinOt = np.sin(self.parameters.Omega*t)/(2*self.parameters.Omega)
        expSt = np.exp( -1j*self.parameters.Sigma*t)

        A0S  = self.A_KS_0
        A0L  = self.A_KL_0


        # phase space loop handled by numpy!
        self.A_KS[t] = expSt*(A0S*cosOt+1*1j*(self.parameters.dl*A0S-1*self.parameters.dchi*A0L)*sinOt)
        self.A_KL[t] = expSt*(A0L*cosOt-1*1j*(self.parameters.dl*A0L+1*self.parameters.dchi*A0S)*sinOt)

from collections import OrderedDict

class LimitedSizeDict(OrderedDict):
    def __init__(self, *args, **kwds):
        self.size_limit = kwds.pop("size_limit", None)
        OrderedDict.__init__(self, *args, **kwds)
        self._check_size_limit()

    def __setitem__(self, key, value):
        OrderedDict.__setitem__(self, key, value)
        self._check_size_limit()

    def _check_size_limit(self):
        if self.size_limit is not None:
            while len(self) > self.size_limit:
                self.popitem(last=False)

if __name__=="__main__":
    print ("Testing {}".format(os.path.basename(__file__)))
    A_KS = np.array([1])
    A_KL = np.array([1])
    # A_KL = np.array([0])
    pp = PhysicalParameters(False)
    a = Amplitude(A_KS, A_KL=A_KL, eps=0, parameters=pp)
    print pp.m_L-pp.m_S
    for t in range(5):
        a_KS = a.get_A_KS(t*pp.KS_tau)[0]
        a_KL = a.get_A_KL(t*pp.KS_tau)[0]
        a_KS_exp = np.exp(-t*pp.KS_tau*(1j*pp.m_S + 0.5*pp.g_S))
        a_KL_exp = np.exp(-t*pp.KS_tau*(1j*pp.m_L + 0.5*pp.g_L))

        print "{:1.7f}     {:1.7f}     {:1.7f}".format(t, a_KL, a_KL_exp)
    # print a.get_A_KL(9)

