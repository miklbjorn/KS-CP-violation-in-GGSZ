from scipy.optimize import minimize
import autograd.numpy as np
import cmath, math
from DefaultEfficiency import DefaultEfficiency
from PhysicalParameters import PhysicalParameters
from DefaultModel import DefaultModel
from Binning import Binning
from Amplitude import Amplitude
import UtilityFunctions as uf
import sys, os

class FullModel(object):
    """Default Model calculates yields using the default setup: ie by calculating Fi and using the yield equations from
    the MI GGSZ LHCb papers.

    It does not include global asymmetries.

    """
    def __init__(self, 
        KS_amplitude_file, 
        KL_amplitude_file, 
        eps=PhysicalParameters().default_eps, 
        parameters=PhysicalParameters(),
        efficiency=None, 
        binning_file="input/KsPiPi_optimal.pickle",
        KL_to_pipi_factor=1,#factor to multiply on the KL->pipi amplitude
        amplitude_asym_factor=1): # Factor to multiply on the A_KS - A_K1 difference

        self.eps = eps
        self.parameters = parameters

        # save predicted yields, so they do not have to be recalculated if just doing statistical fluctuations
        self.predicted_yields = {}

        self.KS_amplitude_file = "not initialized" # will be set in self.update_amplitudes(), DONT CHANGE
        self.KL_amplitude_file = "not initialized" # will be set in self.update_amplitudes(), DONT CHANGE

        self.KL_to_pipi_factor = KL_to_pipi_factor
        self.amplitude_asym_factor = amplitude_asym_factor

        self.update_amplitudes(KS_amplitude_file, KL_amplitude_file)

        if efficiency==None:
            efficiency = DefaultEfficiency(
                self.s12, self.s13)
        self.efficiency = efficiency

        # The time range over which to do the amplitude integrals
        # Set by efficiency as default, as that is where the time-acceptance is defined
        self.time_range = self.efficiency.default_time_range()

        bin_def, bin_def_s12, bin_def_s13 = uf.load_binning(binning_file)


        self.binning = Binning(bin_def, bin_def_s12, bin_def_s13,
                self.s12, self.s13)


        return

    def update_amplitudes(self, KS_amplitude_file, KL_amplitude_file):
        if not (KS_amplitude_file == self.KS_amplitude_file):
            self.A_KS, self.s12,  self.s13  = uf.load_amplitude(KS_amplitude_file)
            self.KS_amplitude_file = KS_amplitude_file
        if not (KL_amplitude_file == self.KL_amplitude_file):
            self.A_KL, self.s12L, self.s13L = uf.load_amplitude(KL_amplitude_file)
            self.KL_amplitude_file = KL_amplitude_file

        if not (np.array_equal(self.s12, self.s12L) and np.array_equal(self.s13, self.s13L)):
            print "Exiting: Different Dalitz coordinates in files:"
            print "  KS: {}".format(KS_amplitude_file)
            print "  KL: {}".format(KL_amplitude_file)
            for i, x in enumerate(self.s12):
                if self.s12L[i] != x:
                    print i, self.s12[i], self.s12L[i], self.s13[i], self.s13L[i]
            sys.exit(-1)

        # Define the amplitudes for a D0 and D0bar decay
        # Here it is assumed that the input A_KS is actually for K1, the CP-even state. It is A(K1 -> D0 h+h-)
        #  and it is assumed that the input A_KL is actually for K2, the CP-odd  state. It is A(K2 -> D0 h+h-)
        # The eps appearing here is due to EFFECT 1: the NON-EXACT MIRROR ASYMMETRY of the D(bar)hh amplitude
        # (EFFECT 2 being the presence of direct KL->pipi)
        A_KS_D0    =  self.A_KS             - self.amplitude_asym_factor*self.eps*self.A_KL # KS in K1, K2 basis
        A_KL_D0    =  self.A_KL             - self.amplitude_asym_factor*self.eps*self.A_KS # KL in K1, K2 basis
        A_KS_D0bar =  self.A_KS.transpose() + self.amplitude_asym_factor*self.eps*self.A_KL.transpose() # Use CP sign when transposing for conjugate
        A_KL_D0bar = -self.A_KL.transpose() - self.amplitude_asym_factor*self.eps*self.A_KS.transpose()

        # The 2D dalitz coordinate dependent 
        self.A_D0    = Amplitude(A_KS_D0   , self.s12, self.s13, self.KL_to_pipi_factor*A_KL_D0   , self.eps, self.parameters)
        self.A_D0bar = Amplitude(A_KS_D0bar, self.s12, self.s13, self.KL_to_pipi_factor*A_KL_D0bar, self.eps, self.parameters)

        # delete the cached predicted yields, as they may change
        self.predicted_yields = {}

    def load_KL_amplitude(self, KL_amplitude_file):
        self.update_amplitudes(self.KS_amplitude_file, KL_amplitude_file)

    def load_KS_amplitude(self, KS_amplitude_file):
        self.update_amplitudes(KS_amplitude_file, self.KL_amplitude_file)

    def set_time_range(self, time_range):
        # set the time range over which to integrate
        # should be a numpy 1D array
        self.time_range = time_range
        # delete the cached predicted yields, as they may change
        self.predicted_yields = {}

    def predict_yields(self, physics_param, Ntot):
        gamma  = physics_param[0] # in RADIANS!
        rB     = physics_param[1]
        deltaB = physics_param[2] # in RADIANS!

        key = "{}_{}_{}_{}".format(gamma, rB, deltaB, Ntot)

        # If already calculated, no reason to do it again
        if self.predicted_yields.has_key(key):
            # print "ALREADY CALCULATED"
            return self.predicted_yields[key]
        
        # print "CALCULATING"
        # If not, we'll get to it
        bin_num = self.binning.get_number_of_bins()
        Np_tot = np.zeros(bin_num*2)
        Nm_tot = np.zeros(bin_num*2)

        dts = self.time_range[1:] - self.time_range[:-1]
        dts = np.append(dts[0], dts)
        for i, t in enumerate(self.time_range): # units etc handled when setting time range
            dt = dts[i]

            A_D0_t    = self.A_D0   .get_A(t)
            A_D0bar_t = self.A_D0bar.get_A(t)

            A_Bp_t = A_D0bar_t + cmath.rect(rB, (deltaB+gamma))*A_D0_t
            A_Bm_t = A_D0_t    + cmath.rect(rB, (deltaB-gamma))*A_D0bar_t

            # A_Bp_sq = abs(A_Bp_t)**2*self.efficiency.get_time_averaged_eff()
            # A_Bm_sq = abs(A_Bm_t)**2*self.efficiency.get_time_averaged_eff()
            A_Bp_sq = abs(A_Bp_t)**2*self.efficiency.get_eff(t)
            A_Bm_sq = abs(A_Bm_t)**2*self.efficiency.get_eff(t)

            Np = np.array([])
            Nm = np.array([])
            for i in range(-bin_num, bin_num+1):
                if i==0: continue
                bin_idx = self.binning.get_bin_indices(i)
                Np = np.append(Np, np.sum(A_Bp_sq[bin_idx]))
                Nm = np.append(Nm, np.sum(A_Bm_sq[bin_idx]))
            Np_tot += Np*dt # add time-slide array to totals array
            Nm_tot += Nm*dt

        debug_normalize_plus_minus_separately = False
        # for debugging purposes it can be useful to generate an equal amount of B+ and B-
        # this can be enforced by setting the flag above to True
        if debug_normalize_plus_minus_separately:
            Np_tot = Np_tot/sum(Np_tot)
            Nm_tot = Nm_tot/sum(Nm_tot)

        yields = np.append(Np_tot, Nm_tot)
        yields = Ntot*yields/sum(yields) # normalize

        self.predicted_yields[key] = yields # save so it will not have to be recalculated     
        
        return yields

    def get_Fi_hat(self):
        '''Returns Fi as they would be measured in a control channel, for the physics
        situation and parameters relevant for the FullModel'''
        yields = self.predict_yields([0., 0., 0.], 1000) # no cabibbo suppressed part, ~ flavour tagged D
        Fi_hat = (yields[16:] + yields[15::-1])/sum(yields) #B- yields are in right order, B+ yields should be reversed
        return Fi_hat

if __name__ == "__main__":
    print ("Testing {}".format(os.path.basename(__file__)))
    parameters = PhysicalParameters(False)
    fm = FullModel(
        "../amplitude_calculations/output/KS_default.pickle",
        "../amplitude_calculations/output/KL_default.pickle",
        1.*cmath.rect(2.2e-3, math.pi/4),
        parameters)
    dm = DefaultModel("../amplitude_calculations/output/KS_default.pickle")
    physics_param = [uf.deg_to_rad(75.), 0.1, uf.deg_to_rad(90.)]
    
    data = fm.predict_yields(physics_param, 4000)
    print "DATA ASYM:", (sum(data[16:])-sum(data[:16]))/sum(data)
    print "B+:", (sum(data[:16]))
    print "B-:", (sum(data[16:]))
    # res = dm.fit(data)
    # print res

    # fm.load_KL_amplitude("../amplitude_calculations/output/KL_seed_13.pickle")
    # data = fm.predict_yields(physics_param, 4000)
    # res = dm.fit(data)

    res = dm.fit_gamma_from_asymmetry(data, 0.1, uf.deg_to_rad(90.))
    print res
    print uf.rad_to_deg(res.x[0])
    data_asym = (sum(data[16:])-sum(data[:16]))/sum(data)

    asym_dm = []
    asym_fm = []
    asym_fm2 = []
    plot_end = 90
    gms = np.linspace(0,plot_end, 6)
    fit_gms=[]

    from matplotlib import pyplot as plt

    for g in gms:
        # print g
        asym_dm.append(dm.predict_asym([uf.deg_to_rad(g), 0*0.1, uf.deg_to_rad(130.)]))
        data = fm.predict_yields([uf.deg_to_rad(g), 0*0.1, uf.deg_to_rad(130.)], 4000)
        asym_fm.append((sum(data[16:])-sum(data[:16]))/sum(data))
        res = dm.fit_gamma_from_asymmetry(data, 0*0.1, uf.deg_to_rad(130.), gamma0=g-5)
        fit_gms.append(uf.rad_to_deg(res.x[0]))
        # print res.message
    # plt.plot(gms, asym_dm, gms, asym_fm, [0, plot_end], [data_asym, data_asym], [uf.rad_to_deg(res.x[0]), uf.rad_to_deg(res.x[0])], [0, data_asym])
    plt.plot(gms, fit_gms-gms)
    plt.show()


    for i in range(len(asym_dm)):
        print "FULL/DEF :  {:1.8f}   {:1.5f}   {:1.5f}".format(asym_fm[i], asym_dm[i], asym_fm[i]/asym_dm[i])

    print asym_dm[-1]
    print asym_fm[-1]

    print sum(dm.sqrt_Fi_Fi_inv)


