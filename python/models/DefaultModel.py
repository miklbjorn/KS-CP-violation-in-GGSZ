from scipy.optimize import minimize
from scipy import stats
import autograd.numpy as np
import cmath, math
from DefaultEfficiency import DefaultEfficiency
from Binning import Binning
from Amplitude import Amplitude
import UtilityFunctions as uf
import os
from autograd import hessian
import autograd 
class DefaultModel(object):
    """Default Model calculates yields using the default setup: ie by calculating Fi and using the yield equations from
    the MI GGSZ LHCb papers.

    It does not include CPV or any global asymmetries.

    Fi and ci, si should be calculated with FullModel.py and passed to DefaultModel.py in the full studies
    However DefaultModel.py can calculate these values for the no-CPV scenario, for testing purposes

    """
    def __init__(self, KS_amplitude_file, efficiency=None, binning_file="input/KsPiPi_optimal.pickle"):
        A_KS, s12, s13 = uf.load_amplitude(KS_amplitude_file)
        self.amplitude = Amplitude(A_KS, s12, s13)
        self.s12 = s12
        self.s13 = s13
        
        bin_def, bin_def_s12, bin_def_s13 = uf.load_binning(binning_file)
        self.binning = Binning(bin_def, bin_def_s12, bin_def_s13,
                s12, s13)
        
        if efficiency==None:
            efficiency = DefaultEfficiency(
                s12, s13)

        self.efficiency = efficiency

        self.Fi, self.ci, self.si = self.calculate_Fi_ci_si()
        self.Fi_inv = np.flip(self.Fi, 0)
        self.sqrt_Fi_Fi_inv = np.sqrt(self.Fi*self.Fi_inv)

        self.bin_num = self.binning.get_number_of_bins()

        self.channel_num = 1 # only fit one Dh channel

        return

    def set_Fi(self, Fi):
        self.Fi = Fi
        self.Fi_inv = np.flip(self.Fi, 0)
        self.sqrt_Fi_Fi_inv = np.sqrt(self.Fi*self.Fi_inv)
        self.bin_num = len(Fi)/2

    def set_ci_si(self, ci,  si):
        self.ci = ci
        self.si = si

    def calculate_Fi_ci_si(self):
        ''' Calculate simple calculation of Fi, ci, si: does not include any CPV effects, not time dependece '''
        bin_num = self.binning.get_number_of_bins()
        Fi = np.array([])
        ci = np.array([])
        si = np.array([])

        A_mag = abs(self.amplitude.get_A(0)) # Just make simple calculation in this class
        A_ph  = np.angle(self.amplitude.get_A(0))
        A_mag_inv = np.transpose(A_mag)
        A_ph_inv  = np.transpose(A_ph)

        
        avg_eff_over_phsp = self.efficiency.get_time_averaged_eff()
        for i in range(-bin_num, bin_num+1):
            if i==0: continue
            bin_idx = self.binning.get_bin_indices(i)
            inv_bin_idx = self.binning.get_bin_indices(-i)
            avg_eff = avg_eff_over_phsp[bin_idx]
            Fi = np.append(Fi, np.sum(
                avg_eff*A_mag[bin_idx]**2))
            ci = np.append(ci, np.sum(
                avg_eff*A_mag[bin_idx]*A_mag_inv[bin_idx]*np.cos(A_ph[bin_idx] - A_ph_inv[bin_idx]) ))
            si = np.append(si, np.sum(
                avg_eff*A_mag[bin_idx]*A_mag_inv[bin_idx]*np.sin(A_ph[bin_idx] - A_ph_inv[bin_idx]) ))

        Fi_inv = np.flip(Fi, 0)
        ci = ci / np.sqrt(Fi*np.flip(Fi, 0))
        si = si / np.sqrt(Fi*np.flip(Fi, 0))
        Fi = Fi/sum(Fi)

        return Fi, ci, si

    def print_Fi_ci_si(self):
        print "Bin  Fi     Fi_inv ci     si"
        for i, f in enumerate(self.Fi):
            print "{: 3d} {: 1.3f} {: 1.3f}  {: 1.3f} {: 1.3f}".format(-8+i + int(i>7), self.Fi[i],self.Fi_inv[i], self.ci[i], self.si[i])


    def fit(self, data, quiet=False, method='BFGS', reserve_method='Powell'):
        if not quiet:
            print "Fitting data using DefaultModel()"

        start_guess_xy = np.array([0.0001, 0.0002, -0.0002, 0.0001]) # offset start guess because that sometimes helps convergence 

        res = minimize(self.yields_chi_square, 
            # [0.01, 0.02, 0.03, 0.04], # offset start guess because that sometimes helps convergence 
            start_guess_xy,
            (data, "string just here to force python to have correct tuple handling"),
            method = method)
        res.method = method

        H = hessian(self.yields_chi_square)
        res = self.add_fit_result_details(res, data, H)

        if not quiet:
            print " Succes:", res.success
            print " Method:", res.method
            print " Msg   :", res.message
            print " Fun   :", res.fun
            print " p-val :", res.p_value

        if not res.success:
            res = minimize(self.yields_chi_square, 
                # [0.01, 0.02, 0.03, 0.04], # offset start guess because that sometimes helps convergence 
                start_guess_xy,
                (data, "string just here to force python to have correct tuple handling"),
                method = reserve_method)
            res.method = reserve_method

            res = self.add_fit_result_details(res, data, H)

            if not quiet:
                print " 2nd try: Succes:", res.success
                print " 2nd try: Method:", res.method
                print " 2nd try: Msg   :", res.message
                print " 2nd try: Fun   :", res.fun
                print " 2nd try: p-val :", res.p_value

        return res

    def add_fit_result_details(self, res, data, Hessian_function):
        ''' adds extra details to the fit result: p value, Hessian from autograd,
        x_err vector, covariance and correlation matrices'''
        
        ndf = len(data)-4
        res.p_value = 1 - stats.chi2.cdf(2*res.fun, ndf)

        # H = hessian(chi2_func)
        res.H = Hessian_function(np.array(res.x), data)

        res.cov_mat = np.linalg.inv(res.H)
        res.cor_mat = uf.get_correlation_matrix(res.cov_mat)
        res.x_unc = np.sqrt(np.diag(res.cov_mat))

        return res



    def yields_chi_square(self, xy, *args):
        ''' The function to me minimized when fitting x, y, is the chi-square
            xy: a list with elements [xm, ym, xp, yp]
            data: The yields to which the fit is being made
        '''
        data = args[0]
        Nplus  = sum(data[0:len(data)/2]) # B+
        Nminus =  sum(data[len(data)/2:]) # B-
        predictions = self.predict_yields(xy, Nplus, Nminus)
        uncertainty_squared = data
        for i, u in enumerate(uncertainty_squared):
            if u==0: uncertainty_squared[i]=1 # make sure there are no divisions by zero
        minLL = np.sum(0.5*(data - predictions)**2/uncertainty_squared) # least squares fit
        return minLL


    def predict_yields(self, xy, Nplus, Nminus, Fi = np.array([0]), normalize=True):
        xm      = xy[0]
        ym      = xy[1]
        xp      = xy[2]
        yp      = xy[3]

        Fi = np.array(Fi)
        Fi_provided = Fi.any()

        if Fi_provided:
            Fi_inv = Fi[-1::-1]
            sqrt_Fi_Fi_inv = np.sqrt(Fi*Fi_inv)
        else:
            Fi = self.Fi
            Fi_inv = self.Fi_inv
            sqrt_Fi_Fi_inv = self.sqrt_Fi_Fi_inv


        # unnormalized distribution cf yield equations
        yp_distribution = (   Fi_inv +(xp**2+yp**2)*Fi      + 2*sqrt_Fi_Fi_inv*(xp*self.ci-yp*self.si))
        ym_distribution = (   Fi     +(xm**2+ym**2)*Fi_inv  + 2*sqrt_Fi_Fi_inv*(xm*self.ci+ym*self.si))
        # predicted yields
        yield_p = Nplus *yp_distribution
        yield_m = Nminus*ym_distribution
        if normalize:
            yield_p = yield_p/sum(yp_distribution)
            yield_m = yield_m/sum(ym_distribution)
        yields = np.append(yield_p, yield_m)

        return yields

    def fit_gamma_from_asymmetry(self, data, rB, dB, gamma0=80, quiet=False):
        if not quiet:
            print "Fitting gamma from asymmetry - fixed (rB, dB)=({:1.4f}, {:3.4f})".format(rB, dB)
        res = minimize(self.asym_chi_square, 
            [uf.deg_to_rad(gamma0)], 
            (data, rB, dB),
            method = 'L-BFGS-B',
            bounds=[(-math.pi/2,math.pi/2)])
        if not quiet:
            print " Succes:", res.success
            print " Msg   :", res.message
            print " Fun   :", res.fun
        return res

    def asym_chi_square(self, gamma, *args):
        ''' The function to me minimized when fitting x, y, is the chi-square
            xy: a list with elements [xm, ym, xp, yp]
            data: The yields to which the fit is being made
        '''
        data = args[0]
        rB = args[1]
        dB = args[2]

        data_asym = (sum(data[16:])-sum(data[:16]))/sum(data)
        pred_asym = self.predict_asym([gamma, rB, dB])
        uncertainty_squared = 0.00000001 #??
        minLL = np.sum(0.5*(data_asym - pred_asym)**2/uncertainty_squared) # least squares fit
        return minLL

    def predict_asym(self, physics_param):
        gamma = physics_param[0]
        rB = physics_param[1]
        dB = physics_param[2]
        xp = rB*math.cos(dB+gamma)
        xm = rB*math.cos(dB-gamma)
        # There is *not* missing a factor 2 on the sums: note that in the ANA-note, the sum is over positive bins only
        pred_asym = sum(self.sqrt_Fi_Fi_inv*self.ci)*(xm-xp)/(1+rB*rB+sum(self.sqrt_Fi_Fi_inv*self.ci)*(xp+xm))
        return pred_asym


if __name__ == "__main__":
    print ("Testing {}".format(os.path.basename(__file__)))
    dm = DefaultModel("../../amplitude_calculations/output/KS_default.pickle",
        binning_file="../input/KsPiPi_optimal.pickle")
    gamma = uf.deg_to_rad(75.)
    r_dk = 0.1
    d_dk = uf.deg_to_rad(130.)
    r_dpi = 0.03
    d_dpi = uf.deg_to_rad(300.)


    import cmath
    xi = (r_dpi/r_dk)*cmath.exp(1j*(d_dpi-d_dk))
    x_xi = xi.real
    y_xi = xi.imag

    xm_dk, ym_dk, xp_dk, yp_dk = uf.get_xy([gamma, r_dk, d_dk])
    xm_dpi, ym_dpi, xp_dpi, yp_dpi = uf.get_xy([gamma, r_dpi, d_dpi])

    N_dk = 6.5e3
    N_dpi = N_dk/0.075 # approximate ratio

    n_dh = dm.predict_yields([xm_dk, ym_dk, xp_dk, yp_dk], N_dk, N_dk)
    res = dm.fit(n_dh)


    import UtilityFunctions as uf
    uf.print_xy_result(res)

    from GammaFitter import GammaFitter
    gf = GammaFitter()
    p, pe, pr = gf.fit_gamma_from_xy_result(res)


