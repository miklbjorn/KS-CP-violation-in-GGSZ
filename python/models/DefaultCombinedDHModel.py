from scipy.optimize import minimize
from scipy import stats
import autograd.numpy as np
import cmath, math
from DefaultEfficiency import DefaultEfficiency
from Binning import Binning
from Amplitude import Amplitude
import UtilityFunctions as uf
import os
from models import DefaultModel
from autograd import hessian

class DefaultCombinedDHModel(DefaultModel):
    """DefaultModelFloatFi uses same model as DefaultModel(), ie. it calculates yields 
    using the default setup: ie by calculating Fi and using the yield equations from
    the MI GGSZ LHCb papers, and does not include CPV or any global asymmetries.

    However it ALSO floats the Fi parameters (so far unconstrained) by overloading the
    chi-square and fit functions

    Fi and ci, si should be calculated with FullModel.py and passed to DefaultModel.py in the full studies
    However DefaultModel.py can calculate these values for the no-CPV scenario

    """
    def __init__(self, KS_amplitude_file, efficiency=None, binning_file="input/KsPiPi_optimal.pickle"):
        super(DefaultCombinedDHModel, self).__init__(
            KS_amplitude_file, efficiency, binning_file)

        self.channel_num = 2 # fits two Dh channels


    def fit(self, data, quiet=False):
        if not quiet:
            print "Fitting data using DefaultCombinedDHModel()x"

        start_guess_xy_xi = np.array([0.0001, 0.0002, -0.0002, 0.0001, 0.0001, -0.0001]) # offset start guess because that sometimes helps convergence 

        res = minimize(self.yields_chi_square, 
            # [0.01, 0.02, 0.03, 0.04], # offset start guess because that sometimes helps convergence 
            start_guess_xy_xi,
            (data, "string just here to force python to have correct tuple handling"),
            method = 'BFGS',
            tol=1e-4)

        H = hessian(self.yields_chi_square)
        res = self.add_fit_result_details(res, data, H)


        if not quiet:
            print " Succes      :", res.success
            print " Msg         :", res.message
            print " Fun         :", res.fun
            print " p-val       :", res.p_value

        return res

    def yields_chi_square(self, xy_xi, *args):
        ''' The function to me minimized when fitting x, y, is the chi-square
            xy: a list with elements [xm, ym, xp, yp]
            data: The yields to which the fit is being made
        '''
        data = args[0]
        N_ch = len(data)/4 # yields per charge in a single channel
        Nplus     = sum(data[0:N_ch]) # B+
        Nminus    =  sum(data[N_ch:2*N_ch]) # B-
        Nplus_xi  = sum(data[2*N_ch:3*N_ch]) # B+, xi channel
        Nminus_xi =  sum(data[3*N_ch:]) # B-, xi channel
        predictions = self.predict_yields(xy_xi, Nplus, Nminus, Nplus_xi, Nminus_xi)
        uncertainty_squared = data
        for i, u in enumerate(uncertainty_squared):
            if u==0: uncertainty_squared[i]=1 # make sure there are no divisions by zero
        minLL = np.sum(0.5*(data - predictions)**2/uncertainty_squared) # least squares fit
        return minLL


    def predict_yields(self, xy_xi, Nplus, Nminus, Nplus_xi, Nminus_xi, Fi = np.array([0]), normalize=True):
        xm      = xy_xi[0]
        ym      = xy_xi[1]
        xp      = xy_xi[2]
        yp      = xy_xi[3]
        x_xi    = xy_xi[4]
        y_xi    = xy_xi[5]

        # calculate xi channel x, y parameters
        xm_xi = x_xi*xm - y_xi*ym;
        ym_xi = x_xi*ym + y_xi*xm;
        xp_xi = x_xi*xp - y_xi*yp;
        yp_xi = x_xi*yp + y_xi*xp;

        yields_base = super(DefaultCombinedDHModel, self).predict_yields(
            [xm, ym, xp, yp], Nplus, Nminus, Fi, normalize)
        
        yields_xi = super(DefaultCombinedDHModel, self).predict_yields(
            [xm_xi, ym_xi, xp_xi, yp_xi], Nplus_xi, Nminus_xi, Fi, normalize)

        yields = np.append(yields_base, yields_xi)
        return yields




if __name__ == "__main__":
    print ("Testing {}".format(os.path.basename(__file__)))
    dcm = DefaultCombinedDHModel("../../amplitude_calculations/output/KS_default.pickle",
        binning_file="../input/KsPiPi_optimal.pickle")
    gamma = uf.deg_to_rad(75.)
    r_dk = 0.1
    d_dk = uf.deg_to_rad(130.)
    r_dpi = 0.005
    d_dpi = uf.deg_to_rad(300.)


    import cmath
    xi = (r_dpi/r_dk)*cmath.exp(1j*(d_dpi-d_dk))
    x_xi = xi.real
    y_xi = xi.imag

    xm_dk, ym_dk, xp_dk, yp_dk = uf.get_xy([gamma, r_dk, d_dk])
    xm_dpi, ym_dpi, xp_dpi, yp_dpi = uf.get_xy([gamma, r_dpi, d_dpi])

    N_dk = 6.5e3
    N_dpi = N_dk/0.075 # approximate ratio

    n_dh = dcm.predict_yields([xm_dk, ym_dk, xp_dk, yp_dk, x_xi, y_xi], N_dk, N_dk, N_dpi, N_dpi)
    res = dcm.fit(n_dh)

    import UtilityFunctions as uf
    uf.print_xy_xi_result(res)

    from GammaFitter import GammaFitter
    gf = GammaFitter()
    p, pe, pr = gf.fit_gamma_from_xy_xi_result(res)





