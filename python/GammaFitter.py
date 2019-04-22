import UtilityFunctions as uf
from scipy.optimize import minimize
from scipy import stats
import autograd.numpy as np
from models import DefaultModel
import math, os
from autograd import hessian

class GammaFitter(object):

    def __init__(self):
        self.has_made_fit = False



    def add_fit_result_details(self, res, data, H):
        
        # add the same results as in the yield fits
        ndf = len(data)-4
        res.p_value = 1 - stats.chi2.cdf(2*res.fun, ndf)

        res.H = H

        res.cov_mat = np.linalg.inv(res.H)
        try:
            res.cor_mat = uf.get_correlation_matrix(res.cov_mat)
        except ValueError:
            # something went wrong when inverting
            # negative elements in covariance matrix diag
            # don't trust results
            res.success = False
            res.cor_mat = 0 * H
        res.x_unc = np.sqrt(np.diag(res.cov_mat))

        # translate to degrees; 1st save raw results
        res.raw_x = np.copy(res.x)
        res.raw_x_unc = np.copy(res.x_unc)

        # ensure angles are positive
        while res.x[0] < 0: res.x[0] += 2*math.pi
        while res.x[2] < 0: res.x[2] += 2*math.pi

        for i in range(len(res.x)):
            if (i % 2) == 0:
                # every other result is an angle, so apply rad_to_deg
                res.x[i] = uf.rad_to_deg(res.x[i])
                res.x_unc[i] = uf.rad_to_deg(res.x_unc[i])

        return res

    ######################################################
    #
    # Functions for fit from xy
    #
    ######################################################

    def fit_from_result(self, fit_result, start_guess = [uf.deg_to_rad(75.), 0.1, uf.deg_to_rad(130.)]):
        xy_vector = np.array(fit_result.x[0:4]) # xm ym xp yp
        xy_cov_mat = np.array(fit_result.cov_mat[0:4,0:4]) # covariance matrix
        return self.fit(xy_vector, xy_cov_mat, start_guess)

    def fit(self, xy_vector, xy_cov_mat, start_guess = [uf.deg_to_rad(75.), 0.1, uf.deg_to_rad(130.)]):

        xy_cov_mat_inv = np.linalg.inv(xy_cov_mat) # covariance matrix

        # bound r and angles to be positive
        bounds = [(0, math.pi), (1e-5, None), (None, 4*math.pi)]

        res = minimize(self.chi_square_func_for_xy, 
            start_guess, 
            (xy_vector, xy_cov_mat_inv),
            method = 'L-BFGS-B',
            bounds=bounds)

        H_func = hessian(self.chi_square_func_for_xy)

        H = H_func(np.array(res.x), xy_vector, xy_cov_mat_inv)
        res = self.add_fit_result_details(res, xy_vector, H)

        return res.x, res.x_unc, res
        return res.x, 1, res

    def chi_square_func_for_xy(self, physics_param, *args):
        ''' chi square for physics param = [xm, ym, xp, yp]'''
        xy_vector = args[0]
        xy_cov_mat_inv = args[1]
        fit_xy_vector = np.array(uf.get_xy(physics_param))
        dx = xy_vector - fit_xy_vector
        chi_square = 0.5*np.dot(
            np.dot(dx.transpose(),xy_cov_mat_inv),
            dx)
        return chi_square



    ######################################################
    #
    # Functions for fit from global aymmetry
    #
    ######################################################

    def fit_gamma_from_global_asym(self, data, eps, rB, dB, sum_Ki_Ki_inv_ci, parameters):
        asym_data = (sum(data[16:])-sum(data[:16]))/sum(data)
        Gamma = 0.5*(parameters.g_S + parameters.g_L)
        mu = parameters.dm/Gamma
        eps_term = 2*eps.real-2*parameters.g_S/(Gamma*(1+mu**2))*(eps.real+mu*eps.imag)
        def chi2(gamma):
            denom = 1 * rB**2+2*sum_Ki_Ki_inv_ci*math.cos(dB)*math.sin(gamma)
            asym_pred =  (2*sum_Ki_Ki_inv_ci*rB*math.sin(dB)*math.sin(gamma)-eps_term)/denom
            chi_square = 0.5*(asym_data-asym_pred)**2
            return chi_square
        result = minimize(chi2, [uf.deg_to_rad(70)])
        return result


if __name__ == "__main__":
    print ("Testing {}".format(os.path.basename(__file__)))
    test_xy_param = uf.get_xy([uf.deg_to_rad(88), 0.09, uf.deg_to_rad(100)])
    gf = GammaFitter()
    dm = DefaultModel("../amplitude_calculations/output/KS_default.pickle")
    data = dm.predict_yields(test_xy_param, Nplus=2000, Nminus=2000)
    result = dm.fit(data)
    print "Fitted xy"
    for i in range(0, 4):
        print "{: 2.3f} +/- {: 2.3f}".format(result.x[i], math.sqrt(result.hess_inv[i,i]))
    
    print "Fitted xy correlations"
    print uf.get_correlation_matrix(result.hess_inv[0:4,0:4])

    param, err, res = gf.fit_gamma_from_xy(result)
    print "Fitted physics param"
    for i in range(0, 3):
        print "{: 2.3f} +/- {: 2.3f}".format(param[i], err[i])