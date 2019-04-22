import UtilityFunctions as uf
from scipy.optimize import minimize
from scipy import stats
import autograd.numpy as np
from models import DefaultCombinedDHModel
import math, os
from autograd import hessian
from GammaFitter import GammaFitter

class GammaFitterXiSetup(GammaFitter):

    def __init__(self):
        self.has_made_fit = False



    def fit_from_result(self, fit_result):
        xy_xi_vector = np.array(fit_result.x[0:6]) # xm ym xp yp
        xy_xi_cov_mat = np.array(fit_result.cov_mat[0:6,0:6]) # covariance matrix
        return self.fit(xy_xi_vector, xy_xi_cov_mat)

    def fit(self, xy_xi_vector, xy_xi_cov_mat):

        xy_xi_cov_mat_inv = np.linalg.inv(xy_xi_cov_mat) # covariance matrix

        # bound r's to be positive
        bounds = [(0, math.pi), 
            (1e-5, None), (None, None),
            (1e-5, None), (None, None)]

        res = minimize(self.chi_square_func_for_xy_xi, 
            [uf.deg_to_rad(70), 0.1, uf.deg_to_rad(130), 0.005, uf.deg_to_rad(300)], 
            (xy_xi_vector, xy_xi_cov_mat_inv),
            method = 'L-BFGS-B',
            bounds=bounds)

        H_func = hessian(self.chi_square_func_for_xy_xi)
        H = H_func(np.array(res.x), xy_xi_vector, xy_xi_cov_mat_inv)

        res = self.add_fit_result_details(res, xy_xi_vector,  H)

        return res.x, res.x_unc, res

    def chi_square_func_for_xy_xi(self, physics_param, *args):
        ''' chi square for physics param = [xm, ym, xp, yp, x_xi, y_xi]'''
        xy_xi_vector = args[0]
        xy_xi_cov_mat_inv = args[1]
        fit_xy_xi_vector = np.array(uf.get_xy_xi(physics_param))
        dx = xy_xi_vector - fit_xy_xi_vector
        chi_square = 0.5*np.dot(
            np.dot(dx.transpose(),xy_xi_cov_mat_inv),
            dx)
        return chi_square    

if __name__ == "__main__":
    print ("Testing {}".format(os.path.basename(__file__)))
    test_xy_param = uf.get_xy_xi([uf.deg_to_rad(75), 0.09, uf.deg_to_rad(130), 0.005, uf.deg_to_rad(300)])
    gf = GammaFitterXiSetup()
    dm = DefaultCombinedDHModel("../amplitude_calculations/output/KS_default.pickle")
    data = dm.predict_yields(test_xy_param, Nplus=7000, Nminus=7000, Nplus_xi=20000, Nminus_xi=20000)
    result = dm.fit(data)


    param, err, res = gf.fit_from_result(result)
    print "Fitted physics param"
    for i in range(0, 5):
        print "{: 2.3f} +/- {: 2.3f}".format(param[i], err[i])