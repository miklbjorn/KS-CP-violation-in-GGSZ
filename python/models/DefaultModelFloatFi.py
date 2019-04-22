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

class DefaultModelFloatFi(DefaultModel):
    """DefaultModelFloatFi uses same model as DefaultModel(), ie. it calculates yields 
    using the default setup: ie by calculating Fi and using the yield equations from
    the MI GGSZ LHCb papers, and does not include CPV or any global asymmetries.

    However it ALSO floats the Fi parameters (so far unconstrained) by overloading the
    chi-square and fit functions

    Fi and ci, si should be calculated with FullModel.py and passed to DefaultModel.py in the full studies
    However DefaultModel.py can calculate these values for the no-CPV scenario

    """
    def __init__(self, KS_amplitude_file, efficiency=None, binning_file="input/KsPiPi_optimal.pickle"):
        super(DefaultModelFloatFi, self).__init__(
            KS_amplitude_file, efficiency, binning_file)


    def yields_chi_square(self, xy_and_Fi, *args):
        ''' The function to me minimized when fitting x, y, is the chi-square
            xy_Fi: a list with elements [{Fi}, xm, ym, xp, yp]
            data: The yields to which the fit is being made
        '''
        data = args[0]
        Nplus  = sum(data[0:len(data)/2]) # B+
        Nminus =  sum(data[len(data)/2:]) # B-
        xy = xy_and_Fi[0:4]
        Fi = xy_and_Fi[4:]
        Fi = np.append(Fi, 1-sum(Fi)) 


        predictions = self.predict_yields(xy, Nplus, Nminus, Fi=Fi)
        uncertainty_squared = data
        for i, u in enumerate(uncertainty_squared):
            if u==0: uncertainty_squared[i]=1 # make sure there are no divisions by zero
        minLL = np.sum(0.5*(data - predictions)**2/uncertainty_squared) # least squares fit
        return minLL


    def fit(self, data, quiet=False):
        if not quiet:
            print "Fitting data using DefaultModelFloatFi()"

        start_guess_xy = np.array([0.0001, 0.0002, -0.0002, 0.0001]) # offset start guess because that sometimes helps convergence 
        start_Fi = self.Fi[:-1] # start with all Fi equal to 1 (normalisation does not matter, handled elsewhere)

        start_guess = np.append(start_guess_xy, start_Fi)

        bounds = []
        for i, x in enumerate(start_guess):
            if i < 4 : 
                bounds.append((None, None))
            else :
                bounds.append((0, 1))

        res = minimize(self.yields_chi_square, 
            start_guess,
            (data, "string just here to force python to have correct tuple handling"),
            method = 'L-BFGS-B',
            bounds=bounds)

        H = hessian(self.yields_chi_square)
        res = self.add_fit_result_details(res, data, H)

        if not quiet:
            print " Succes:", res.success
            print " Msg   :", res.message
            print " Fun   :", res.fun
            print " p-val :", res.p_value

        return res


if __name__ == "__main__":

    dm0 = DefaultModel("../amplitude_calculations/output/KS_default.pickle")

    print (issubclass(DefaultModel, object))

    print ("Testing {}".format(os.path.basename(__file__)))
    dm = DefaultModelFloatFi("../amplitude_calculations/output/KS_default.pickle")
    dm.print_Fi_ci_si()
    dm.lolzor()

