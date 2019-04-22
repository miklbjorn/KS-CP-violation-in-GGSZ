import autograd.numpy as np
from scipy import stats
from scipy.optimize import minimize
import sys, math, cPickle
from collections import defaultdict
import UtilityFunctions as uf
from GammaFitter import GammaFitter

class ToyStudy():
    """ToyStudy runs a standard toy study: given an input FullModel and a fit DefaultModel, make N fits
    and save biases, pulls etc"""
    def __init__(self, 
        fit_model,
        generator_model,
        input_phys_param,
        N_signal, 
        N_toys,
        seed=1):

        self.fit_model = fit_model
        self.generator_model = generator_model
        self.input_phys_param = input_phys_param
        self.N_signal = N_signal 
        self.N_toys = N_toys
        self.seed = seed

        self.saved_outputs = defaultdict(list)
        self.gaus_fits = {}
        self.input_xy = uf.get_xy(input_phys_param)
        self.xy = ["xm", "ym", "xp", "yp"]
        self.phys_param = ["gamma", "rB", "dB"] # physics parameters

        self.gf = GammaFitter()

    def run(self, quiet=False, accept_fails_if_p_more_than=1.0):
        np.random.seed(self.seed)
        self.saved_outputs = defaultdict(list)

        if not quiet: print "Running toy study"
        if not quiet: print "-" # line that get deleted
        self.fit_fails = 0.
        for n in range(0, self.N_toys):
            sys.stdout.write("\033[F") #back to previous line
            sys.stdout.write("\033[K") #clear line
            if not quiet: print ("  Toy: {}/{}".format(n+1, self.N_toys))

            yields = self.generator_model.predict_yields(self.input_phys_param, self.N_signal)
            fluctuated_yields = uf.fluctuate_yields(yields)
            result = self.fit_model.fit(fluctuated_yields, quiet=False)

            if not result.success:
                if not result.p_value > accept_fails_if_p_more_than:
                    self.fit_fails += 1.
                    continue
            
            result_corr = uf.get_correlation_matrix(result.hess_inv)

            for i in range(0, 4):
                self.saved_outputs["{}_val".format(self.xy[i])]  .append( result.x[i])
                self.saved_outputs["{}_bias".format(self.xy[i])] .append( result.x[i] - self.input_xy[i])
                self.saved_outputs["{}_err".format(self.xy[i])]  .append( math.sqrt(result.hess_inv[i, i]))
                self.saved_outputs["{}_pull".format(self.xy[i])] .append( (result.x[i] - self.input_xy[i])/math.sqrt(result.hess_inv[i, i]))
                for j in range(0, 4):
                    self.saved_outputs["{}_{}_corr".format(self.xy[i], self.xy[j])]  .append( result_corr[i, j])


            result_phys_param, result_phys_param_err, result_phys_param = self.gf.fit_gamma_from_xy_result(result)
            result_phys_param_corr = uf.get_correlation_matrix(result_phys_param.hess_inv)
            for i in range(0, 3):
                self.saved_outputs["{}_val".format(self.phys_param[i])]  .append( result_phys_param.x[i])
                self.saved_outputs["{}_bias".format(self.phys_param[i])] .append( result_phys_param.x[i] - self.input_phys_param[i])
                self.saved_outputs["{}_err".format(self.phys_param[i])]  .append( math.sqrt(result_phys_param.hess_inv[i, i]))
                self.saved_outputs["{}_pull".format(self.phys_param[i])] .append( (result_phys_param.x[i] - self.input_phys_param[i])/math.sqrt(result_phys_param.hess_inv[i, i]))
                for j in range(0, 3):
                    self.saved_outputs["{}_{}_corr".format(self.phys_param[i], self.phys_param[j])]  .append( result_phys_param_corr[i, j])

            self.saved_outputs["p_val"].append(1 - stats.chi2.cdf(2*result.fun, 28)) # n dof = 32 yields - 4 param
            self.saved_outputs["fit_chi2"].append(2*result.fun)


        self.fit_bias_and_pull()

        print "  Fit success rate: {:1.2%}  ({:1.0f}/{} failed)".format((self.N_toys-self.fit_fails)/float(self.N_toys), self.fit_fails, self.N_toys)
        return self.saved_outputs["p_val"]

    def fit_bias_and_pull(self):
        for i in range(0, 4):
            self.gaus_fits["{}_val".format(self.xy[i])]  = self.fit_gaussian(self.saved_outputs["{}_val".format(self.xy[i])] )
            self.gaus_fits["{}_bias".format(self.xy[i])] = self.fit_gaussian(self.saved_outputs["{}_bias".format(self.xy[i])] )
            self.gaus_fits["{}_pull".format(self.xy[i])] = self.fit_gaussian(self.saved_outputs["{}_pull".format(self.xy[i])] )
        for i in range(0, 3):
            self.gaus_fits["{}_val".format(self.phys_param[i])]  = self.fit_gaussian(self.saved_outputs["{}_val".format(self.phys_param[i])] )
            self.gaus_fits["{}_bias".format(self.phys_param[i])] = self.fit_gaussian(self.saved_outputs["{}_bias".format(self.phys_param[i])] )
            self.gaus_fits["{}_pull".format(self.phys_param[i])] = self.fit_gaussian(self.saved_outputs["{}_pull".format(self.phys_param[i])] )

    def save_results(self, file_name):
        file = open(file_name, "wb")
        cPickle.dump([
            self.saved_outputs, 
            self.fit_fails,
            self.gaus_fits
            ], file, protocol=2)
        file.close()

    def load_results(self, file_name):
        file = open(file_name, "r")
        list_to_load = cPickle.load(file)
        file.close()
        self.saved_outputs = list_to_load[0]
        self.fit_fails = list_to_load[1]
        self.gaus_fits = list_to_load[2]

    def gaussian_log_likelihood(self, p, *args):
        data = args[0]
        minus_ll = 0.5 * sum( (p[0]-data)**2 / p[1]**2  + 2*math.log(abs(p[1])) )
        return minus_ll

    def fit_gaussian(self, data):
        data = np.array(data)
        res = minimize(self.gaussian_log_likelihood, 
            [data.mean()+0.001, data.std()+0.001], # slight offset avoids exact zero jacobian
            (data, "string just here to force python to have correct tuple handling"),
            method = 'BFGS')
        mu = res.x[0]
        sigma = res.x[1]
        mu_err = math.sqrt(res.hess_inv[0,0])
        sigma_err = math.sqrt(res.hess_inv[1,1])
        return mu, sigma, mu_err, sigma_err

    def get_list(self, name):
        if self.saved_outputs.has_key(name):
            return self.saved_outputs[name]
        else:
            print "Unknown key:", name
            return None

    def get_gaus(self, name):
        if self.gaus_fits.has_key(name):
            return self.gaus_fits[name]
        else:
            print "Unknown key:", name
            return None


if __name__=="__main__":
    from models import DefaultModel
    from models import FullModel
    from PhysicalParameters import PhysicalParameters
    dm = DefaultModel("../amplitude_calculations/output/KS_default.pickle")
    fm = FullModel(
        "../amplitude_calculations/output/KS_default.pickle",
        "../amplitude_calculations/output/KL_default.pickle",
        0,
        PhysicalParameters(False))
    ts = ToyStudy(dm, fm, [uf.deg_to_rad(85), 0.1, uf.deg_to_rad(100)], 4000, 10)
    
    test = np.random.normal(0, 1, size=1000)
    print test.mean()
    print test.std()
    # print test
    print ts.fit_gaussian(test)

    print ts.run()
    ts.seed = 3
    print ts.run()
    ts.seed = 1
    print ts.run()
