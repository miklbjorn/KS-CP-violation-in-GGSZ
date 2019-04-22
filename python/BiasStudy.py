#coding=UTF-8

import UtilityFunctions as uf
import autograd.numpy as np
import os, sys, cPickle
from collections import defaultdict
from GammaFitter import GammaFitter

class BiasStudy():
    """docstring for BiasStudy

    - can use a list of KL amplitudes, and report 
    - can vary input physics parameters
    - can use both "fit avg yields" stup and proper "toy studies" setup
    - can use both default gamma etc, but also global yields based

    - save all biases, mean bias etc


    """
    def __init__(self, 
        generator_model=None, # Model used to generate yields (with efficiency/mat interactions set, etc)
        param_sets=[],# A Nx3 or Nx5 with physics param to use, angles in RAD
        N_signal= [14000, 14000/0.075], 
        KL_amplitude_files = [], # Can specify more than one KL-amplitude file
        empty=False # make an empty study, for loading
        ):
        
        self.no_printed_warning = True

        if empty:
            return

        if generator_model==None:
            raise ValueError("You MUST supply generator model to non-empty BiasStudy")
        if param_sets==[]:
            raise ValueError("You MUST supply param_sets to non-empty BiasStudy")

        self.generator_model = generator_model
        self.param_sets = param_sets

        self.param_sets_deg = []
        for p in param_sets:
            p_deg = []
            for i in range(len(p)):
                if i % 2 == 0:
                    p_deg.append(uf.rad_to_deg(p[i]))
                else:
                    p_deg.append(p[i])
            self.param_sets_deg.append(p_deg)

        self.yields = {}

        self.N_signal = N_signal
        self.KL_amplitude_files = KL_amplitude_files 

        if KL_amplitude_files == []:
            self.KL_amplitude_files = [generator_model.KL_amplitude_file]

        self.yields_init = False

        self.results = {}



    def init_yields(self, quiet=False):

        nf = len(self.KL_amplitude_files)
        npar = len(self.param_sets)
        total_yields= nf*npar

        if not quiet:
            print "Making yields for bias study:"
            print "  # of params  : {}".format(npar)
            print "  # of A(KL)'s : {}".format(nf)
            print "  # TOTAL yield sets : {}".format(total_yields)

        yield_n = 0
        if not quiet: print ("  Yield gen: {}/{}".format(yield_n, total_yields))
        for fi, f in enumerate(self.KL_amplitude_files):
            self.generator_model.load_KL_amplitude(f)
            for pi, p in enumerate(self.param_sets):
                physics_param = p

                sys.stdout.write("\033[F") #back to previous line
                sys.stdout.write("\033[K") #clear line
                yield_n += 1
                if not quiet: print ("  Yield gen: {}/{}".format(yield_n, total_yields))

                if len(p) == 3:
                    # single channel setup
                    avg_yields = self.generator_model.predict_yields(p, self.N_signal[0]) # Yields calculated using FullModel.py
                elif len(p) == 5:
                    # two channel channel setup
                    avg_yields_ch1 = self.generator_model.predict_yields(p[0:3], self.N_signal[0]) # Yields calculated using FullModel.py
                    avg_yields_ch2 = self.generator_model.predict_yields([p[0], p[3], p[4]], self.N_signal[1]) # Yields calculated using FullModel.py
                    avg_yields = np.append(avg_yields_ch1, avg_yields_ch2)

                self.yields["{}_{}".format(fi, pi)] = avg_yields

        self.yields_init = True



    def run(self, fit_model, gamma_fitter, name="default", quiet=False):

        if not self.yields_init:
            self.init_yields()


        self.fit_observables(  fit_model,    fit_model.channel_num, name, quiet)
        self.fit_physics_param(gamma_fitter, fit_model.channel_num, name, quiet)


    def fit_observables(self, fit_model, channel_num, name="default", quiet=False):

        nf = len(self.KL_amplitude_files)
        npar = len(self.param_sets)
        total_fits= nf*npar
        n_xy = 6
        if channel_num == 1:
            n_xy=4    

        self.results["{}-fit-success".format(name)] = np.zeros((nf, npar))
        self.results["{}-xy-val".format(name)] = np.zeros((nf, npar, n_xy))
        self.results["{}-xy-bias".format(name)] = np.zeros((nf, npar, n_xy))
        self.results["{}-xy-res".format(name)] = {}

        fit_n = 0
        print "Fitting yields for observables for :", name
        print " Fit:", fit_n
        for fi, f in enumerate(self.KL_amplitude_files):
            for pi, p in enumerate(self.param_sets):
                
                if channel_num == 1:
                    input_xy = uf.get_xy(p[0:3])
                    avg_yields = self.yields["{}_{}".format(fi, pi)][0:32]
                else:
                    if len(p)<5:
                        raise ValueError("Cannot fit with 2 channel Model, only generated 1 channel yields")
                    input_xy = uf.get_xy_xi(p)
                    avg_yields = self.yields["{}_{}".format(fi, pi)]
                    
                sys.stdout.write("\033[F\033[K") #back to previous line, clear line
                fit_n += 1
                if not quiet: print ("  Fit: {}/{}".format(fit_n, total_fits))

                 # Yields calculated using FullModel.py
                direct_result = fit_model.fit(avg_yields, quiet=True)
                self.results["{}-xy-res".format(name)]["{}_{}".format(fi, pi)] = direct_result
                self.results["{}-xy-val".format(name)][fi, pi, :] = direct_result.x[0:n_xy]
                self.results["{}-xy-bias".format(name)][fi, pi, :] = np.array(direct_result.x[0:n_xy]) - np.array(input_xy)
                self.results["{}-fit-success".format(name)][fi, pi] = direct_result.success
                
        return

    def fit_physics_param(self, gamma_fitter, channel_num, name="default", quiet=False):

        nf = len(self.KL_amplitude_files)
        npar = len(self.param_sets)
        total_fits= nf*npar

        n_phys = 5
        if channel_num == 1:
            n_phys=3


        self.results["{}-phys-fit-success".format(name)] = np.zeros((nf, npar))
        self.results["{}-physics_param-val".format(name)] = np.zeros((nf, npar, n_phys))
        self.results["{}-physics_param-bias".format(name)] = np.zeros((nf, npar, n_phys))
        self.results["{}-phys-fit-res".format(name)] = {}

        fit_n = 0
        print "Fitting observables for physics_param:", name
        print " Fit:", fit_n
        for fi, f in enumerate(self.KL_amplitude_files):
            for pi, p in enumerate(self.param_sets):
                physics_param_deg = self.param_sets_deg[pi]
                
                if channel_num == 1:
                    physics_param_deg = physics_param_deg[0:3]
                    
                direct_result = self.results["{}-xy-res".format(name)]["{}_{}".format(fi, pi)]

                sys.stdout.write("\033[F\033[K") #back to previous line, clear line
                fit_n += 1
                if not quiet: print ("  Fit: {}/{}".format(fit_n, total_fits))

                # start guess close to actual values for better convergence
                phys_fit = gamma_fitter.fit_from_result(direct_result, start_guess = np.array(p)*0.95)
                # the fit is a bit sensitive to initial parameters so try a default set, as well as the one close to inputs
                # use it if it converged and the solution is better
                phys_fit_default_start = gamma_fitter.fit_from_result(direct_result)
                if phys_fit_default_start[2].fun < phys_fit[2].fun and phys_fit_default_start[2].success:
                    phys_fit = phys_fit_default_start

                direct_physics_param, direct_physics_param_err, direct_physics_param_result = phys_fit
                self.results["{}-physics_param-val".format(name)][fi, pi, :] = direct_physics_param
                self.results["{}-physics_param-bias".format(name)][fi, pi, :] = np.array(direct_physics_param) - np.array(physics_param_deg)
                self.results["{}-phys-fit-success".format(name)][fi, pi] = direct_physics_param_result.success
                        # Add to the results 

    def get_results(self, var_name, result_name="default", result="bias", KL_results="all"):
        xy_vars = ["xm", "ym", "xp", "yp", "x_xi", "y_xi"]
        physics_param_vars = ["gamma", "r", "d", "r_xi", "d_xi"]
        if (not (var_name in xy_vars)) and (not (var_name in physics_param_vars)):
            raise ValueError("'{}' not an xy parameter (xm, ym, xp, yp, x_xi, y_xi) or physics parameter (gamma, r, d, r_xi, d_xi)".format(var_name))
        
        # Pick out the right parameters, and choose fit values or bais
        if var_name in xy_vars:
            idx = xy_vars.index(var_name)
            var_type = "xy"
        else:
            idx = physics_param_vars .index(var_name)
            var_type = "physics_param"
        
        # select results where both fits succeeded
        array = self.results["{}-{}-{}".format(result_name, var_type, result)][:,:,idx]
        success = self.results["{}-fit-success".format(result_name)]
        phys_success = self.results["{}-phys-fit-success".format(result_name)]
        weights = np.minimum(success, phys_success)

        avg_weights = np.average(weights,axis=0)
        avg_obs_suc = np.average(success,axis=0)
        avg_phy_suc = np.average(phys_success,axis=0)

        for i, aw in enumerate(avg_weights):
            if aw == 0:
                weights[1, i] =1
                if self.no_printed_warning:
                    print ("NO FITS converged for param set number '", i, "' (zero indexed)",
                    " name:", result_name)
                    print "setting 1 weight to 1 (but fix!"
                    print "weight   :", avg_weights
                    print "obs fits :", avg_obs_suc
                    print "phys fits:", avg_phy_suc
                    self.no_printed_warning = False


        # The determine what should be done with the different KL amplitudes
        if KL_results=="all":
            # Return a (# KL amplitudes) x (#input params) array with all results
            return array
        elif KL_results=="weights":
            return weights
        elif KL_results=="avg":
            return np.average(array, axis=0, weights=weights)
        elif KL_results=="std":
            # calculate weigthed std (not in numpy)
            avg = np.average(array, axis=0, weights=weights)
            var = np.average((array-avg)**2, axis=0, weights=weights)
            return np.sqrt(var)

    def save_params(self, file_name):
        file = open(file_name, "wb")
        cPickle.dump([
            self.param_sets, self.param_sets_deg
            ], file, protocol=2)
        file.close()

    def load_params(self, file_name):
        file = open(file_name, "r")
        list_to_load = cPickle.load(file)
        file.close()
        self.param_sets = list_to_load[0]
        self.param_sets_deg = list_to_load[1]

    def save_results(self, file_name):
        file = open(file_name, "wb")
        cPickle.dump([
            self.results
            ], file, protocol=2)
        file.close()

    def load_results(self, file_name):
        file = open(file_name, "r")
        list_to_load = cPickle.load(file)
        file.close()
        self.results = list_to_load[0]

    def save_yields(self, file_name):
        file = open(file_name, "wb")
        cPickle.dump([
            self.yields
            ], file, protocol=2)
        file.close()

    def load_yields(self, file_name):
        file = open(file_name, "r")
        list_to_load = cPickle.load(file)
        file.close()
        self.yields = list_to_load[0]
        self.yields_init = True


if __name__=="__main__":
    print ("Testing {}".format(os.path.basename(__file__)))
    from models import DefaultModel, DefaultCombinedDHModel
    from models import FullModel
    from PhysicalParameters import PhysicalParameters
    from GammaFitter import GammaFitter
    from GammaFitterXiSetup import GammaFitterXiSetup

    dm = DefaultModel("../amplitude_calculations/output/KS_default.pickle")
    dcm = DefaultCombinedDHModel("../amplitude_calculations/output/KS_default.pickle")
    fm = FullModel(
        "../amplitude_calculations/output/KS_default.pickle",
        "../amplitude_calculations/output/KL_default.pickle",
        PhysicalParameters().default_eps,
        PhysicalParameters(True))
    bs = BiasStudy(fm, 
        param_sets = np.array([[uf.deg_to_rad(x), 0.1, uf.deg_to_rad(130)] for x in [ 80]]),
        KL_amplitude_files = ["../amplitude_calculations/output/KL_seed_1.pickle", "../amplitude_calculations/output/KL_seed_2.pickle"])
    

    gf = GammaFitter()
    gf_xi = GammaFitterXiSetup()

    bs.init_yields()
    # bs.load_yields("tmp_n.pickle")
    bs.run(dm, gf)
    bs.run(dcm, gf_xi, "com")


    # print bs.get_results("xm", KL_results="all")
    # print bs.get_results("xm", KL_results="avg")
    # print bs.get_results("xm", KL_results="std")

    # # print bs.results["xy-val"]
    # print bs.results["xy-bias"]
