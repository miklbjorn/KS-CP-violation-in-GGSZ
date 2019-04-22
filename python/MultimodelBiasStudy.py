#coding=UTF-8

import UtilityFunctions as uf
import autograd.numpy as np
import pandas as pd
import os, sys, cPickle
from collections import defaultdict
from GammaFitter import GammaFitter

class MultimodelBiasStudy():
    """docstring for BiasStudy

    - can use a list of KL amplitudes, and report 
    - can vary input physics parameters
    - can use both "fit avg yields" stup and proper "toy studies" setup
    - can use both default gamma etc, but also global yields based

    - save all biases, mean bias etc


    """
    def __init__(self, 
        generator_models=None, # Model used to generate yields (with efficiency/mat interactions set, etc)
        param_sets=None,# A Nx3 or Nx5 with physics param to use, angles in RAD
        N_signal= [14000, 14000/0.075], 
        KL_amplitude_files = [], # Can specify more than one KL-amplitude file
        empty=False):
        
        if empty:
            return
        if generator_models==None:
            raise ValueError("You MUST supply list od generator models to non-empty BiasStudy")
        if param_sets==[]:
            raise ValueError("You MUST supply param_sets to non-empty BiasStudy")


        self.no_printed_warning = True

        self.generator_models = generator_models
        self.param_sets = param_sets

        self.param_sets_deg = []

        p_dicts = []
        for p in param_sets:
            p_deg = []
            p_dict = {}
            p_dict['g_val'] = p[0]
            p_dict['g_val_deg'] = uf.rad_to_deg(p[0])
            p_dict['r_val'] = p[1]
            p_dict['d_val'] = p[2]
            p_dict['d_val_deg'] = uf.rad_to_deg(p[2])
            p_dicts.append(p_dict)

        self.param = pd.DataFrame(p_dicts)

        self.yields = pd.DataFrame(columns=[
            'param_index', 'gen_model_index', 'KL_amp_index', 'yield_set'])

        self.fits = pd.DataFrame(columns=[
            'param_index', 'gen_model_index', 'KL_amp_index', 'type', 'result', 'success'])

        self.results = pd.DataFrame(columns=[
            'param_index', 'gen_model_index', 'KL_amp_index', 'param', 'val', 'bias', 'uncertainty'])

        self.N_signal = N_signal
        self.KL_amplitude_files = KL_amplitude_files 

        if KL_amplitude_files == []:
            self.KL_amplitude_files = [generator_model.KL_amplitude_file]

        self.yields_init = False





    def init_yields(self, quiet=False):
        nm = len(self.generator_models)
        nf = len(self.KL_amplitude_files)
        npar = len(self.param_sets)
        total_yields= nf*npar*nm

        if not quiet:
            print "Making yields for bias study:"
            print "  # of models  : {}".format(nm)
            print "  # of params  : {}".format(npar)
            print "  # of A(KL)'s : {}".format(nf)
            print "  # TOTAL yield sets : {}".format(total_yields)

        yield_n = 0
        yield_dfs = []
        if not quiet: print ("  Yield gen: {}/{}".format(yield_n, total_yields))
        for mi, m in enumerate(self.generator_models):
            for fi, f in enumerate(self.KL_amplitude_files):
                m.load_KL_amplitude(f)
                for param_index, param_row in self.param.iterrows():
                    physics_param = [param_row.g_val, param_row.r_val, param_row.d_val]

                    sys.stdout.write("\033[F") #back to previous line
                    sys.stdout.write("\033[K") #clear line
                    yield_n += 1
                    if not quiet: print ("  Yield gen: {}/{}".format(yield_n, total_yields))

                    
                    avg_yields = m.predict_yields(physics_param, float(self.N_signal[0])/nm) # Yields calculated using FullModel.py
                    yield_dfs.append(pd.DataFrame([[
                        param_index, mi, fi, avg_yields]],
                        columns=self.yields.columns
                        ))
                    

        self.yields = pd.concat(yield_dfs)
        self.yields_init = True



    def run(self, fit_models, avg_fit_model, gamma_fitter, name="default", quiet=False):

        if not self.yields_init:
            self.init_yields()


        self.fit_observables(  fit_models, avg_fit_model, avg_fit_model.channel_num, name, quiet)
        self.fit_physics_param(gamma_fitter, avg_fit_model.channel_num, name, quiet)


    def fit_observables(self, fit_models, avg_fit_model, channel_num, name="default", quiet=False):

        nm = len(self.generator_models)
        nf = len(self.KL_amplitude_files)
        npar = len(self.param_sets)
        total_fits= nf*npar
        n_xy = 6
        if channel_num == 1:
            n_xy=4    

   
        fit_n = 0
        print "Fitting yields for observables for :", name
        print " Fit:", fit_n



        for fi, f in enumerate(self.KL_amplitude_files):
            for param_index, param_row in self.param.iterrows():
                input_xy = uf.get_xy([param_row.g_val, param_row.r_val, param_row.d_val])
                
                total_yields = np.zeros(32)

                def add_result(model_index, direct_result):
                    self.fits = self.fits.append({'param_index': param_index, 'gen_model_index': model_index, 'KL_amp_index': fi,
                        'type': 'observable_fit', 'result': direct_result, 'success': direct_result.success}, ignore_index=True)
                    
                    xy = np.array(direct_result.x[0:n_xy])
                    xy_unc = np.array(direct_result.x_unc[0:n_xy])
                    bias = xy - np.array(input_xy)
                    xy_names = ['xm', 'ym', 'xp', 'yp']
                    for i, xy_val in enumerate(xy):
                        self.results = self.results.append({'param_index': param_index, 'gen_model_index': model_index, 'KL_amp_index': fi,
                            'param': xy_names[i], 'val': xy_val, 'bias': bias[i], 'uncertainty': xy_unc[i]}, ignore_index=True)

                for mi, m in enumerate(self.generator_models):
                    y = self.yields
                    avg_yields = y[
                        (y.param_index == param_index) & (y.gen_model_index == mi) & (y.KL_amp_index == fi)].yield_set.values[0]
                    total_yields += avg_yields
                 
                    direct_result = fit_models[mi].fit(avg_yields, quiet=True)
                    add_result(mi, direct_result)
                
                direct_result = avg_fit_model.fit(total_yields, quiet=True)
                add_result(-1, direct_result)
        return

    def fit_physics_param(self, gamma_fitter, channel_num, name="default", quiet=False):
        n_phys = 5
        if channel_num == 1:
            n_phys=3

        def add_result(fit_row, phys_fit, input_param):
            direct_physics_param, direct_physics_param_err, direct_result = phys_fit
            self.fits = self.fits.append({'param_index': fit_row.param_index, 'gen_model_index': fit_row.gen_model_index, 'KL_amp_index': fit_row.KL_amp_index,
                'type': 'physics_param_fit', 'result': direct_result, 'success': direct_result.success}, ignore_index=True)
            
 
            bias =  np.array(direct_physics_param) - np.array(input_param)
            pp_names = ['g', 'r', 'd']
            for i, pp_val in enumerate(direct_physics_param):
                self.results = self.results.append({'param_index': fit_row.param_index, 'gen_model_index': fit_row.gen_model_index, 'KL_amp_index': fit_row.KL_amp_index,
                    'param': pp_names[i], 'val': pp_val, 'bias': bias[i], 'uncertainty': direct_physics_param_err[i]}, ignore_index=True)


        for fit_index, fit_row in self.fits[self.fits.type == 'observable_fit'].iterrows():
            physics_param = self.param.loc[fit_row.param_index]
            physics_param_val = np.array([physics_param.g_val, physics_param.r_val, physics_param.d_val])
            physics_param_deg = [physics_param.g_val_deg, physics_param.r_val, physics_param.d_val_deg]
            
            # the fit is a bit sensitive to initial parameters so try a default set, as well as something close to inputs
            phys_fit = gamma_fitter.fit_from_result(fit_row.result, start_guess = physics_param_val*0.95)
            phys_fit_default_start = gamma_fitter.fit_from_result(fit_row.result)
            if phys_fit_default_start[2].fun < phys_fit[2].fun and phys_fit_default_start[2].success:
                phys_fit = phys_fit_default_start
            add_result(fit_row, phys_fit, physics_param_deg)

                            # Add to the results 

    def get_results(self, var_name, result_name="default", result="bias", KL_results="all"):
        xy_vars = ["xm", "ym", "xp", "yp", "x_xi", "y_xi"]
        physics_param_vars = ["g", "r", "d", "r_xi", "d_xi"]
        if (not (var_name in xy_vars)) and (not (var_name in physics_param_vars)):
            raise ValueError("'{}' not an xy parameter (xm, ym, xp, yp, x_xi, y_xi) or physics parameter (g, r, d, r_xi, d_xi)".format(var_name))
        
        # Pick out the right parameter and return results for the full average fit
        res = self.results[(self.results.param == var_name) & (self.results.gen_model_index == -1)]
        obs_fits = self.fits[(self.fits.gen_model_index==-1) & (self.fits.type == 'observable_fit')]
        phys_fits = self.fits[(self.fits.gen_model_index==-1) & (self.fits.type == 'physics_param_fit')]
        res = res.merge(obs_fits, on=['gen_model_index', 'param_index', 'KL_amp_index'], how='outer')
        res = res.merge(phys_fits, on=['gen_model_index', 'param_index', 'KL_amp_index'], how='outer',
            suffixes=('_obs', '_phys'))
        res.success_obs = pd.to_numeric(res.success_obs)
        res.success_phys = pd.to_numeric(res.success_phys)
        
        # Print all failed fits
        mean_success = res.groupby(['gen_model_index', 'param_index']).mean()
        for p, row in mean_success[(mean_success.success_obs < 1) &( mean_success.success_obs < 1)].iterrows():
            print 'some fits did not converge for   (gen_model, param) =', p
        for p, row in mean_success[(mean_success.success_obs == 0) & (mean_success.success_obs == 0 )].iterrows():
            print 'NO fits converged for   (gen_model, param) =', p

        # only use fully converged fits
        res = res[(res.success_obs) & (res.success_phys)]
        
        res = res[['param', 'param_index', 'KL_amp_index', 'bias']]

        # calc mean bias and KL spread
        calc_group = res.groupby(['param', 'param_index'])
        res_mean = calc_group.mean()
        res_std = calc_group.std()
        res = res_mean.merge(res_std, on=['param', 'param_index'], suffixes=('_mean', '_std'))

        # finally add in the parameter values
        # also add in the xm, ym initial values
        res = res.join(self.param, on=['param_index']).reset_index()
        res['xm_val'] = res.r_val*np.cos((res.d_val_deg-res.g_val_deg)/180.*np.pi)
        res['ym_val'] = res.r_val*np.sin((res.d_val_deg-res.g_val_deg)/180.*np.pi)
        res['xp_val'] = res.r_val*np.cos((res.d_val_deg+res.g_val_deg)/180.*np.pi)
        res['yp_val'] = res.r_val*np.sin((res.d_val_deg+res.g_val_deg)/180.*np.pi)
        return res
        

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
            self.results, self.fits
            ], file, protocol=2)
        file.close()

    def load_results(self, file_name):
        file = open(file_name, "r")
        list_to_load = cPickle.load(file)
        file.close()
        self.results = list_to_load[0]
        self.fits = list_to_load[1]

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

    def save(self, file_name):
        file = open(file_name, "wb")
        cPickle.dump([
            self.yields, self.param, self.results, self.fits
            ], file, protocol=2)
        file.close()

    def load(self, file_name):
        file = open(file_name, "r")
        list_to_load = cPickle.load(file)
        file.close()
        self.yields = list_to_load[0]
        self.param = list_to_load[1]
        self.results = list_to_load[2]
        self.fits = list_to_load[3]

if __name__=="__main__":
    print ("Testing {}".format(os.path.basename(__file__)))
    from models import DefaultModel, DefaultCombinedDHModel
    from models import FullModel
    from PhysicalParameters import PhysicalParameters
    from GammaFitter import GammaFitter
    from GammaFitterXiSetup import GammaFitterXiSetup


    from datetime import datetime
    startTime = datetime.now()

    #do something


    # dcm = DefaultCombinedDHModel("../amplitude_calculations/output/KS_default.pickle")
    amp_dir = '/data/lhcb/users/bjoern/Analyses/B2DPi/data/various_studies/ks_cpv/amplitudes'
    
    dm = 1#DefaultModel("{}/Belle2018_KS_default.pickle".format(amp_dir))
    
    fm = 2#FullModel(
        # "{}/Belle2018_KS_default.pickle".format(amp_dir),
        # "{}/Belle2018_KL_default.pickle".format(amp_dir),
        # PhysicalParameters().default_eps,
        # PhysicalParameters(True))
    
    fm2 = 3#FullModel(
        # "{}/Belle2018_KS_default.pickle".format(amp_dir),
        # "{}/Belle2018_KL_default.pickle".format(amp_dir),
        # PhysicalParameters().default_eps*0,
        # PhysicalParameters(False))
    
    bs = MultimodelBiasStudy([fm, fm2], 
        param_sets = np.array([[uf.deg_to_rad(x), 0.1, uf.deg_to_rad(130)] for x in [75]]),
        KL_amplitude_files = ["{}/Belle2018_KL_seed_1.pickle".format(amp_dir)])
        # KL_amplitude_files = ["{}/Belle2018_KL_seed_1.pickle".format(amp_dir), "{}/Belle2018_KL_seed_2.pickle".format(amp_dir)])






    
    # bs.init_yields()
    # print bs.yields.head()
    # bs.save_yields('test.pickle')
    # bs.fit_observables(dm, 4)
    bs.load_results('test2.pickle')
    print bs.fits.head()
    # print bs.results.head()
    # print bs.results.tail()
    gf = GammaFitter()
    bs.fit_physics_param(gf, 1)
    print "THIS TOOK: ", datetime.now() - startTime 
    print bs.results.head()
    print bs.results.tail()
    print "THIS TOOK: ", datetime.now() - startTime