# This script runs a complete bias study and
# saves all relevant results for further processing
# eg. in jupyter
# it can be run on the batch system

import sys, os, cPickle
import argparse
from DefaultEfficiency import DefaultEfficiency
from PhysicalParameters import PhysicalParameters
from GammaFitter import GammaFitter
from models import FullModel, DefaultModel
from PhaseSpaceDependentMomentum import PhaseSpaceDependentMomentum
from TimeDependentEfficiency import TimeDependentEfficiency
from BiasStudy import BiasStudy
import UtilityFunctions as uf
import numpy as np
import subprocess

# pass args

parser = argparse.ArgumentParser()
parser.add_argument(
        '-s',
        '--setup',
        type=str,
        required=True)
parser.add_argument(
        '-n',
        '--n_kl',
        type=int,
        default=100)
parser.add_argument(
        '-m',
        '--model',
        type=str,
        default='Belle2018',
        choices = ['EvtGen', 'Belle2018', 'Belle2018AltPhase', 'Belle2010', 'Belle2006', 'BaBar2005', 'Cleo2002'],
)
parser.add_argument(
        '-mode',
        '--mode',
        type=str,
        default='iv_scan',
        choices = ['iv_scan', 'gamma_scan', 'r_scan', 'delta_scan'],
)
parser.add_argument(
        '-rv',
        '--r_val',
        type=float,
        default=0.1,
)
parser.add_argument(
        '-dv',
        '--delta_val',
        type=float,
        default=130,
)
parser.add_argument(
        '-mk',
        '--model_Ki',
        action='store_true'
)
parser.add_argument(
        '-iv',
        '--iv_num',
        type=int,
        default=-1
)
parser.add_argument(
        '-pek',
        '--perfect_eff_Ki',
        action='store_true'
)
parser.add_argument(
        '-p',
        '--points',
        type=int,
        default=500,
)
parser.add_argument(
        '-tp',
        '--time_range_points',
        type=int,
        default=100,
)
# multipliers to vary certain parameters as systematics
parser.add_argument(
        '-pm',
        '--p_multiplier',
        type=float,
        default=1.0,
)
parser.add_argument(
        '-rm',
        '--rho_multiplier',
        type=float,
        default=1.0,
)
parser.add_argument(
        '-zm',
        '--zmax_multiplier',
        type=float,
        default=1.0,
)
parser.add_argument(
        '-a',
        '--absorbtion',
        type=int,
        default=0,
)
parser.add_argument(
        '-q',
        '--quantile',
        type=float,
        default=0
)
args = parser.parse_args()

points_name = ''
points_file_name = ''
include_absorbtion = None
quantile_suffix = ''
if not args.mode == 'iv_scan':
    points_name += '_'.format(args.mode)
if args.mode == 'r_scan':
    points_name += 'r_scan_dv{}'.format(args.delta_val)
if args.mode == 'delta_scan':
    points_name += 'delta_scan_rv{}'.format(args.r_val)
if not args.points==500:
    points_name += '_p{}'.format(args.points)
    points_file_name = '_n{}'.format(args.points)
if not args.time_range_points == 100:
    points_name += '_tp{}'.format(args.time_range_points)
if not args.p_multiplier == 1.0:
    points_name += '_pm{}'.format(args.p_multiplier)
if not args.rho_multiplier == 1.0:
    points_name += '_rm{}'.format(args.rho_multiplier)
if not args.zmax_multiplier == 1.0:
    points_name += '_zm{}'.format(args.zmax_multiplier)
if not args.absorbtion == 0:
    points_name += '_a{}'.format(args.absorbtion)
    if args.absorbtion == 1:
        include_absorbtion = 'high_p'
    if args.absorbtion == 2:
        include_absorbtion = 'low_p'
if not args.quantile == 0:
    points_name += '_q{:1.3f}'.format(args.quantile)
    quantile_suffix = '_q{:1.3f}'.format(args.quantile)
if args.model_Ki: 
    points_name += '_modelKi'
if args.perfect_eff_Ki: 
    points_name += '_perfecteffKi'

cases = ['full', 'eps_only', 'mat_only']

amp_dir = '/data/lhcb/users/bjoern/Analyses/B2DPi/data/various_studies/ks_cpv/amplitudes'
output_dir = '/data/lhcb/users/bjoern/Analyses/B2DPi/data/various_studies/ks_cpv/simple_full_studies'
output_dir = '{}/setup_{}_{}{}_n_{}'.format(output_dir, args.model, args.setup, points_name, args.n_kl)
subprocess.call(['mkdir', '-p', output_dir])
print 'output_dir:', output_dir

pp = PhysicalParameters(include_mixing=True, sum_chi=include_absorbtion)
pp_no_mat = PhysicalParameters(include_mixing=False, sum_chi=include_absorbtion)

# A default model
fm_default = FullModel(
        "{}/{}_KS_default{}.pickle".format(amp_dir, args.model, points_file_name), 
        "{}/{}_KL_default{}.pickle".format(amp_dir, args.model, points_file_name), 
        eps = PhysicalParameters().default_eps, 
        parameters = PhysicalParameters(), 
        binning_file="../python/input/KsPiPi_optimal.pickle") 

# extract the phase space coordinates from default model
s12 = fm_default.s12
s13 = fm_default.s13

if args.setup=='default':
    eff = DefaultEfficiency(s13, s13)
elif args.setup=='LL_lhcb':
    P = PhaseSpaceDependentMomentum(s12, s13, setup='lhcb_LL_p'+quantile_suffix)
    PZ = PhaseSpaceDependentMomentum(s12, s13, setup='lhcb_LL_pz'+quantile_suffix)
    pp.p = P.p*args.p_multiplier
    pp.rho = 2.3/2.0*pp.rho*args.rho_multiplier
    pp_no_mat.p = P.p*args.p_multiplier
    eff = TimeDependentEfficiency(s12, s13, p= PZ.p*args.p_multiplier, 
        z_min=0, z_max = 280*args.zmax_multiplier)
elif args.setup=='DD_lhcb':
    P = PhaseSpaceDependentMomentum(s12, s13, setup='lhcb_DD_p'+quantile_suffix)
    PZ = PhaseSpaceDependentMomentum(s12, s13, setup='lhcb_DD_pz'+quantile_suffix)
    pp.p = P.p*args.p_multiplier
    pp.rho = 1.6/2.0*pp.rho*args.rho_multiplier
    pp_no_mat.p = P.p*args.p_multiplier
    eff = TimeDependentEfficiency(s12, s13, p= PZ.p*args.p_multiplier, 
        z_min=280, z_max = 2350*args.zmax_multiplier)
elif args.setup=='VDX_belle':
    P = PhaseSpaceDependentMomentum(s12, s13, setup='belle_VDX_p'+quantile_suffix)
    PT = PhaseSpaceDependentMomentum(s12, s13, setup='belle_VDX_pt'+quantile_suffix)
    pp.p = P.p*args.p_multiplier
    pp.rho = 4.0/2.0*pp.rho*args.rho_multiplier
    pp_no_mat.p = P.p*args.p_multiplier
    eff = TimeDependentEfficiency(s12, s13, p= PT.p*args.p_multiplier, 
        z_min=0, z_max = 80*args.zmax_multiplier)
elif args.setup=='VDX_belle_altSig':
    P = PhaseSpaceDependentMomentum(s12, s13, setup='belle_VDX_p'+quantile_suffix)
    PT = PhaseSpaceDependentMomentum(s12, s13, setup='belle_VDX_pt'+quantile_suffix)
    pp = PhysicalParameters(use_low_p_dchi=True, sum_chi=include_absorbtion)
    pp_no_mat = PhysicalParameters(include_mixing=False, use_low_p_dchi=True, sum_chi=include_absorbtion)
    pp.p = P.p*args.p_multiplier
    pp.rho = 4.0/2.0*pp.rho*args.rho_multiplier
    pp_no_mat.p = P.p*args.p_multiplier
    eff = TimeDependentEfficiency(s12, s13, p= PT.p*args.p_multiplier, 
        z_min=0, z_max = 80*args.zmax_multiplier)
elif args.setup=='all_belle':
    P = PhaseSpaceDependentMomentum(s12, s13, setup='belle_all_p'+quantile_suffix)
    PT = PhaseSpaceDependentMomentum(s12, s13, setup='belle_all_pt'+quantile_suffix)
    pp.p = P.p*args.p_multiplier
    pp.rho = 2.5/2.0*pp.rho*args.rho_multiplier
    pp_no_mat.p = P.p*args.p_multiplier
    eff = TimeDependentEfficiency(s12, s13, p= PT.p*args.p_multiplier, 
        z_min=0, z_max = 1130*args.zmax_multiplier)
else:
    print "Unkown setup: {}. Exiting!".format(args.setup)
    sys.exit(-1)

pp.update()
pp_no_mat.update()

t_range = eff.default_time_range()
if not args.time_range_points == 100:
    t_min = np.min(t_range)
    t_max = np.max(t_range)
    t_range = np.linspace(t_min, t_max, args.time_range_points)

print 'Making FullModels and calculating Fis!'
full_models = {}
Fis = {}
for c in cases:
    print "  ",c
    eps = PhysicalParameters().default_eps
    param = pp
    if c=='mat_only':
        eps = 0
    if c=='eps_only':
        param = pp_no_mat
    full_models[c] = FullModel(
        "{}/{}_KS_default{}.pickle".format(amp_dir, args.model, points_file_name), 
        "{}/{}_KL_default{}.pickle".format(amp_dir, args.model, points_file_name), 
        eps = eps, 
        parameters = param, 
        efficiency = eff,
        binning_file="../python/input/KsPiPi_optimal.pickle")
    full_models[c].set_time_range(t_range)
    Fis[c] = full_models[c].get_Fi_hat()

print 'Making fit models!'
fit_models = {}
for c in cases:
    print "  ",c
    fit_eff = eff
    if args.perfect_eff_Ki:
        fit_eff = DefaultEfficiency(s13, s13)
    fit_models[c] = DefaultModel("{}/{}_KS_default{}.pickle".format(amp_dir, args.model, points_file_name),
                 binning_file="../python/input/KsPiPi_optimal.pickle",
                 efficiency = fit_eff)
    if not args.model_Ki and not args.perfect_eff_Ki:
        # set data driven Fi UNLESS -mk flag is parsed
        fit_models[c].set_Fi(Fis[c])


KL_input_files = []
for i in range(0, args.n_kl):
    KL_input_files.append("{}/{}_KL_seed_{}{}.pickle".format(amp_dir, args.model, i+1, points_file_name))

# Input values for the single point studies
input_values = [
    [uf.deg_to_rad(75.), 0.1  , uf.deg_to_rad(130.)],#, # DK
    [uf.deg_to_rad(75.), 0.005, uf.deg_to_rad(300.)] # Dpi (1)
    # [uf.deg_to_rad(75.), 0.03 , uf.deg_to_rad(350.)], # Dpi (2)
    # [uf.deg_to_rad(75.), 0.1  , uf.deg_to_rad( 40.)], # DK*
    # [uf.deg_to_rad(75.), 0.2  , uf.deg_to_rad(330.)], # D*K
    # [uf.deg_to_rad(75.), 0.2  , uf.deg_to_rad(190.)] # DK0*
]

gf = GammaFitter()

if args.mode == 'iv_scan':
    print "making single studies for different input vals"
    bs_iv = {}
    for i, iv in enumerate(input_values):
        if args.iv_num >= 0 and not (i==args.iv_num):
            continue
        for c in cases:
            print "  study:", i, c
            key = '{}_{}'.format(c, i)
            bs_iv[key] = BiasStudy(full_models[c], # FullModel defined above
                              param_sets = [iv],
                              KL_amplitude_files = KL_input_files)
            bs_iv[key].init_yields()
            bs_iv[key].run(fit_models[c], gf)

            file_name = 'iv_{}_{}'.format(i, c)
            bs_iv[key].save_yields('{}/iv_yields_{}.pickle'.format(output_dir, file_name))
            bs_iv[key].save_results('{}/iv_results_{}.pickle'.format(output_dir, file_name))
            bs_iv[key].save_params('{}/iv_params_{}.pickle'.format(output_dir, file_name))
            print 'Saving results to {}/iv_params_{}.pickle'.format(output_dir, file_name)

if args.mode == 'gamma_scan':
    input_gammas = np.linspace(45,90,10)
    params = np.array([[uf.deg_to_rad(x), 0.1, uf.deg_to_rad(130)] for x in input_gammas])
    print 'making gamma scans'
    bs_gamma = {}
    for c in cases:
        print '  study:', c
        bs_gamma[c] = BiasStudy(full_models[c], # FullModel defined above
                          param_sets=params,
                          KL_amplitude_files = KL_input_files)
        bs_gamma[c].init_yields()
        bs_gamma[c].run(fit_models[c], gf)

        file_name = 'g_scan_{}'.format( c)
        bs_gamma[c].save_yields('{}/g_scan_yields_{}.pickle'.format(output_dir, file_name))
        bs_gamma[c].save_results('{}/g_scan_results_{}.pickle'.format(output_dir, file_name))
        bs_gamma[c].save_params('{}/g_scan_params_{}.pickle'.format(output_dir, file_name))

if args.mode == 'r_scan':
    input_rs = np.linspace(0.005,0.3,10)
    params = np.array([[uf.deg_to_rad(75), x, uf.deg_to_rad(args.delta_val)] for x in input_rs])
    print 'making r scans'
    bs_gamma = {}
    for c in cases:
        print '  study:', c
        bs_gamma[c] = BiasStudy(full_models[c], # FullModel defined above
                          param_sets=params,
                          KL_amplitude_files = KL_input_files)
        bs_gamma[c].init_yields()
        bs_gamma[c].run(fit_models[c], gf)

        file_name = 'r_scan_{}'.format( c)
        bs_gamma[c].save_yields('{}/r_scan_yields_{}.pickle'.format(output_dir, file_name))
        bs_gamma[c].save_results('{}/r_scan_results_{}.pickle'.format(output_dir, file_name))
        bs_gamma[c].save_params('{}/r_scan_params_{}.pickle'.format(output_dir, file_name))

if args.mode == 'delta_scan':
    input_ds = np.linspace(0,360,10)
    params = np.array([[uf.deg_to_rad(75), args.r_val, uf.deg_to_rad(x)] for x in input_ds])
    print 'making delta scans'
    bs_gamma = {}
    for c in cases:
        print '  study:', c
        bs_gamma[c] = BiasStudy(full_models[c], # FullModel defined above
                          param_sets=params,
                          KL_amplitude_files = KL_input_files)
        bs_gamma[c].init_yields()
        bs_gamma[c].run(fit_models[c], gf)

        file_name = 'd_scan_{}'.format( c)
        bs_gamma[c].save_yields('{}/d_scan_yields_{}.pickle'.format(output_dir, file_name))
        bs_gamma[c].save_results('{}/d_scan_results_{}.pickle'.format(output_dir, file_name))
        bs_gamma[c].save_params('{}/d_scan_params_{}.pickle'.format(output_dir, file_name))

