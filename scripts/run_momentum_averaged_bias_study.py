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
from MultimodelBiasStudy import MultimodelBiasStudy
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
        '-ch',
        '--channels',
        type=str,
        choices = ['1', '2'],
        default='1')
parser.add_argument(
        '-c',
        '--case',
        type=str,
        default='full',
        choices=['full', 'eps_only', 'mat_only'])
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
        '-iv',
        '--iv_num',
        type=int,
        default=-1
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
        '-sp',
        '--scan_point',
        type=int,
        default=1,
)
parser.add_argument(
        '-tp',
        '--time_range_points',
        type=int,
        default=100,
)
parser.add_argument(
        '-p',
        '--points',
        type=int,
        default=500,
)
parser.add_argument(
        '-mk',
        '--model_Ki',
        action='store_true'
)
parser.add_argument(
        '-pek',
        '--perfect_eff_Ki',
        action='store_true'
)
parser.add_argument(
        '-qs',
        '--quantiles',
        type=int,
        default=5,
        choices=[1, 5, 9, 10, 20],
)

args = parser.parse_args()

points_name = ''
points_file_name = ''
include_absorbtion = None
quantile_suffix = ''
if not args.mode == 'iv_scan':
    points_name += '_{}'.format(args.mode)
if args.mode == 'r_scan':
    points_name += '_dv{}'.format(args.delta_val)
if args.mode == 'delta_scan':
    points_name += '_rv{}'.format(args.r_val)
if not args.quantiles == 1:
    points_name += '_q{}'.format(args.quantiles)
    quantile_suffix = '_q{}'.format(args.quantiles)
if args.model_Ki: 
    points_name += '_modelKi'
if args.perfect_eff_Ki: 
    points_name += '_perfecteffKi'
if not args.time_range_points == 100:
    points_name += '_tp{}'.format(args.time_range_points)
if not args.points==500:
    points_name += '_p{}'.format(args.points)
    points_file_name = '_n{}'.format(args.points)

amp_dir = '/data/lhcb/users/bjoern/Analyses/B2DPi/data/various_studies/ks_cpv/amplitudes'
output_dir = '/data/lhcb/users/bjoern/Analyses/B2DPi/data/various_studies/ks_cpv/simple_full_studies'
output_dir = '{}/setup_{}_{}{}_n_{}'.format(output_dir, args.model, args.setup, points_name, args.n_kl)
subprocess.call(['mkdir', '-p', output_dir])



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

quantile_suffixes = ['']
if args.quantiles == 5:
    quantile_suffixes = ['_q{:1.3f}'.format(q) for q in np.linspace(0.1, 0.9, 5)]
if args.quantiles == 9:
    quantile_suffixes = ['_q{:1.3f}'.format(q) for q in np.linspace(0.1, 0.9, 9)]
if args.quantiles == 10:
    quantile_suffixes = ['_q{:1.3f}'.format(q) for q in np.linspace(0.05, 0.95, 10)]
if args.quantiles == 20:
    quantile_suffixes = ['_q{:1.3f}'.format(q) for q in np.linspace(0.025, 0.975, 20)]

cases = [args.case]
full_models = []
eff = {}
Fi_hat = {}

pp = {}
pp_no_mat = {}
P = {}
PZ = {}
PT = {}

for q in quantile_suffixes:
    pp[q] = PhysicalParameters(include_mixing=True, sum_chi=include_absorbtion)
    pp_no_mat[q] = PhysicalParameters(include_mixing=False, sum_chi=include_absorbtion)

    if args.setup=='default':
        eff[q] = DefaultEfficiency(s13, s13)
    
    elif args.setup=='LL_lhcb':
        P[q] = PhaseSpaceDependentMomentum(s12, s13, setup='lhcb_LL_p'+q)
        PZ[q] = PhaseSpaceDependentMomentum(s12, s13, setup='lhcb_LL_pz'+q)
        pp[q].p = P[q].p
        pp[q].rho = 2.3/2.0*pp[q].rho
        pp_no_mat[q].p = P[q].p
        eff[q] = TimeDependentEfficiency(s12, s13, p= PZ[q].p, 
            z_min=0, z_max = 280)
    
    elif args.setup=='DD_lhcb':
        P[q] = PhaseSpaceDependentMomentum(s12, s13, setup='lhcb_DD_p'+q)
        PZ[q] = PhaseSpaceDependentMomentum(s12, s13, setup='lhcb_DD_pz'+q)
        pp[q].p = P[q].p
        pp[q].rho = 1.6/2.0*pp[q].rho
        pp_no_mat[q].p = P[q].p
        eff[q] = TimeDependentEfficiency(s12, s13, p= PZ[q].p, 
            z_min=280, z_max = 2350)
    
    elif args.setup=='VDX_belle':
        P[q] = PhaseSpaceDependentMomentum(s12, s13, setup='belle_VDX_p'+q)
        PT[q] = PhaseSpaceDependentMomentum(s12, s13, setup='belle_VDX_pt'+q)
        pp[q].p = P[q].p
        pp[q].rho = 4.0/2.0*pp[q].rho
        pp_no_mat[q].p = P[q].p
        eff[q] = TimeDependentEfficiency(s12, s13, p= PT[q].p, 
            z_min=0, z_max = 80)
    else:
        print "Unkown setup: {}. Exiting!".format(args.setup)
        sys.exit(-1)

    pp[q].update()
    pp_no_mat[q].update()


    eps = PhysicalParameters().default_eps
    param = pp[q]
    if args.case == 'mat_only':
        eps = 0
    if args.case == 'eps_only':
        param = pp_no_mat[q]
    full_models.append(FullModel(
        "{}/{}_KS_default{}.pickle".format(amp_dir, args.model, points_file_name), 
        "{}/{}_KL_default{}.pickle".format(amp_dir, args.model, points_file_name), 
        eps = eps, 
        parameters = param, 
        efficiency = eff[q],
        binning_file="../python/input/KsPiPi_optimal.pickle"))
    full_models[-1].set_time_range(eff[q].default_time_range())
    Fi_hat[q] = full_models[-1].get_Fi_hat()

print 'Calculating average Fi'
yields = np.zeros(32)
for fm in full_models:
    yields += fm.predict_yields([0., 0., 0.], 1000)
p_avg_Fi_hat = (yields[16:] + yields[15::-1])/sum(yields)

print 'Making fit models!'
fit_models = []
for q in quantile_suffixes:
    print "  ", q
    fit_eff = eff[q]
    if args.perfect_eff_Ki:
        fit_eff = DefaultEfficiency(s13, s13)
    fit_models.append(DefaultModel("{}/{}_KS_default{}.pickle".format(amp_dir, args.model, points_file_name),
                 binning_file="../python/input/KsPiPi_optimal.pickle",
                 efficiency = fit_eff))
    if (not args.model_Ki) and (not args.perfect_eff_Ki):
        # set data driven Fi UNLESS -mk flag is parsed
        fit_models[-1].set_Fi(Fi_hat[q])

print 'Making average fit model:'
p_avg_fit_model =  DefaultModel("{}/{}_KS_default{}.pickle".format(amp_dir, args.model, points_file_name),
                 binning_file="../python/input/KsPiPi_optimal.pickle",
                 efficiency = DefaultEfficiency(s13, s13)) # Does not matter

# avg the ci and si
print 'Setting average ci, si and Fi'
ci_sum = np.zeros(16)
si_sum = np.zeros(16)
for fm in fit_models:
    ci_sum += fm.ci
    si_sum += fm.si
ci_avg = ci_sum/len(fit_models)
si_avg = si_sum/len(fit_models)

p_avg_fit_model.set_ci_si(ci_avg, si_avg)
p_avg_fit_model.set_Fi(p_avg_Fi_hat)


KL_input_files = []
for i in range(0, args.n_kl):
    KL_input_files.append("{}/{}_KL_seed_{}{}.pickle".format(amp_dir, args.model, i+1, points_file_name))

# Input values for the single point studies
input_values = [
    [uf.deg_to_rad(75.), 0.1  , uf.deg_to_rad(130.)], # DK
    [uf.deg_to_rad(75.), 0.005, uf.deg_to_rad(300.)] # Dpi (1)
    # [uf.deg_to_rad(75.), 0.03 , uf.deg_to_rad(350.)], # Dpi (2)
    # [uf.deg_to_rad(75.), 0.1  , uf.deg_to_rad( 40.)], # DK*
    # [uf.deg_to_rad(75.), 0.2  , uf.deg_to_rad(330.)], # D*K
    # [uf.deg_to_rad(75.), 0.2  , uf.deg_to_rad(190.)] # DK0*
]

gf = GammaFitter()

if args.mode == 'iv_scan':
    print "making single studies for different input vals"
    for i, iv in enumerate(input_values):
        
        if args.iv_num >= 0 and not (i==args.iv_num):
            continue
        
        print "  study:", i
        bs_iv = MultimodelBiasStudy(full_models, # FullModel defined above
                          param_sets = [iv],
                          KL_amplitude_files = KL_input_files)
        bs_iv.init_yields()
        bs_iv.run(fit_models, p_avg_fit_model, gf)

        file_name = 'iv_{}_{}'.format(i, args.case)
        bs_iv.save('{}/momentum_avg_{}.pickle'.format(output_dir, file_name))


if args.mode == 'gamma_scan':
    input_gammas = np.linspace(45,90,10)
    params = np.array([[uf.deg_to_rad(input_gammas[args.scan_point]), 0.1, uf.deg_to_rad(130)]])
    print 'making gamma scans'
    bs_gamma = MultimodelBiasStudy(full_models, # FullModel defined above
                          param_sets=params,
                          KL_amplitude_files = KL_input_files)
    bs_gamma.init_yields()
    bs_gamma.run(fit_models, p_avg_fit_model, gf)

    file_name = 'g_scan_{}_sp{}'.format(args.case, args.scan_point)
    bs_gamma.save('{}/momentum_avg_{}.pickle'.format(output_dir, file_name))

if args.mode == 'r_scan':
    input_rs = np.linspace(0.005,0.05,8)
    input_rs = np.append(input_rs, [0.75, 0.1])
    params = np.array([[uf.deg_to_rad(75), input_rs[args.scan_point], uf.deg_to_rad(args.delta_val)]])
    print 'making r scans'
    bs_gamma = MultimodelBiasStudy(full_models, # FullModel defined above
                      param_sets=params,
                      KL_amplitude_files = KL_input_files)
    bs_gamma.init_yields()
    bs_gamma.run(fit_models, p_avg_fit_model, gf)

    file_name = 'r_scan_{}_sp{}'.format(args.case, args.scan_point)
    bs_gamma.save('{}/momentum_avg_{}.pickle'.format(output_dir, file_name))

if args.mode == 'delta_scan':
    input_ds = np.linspace(0,360,25)
    params = np.array([[uf.deg_to_rad(75), args.r_val, uf.deg_to_rad(input_ds[args.scan_point])]])
    print 'making delta scans'
    bs_gamma = MultimodelBiasStudy(full_models, # FullModel defined above
                      param_sets=params,
                      KL_amplitude_files = KL_input_files)
    bs_gamma.init_yields()
    bs_gamma.run(fit_models, p_avg_fit_model, gf)

    file_name = 'd_scan_{}_sp{}'.format(args.case, args.scan_point)
    bs_gamma.save('{}/momentum_avg_{}.pickle'.format(output_dir, file_name))

