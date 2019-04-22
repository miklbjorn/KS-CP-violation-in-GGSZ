# Submit complete study with either single or double channel fits

import os, sys, argparse, subprocess
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument(
        '-s',
        '--setup',
        type=str,
        choices = ['default', 'LL_lhcb', 'DD_lhcb', 'VDX_belle', 'VDX_belle_altSig', 'all_belle'],
        required=True)
parser.add_argument(
        '-c',
        '--channels',
        type=str,
        choices = ['1', '2'],
        default='1')
parser.add_argument(
        '-n',
        '--n_kl',
        type=int,
        default=50)
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
        '-pek',
        '--perfect_eff_Ki',
        action='store_true'
)
parser.add_argument(
        '-p',
        '--points',
        type=int,
        default=500
)
parser.add_argument(
        '-tp',
        '--time_range_points',
        type=int,
        default=100
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
parser.add_argument(
        '-iv',
        '--iv_num',
        type=int,
        default=-1
)

args = parser.parse_args()


points_name = ''
if not args.mode == 'iv_scan':
    points_name += '_'.format(args.mode)
if args.mode == 'r_scan':
    points_name += 'r_scan_dv{}'.format(args.delta_val)
if args.mode == 'delta_scan':
    points_name += 'delta_scan_rv{}'.format(args.r_val)
if not args.points==500:
    points_name += '_p{}'.format(args.points)
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
if not args.quantile == 0:
    points_name += '_q{:1.3f}'.format(args.quantile)
if args.model_Ki: 
    points_name += '_modelKi'
if args.perfect_eff_Ki: 
    points_name += '_perfecteffKi'
if args.iv_num >= 0:
    points_name += '_ivnum{}'.format(args.iv_num)

input_dir = os.environ["COMMON_INPUTS_PATH_GGSZ"]
with open("{}/batch_templates/oxford_batch_template.sh".format(input_dir), "r") as file:
    template_text = file.read()

executable = 'run_complete_simple_bias_study.py'
if args.channels == '2':
    print "NOT SET UP TO WORK WITH 2 CHANNELS YET!"
    sys.exit(-1)

subdir = os.getcwd()
job_content = "cd {}\n".format(subdir)
job_content += 'python run_complete_simple_bias_study.py'
for a in sys.argv[1:]:
    # pass all arguments on to the job
    job_content += ' {}'.format(a)

name = 'bias_{}_{}_{}_{}{}'.format(args.model, args.setup, args.channels, args.n_kl, points_name)


output_dir = '/data/lhcb/users/bjoern/Analyses/B2DPi/data/various_studies/ks_cpv/simple_full_studies'
template_text = template_text.format(
    name = name,
    log = "{}/batch_logs/{}.log".format(output_dir, name),
    queue = "short",
    exec_job = job_content)

batch_script = '{}/batch_scripts/{}.sh'.format(subdir, name)
print batch_script
script_file = open(batch_script, 'wb')
script_file.write(template_text)
script_file.close()
subprocess.call(["qsub", batch_script])
