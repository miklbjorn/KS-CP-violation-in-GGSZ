# Submit complete study with either single or double channel fits

import os, sys, argparse, subprocess

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
        '-sp',
        '--scan_point',
        type=int,
        default=1,
)
parser.add_argument(
        '-pek',
        '--perfect_eff_Ki',
        action='store_true'
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
        '-qs',
        '--quantiles',
        type=int,
        default=1,
        choices=[1, 5, 9, 10, 20],
)
parser.add_argument(
        '-iv',
        '--iv_num',
        type=int,
        default=-1
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
if args.iv_num >= 0:
    points_name += '_ivnum{}'.format(args.iv_num)

input_dir = os.environ["COMMON_INPUTS_PATH_GGSZ"]
with open("{}/batch_templates/oxford_batch_template.sh".format(input_dir), "r") as file:
    template_text = file.read()

executable = 'run_momentum_averaged_bias_study.py'


subdir = os.getcwd()
job_content = "cd {}\n".format(subdir)
job_content += 'python {}'.format(executable)
for a in sys.argv[1:]:
    # pass all arguments on to the job
    job_content += ' {}'.format(a)

scan_point_name = '' if args.mode=='iv_scan' else '_sp{}'.format(args.scan_point)

name = 'p_avg_bias_{}_{}_{}_{}{}{}'.format(args.model, args.setup, args.channels, 
    args.n_kl, points_name, scan_point_name)


output_dir = '/data/lhcb/users/bjoern/Analyses/B2DPi/data/various_studies/ks_cpv/simple_full_studies'
template_text = template_text.format(
    name = name,
    log = "{}/batch_logs/{}.log".format(output_dir, name),
    queue = "normal",
    exec_job = job_content)

batch_script = '{}/batch_scripts/{}.sh'.format(subdir, name)
print batch_script
script_file = open(batch_script, 'wb')
script_file.write(template_text)
script_file.close()
subprocess.call(["qsub", batch_script])
