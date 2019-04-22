submit_gamma = False
submit_delta = False
submit_delta_dpi = True
submit_r = False

import subprocess, os
from itertools import product

qs = '5' # quantiles
n_KL = '5' # KL files
n_KL_dpi = '20' # fewer KL files in high resolution Dpi (20 was plenty anyway)
model = 'Belle2018'
experiments = ['LL_lhcb', 'DD_lhcb', 'VDX_belle']

if submit_gamma:
    cases = ['full', 'eps_only', 'mat_only']
    for i, exp, case in product(range(10), experiments, cases):
        subprocess.call(['python', 'submit_momentum_averaged_bias_study.py',
            '-s', exp,
            '-mode', 'gamma_scan',
            '-sp', str(i),
            '-m', model,
            '-n', n_KL,
            '-qs', qs,
            '-c', case])

if submit_delta:
    r_vals = [0.05, 0.1, 0.25]
    counter = 0
    for i, exp, r in product(range(25), experiments, r_vals):

        output_dir = '/data/lhcb/users/bjoern/Analyses/B2DPi/data/various_studies/ks_cpv/simple_full_studies'
        
        points_name = '_delta_scan_rv{}_q{}'.format(r, qs)

        output_dir = '{}/setup_{}_{}{}_n_{}'.format(output_dir, model, exp, points_name, n_KL)

        file_name = 'd_scan_{}_sp{}'.format('full', i)
        file_name = '{}/momentum_avg_{}.pickle'.format(output_dir, file_name)

        
        if not os.path.exists(file_name) or True:
            print ' ', file_name
            counter+=1
            subprocess.call(['python', 'submit_momentum_averaged_bias_study.py',
                '-s', exp,
                '-mode', 'delta_scan',
                '-sp', str(i),
                '-m', model,
                '-n', n_KL,
                '-qs', qs,
                '-p', '1000',
                '-rv', str(r)])
    print 'submitted {} jobs'.format(counter)

if submit_delta_dpi:
    r_vals = [0.005]
    counter = 0
    for i, exp, r in product(range(25), experiments, r_vals):

        output_dir = '/data/lhcb/users/bjoern/Analyses/B2DPi/data/various_studies/ks_cpv/simple_full_studies'
        
        points_name = '_delta_scan_rv{}_q{}_p2000'.format(r, qs)

        output_dir = '{}/setup_{}_{}{}_n_{}'.format(output_dir, model, exp, points_name, n_KL_dpi)

        file_name = 'd_scan_{}_sp{}'.format('full', i)
        file_name = '{}/momentum_avg_{}.pickle'.format(output_dir, file_name)

        
        if not os.path.exists(file_name) or True:
            print ' ', file_name
            counter+=1
            subprocess.call(['python', 'submit_momentum_averaged_bias_study.py',
                '-s', exp,
                '-mode', 'delta_scan',
                '-sp', str(i),
                '-m', model,
                '-n', n_KL_dpi,
                '-qs', qs,
                # '-p', '2000', # higher phase-space resolution needed for small r
                '-rv', str(r)])
    print 'submitted {} jobs'.format(counter)


if submit_r:
    d_vals = [130, 300]
    counter = 0
    for i, exp, d in product(range(10), experiments, d_vals):

        output_dir = '/data/lhcb/users/bjoern/Analyses/B2DPi/data/various_studies/ks_cpv/simple_full_studies'
        
        points_name = '_r_scan_rv{}_q5'.format(d)

        output_dir = '{}/setup_{}_{}{}_n_{}'.format(output_dir, model, exp, points_name, n_KL)

        file_name = 'r_scan_{}_sp{}'.format('full', i)
        file_name = '{}/momentum_avg_{}.pickle'.format(output_dir, file_name)

        if True:#not os.path.exists(file_name):
            counter+=1
            subprocess.call(['python', 'submit_momentum_averaged_bias_study.py',
                '-s', exp,
                '-mode', 'r_scan',
                '-sp', str(i),
                '-m', model,
                '-n', n_KL,
                '-qs', qs,
                '-dv', str(d)])
    print 'submitted {} jobs'.format(counter)
