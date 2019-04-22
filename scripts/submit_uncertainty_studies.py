import subprocess, os, time

experiments = ['LL_lhcb', 'DD_lhcb', 'VDX_belle']

n_KL = '5'
n_KL_momentum = '5'

iv = '1'
p = '500'

num_jobs = 0

for exp in experiments:

    # Default values for non-momentum averaged
    subprocess.call(['python', 'submit_simple_bias_study.py',
        '-s', exp, '-m', 'Belle2018', '-n', n_KL, '-p', p, '-iv', iv])
    num_jobs += 1

    # Model uncertainties
    subprocess.call(['python', 'submit_simple_bias_study.py',
        '-s', exp, '-m', 'Belle2010', '-n', n_KL, '-p', p, '-iv', iv])
    num_jobs += 1

    subprocess.call(['python', 'submit_simple_bias_study.py',
        '-s', exp, '-m', 'EvtGen', '-n', n_KL, '-p', p, '-iv', iv])
    num_jobs += 1

    # rho uncertainties
    subprocess.call(['python', 'submit_simple_bias_study.py',
        '-s', exp, '-rm', '0.9', '-n', n_KL, '-p', p, '-iv', iv])
    num_jobs += 1

    subprocess.call(['python', 'submit_simple_bias_study.py',
        '-s', exp, '-rm', '1.1', '-n', n_KL, '-p', p, '-iv', iv])
    num_jobs += 1

    # z uncertainties
    subprocess.call(['python', 'submit_simple_bias_study.py',
        '-s', exp, '-zm', '0.9', '-n', n_KL, '-p', p, '-iv', iv])
    num_jobs += 1

    subprocess.call(['python', 'submit_simple_bias_study.py',
        '-s', exp, '-zm', '1.1', '-n', n_KL, '-p', p, '-iv', iv])


    # momentum uncertainties


    for case in ['full', 'eps_only', 'mat_only']:
        # Default values for momentum averaged
        subprocess.call(['python', 'submit_momentum_averaged_bias_study.py',
            '-s', exp, '-c', case,'-m', 'Belle2018', '-n', n_KL_momentum, '-p', p, '-iv', iv])
        num_jobs += 1

        # Default values for momentum averaged
        subprocess.call(['python', 'submit_momentum_averaged_bias_study.py',
            '-s', exp, '-c', case,'-m', 'Belle2018', '-n', n_KL_momentum, '-p', p, '-iv', iv, '-qs', '10'])
        num_jobs += 1

        # Default values for momentum averaged
        subprocess.call(['python', 'submit_momentum_averaged_bias_study.py',
            '-s', exp, '-c', case,'-m', 'Belle2018', '-n', n_KL_momentum, '-p', p, '-iv', iv, '-qs', '20'])
        num_jobs += 1

print "submitted {} jobs".format(num_jobs)
