from memory_profiler import profile
from models import DefaultModel, FullModel

# specific to Oxford cluster setup ...
amp_dir = '/data/lhcb/users/bjoern/Analyses/B2DPi/data/various_studies/ks_cpv/amplitudes'

@profile
def default_model():
    amplitude_file = '{}/EvtGen_KS_default.pickle'.format(amp_dir)
    binning_file = '../python/input/KsPiPi_optimal.pickle'
    dm = DefaultModel(amplitude_file, binning_file=binning_file)
    dm.print_Fi_ci_si()

@profile
def full_model():
    amplitude_file = '{}/EvtGen_KS_default.pickle'.format(amp_dir)
    amplitude_file2= '{}/EvtGen_KL_seed_3.pickle'.format(amp_dir)
    binning_file = '../python/input/KsPiPi_optimal.pickle'
    fm = FullModel(amplitude_file, amplitude_file2, binning_file=binning_file)
    n = fm.predict_yields([0.2, 0.1, 0.5], 2000)

if __name__ == '__main__':
    # default_model()
    full_model()
