# This script calculates and prints Fi, ci, si 
# for an input K amplitude and binning

amplitude_file = '../amplitude_calculations/output/EvtGen_KS_default.pickle'
binning_file = '../python/input/KsPiPi_optimal.pickle'

from models import DefaultModel

dm = DefaultModel(amplitude_file, binning_file=binning_file)
dm.print_Fi_ci_si()

# access them as arrays as: dm.Fi, dm.ci, dm.si
# indexing from -8 to 8 (as far as I remember)
