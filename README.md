# Neutral kaon CPV and material interactions in GGSZ measurements of the CKM angle gamma
This repository holds (part of) the code base used to make the calculations of the paper [CP violation and material interaction of neutral kaons in measurements of the CKM angle γ using B±→DK± decays where D→K0Sπ+π−](https://arxiv.org/abs/1904.01129). There is only part of the code because
1. some of the amplitude model implementations used in the paper are not made by me, and not public
2. a lot of random scripts have not been included.

The code has four parts:
1. a c++ program that calculates the D→K0Sπ+π− decay amplitude at a specified set of grid points in phase space
2. a set of python classes that can calculate experimental yields and fit them (among other things) given the output of step 1
3. a set of python scripts that use these classes to run all studies going into the paper
4. a set of Jupyter notebooks that make the plots of the paper, based on these studies

## 1. The amplitude calculations
This program is described in the [README in the amplitudes_calculations subdirectory](amplitude_calculations/README.md)

## 2. The python classes
Are in the [python subdir](python). They are reasonably self-explanatory - for their use, see the scripts described below.

## 3. The scripts
The scripts [scripts/run_complete_simple_bias_study.py](scripts/run_complete_simple_bias_study.py) and [scripts/run_momentum_averaged_bias_study.py](scripts/run_momentum_averaged_bias_study.py) run a well-defined bias study for a single momentum input, or averaging over a momentum distribution, respectively. There are also scripts that submit them to the Oxford batch system in the same folder (easy to modify for other batch systems, I imagine).

The scripts [scripts/submit_uncertainty_studies.py](scripts/submit_uncertainty_studies.py) and [scripts/submit_scans.py](scripts/submit_scans.py) runs the above scripts with a series of different inputs, in order to produce all of the calculations needed for the bias values and uncertainties in the paper.

## 4. The notebooks
There are various bits and pieces in the notebooks. The important one is [this one](notebooks/MakePaperBiasPlots.ipynb), because it makes the plots of the paper, including all uncertainty calculations etc.
