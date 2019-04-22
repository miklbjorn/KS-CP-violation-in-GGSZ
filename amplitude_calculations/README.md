# Calculating the KS and KL amplitudes
This folder holds the c++ program that calculates KS and KL amplitudes. They are subsequently loaded by the python program that runs the CPV studies. 

This version only includes the D decay amplitude model from EvtGen (see below).

## The calculated amplitudes and how to use them

## Recompiling the program
The program depends on [EvtGen](https://evtgen.hepforge.org/), which in turn depends on [HepMC](http://lcgapp.cern.ch/project/simu/HepMC/) (built on top of the old HepMC version 2.06.08, which works fine for the purposes here), as well as having a working [ROOT](https://root.cern.ch/) installation. The HepMC-source and EvtGen files are included in this git repository, so everything compiles out of the box.


ROOT is assumed to be installed on the system (else none of the python programs will work either).

Then you can build the main program (which automatically compiles the needed HepMC/EvtGenBase library). From the *amplitude_calculations* directory
```
mkdir build
cd build
cmake ..
make
make install
```

The c++ program uses the header-only [Templatized C++ Command Line Parser Library](http://tclap.sourceforge.net/) for argument passing. 

## Running the program
To make a histogram with a default KS amplitude call, with a different mode
```
cd .. # (only if you are still in build directory)

# default
./k_amplitudes


# generate KL amplitude with seed 33 for unknown, randomized parameters
./k_amplitude -k KL -s 33
```

The program has the following extra options:

```
USAGE: 

   ./k_amplitudes  [-r <double>] [-s <int>] [-n <int>] [-o <string>] [--]
                   [--version] [-h]


Where: 

   -k <string={KS/KL}>,  --kaon_type <string={KS/KL}>
     Specify wether you want a 'KS' (default) or 'KL' amplitude.

   -r <double>,  --r_value <double>
     r value used when transforming KS to KL amplitude. (default =
     0.236*0.236 = tan^2 \theta_cabibbo).

   -s <int>,  --seed <int>
     Seed used in KL amplitude variation. 0 (default) gives non-random r
     and delta's).

   -n <int>,  --num_points <int>
     Number of points along each of the s(K, h+) and s(K, pi-) axes for
     which amp should be calculated. (default: 500).

   -o <string>,  --output_dir <string>
     Directory in which to put output. (Default './output').

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.


   Enjoy your amplitudes!
```
