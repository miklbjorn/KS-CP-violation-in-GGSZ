# This module holds all interaction with ROOT
# It holds functions that load necessary input from ROOT
# files and saves them as python objects
# This allows for all other python code to not rely on the ROOT
# library, which makes much, much easier to use conda environemnts,
# jupyter notebooks, etc. (at least w/o somewhat cumbersome setup maneuvers)
#
# ROOTFunctions.py can also be run as a stand-alone script, doing the pickling
# call : python ROOTFunctions.py -amp/bin <root file>
# to pickle either amplitude or binning file

import ROOT, sys, cPickle, math, cmath
import numpy as np
from Amplitude import Amplitude
from Binning import Binning

def amplitude_from_root_to_python(amplitude_file):
    ''' Saves an amplitude definition from a ROOT file into a pickled
    numpy 2D array that can subsequently be loaded and passed to an Amplitude() instance.
    Also saves the s12 and s13 1D numpy arrays in the pickled file

    A function that loads the file correctly is defined in UtilityFunctions.py
    (in a separate file so that it does not need to import ROOT)
    '''
    print "Pickling amplitude from ROOT file:", amplitude_file

    if ".root" not in amplitude_file:
        print "Input file must be a ROOT file (no '.root' in '{}')".format(amplitude_file)
        sys.exit(-1)

    root_file = ROOT.TFile(amplitude_file, "read")
    amp_hist = root_file.Get("h_amp")
    ph_hist = root_file.Get("h_ph")
    
    nbin_x = amp_hist.GetNbinsX()
    nbin_y = amp_hist.GetNbinsY()

    A = np.empty((nbin_x, nbin_y), dtype=complex)
    s12 = np.empty((nbin_x), dtype=float)
    s13 = np.empty((nbin_y), dtype=float)

    for x in range(0, nbin_x):
        for y in range(0, nbin_y):
            s12ij = amp_hist.GetXaxis().GetBinCenter(x+1) # bin index starts at 1 (0=underflow)
            s13ij = amp_hist.GetYaxis().GetBinCenter(y+1)
            Aij   = amp_hist.GetBinContent(x+1, y+1)
            pij   = ph_hist.GetBinContent(x+1, y+1)
            A[x, y] = cmath.rect(Aij, pij)
            s12[x] = s12ij
            s13[y] = s13ij

    pickle_file = amplitude_file.replace(".root", ".pickle")
    print "Output file:", pickle_file

    file_object = open(pickle_file, "wb")
    cPickle.dump([A, s12, s13], file_object, protocol=2) 

    root_file.Close()
    file_object.close()

    return

def binning_from_root_to_python(binning_file):
    ''' Saves a binning definition from a ROOT file into a pickled
    numpy 2D array that can subsequently be loaded and passed to an Binning() instance.
    Also saves the s12 and s13 1D numpy arrays in the pickled file

    A function that loads the file correctly is defined in UtilityFunctions.py
    (in a separate file so that it does not need to import ROOT)
    '''

    print "Pickling binning from ROOT file:", binning_file

    if ".root" not in binning_file:
        print "Input file must be a ROOT file (no '.root' in '{}')".format(binning_file)
        sys.exit(-1)    
    bin_file = ROOT.TFile(binning_file, "r")
    bin_hist = bin_file.Get("dkpp_bin_h")
    
    bin_num_x = bin_hist.GetNbinsX()
    bin_num_y = bin_num_x
    s12 = np.array([])
    for x in range(1, bin_num_x + 1):
        s12 = np.append(s12, bin_hist.GetXaxis().GetBinCenter(x))
    s13 = s12

    bins = np.empty((bin_num_x, bin_num_y))
    for x, sx in enumerate(s12):
        for y, sy in enumerate(s13):
            bin_x = bin_hist.GetXaxis().FindBin(sx)
            bin_y = bin_hist.GetYaxis().FindBin(sy)
            raw_bin_no = bin_hist.GetBinContent(bin_x, bin_y)
            if raw_bin_no != 0: raw_bin_no = 9 - raw_bin_no # much switch ordering for KsPiPi bins (so far, all bins)
            bins[x, y] = raw_bin_no
            if x < y: bins[x, y] *= -1
            if x == y: bins[x, y] *= 0 # Do not assign the diagonal to any bins

    pickle_file = binning_file.replace(".root", ".pickle")
    print "Output file:", pickle_file

    file_object = open(pickle_file, "wb")
    cPickle.dump([bins, s12, s13], file_object, protocol=2)   

    file_object.close()  
    bin_file.Close()
    return

def babar_model_from_root_to_python(amplitude_file):
    ''' Saves an amplitude definition from a ROOT file into a pickled set of two
    numpy 2D array that can subsequently be loaded
    saves one for the abs value and one for the delta phi (no phi in the ROOT files)
    Also saves the s12 and s13 1D numpy arrays in the pickled file

    A function that loads the file correctly is defined in UtilityFunctions.py
    (in a separate file so that it does not need to import ROOT)
    '''
    print "Pickling amplitude from ROOT file:", amplitude_file

    if ".root" not in amplitude_file:
        print "Input file must be a ROOT file (no '.root' in '{}')".format(amplitude_file)
        sys.exit(-1)

    root_file = ROOT.TFile(amplitude_file, "read")
    amp_hist = root_file.Get("dkpp_2008_amp")
    ph_hist = root_file.Get("dkpp_2008_dph")
    
    nbin_x = amp_hist.GetNbinsX()
    nbin_y = amp_hist.GetNbinsY()

    A = np.empty((nbin_x, nbin_y), dtype=float)
    dph = np.empty((nbin_x, nbin_y), dtype=float)
    s12 = np.empty((nbin_x), dtype=float)
    s13 = np.empty((nbin_y), dtype=float)

    for x in range(0, nbin_x):
        for y in range(0, nbin_y):
            s12ij = amp_hist.GetXaxis().GetBinCenter(x+1) # bin index starts at 1 (0=underflow)
            s13ij = amp_hist.GetYaxis().GetBinCenter(y+1)
            Aij   = amp_hist.GetBinContent(x+1, y+1)
            pij   = ph_hist.GetBinContent(x+1, y+1)
            A[x, y] = Aij
            dph[x, y] = pij
            s12[x] = s12ij
            s13[y] = s13ij

    pickle_file = amplitude_file.replace(".root", ".pickle")
    print "Output file:", pickle_file

    file_object = open(pickle_file, "wb")
    cPickle.dump([A, dph, s12, s13], file_object, protocol=2) 

    root_file.Close()
    file_object.close()

    return
if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1]=="-amp":
        amplitude_from_root_to_python(sys.argv[2])
    elif len(sys.argv) > 2 and sys.argv[1]=="-bin":
        binning_from_root_to_python(sys.argv[2])
    elif len(sys.argv) > 2 and sys.argv[1]=="-babar":
        babar_model_from_root_to_python(sys.argv[2])
    else:
        print "Call using either option below:"
        print "  '-amp <root-file>' to pickle amplitude file"
        print "  '-bin <root-file>' to pickle binning file"




