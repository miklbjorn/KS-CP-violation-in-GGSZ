''' Utility functions that are needed in multiple places, but need not be in a specific class'''

import cmath, math
import autograd.numpy as np
import cPickle

cached_amplitudes = {}

def load_amplitude(file, get_s_array=True):
    '''Loads an amplitude from a .pickle file (created from ROOT file via ROOTFunctions.py)
    if get_s_array==true: also return 1D numpy arrays with the Dalitz plot coordinates
    corresponding to the amplitude values, for each axis in the Dalitz plot
    '''
    
    if file in cached_amplitudes:
        print 'cached', file
        packed_list = cached_amplitudes[file]
    else:
        print 'loading', file

        pickle_file = open(file, "r")
        packed_list = cPickle.load(pickle_file)
        cached_amplitudes[file] = packed_list
    
    A = packed_list[0]

    if get_s_array:
        s12 = packed_list[1]
        s13 = packed_list[2]
        return A, s12, s13
    return A

def load_binning(file, get_s_array=True):
    '''Loads an binning from a .pickle file (created from ROOT file via ROOTFunctions.py)
    if get_s_array==true: also return 1D numpy arrays with the Dalitz plot coordinates
    corresponding to the bin values, for each axis in the Dalitz plot
    '''
    pickle_file = open(file, "r")
    packed_list = cPickle.load(pickle_file)
    
    bins = packed_list[0]

    if get_s_array:
        s12 = packed_list[1]
        s13 = packed_list[2]
        return bins, s12, s13
    return bins

def load_babar_model(file, get_s_array=True):
    '''Loads an amplitude from a .pickle file (created from ROOT file via ROOTFunctions.py)
    if get_s_array==true: also return 1D numpy arrays with the Dalitz plot coordinates
    corresponding to the amplitude values, for each axis in the Dalitz plot
    '''
    pickle_file = open(file, "r")
    packed_list = cPickle.load(pickle_file)
    
    A = packed_list[0]
    dph = packed_list[1]

    if get_s_array:
        s12 = packed_list[2]
        s13 = packed_list[3]
        return A, dph, s12, s13
    return A

def get_spm(s12, s13):
    S12, S13 = np.meshgrid(s12, s13)
    mD = 1.865
    mK = 0.497
    mPi = 0.139
    minimum_possible_spm = (2*mPi)**2
    maximum_possible_spm = (mD-mK)**2
    spm = np.ones((len(s12),len(s13)))*mD**2 + mK**2 + 2*mPi**2 - S12 - S13
    spm[spm<minimum_possible_spm] = minimum_possible_spm
    spm[spm>maximum_possible_spm] = maximum_possible_spm
    return spm

    
def get_xy(physics_param):
    ''' takes an input vector with [gamma, rB, deltaB] and returns x, y
    angles in RADIANS
    '''
    gamma  = physics_param[0]
    rB     = physics_param[1]
    deltaB = physics_param[2]

    xm = rB * np.cos(deltaB - gamma)
    xp = rB * np.cos(deltaB + gamma)
    ym = rB * np.sin(deltaB - gamma)
    yp = rB * np.sin(deltaB + gamma)

    return [xm, ym, xp, yp]

def get_xy_xi(physics_param):
    ''' takes an input vector with [gamma, rB, deltaB] and returns x, y
    angles in RADIANS
    '''
    gamma  = physics_param[0]
    r_dk   = physics_param[1]
    d_dk   = physics_param[2]
    r_dpi  = physics_param[3]
    d_dpi  = physics_param[4]

    xm = r_dk * np.cos(d_dk - gamma)
    xp = r_dk * np.cos(d_dk + gamma)
    ym = r_dk * np.sin(d_dk - gamma)
    yp = r_dk * np.sin(d_dk + gamma)

    x_xi = (r_dpi/r_dk)*np.cos(d_dpi-d_dk)
    y_xi = (r_dpi/r_dk)*np.sin(d_dpi-d_dk)

    return [xm, ym, xp, yp, x_xi, y_xi]


def deg_to_rad(angle):
    return float(angle)/180.*math.pi

def rad_to_deg(angle):
    return float(angle)/math.pi*180.

def get_correlation_matrix(cov_matrix):
    corr_matrix = cov_matrix*0
    n = np.shape(cov_matrix)[0]
    for i in range(0, n):
        for j in range(0, n):
            corr_matrix[i, j] = cov_matrix[i, j]/math.sqrt(cov_matrix[i, i]*cov_matrix[j, j])
    return corr_matrix

def get_covariance_matrix(corr_matrix, uncertainty):
    cov_matrix = corr_matrix*0
    n = np.shape(cov_matrix)[0]
    for i in range(0, n):
        for j in range(0, n):
            cov_matrix[i, j] = corr_matrix[i, j]*uncertainty[i]*uncertainty[j]
    return cov_matrix

def print_xy_result(result):
    print "Fit results for x and y parameters:"
    print "xm = {: 1.4f} +/- {: 1.4f}".format(result.x[0], result.x_unc[0])
    print "ym = {: 1.4f} +/- {: 1.4f}".format(result.x[1], result.x_unc[1])
    print "xp = {: 1.4f} +/- {: 1.4f}".format(result.x[2], result.x_unc[2])
    print "yp = {: 1.4f} +/- {: 1.4f}".format(result.x[3], result.x_unc[3])
    print "Correlation matrix"
    print result.cor_mat[0:4, 0:4]

def print_xy_xi_result(result):
    print "Fit results for x and y parameters:"
    print "xm   = {: 1.4f} +/- {: 1.4f}".format(result.x[0], result.x_unc[0])
    print "ym   = {: 1.4f} +/- {: 1.4f}".format(result.x[1], result.x_unc[1])
    print "xp   = {: 1.4f} +/- {: 1.4f}".format(result.x[2], result.x_unc[2])
    print "yp   = {: 1.4f} +/- {: 1.4f}".format(result.x[3], result.x_unc[3])
    print "x_xi = {: 1.4f} +/- {: 1.4f}".format(result.x[4], result.x_unc[4])
    print "y_xi = {: 1.4f} +/- {: 1.4f}".format(result.x[5], result.x_unc[5])
    print "Correlation matrix"
    print result.cor_mat[0:6, 0:6]

def print_yields(N):
    bin_num = len(N)/2
    print " i : B+(i)  | B-(-i)"
    print "--------------------"
    for i in range(0, bin_num):
        this_bin = bin_num/2 - i
        if this_bin <= 0: this_bin-=1
        print "{:2d} : {:5.1f}  |  {:5.1f}".format(this_bin, N[bin_num-1-i], N[bin_num+i])
    print "--------------------"


def compare_yields(N1, N2):
    bin_num = min(len(N1), len(N2))/2

    if len(N1)!=len(N2):
        print "COMPARING YIELDS WItH DIFFERENT BIN NUMBERS, WILL NOT MAKE SENSE"

    print " i :  N1     |    N2    ||    N1    |    N2  "
    print " i : B+(i)   |  B+(i)   ||   B-(-i) |  B-(-i)"
    print "--------------------"
    for i in range(0, bin_num):
        this_bin = bin_num/2 - i
        if this_bin <= 0: this_bin-=1
        print "{:2d} : {:6.1f}  |  {:6.1f}  ||  {:6.1f}  |  {:6.1f}".format(
            this_bin, N1[bin_num-1-i], N2[bin_num-1-i], N1[bin_num+i], N2[bin_num+i])
    print "--------------------"



def fluctuate_yields(yields):
    fluctuated_yields = []
    for i, y in enumerate(yields):
        fluctuated_yields.append(np.random.poisson(y)) # numpy.poisson does work even if the input-average is a non-integer
    return fluctuated_yields





