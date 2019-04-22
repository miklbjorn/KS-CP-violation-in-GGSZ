import os
import numpy as np
import UtilityFunctions as uf

class Binning:
    """Binning() holds the binning defintions

    Definitions are input as numpy arrays, retrieved from pickle files by UtilityFunctions.load_binning

    One must pass
    bin_defintions      : a 2D numpy array with bin definitions
    definition_s12/s13  : 1D arrays with s12/s13 values corresponding to the bins in bin_definitions
    application_s12/s13 : 1D arrays with s12/s13 values for which the bin_definitions are wanted

    if definition_{s12,s13} = application_{s12,s13} the Binning.get_bins() = bin_defintions
    But having separate definition/application arrays allow for more freedom

    """
    def __init__(self, bin_definitions, definition_s12, definition_s13, application_s12, application_s13):
        
        bin_num_x = len(application_s12)
        bin_num_y = len(application_s13)
        self.bins = np.empty((bin_num_x, bin_num_y))    # the array that will hold the bin definitions for the DP coords
                                                        # in application_s1k
        for x, sx in enumerate(application_s12):
            for y, sy in enumerate(application_s13):
                if sx > max(definition_s12):
                    bin_x = len(definition_s12) - 1
                else:
                    bin_x = np.argwhere(definition_s12 >= sx)[0][0] # Get the x-bin number in bin_definitions that correspond to sx in application_s12
                if sy > max(definition_s13):
                    bin_y = len(definition_s13) - 1
                else:
                    bin_y = np.argwhere(definition_s13 >= sy)[0][0] # Get the y-bin number in bin_definitions that correspond to sy in application_s13
            
                self.bins[x][y] = bin_definitions[bin_x][bin_y] # set correct bin definition

        self.s12 = application_s12
        self.s13 = application_s13
        self.bin_num = int(np.max(self.bins))

    def get_number_of_bins(self):
        return self.bin_num

    def get_bin_indices(self, i):
        return self.bins==i

    def get_bins(self):
        return self.bins

    def plot_bins(self):
        return "to be implemented"

if __name__ == "__main__":
    print ("Testing {}".format(os.path.basename(__file__)))

    bin_def, bin_def_s12, bin_def_s13 = uf.load_binning("input/KsPiPi_optimal.pickle")
    # b = Binning(bin_def, bin_def_s12, bin_def_s13, bin_def_s12, bin_def_s13)
    b = Binning(bin_def, bin_def_s12, bin_def_s13, [1.2, 2.3], [1.2, 2.3])
    print np.array_equal(bin_def, b.get_bins())
    print b.get_bins()