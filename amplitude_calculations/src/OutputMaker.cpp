#include "OutputMaker.h"
#include <iostream>
#include <fstream>
#include "TH2F.h"
#include "TFile.h"

OutputMaker::OutputMaker(){};
OutputMaker::~OutputMaker(){};

void OutputMaker::write_output(
    const DtoKpipiAmplitude& amp,
    const std::string& output_dir,
    const int numpoints)
    const{

    std::string output_file_name = get_file_name(amp, numpoints);
    std::string output_destination = output_dir + "/" + output_file_name;

    // double min_s = amp.get_min_s();
    // double max_s = amp.get_max_s();
    double min_s = 0.3;
    double max_s = 3;
    double delta_s = (max_s - min_s)/(numpoints - 1);
    std::cout << "Calculating (" << numpoints << ")^2 points, " <<
    "in range [" << min_s << " - " << max_s << " (Gev/c^2)^2] for (Ks, pi+) and (Ks, pi-)" << std::endl;
    
    std::cout << "Writing amplitudes to " << output_destination << std::endl;

    TFile * file = new TFile(output_destination.c_str(), "recreate");
    TH2F * h_amp = new TH2F("h_amp", "", numpoints, min_s, max_s, numpoints, min_s, max_s);
    TH2F * h_ph  = new TH2F("h_ph" , "", numpoints, min_s, max_s, numpoints, min_s, max_s);

    for (int i = 0; i < numpoints; i++){
        for (int j = 0; j < numpoints; j++){
            double sKp = min_s + i*delta_s;
            double sKm = min_s + j*delta_s;
            dcomplex amp_val = amp.get_amplitude(sKp, sKm);
            h_amp->SetBinContent(j+1, i+1, abs(amp_val));
            h_ph ->SetBinContent(j+1, i+1, arg(amp_val));
        }
    }
    h_amp->Write("", 1);
    h_ph->Write("", 1);
    file->Close();

    // Finally, pickle the output file using python!
    std::cout << "Using python/ROOTFunctions.py to pickle ROOT output to python format:\n";
    std::string command = Form("python2.7 ../python/ROOTFunctions.py -amp %s", output_destination.c_str());
    std::cout << "command: " << command << std::endl;
    std::system(command.c_str());
};

std::string OutputMaker::get_file_name(
    const DtoKpipiAmplitude& amp,
    const int numpoints
    ) const {
    std::string num_point_name = "";
    if (numpoints != 500){
        num_point_name = "_n" + std::to_string(numpoints);
    }
    if (amp.get_kaon_type() == "KS" || amp.get_seed() == 0){
        return amp.get_amp_prefix() + amp.get_kaon_type() + "_default" + num_point_name  + ".root";
    } else {
        return amp.get_amp_prefix() + "KL_seed_" + std::to_string(amp.get_seed()) + num_point_name +".root";
    }
}