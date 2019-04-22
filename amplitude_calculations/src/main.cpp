#include <iostream>
#include "DtoKpipiAmplitude.h"
#include "EvtGenAmplitude.h"
#include "Cleo2002Amplitude.h"
#include "BaBar2005Amplitude.h"
#include "Belle2006Amplitude.h"
#include "Belle2010Amplitude.h"
#include "Belle2018AmplitudeWithAlternativePhaseConvention.h"
#include "Belle2018Amplitude.h"
#include "OutputMaker.h"
#include "tclap/CmdLine.h"

int main(int argc , char* argv[] ){

    int seed, numpoints;
    double r_value;
    std::string kaon_type, output_dir, model;

    //////////////////////////////////////////////////////////////
    // ---- PROCESS ARGUMENTS ------------------------------------
    //////////////////////////////////////////////////////////////

    try { // The command line parsing
        TCLAP::CmdLine cmd("Enjoy your amplitudes!", ' ', "0.1");

        TCLAP::ValueArg<std::string> output_dir_arg(
            "o", // short name
            "output_dir", // long name
            "Directory in which to put output. (Default './output').", // description
            false, // is it required (no there is a default)
            "./output", // default
            "string");
        cmd.add(output_dir_arg);

        TCLAP::ValueArg<int> numpoints_arg(
            "n", // short name
            "num_points", // long name
            "Number of points along each of the s(K, h+) and s(K, pi-) axes for which amp should be calculated. (default: 500).", // description
            false, // is it required (no there is a default)
            500, // default
            "int");
        cmd.add(numpoints_arg);

        TCLAP::ValueArg<int> seed_arg(
            "s", // short name
            "seed", // long name
            "Seed used in KL amplitude variation. 0 (default) gives non-random r and delta's).", // description
            false, // is it required (no there is a default)
            0, // default
            "int");
        cmd.add(seed_arg);

        TCLAP::ValueArg<double> r_val_arg(
            "r", // short name
            "r_value", // long name
            "r value used when transforming KS to KL amplitude. (default = 0.236*0.236 = tan^2 \\theta_cabibbo).", // description
            false, // is it required (no there is a default)
            0.236*0.236, // default
            "double");
        cmd.add(r_val_arg);


        TCLAP::ValueArg<std::string> kaon_type_arg(
            "k", // short name
            "kaon_type", // long name
            "Specify wether you want a 'KS' (default) or 'KL' amplitude.", // description
            false, // is it required (no there is a default)
            "KS", // default
            "string={KS/KL}");
        cmd.add(kaon_type_arg);

        TCLAP::ValueArg<std::string> model_arg(
            "m",
            "model",
            "which model to use, from EvtGen (default), Belle2018, Belle2010, Belle2006, BaBar2005, Cleo2002",
            false,
            "EvtGen",
            "string={EvtGen (default), Belle2018, Belle2010, Belle2006, BaBar2005, Cleo2002}");
        cmd.add(model_arg);

        cmd.parse(argc, argv);

        output_dir = output_dir_arg.getValue();
        numpoints = numpoints_arg.getValue();
        seed = seed_arg.getValue();
        kaon_type = kaon_type_arg.getValue();
        model = model_arg.getValue();
        r_value = r_val_arg.getValue();

    } catch (TCLAP::ArgException &e){
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    }

    //////////////////////////////////////////////////////////////
    // ---- Calculate amplitude ----------------------------------
    //////////////////////////////////////////////////////////////

    std::cout << "Running KS-amplitude calculator\n";

    DtoKpipiAmplitude* amp;
    if (model=="EvtGen"){
        amp = new EvtGenAmplitude(kaon_type, seed, r_value);
    } else if (model=="Belle2018"){
        amp = new Belle2018Amplitude("./belle2018/belle2018_param.txt", kaon_type, seed, r_value);
    } else if (model=="Belle2018AltPhase"){
        amp = new Belle2018AmplitudeWithAlternativePhaseConvention("./belle2018/belle2018_param.txt", kaon_type, seed, r_value);
    } else if (model=="Belle2010"){
        amp = new Belle2010Amplitude(kaon_type, seed, r_value);
    } else if (model=="Belle2006"){
        amp = new Belle2006Amplitude(kaon_type, seed, r_value);
    } else if (model=="BaBar2005"){
        amp = new BaBar2005Amplitude(kaon_type, seed, r_value);
    } else if (model=="Cleo2002"){
        amp = new Cleo2002Amplitude(kaon_type, seed, r_value);
    } else {
        std::cout << "Unknown model: " << model << "! Exiting!\n";
        exit(-1);
    }


    OutputMaker out = OutputMaker();
    out.write_output(*amp, output_dir, numpoints);

    delete amp;

    return 0;
}