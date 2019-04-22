#include "EvtGenAmplitude.h"
#include "EvtGenBase/EvtDalitzPoint.hh"
#include "EvtGenBase/EvtDalitzPlot.hh"
#include "EvtGenBase/EvtDalitzReso.hh"
#include "EvtGenBase/EvtCyclic3.hh"
#include <iostream>
#include <random>

EvtGenAmplitude::EvtGenAmplitude(std::string kaon_type, int seed, double r_value) :
DtoKpipiAmplitude(kaon_type, seed, r_value)
{
    // Calculate the amplitude factors relevant for _kaon_type and _seed.
    set_amplitude_factors();
};


EvtGenAmplitude::~EvtGenAmplitude(){}


void EvtGenAmplitude::set_amplitude_factors(){
    if (_kaon_type == "KS"){
        std::cout << "Using default KS amplitude factors\n";
        if (_seed) std::cout << "seed is non-zero (" << _seed << ") but IGNORED for KS\n";
        _CF_factor = EvtComplex(1.);
        _DCS_factor = EvtComplex(1.);
        _nonRes_factor = EvtComplex(1.);
        for (int i = 0; i < 7; i++){
            _pipi_factors.push_back(EvtComplex(1.));
        }

    } else if (_seed == 0){
        std::cout << "Using default KL amplitude factors\n";
        
        double r = _r_value;
        double delta = 0.;

        _CF_factor = EvtComplex(-1.);
        _DCS_factor = EvtComplex(1.);
        _nonRes_factor = EvtComplex(-1.);
        // _nonRes_factor = get_amplitude_factor(r, delta);
        for (int i = 0; i < 7; i++){
            _pipi_factors.push_back(get_amplitude_factor(r, delta));
        } 
    } else {
        std::cout << "Varying KL amplitude factors using seed: " << _seed << std::endl;
        double r_mean = _r_value;
        double r_std = 0.5*_r_value;
        double delta_min = 0.;
        double delta_max = 2*3.1415;

        std::default_random_engine generator(_seed);
        std::normal_distribution<double> r(r_mean, r_std);
        std::uniform_real_distribution<double> delta(delta_min, delta_max);

        _CF_factor = EvtComplex(-1.);
        _DCS_factor = EvtComplex(1.);
        _nonRes_factor = get_amplitude_factor(r(generator), delta(generator));
        for (int i = 0; i < 7; i++){
            _pipi_factors.push_back(get_amplitude_factor(r(generator), delta(generator)));
        }  
    }

    std::cout << "Factors to be multiplied onto D->KsPiPi model amplitudes:\n";
    std::cout << "  CF       factor: " << _CF_factor            << std::endl;
    std::cout << "  DSC      factor: " << _DCS_factor           << std::endl;
    std::cout << "  nonRes   factor: " << _nonRes_factor            << std::endl;
    std::cout << "  rho0     factor: " << _pipi_factors[0]            << std::endl;
    std::cout << "  omega    factor: " << _pipi_factors[1]            << std::endl;
    std::cout << "  f0_980   factor: " << _pipi_factors[2]            << std::endl;
    std::cout << "  f0_1370  factor: " << _pipi_factors[3]            << std::endl;
    std::cout << "  f2_1270  factor: " << _pipi_factors[4]            << std::endl;
    std::cout << "  sigma    factor: " << _pipi_factors[5]            << std::endl;
    std::cout << "  sigma2   factor: " << _pipi_factors[6]            << std::endl;

}

EvtComplex EvtGenAmplitude::get_amplitude_factor(double r, double delta) const{
    dcomplex factor = DtoKpipiAmplitude::get_amplitude_factor(r, delta);
    return EvtComplex(factor.real(), factor.imag());
}

dcomplex EvtGenAmplitude::get_amplitude(double sKp, double sKm) const {

    // If outside kinematic limits, return zero
    if (!is_valid_point(sKp, sKm)){
        return 0;
    }
    
    // Calculate the 3rd Dalitz
    double spm = _mD*_mD + _mKs*_mKs + 2*_mPi*_mPi- sKp - sKm;
    EvtDalitzPoint point( _mKs, _mPi, _mPi, sKp, spm, sKm );

    static const EvtDalitzPlot plot( _mKs, _mPi, _mPi, _mD );

    // This corresponds to relativistic Breit-Wigner distributions. Not K-matrix.
    // Defining resonances.
        // Useful constants.
    const EvtSpinType::spintype& SCALAR = EvtSpinType::SCALAR;
    const EvtSpinType::spintype& VECTOR = EvtSpinType::VECTOR;
    const EvtSpinType::spintype& TENSOR = EvtSpinType::TENSOR;

    const EvtDalitzReso::NumType& RBW   = EvtDalitzReso::RBW_CLEO_ZEMACH;
    const EvtDalitzReso::NumType& GS    = EvtDalitzReso::GS_CLEO_ZEMACH;

    const EvtCyclic3::Pair& AB = EvtCyclic3::AB;
    const EvtCyclic3::Pair& AC = EvtCyclic3::AC;
    const EvtCyclic3::Pair& BC = EvtCyclic3::BC;

    static EvtDalitzReso KStarm      ( plot, BC, AC, VECTOR, 0.893606, 0.0463407, RBW );
    static EvtDalitzReso KStarp      ( plot, BC, AB, VECTOR, 0.893606, 0.0463407, RBW );
    static EvtDalitzReso rho0        ( plot, AC, BC, VECTOR, 0.7758  , 0.1464   , GS  );
    static EvtDalitzReso omega       ( plot, AC, BC, VECTOR, 0.78259 , 0.00849  , RBW );
    static EvtDalitzReso f0_980      ( plot, AC, BC, SCALAR, 0.975   , 0.044    , RBW );
    static EvtDalitzReso f0_1370     ( plot, AC, BC, SCALAR, 1.434   , 0.173    , RBW );
    static EvtDalitzReso f2_1270     ( plot, AC, BC, TENSOR, 1.2754  , 0.1851   , RBW );
    static EvtDalitzReso K0Starm_1430( plot, BC, AC, SCALAR, 1.459   , 0.175    , RBW );
    static EvtDalitzReso K0Starp_1430( plot, BC, AB, SCALAR, 1.459   , 0.175    , RBW );
    static EvtDalitzReso K2Starm_1430( plot, BC, AC, TENSOR, 1.4256  , 0.0985   , RBW );
    static EvtDalitzReso K2Starp_1430( plot, BC, AB, TENSOR, 1.4256  , 0.0985   , RBW );
    static EvtDalitzReso sigma       ( plot, AC, BC, SCALAR, 0.527699, 0.511861 , RBW );
    static EvtDalitzReso sigma2      ( plot, AC, BC, SCALAR, 1.03327 , 0.0987890, RBW );
    static EvtDalitzReso KStarm_1680 ( plot, BC, AC, VECTOR, 1.677   , 0.205    , RBW );

    EvtComplex amp = 0.;



    // Adding terms to the amplitude with their corresponding amplitude and phase terms.
    // Multiply by the relevant correction factors (calculated during construction)
    amp += EvtComplex(   .848984 ,   .893618  )                                   * _nonRes_factor;
    amp += EvtComplex( -1.16356  ,  1.19933   ) * KStarm      .evaluate( point )  * _CF_factor;
    amp += EvtComplex(   .106051 , - .118513  ) * KStarp      .evaluate( point )  *_DCS_factor; // DCS
    amp += EvtComplex(  1.0      ,  0.0       ) * rho0        .evaluate( point )  *_pipi_factors[0]; // pi_res
    amp += EvtComplex( - .0249569,   .0388072 ) * omega       .evaluate( point )  *_pipi_factors[1]; // pi_res
    amp += EvtComplex( - .423586 , - .236099  ) * f0_980      .evaluate( point )  *_pipi_factors[2]; // pi_res
    amp += EvtComplex( -2.16486  ,  3.62385   ) * f0_1370     .evaluate( point )  *_pipi_factors[3]; // pi_res
    amp += EvtComplex(   .217748 , - .133327  ) * f2_1270     .evaluate( point )  *_pipi_factors[4]; // pi_res
    amp += EvtComplex(  1.62128  ,  1.06816   ) * K0Starm_1430.evaluate( point )  * _CF_factor;
    amp += EvtComplex(   .148802 ,   .0897144 ) * K0Starp_1430.evaluate( point )  *_DCS_factor; // DCS
    amp += EvtComplex(  1.15489  , - .773363  ) * K2Starm_1430.evaluate( point )  * _CF_factor;
    amp += EvtComplex(   .140865 , - .165378  ) * K2Starp_1430.evaluate( point )  *_DCS_factor; // DCS
    amp += EvtComplex( -1.55556  , - .931685  ) * sigma       .evaluate( point )  *_pipi_factors[5]; // pi_res
    amp += EvtComplex( - .273791 , - .0535596 ) * sigma2      .evaluate( point )  *_pipi_factors[6]; // pi_res
    amp += EvtComplex( -1.69720  ,   .128038  ) * KStarm_1680 .evaluate( point )  * _CF_factor;

    return dcomplex(real(amp), imag(amp));
}

