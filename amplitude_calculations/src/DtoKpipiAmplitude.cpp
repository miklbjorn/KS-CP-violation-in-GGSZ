#include "DtoKpipiAmplitude.h"
#include "EvtGenBase/EvtDalitzPoint.hh"
#include <iostream>
#include <random>

DtoKpipiAmplitude::DtoKpipiAmplitude(std::string kaon_type, int seed, double r_value) :
_kaon_type(kaon_type),
_seed(seed),
_r_value(r_value),
_mD(1.86483),
_mKs(0.497611),
_mPi(0.13957061)
{
    if (_kaon_type == "KL" || _kaon_type == "KS" ) {
        std::cout << "Initialising amplitude for : " << _kaon_type << std::endl;
    } else {
        std::cout << "UNKNOWN KAON TYPE: " << _kaon_type << "! EXITING!\n";
        exit(-1);
    }


};

DtoKpipiAmplitude::~DtoKpipiAmplitude(){}

dcomplex DtoKpipiAmplitude::get_amplitude_factor(double r, double delta) const {
    if (r<0) r=0;
    dcomplex tmp = dcomplex(1);
    tmp -= dcomplex(2*r*cos(delta), 2*r*sin(delta));
    tmp *= dcomplex(-1.0, 0.);
    return tmp;
}

bool DtoKpipiAmplitude::is_valid_point(double sKp, double sKm) const {
    // Calculate the 3rd Dalitz
    double spm = _mD*_mD + _mKs*_mKs + 2*_mPi*_mPi- sKp - sKm;
    EvtDalitzPoint point( _mKs, _mPi, _mPi, sKp, spm, sKm );

    return point.isValid();

}



