#ifndef __EvtGenAmplitude__
#define __EvtGenAmplitude__

#include <string>
#include "EvtGenBase/EvtComplex.hh"
#include "DtoKpipiAmplitude.h"

class EvtGenAmplitude : public DtoKpipiAmplitude {

    public:
        EvtGenAmplitude(std::string kaon_type, int seed=0, double r_value=0.236*0.236);
        ~EvtGenAmplitude(); // non-virtual: should not be a base class

        virtual dcomplex get_amplitude(double sKp, double sKm) const;
        virtual std::string get_amp_prefix() const {return "EvtGen_";};

    protected:
        void set_amplitude_factors();
        EvtComplex get_amplitude_factor(double r, double delta) const;

        // Factors between KS and KL amplitudes, A(KL)/A(KS) for the resonances
        EvtComplex _CF_factor; // Cabibbo favoured
        EvtComplex _DCS_factor; // Doube Cabibbo suppressed
        EvtComplex _nonRes_factor; // Non-resonant part
        std::vector<EvtComplex> _pipi_factors; // vector holding the individual pipi-resonance values
};


#endif
