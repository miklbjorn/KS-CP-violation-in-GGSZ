//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtBCSFF.cc
//
// Description: form factors for Bc->Slnu
//
// Modification history:
//
//    AVL     Jul 6, 2012             Module created
//    AVL     Feb 5, 2018            D0, D*0 modes added
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtBCSFF.hh"

#include "EvtGenBase/EvtId.hh"
#include <string>
#include <math.h>
#include "EvtGenBase/EvtPDL.hh"
#include <stdlib.h>
#include <iostream>

using namespace std;

EvtBCSFF::EvtBCSFF(int idS, int fit) {

  idScalar = idS;
  whichfit = fit;
  MBc = EvtPDL::getMeanMass(EvtPDL::getId("B_c+"));
  MD0 = EvtPDL::getMeanMass(EvtPDL::getId("D0"));
  return;
}

void EvtBCSFF::getscalarff(EvtId /*p*/, EvtId /*d*/, double t, double /*mass*/,
                           double *fpf, double *f0f) {
  
  double q2 = t;

  if (whichfit == 0) {
    *fpf = 1;
    *f0f = 0;
    return;
  }

  if (idScalar == EvtPDL::getId("chi_c0").getId() ) { // Bc -> chi_c0
    if (whichfit == 3) { // FF from Wang et al 10.1103/PhysRevD.79.114018
      double ratio = q2 / (MBc*MBc);
      
      double fpf_0 = 0.47;
      double fpf_c1 = 2.03;
      double fpf_c2 = 0.43;

      double f0f_0 = 0.47;
      double f0f_c1 = -0.45;
      double f0f_c2 = -1.31;
      
      *fpf = fpf_0*exp(fpf_c1*ratio+fpf_c2*ratio*ratio);
      *f0f = f0f_0*exp(f0f_c1*ratio+f0f_c2*ratio*ratio);
      return;

    } else {
      EvtGenReport(EVTGEN_ERROR,"EvtBCSFF") << "Must choose 0 (fpf = 1) or 3 (Wang).\n";
      ::abort();
    }
  } else if (idScalar == EvtPDL::getId("D0").getId() || idScalar == EvtPDL::getId("anti-D0").getId()) { // Bc -> D0
    if (whichfit == 1) {// FF from Kiselev:2002vz, tables III, IV
      double q2invmass = q2 / (MBc * MBc - MD0 * MD0);
      double den = 1 - q2 / (5.0 * 5.0);
      if (fabs(den) < 1e-10) {
	*fpf = 0;
	*f0f = 0;
      } else {
	double fPlus = 0.32 / den;
	double fMinus = -0.34 / den;
	*fpf = fPlus;
	*f0f = q2invmass * fMinus + fPlus;
      }
    } else if (whichfit == 2) {// FF from Ebert:2003cn, Fig 9
      double ratio=q2/MBc/MBc;
      double const fPlus_0 = 0.143, fPlus_a = 0.7, fPlus_b = 2.13;
      double const f0_0 = 0.136, f0_a = 1.63, f0_b = -0.139;
      *fpf = fPlus_0/ (1 - fPlus_a * ratio - fPlus_b * ratio * ratio);
      *f0f = f0_0/ (1 - f0_a * ratio - f0_b * ratio * ratio);
      return;
    } else {
      EvtGenReport(EVTGEN_ERROR,"EvtBCSFF") << "Should choose 1 (Kiselev:2002vz) or 2 (Ebert:2003cn).\n";
    }
  } else {
    EvtGenReport(EVTGEN_ERROR,"EvtBCSFF") << "Only chi_c0 and D0/anti-D0 implemented.\n";
    ::abort();
  }
}


void EvtBCSFF::getvectorff(EvtId, EvtId, double, double,
                           double*, double*, double*, double*) {

  EvtGenReport(EVTGEN_ERROR,"EvtGen") << "Not implemented :getvectorff in EvtBCSFF.\n";  
  ::abort();

}

void EvtBCSFF::gettensorff(EvtId /*p*/, EvtId /*d*/, double /*t*/, double /*mass*/,
                              double*, double*, double*, double*) {
  

  EvtGenReport(EVTGEN_ERROR,"EvtGen") << "Not implemented :gettensorff in EvtBCSFF.\n";  
  ::abort();

}

void EvtBCSFF::getbaryonff(EvtId, EvtId, double, double,
                           double*, double*, double*, double*) {
  
  EvtGenReport(EVTGEN_ERROR,"EvtGen") << "Not implemented :getbaryonff in EvtBCSFF.\n";  
  ::abort();

}

void EvtBCSFF::getdiracff(EvtId, EvtId, double, double,
                          double*, double*, double*, double*, double*, double*) {
  
  EvtGenReport(EVTGEN_ERROR,"EvtGen") << "Not implemented :getdiracff in EvtBCSFF.\n";
  ::abort();

}

void EvtBCSFF::getraritaff(EvtId, EvtId, double, double,
                           double*, double*, double*, double*, double*, double*, double*, double*) {
  
  EvtGenReport(EVTGEN_ERROR,"EvtGen") << "Not implemented :getraritaff in EvtBCSFF.\n";
  ::abort();

}

