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
// Module: EvtBcSMuNu.cc
//
// Description: Routine to implement semileptonic Bc -> scalar l nu decays
//
// Modification history:
//
//    AVL     Jul 6, 2012             Module created
//    AVL     Feb 5, 2018            D0, D*0 modes added
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSemiLeptonicScalarAmp.hh"
#include <string>
#include <iostream>

#include "EvtGenModels/EvtBcSMuNu.hh"
#include "EvtGenModels/EvtBCSFF.hh"


using namespace std;



EvtBcSMuNu::~EvtBcSMuNu() {
}

std::string EvtBcSMuNu::getName(){
  return "BC_SMN";     
}


EvtDecayBase* EvtBcSMuNu::clone(){
  return new EvtBcSMuNu;

}

void EvtBcSMuNu::decay( EvtParticle *p ){
  p->initializePhaseSpace(getNDaug(),getDaugs());
  calcamp->CalcAmp(p,_amp2,ffmodel);
}


void EvtBcSMuNu::init(){
  checkNArg(1);
  checkNDaug(3);

  //We expect the parent to be a scalar 
  //and the daughters to be X lepton neutrino

  checkSpinParent(EvtSpinType::SCALAR);

  checkSpinDaughter(0,EvtSpinType::SCALAR);
  checkSpinDaughter(1,EvtSpinType::DIRAC);
  checkSpinDaughter(2,EvtSpinType::NEUTRINO);

  idScalar = getDaug(0).getId();
  whichfit = int(getArg(0)+0.1);

  ffmodel = new EvtBCSFF(idScalar,whichfit);
  
  calcamp = new EvtSemiLeptonicScalarAmp; 
 
}

void EvtBcSMuNu::initProbMax() {

  EvtId parId = getParentId();
  EvtId mesonId = getDaug(0);
  EvtId lepId = getDaug(1);
  EvtId nuId = getDaug(2);

  int nQ2Bins = 200;
  double maxProb = calcamp->CalcMaxProb(parId, mesonId, lepId, nuId, ffmodel, nQ2Bins);

  if (verbose()) {
    EvtGenReport(EVTGEN_INFO,"EvtBcSMuNu") << "Max prob = " << maxProb << endl;
  }

  setProbMax(maxProb);
  

}
