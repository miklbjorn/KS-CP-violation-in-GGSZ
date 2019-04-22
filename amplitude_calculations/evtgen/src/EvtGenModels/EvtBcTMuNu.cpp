
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
// Module: EvtBcTMuNu.cc
//
// Description: Routine to implement semileptonic B->T lnu decays
//
// Modification history:
//
//    AVL     July 6, 2012        Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSemiLeptonicTensorAmp.hh"
#include <string>
#include <iostream>

#include "EvtGenModels/EvtBcTMuNu.hh"
#include "EvtGenModels/EvtBCTFF.hh"


using namespace std;



EvtBcTMuNu::~EvtBcTMuNu() {
}

std::string EvtBcTMuNu::getName(){
  return "BC_TMN";     
}


EvtDecayBase* EvtBcTMuNu::clone(){
  return new EvtBcTMuNu;

}

void EvtBcTMuNu::decay( EvtParticle *p ){

  p->initializePhaseSpace(getNDaug(),getDaugs());
  calcamp->CalcAmp(p,_amp2,ffmodel);

}


void EvtBcTMuNu::init(){
  
  checkNArg(1);
  checkNDaug(3);

  //We expect the parent to be a scalar
  //and the daughters to be X lepton neutrino

  checkSpinParent(EvtSpinType::SCALAR);

  checkSpinDaughter(0,EvtSpinType::TENSOR);
  checkSpinDaughter(1,EvtSpinType::DIRAC);
  checkSpinDaughter(2,EvtSpinType::NEUTRINO);
  
  idTensor = getDaug(0).getId();
  whichfit = int(getArg(0)+0.1);

  ffmodel = new EvtBCTFF(idTensor,whichfit);
  
  calcamp = new EvtSemiLeptonicTensorAmp; 
 
}

void EvtBcTMuNu::initProbMax() {

  EvtId parId = getParentId();
  EvtId mesonId = getDaug(0);
  EvtId lepId = getDaug(1);
  EvtId nuId = getDaug(2);

  int nQ2Bins = 200;
  double maxProb = calcamp->CalcMaxProb(parId, mesonId, lepId, nuId, ffmodel, nQ2Bins);

  if (verbose()) {
    EvtGenReport(EVTGEN_INFO,"EvtBcTMuNu") << "Max prob = " << maxProb << endl;
  }

  setProbMax(maxProb);

}
