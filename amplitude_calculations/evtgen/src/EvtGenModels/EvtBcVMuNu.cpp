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
// Module: EvtBcVMuNu.cc
//
// Description: Routine to implement semileptonic Bc -> vector l nu decays
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
#include "EvtGenBase/EvtSemiLeptonicVectorAmp.hh"
#include <string>
#include <iostream>

#include "EvtGenModels/EvtBcVMuNu.hh"
#include "EvtGenModels/EvtBCVFF.hh"


using namespace std;



EvtBcVMuNu::~EvtBcVMuNu() {
}

std::string EvtBcVMuNu::getName(){
  return "BC_VMN";     
}


EvtDecayBase* EvtBcVMuNu::clone(){
  return new EvtBcVMuNu;

}

void EvtBcVMuNu::decay( EvtParticle *p ){

  p->initializePhaseSpace(getNDaug(),getDaugs());
  calcamp->CalcAmp(p,_amp2,ffmodel);
}


void EvtBcVMuNu::init(){
  checkNArg(1);
  checkNDaug(3);

  //We expect the parent to be a scalar 
  //and the daughters to be X lepton neutrino

  checkSpinParent(EvtSpinType::SCALAR);

  checkSpinDaughter(0,EvtSpinType::VECTOR);
  checkSpinDaughter(1,EvtSpinType::DIRAC);
  checkSpinDaughter(2,EvtSpinType::NEUTRINO);

  idVector = getDaug(0).getId();
  whichfit = int(getArg(0)+0.1);

  ffmodel = new EvtBCVFF(idVector,whichfit);

  calcamp = new EvtSemiLeptonicVectorAmp; 
 
}

void EvtBcVMuNu::initProbMax() {

  EvtId parId = getParentId();
  EvtId mesonId = getDaug(0);
  EvtId lepId = getDaug(1);
  EvtId nuId = getDaug(2);

  int nQ2Bins = 200;
  double maxProb = calcamp->CalcMaxProb(parId, mesonId, lepId, nuId, ffmodel, nQ2Bins);

  if (verbose()) {
    EvtGenReport(EVTGEN_INFO,"EvtBcVMuNu") << "Max prob = " << maxProb << endl;
  }

  setProbMax(maxProb);

}
