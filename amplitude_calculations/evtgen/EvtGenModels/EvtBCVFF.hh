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
// Module: EvtGen/EvtBCVFF.hh
//
// Description: Form factors for EvtBCVFF model
//
// Modification history:
//
//    AVL     Jul 6, 2012             Module created
//    AVL     Feb 5, 2018            D0, D*0 modes added
//
//------------------------------------------------------------------------

#ifndef EVTBCVFF_HH
#define EVTBCVFF_HH

#include "EvtGenBase/EvtSemiLeptonicFF.hh"

class EvtId;

class EvtBCVFF : public EvtSemiLeptonicFF {

public:
  EvtBCVFF(int idV, int fit);
  void getvectorff( EvtId parent, EvtId daught,
                       double t, double mass, double *a1f,
                       double *a2f, double *vf, double *a0f );

  void getscalarff(EvtId, EvtId, double, double, double*, 
		   double*);

  void gettensorff(EvtId, EvtId, double, double, double*, 
		   double*, double*, double*);

  void getbaryonff(EvtId, EvtId, double, double, double*, 
		   double*, double*, double*);

  void getdiracff(EvtId, EvtId, double, double, double*, double*,
                  double*, double*, double*, double*);

  void getraritaff(EvtId, EvtId, double, double, double*, double*, 
		   double*, double*, double*, double*, double*, double*);

private:
  int idVector, whichfit;
  double MBc, MD0, Mpsi, Mpsi2S, kappa, Mchi;
};

#endif

