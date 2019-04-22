//
//  Sample test program for running EvtGen
//  

#include "EvtGen/EvtGen.hh"

#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtHepMCEvent.hh"
#include "EvtGenBase/EvtSimpleRandomEngine.hh"
#include "EvtGenBase/EvtMTRandomEngine.hh"
#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtDecayBase.hh"

#ifdef EVTGEN_EXTERNAL
#include "EvtGenExternal/EvtExternalGenList.hh"
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <list>

bool filter(HepMC::GenEvent* event)
{
  bool hasLepton = false;
  bool hasCharm = false;

  for (HepMC::GenEvent::particle_iterator it = event->particles_begin(); 
        it != event->particles_end(); ++it)
  {
    if ( abs((*it)->pdg_id()) == 11 || abs((*it)->pdg_id()) == 13 )
    {
      hasLepton = true;
    }
    int id = abs((*it)->pdg_id());
    if ( id > 400 && id < 500 )
    {
      hasCharm = true;
    }
  }
  return (hasLepton && (!hasCharm)); 
}

int main(int argc, char** argv) {

  EvtParticle* parent(0);

  // Define the random number generator
  EvtRandomEngine* eng = 0;

#ifdef EVTGEN_CPP11
  // Use the Mersenne-Twister generator (C++11 only)
  eng = new EvtMTRandomEngine();
#else
  eng = new EvtSimpleRandomEngine();
#endif

  EvtRandom::setRandomEngine(eng);

  EvtAbsRadCorr* radCorrEngine = 0;
  std::list<EvtDecayBase*> extraModels;

#ifdef EVTGEN_EXTERNAL
  bool convertPythiaCodes(false);
  bool useEvtGenRandom(true);
  EvtExternalGenList genList(convertPythiaCodes, "", "gamma", useEvtGenRandom);
  radCorrEngine = genList.getPhotosModel();
  extraModels = genList.getListOfModels();
#endif

  //Initialize the generator - read in the decay table and particle properties
  EvtGen myGenerator("../DECAY.DEC","../evt.pdl", eng,
  		     radCorrEngine, &extraModels);

  //If I wanted a user decay file, I would read it in now.
  myGenerator.readUDecay("exampleFiles/Btousemileptonic.dec");

  static EvtId UPS4 = EvtPDL::getId(std::string("Upsilon(4S)"));

  int nEvents(40000);

  std::ofstream hepmcFile("hepMCtest");

  // Loop to create nEvents, starting from an Upsilon(4S)
  int i;
  for (i = 0; i < nEvents; i++) {

    std::cout<<"Event number "<<i<<std::endl;

    // Set up the parent particle
    EvtVector4R pInit(EvtPDL::getMass(UPS4), 0.0, 0.0, 0.0);
    parent = EvtParticleFactory::particleFactory(UPS4, pInit);
    parent->setVectorSpinDensity();      

    // Generate the event
    myGenerator.generateDecay(parent);    
    
    // Write out the results
    EvtHepMCEvent theEvent;
    theEvent.constructEvent(parent);
    HepMC::GenEvent* genEvent = theEvent.getEvent();

    if ( filter(genEvent) )
    {
//      genEvent->print(hepmcFile);
      hepmcFile<<(*genEvent);
    }

    parent->deleteTree();

  }

  hepmcFile.close();

  delete eng;
  return 0;

}