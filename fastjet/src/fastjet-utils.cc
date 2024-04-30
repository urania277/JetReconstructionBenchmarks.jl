// fastjet-lib.cc

// MIT Licenced, Copyright (c) 2023-2024 CERN
//
// Utilities common to all fastjet codes used for benchmarking and validation

// Original version of this code Philippe Gras, IRFU
// Modified by Graeme A Stewart, CERN

#include "fastjet/PseudoJet.hh"
#include <iostream> // needed for io
#include <cstdio>   // needed for io
#include <string>
#include <vector>
#include <chrono>

#include <unistd.h>
#include <stdlib.h>

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/ReaderAscii.h"

using namespace std;

vector<vector<fastjet::PseudoJet>> read_input_events(const char* fname, long maxevents = -1) {
  // Read input events from a HepMC3 file, return the events in a vector
  // Each event is a vector of initial particles
  
  HepMC3::ReaderAscii input_file (fname);

  int events_parsed = 0;

  std::vector<std::vector<fastjet::PseudoJet>> events;

  while(!input_file.failed()) {
    if (maxevents >= 0 && events_parsed >= maxevents) break;

    std::vector<fastjet::PseudoJet> input_particles;

    HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::MM);

    // Read event from input file
    input_file.read_event(evt);

    // If reading failed - exit loop
    if (input_file.failed()) break;

    ++events_parsed;
    input_particles.clear();
    input_particles.reserve(evt.particles().size());
    for(auto p: evt.particles()){
      if(p->status() == 1){
	      input_particles.emplace_back(p->momentum().px(),
				     p->momentum().py(),
				     p->momentum().pz(),
				     p->momentum().e());
      }
    }
    events.push_back(input_particles);
  }

  cout << "Read " << events_parsed << " events from " << fname << endl;
  return events;
}