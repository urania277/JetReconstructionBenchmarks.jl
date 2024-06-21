// Copyright (C) 2021 Torbjorn Sjostrand.
// Copyright (C) 2021 Ph. Gras.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Mikhail Kirsanov <Mikhail.Kirsanov@cern.ch>.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"

using namespace Pythia8;

int main() {

  // Interface for conversion from Pythia8::Event to HepMC
  // event. Specify file where HepMC events will be stored.
  Pythia8::Pythia8ToHepMC topHepMC("events-AuAu.hepmc3");

  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;

  // Setup the beams.
  pythia.readString("Beams:idA = 1000791970");
  pythia.readString("Beams:idB = 1000791970"); // The lead ions.
  pythia.readString("Beams:eCM = 2760.0");
  pythia.readString("Beams:frameType = 1");

  // Initialize the Angantyr model to fit the total and semi-inclusive
  // cross sections in Pythia within some tolerance.
  pythia.readString("HeavyIon:SigFitErr = "
                    "0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0");
  // These parameters are typicall suitable for sqrt(S_NN)=5TeV
  pythia.readString("HeavyIon:SigFitDefPar = "
                    "17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0");
  // A simple genetic algorithm is run for 20 generations to fit the
  // parameters.
  pythia.readString("HeavyIon:SigFitNGen = 20");

  pythia.init();
  Hist mult("charged multiplicity", 100, -0.5, 799.5);

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < 100; ++iEvent) {
    if (!pythia.next()) continue;

    // Find number of all final charged particles and fill histogram.
    int nCharged = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged())
        ++nCharged;
    mult.fill( nCharged );

    // Construct new empty HepMC event, fill it and write it out.
    topHepMC.writeNextEvent( pythia );


  // End of event loop. Statistics. Histogram.
  }
  pythia.stat();
  cout << mult;

  // Done.
  return 0;
}
