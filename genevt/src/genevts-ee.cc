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
  Pythia8::Pythia8ToHepMC topHepMC("events-ee-120.hepmc3");

  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;

  // Allow no substructure in e+- beams: normal for corrected LEP data.
  pythia.readString("PDF:lepton = off");
  // Process selection.
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  // Switch off all Z0 decays and then switch back on those to quarks.
  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfAny = 1 2 3 4 5");

  // e+e- beams
  pythia.readString("Beams:idA =  11");
  pythia.readString("Beams:idB = -11");
  
  // Can set beam energy here, like this:
  // double mZ = pythia.particleData.m0(23);
  // pythia.settings.parm("Beams:eCM", mZ);
  // Or with a manual valus
  pythia.readString("Beams:eCM = 120.");

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
