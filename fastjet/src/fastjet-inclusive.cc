// fastjet-finder.cc
// MIT Licenced, Copyright (c) 2023-2024 CERN
//
// Code to run and time the jet finding of against various
// HepMC3 input files

// Original version of this code Philippe Gras, IRFU
// Modified by Graeme A Stewart, CERN

#include "fastjet/ClusterSequence.hh"
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

#include "fastjet-utils.hh"

using namespace std;
using Time = std::chrono::high_resolution_clock;
using us = std::chrono::microseconds;

vector<fastjet::PseudoJet> run_fastjet_clustering(std::vector<fastjet::PseudoJet> input_particles,
  fastjet::Strategy strategy, fastjet::JetAlgorithm algorithm, double R, double ptmin) {
  // create a jet definition: a jet algorithm with a given radius parameter
  fastjet::RecombinationScheme recomb_scheme=fastjet::E_scheme;
  fastjet::JetDefinition jet_def(algorithm, R, recomb_scheme, strategy);


  // run the jet clustering with the above jet definition
  fastjet::ClusterSequence clust_seq(input_particles, jet_def);

  // // get the resulting jets ordered in pt
  vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));

  return inclusive_jets;
}

int main(int argc, char* argv[]) {
  // Default values
  int maxevents = -1;
  int trials = 8;
  string mystrategy = "Best";
  bool dump = false;
  int power = -1;
  double R = 0.4;
  double ptmin = 0.5;
  int njets = -1;
  double dmin = 0.0;
  string dump_file = "";

  string usage = " [-h] [-m max_events] [-n trials] [-s strategy] [-p power] [-d dump_file] <HepMC3_input_file>";

  int opt;
  while ((opt = getopt(argc, argv, "m:n:s:p:d:R:P:h")) != -1) {
    switch(opt) {
      case 'm':
        maxevents = stoi(optarg);
        break;
      case 'n':
        trials = stoi(optarg);
        break;
      case 's':
        mystrategy = string(optarg);
        break;
      case 'p':
        power = stoi(optarg);
        break;
      case 'R':
        R = stod(optarg);
        break;
      case 'P':
        ptmin = stod(optarg);
        break;
      case 'd':
        dump_file = string(optarg);
        dump = true;
        break;
      case 'h':
        std::cout << "Usage: " << argv[0] << usage << std::endl;
        std::cout << " HepMC3_input_file: File with input events in HepMC3 format"  << std::endl;
        std::cout << " -m max_events: default is " << maxevents << ", which is all the events in the file" << std::endl;
        std::cout << " -n trials: default is " << trials << ", which is the number of repeats to do" << std::endl;
        std::cout << " -s strategy: valid values are 'Best' (default), 'N2Plain', 'N2Tiled'" << std::endl;
        std::cout << " -p power: -1=antikt, 0=cambridge_achen, 1=inclusive kt" << std::endl;
        std::cout << " -R size: R parameter, cone size (default = 0.4)" << std::endl;
        std::cout << " -P pt_min: minimum pt for inclusive jet output (default = 5.0)" << std::endl;
        std::cout << " -d dump_file: output jets are printed to here (use '-' for stdout)" << std::endl;
        std::cout << " -h: print this message"  << std::endl;
        exit(-1);
        break;
      default:
        std::cout << "Usage: " << argv[0] << usage << std::endl;
        exit(EXIT_FAILURE);
    }
  }

  if (optind >= argc) {
    std::cerr << "No <HepMC3_input_file> argument after options" << std::endl;
    std::cout << "Usage: " << argv[0] << usage << std::endl;
    exit(EXIT_FAILURE);
  }
  if (optind < argc-1) {
    std::cerr << "Unexpected arguments after HepMC3 file (which must be the last argument):" << std::endl;
    for (auto arg = optind+1; arg != argc; arg++) cout << " " << argv[arg];
    cout << endl;
    std::cout << "Usage: " << argv[0] << usage << std::endl;
    exit(EXIT_FAILURE);
  }

  string input_file = argv[optind];

  // read in input events
  //----------------------------------------------------------
  auto events = read_input_events(input_file.c_str(), maxevents);
  
  // Set strategy
  fastjet::Strategy strategy = fastjet::Best;
  if (mystrategy == string("N2Plain")) {
    strategy = fastjet::N2Plain;
  } else if (mystrategy == string("N2Tiled")) {
    strategy = fastjet::N2Tiled;
  }

  auto algorithm = fastjet::antikt_algorithm;
  if (power == 0) {
    algorithm = fastjet::cambridge_aachen_algorithm;
  } else if (power == 1) {
    algorithm = fastjet::kt_algorithm;
  }

  std::cout << "Strategy: " << mystrategy << "; Alg: " << power << endl;

  auto dump_fh = stdout;
  if (dump) {
    if (dump_file != "-") {
      dump_fh = fopen(dump_file.c_str(), "w");
    }
  }

  double time_total = 0.0;
  double time_total2 = 0.0;
  double sigma = 0.0;
  double time_lowest = 1.0e20;
  for (long trial = 0; trial < trials; ++trial) {
    std::cout << "Trial " << trial << " ";
    auto start_t = std::chrono::steady_clock::now();
    for (size_t ievt = 0; ievt < events.size(); ++ievt) {
      auto inclusive_jets = run_fastjet_clustering(events[ievt], strategy, algorithm, R, ptmin);

      if (dump) {
         fprintf(dump_fh, "Jets in processed event %zu\n", ievt+1);
    
        // print out the details for each jet
        for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
          fprintf(dump_fh, "%5u %15.10f %15.10f %15.10f\n",
          i, inclusive_jets[i].rap(), inclusive_jets[i].phi(),
          inclusive_jets[i].perp());
        }
      }
    }
    auto stop_t = std::chrono::steady_clock::now();
    auto elapsed = stop_t - start_t;
    auto us_elapsed = double(chrono::duration_cast<chrono::microseconds>(elapsed).count());
    std::cout << us_elapsed << " us" << endl;
    time_total += us_elapsed;
    time_total2 += us_elapsed*us_elapsed;
    if (us_elapsed < time_lowest) time_lowest = us_elapsed;
  }
  time_total /= trials;
  time_total2 /= trials;
  if (trials > 1) {
    sigma = std::sqrt(double(trials)/(trials-1) * (time_total2 - time_total*time_total));
  } else {
    sigma = 0.0;
  }
  double mean_per_event = time_total / events.size();
  double sigma_per_event = sigma / events.size();
  time_lowest /= events.size();
  std::cout << "Processed " << events.size() << " events, " << trials << " times" << endl;
  std::cout << "Total time " << time_total << " us" << endl;
  std::cout << "Time per event " << mean_per_event << " +- " << sigma_per_event << " us" << endl;
  std::cout << "Lowest time per event " << time_lowest << " us" << endl;

  return 0;
}
