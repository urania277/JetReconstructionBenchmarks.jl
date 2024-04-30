#include "fastjet/ClusterSequence.hh"
#include <vector>

std::vector<std::vector<fastjet::PseudoJet>> read_input_events(const char* fname, long maxevents = -1);
