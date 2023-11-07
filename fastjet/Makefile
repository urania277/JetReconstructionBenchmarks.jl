# Makefile for fastjet timing executable
# Adjust these options for your system, in particulat to point to
# the FastJet and HepMC3 install areas
HEPMC3_DIR=$(HOME)/.local/fastjet
FASTJET_DIR=$(HOME)/.local/fastjet

.PHONY: all clean

CXXFLAGS+=-MMD
CXXFLAGS+=-O2 -ggdb --std=c++17

all: fastjet-finder

clean:
	$(RM) -r fastjet-finder *.o *.d *.dSYM

SRC = fastjet-finder.cc

DEPS:=$(addsuffix .d,$(basename $(SRC)))

fastjet-finder: fastjet-finder.cc
	$(LINK.cc) -o $@ -Wno-deprecated-declarations  -I $(FASTJET_DIR)/include -L $(FASTJET_DIR)/lib \
-I $(HEPMC3_DIR)/include -L $(HEPMC3_DIR)/lib $< -l fastjet -lHepMC3 -lm

-include $(DEPS)
