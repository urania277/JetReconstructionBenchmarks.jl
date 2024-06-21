.PHONY: all

# include ../Makefile.def

CXXFLAGS+=-MMD
CXXFLAGS+=-O2 -std=c++17
#CXXFLAGS+=-g -Og

PYTHIA_DIR=$(HOME)/.local/pythia
HEPMC3_DIR=/opt

all: genevts-pp genevts-ee genevts-AA

clean:
	$(RM) genevts *.o *.d

SRC = genevts.cc

DEPS:=$(addsuffix .d,$(basename $(SRC)))

genevts-pp: genevts-pp.cc
	$(LINK.cc) -o $@ -I $(PYTHIA_DIR)/include -I $(HEPMC3_DIR)/include -L $(PYTHIA_DIR)/lib -L $(HEPMC3_DIR)/lib $<  -lpythia8 -lHepMC3 

genevts-ee: genevts-ee.cc
	$(LINK.cc) -o $@ -I $(PYTHIA_DIR)/include -I $(HEPMC3_DIR)/include -L $(PYTHIA_DIR)/lib -L $(HEPMC3_DIR)/lib $<  -lpythia8 -lHepMC3 

genevts-AA: genevts-AA.cc
	$(LINK.cc) -o $@ -I $(PYTHIA_DIR)/include -I $(HEPMC3_DIR)/include -L $(PYTHIA_DIR)/lib -L $(HEPMC3_DIR)/lib $<  -lpythia8 -lHepMC3 


-include $(DEPS)
