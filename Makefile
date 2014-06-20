# C++ flags 
#-----------------------------------
CXX = g++

# set library flags for compiler
#-----------------------------------
FJFLAGS	:= $(shell fastjet-config --cxxflags)
FJLIBS	:= $(shell fastjet-config --libs )
ROOTCFLAGS	:= $(shell root-config --cflags)
ROOTGLIBS	:= $(shell root-config --glibs)
HepMClib        = -L$(HepMCSYS)/lib -lHepMC

INCLUDES = -I$(HepMCSYS)/include/ -I/usr/local/include/  -I$(ROOTSYS)/include -I$(DELPHESPATH)

# set platform specific flags
#-----------------------------------
PLATFORM_TYPE := $(shell uname -s)
ifeq ($(PLATFORM_TYPE), Linux)
	SOFLAGS := -shared
else
	ifeq ($(PLATFORM_TYPE), Darwin)
		SOFLAGS := -dynamiclib -undefined dynamic_lookup
	endif
endif

CXXFLAGS = -ansi -std=c++11 -g -pedantic -Wall -O3 -Wno-long-long $(FJFLAGS) $(ROOTCFLAGS) $(INCLUDES)
CXXFLAGSLINK = -Wall -g -O3 -Wno-long-long $(FJFLAGS) $(ROOTCFLAGS) $(INCLUDES)
# rules
# ----------------------------------------------------------------------------

.SUFFIXES:      .o .cxx .f .exe .C


# instructions for building a .o file from a .cxx file
# ----------------------------------------------------------------------------


FILES = Qjets.o QjetsPlugin.o 
FILEROOT  = HelperClasses.o LocalSettings.o HiggsAnalysis.o JetSubstructure.o Run.o 
MINIFILES = DelphesNTuple.o miniha.o HiggsHist.o HelperClasses.o LocalSettings.o AtlasStyle.o
EXTRAFILES = Extras.o LocalSettings.o HelperClasses.o AtlasStyle.o

all:  $(FILES) $(FILEROOT)  lib/libQjets.a  ha sample miniha extras

lib/libQjets.a: $(FILES) $(FILES:.cc=.o)
	ar cq lib/libQjets.a $(FILES)

lib/libQjets.so: $(FILES) $(FILES:.cc=.o)
	$(CXX) $(CXXFLAGS) $(SOFLAGS) -o $@ $(FILES)


ha: $(FILES) $(FILEROOT)
	$(CXX) $(CXXFLAGSLINK) $(FILEROOT) -lQjets  -lDelphes -L./lib -L$(DELPHESPATH) $(ROOTGLIBS) $(FJLIBS) -o $@

miniha: $(MINIFILES)
	$(CXX) $(CXXFLAGSLINK) $(MINIFILES) -lQjets  -lDelphes -L./lib -L$(DELPHESPATH) $(ROOTGLIBS) $(FJLIBS) -o $@

sample: sample.o
	$(CXX) $(CXXFFLAGSLINK) sample.o   -lDelphes -L./lib -L$(DELPHESPATH) $(ROOTGLIBS) $(FJLIBS) -o $@

extras: $(EXTRAFILES)
	$(CXX) $(CXXFFLAGSLINK) $(EXTRAFILES) -lDelphes -L$(DELPHESPATH) $(ROOTGLIBS) -o $@

.o: %.C %.h
	rm  lib/*; g++ -fPIC  -O3 -c $(ROOTCFLAGS)  $(FJFLAGS) $< -o $@ 

clean:
	rm $(FILES) $(FILEROOT)  lib/*
