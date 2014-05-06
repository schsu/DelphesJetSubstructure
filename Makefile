PLATFORM_TYPE := $(shell uname -s)
FJFLAGS := $(shell fastjet-config --cxxflags)
FJLIBS 	:= $(shell fastjet-config --libs)
ROOTCFLAGS	:= $(shell root-config --cflags)
ROOTGLIBS	:= $(shell root-config --glibs)
#----------------------------------------------
DELPHESPATH := /MG5/Delphes/
#-----------------------------------------------

ifeq ($(PLATFORM_TYPE), Linux)
	SOFLAGS := -shared
else
	ifeq ($(PLATFORM_TYPE), Darwin)
		SOFLAGS := -dynamiclib -undefined dynamic_lookup
	endif
endif

FILES = JetSub.o TreeStruct.o
LINKS = JetSubDict.o
LIBS = Qjets/Qjets.o Qjets/QjetsPlugin.o $(DELPHESPATH)/libDelphes.so

all: JetSubDict.C JetSubDict.o Qjets/lib/libQjets.a JetSub.o TreeStruct.o JetSub.so
	
JetSub.a: $(LIBS) $(FILES) $(LINKS)
	ar cq JetSub.a $(FILES)
	
JetSub.so: $(LIBS) $(FILES) $(LINKS)
	gcc $(SOFLAGS) -o $@ $^
	
JetSubDict.C: $(FILES:.o=.h) JetSub_LinkDef.h
	rootcint -f $@ -c -p -I$(DELPHESPATH) $^
	
JetSubDict.o: JetSubDict.C JetSubDict.h
	gcc -fPIC -O3 -c -I$(DELPHESPATH) $(ROOTCFLAGS) $(FJFLAGS) $< -o $@
	
Qjets/lib/libQjets.a: Qjets/Qjets.C Qjets/QjetsPlugin.C
	cd Qjets && make
	
%.o: %.C %.h
	gcc -fPIC -O3 -c -I$(DELPHESPATH) $(ROOTCFLAGS) $(FJFLAGS) $< -o $@ 
	
clean:
	rm $(FILES) $(LINKS:.o=.*) JetSub.so