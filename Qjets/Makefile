PLATFORM_TYPE := $(shell uname -s)
FJFLAGS := $(shell fastjet-config --cxxflags)
FJLIBS := $(shell fastjet-config --libs)

ifeq ($(PLATFORM_TYPE), Linux)
	SOFLAGS := -shared
else
	ifeq ($(PLATFORM_TYPE), Darwin)
		SOFLAGS := -dynamiclib -undefined dynamic_lookup
	endif
endif

FILES = Qjets.o QjetsPlugin.o

all: $(FILES) lib/libQjets.a 

lib/libQjets.a: $(FILES) $(FILES:.cc=.o)
	ar cq lib/libQjets.a $(FILES)

lib/libQjets.so: $(FILES) $(FILES:.cc=.o)
	gcc $(SOFLAGS) -o $@ $(FILES)

%.o: %.C %.h
	rm  lib/*; gcc -fPIC -O3 -c $(FJFLAGS) $< -o $@ 

clean:
	rm $(FILES) lib/*