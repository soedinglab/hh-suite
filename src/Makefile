CS_DIR = ../cs/src

# Convert BITS=?? to a -m flag
BITS_FLAG =
ifeq (32,$(BITS))
BITS_FLAG = -m32
endif
ifeq (64,$(BITS))
BITS_FLAG = -m64
endif

ifdef ICC
CXX = env LANG=C icc
CXXFLAGS = -fast -Wall -fno-strict-aliasing -finline-functions -funroll-loops -fp-model fast=2 -fast-transcendentals -wd981 -wd383 -wd2259 -wd1572 -D HH_SSE3 -I$(CS_DIR)
HHSEARCH_LIBS = -lpthread -lrt
HHBLITS_LIBS = -lpthread -lrt -openmp
OBJECT_LIBS = -openmp
else
CXX = g++
CXXFLAGS = -O3 -Wall -Wno-deprecated -Wno-char-subscripts -fno-strict-aliasing -DHH_SSE3 -march=core2 -I$(CS_DIR) $(BITS_FLAG)
HHSEARCH_LIBS = -lpthread -lrt
HHBLITS_LIBS = -lpthread -lrt -fopenmp
OBJECT_LIBS = -fopenmp
endif

CPP = $(CXX) $(CXXFLAGS)
CPP_VALGRIND = $(CXX) -g $(CXXFLAGS)

ifeq (32,$(BITS))
HHALIGN_LIBS = 
else
HHALIGN_LIBS = -DHH_PNG -lpng -lz -I/usr/include/ -L/usr/lib
endif

SOURCES = hhsearch.C hhhitlist.C hhhit.C hhalignment.C hhhmm.C hhmatrices.C hhdecl.C hhfunc.C hhfullalignment.C hhfullalignment.C hash.C hhutil.C util.C list.C thread.C 

CS_SOURCES_WITH_APPS = $(filter-out $(wildcard $(CS_DIR)/*test*.cc), $(wildcard $(CS_DIR)/*.cc))
CS_SOURCES = $(filter-out $(wildcard $(CS_DIR)/*app*.cc), $(CS_SOURCES_WITH_APPS))
CS_OBJECTS = $(patsubst %.cc, %.o, $(CS_SOURCES))

HEADERS = hhhitlist.h hhhit.h hhalignment.h hhhmm.h hhfullalignment.h hhfullalignment.h hash.h list.h

TARGETS = hhblits hhsearch hhmake hhalign hhfilter hhconsensus
TARGETS_STATIC = hhblits_static hhsearch_static hhmake_static hhalign_static hhfilter_static hhconsensus_static
TARGETS_VALGRIND = hhblits_valgrind hhsearch_valgrind hhalign_valgrind

all: $(TARGETS)

all_static: $(TARGETS_STATIC)

all_valgrind: $(TARGETS_VALGRIND)

$(CS_OBJECTS): $(CS_DIR)/%.o: $(CS_DIR)/%.cc
	$(CPP) $(OBJECT_LIBS) -c $< -o $@

hhsearch: hhsearch.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	 $(CPP) hhsearch.C -o ../bin/hhsearch $(CS_OBJECTS) $(HHSEARCH_LIBS)

hhsearch_static: hhsearch.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	 $(CPP) -static hhsearch.C -o ../bin/hhsearch $(CS_OBJECTS) $(HHSEARCH_LIBS)

hhsearch_valgrind: hhsearch.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	  $(CPP_VALGRIND) hhsearch.C -o ../bin/hhsearch_valgrind $(CS_OBJECTS) $(HHSEARCH_LIBS)

hhblits: hhblits.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	 $(CPP) hhblits.C -o ../bin/hhblits $(CS_OBJECTS) $(HHBLITS_LIBS)

hhblits_static: hhblits.C $(SOURCES) $(HEADERS)$(CS_OBJECTS) 
	 $(CPP) -static hhblits.C -o ../bin/hhblits $(CS_OBJECTS) $(HHBLITS_LIBS)

hhblits_valgrind: hhblits.C $(SOURCES) $(HEADERS)$(CS_OBJECTS) 
	  $(CPP_VALGRIND) hhblits.C -o ../bin/hhblits_valgrind $(CS_OBJECTS) $(HHBLITS_LIBS)

hhmake: hhmake.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	$(CPP) hhmake.C -o ../bin/hhmake $(CS_OBJECTS)

hhmake_static: hhmake.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	$(CPP) -static hhmake.C -o ../bin/hhmake $(CS_OBJECTS)

hhalign: hhalign.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	 $(CPP) hhalign.C -o ../bin/hhalign $(CS_OBJECTS) $(HHALIGN_LIBS)

hhalign_static: hhalign.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	 $(CPP) -static hhalign.C -o ../bin/hhalign $(CS_OBJECTS) $(HHALIGN_LIBS)

hhalign_valgrind: hhalign.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	  $(CPP_VALGRIND) hhalign.C -o ../bin/hhalign_valgrind $(CS_OBJECTS) $(HHALIGN_LIBS)

hhfilter: hhfilter.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	  $(CPP) hhfilter.C -o ../bin/hhfilter $(CS_OBJECTS)

hhfilter_static: hhfilter.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	  $(CPP) -static hhfilter.C -o ../bin/hhfilter $(CS_OBJECTS) 

hhconsensus: hhconsensus.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	     $(CPP) hhconsensus.C -o ../bin/hhconsensus $(CS_OBJECTS)

hhconsensus_static: hhconsensus.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	     $(CPP) -static hhconsensus.C -o ../bin/hhconsensus $(CS_OBJECTS)

hhconsensus_valgrind: hhconsensus.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	     $(CPP_VALGRIND) hhconsensus.C -o ../bin/hhconsensus_valgrind $(CS_OBJECTS)

clean: 	
	-rm -f *.o ../cs/src/*.o
	-rm -f ../bin/*
	-cp .hhdefaults ../bin/
