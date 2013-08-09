CS_DIR = cs
FFINDEX_DIR?=../lib/ffindex/src

OS:= $(shell uname)

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
  CXXFLAGS += -fast -Wall -fno-strict-aliasing -finline-functions -funroll-loops -fp-model fast=2 -fast-transcendentals -wd981 -wd383 -wd2259 -wd1572 -D HH_SSE3 -I$(CS_DIR) -I$(FFINDEX_DIR)
  HHSEARCH_LIBS = -lpthread -lrt
  HHBLITS_LIBS = -lpthread -lrt -openmp
  CS_LIBS = -openmp
  OBJECT_LIBS = -openmp
else # ifndef ICC
  CXX = g++
  CXXFLAGS += -g -O3 -Wall -Wno-vla -Wno-deprecated -Wno-char-subscripts -fno-strict-aliasing -I$(CS_DIR) $(BITS_FLAG) -I$(FFINDEX_DIR) -fopenmp
  HHSEARCH_LIBS = -lpthread
  HHBLITS_LIBS = -lpthread -fopenmp
  CS_LIBS = -fopenmp
  HHALIGN_LIBS = 
  OBJECT_LIBS = -fopenmp
  ifdef WITH_NATIVE
    CXXFLAGS += -march=native
  endif
  ifdef WITH_SSE41
    CXXFLAGS += -DHH_SSE41 -msse4.1 -DHH_SSSE3 -mssse3
  else
    ifndef NO_SSE2
      ifdef NO_SSE3
        CXXFLAGS += -DHH_SSE2 -msse2
      else
        CXXFLAGS += -DHH_SSE3 -msse3
      endif
    endif
  endif
  ifeq ($(OS), Darwin)
    CXXFLAGS+=-DHH_MAC
  endif
  ifeq ($(OS), Linux)
    HHSEARCH_LIBS = -lpthread -lrt
    HHBLITS_LIBS +=-lrt
    ifndef NO_PNG
      HHALIGN_LIBS = -DNO_FREETYPE -DHH_PNG -lpng -lz
    endif
  endif
endif # ifndef ICC

CPP = $(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS)
CPP_VALGRIND = $(CXX) -g $(CXXFLAGS)


SOURCES = hhsearch.C hhhitlist.C hhhit.C hhalignment.C hhhmm.C hhmatrices.C hhdecl.C hhfunc.C hhprefilter.C hhworker.C hhfullalignment.C hhhalfalignment.C hhutil.C util.C hash.C list.C

CS_SOURCES_WITH_APPS = $(filter-out $(wildcard $(CS_DIR)/*test*.cc), $(wildcard $(CS_DIR)/*.cc))
CS_SOURCES = $(filter-out $(wildcard $(CS_DIR)/*app.cc), $(CS_SOURCES_WITH_APPS))
CS_OBJECTS = $(patsubst %.cc, %.o, $(CS_SOURCES))

HEADERS = hhhitlist.h hhhit.h hhalignment.h hhhmm.h hhfullalignment.h hhfullalignment.h hash.h list.h

TARGETS = hhblits hhsearch hhmake hhalign hhfilter hhconsensus cstranslate
TARGETS_STATIC = hhblits_static hhsearch_static hhmake_static hhalign_static hhfilter_static hhconsensus_static cstranslate_static
TARGETS_VALGRIND = hhblits_valgrind hhsearch_valgrind hhalign_valgrind

all: $(TARGETS)

all_static: $(TARGETS_STATIC)

all_valgrind: $(TARGETS_VALGRIND)

$(CS_OBJECTS): $(CS_DIR)/%.o: $(CS_DIR)/%.cc
	$(CPP) $(OBJECT_LIBS) -c $< -o $@

cstranslate: $(CS_DIR)/cstranslate_app.cc $(CS_DIR)/*.cc $(CS_DIR)/*.h $(CS_OBJECTS)
	$(CPP) $(CS_LIBS) -lm -o cstranslate $(CS_DIR)/cstranslate_app.cc $(CS_OBJECTS)

cstranslate_static: $(CS_DIR)/cstranslate_app.cc $(CS_DIR)/*.cc $(CS_DIR)/*.h $(CS_OBJECTS)
	$(CPP) $(CS_LIBS) -static -lm -o cstranslate $(CS_DIR)/cstranslate_app.cc $(CS_OBJECTS)

hhsearch: hhsearch.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	 $(CPP) hhsearch.C -o hhsearch $(CS_OBJECTS) $(HHSEARCH_LIBS) $(OBJECT_LIBS)

hhsearch_static: hhsearch.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	 $(CPP) -static hhsearch.C -o hhsearch $(CS_OBJECTS) $(HHSEARCH_LIBS) $(OBJECT_LIBS)

hhsearch_valgrind: hhsearch.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	  $(CPP_VALGRIND) hhsearch.C -o hhsearch_valgrind $(CS_OBJECTS) $(HHSEARCH_LIBS) $(OBJECT_LIBS)

hhblits: hhblits.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	 $(CPP) hhblits.C -o hhblits $(CS_OBJECTS) $(HHBLITS_LIBS)  $(OBJECT_LIBS) $(FFINDEX_DIR)/libffindex.a

hhblits_static: hhblits.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	 $(CPP) -static hhblits.C -o hhblits $(CS_OBJECTS) $(HHBLITS_LIBS)  $(OBJECT_LIBS) $(FFINDEX_DIR)/libffindex.a

hhblits_valgrind: hhblits.C $(SOURCES) $(HEADERS)$(CS_OBJECTS)
	  $(CPP_VALGRIND) hhblits.C -o hhblits_valgrind $(CS_OBJECTS) $(HHBLITS_LIBS)  $(OBJECT_LIBS) $(FFINDEX_DIR)/libffindex.a

hhmake: hhmake.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	$(CPP) hhmake.C -o hhmake $(CS_OBJECTS) $(OBJECT_LIBS)

hhmake_static: hhmake.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	$(CPP) -static hhmake.C -o hhmake $(CS_OBJECTS) $(OBJECT_LIBS)

hhalign: hhalign.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	 $(CPP) hhalign.C -o hhalign $(CS_OBJECTS) $(HHALIGN_LIBS) $(OBJECT_LIBS)

hhalign_static: hhalign.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	 $(CPP) -static hhalign.C -o hhalign $(CS_OBJECTS) $(HHALIGN_LIBS) $(OBJECT_LIBS)

hhalign_valgrind: hhalign.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	  $(CPP_VALGRIND) hhalign.C -o hhalign_valgrind $(CS_OBJECTS) $(HHALIGN_LIBS) $(OBJECT_LIBS)

hhfilter: hhfilter.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	  $(CPP) hhfilter.C -o hhfilter $(CS_OBJECTS) $(OBJECT_LIBS)

hhfilter_static: hhfilter.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	  $(CPP) -static hhfilter.C -o hhfilter $(CS_OBJECTS) $(OBJECT_LIBS)

hhfilter_valgrind: hhfilter.C $(SOURCES) $(HEADERS)$(CS_OBJECTS)
	  $(CPP_VALGRIND) hhfilter.C -o hhfilter_valgrind $(CS_OBJECTS) $(OBJECT_LIBS)

hhconsensus: hhconsensus.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	     $(CPP) hhconsensus.C -o hhconsensus $(CS_OBJECTS) $(OBJECT_LIBS)

hhconsensus_static: hhconsensus.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	     $(CPP) -static hhconsensus.C -o hhconsensus $(CS_OBJECTS) $(OBJECT_LIBS)

hhconsensus_valgrind: hhconsensus.C $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	     $(CPP_VALGRIND) hhconsensus.C -o hhconsensus_valgrind $(CS_OBJECTS) $(OBJECT_LIBS)

hhsuite-userguide.pdf: hhsuite-userguide.tex
	pdflatex hhsuite-userguide.tex
	bibtex hhsuite-userguide
	pdflatex hhsuite-userguide.tex
	pdflatex hhsuite-userguide.tex
	cp hhsuite-userguide.pdf ../

clean: 	
	-rm -f *.o cs/*.o
	-rm -f hhconsensus hhconsensus_valgrind hhfilter_valgrind hhfilter hhalign_valgrind hhalign hhmake hhblits_valgrind hhblits hhsearch_valgrind hhsearch cstranslate
