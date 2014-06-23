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

CC=gcc -fopenmp
INC=-I$(FFINDEX_DIR)
CFLAGS+=-std=c99 -O2 -g -fPIC -Wall -pedantic $(INC)

ifdef ICC
  CXX = env LANG=C icc
  CXXFLAGS += -fast -Wall -fno-strict-aliasing -finline-functions -funroll-loops -fp-model fast=2 -fast-transcendentals -wd981 -wd383 -wd2259 -wd1572 -I$(CS_DIR) -I$(FFINDEX_DIR)
  HHSEARCH_LIBS = -lpthread -lrt -openmp
  HHBLITS_LIBS = -lpthread -lrt -openmp
  CS_LIBS = -openmp
  OBJECT_LIBS = -openmp
else # ifndef ICC
  CXX = g++
  CXXFLAGS += -g -std=c++11 -ffast-math -O3 -Wall -Wno-vla -Wno-deprecated -Wno-char-subscripts -fno-strict-aliasing -I$(CS_DIR) $(BITS_FLAG) -I$(FFINDEX_DIR) -fopenmp
  HHSEARCH_LIBS = -lpthread
  HHBLITS_LIBS = -lpthread -fopenmp
  CS_LIBS = -fopenmp
  HHALIGN_LIBS = 
  OBJECT_LIBS = -fopenmp
	ifdef AVX2
	CXXFLAGS += -DAVX2_SUPPORT -mavx2
	else
	CXXFLAGS += -DSSE -msse3
	endif
  ifeq ($(OS), Darwin)
    CXXFLAGS+=-DHH_MAC -Wa,-q
  endif
  ifeq ($(OS), Linux)
    HHSEARCH_LIBS = -lrt
    HHBLITS_LIBS +=-lrt
    ifndef NO_PNG
      HHALIGN_LIBS = -DNO_FREETYPE -DHH_PNG -lpng -lz
    endif
  endif
endif # ifndef ICC

CPP = $(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) 
CPP_VALGRIND = $(CXX) -g $(CXXFLAGS)


SOURCES = hhsearch.cpp hhhitlist.cpp hhhit.cpp hhalignment.cpp hhhmm.cpp hhmatrices.cpp hhdecl.cpp hhfunc.cpp hhfullalignment.cpp hhhalfalignment.cpp hhutil.cpp util.cpp hash.cpp list.cpp hhviterbi.cpp hhviterbialgorithm.cpp hhviterbimatrix.cpp hhhmmsimd.cpp hhviterbirunner.cpp
HH_OBJECTS = hhalignment.o hhfullalignment.o hhhalfalignment.o hhhitlist.o hhhmm.o hhdecl.o hhhit.o hhutil.o util.o hhblits.o hhfunc.o hhmatrices.o hhdatabase.o hhprefilter.o hhviterbi.o hhviterbimatrix.o hhhmmsimd.o hhviterbirunner.o hhviterbialgorithm_with_celloff.o hhviterbialgorithm_without_celloff.o hhposteriormatrix.o hhposteriordecoderrunner.o hhposteriordecoder.o hhmacalgorithm.o hhforwardalgorithm.o hhbackwardalgorithm.o hhbacktracemac.o hhsearch.o hhalign_ori.o

CS_SOURCES_WITH_APPS = $(filter-out $(wildcard $(CS_DIR)/*test*.cc), $(wildcard $(CS_DIR)/*.cc))
CS_SOURCES = $(filter-out $(wildcard $(CS_DIR)/*app.cc), $(CS_SOURCES_WITH_APPS))
CS_OBJECTS = $(patsubst %.cc, %.o, $(CS_SOURCES))

HEADERS = hhhitlist.h hhhit.h hhalignment.h hhhmm.h hhfullalignment.h hhfullalignment.h hash.h list.h hhfunc.h hhmatrices.h

TARGETS = hhblits hhblits_omp hhblits_mpi hhsearch hhmake hhalign hhfilter hhconsensus cstranslate a3m_extractor a3m_reducer a3m_database_filter ffindex_from_fasta ffindex_from_fasta_with_split 
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

hhviterbialgorithm_with_celloff.o: hhviterbialgorithm.cpp
	$(CPP) -c -o hhviterbialgorithm_with_celloff.o hhviterbialgorithm.cpp -DVITERBI_CELLOFF=1

hhviterbialgorithm_without_celloff.o: hhviterbialgorithm.cpp 
	$(CPP) -c -o hhviterbialgorithm_without_celloff.o hhviterbialgorithm.cpp

hhalign_ori: hhalign_app.cpp $(HH_OBJECTS) $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	$(CPP) hhalign_app.cpp -o hhalign_ori $(HH_OBJECTS) $(CS_OBJECTS) $(HHSEARCH_LIBS) $(OBJECT_LIBS) $(FFINDEX_DIR)/libffindex.a

hhsearch: hhsearch_app.cpp $(HH_OBJECTS) $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	$(CPP) hhsearch_app.cpp -o hhsearch $(HH_OBJECTS) $(CS_OBJECTS) $(HHSEARCH_LIBS) $(OBJECT_LIBS) $(FFINDEX_DIR)/libffindex.a

hhsearch_static: hhsearch.cpp $(HH_OBJECTS) $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	$(CPP) -static hhsearch_app.cpp -o hhsearch $(HH_OBJECTS) $(CS_OBJECTS) $(HHSEARCH_LIBS) $(OBJECT_LIBS) $(FFINDEX_DIR)/libffindex.a

hhsearch_valgrind: hhsearch.cpp $(HH_OBJECTS) $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	$(CPP_VALGRIND) hhsearch.cpp -o hhsearch_valgrind $(HH_OBJECTS) $(CS_OBJECTS) $(HHSEARCH_LIBS) $(OBJECT_LIBS)

hhblits: hhblits_app.cpp $(HH_OBJECTS) $(CS_OBJECTS)
	$(CPP) -g hhblits_app.cpp -o hhblits $(HH_OBJECTS) $(CS_OBJECTS) $(HHBLITS_LIBS) $(OBJECT_LIBS) $(FFINDEX_DIR)/libffindex.a

hhblits_mpi: hhblits_mpi.cpp $(HH_OBJECTS) $(CS_OBJECTS)
	mpic++ $(CXXFLAGS) -g hhblits_mpi.cpp -o hhblits_mpi $(HH_OBJECTS) $(CS_OBJECTS) $(HHBLITS_LIBS) $(OBJECT_LIBS) $(FFINDEX_DIR)/libffindex.a
	
hhblits_omp: hhblits_omp.cpp $(HH_OBJECTS) $(CS_OBJECTS)
	mpic++ $(CXXFLAGS) -g hhblits_omp.cpp -o hhblits_omp $(HH_OBJECTS) $(CS_OBJECTS) $(HHBLITS_LIBS) $(OBJECT_LIBS) $(FFINDEX_DIR)/libffindex.a
	
	
hhblits_static: hhblits_app.cpp $(HH_OBJECTS) $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	$(CPP) -static hhblits_app.cpp -o hhblits $(HH_OBJECTS) $(CS_OBJECTS) $(HHBLITS_LIBS)  $(OBJECT_LIBS) $(FFINDEX_DIR)/libffindex.a

hhblits_valgrind: hhblits_app.cpp $(HH_OBJECTS) $(SOURCES) $(HEADERS)$(CS_OBJECTS)
	$(CPP_VALGRIND) hhblits_app.cpp -o hhblits_valgrind $(HH_OBJECTS) $(CS_OBJECTS) $(HHBLITS_LIBS)  $(OBJECT_LIBS) $(FFINDEX_DIR)/libffindex.a

hhmake: hhmake.cpp $(HH_OBJECTS) $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	$(CPP) hhmake.cpp -o hhmake $(HH_OBJECTS) $(CS_OBJECTS) $(OBJECT_LIBS) $(FFINDEX_DIR)/libffindex.a

hhmake_static: hhmake.cpp $(HH_OBJECTS) $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	$(CPP) -static hhmake.cpp -o hhmake $(HH_OBJECTS) $(CS_OBJECTS) $(OBJECT_LIBS) $(FFINDEX_DIR)/libffindex.a

hhalign: hhalign.cpp $(HH_OBJECTS) $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	$(CPP) hhalign.cpp -o hhalign $(HH_OBJECTS) $(CS_OBJECTS) $(HHALIGN_LIBS) $(OBJECT_LIBS) $(FFINDEX_DIR)/libffindex.a

hhalign_static: hhalign.cpp $(HH_OBJECTS) $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	$(CPP) -static hhalign.cpp -o hhalign $(HH_OBJECTS) $(CS_OBJECTS) $(HHALIGN_LIBS) $(OBJECT_LIBS) $(FFINDEX_DIR)/libffindex.a

hhalign_valgrind: hhalign.cpp $(HH_OBJECTS) $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	$(CPP_VALGRIND) hhalign.cpp -o hhalign_valgrind $(HH_OBJECTS) $(CS_OBJECTS) $(HHALIGN_LIBS) $(OBJECT_LIBS)

hhfilter: hhfilter.cpp $(HH_OBJECTS) $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	$(CPP) hhfilter.cpp -o hhfilter $(HH_OBJECTS) $(CS_OBJECTS) $(OBJECT_LIBS) $(FFINDEX_DIR)/libffindex.a

hhfilter_static: hhfilter.cpp $(HH_OBJECTS) $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	$(CPP) -static hhfilter.cpp -o hhfilter $(HH_OBJECTS) $(CS_OBJECTS) $(OBJECT_LIBS) $(FFINDEX_DIR)/libffindex.a

hhfilter_valgrind: hhfilter.cpp $(HH_OBJECTS) $(SOURCES) $(HEADERS)$(CS_OBJECTS)
	$(CPP_VALGRIND) hhfilter.cpp -o hhfilter_valgrind $(HH_OBJECTS) $(CS_OBJECTS) $(OBJECT_LIBS)

hhconsensus: hhconsensus.cpp $(HH_OBJECTS) $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	$(CPP) hhconsensus.cpp -o hhconsensus $(HH_OBJECTS) $(CS_OBJECTS) $(OBJECT_LIBS) $(FFINDEX_DIR)/libffindex.a

hhconsensus_static: hhconsensus.cpp $(HH_OBJECTS) $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	$(CPP) -static hhconsensus.cpp -o hhconsensus $(HH_OBJECTS) $(CS_OBJECTS) $(OBJECT_LIBS) $(FFINDEX_DIR)/libffindex.a

hhconsensus_valgrind: hhconsensus.cpp $(HH_OBJECTS) $(SOURCES) $(HEADERS) $(CS_OBJECTS)
	$(CPP_VALGRIND) hhconsensus.cpp -o hhconsensus_valgrind $(HH_OBJECTS) $(CS_OBJECTS) $(OBJECT_LIBS)

a3m_extractor: a3m_extractor.cpp a3m_compress.cpp a3m_compress.h	     
	$(CPP) a3m_extractor.cpp -o a3m_extractor a3m_compress.cpp a3m_compress.h $(FFINDEX_DIR)/libffindex.a
	
a3m_reducer: a3m_reducer.cpp a3m_compress.cpp a3m_compress.h	     
	$(CPP) a3m_reducer.cpp -o a3m_reducer a3m_compress.cpp a3m_compress.h $(FFINDEX_DIR)/libffindex.a
	
a3m_database_filter: a3m_database_filter.cpp	     
	$(CPP) a3m_database_filter.cpp -o a3m_database_filter a3m_compress.cpp a3m_compress.h $(FFINDEX_DIR)/libffindex.a

ffindex_from_fasta: ffindex_from_fasta.o
	$(CC) -o $@ ffindex_from_fasta.o $(FFINDEX_DIR)/libffindex.a

ffindex_from_fasta_with_split: ffindex_from_fasta_with_split.o
	$(CC) -o $@ ffindex_from_fasta_with_split.o $(FFINDEX_DIR)/libffindex.a

hhsuite-userguide.pdf: hhsuite-userguide.tex
	pdflatex hhsuite-userguide.tex
	bibtex hhsuite-userguide
	pdflatex hhsuite-userguide.tex
	pdflatex hhsuite-userguide.tex
	cp hhsuite-userguide.pdf ../

clean: 	
	-rm -f *.o cs/*.o
	-rm -f hhconsensus hhconsensus_valgrind hhfilter_valgrind hhfilter hhalign_valgrind hhalign hhmake hhblits_valgrind hhblits hhsearch_valgrind hhsearch cstranslate
