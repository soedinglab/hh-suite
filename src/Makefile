CXX = g++
CXXFLAGS = -O3 -Wall -Wno-char-subscripts -fno-strict-aliasing
CPP = $(CXX) $(CXXFLAGS)

SOURCES = hhsearch.C hhhitlist.C hhhit.C hhalignment.C hhhmm.C hhmatrices.C hhdecl.C hhfunc.C hhfullalignment.C hhfullalignment.C hash.C hhutil.C util.C list.C thread.C cs_counts.cpp amino_acid.cpp simple_cluster.cpp cluster.cpp training_profile.cpp profile.cpp sequence.cpp matrix.cpp

HEADERS = hhhitlist.h hhhit.h hhalignment.h hhhmm.h hhfullalignment.h hhfullalignment.h hash.h list.h cs_counts.h amino_acid.h simple_cluster.h cluster.h training_profile.h profile.h sequence.h matrix.h

TARGETS = hhsearch hhmake hhalign hhfilter hhconsensus

all: $(TARGETS)

hhsearch: hhsearch.C $(SOURCES) $(HEADERS)
	  $(CPP) -lpthread -lrt hhsearch.C -o ../bin/hhsearch

hhsearch_valgrind: hhsearch.C $(SOURCES) $(HEADERS)
	  $(CXX) -g -O -Wall -Wno-char-subscripts -fno-strict-aliasing -lpthread -lrt hhsearch.C -o ../bin/hhsearch_valgrind

hhmake: hhmake.C $(SOURCES) $(HEADERS)
	$(CPP) hhmake.C -o ../bin/hhmake

hhalign: hhalign.C $(SOURCES) $(HEADERS)
	 $(CPP) -I/usr/include/ -L/usr/lib -lpng -lz hhalign.C -o ../bin/hhalign

hhalign_valgrind: hhalign.C $(SOURCES) $(HEADERS)
	  $(CXX) -g -O -Wall -Wno-char-subscripts -fno-strict-aliasing -I/usr/include/ -L/usr/lib -lpng -lz hhalign.C -o ../bin/hhalign_valgrind

hhfilter: hhfilter.C $(SOURCES) $(HEADERS)
	  $(CPP) hhfilter.C -o ../bin/hhfilter

hhconsensus: hhconsensus.C $(SOURCES) $(HEADERS)
	     $(CPP) hhconsensus.C -o ../bin/hhconsensus

clean: 	
	-rm -f *.o 
	-rm -f ../bin/*
