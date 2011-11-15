INSTALL_DIR?=/usr/local

all_static: ffindex_static
	cd src && make all_static

all: ffindex_static
	cd src && make all

hhblits_static: hhblits_static
	cd src && make hhblits_static

hhblits: cs ffindex
	cd src && make all

cs:
	cd lib/cs/src && make OPENMP=1 cssgd

ffindex:
	cd lib/ffindex && make

ffindex_static:
	cd lib/ffindex && make FFINDEX_STATIC=1
	
install:
	cd lib/ffindex && make install INSTALL_DIR=$(INSTALL_DIR)
	mkdir -p $(INSTALL_DIR)/bin
	install bin/hhblits $(INSTALL_DIR)/bin/hhblits
	install bin/hhblits $(INSTALL_DIR)/bin/cstranslate
	install bin/hhblits $(INSTALL_DIR)/bin/hhalign
	install bin/hhblits $(INSTALL_DIR)/bin/hhconsensus
	install bin/hhblits $(INSTALL_DIR)/bin/hhfilter
	install bin/hhblits $(INSTALL_DIR)/bin/hhmake
	install bin/hhblits $(INSTALL_DIR)/bin/hhsearch
	mkdir -p $(INSTALL_DIR)/lib/hh
	install data/context_data.lib $(INSTALL_DIR)/lib/hh/context_data.lib
	install data/cs219.lib $(INSTALL_DIR)/lib/hh/cs219.lib


clean:
	cd lib/cs/src && make clean
	cd lib/ffindex && make clean
	cd src && make clean
