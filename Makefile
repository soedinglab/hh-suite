INSTALL_DIR?=/usr/local

hhblits: cs ffindex
	cd src && make hhblits

hhblits_static: ffindex_static
	cd src && make hhblits_static

cs:
	cd lib/cs/src && make OPENMP=1 cssgd

ffindex:
	cd lib/ffindex && make

ffindex_static:
	cd lib/ffindex && make FFINDEX_STATIC=1
	
install:
	cd lib/ffindex && make install INSTALL_DIR=$(INSTALL_DIR)
	install bin/hhblits $(INSTALL_DIR)/bin/hhblits

clean:
	cd lib/cs/src && make clean
	cd lib/ffindex && make clean
	cd src && make clean
