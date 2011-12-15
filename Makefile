INSTALL_DIR?=$(PWD)
libdir=`([ -d /usr/lib64 ] && echo lib64) || echo lib`
# Overriding this is currently not fully supported as the code won't know
# to what this is set then.
INSTALL_LIB_DIR?=$(INSTALL_DIR)/$(libdir)/hh

dist_name=hh-suite-2.2.22

all_static: ffindex_static
	cd src && make all_static

all: ffindex_static
	cd src && make all

hhblits_static: hhblits_static
	cd src && make hhblits_static

hhblits: ffindex
	cd src && make all

#cs:
#	cd lib/cs/src && make OPENMP=1 cssgd

ffindex:
	cd lib/ffindex && make

ffindex_static:
	cd lib/ffindex && make FFINDEX_STATIC=1
	
install:
	cd lib/ffindex && make install INSTALL_DIR=$(INSTALL_DIR)
	mkdir -p $(INSTALL_DIR)/bin
	install src/hhblits     $(INSTALL_DIR)/bin/hhblits
	install src/cstranslate $(INSTALL_DIR)/bin/cstranslate
	install src/hhalign     $(INSTALL_DIR)/bin/hhalign
	install src/hhconsensus $(INSTALL_DIR)/bin/hhconsensus
	install src/hhfilter    $(INSTALL_DIR)/bin/hhfilter
	install src/hhmake      $(INSTALL_DIR)/bin/hhmake
	install src/hhsearch    $(INSTALL_DIR)/bin/hhsearch
	mkdir -p $(INSTALL_LIB_DIR)
	install data/context_data.lib $(INSTALL_LIB_DIR)/context_data.lib
	install data/cs219.lib        $(INSTALL_LIB_DIR)/cs219.lib
	install src/.hhdefaults       $(INSTALL_LIB_DIR)/hhdefaults

deinstall:
	rm -f $(INSTALL_DIR)/bin/hhblits $(INSTALL_DIR)/bin/cstranslate $(INSTALL_DIR)/bin/hhalign \
		$(INSTALL_DIR)/bin/hhconsensus $(INSTALL_DIR)/bin/hhfilter $(INSTALL_DIR)/bin/hhmake $(INSTALL_DIR)/bin/hhsearch
	rmdir $(INSTALL_DIR)/bin || true
	rm -f $(INSTALL_LIB_DIR)/context_data.lib $(INSTALL_LIB_DIR)/cs219.lib $(INSTALL_LIB_DIR)/hhdefaults
	rmdir $(INSTALL_LIB_DIR) || true
	cd lib/ffindex && make deinstall INSTALL_DIR=$(INSTALL_DIR)

clean:
	#cd lib/cs/src && make clean
	cd lib/ffindex && make clean
	cd src && make clean

dist/$(dist_name).tar.gz:
	mkdir -p dist
	git archive --prefix=$(dist_name)/ -o dist/$(dist_name).tar.gz HEAD
	cd dist && tar xf $(dist_name).tar.gz
	mkdir -p dist/$(dist_name)/bin
	cd dist/$(dist_name) && rsync --exclude .git --exclude .hg -av ../../lib .
	cd dist && tar czf $(dist_name).tar.gz $(dist_name)
