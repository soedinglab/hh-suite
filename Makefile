# This can be overridden e.g.: make install INSTALL_DIR=...
INSTALL_DIR?=$(PWD)

# Guess wether to use lib or lib64
libdir=`([ -d /usr/lib64 ] && echo lib64) || echo lib`

# Overriding this is currently not fully supported as the code won't know
# to what this is set then.
INSTALL_LIB_DIR?=$(INSTALL_DIR)/$(libdir)/hh
INSTALL_SCRIPT_DIR?=$(INSTALL_LIB_DIR)/scripts
INSTALL_DATA_DIR?=$(INSTALL_LIB_DIR)/scripts

dist_name=hh-suite-2.2.23

all_static: ffindex_static
	cd src && make all_static

all: ffindex
	cd src && make all

hhblits_static: hhblits_static
	cd src && make hhblits_static

hhblits: ffindex
	cd src && make all

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
	install src/.hhdefaults       $(INSTALL_LIB_DIR)/hhdefaults
	mkdir -p $(INSTALL_DATA_DIR)
	install data/context_data.lib $(INSTALL_DATA_DIR)/context_data.lib
	install data/cs219.lib        $(INSTALL_DATA_DIR)/cs219.lib
	mkdir -p $(INSTALL_SCRIPTS_DIR)
	install dist-scripts/Align.pm        $(INSTALL_SCRIPTS_DIR)/Align.pm
	install dist-scripts/MyPaths.pm      $(INSTALL_SCRIPTS_DIR)/MyPaths.pm
	install dist-scripts/addss.pl        $(INSTALL_SCRIPTS_DIR)/addss.pl
	install dist-scripts/create_cs_db.pl $(INSTALL_SCRIPTS_DIR)/create_cs_db.pl
	install dist-scripts/create_db.pl    $(INSTALL_SCRIPTS_DIR)/create_db.pl
	install dist-scripts/create_profile_from_hhm.pl   $(INSTALL_SCRIPTS_DIR)/create_profile_from_hhm.pl
	install dist-scripts/create_profile_from_hmmer.pl $(INSTALL_SCRIPTS_DIR)/create_profile_from_hmmer.pl
	install dist-scripts/hhmakemodel.pl $(INSTALL_SCRIPTS_DIR)/hhmakemodel.pl
	install dist-scripts/reformat.pl    $(INSTALL_SCRIPTS_DIR)/reformat.pl

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
