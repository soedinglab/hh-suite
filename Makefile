# This can be overridden e.g.: make install INSTALL_DIR=...
INSTALL_DIR?=$(PWD)

# Overriding this is currently not fully supported as the code won't know
# to what this is set then. You can try setting HHLIB.
INSTALL_BIN_DIR?=$(INSTALL_DIR)/bin
INSTALL_LIB_DIR?=$(INSTALL_DIR)/lib
INSTALL_SCRIPTS_DIR?=$(INSTALL_DIR)/scripts
INSTALL_DATA_DIR?=$(INSTALL_DIR)/data

dist_name=hhsuite-2.0.16

#all_static: ffindex_static
#	$(MAKE) -C src all_static

all: ffindex
	$(MAKE) -C src all

doc:
	$(MAKE) -C src hhsuite-userguide.pdf

hhblits_static: hhblits_static
	$(MAKE) -C src hhblits_static

hhblits: ffindex
	$(MAKE) -C src all

ffindex:
	$(MAKE) -C lib/ffindex

#ffindex_static:
#	$(MAKE) -C lib/ffindex FFINDEX_STATIC=1

install:
	$(MAKE) -C lib/ffindex install INSTALL_DIR=$(INSTALL_DIR)
	mkdir -p $(INSTALL_BIN_DIR)
	install src/hhblits                           $(INSTALL_BIN_DIR)/hhblits
	install src/hhblits_omp                       $(INSTALL_BIN_DIR)/hhblits_omp
	install src/hhalign                           $(INSTALL_BIN_DIR)/hhalign
	install src/hhconsensus                       $(INSTALL_BIN_DIR)/hhconsensus
	install src/hhfilter                          $(INSTALL_BIN_DIR)/hhfilter
	install src/hhmake                            $(INSTALL_BIN_DIR)/hhmake
	install src/hhsearch                          $(INSTALL_BIN_DIR)/hhsearch
	install src/ffindex_from_fasta                $(INSTALL_BIN_DIR)/ffindex_from_fasta
	install src/cstranslate                       $(INSTALL_BIN_DIR)/cstranslate
	mkdir -p $(INSTALL_LIB_DIR)
	mkdir -p $(INSTALL_DATA_DIR)
	ifneq ($(INSTALL_DATA_DIR),data)
		install -m 0644 data/context_data.lib         $(INSTALL_DATA_DIR)/context_data.lib
		install -m 0644 data/context_data.crf         $(INSTALL_DATA_DIR)/context_data.crf
		install -m 0644 data/cs219.lib                $(INSTALL_DATA_DIR)/cs219.lib
		install -m 0644 data/do_not_delete            $(INSTALL_DATA_DIR)/do_not_delete
		install -m 0644 data/do_not_delete.phr        $(INSTALL_DATA_DIR)/do_not_delete.phr
		install -m 0644 data/do_not_delete.pin        $(INSTALL_DATA_DIR)/do_not_delete.pin
		install -m 0644 data/do_not_delete.psq        $(INSTALL_DATA_DIR)/do_not_delete.psq
	endif
	mkdir -p $(INSTALL_SCRIPTS_DIR)
	ifneq ($(INSTALL_SCRIPTS_DIR),scripts)
		install -m 0644 scripts/Align.pm              $(INSTALL_SCRIPTS_DIR)/Align.pm
		install -m 0644 scripts/HHPaths.pm            $(INSTALL_SCRIPTS_DIR)/HHPaths.pm
		install scripts/splitfasta.pl                 $(INSTALL_SCRIPTS_DIR)/splitfasta.pl
		install scripts/addss.pl                      $(INSTALL_SCRIPTS_DIR)/addss.pl
		install scripts/create_profile_from_hhm.pl    $(INSTALL_SCRIPTS_DIR)/create_profile_from_hhm.pl
		install scripts/create_profile_from_hmmer.pl  $(INSTALL_SCRIPTS_DIR)/create_profile_from_hmmer.pl
		install scripts/hhmakemodel.pl                $(INSTALL_SCRIPTS_DIR)/hhmakemodel.pl
		install scripts/reformat.pl                   $(INSTALL_SCRIPTS_DIR)/reformat.pl
		install scripts/multithread.pl                $(INSTALL_SCRIPTS_DIR)/multithread.pl
		install scripts/hhsuitedb.pl                  $(INSTALL_SCRIPTS_DIR)/hhsuitedb.pl
		install scripts/checkA3M.pl                   $(INSTALL_SCRIPTS_DIR)/checkA3M.pl
		install scripts/pdb2fasta.pl                  $(INSTALL_SCRIPTS_DIR)/pdb2fasta.pl
		install scripts/pdbfilter.pl                  $(INSTALL_SCRIPTS_DIR)/pdbfilter.pl
	endif

deinstall:
	$(MAKE) -C lib/ffindex deinstall INSTALL_DIR=$(INSTALL_DIR)
	rm -f $(INSTALL_BIN_DIR)/hhblits 
	rm -f $(INSTALL_BIN_DIR)/hhblits_omp
	rm -f $(INSTALL_BIN_DIR)/hhalign
	rm -f $(INSTALL_BIN_DIR)/hhconsensus 
	rm -f $(INSTALL_BIN_DIR)/hhfilter
	rm -f $(INSTALL_BIN_DIR)/hhmake
	rm -f $(INSTALL_BIN_DIR)/hhsearch
	rm -f $(INSTALL_BIN_DIR)/ffindex_from_fasta
	rm -f $(INSTALL_BIN_DIR)/cstranslate
	
	rm -f $(INSTALL_DATA_DIR)/context_data.lib
	rm -f $(INSTALL_DATA_DIR)/context_data.crf 
	rm -f $(INSTALL_DATA_DIR)/cs219.lib
	rm -f $(INSTALL_DATA_DIR)/do_not_delete
	rm -f $(INSTALL_DATA_DIR)/do_not_delete.phr
	rm -f $(INSTALL_DATA_DIR)/do_not_delete.pin
	rm -f $(INSTALL_DATA_DIR)/do_not_delete.psq
	
	rm -f $(INSTALL_SCRIPTS_DIR)/Align.pm
	rm -f $(INSTALL_SCRIPTS_DIR)/HHPaths.pm
	rm -f $(INSTALL_SCRIPTS_DIR)/splitfasta.pl
	rm -f $(INSTALL_SCRIPTS_DIR)/addss.pl
	rm -f $(INSTALL_SCRIPTS_DIR)/create_profile_from_hhm.pl
	rm -f $(INSTALL_SCRIPTS_DIR)/create_profile_from_hmmer.pl
	rm -f $(INSTALL_SCRIPTS_DIR)/hhmakemodel.pl
	rm -f $(INSTALL_SCRIPTS_DIR)/reformat.pl
	rm -f $(INSTALL_SCRIPTS_DIR)/multithread.pl
	rm -f $(INSTALL_SCRIPTS_DIR)/hhsuitedb.pl
	rm -f $(INSTALL_SCRIPTS_DIR)/checkA3M.pl
	rm -f $(INSTALL_SCRIPTS_DIR)/pdb2fasta.pl
	rm -f $(INSTALL_SCRIPTS_DIR)/pdbfilter.pl
	
	rmdir $(INSTALL_BIN_DIR) || true
	rmdir $(INSTALL_DATA_DIR) || true
	rmdir $(INSTALL_SCRIPTS_DIR) || true
	rmdir $(INSTALL_LIB_DIR) || true

.PHONY: clean
clean:
	cd lib/ffindex && $(MAKE) clean
	$(MAKE) -C src clean

dist/$(dist_name).tar.gz:
	make clean
	mkdir -p dist
	git archive --prefix=$(dist_name)/ -o dist/$(dist_name).tar.gz HEAD
	cd dist && tar xf $(dist_name).tar.gz
	mkdir -p dist/$(dist_name)/bin
	cd dist/$(dist_name) && rsync --exclude .git --exclude .hg -av ../../lib .
	cd dist && tar czf $(dist_name).tar.gz $(dist_name)
