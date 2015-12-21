all: TMalign TMscore modeller install clean

TMalign.f:
	wget http://zhanglab.ccmb.med.umich.edu/TM-align/TMalign.f

TMalign: TMalign.f
	gfortran -static -O3 -ffast-math -lm -o TMalign TMalign.f

TMscore.f:
	wget http://zhanglab.ccmb.med.umich.edu/TM-score/TMscore.f

TMscore: TMscore.f
	gfortran -static -O3 -ffast-math -lm -o TMscore TMscore.f

modeller-9.13.tar.gz:
	wget http://salilab.org/modeller/9.13/modeller-9.13.tar.gz

modeller: modeller-9.13.tar.gz
	tar xzvf modeller-9.13.tar.gz
	echo -e "3\n$(abspath $(CURDIR)/../bin/modeller9.13/)\nXXXXXXXXXX\n\n\n\n\n\n\n\n\n" > opt_modeller.opt
	mkdir -p $(abspath $(CURDIR)/../bin/modeller9.13/)
	cd modeller-9.13 && ./Install < ../opt_modeller.opt
	rm opt_modeller.opt
	-patch --forward $(abspath $(CURDIR)/../bin/modeller9.13/)/src/include/modeller.h < patch_modeller_h_missing_files.patch
	-cd modeller_cuser && ./buildall.sh $(abspath $(CURDIR)/../bin/modeller9.13/)

install:
	cp TMalign ../bin/
	cp TMscore ../bin/

clean:
	rm -f TMalign*
	rm -f TMscore*
	rm -rf modeller-9.13
	rm -f modeller-9.13.tar.gz
