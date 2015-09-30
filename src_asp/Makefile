CC=gcc
CFLAGS= -O3 -Wall -Wshadow -pedantic -g -fsanitize=address -fno-omit-frame-pointer #-m64
# The source files, object files, libraries and executable name.
SRCFILES= hotPrimers.c interface.c eprintf.c sequenceData.c primers.c data.c tab.c stringUtil.c
OBJFILES= hotPrimers.o interface.o eprintf.o sequenceData.o primers.o data.o tab.o stringUtil.o
LIBS= -lm
EXECFILE=hotPrimers
DIRECTORY=HotPrimers
VERSION=0.2
# The make rule for the executable
.PHONY : all
all : $(EXECFILE)
$(EXECFILE) : $(OBJFILES)

	$(CC) $(CFLAGS) -o $(EXECFILE) $(OBJFILES) $(LIBS)
$(EXECFILE).o: interface.h eprintf.h
interface.o: interface.h
eprintf.o: eprintf.h

doc:
	cd ../Doc; make clean; make pdf; cd ../$(DIRECTORY)_$(VERSION)

# Other Standard make rules
lint : 
	lint $(SRCFILES) | more
clean:
	rm -f *.o *~
tarfile:
	cd ../Doc; make clean; make pdf; cd ../$(DIRECTORY)_$(VERSION)
	mkdir $(DIRECTORY)_$(VERSION)
	mkdir $(DIRECTORY)_$(VERSION)/Data/
	mkdir $(DIRECTORY)_$(VERSION)/Data/Genome/
	mkdir $(DIRECTORY)_$(VERSION)/Data/VCF/
	cp ../../Data/Genome/wget.sh $(DIRECTORY)_$(VERSION)/Data/Genome/
	cp ../../Data/VCF/wget.sh $(DIRECTORY)_$(VERSION)/Data/VCF/
	cp ../../Data/Hotspot_regions.txt $(DIRECTORY)_$(VERSION)/Data/hotspotRegions.txt
	cp -rf $(SRCFILES) *.h  README COPYING Makefile ../Doc/hotPrimersDoc.pdf $(DIRECTORY)_$(VERSION)
	tar cvzfh $(EXECFILE)_$(VERSION).tgz $(DIRECTORY)_$(VERSION)
	mv $(EXECFILE)_$(VERSION).tgz ../
	/bin/rm -r $(DIRECTORY)_$(VERSION)

