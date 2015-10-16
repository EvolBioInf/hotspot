CFLAGS=  -O3 -W -Wall -g #-fsanitize=address
CPPFLAGS= -I common
LIBS= -lm -ldivsufsort -lgsl -lblas

all: asp aso xov six
asp aso xov six: common/eprintf.o common/sequenceData.o common/stringUtil.o common/tab.o

aso: LIBS+= -lz
aso: src_aso/aso.o src_aso/complexity.o src_aso/esa.o src_aso/data.o src_aso/interface.o src_aso/oligos.o src_aso/shulen.o
	$(CC) $^ -o $@ $(LIBS)

asp: src_asp/asp.o src_asp/data.o src_asp/interface.o src_asp/primers.o src_asp/asPrimers.o src_asp/univPrimers.o
	$(CC) $^ -o $@

xov: src_xov/xov.o src_xov/interface.o src_xov/ml.o
	$(CC) $^ -o $@ $(LIBS)

six: src_six/six.o src_six/interface.o src_six/gsl_rng.o
	$(CC) $^ -o $@ $(LIBS)

src_asp/%.o: CPPFLAGS+= -Isrc_asp -I.
src_aso/%.o: CPPFLAGS+= -Isrc_aso


.PHONY: clean
clean:
	rm -rf *.o asp aso xov six src_asp/*.o src_aso/*.o src_xov/*.o src_six/*.o common/*.o


# Helpers for downloading data
.PHONY: get-data-mus

get-all-data-mus:
	mkdir -p data/mus/genome;
	mkdir -p data/mus/vcf;
	cd data/mus/genome && wget 'ftp://ftp.ncbi.nlm.nih.gov/genomes/M_musculus/Assembled_chromosomes/seq/*_ref_*.fa.gz'
	cd data/mus/genome && gunzip *.gz
	cd data/mus/vcf && wget 'ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/*.*'
get-chr19-data-mus:
	mkdir -p data/mus/genome;
	mkdir -p data/mus/vcf;
	cd data/mus/genome && wget 'ftp://ftp.ncbi.nlm.nih.gov/genomes/M_musculus/Assembled_chromosomes/seq/*_ref_*chr19.fa.gz'
	cd data/mus/genome && gunzip *.gz
	cd data/mus/vcf && wget 'ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/mouse_10090/VCF/*chr_19*.*'
