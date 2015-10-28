SUBDIRS= common src_xov src_aso src_asp src_six

# Helpers for downloading data
.PHONY: get-all-data-mus get-chr19-data-mus

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
