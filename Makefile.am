SUBDIRS= common src_xov src_aso src_asp src_six

dist_noinst_DATA= COPYING README.md hotspotDoc.pdf

# Helpers for downloading data
.PHONY: get-all-data-mus get-chr19-data-mus data_dirs

data_dirs:
	mkdir -p data/mus/genome
	mkdir -p data/mus/vcf

data/mus/mouseAll.tgz:
	cd data/mus/ && wget 'http://guanine.evolbio.mpg.de/hotspot/mouseAll.tgz'

data/mus/mouse19.tgz:
	cd data/mus/ && wget 'http://guanine.evolbio.mpg.de/hotspot/mouse19.tgz'

get-all-data-mus: data_dirs data/mus/mouseAll.tgz
	cd data/mus/ && tar xzf mouseAll.tgz

get-chr19-data-mus: data_dirs data/mus/mouse19.tgz
	cd data/mus/ && tar xvf mouse19.tgz
