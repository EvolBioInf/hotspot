# hotspot - Software to Support Sperm-Typing for Detecting Recombination Hotspots
Bernhard Haubold & Fabian Klötzl

## GETTING STARTED

To build and install all programs execute the following steps.

    autoreconf -i
    ./configure
    make
    make install

## Usage

Download the example data.

    make get-chr19-data-mus

### Executing `asp`
    asp -g ./data/mus/genome/ -s ./data/mus/vcf/ ./data/mus/exampleHotSpots19.txt
### Executing `aso`
    aso -g ./data/mus/genome/ -s ./data/mus/vcf/ ./data/mus/exampleHotSpots19.txt
### Executing Universals in `aso`
    aso -u -g ./data/mus/genome/ -s ./data/mus/vcf/ ./data/mus/exampleHotSpots19.txt
### Executing `xov`
    xov ./data/mus/exampleResults.txt
### Executing `six`
    six
### Testing `xov` using `six`
    six | xov
