# hotspot - Software to Support Sperm-Typing for Detecting Recombination Hotspots
Linda Odenthal-Hesse, Julien Dutheil, Fabian Klötzl, Bernhard Haubold

## GETTING STARTED

To build and install all programs execute the following steps. The first step is not required, if you download the latest release. Please note that `libdivsufsort`, `gsl` and `tabix` are required prerequisites.

    autoreconf -i
    ./configure
    make
    make install

For extensive installation instructions see the [hotspotDoc.pdf](hotspotDoc.pdf).

## Usage

Download the example data from [our server](http://guanine.evolbio.mpg.de/hotspot/) or via the commandline.

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
    
## License

Copyright (C) 2015,  Bernhard Haubold

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

## Contact

In case of problems, file a bug or contact us via haubold@evolbio.mpg.de
