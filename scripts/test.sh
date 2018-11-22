./src_aso/aso -g data/mus/genome/ -s data/mus/vcf/ data/mus/exampleHotSpots19.txt > tmp.out 
DIFF=$(diff tmp.out data/test/aso1.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(aso, allele-specific)\tpass\n"
else
    printf "Test(aso, allele-specific)\tfail: %s\n" ${DIFF}
fi

./src_aso/aso -u -l 100 -g data/mus/genome/ -s data/mus/vcf/ data/mus/exampleHotSpots19.txt > tmp.out
DIFF=$(diff tmp.out data/test/aso2.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(aso, universal)\t\tpass\n"
else
    printf "Test(aso, universal)\t\tfail: %s\n" ${DIFF}
fi

./src_asp/asp -g data/mus/genome/ -s data/mus/vcf/ data/mus/exampleHotSpots19.txt > tmp.out 
DIFF=$(diff tmp.out data/test/asp1.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(asp, forward)\t\tpass\n"
else
    printf "Test(asp, forward)\t\tfail: %s\n" ${DIFF}
fi

./src_asp/asp -r -g data/mus/genome/ -s data/mus/vcf/ data/mus/exampleHotSpots19.txt > tmp.out
DIFF=$(diff tmp.out data/test/asp2.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(asp, reverse)\t\tpass\n"
else
    printf "Test(asp, reverse)\t\tfail: %s\n" ${DIFF}
fi

./src_xov/xov data/mus/exampleResults.txt > tmp.out
DIFF=$(diff tmp.out data/test/xov.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(xov)\t\t\tpass\n"
else
    printf "Test(xov)\t\t\tfail: %s\n" ${DIFF}
fi

./src_six/six -s 13 > tmp.out
DIFF=$(diff tmp.out data/test/six.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(six)\t\t\tpass\n"
else
    printf "Test(six)\t\t\tfail: %s\n" ${DIFF}
fi

rm tmp.out
