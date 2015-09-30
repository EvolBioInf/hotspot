/***** interface.c ************************************************
 * Description: Routine for gathering arguments from the command
 *              line.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Jun 20 13:12:10 2004.
 *****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include "interface.h"
#include "eprintf.h"

Args *args;

Args *getArgs(int argc, char *argv[]) {
  char c;
  char *optString = "hvmug:s:l:d:";

  args = (Args *)emalloc(sizeof(Args));
  args->l = DEFAULT_L;
  args->g = NULL;
  args->s = NULL;
  args->d = -1;
  args->u = 0;
  args->h = 0;
  args->v = 0;
  args->e = 0;

  c = getopt(argc, argv, optString);
  while (c != -1) {
    switch (c) {
    case 'u': /* universals */
      args->u = 1;
      args->l = DEFAULT_LU;
      args->d = DEFAULT_D;
      break;
    case 'd': /* distance between universals */
      args->d = atoi(optarg);
      break;
    case 'l': /* oligo length */
      args->l = atoi(optarg);
      break;
    case 'g': /* path to genome files */
      args->g = estrdup(optarg);
      break;
    case 's': /* path to vcf SNP files */
      args->s = estrdup(optarg);
      break;
    case '?': /* fall-through is intentional */
    case 'h': /* print help */
      args->h = 1;
      break;
    case 'v': /* print version */
      args->v = 1;
      break;
    default:
      printf("# unknown argument: %c\n", c);
      args->e = 1;
      return args;
    }
    c = getopt(argc, argv, optString);
  }
  args->inputFiles = argv + optind;
  args->numInputFiles = argc - optind;
  if(!args->g){
    printf("ERROR: Please enter location of genome files using the -g option\n");
    args->e = 1;
  }
  if(!args->s){
    printf("ERROR: Please enter location of SNP files using the -s option\n");
    args->e = 1;
  }
  if(args->d > 0 && !args->u)
    args->u = 1;
  return args;
}

void freeArgs(Args *args){
  free(args->s);
  free(args->g);
  free(args);
}

void printUsage() {
  printf("Usage: %s [options] [inputFiles]\n", progname());
  printf("Purpose: Designing allele-specific oligos for detecting "
	 "recombination hotspots\n");
  printf("Example: aso -g data/mus/genome -s data/mus/vcf "
	 "data/mus/exampleHotSpots19.txt\n");
  printf("Options:\n");
  printf("\t-g STR path to folder containing genome files in FASTA format\n");
  printf("\t-s STR path to folder containing SNPs in VCF format\n");
  printf("\t[-l NUM oligo length; default aso: %d; default universal: %d]\n", DEFAULT_L, DEFAULT_LU);
  printf("\t[-u generate universals; default: allele-specific oligos]\n");
  printf("\t[-d distance between universals; default: %d]\n", DEFAULT_D);
  printf("\t[-h print this help message and exit]\n");
  printf("\t[-v print program information and exit]\n");
  exit(0);
}

void printSplash(char *version) {
  printf("%s %s\n", progname(), version);
  printf("Written by Bernhard Haubold.\n");
  printf("Distributed under the GNU General Public License.\n");
  printf("Please send bug reports to haubold@evolbio.mpg.de\n");
  exit(0);
}
