/***** interface.c ************************************************
 * Description: Routine for gathering arguments from the command
 *              line.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Jun 20 13:12:10 2004.
 *****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "interface.h"
#include "eprintf.h"

Args *args;

Args *getArgs(int argc, char *argv[]) {
  char c;
  char *optString = "hvubm:M:o:rg:s:f:d:D:";
  short minLenIsSet;
  short maxLenIsSet;

  minLenIsSet = 0;
  maxLenIsSet = 0;
  args = (Args *)emalloc(sizeof(Args));
  args->f = DEFAULT_F;
  args->m = DEFAULT_MI;
  args->M = DEFAULT_MA;
  args->g = NULL;
  args->s = NULL;
  args->d = -1;
  args->D = -1;
  args->u = 0;
  args->o = DEFAULT_O;
  args->r = 0;
  args->b = 0;
  args->h = 0;
  args->v = 0;
  args->e = 0;

  while ((c = getopt(argc, argv, optString)) != -1) {
    switch (c) {
    case 'u': /* universals? */
      args->u = 1;
      break;
    case 'd': /* distance between forward universals */
      args->d = atoi(optarg);
      break;
    case 'D': /* distance between reverse universals */
      args->D = atoi(optarg);
      break;
    case 'f': /* length of flanking region */
      args->f = atoi(optarg);
      break;
    case 'm': /* minimum length */
      args->m = atoi(optarg);
      minLenIsSet = 1;
      break;
    case 'b': /* print debugging information? */
      args->b = 1;
      break;
    case 'M': /* maximum length */
      args->M = atoi(optarg);
      maxLenIsSet = 1;
      break;
    case 'r': /* reverse primer */
      args->r = 1;
      break;
    case 'o': /* optimal GC content */
      args->o = atof(optarg);
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
  }
  if(!args->g){
    printf("ERROR: Please enter the name of the folder containing the chromosome sequences using -g\n");
    args->e = 1;
  }
  if(!args->s){
    printf("ERROR: Please enter the name of the folder containing the SNP data using -s\n");
    args->e = 1;
  }
  if(args->u){
    if(args->d == -1)
      args->d = DEFAULT_DF;
    if(args->D == -1)
      args->D = DEFAULT_DR;
    if(!minLenIsSet)
      args->m = DEFAULT_UI;
    if(!maxLenIsSet)
      args->M = DEFAULT_UA;
  }
  args->inputFiles = argv + optind;
  args->numInputFiles = argc - optind;
  return args;
}

void freeArgs(Args *a) {
  free(a->g);
  free(a->s);
  free(a);
  a = NULL;
}

void printUsage() {
  printf("Usage: %s -g genomeDirectory -s snpDirectory [options] [inputFiles]\n", progname());
  printf("Purpose: Designing allele-specific primers for detecting recombination hotspots\n");
  printf("Example: asp -g data/mus/genome -s data/mus/vcf data/mus/exampleHotSpots19.txt\n");
  printf("Options:\n");
  printf("\t-g STR path to folder containing genome files\n");
  printf("\t-s STR path to folder containing SNPs in VCF format\n");
  printf("\t[-m NUM minimum primer length; default: %d]\n", DEFAULT_MI);
  printf("\t[-M NUM maximum primer length; default: %d]\n", DEFAULT_MA);
  printf("\t[-o NUM optimal GC content; default: %.2f]\n", DEFAULT_O);
  printf("\t[-f NUM length of flanking region; default: %d]\n", DEFAULT_F);
  printf("\t[-r reverse primer; default: forward]\n");
  printf("\t[-u universals of default lengths %d-%d; default: allele-specific primers]\n",DEFAULT_UI,DEFAULT_UA);
  printf("\t[-d NUM distance between forward universals; default: %d]\n",DEFAULT_DF);
  printf("\t[-D NUM distance between reverse universals; default: %d]\n",DEFAULT_DR);
  printf("\t[-b print debugging information]\n");
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
