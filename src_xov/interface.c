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

Args *getArgs(int argc, char *argv[]){
  char c;
  char *optString = "hvpdoc:m:M:i:";

  args = (Args *)emalloc(sizeof(Args));
  args->c = DEFAULT_C;
  args->m = DEFAULT_M;
  args->M = DEFAULT_MM;
  args->i = DEFAULT_I;
  args->o = 0;
  args->d = 0;
  args->p = 0;
  args->h = 0;
  args->v = 0;
  args->e = 0;

  c = getopt(argc, argv, optString);
  while(c != -1){
    switch(c){
    case 'c':
      args->c = atof(optarg);
      break;
    case 'm': /* minimum x-value for log-likelihood curve */
      args->m = atof(optarg);
      break;
    case 'M': /* maximum x-value for log-likelihood curve */
      args->M = atof(optarg);
      break;
    case 'p': /* print log-likelihood curve? */
      args->p = 1;
      break;
    case 'o': /* use Poisson approximation */
      args->o = 1;
      break;
    case 'd': /* print data? */
      args->d = 1;
      break;
    case 'i': /* number of intervals for log-likelihood curve */
      args->i = atoi(optarg);
      break;
    case '?':                           /* fall-through is intentional */
    case 'h':                           /* print help */
      args->h = 1;
      break;
    case 'v':                           /* print version */
      args->v = 1;
      break;
    default:
      printf("# unknown argument: %c\n",c);
      args->e = 1;
      return args;
    }
    c = getopt(argc, argv, optString);
  }
  args->inputFiles = argv + optind;
  args->numInputFiles = argc - optind;
  return args;
}


void printUsage(){
  printf("Usage: %s [options] [inputFiles]\n",progname());
  printf("Purpose: Compute rate of crossover from experimental data\n");
  printf("Example: xov foo.dat\n");
  printf("Options:\n");
  printf("\t[-c <NUM> confidence interval; default: %.2f]\n",DEFAULT_C);
  printf("\t[-p print log-likelihood curve; default: print ML-estimates]\n");
  printf("\t[-m <NUM> minimum x-value for log-likelihood curve; default %.3f]\n",DEFAULT_M);
  printf("\t[-M <NUM> maximum x-value for log-likelihood curve; default: %.3f]\n",DEFAULT_MM);
  printf("\t[-i <NUM> number of intervals for log-likelihood curve; default: %d]\n",DEFAULT_I);
  printf("\t[-o use Poisson approximation; default: use exact likelihood function]\n");
  printf("\t[-d print data and exit (debug)]\n");
  printf("\t[-h print this help message and exit]\n");
  printf("\t[-v print program information and exit]\n");
  exit(0);
}

void printSplash(char *version){
  printf("%s %s\n",progname(),version);
  printf("Written by Bernhard Haubold.\n");
  printf("Distributed under the GNU General Public License.\n");
  printf("Please send bug reports to haubold@evolbio.mpg.de\n");
  exit(0);
}
