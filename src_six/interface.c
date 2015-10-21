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
  char *optString = "hvr:m:n:x:s:";
  int i;

  args = (Args *)emalloc(sizeof(Args));
  args->m = DEFAULT_M;
  args->n = DEFAULT_N;
  args->x = DEFAULT_X;
  args->r = DEFAULT_R;
  args->s = 0;
  args->mol = (int *)emalloc(DEFAULT_M*sizeof(int));
  args->mol[0] = DEFAULT_MM0;
  args->mol[1] = DEFAULT_MM1;
  args->mol[2] = DEFAULT_MM2;
  args->h = 0;
  args->v = 0;
  args->e = 0;

  c = getopt(argc, argv, optString);
  while(c != -1){
    switch(c){
    case 'm':                           /* number of experiments */
      args->m = atoi(optarg);
      args->mol = (int *)erealloc(args->mol,sizeof(int)*args->m);
      for(i=0;i<args->m;i++)
	args->mol[i] = atoi(argv[optind+i]);
      break;
    case 's':                           /* seed for random number generator */
      args->s = atoi(optarg);
      break;
    case 'n':                           /* replicates per experiment */
      args->n = atoi(optarg);
      break;
    case 'r':                           /* replicate experiments */
      args->r = atoi(optarg);
      break;
    case 'x':                           /* crossover rate */
      args->x = atof(optarg);
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
  printf("Purpose: Simulate crossover data for analysis with xov\n");
  printf("Example: six\n");
  printf("Options:\n");
  printf("\t[-r <NUM> number of replicates; default: %d\n",DEFAULT_R);
  printf("\t[-m <NUM> number of experiments per replicate, followed by number of molecules per exp.; default: %d %d %d %d]\n", DEFAULT_M,DEFAULT_MM0,DEFAULT_MM1,DEFAULT_MM2);
  printf("\t[-n <NUM> number of replicates per experiment; default: %d]\n",DEFAULT_N);
  printf("\t[-x <NUM> crossover rate in cM; default: %.1f\n",DEFAULT_X);
  printf("\t[-s <NUM> seed for random number generator; default: generated internally]\n");
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
