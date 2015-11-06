/***** six.c **************************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sat Oct 10 12:30:06 2015
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "gsl_rng.h"
#include "interface.h"
#include "eprintf.h"
#include "hotspot.h"

void runSim(Args *args, gsl_rng *r){
  int i, j, k, pos;
  double p;
  p = (double)args->x / 100.;
  for(i=0;i<args->m;i++){
    pos = 0;
    for(j=0;j<args->n;j++){
      for(k=0;k<args->mol[i];k++)
	if(gsl_rng_uniform(r) < p){
	  pos++;
	  break;
	}
    }
    printf("\t%d|%d|%d",args->mol[i],pos,args->n-pos);
  }
  printf("\n");
}

int main(int argc, char *argv[]){
  int i;
  char *version;
  Args *args;
  gsl_rng *r;

  version = VERSION;
  setprogname2("six");
  args = getArgs(argc, argv);
  if(args->v)
    printSplash(version);
  if(args->h || args->e)
    printUsage(version);
  r = ini_gsl_rng(args);
  printf("# Int\tChr\tStart\tEnd");
  for(i=0;i<args->m;i++)
    printf("\td|n-k|k");
  printf("\n");
  for(i=0;i<args->r;i++){
    printf("Int%d\t1\t15,000,001\t15,002,000",i+1); 
    runSim(args,r);
  }
  free_gsl_rng(r,args);
  free(args);
  free(progname());
  return 0;
}
