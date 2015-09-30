/***** xov.c **************************************
 * Description: Compute rate of crossover from
 *   experimental data.
 * Reference: Kauppi et al. (2009). Analysis of human
 *   recombination products from human sperm.
 *   In: Keeney, S. (ed.), Meiosis, Volume 1,
 *   Molecular and Genetic Methods, Volume 557.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Sep 16 12:06:33 2015
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_cdf.h>
#include "interface.h"
#include "eprintf.h"
#include "ml.h"

void scanFile(Args *args, FILE *fp){
  char *line, *name, *token;
  size_t n;
  int status, len;
  int pos[2], p, i, l, tp;
  Data *data;
  double ll, ul, cf;

  line = NULL;
  name = (char *)emalloc(256);
  while((status = getline(&line,&n,fp)) != -1){
    if(line[0] != '#')
      break;
  }
  if(!args->p)
    printf("# SNP\tbp\t[cM\tcM\tcM]\n");
  data = newData(1);
  while(status != -1){
    data->n = 0;
    name = strtok(line,"\t");     /* read name of region */
    token = strtok(NULL,"\t");    
    len = atoi(token);            /* length of region in bp */
    tp = 0;
    while((token = strtok(NULL,"\t")) != NULL){ /* iterating across experiments */
      if(data->n + 1 == data->maxN)
	expandData(data);
      l = strlen(token);
      p = 0;
      for(i=0;i<l;i++){  /* determine positions of commata */
	if(token[i] == ','){
	  token[i] = '\0';
	  pos[p++] = i+1;
	}
      }
      data->numMol[data->n] = atoi(token); /* number of molecules */
      token += pos[0];
      data->numPos[data->n] = atoi(token);
      token -= pos[0];
      token += pos[1];
      data->numNeg[data->n] = atoi(token);
      token -= pos[1];
      tp += data->numPos[data->n];
      data->n++;
    }
    if(args->p && tp){
      printLikFun(data, args);
      return;
    }
    if(tp)
      cf = crossoverFreq(data);
    else
      cf = 0.0;
    data->mlXover = cf;
    data->ci = args->c;
    if(cf){
      ul = upperCL(data,cf);
      ll = lowerCL(data,cf);
    }else{
      ul = 0.;
      ll = 0.;
    }
    /* convert to centi-Morgans */
    ll *= 100;
    cf *= 100;
    ul *= 100;
    printf("%s\t%*d\t%.3f\t%.3f\t%.3f\n",name,4,len,ll,cf,ul);
    status = getline(&line,&n,fp);
  } 
  free(line);
}

int main(int argc, char *argv[]){
  int i;
  char *version;
  Args *args;
  FILE *fp;

  version = "0.1";
  setprogname2("xov");
  args = getArgs(argc, argv);
  if(args->v)
    printSplash(version);
  if(args->h || args->e)
    printUsage(version);
  if(args->numInputFiles == 0){
    fp = stdin;
    scanFile(args,fp);
  }else{
    for(i=0;i<args->numInputFiles;i++){
      fp = efopen(args->inputFiles[i],"r");
      scanFile(args,fp);
      fclose(fp);
    }
  }
  free(args);
  free(progname());
  return 0;
}
