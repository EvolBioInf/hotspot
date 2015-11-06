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
#include <ctype.h>
#include <gsl/gsl_cdf.h>
#include "interface.h"
#include "eprintf.h"
#include "ml.h"
#include "hotspot.h"

/* intVal: return integer value of string */
int intVal(char *str){
  int l, i, j;
  char *num;

  num = (char *)emalloc(256);
  l = strlen(str);
  j = 0;
  for(i=0;i<l;i++)
    if(isdigit(str[i]))
      num[j++] = str[i];
  num[j] = '\0';
  
  free(num);
  return atoi(num);
}

void scanFile(Args *args, FILE *fp){
  char *line, *name, *token, *chr;
  size_t n;
  int status, len, l, start, end;
  int pos[2], p, i, tp, tn;
  Data *data;
  double ll, ul, cf;

  line = NULL;
  while((status = getline(&line,&n,fp)) != -1){
    if(line[0] != '#')
      break;
  }
  if(!args->p)
    printf("# Int\tChr\tStart\tEnd\tLen\t[%%\t%%\t%%]\t[cM/Mb\tcM/Mb\tcM/Mb]\n");
  data = newData(1,args->o);
  while(status != -1){
    data->n = 0;
    name = estrdup(strtok(line,"\t")); /* read name of region */
    chr = estrdup(strtok(NULL,"\t"));  /* read chromosome */
    start = intVal(strtok(NULL,"\t")); /* start position */
    end = intVal(strtok(NULL,"\t"));   /* end position */
    len = end - start + 1;             /* length of region in bp */
    tp = 0; /* total positives */
    tn = 0; /* total negatives */
    while((token = strtok(NULL,"\t")) != NULL){ /* iterating across experiments */
      if(data->n + 1 == data->maxN)
	expandData(data);
      l = strlen(token);
      p = 0;
      for(i=0;i<l;i++){  /* determine positions of commata */
	if(token[i] == '|'){
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
      tn += data->numNeg[data->n];
      data->n++;
    }
    if(args->d){
      printData(name,len,data);
      exit(0);
    }
    if(args->p && tp){
      printLikFun(data, args);
      return;
    }
    cf = crossoverFreq(data);
    data->mlXover = cf;
    data->ci = args->c;
    ul = upperCL(data,cf);
    ll = lowerCL(data,cf);
    /* convert to centi-Morgans */
    ll *= 100;
    cf *= 100;
    ul *= 100;
    if(cf < 100){
      printf("%s\t%s\t%d\t%d\t%*d\t%.3f\t%.3f\t%.3f",name,chr,start,end,4,len,ll,cf,ul);
      /* convert to cM per Mb */
      ll = ll / len * 1000000;
      cf = cf / len * 1000000;
      ul = ul / len * 1000000;
      printf("\t%.3f\t%.3f\t%.3f\n",ll,cf,ul);
    }else
      printf("%s\t%d\t%d\t%*d\tNA\tNA\tNA\tNA\tNA\tNA\n",name,start,end,4,len);
    status = getline(&line,&n,fp);
    free(name);
  } 
  free(line);
  freeData(data);
}

int main(int argc, char *argv[]){
  int i;
  char *version;
  Args *args;
  FILE *fp;

  version = VERSION;
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
