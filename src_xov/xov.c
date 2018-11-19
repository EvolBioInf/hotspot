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
#include "tab.h"
#include "interface.h"
#include "eprintf.h"
#include "ml.h"
#include "config.h"

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
  i = atoi(num);
  free(num);

  return i;
}

void scanFile(Args *args, FILE *fp){
  char *line, *name, *token, *chr;
  int len, l, start, end;
  int pos[2], p, i, j, tp, tn;
  Data *data;
  double ll, ul, cf;

  if(!args->p)
    printf("# Int\tChr\tStart\tEnd\tLen\t[%%\t%%\t%%]\t[cM/Mb\tcM/Mb\tcM/Mb]\n");
  data = newData(1,args->o);
  tabSetFieldSep(args->D);
  while((line = tabGetLine(fp)) != NULL){
    if(line[0] == '#' || strlen(line) == 0)
      continue;
    if(tabNfield() < 5) {
      fprintf(stderr, "ERROR[xov]: There is a problem scanning the columns in your input file.\n");
      fprintf(stderr, "ERROR[xov]: Please check the file formatting and consider adjusting the delimiter using -D.\n");
      exit(-1);
    }
    data->n = 0;
    name = estrdup(tabField(0)); /* read name of region */
    chr = estrdup(tabField(1));  /* read chromosome */
    start = intVal(tabField(2)); /* start position */
    end = intVal(tabField(3));   /* end position */
    len = end - start + 1;          /* length of region in bp */
    tp = 0; /* total positives */
    tn = 0; /* total negatives */
    for(i = 4; i < tabNfield(); i++) {
      if(data->n + 1 == data->maxN)
	expandData(data);
      token = tabField(i);
      l = strlen(token);
      p = 0;
      for(j=0;j<l;j++){  /* determine positions of commata */
	if(token[j] == '|'){
	  token[j] = '\0';
	  pos[p++] = j+1;
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
