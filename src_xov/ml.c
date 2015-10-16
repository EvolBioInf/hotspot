/***** ml.c ***************************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Fri Sep 18 14:59:59 2015
 **************************************************/
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_cdf.h>
#include "ml.h"
#include "eprintf.h"

double logLik(double r, void *param){
  double ll,x, k, pnr;
  int i, n;
  Data *data;

  ll = 0;
  data = (Data *)param;
  for(i=0;i<data->n;i++){
    x = data->numMol[i];
    k = data->numNeg[i];
    n = k + data->numPos[i];
    pnr = pow(1.-r,x);
    ll += gsl_sf_lnchoose(n,k) + k*x*log(1-r) + (n-k)*log(1-pnr);
  }
  if(ll == 0.)
    ll = DBL_MAX;
  else
    ll *= -1;
  return ll;
}

double logLikPoi(double lambda, void *param){
  double ll, x, k;
  int i, n;
  Data *data;

  ll = 0;
  data = (Data *)param;
  for(i=0;i<data->n;i++){
    x = data->numMol[i];
    k = data->numNeg[i];
    n = k + data->numPos[i];
    if(lambda > 0.)
      ll += gsl_sf_lnchoose(n,k) - k*x*lambda+(n-k)*log(1-exp(-x*lambda));
  }
  if(ll == 0.)
    ll = DBL_MAX;
  else
    ll *= -1;
  return ll;
}

/* printLikFun: print the likelihood function
 *   NB: cM = crossoverProb * 100
 *   the user thinks in terms of centi-Morgans,
 *   the program computes crossover-Probabilities
 */
void printLikFun(Data *data, Args *args){
  double step, x, ub, l;
  x = args->m/100.;
  ub = args->M/100.;
  step = (ub - x)/(double)args->i;
  while(x <= ub){
    if(x>0){
      if(args->o)
	l = -logLikPoi(x, data);
      else
	l = -logLik(x, data);
    }else{
      if(args->o)
	l = -logLikPoi(DBL_EPSILON, data);
      else
	l = -logLik(DBL_EPSILON, data);
    }
    printf("%.3f\t%.3e\n",x*100.,l);
    x += step;  
  }
}

double crossoverFreq(Data *data){
  int status;
  int iter = 0, max_iter = 100;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
  double m = 0.5;
  double a = 0.0+DBL_EPSILON, b = 0.9;
  double la, lm, lb;
  gsl_function F;

  if(data->o)
    F.function = &logLikPoi;
  else
    F.function = &logLik;
  F.params = data;

  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (T);
  la = F.function(a,data);
  lb = logLik(b,data);
  lm = logLik(m,data);
  /* printf("la: %e; lm: %e; lb: %e\n",la,lm,lb); */
  while(!(lm<la && lm<lb)){
    b = m;
    m = b/2.;
    la = logLik(a,data);
    lb = logLik(b,data);
    lm = logLik(m,data);
  }
  /* printf("la: %e; lm: %e; lb: %e\n",la,lm,lb); */
  gsl_min_fminimizer_set (s, &F, m, a, b);
  do {
      iter++;
      status = gsl_min_fminimizer_iterate (s);

      m = gsl_min_fminimizer_x_minimum (s);
      a = gsl_min_fminimizer_x_lower (s);
      b = gsl_min_fminimizer_x_upper (s);

      status = gsl_min_test_interval (a, b, DBL_MIN, 0.0);
  }  while (status == GSL_CONTINUE && iter < max_iter);
  gsl_min_fminimizer_free(s);
  return (a+b)/2.;
}

Data *newData(int maxN, int poi){
  Data *d;
  d = (Data *)emalloc(sizeof(Data));
  d->n = 0;
  d->numPos = (int *)emalloc(maxN * sizeof(int));
  d->numNeg = (int *)emalloc(maxN * sizeof(int));
  d->numMol = (int *)emalloc(maxN * sizeof(int));
  d->maxN = maxN;
  if(poi)
    d->o = 1;
  else
    d->o = 0;
  return d;
}

void expandData(Data *d){
  d->maxN *= 2;
  d->numPos = (int *)erealloc(d->numPos,d->maxN * sizeof(int));
  d->numNeg = (int *)erealloc(d->numNeg,d->maxN * sizeof(int));
  d->numMol = (int *)erealloc(d->numMol,d->maxN * sizeof(int));
}


double likConf(double lambda, void *param){
  double ll, ml, l, x;
  Data *data;

  data = (Data *)param;
  if(data->o){
    ml = logLikPoi(data->mlXover, param);
    ll = logLikPoi(lambda, param);
  }else{
    ml = logLik(data->mlXover, param);
    ll = logLik(lambda, param);
  }
  x = gsl_cdf_chisq_Pinv(data->ci,1)/2.;
  l = ll - ml - x;

  return l;
}

double upperCL(Data *data, double mlEst){
    int status;
    int iter = 0;
    int max_iter = 100;
    const gsl_root_fsolver_type *solverType;
    gsl_root_fsolver *s;
    gsl_function fun;
    double xLo, xHi;

    
    fun.function = &likConf;
    fun.params = data;

    solverType = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(solverType);
    /* search for upper bound */
    xLo = mlEst;
    if(xLo < 1.0)
      xHi = 1.0;
    else
      return 1.0;
    xHi = 1.0 - DBL_EPSILON;
    status = gsl_root_fsolver_set(s,&fun,xLo,xHi);
    do{
	iter++;
	status = gsl_root_fsolver_iterate(s);
	xLo = gsl_root_fsolver_x_lower(s);
	xHi = gsl_root_fsolver_x_upper(s);
	status = gsl_root_test_interval(xLo, xHi, 0, 0.001);
    }while(status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free(s);
    return (xLo + xHi) / 2.0;
}

double lowerCL(Data *data, double mlEst){
    int status;
    int iter = 0;
    int max_iter = 100;
    const gsl_root_fsolver_type *solverType;
    gsl_root_fsolver *s;
    gsl_function fun;
    double xLo, xHi;

    fun.function = &likConf;
    fun.params = data;

    solverType = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(solverType);
    /* search for lower bound */
    xLo = 0. + DBL_EPSILON;
    xHi = mlEst;
    if(xHi == 0.)
      return 0.;
    status = gsl_root_fsolver_set(s,&fun,xLo,xHi);
    do{
	iter++;
	status = gsl_root_fsolver_iterate(s);
	xLo = gsl_root_fsolver_x_lower(s);
	xHi = gsl_root_fsolver_x_upper(s);
	status = gsl_root_test_interval(xLo, xHi, 0, 0.001);
    }while(status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free(s);

    return (xLo + xHi) / 2.0;
}

void printData(char *name, int len, Data *data){
  int i;

  printf("%s\t%d\t%d,%d,%d",name,len,data->numMol[0],data->numPos[0],data->numNeg[0]);
  for(i=1;i<data->n;i++)
    printf("\t%d,%d,%d",data->numMol[i],data->numPos[i],data->numNeg[i]);
  printf("\n");
}

void freeData(Data *data){
  free(data->numPos);
  free(data->numNeg);
  free(data->numMol);
  free(data);
}
