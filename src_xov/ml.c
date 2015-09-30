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

double logLik(double lambda, void *param){
  double ll, x, k;
  int i, n;
  Data *data;

  ll = 0;
  data = (Data *)param;
  for(i=0;i<data->n;i++){
    x = data->numMol[i];
    k = data->numPos[i];
    n = k + data->numNeg[i];
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
    if(x>0)
      l = -logLik(x, data);
    else
      l = -logLik(DBL_EPSILON, data);
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
  double a = 0.0, b = 1.0;
  gsl_function F;

  F.function = &logLik;
  F.params = data;

  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (T);
  gsl_min_fminimizer_set (s, &F, m, a, b);

  do {
      iter++;
      status = gsl_min_fminimizer_iterate (s);

      m = gsl_min_fminimizer_x_minimum (s);
      a = gsl_min_fminimizer_x_lower (s);
      b = gsl_min_fminimizer_x_upper (s);

      status = gsl_min_test_interval (a, b, 0.00000001, 0.0);
    }  while (status == GSL_CONTINUE && iter < max_iter);
  return (a+b)/2.;
}

Data *newData(int maxN){
  Data *d;
  d = (Data *)emalloc(sizeof(Data));
  d->n = 0;
  d->numPos = (int *)emalloc(maxN * sizeof(int));
  d->numNeg = (int *)emalloc(maxN * sizeof(int));
  d->numMol = (int *)emalloc(maxN * sizeof(int));
  return d;
}

void expandData(Data *d){
  d->maxN *= 2;
  d->numPos = (int *)emalloc(d->maxN * sizeof(int));
  d->numNeg = (int *)emalloc(d->maxN * sizeof(int));
  d->numMol = (int *)emalloc(d->maxN * sizeof(int));
}


double likConf(double lambda, void *param){
  double ll, ml, l, x;
  Data *data;

  data = (Data *)param;
  ml = logLik(data->mlXover, param);
  ll = logLik(lambda, param);
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
    xHi = 1.0;
    status = gsl_root_fsolver_set(s,&fun,xLo,xHi);
    do{
	iter++;
	status = gsl_root_fsolver_iterate(s);
	xLo = gsl_root_fsolver_x_lower(s);
	xHi = gsl_root_fsolver_x_upper(s);
	status = gsl_root_test_interval(xLo, xHi, 0, 0.0000001);
    }while(status == GSL_CONTINUE && iter < max_iter);
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
    xLo = 0.;
    xHi = mlEst;
    status = gsl_root_fsolver_set(s,&fun,xLo,xHi);
    do{
	iter++;
	status = gsl_root_fsolver_iterate(s);
	xLo = gsl_root_fsolver_x_lower(s);
	xHi = gsl_root_fsolver_x_upper(s);
	status = gsl_root_test_interval(xLo, xHi, 0, DBL_EPSILON);
    }while(status == GSL_CONTINUE && iter < max_iter);
    return (xLo + xHi) / 2.0;
}
