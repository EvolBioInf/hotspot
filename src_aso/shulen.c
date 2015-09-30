/***** shulen.c ***************************************************
 * Description: Compute the null distribution of shustring lengths
 *              in a query/sbjct comparison given the G/C content
 *              of query and sbjct.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Jun 20 13:03:03 2004.
 * License: GNU General Public
 *****************************************************************/
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_sf_gamma.h>
#include "shulen.h"
#include "interface.h" 

int thresholdReached = 0;

double sum(double x, double p, double l){
  double s;
  double k;

  s = 0.0;
  if(!thresholdReached){
    for(k=0; k<=x; k++){
      s += exp(log(pow(2,x)*pow(p,k)*pow(0.5-p,x-k)*pow(1-pow(p,k)*pow(0.5-p,x-k),l)) + gsl_sf_lnchoose(x,k));
      if(s >= 1.0-DBL_EPSILON){
	thresholdReached = 1.0;
	s = 1.0;
      }
    }
  }else{
    s = 1.0;
  }
  return s;
}

double expShulen(double gc, double l){
  double cp;      /* cumulative probability */
  double p;       /* G/C-content of query */
  double d;       /* maximum divergence */
  double q;       /* G/C-content of sbjct */
  double prob;    /* probability */
  double m;       /* mean shustring length */
  double x;       /* current shustring length */
  double prevP1, curP1;

  p = q = gc;
  cp = 0.0;
  x = 0;
  m = 0.0;
  prevP1 = 0.0;
  d = 1.0 - gc*gc;
  while(cp < 1.0-DBL_EPSILON){
    x++;
    curP1 = sum(x,p/2,l);              /* exact formula */
    curP1 *= 1.0 - pow(1.0 - d,x);
    prob = curP1 - prevP1;             /* exact probability */
    prevP1 = curP1;
    if(prob < 0)
      prob = 0;
    cp += prob;
    m += (double)x * prob;
  }
  thresholdReached = 0;
  return m;
}
