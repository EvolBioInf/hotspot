/***** complexity.c *******************************
 * Description:
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Jul 15 10:49:56 2015
 **************************************************/
#include <math.h>
#include <limits.h>
#include "eprintf.h"
#include "esa.h"
#include "factor.h"
#include "shulen.h"
#include "sequenceData.h"

long max3(long x1, long x2, long x3){
  long m;

  m = LONG_MIN;
  if(m < x1)
    m = x1;
  if(m < x2)
    m = x2;
  if(m < x3)
    m = x3;

  return m;
}

/* complexity: compute match complexity from the number of match factors*/
double complexity(Sequence *seq) {
  long i, l1, l2;
  double l, c, esl, cObs, cMax, cMin, cNor, gc;
  Esa *esa;
  
  esa = getEsa(seq);
  /* prevent out of bounds error */
  esa->isa = (long *)erealloc(esa->isa,(esa->n+1)*sizeof(long));
  esa->isa[esa->n] = 0;
  /* count match factors */
  c = 0;
  i = 0;
  while(i < esa->n){
    l1 = esa->lcp[esa->isa[i]];
    l2 = esa->lcp[esa->isa[i]+1];
    i += max3(l1,l2,1);
    c++;
  }
  l = esa->n;
  cObs = (c-1)/l;
  cMin = 2./l;
  gc = gcContent(seq);
  esl = expShulen(gc,seq->len);
  cMax = 1./(esl - 1.);
  cNor = (cObs - cMin)/(cMax - cMin);

  freeEsa();

  return cNor;
}


