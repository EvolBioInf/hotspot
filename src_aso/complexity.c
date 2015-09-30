/***** complexity.c *******************************
 * Description:
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Jul 15 10:49:56 2015
 **************************************************/
#include <math.h>
#include "eprintf.h"
#include "esa.h"
#include "factor.h"
#include "shulen.h"
#include "sequenceData.h"

/* complexity: compute match complexity */
double complexity(Sequence *seq) {
  long i, *ml;
  double l, c, esl, cObs, cMax, cMin, cNor, gc;
  Esa *esa;
  
  esa = getEsa(seq);
  /* construct and fill array of match lengths */
  ml = (long *)emalloc(esa->n*sizeof(long));
  for(i=0;i<esa->n-1;i++){
    if(esa->lcp[i] < esa->lcp[i+1])
      ml[esa->sa[i]] = esa->lcp[i+1];
    else
      ml[esa->sa[i]] = esa->lcp[i];
  }
  ml[esa->sa[i]] = esa->lcp[i];
  /* colmpute observed number of matches */
  c = 0;
  i = 0;
  while(i < esa->n){
    c++;
    i+= ml[i];
  }
  l = esa->n;
  cObs = (c-1)/l;
  cMin = 2./l;
  gc = gcContent(seq);
  esl = expShulen(gc,seq->len);
  cMax = 1./(esl - 1.);
  cNor = (cObs - cMin)/(cMax - cMin);

  free(ml);
  freeEsa();

  return cNor;
}
