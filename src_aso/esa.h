/***** esa.h **************************************
 * Description: Header file for computation
 *   of enhance suffix array implemented
 *   in esa.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Jul 15 11:17:08 2013
 **************************************************/
#ifndef ESA
#define ESA

#include <divsufsort.h>
#include "sequenceData.h"
#include "stringUtil.h"

/* define data container */
typedef struct esa{
  long *sa;                /* suffix array */
  long *lcp;               /* longest common prefix array */
  long n;                  /* length of sa and lcp */
} Esa;
/* in esa.c */
Esa *getEsa(Sequence *seq);
void freeEsa();

#endif
