/***** oligos.h ***********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Jun  7 11:44:00 2015
 **************************************************/
#ifndef OLIGOS
#define OLIGOS

#include <stdio.h>
#include "sequenceData.h"
#include "interface.h"

void snpOligos(Args *args, Sequence *seq, FILE *snp);
void universals(Args *args, Sequence *seq, FILE *snp);

#endif

