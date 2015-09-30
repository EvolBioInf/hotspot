/***** primers.h **********************************
 * Description:
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Jun  7 11:44:04 2015
 **************************************************/
#ifndef PRIMERS
#define PRIMERS

#include <stdio.h>
#include "sequenceData.h"
#include "interface.h"

void forwardP(Args *args, Sequence *seq, FILE *snp);
void reverseP(Args *args, Sequence *seq, FILE *snp);
void forwardPuniv(Args *args, Sequence *seq, FILE *snp, int start, int end);
void reversePuniv(Args *args, Sequence *seq, FILE *snp, int start, int end);
int optimalStartPos(Args *args, char *primer, int minStart);
float gc(const char *primer, size_t startPos);
void complement(char *primer);
void reverse(char *primer);
double tm(char *oligo, size_t startPos);
int intersects(int *pos, int n, int s, int e);

#endif
