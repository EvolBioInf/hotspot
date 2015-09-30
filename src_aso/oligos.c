/***** oligos.c ***********************************
 * Description:
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Jun  7 11:43:57 2015
 **************************************************/
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include "eprintf.h"
#include "tab.h"
#include "oligos.h"
#include "complexity.h"
#include "data.h"

void snpOligos(Args *args, Sequence *seq, FILE *snp) {
  int i, j, pos, len, offset;
  char *snpId, *intName, *alt, *chr;

  len = args->l;
  intName = estrdup(tabField(INT_FIELD));
  while (tabGetLine(snp) != NULL) {
    chr = tabField(0);
    pos = atoi(tabField(1));
    pos--; /* zero-based counting within progr. */
    snpId = tabField(2);
    alt = tabField(4);
    if (strchr(alt, ','))
      continue;
    /* print first oligo */
    offset = len / 2;
    printf("%s\t%s\t%s\t%d\t", intName, snpId, chr, pos + 1);
    for (i = 0; i < len; i++)
      printf("%c", seq->seq[pos + i - offset]);
    printf("\t");
    /* print second oligo */
    for (i = 0; i < len / 2; i++)
      printf("%c", seq->seq[pos + i - offset]);
    printf("%c", alt[0]);
    for (j = i + 1; j < len; j++)
      printf("%c", seq->seq[pos + j - offset]);
    printf("\n");
  }
  free(intName);
}

void universals(Args *args, Sequence *seq, FILE *snp) {
  int i, snpPos, hotSpotEnd;
  char *intName, *chr;
  double cp;
  long origLen, oligoStart;

  intName = estrdup(tabField(INT_FIELD)); /* this comes from the query on the regions file */
  origLen = seq->len;
  seq->len = args->l;
  oligoStart = getHotSpotStart();
  hotSpotEnd = getHotSpotEnd();
  /* get the fist SNP (if there are any) */
  chr = emalloc(256*sizeof(char));
  snpPos = -1;
  if(fscanf(snp,"%s %d %*s %*s %*s %*s %*s %*s",chr,&snpPos) == EOF)
    snpPos = hotSpotEnd;

  while(oligoStart + args->l < hotSpotEnd){
    while(oligoStart + args->l < snpPos){
      seq->seq += oligoStart;
      cp = complexity(seq);
      if(cp > 1.)
	cp = 1.;
      printf("%s\t%s\t%ld\t%ld\t%.3f\t",intName,chr,oligoStart+1,oligoStart+args->l,cp);
      /* print forward oligo */
      for(i=0;i<args->l;i++)
	printf("%c",seq->seq[i]);
      printf("\t");
      /* print reverse oligo */
      for(i=args->l-1;i>=0;i--){
	if(seq->seq[i] == 'A')
	  printf("T");
	else if(seq->seq[i] == 'C')
	  printf("G");
	else if(seq->seq[i] == 'G')
	  printf("C");
	else if(seq->seq[i] == 'T')
	  printf("A");
      }
      printf("\n");
      seq->seq -= oligoStart;
      oligoStart += (args->d + args->l - 1);
    }
    if(fscanf(snp,"%s %d %*s %*s %*s %*s %*s %*s",chr,&snpPos) == EOF)
      snpPos = hotSpotEnd;
  }
  seq->len = origLen;
  free(chr);
  free(intName);
}
