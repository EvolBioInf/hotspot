/***** asPrimers.c ********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Sep 21 11:27:46 2015
 **************************************************/
#include <stdio.h>
#include <string.h>
#include "data.h"
#include "tab.h"
#include "primers.h"
#include "eprintf.h"

void forwardP(Args *args, Sequence *seq, FILE *snp) {
  int i, pos, maxStart, optStart;
  char *snpId, *alt, *hsName, *chr, *primer;

  hsName = (char *)emalloc(256);
  primer = (char *)emalloc(args->M + 1);
  alt = (char *)emalloc(256);
  snpId = (char *)emalloc(256);
  chr = (char *)emalloc(256);
  /* get name of hotspot */
  strcpy(hsName, tabField(INT_FIELD));
  while(fscanf(snp,"%s %d %s %*s %s %*s %*s %*s",chr,&pos,snpId,alt) != EOF){
    pos--; /* zero-based counting within progr. */
    maxStart = pos - args->M + 1;
    if (strchr(alt, ','))
      continue;
    for (i = 0; i < args->M; i++)
      primer[i] = seq->seq[maxStart + i];
    primer[i] = '\0';
    /* print first primer */
    optStart = optimalStartPos(args, primer, args->M - args->m);
    printf("%s\t%s\t%s\t%d\t", hsName, snpId, chr, pos + 1);
    for (i = optStart; i < args->M; i++)
      printf("%c", primer[i]);
    printf("\t%.2f\t%.1f\t", gc(primer, optStart), tm(primer, optStart));
    primer[args->M - 1] = alt[0];
    /* print second primer */
    optStart = optimalStartPos(args, primer, args->M - args->m);
    for (i = optStart; i < args->M; i++)
      printf("%c", primer[i]);
    printf("\t%.2f\t%.1f\n", gc(primer, optStart), tm(primer, optStart));
  }
  free(hsName);
  free(primer);
  free(alt);
  free(snpId);
  free(chr);
}

void reverseP(Args *args, Sequence *seq, FILE *snp) {
  int i, pos, optStart;
  char *snpId, *alt, *hsName, *chr, *primer;

  hsName = emalloc(256);
  chr = emalloc(256);
  snpId = emalloc(256);
  alt = emalloc(256);
  
  primer = emalloc(args->M + 1);
  /* get name of hotspot */
  strcpy(hsName, tabField(0));


  while(fscanf(snp,"%s %d %s %*s %s %*s %*s %*s",chr,&pos,snpId,alt) != EOF){
    /* zero-based counting within progr. */
    pos--;
    if (strchr(alt, ',')) 
      continue;
    for (i = 0; i < args->M; i++)
      primer[i] = seq->seq[pos + i];

    primer[args->M] = '\0';
    reverse(primer);
    complement(primer);

    /* print first primer */
    optStart = optimalStartPos(args, primer, args->M - args->m);
    printf("%s\t%s\t%s\t%d\t", hsName, snpId, chr, pos + 1);

    for (i = optStart; i < args->M; i++)
      printf("%c", primer[i]);

    printf("\t%.2f\t%.1f\t", gc(primer, optStart), tm(primer, optStart));
    complement(alt);
    primer[args->M - 1] = alt[0];
		
    /* print second primer */
    optStart = optimalStartPos(args, primer, args->M - args->m);

    for (i = optStart; i < args->M; i++)
      printf("%c", primer[i]);

    printf("\t%.2f\t%.1f\n", gc(primer, optStart), tm(primer, optStart));
  }
  free(chr);
  free(snpId);
  free(alt);
  free(hsName);
  free(primer);
}
