/***** univPrimers.c ******************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Sep 21 11:27:55 2015
 **************************************************/
#include <stdio.h>
#include <string.h>
#include "tab.h"
#include "data.h"
#include "eprintf.h"
#include "primers.h"

void forwardPuniv(Args *args, Sequence *seq, FILE *snp, int hsStart, int hsEnd){
  char *hsName, *chr, *primer;
  int primerStart, pos, *snpPos, i, optStart;
  int p, s, e, numSnp;

  hsName = (char *)emalloc(256);
  primer = (char *)emalloc(args->M + 1);
  chr = (char *)emalloc(256);
  /* get all SNPs */
  numSnp = 0;
  snpPos = NULL;
  while(fscanf(snp,"%s %d %*s %*s %*s %*s %*s %*s",chr,&pos) != EOF){
    snpPos = (int *)erealloc(snpPos,sizeof(int)*(numSnp+1));
    snpPos[numSnp++] = pos;
  }
  /* get name of hotspot */
  strcpy(hsName, tabField(INT_FIELD));
  primerStart = hsStart;
  while(primerStart + args->d + args->M < hsEnd){
    /* get primer */
    for(i=0;i<args->M;i++)
      primer[i] = seq->seq[primerStart+i];
    primer[i] = '\0';
    optStart = optimalStartPos(args, primer, args->M - args->m);
    /* check whether primer intersects SNP */
    s = primerStart + optStart;
    e = primerStart + optStart + args->M - 1;
    while((p = intersects(snpPos, numSnp, s, e))){
      primerStart = p + 1;
      for(i=0;i<args->M;i++)
	primer[i] = seq->seq[primerStart+i];
      primer[i] = '\0';
      optStart = optimalStartPos(args, primer, args->M - args->m);
      s = primerStart + optStart;
      e = primerStart + optStart + args->M - 1;
    }
    /* print primer */
    printf("%s\t%s\t%d\t",hsName,chr,primerStart+optStart+1);
    for(i=optStart; i<args->M; i++)
      printf("%c",primer[i]);
    printf("\t%.2f",gc(primer, optStart));
    printf("\t%.2f\n",tm(primer, optStart));
    /* advance primer start */
    primerStart += (args->d + args->M - 1);
  }
  free(hsName);
  free(primer);
  free(chr);
}

void reversePuniv(Args *args, Sequence *seq, FILE *snp, int hsStart, int hsEnd){
  char *hsName, *chr, *primer;
  int primerStart, pos, i, optStart;
  int *snpPos, numSnp, s, e, p;

  hsName = (char *)emalloc(256);
  primer = (char *)emalloc(args->M + 1);
  chr = (char *)emalloc(256);
  /* get name of hotspot */
  strcpy(hsName, tabField(INT_FIELD));

  /* get all SNPs */
  numSnp = 0;
  snpPos = NULL;
  while(fscanf(snp,"%s %d %*s %*s %*s %*s %*s %*s",chr,&pos) != EOF){
    snpPos = (int *)erealloc(snpPos,sizeof(int)*(numSnp+1));
    snpPos[numSnp++] = pos;
  }
  /* walk from right to left while avoiding SNPs */
  primerStart = hsEnd;
  while(primerStart > hsStart + args->D + args->M){
    /* get primer */
    for(i=0;i<args->M;i++)
      primer[i] = seq->seq[primerStart+i];
    primer[i] = '\0';
    reverse(primer);
    complement(primer);
    optStart = optimalStartPos(args, primer, args->M - args->m);
    s = primerStart + optStart;
    e = primerStart + optStart + args->M - 1;
    /* check whether primer intersects SNP */
    while((p = intersects(snpPos, numSnp, s, e)) && p - hsStart > args->D){
      primerStart = p + 1;
      for(i=0;i<args->M;i++)
	primer[i] = seq->seq[primerStart+i];
      primer[i] = '\0';
      reverse(primer);
      complement(primer);
      optStart = optimalStartPos(args, primer, args->M - args->m);
      s = primerStart + optStart;
      e = primerStart + optStart + args->M - 1;
    }
    printf("%s\t%s\t%d\t",hsName,chr,primerStart+1);
    for(i=optStart; i<args->M; i++)
      printf("%c",primer[i]);
    printf("\t%.2f",gc(primer, optStart));
    printf("\t%.2f\n",tm(primer, optStart));
    primerStart -= (args->D + args->M);
  }
  free(hsName);
  free(primer);
  free(chr);
}
