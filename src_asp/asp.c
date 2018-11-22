/***** asp.c **************************************
 * Description:
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Jun  7 21:51:05 2015
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "interface.h"
#include "data.h"
#include "tab.h"
#include "eprintf.h"
#include "primers.h"
#include "config.h"

void scanFile(FILE *fp, Args *args);

int main(int argc, char *argv[]) {
  int i;
  char *version, *cmd;
  Args *args;
  FILE *fp;

  cmd = (char *)emalloc(256 * sizeof(char));
  version = VERSION;
  setprogname2("asp");
  args = getArgs(argc, argv);
  if (args->v)
    printSplash(version);
  if (args->h || args->e)
    printUsage(version);
  if (args->numInputFiles == 0) {
    fp = stdin;
    scanFile(fp, args);
  } else {
    for (i = 0; i < args->numInputFiles; i++) {
      sprintf(cmd, "sort -k 2 %s", args->inputFiles[i]);
      fp = popen(cmd, "r");
      scanFile(fp, args);
      pclose(fp);
    }
  }
  tabReset();
  freeArgs(args);
  free(progname());
  free(cmd);
  return 0;
}

void scanFile(FILE *fp, Args *args) {
  Sequence *seq = NULL;
  char *oldChr;
  int start, end;
  FILE *snp;

  oldChr = (char *)emalloc(sizeof(char));
  oldChr[0] = '\0';
  while (tabGetLine(fp) != NULL) {
    /* read chromosome, start, and end of hotspot from line */
    if (tabField(0)[0] == '#')
      continue;
    if (strcmp(oldChr, hotSpotGetChr()) != 0) {
      if (seq)
	freeSequence(seq);
      seq = hotSpotSeq(args);
      free(oldChr);
      oldChr = estrdup(hotSpotGetChr());
    }
    start = atoi(tabField(START_FIELD));
    end = atoi(tabField(END_FIELD));
    if (args->r) {
      /* reverse primers are downstream of hotspot */
      start = end;
      end += args->f;
      snp = hotSpotSnp(args, start, end);
      if(args->u)
	reversePuniv(args, seq, snp, start, end);
      else
	reverseP(args, seq, snp);
    } else {
      /* forward primers are upstream of hotspot */
      end = start;
      start -= args->f;
      snp = hotSpotSnp(args, start, end);
      if(args->u)
	forwardPuniv(args, seq, snp, start, end);
      else
	forwardP(args, seq, snp);
    }
    pclose(snp);
  }
  freeSequence(seq);
  free(oldChr);
}
