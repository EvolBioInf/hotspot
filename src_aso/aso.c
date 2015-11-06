/***** hotOligos.c ********************************
 * Description:
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Jun  7 16:35:35 2015
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tab.h"
#include "interface.h"
#include "eprintf.h"
#include "oligos.h"
#include "data.h"
#include "hotspot.h"

void scanFile(FILE *fp, Args *args) {
  Sequence *seq = NULL;
  char *oldChr;
  FILE *snp;
  int c;

  oldChr = (char *)emalloc(sizeof(char));
  oldChr[0] = '\0';
  while (tabGetLine(fp) != NULL) {
    /* read chromosome, start, and end of hotspot from line */
    if (tabField(0)[0] == '#')
      continue;
    if (strcmp(oldChr, hotSpotGetChr()) != 0) {
      if (seq)
	free(seq);
      seq = hotSpotSeq(args);
      free(oldChr);
      oldChr = estrdup(hotSpotGetChr());
    }
    c++;
    snp = hotSpotSnp(args);
    if (args->u)
      universals(args, seq, snp);
    else
      snpOligos(args, seq, snp);
    pclose(snp);
  }
  freeSequence(seq);
  free(oldChr);
}

int main(int argc, char *argv[]) {
  int i;
  char *version, *cmd;
  Args *args;
  FILE *fp;

  cmd = (char *)emalloc(256 * sizeof(char));
  version = VERSION;
  setprogname2("aso");
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
      sprintf(cmd, "sort %s", args->inputFiles[i]);
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

