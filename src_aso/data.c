/***** data.c *************************************
 * Description:
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Jun  8 12:19:13 2015
 **************************************************/
#include <stdio.h>
#include <glob.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include "eprintf.h"
#include "data.h"
#include "tab.h"

int hotSpotStart = -1;
int hotSpotEnd = -1;
char *hotSpotChr = NULL;

char *hotSpotGetChr() {
  char *chr;

  chr = tabField(CHR_FIELD);
  chr += 3;
  if (chr[0] == '0')
    chr++;
  return chr;
}

Sequence *hotSpotSeq(Args *args) {
  char *chr, *name;
  int fd, status;
  Sequence *seq;
  glob_t *pglob;

  pglob = NULL;
  pglob = (glob_t *)emalloc(sizeof(glob_t));
  name = (char *)emalloc(256 * sizeof(char));
  name[0] = '\0';
  chr = hotSpotGetChr();

  /* construct name of genome file */
  strcat(name, args->g);
  strcat(name, "/*_chr");
  strcat(name, chr);
  strcat(name, ".fa*");
  /* open genome file and read genome */
  status = glob(name, 0, NULL, pglob);
  if (status == GLOB_NOMATCH) {
    fprintf(stderr, "ERROR: Could not open genome file %s\n", name);
    exit(-1);
  }
  fd = open(pglob->gl_pathv[0], O_RDONLY, 0);
  if (fd < 0) {
    fprintf(stderr, "ERROR: Could not open genome file %s\n", name);
    exit(-1);
  }
  seq = readFasta(fd);
  close(fd);
  free(name);
  globfree(pglob);
  free(pglob);
  return seq;
}

FILE *hotSpotSnp(Args *args) {
  char *name, *chr, *cmd;
  int start, end;
  FILE *sp;

  name = (char *)emalloc(256 * sizeof(char));
  cmd = (char *)emalloc(256 * sizeof(char));
  /* construct name of snp file */
  name[0] = '\0';
  chr = hotSpotGetChr();
  hotSpotChr = chr;
  start = atoi(tabField(START_FIELD));
  end = atoi(tabField(END_FIELD));
  hotSpotStart = start;
  hotSpotEnd = end;
  strcat(name, args->s);
  strcat(name, "/vcf_chr_");
  strcat(name, chr);
  strcat(name, ".vcf.gz");
  if (args->u)
    sprintf(cmd, "tabix %s %s:%d-%d", name, chr, start, end);
  else
    sprintf(cmd, "tabix %s %s:%d-%d | grep snp", name, chr, start, end);
  /* open pipeline for retrieving SNPs */
  sp = popen(cmd, "r");
  if (sp == NULL) {
    fprintf(stderr, "ERROR: Could not retrieve SNPs via %s\n", cmd);
    exit(-1);
  }
  free(name);
  free(cmd);
  return sp;
}

int getHotSpotStart(){
  return hotSpotStart;
}

int getHotSpotEnd(){
  return hotSpotEnd;
}

char *getHotSpotChr(){
  return hotSpotChr;
}
