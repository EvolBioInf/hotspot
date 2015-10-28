/***** data.c *************************************
 * Description:
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Jun  8 12:29:16 2015
 **************************************************/
#define _GNU_SOURCE
#include <glob.h>
#include <string.h>
#include <stdio.h>
#include <err.h>
#include <fcntl.h>
#include <unistd.h>
#include "eprintf.h"
#include "data.h"
#include "tab.h"

char *hotSpotChr = NULL;
int hotSpotStart = -1;
int hotSpotEnd = -1;

char *hotSpotGetChr() {
  char *chr;

  chr = tabField(CHR_FIELD);
  chr += 3;
  if (chr[0] == '0')
    chr++;
  return chr;
}

Sequence *hotSpotSeq(Args *args) {
  char *chr, *name = NULL;
  int fd, status;
  Sequence *seq;
  glob_t *pglob;

  pglob = (glob_t *)emalloc(sizeof(glob_t));
  chr = hotSpotGetChr();
  hotSpotChr = chr;
  /* construct name of genome file */
  int check = asprintf( &name, "%s/*_chr%s.fa", args->g, chr);
  if(check == -1){
    err(1, "Out of mem?");
  }
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

FILE *hotSpotSnp(Args *args, int start, int end) {
  char *name = NULL, *cmd = NULL, *chr;
  FILE *sp;

  chr = hotSpotGetChr();

  /* construct name of snp file */
  int check = asprintf( &name, "%s/vcf_chr_%s.vcf.gz", args->s, chr);
  if(check == -1){
    err(1, "Out of mem?");
  }

  check = asprintf( &cmd, "tabix %s %s:%d-%d | grep snp", name, chr, start, end);
  if(check == -1){
    err(1, "Out of mem?");
  }
  if(args->b)
    printf("%s\n",cmd);
  /* open pipeline for retrieving SNPs */
  sp = popen(cmd, "r");
  if (sp == NULL) {
    err( -1, "ERROR: Could not retrieve SNPs via %s\n", cmd);
  }

  free(name);
  free(cmd);
  return sp;
}
