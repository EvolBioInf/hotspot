/***** sequencedata.c *********************************************
 * Description: Collection of routines for reading and
 * manipulating sequence data.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Sun Jun  6 10:34:31 2004.
 *****************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <ctype.h>
#include <glob.h>
#include "sequenceData.h"
#include "stringUtil.h"
#include "eprintf.h"
#include "tab.h"

static int lastSequence = 0;
static char *line = NULL;

/* convertToACGT: convert nucleotide data to acgt alphabet.
 */
void convertToAcgt(Sequence *seq) {
  replace(seq->seq, 'r', 'g');
  replace(seq->seq, 'y', 't');
  replace(seq->seq, 'm', 'a');
  replace(seq->seq, 'k', 'g');
  replace(seq->seq, 's', 'g');
  replace(seq->seq, 'w', 'a');
  replace(seq->seq, 'h', 'a');
  replace(seq->seq, 'b', 'g');
  replace(seq->seq, 'v', 'g');
  replace(seq->seq, 'd', 'g');
  replace(seq->seq, 'n', 'g');
  replace(seq->seq, 'u', 't');

  replace(seq->seq, 'R', 'G');
  replace(seq->seq, 'Y', 'T');
  replace(seq->seq, 'M', 'A');
  replace(seq->seq, 'K', 'G');
  replace(seq->seq, 'S', 'G');
  replace(seq->seq, 'W', 'A');
  replace(seq->seq, 'H', 'A');
  replace(seq->seq, 'B', 'G');
  replace(seq->seq, 'V', 'G');
  replace(seq->seq, 'D', 'G');
  replace(seq->seq, 'N', 'G');
  replace(seq->seq, 'U', 'T');
}

/* everything that is not [acgtACGT] is flagged by a -1 */
int *getRestrictedDnaDictionary(int *dic) {
  int i;

  if (dic == NULL)
    dic = (int *)emalloc((DICSIZE + 1) * sizeof(int));

  for (i = 0; i < DICSIZE; i++)
    dic[i] = 0;

  dic['a'] = 1; /* a */
  dic['c'] = 1; /* c */
  dic['g'] = 1; /* g */
  dic['t'] = 1; /* t */
  dic['A'] = 1; /* A */
  dic['C'] = 1; /* C */
  dic['G'] = 1; /* G */
  dic['T'] = 1; /* T */

  return dic;
}

/* getDnaDictionary: create DNA dictionary */
int *getDnaDictionary(int *dic) {
  int i;

  if (dic == NULL)
    dic = (int *)malloc((DICSIZE + 1) * sizeof(int));

  for (i = 0; i < DICSIZE; i++)
    dic[i] = 0;

  dic['a'] = 0; /* a */
  dic['c'] = 1; /* c */
  dic['g'] = 2; /* g */
  dic['t'] = 3; /* t */
  dic['A'] = 0; /* A */
  dic['C'] = 1; /* C */
  dic['G'] = 2; /* G */
  dic['T'] = 3; /* T */
  dic['r'] = dic['g'];
  dic['R'] = dic['g'];
  dic['y'] = dic['t'];
  dic['Y'] = dic['t'];
  dic['m'] = dic['a'];
  dic['M'] = dic['a'];
  dic['k'] = dic['g'];
  dic['K'] = dic['g'];
  dic['s'] = dic['g'];
  dic['S'] = dic['g'];
  dic['w'] = dic['a'];
  dic['W'] = dic['a'];
  dic['h'] = dic['a'];
  dic['H'] = dic['a'];
  dic['b'] = dic['g'];
  dic['B'] = dic['g'];
  dic['v'] = dic['g'];
  dic['V'] = dic['g'];
  dic['d'] = dic['g'];
  dic['D'] = dic['g'];
  dic['n'] = dic['g'];
  dic['N'] = dic['g'];
  dic['u'] = dic['t'];
  dic['U'] = dic['t'];

  return dic;
}

/* reverse and complement a sequence */
Sequence *revcomp(Sequence *seq) {
  long i, j, n;
  char c;
  Sequence *newSeq;
  newSeq = (Sequence *)emalloc(sizeof(Sequence));

  /*   n = strlen(seq->seq); */
  n = seq->len;
  newSeq->seq = (char *)emalloc((n + 1) * sizeof(char));
  newSeq->freqTab = NULL;
  newSeq->sbjctId = NULL;

  newSeq->id = strdup2(seq->id);
  j = 0;
  for (i = n - 1; i >= 0; i--) {
    c = seq->seq[i];
    switch (c) {
    case BORDER:
      newSeq->seq[j++] = BORDER;
      break;
    case 'A':
      newSeq->seq[j++] = 'T';
      break;
    case 'C':
      newSeq->seq[j++] = 'G';
      break;
    case 'G':
      newSeq->seq[j++] = 'C';
      break;
    case 'T':
      newSeq->seq[j++] = 'A';
      break;
    default:
      newSeq->seq[j++] = c;
      break;
    }
  }
  newSeq->seq[n] = '\0';
  return newSeq;
}
/* Get next sequence from an open data stream in FASTA format; this stream may
 * be the stdin */
Sequence *getPermanentNextSequence(FILE *fp) {
  Sequence *sequence;
  int seqlen, seqi, i, l;
  int currentBuffer;

  if (lastSequence) {
    return NULL;
  }
  if (line == NULL) {
    line = (char *)emalloc((SEQLINE + 2) * sizeof(char));
    line = fgets(line, SEQLINE, fp);
  }
  /* make a sequence object */
  sequence = (Sequence *)emalloc(sizeof(Sequence));
  /* allocate memory for sequence id */
  sequence->id = (char *)emalloc((strlen(line) + 1) * sizeof(char));
  /* copy sequence id */
  strcpy(sequence->id, chomp(line));
  /* allocate memory for sequence string */
  sequence->seq = (char *)emalloc((SEQBUFFER + 1) * sizeof(char));
  seqlen = 0;
  currentBuffer = SEQBUFFER;
  seqi = 0;
  while ((line = fgets(line, SEQLINE, fp)) != NULL) {
    if (strstr(line, ">") != NULL) {
      sequence->seq[seqi++] = '\0';
      sequence->seq = (char *)realloc(sequence->seq, seqi * sizeof(char));
      return sequence;
    }
    if (strlen(line) > SEQLINE) {
      printf("error in getNextSequence: cannot deal with lines longer "
	     "than %d bp.\n",
	     SEQLINE);
      printf("  change the SEQLINE parameter in file sequenceData.h and "
	     "recompile.\n");
      exit(2);
    }
    l = strlen(line);
    /* disregard the final carriage return */
    if (line[l - 1] == '\n')
      l--;
    seqlen += l;
    if (seqlen > currentBuffer) {
      currentBuffer += SEQBUFFER;
      sequence->seq = (char *)erealloc(sequence->seq, currentBuffer);
    }
    for (i = 0; i < l; i++) {
      sequence->seq[seqi++] = line[i];
    }
    /* sequence->seq = strncat(sequence->seq,line,strlen(line)-1); */
  }
  sequence->seq[seqi++] = '\0';
  sequence->seq = (char *)realloc(sequence->seq, seqi * sizeof(char));
  sequence->len = seqi - 1;
  lastSequence = 1;
  return sequence;
}

void resetSequenceReader() {
  line = NULL;
  lastSequence = 0;
}

/* convert multiple sequences contained in seq into an
 * array of sequences each representing a single sequence
 */
Sequence **sequence2array(Sequence *seq) {
  Sequence **seqs;
  size_t i, j, k, len;

  /* allocate space for sequences */
  seqs = (Sequence **)emalloc(seq->numSeq * sizeof(Sequence *));
  for (i = 0; i < seq->numSeq; i++) {
    seqs[i] = (Sequence *)emalloc(sizeof(Sequence));
    seqs[i]->freqTab = (long *)emalloc(DICSIZE * sizeof(long));
    seqs[i]->numNuc = 0;
    for (j = 0; j < DICSIZE; j++)
      seqs[i]->freqTab[j] = 0;
  }
  /* deal with first sequence */
  len = seq->borders[0] + 1;
  seqs[0]->seq = (char *)emalloc((len + 1) * sizeof(char));
  seqs[0]->len = len;
  for (i = 0; i < len; i++) {
    seqs[0]->seq[i] = seq->seq[i];
    seqs[0]->freqTab[(int)seq->seq[i]]++;
    seqs[0]->numNuc++;
  }
  seqs[0]->numNuc *= 2;
  seqs[0]->seq[len - 1] = BORDER;
  seqs[0]->seq[len] = '\0';
  seqs[0]->id = emalloc(6 * sizeof(char));
  seqs[0]->id[0] = '\0';
  strcat(seqs[0]->id, "strId");
  seqs[0]->numSeq = 1;
  seqs[0]->borders = (size_t *)emalloc(sizeof(size_t));
  seqs[0]->borders[0] = seq->borders[0];
  seqs[0]->headers = (char **)emalloc(sizeof(char *));
  seqs[0]->headers[0] =
    (char *)emalloc((strlen(seq->headers[0]) + 1) * sizeof(char));
  seqs[0]->headers[0] = strcpy(seqs[0]->headers[0], seq->headers[0]);
  /* deal with remaining sequences */
  for (i = 1; i < seq->numSeq; i++) {
    len = seq->borders[i] - seq->borders[i - 1];
    seqs[i]->len = len;
    seqs[i]->seq = (char *)emalloc((len + 1) * sizeof(char));
    k = 0;
    for (j = seq->borders[i - 1] + 1; j < seq->borders[i]; j++) {
      seqs[i]->seq[k++] = seq->seq[j];
      seqs[i]->freqTab[(int)seq->seq[j]]++;
      seqs[i]->numNuc++;
    }
    seqs[i]->seq[len - 1] = BORDER;
    seqs[i]->seq[len] = '\0';
    seqs[i]->id = emalloc(6 * sizeof(char));
    seqs[i]->id[0] = '\0';
    strcat(seqs[i]->id, "strId");
    seqs[i]->numSeq = 1;
    seqs[i]->borders = (size_t *)emalloc(sizeof(size_t));
    seqs[i]->borders[0] = len - 1;
    seqs[i]->headers = (char **)emalloc(sizeof(char *));
    seqs[i]->headers[0] =
      (char *)emalloc((strlen(seq->headers[i]) + 1) * sizeof(char));
    seqs[i]->headers[0][0] = '\0';
    seqs[i]->headers[0] = strcpy(seqs[i]->headers[0], seq->headers[i]);
    seqs[i]->numNuc *= 2;
  }
  return seqs;
}

/* read FASTA-formatted sequence data from an open file descriptor
 * into single sequence string
 */
Sequence *readFasta(int fd) {
  Sequence *s;
  char buf[BUFSIZ];
  int headerOpen;
  int headerLen;
  int i;
  size_t maxLen;
  int c;
  /* size_t sl; */

  if (fd < 0)
    return NULL;

  s = (Sequence *)emalloc(sizeof(Sequence));
  s->freqTab = (long *)emalloc(DICSIZE * sizeof(long));
  for (i = 0; i < DICSIZE; i++)
    s->freqTab[i] = 0;
  s->borders = (size_t *)emalloc(sizeof(size_t));
  s->headers = (char **)emalloc(sizeof(char *));
  s->id = (char *)emalloc(6 * sizeof(char));
  s->id[0] = '\0';
  strcat(s->id, "strId");
  headerOpen = 0;
  s->len = 0;
  s->numSeq = 0;
  s->sbjctId = NULL;
  maxLen = 0;
  headerLen = 0;

  while ((c = read(fd, buf, BUFSIZ)) > 0) {
    if (s->len + c + 1 > maxLen) {
      if (maxLen >= BUFSIZ) {
	maxLen *= 2;
	s->seq = (char *)erealloc(s->seq, (maxLen + 2) * sizeof(char));
      } else {
	maxLen = BUFSIZ;
	s->seq = (char *)emalloc((maxLen + 2) * sizeof(char));
      }
    }
    /* sl = maxLen; */
    /* go through the c characters just read into buf */
    for (i = 0; i < c; i++) {
      if (buf[i] == '>') {
	headerOpen = 1;
	/* take care of sequence borders */
	s->borders = (size_t *)erealloc(s->borders, (s->numSeq + 1) *
					sizeof(size_t));
	if (s->numSeq >= 1) {
	  s->seq[s->len] = BORDER; /* unique character to terminate
				      each sequence */
	  s->borders[s->numSeq - 1] = s->len;
	  s->len++;
	}
	/* take care of sequence headers */
	s->headers = (char **)erealloc(s->headers, (s->numSeq + 1) *
				       sizeof(char *));
	s->headers[s->numSeq] =
	  (char *)emalloc((BUFSIZ + 1) * sizeof(char));
	headerLen = 0;
	s->numSeq++;
      }
      if (headerOpen) {
	if (buf[i] == '\n') {
	  headerOpen = 0;
	  s->headers[s->numSeq - 1][headerLen] = '\0';
	  /* trim header to actual length */
	  s->headers[s->numSeq - 1] =
	    (char *)erealloc(s->headers[s->numSeq - 1],
			     (headerLen + 1) * sizeof(char));
	} else if (headerLen < BUFSIZ && isprint(buf[i]))
	  s->headers[s->numSeq - 1][headerLen++] = buf[i];
      } else if (buf[i] != '\n') {
	s->seq[s->len] = buf[i];
	s->freqTab[(int)buf[i]]++;
	s->len++;
      }
    }
  }
  /* add last border */
  if (s->len < maxLen)
    s->seq[s->len] = BORDER;
  else {
    printf("ERROR [readFasta]: s->len: %d; maxLen: %d\n", (int)s->len,
	   (int)maxLen);
    exit(0);
  }
  s->len++;
  /* set end of last sequence read */
  s->borders[s->numSeq - 1] = s->len - 1;
  /* trim sequence to actual size */
  s->seq = (char *)erealloc(s->seq, (s->len + 1) * sizeof(char));
  return s;
}

/* Get next sequence from an open data stream in FASTA format; this stream may
 * be the stdin */
Sequence *getNextSequence(FILE *fp) {
  Sequence *sequence;
  size_t seqi, i, l;
  size_t currentBuffer;

  if (fp == NULL || lastSequence) {
    return NULL;
  }

  if (line == NULL) {
    line = (char *)malloc((SEQLINE + 2) * sizeof(char));
    line = fgets(line, SEQLINE, fp);
  }

  /* make a sequence object */
  sequence = (Sequence *)malloc(sizeof(Sequence));
  /* allocate memory for sequence id */
  sequence->id = (char *)malloc((strlen(line) + 1) * sizeof(char));
  /* copy sequence id */
  strcpy(sequence->id, line);
  /* allocate memory for sequence string */
  sequence->seq = (char *)malloc((SEQBUFFER + 1) * sizeof(char));
  sequence->numSeq = 1;
  sequence->headers = (char **)emalloc(sizeof(char *));
  sequence->headers[0] = (char *)emalloc(sizeof(char));

  sequence->borders = (size_t *)emalloc(sizeof(size_t));

  sequence->len = 0;
  currentBuffer = SEQBUFFER;
  seqi = 0;
  while ((line = fgets(line, SEQLINE, fp)) != NULL) {
    if (line[0] == '>') {
      sequence->len++;
      sequence->seq =
	(char *)erealloc(sequence->seq, sequence->len * sizeof(char));
      sequence->seq[sequence->len - 1] = BORDER;
      sequence->borders[0] = sequence->len - 1;
      return sequence;
    }
    if (strlen(line) > SEQLINE) {
      printf("error in getNextSequence: cannot deal with lines longer "
	     "than %d bp.\n",
	     SEQLINE);
      printf("  change the SEQLINE parameter in file sequenceData.h and "
	     "recompile.\n");
      exit(2);
    }
    l = strlen(line);
    /* disregard the final carriage return */
    if (line[l - 1] == '\n')
      l--;
    sequence->len += l;
    if (sequence->len > currentBuffer) {
      currentBuffer += SEQBUFFER;
      sequence->seq = (char *)erealloc(sequence->seq, currentBuffer);
    }
    for (i = 0; i < l; i++) {
      sequence->seq[seqi++] = line[i];
    }
    /* sequence->seq = strncat(sequence->seq,line,strlen(line)-1); */
  }
  sequence->len++;
  sequence->seq =
    (char *)erealloc(sequence->seq, sequence->len * sizeof(char));
  sequence->seq[sequence->len - 1] = BORDER;
  sequence->borders[0] = sequence->len - 1;
  lastSequence = 1;
  return sequence;
}

/* freeSequence: free the data structure Sequence */
void freeSequence(Sequence *seq) {
  size_t i;

  for (i = 0; i < seq->numSeq; i++)
    free(seq->headers[i]);
  free(seq->headers);
  free(seq->borders);
  free(seq->id);
  free(seq->seq);
  free(seq->freqTab);
  free(seq->sbjctId);
  free(seq);
}

/* prepareSeq: prepares sequence string for analysis by shustring-type programs.
 * Does the following: 1) set all residues to upper case
 *                     2) generate reverse complement
 *                     3) concatenate reverse complement to end of forward
 * strand
 */
void prepareSeq(Sequence *sequence) {
  Sequence *rstrand;
  size_t i;
  char *nuc = "TCAGtcag";

  strtoupper(sequence->seq, sequence->len);
  /* take care of reverse strand */
  rstrand = revcomp(sequence);
  rstrand->headers = (char **)emalloc(sizeof(char *));
  rstrand->headers[0] = (char *)emalloc(sizeof(char));
  rstrand->borders = (size_t *)emalloc(sizeof(size_t));
  rstrand->numSeq = 1;
  sequence->seq[sequence->len] = '\0';
  sequence->len += sequence->len;
  sequence->seq =
    (char *)erealloc(sequence->seq, (sequence->len + 1) * sizeof(char));
  sequence->borders = (size_t *)erealloc(
					 sequence->borders, 2 * sequence->numSeq * sizeof(size_t));
  for (i = 1; i < sequence->numSeq; i++) {
    sequence->borders[2 * sequence->numSeq - i - 1] =
      sequence->len - sequence->borders[i - 1] - 2;
  }
  sequence->borders[2 * sequence->numSeq - 1] = sequence->len - 1;
  /* move first border of reverted sequences to the end */
  rstrand->seq++;
  strncat(sequence->seq, rstrand->seq, sequence->len);
  rstrand->seq--;
  sequence->seq[sequence->len - 1] = BORDER;
  sequence->seq[sequence->len] = '\0';
  freeSequence(rstrand);
  sequence->numNuc = 0;
  for (i = 0; i < 8; i++)
    sequence->numNuc += sequence->freqTab[(int)nuc[i]];
  sequence->numNuc *= 2;
}

/* catSeq: concatenate the sequences contained in two Sequence objects */
Sequence *catSeq(Sequence *seq1, Sequence *seq2) {
  Sequence *cat;
  size_t i, j, n;

  cat = (Sequence *)emalloc(sizeof(Sequence));
  cat->seq = (char *)emalloc((strlen(seq1->seq) + strlen(seq2->seq) + 1) *
			     sizeof(char));
  cat->seq[0] = '\0';
  cat->seq = strncat(cat->seq, seq1->seq, seq1->len);
  cat->seq = strncat(cat->seq, seq2->seq, seq2->len);
  cat->id = (char *)emalloc(6 * sizeof(char));
  cat->id[0] = '\0';
  strcat(cat->id, "strId");
  n = seq1->numSeq + seq2->numSeq;
  cat->numSeq = n;
  cat->numQuery = seq1->numSeq;
  cat->freqTab = NULL;
  cat->queryStart = 0;
  cat->queryEnd = seq1->len - 1;
  cat->borders = (size_t *)emalloc(2 * n * sizeof(size_t));
  cat->headers = (char **)emalloc(n * sizeof(char *));
  /* take care of the n headers */
  for (i = 0; i < seq1->numSeq; i++) {
    cat->headers[i] =
      (char *)emalloc((strlen(seq1->headers[i]) + 1) * sizeof(char));
    cat->headers[i] = strcpy(cat->headers[i], seq1->headers[i]);
  }
  j = i;
  for (i = 0; i < seq2->numSeq; i++) {
    cat->headers[j] =
      (char *)emalloc((strlen(seq2->headers[i]) + 1) * sizeof(char));
    cat->headers[j] = strcpy(cat->headers[j], seq2->headers[i]);
    j++;
  }
  /* take care of the 2n borders */
  for (i = 0; i < 2 * seq1->numSeq; i++) {
    cat->borders[i] = seq1->borders[i];
  }
  j = i;
  for (i = 0; i < 2 * seq2->numSeq; i++) {
    cat->borders[j + i] =
      seq1->borders[2 * seq1->numSeq - 1] + seq2->borders[i] + 1;
  }
  /* sbjct IDs */
  cat->sbjctId = (int *)emalloc(seq2->len * sizeof(int));
  for (i = 0; i < seq2->len; i++)
    cat->sbjctId[i] = -1;
  for (i = 0; i <= seq2->borders[0]; i++)
    cat->sbjctId[i] = 0;
  for (i = 1; i < seq2->numSeq; i++)
    for (j = seq2->borders[i - 1] + 1; j <= seq2->borders[i]; j++)
      cat->sbjctId[j] = i;
  n = seq2->numSeq - 1;
  for (; i < 2 * seq2->numSeq; i++) {
    for (j = seq2->borders[i - 1] + 1; j <= seq2->borders[i]; j++)
      cat->sbjctId[j] = n;
    n--;
  }

  cat->len = seq1->len + seq2->len;
  cat->numQueryNuc = seq1->numNuc;
  cat->numSbjctNuc = seq2->numNuc;
  cat->numNuc = seq1->numNuc + seq2->numNuc;
  return cat;
}

/* /\* randomizeSeq: in-place randomization of sequence string of query and
 * subject sequence(s) of seq  */
/*  * deals with: */
/*  * - multiple sequences */
/*  * - forward & reverse strands */
/*  *\/ */
/* void randomizeSeq(Sequence *seq){ */
/*   int  i, j, r1, r2, hi, lo; */
/*   int startInd, endInd, x; */
/*   char tmp; */

/*   /\* shuffle query *\/ */
/*   lo = 0; */
/*   for(i=0;i<seq->numQuery;i++){ */
/*     if(i) */
/*       lo = seq->borders[i-1] + 1; */
/*     hi = seq->borders[i] - 1; */
/*     for(j=lo; j < seq->borders[i]; j++){ */
/*       /\* forward strand *\/ */
/*       r1 = getRandMinMaxInt(j, hi); */
/*       tmp = seq->seq[r1]; */
/*       seq->seq[r1] = seq->seq[j]; */
/*       seq->seq[j] = tmp; */
/*       /\* reverse strand *\/ */
/*       r1 = seq->queryEnd - r1 - 1; */
/*       r2 = seq->queryEnd - j - 1; */
/*       tmp = seq->seq[r1]; */
/*       seq->seq[r1] = seq->seq[r2]; */
/*       seq->seq[r2] = tmp; */
/*     } */
/*   } */
/*   /\* shuffle sbjct *\/ */
/*   lo = seq->queryEnd + 1; */
/*   startInd = 2*seq->numQuery; */
/*   x = seq->queryEnd - 1; */
/*   endInd = startInd + seq->numSeq - seq->numQuery; */
/*   for(i=startInd;i<endInd;i++){ */
/*     if(i>startInd) */
/*       lo = seq->borders[i-1] + 1; */
/*     hi = seq->borders[i] - 1; */
/*     for(j=lo; j < seq->borders[i]; j++){ */
/*       /\* forward strand *\/ */
/*       r1 = getRandMinMaxInt(j, hi); */
/*       tmp = seq->seq[r1]; */
/*       seq->seq[r1] = seq->seq[j]; */
/*       seq->seq[j] = tmp; */
/*       /\* reverse strand *\/ */
/*       r1 = seq->len - r1 + x; */
/*       r2 = seq->len - j + x; */
/*       tmp = seq->seq[r1]; */
/*       seq->seq[r1] = seq->seq[r2]; */
/*       seq->seq[r2] = tmp; */
/*     } */
/*   } */
/* } */

/* /\* randomizeSbjct: in-place randomization of sequence string of subject
 * sequence(s) of seq  */
/*  * deals with: */
/*  * - multiple sequences */
/*  * - forward & reverse strands */
/*  *\/ */
/* void randomizeSbjct(Sequence *seq){ */
/*   int  i, j, r1, r2, hi, lo; */
/*   int startInd, endInd, x; */
/*   char tmp; */

/*   /\* shuffle sbjct *\/ */
/*   lo = seq->queryEnd + 1; */
/*   startInd = 2*seq->numQuery; */
/*   x = seq->queryEnd - 1; */
/*   endInd = startInd + seq->numSeq - seq->numQuery; */
/*   for(i=startInd;i<endInd;i++){ */
/*     if(i>startInd) */
/*       lo = seq->borders[i-1] + 1; */
/*     hi = seq->borders[i] - 1; */
/*     for(j=lo; j < seq->borders[i]; j++){ */
/*       /\* forward strand *\/ */
/*       r1 = getRandMinMaxInt(j, hi); */
/*       tmp = seq->seq[r1]; */
/*       seq->seq[r1] = seq->seq[j]; */
/*       seq->seq[j] = tmp; */
/*       /\* reverse strand *\/ */
/*       r1 = seq->len - r1 + x; */
/*       r2 = seq->len - j + x; */
/*       tmp = seq->seq[r1]; */
/*       seq->seq[r1] = seq->seq[r2]; */
/*       seq->seq[r2] = tmp; */
/*     } */
/*   } */
/* } */

/* cloneSeq: make exact copy of Sequence object */
Sequence *cloneSeq(Sequence *ori) {
  Sequence *clone;
  size_t i;

  clone = (Sequence *)emalloc(sizeof(Sequence));
  clone->seq = (char *)emalloc((ori->len + 1) * sizeof(char));
  clone->seq = strncpy(clone->seq, ori->seq, ori->len);
  clone->id = (char *)emalloc(6 * sizeof(char));
  clone->id = strncpy(clone->id, ori->id, 6);
  clone->numSeq = ori->numSeq;
  clone->numQuery = ori->numQuery;
  clone->borders = (size_t *)emalloc(ori->numSeq * sizeof(size_t));
  for (i = 0; i < ori->numSeq; i++)
    clone->borders[i] = ori->borders[i];
  clone->headers = (char **)emalloc(ori->numSeq * sizeof(char *));
  for (i = 0; i < ori->numSeq; i++) {
    clone->headers[i] =
      (char *)emalloc((strlen(ori->headers[i]) + 1) * sizeof(char));
    clone->headers[i] = strcpy(clone->headers[i], ori->headers[i]);
  }
  clone->len = ori->len;
  clone->freqTab = (long *)emalloc(DICSIZE * sizeof(long));
  for (i = 0; i < DICSIZE; i++)
    clone->freqTab[i] = ori->freqTab[i];
  clone->numQueryNuc = ori->numQueryNuc;
  clone->numSbjctNuc = ori->numSbjctNuc;
  clone->numNuc = ori->numNuc;
  clone->queryStart = ori->queryStart;
  clone->queryEnd = ori->queryEnd;

  return clone;
}

double gcContent(Sequence *seq) {
  size_t i, j, numChar;
  size_t min, max, gc;

  min = 0;
  gc = 0;
  numChar = 0;
  for (i = 0; i < seq->numSeq; i++) {
    if(seq->numSeq > 1)
      max = seq->borders[i];
    else
      max = seq->len;
    if (i)
      min = seq->borders[i - 1] + 1;
    for (j = min; j < max; j++) {
      if (seq->seq[j] == 'A' || seq->seq[j] == 'C' ||
	  seq->seq[j] == 'G' || seq->seq[j] == 'T') {
	numChar++;
	if (seq->seq[j] == 'G' || seq->seq[j] == 'C')
	  gc++;
      }
    }
  }
  return (double)gc / (double)(numChar);
}
