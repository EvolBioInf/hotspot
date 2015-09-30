/***** sequenceData.h **************************************
 * Description: Header file for sequence manipulation tasks.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Sun Jun  6 10:34:25 2004.
 ***********************************************************/
#ifndef SEQUENCEDATA
#define SEQUENCEDATA

#include <stdio.h>
#include <stdlib.h>

#define SEQLINE 1000 /* maximal length of one line in FASTA file; hard bound   \
                        */
#define SEQBUFFER 5000000 /* define the size of the sequence buffer */
#define DICSIZE 256
#define BORDER 'Z'
/* #define BORDER '\255' */

/* basic sequence type representing >= 1 entry in FASTA file */
typedef struct sequence {
	char *seq;          /* the sequence */
	char *id;           /* the sequence id */
	size_t numSeq;      /* number of sequences represented */
	int numQuery;       /* number of query sequences */
	size_t *borders;    /* position of last character of sequence in seq */
	char **headers;     /* FASTA header lines */
	size_t len;         /* sequence length */
	long *freqTab;      /* frequency table */
	size_t numQueryNuc; /* number of nucleotides in query sequence */
	size_t numSbjctNuc; /* number of nucleotides in sbjct sequence */
	size_t numNuc;      /* number of nucleotides in sequence */
	size_t queryStart;
	size_t queryEnd;
	size_t effQueryNuc; /* number of nucleotides in query that are the starting
	          * point for shustrings that are longer than expected
	          * by chance alone */
	int *sbjctId;       /* sequence id for each sbjct position */
} Sequence;

Sequence *revcomp(Sequence *seq);
Sequence *getNextSequence(FILE *fp);

int *getDnaDictionary(int *dic);
int *getRestrictedDnaDictionary(int *dic);

void freeSequence(Sequence *seq);
Sequence *getPermanentNextSequence(FILE *fp);
void convertToAcgt(Sequence *seq);
void resetSequenceReader();
Sequence *readFasta(int fd);
Sequence **sequence2array(Sequence *seq);
void prepareSeq(Sequence *sequence);
Sequence *catSeq(Sequence *seq1, Sequence *seq2);
/* void randomizeSeq(Sequence *seq); */
/* void randomizeSbjct(Sequence *seq); */
char *hotSpotGetChr();
Sequence *cloneSeq(Sequence *ori);
double gcContent(Sequence *seq);

#endif
