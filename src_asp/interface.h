#ifndef INTERFACE
#define INTERFACE

#define DEFAULT_MI 14   /* default minimul length */
#define DEFAULT_MA 19   /* default maximal length */
#define DEFAULT_O 0.5   /* optimal GC content */
#define DEFAULT_F 5000  /* length of flanking region */
#define DEFAULT_DF 600  /* distance between forward universals */
#define DEFAULT_DR 3000 /* distance between reverse universals */

/* define argument container */
typedef struct args {
  char h;  /* help message? */
  char *g; /* path to genome */
  char *s; /* snp data in vcf format */
  char r;  /* non-SNP oligos */
  char u;  /* universals? */
  char v;  /* version message? */
  char e;  /* error message? */
  char **inputFiles;
  int m;   /* minimum length */
  int M;   /* maximum length */
  int d;   /* distance between universals, forward */
  int D;   /* distance between universals, reverse */
  int f;   /* length of flanking region */
  float o; /* optimal gc content */
  int numInputFiles;
} Args;

Args *getArgs(int argc, char *argv[]);
void printUsage();
void printSplash(char *version);
void freeArgs(Args *args);

#endif
