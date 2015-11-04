#ifndef INTERFACE
#define INTERFACE

#define DEFAULT_L 17 /* default oligo length */
#define DEFAULT_LU 100                /* default length of universals */
#define DEFAULT_D 200                 /* default distance between universals */

/* define argument container */
typedef struct args{
  char h;   /* help message? */
  char *g;  /* path to genome */
  char *s;  /* snp data in vcf format */
  char u;   /* universals? */
  char v;   /* version message? */
  char e;   /* error message? */
  char **inputFiles;
  int l;    /* length */
  int d;    /* distance between universals */
  int numInputFiles;
} Args;

Args *getArgs(int argc, char *argv[]);
void freeArgs(Args *);
void printUsage();
void printSplash(char *version);

#endif
