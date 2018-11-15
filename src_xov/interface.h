#ifndef INTERFACE
#define INTERFACE

#define DEFAULT_C 0.95
#define DEFAULT_M 0.
#define DEFAULT_MM 1.
#define DEFAULT_I 100
#define DEFAULT_D "\t"

/* define argument container */
typedef struct args{
  char h;   /* help message? */
  char v;   /* version message? */
  char e;   /* error message? */
  char p;   /* print log-likelihood curve */
  char d;   /* print data */
  char *D;  /* delimiter */
  char o;   /* use Poisson approximation? */
  float m;  /* minimum x-value for log-likelihood curve */
  float M;  /* maximum x-value for log-likelihood curve */
  int i;    /* number of intervals for log-likelihood curve */
  float c;  /* confidence interval */
  char **inputFiles;
  int numInputFiles;
} Args;

Args *getArgs(int argc, char *argv[]);
void printUsage();
void printSplash(char *version);

#endif
