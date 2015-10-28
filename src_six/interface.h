#ifndef INTERFACE
#define INTERFACE

#define DEFAULT_M 3 
#define DEFAULT_R 5
#define DEFAULT_MM0 60
#define DEFAULT_MM1 120
#define DEFAULT_MM2 240
#define DEFAULT_N 12
#define DEFAULT_X 1.0

/* define argument container */
typedef struct args{
  char h;   /* help message? */
  char v;   /* version message? */
  char e;   /* error message? */
  char **inputFiles;
  int s;    /* seed for random number generator */
  int m;    /* number of experiments */
  int n;    /* number of replicates per experiment */
  int r;    /* number of replicates */
  int *mol; /* number of molecules per experiment */
  float x;  /* rate of crossover */
  int numInputFiles;
} Args;

Args *getArgs(int argc, char *argv[]);
void printUsage();
void printSplash(char *version);

#endif
