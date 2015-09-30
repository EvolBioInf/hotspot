/***** ml.h ***************************************
 * Description: Header for maximum likelihood 
 *   computation.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Fri Sep 18 15:00:01 2015
 **************************************************/
#ifndef ML
#define ML

#include "interface.h"

typedef struct data{
  int n;
  int maxN;
  int *numPos;
  int *numNeg;
  int *numMol;
  double mlXover;
  double ci;
} Data;

Data *newData(int maxN);
void expandData(Data *d);
double crossoverFreq(Data *data);
double upperCL(Data *data, double mlEst);
double lowerCL(Data *data, double mlEst);
void printLikFun(Data *data, Args *args);

#endif
