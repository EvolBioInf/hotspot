/***** ran.h *****************************************************
 * Description: Header file for the random number generator
 *              contained in the file ran.c. See that
 *              file for copyright information, etc.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Fri Nov 19 12:12:56 2004.
 *****************************************************************/
#ifndef RAN
#define RAN

#include "interface.h"

#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

/* uniform random number generation */
void init_genrand(unsigned long s);
unsigned long genrand_int32(void);
long genrand_int31(void);
double genrand_real1(void);
double genrand_real2(void);
double genrand_real3(void);
double genrand_res53(void);

/* normal distribution */
double rand_Normal (void);

/* initialization/closing of random number generator */
int getSeed(Args *args);
void saveSeed(Args *args);

/* miscellaneous functions */
unsigned long getRandMinMaxInt(long min, long max);

#endif
