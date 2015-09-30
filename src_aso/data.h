/***** data.h *************************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Jun  8 12:19:17 2015
 **************************************************/
#ifndef DATA
#define DATA

#include <stdio.h>
#include "interface.h"
#include "sequenceData.h"

#define INT_FIELD 0
#define CHR_FIELD 1
#define START_FIELD 2
#define END_FIELD 3

Sequence *hotSpotSeq(Args *args);
FILE *hotSpotSnp(Args *args);
char *hotSpotGetChr();
int getHotSpotStart();
int getHotSpotEnd();
char *getHotSpotChr();

#endif
