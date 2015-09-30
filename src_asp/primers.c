/***** primers.c **********************************
 * Description:
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Jun  7 11:44:03 2015
 **************************************************/
#include <string.h>
#include <math.h>
#include "eprintf.h"
#include "primers.h"
#include "tab.h"
#include "stringUtil.h"
#include "data.h"

float gc(const char *primer, size_t startPos) {
  primer += startPos;
  size_t length = 0, gc = 0;
  char c;

  while( (c= *primer++) != '\0'){
    length++;
    if (c == 'c' || c == 'C' || c == 'g' ||c == 'G'){
      gc++;
    }
  }
	
  return gc / (float)length;
}

float absolute(float x) {
  return fabsf(x);
}

int optimalStartPos(Args *args, char *primer, int minStart) {
  int i, l, c, mi;
  float gc, ad, minDiff;

  gc = 0.;
  l = strlen(primer);
  c = 0;
  for (i = minStart; i < l; i++) {
    c++;
    if (primer[i] == 'c' || primer[i] == 'C' || primer[i] == 'g' ||
	primer[i] == 'G')
      gc++;
  }
  minDiff = absolute(args->o - gc / c);
  mi = minStart;
  c = 0;
  for (i = minStart - 1; i >= 0; i--) {
    c++;
    if (primer[i] == 'c' || primer[i] == 'C' || primer[i] == 'g' ||
	primer[i] == 'G')
      gc++;
    ad = absolute(args->o - gc / c);
    if (ad < minDiff) {
      minDiff = ad;
      mi = i;
    }
  }
  return mi;
}

void complement(char *seq) {
  char *p = seq;
  size_t i, n = strlen(seq);
  for (i = 0; i < n; i++, p++) {
    switch(*p){
    case 'a':
    case 'A':
      *p = 'T'; break;
    case 'c':
    case 'C':
      *p = 'G'; break;
    case 'g':
    case 'G':
      *p = 'C'; break;
    case 't':
    case 'T':
      *p = 'A'; break;
    default:
      *p = '!'; break;
    }
  }
}

/* tm: melting temperature
 * source: "Basic Melting Temperature Calculations" on
 * http://www.basic.northwestern.edu/biotools/oligocalc.html
 */
double tm(char *oligo, size_t startPos){
  int i, len, na, nc, ng, nt;
  double tm;

  oligo += startPos;
  len = strlen(oligo);
  na = nc = ng = nt = 0;
  for(i=0;i<len;i++){
    if(oligo[i] == 'A')
      na++;
    else if(oligo[i] == 'C')
      nc++;
    else if(oligo[i] == 'G')
      ng++;
    else if(oligo[i] == 'T')
      nt++;
  }
  if(len < 14)
    tm = (na+nt) * 2. + (ng+nc) * 4.;
  else
    tm = 64.9+41.*(ng+nc-16.4)/len;
  return tm;
}

/* intersects: given an array pos of n positions,
 *   determine by binary search whether interval 
 *   (s,e) intersects any of them; if so, return 
 *   the point of of interesection; otherwise, 
 *   return 0 (SNP positions start at 1)
 */
int intersects(int *pos, int n, int s, int e){
  int l, u, i;

  l = 0;
  u = n - 1;
  
  while(u > l){
    i = (l + u) / 2;
    if(s <= pos[i] && e >= pos[i])
      return pos[i];
    else if(s < pos[i])
      u = i - 1;
    else if(s > pos[i])
      l = i + 1;
    else
      return pos[i];
  }
  return 0;
}
