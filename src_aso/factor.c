/***** factor.c ***********************************
 * Description: Compute the longest previous factor
 *   array using a suffix array and a longest
 *   common prefix array.
 * Reference: Crochemore, M., Ilie, L. and Smyth,
 *   W. F. (2008). A simple algorithm for com-
 *   puting the Lempel Ziv factorization. In:
 *   Data Compression Conference, p. 482-488.
 *   Computing longest previous factor in linear
 *   time and applications.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Jul 15 10:29:09 2013
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "factor.h"
#include "stack.h"
#include "eprintf.h"

long minimum(long a, long b) {
	if (a < b)
		return a;
	else
		return b;
}

long maximum(long a, long b) {
	if (a > b)
		return a;
	else
		return b;
}

/*
 * computeLpf: Compute longest previous factor
 * Reference: M. Crochemore, L. Ilie, W.F. Smyth.
 *   A simple algorithm for computing the Lempel-Ziv
 *   factorization, in: J.A. Storer, M.W. Marcellini
 *   (Eds.), 18th Data Compression Conference, IEEE
 *   Computer Society, Los Alamitos, CA, 2008,
 *   pp. 482-488.
 */
long *computeLpf(long *sa, long *lcp, long n) {
	long i;
	long *lpf;

	lpf = (long *)emalloc(n * sizeof(long));
	sa[n] = -1;
	lcp[n] = 0;

	stackInit(1);
	stackPush(0);

	for (i = 1; i <= n; i++) {
		while (!stackEmpty() &&
		       (sa[i] < sa[stackTop()] ||
		        (sa[i] > sa[stackTop()] && lcp[i] <= lcp[stackTop()]))) {
			if (sa[i] < sa[stackTop()]) {
				lpf[sa[stackTop()]] = maximum(lcp[stackTop()], lcp[i]);
				lcp[i] = minimum(lcp[stackTop()], lcp[i]);
			} else
				lpf[sa[stackTop()]] = lcp[stackTop()];
			stackPop();
		}
		if (i < n)
			stackPush(i);
	}
	freeStack();
	return lpf;
}

long computeLempelZivFact(long *lpf, long n) {
	long *lz, i;

	lz = (long *)emalloc(n * sizeof(long));
	lz[0] = 0;
	i = 0;
	while (lz[i] < n) {
		lz[i + 1] = lz[i] + maximum(1, lpf[lz[i]] + 1);
		i++;
	}

	free(lz);
	free(lpf);

	return i;
}

long numLzFact(long *sa, long *lcp, long n) {
	long *lpf;

	lpf = computeLpf(sa, lcp, n);
	return computeLempelZivFact(lpf, n);
}
