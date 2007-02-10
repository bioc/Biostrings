/****************************************************************************
                      A NAIVE METHOD FOR EXACT MATCHING
		             Author: Herve Pages

 Here is how the "naive" (aka "memcmp", aka "blunt") method for finding exact
 matches is described in Dan Gusfield book "Algorithms on strings, trees, and
 sequences" (slightly modified):
     The naive method aligns the left end of P (the pattern) with the left
     end of S (the subject) and then compares the characters of P and S left
     to right until either two unequal characters are found or until P is
     exhausted, in which case an occurence of P is reported. In either case,
     P is then shifted one place to the right, and the comparisons are
     restarted from the left end of P. This process repeats until the right
     end of P shifts past the right end of S.
 
 Why do we need a "naive" algo?
   - For QC: we can validate other more sophisticated matching algo by
     comparing their results to those obtains with the "naive" algo.
   - To use as a point of reference when comparing performance.
   
 ****************************************************************************/
#include "Biostrings.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/****************************************************************************/
static int debug = 0;

SEXP match_naive_debug()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in 'match_naive.c'\n", debug ? "on" : "off");
#else
	Rprintf("Debug mode not available in 'match_naive.c'\n");
#endif
	return R_NilValue;
}


/****************************************************************************
 * A memcmp-based implementation of the "naive" method
 * ===================================================
 */

/* Returns the number of matches */
static int naive_search(char *P, int nP, char *S, int nS,
		int is_count_only, SEXP *p_matchpos, PROTECT_INDEX matchpos_pi)
{
	int count = 0, *index, n1, n2;

	if (!is_count_only) {
		index = INTEGER(*p_matchpos);
	}
	n1 = 0;
	n2 = n1 + nP;
	for (n1 = 0, n2 = nP; n2 <= nS; n1++, n2++) {
		if (memcmp(P, S + n1, nP) != 0)
			continue;
		if (!is_count_only) {
			if (count >= LENGTH(*p_matchpos)) {
				*p_matchpos = Biostrings_expandMatchIndex(
						*p_matchpos, n2, nS - n2);
				REPROTECT(*p_matchpos, matchpos_pi);
				index = INTEGER(*p_matchpos);
			}
			index[count] = n1;
		}
		count++;
	}
	return count;
}

SEXP match_naive(SEXP p_xp, SEXP p_offset, SEXP p_length,
		SEXP s_xp, SEXP s_offset, SEXP s_length,
		SEXP count_only)
{
	int pat_offset, pat_length, subj_offset, subj_length,
	    is_count_only, count;
	char *pat, *subj;
	SEXP ans;
	int matchpos_length;
	SEXP matchpos = R_NilValue;
	PROTECT_INDEX matchpos_pi;

	pat_offset = INTEGER(p_offset)[0];
	pat_length = INTEGER(p_length)[0];
	pat = CHAR(R_ExternalPtrTag(p_xp)) + pat_offset;
	subj_offset = INTEGER(s_offset)[0];
	subj_length = INTEGER(s_length)[0];
	subj = CHAR(R_ExternalPtrTag(s_xp)) + subj_offset;
	is_count_only = LOGICAL(count_only)[0];

	if (!is_count_only) {
		matchpos_length = Biostrings_estimateExpectedMatchCount(
					pat_length, subj_length, 4);
		PROTECT_WITH_INDEX(matchpos, &matchpos_pi);
		matchpos = allocVector(INTSXP, matchpos_length);
		REPROTECT(matchpos, matchpos_pi);
	}
	count = naive_search(pat, pat_length, subj, subj_length,
				is_count_only, &matchpos, matchpos_pi);
	if (is_count_only) {
		PROTECT(ans = allocVector(INTSXP, 1));
		INTEGER(ans)[0] = count;
		UNPROTECT(1);
	} else {
		PROTECT(ans = allocVector(INTSXP, count));
		memcpy(INTEGER(ans), INTEGER(matchpos), sizeof(int) * count);
		UNPROTECT(2);
	}
	return ans;
}

