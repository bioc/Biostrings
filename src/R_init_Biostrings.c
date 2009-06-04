#include "Biostrings.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

#define REGISTER_CCALLABLE(fun) \
	R_RegisterCCallable("Biostrings", #fun, (DL_FUNC) &fun)

static const R_CallMethodDef callMethods[] = {

/* utils.c */
	CALLMETHOD_DEF(debug_utils, 0),

/* copy_seq.c */
	CALLMETHOD_DEF(debug_copy_seq, 0),

/* RoSeq_utils.c */
	CALLMETHOD_DEF(debug_RoSeq_utils, 0),
	CALLMETHOD_DEF(new_RawPtr_from_STRSXP, 5),

/* XString_class.c */
	CALLMETHOD_DEF(debug_XString_class, 0),
	CALLMETHOD_DEF(init_DNAlkups, 2),
	CALLMETHOD_DEF(init_RNAlkups, 2),
	CALLMETHOD_DEF(new_RawPtr_from_XString, 4),

/* XStringSet_class.c */
	CALLMETHOD_DEF(debug_XStringSet_class, 0),
	CALLMETHOD_DEF(XStringSet_unlist, 1),
	CALLMETHOD_DEF(XStringSet_as_STRSXP, 2),
	CALLMETHOD_DEF(XStringSet_order, 1),
	CALLMETHOD_DEF(XStringSet_duplicated, 1),
	CALLMETHOD_DEF(XStringSet_in_set, 2),
	CALLMETHOD_DEF(XStringSet_match, 3),

/* xscat.c */
	CALLMETHOD_DEF(XString_xscat, 1),
	CALLMETHOD_DEF(XStringSet_xscat, 1),

/* fasta_io.c */
	CALLMETHOD_DEF(debug_fasta_io, 0),
	CALLMETHOD_DEF(fasta_info, 2),
	CALLMETHOD_DEF(RawPtr_loadFASTA, 4),
	CALLMETHOD_DEF(fastq_geometry, 1),
	CALLMETHOD_DEF(read_fastq, 2),

/* letter_frequency.c */
	CALLMETHOD_DEF(XString_letter_frequency, 3),
	CALLMETHOD_DEF(XStringSet_letter_frequency, 4),
	CALLMETHOD_DEF(XString_oligo_frequency, 7),
	CALLMETHOD_DEF(XStringSet_oligo_frequency, 8),
	CALLMETHOD_DEF(XStringSet_nucleotide_frequency_at, 7),
	CALLMETHOD_DEF(XStringSet_consensus_matrix, 5),

/* char_translate.c */
	CALLMETHOD_DEF(XStringSet_char_translate, 3),

/* replace_letter_at.c */
	CALLMETHOD_DEF(XString_replace_letter_at, 6),
	CALLMETHOD_DEF(XString_inplace_replace_letter_at, 4),

/* inject_code.c */
	CALLMETHOD_DEF(inject_code, 4),

/* SparseList_utils.c */
	CALLMETHOD_DEF(debug_SparseList_utils, 0),

/* Dups_utils.c */
	CALLMETHOD_DEF(debug_Dups_utils, 0),
	CALLMETHOD_DEF(Dups_diff, 2),

/* match_reporting.c */
	CALLMETHOD_DEF(debug_match_reporting, 0),

/* MIndex_utils.c */
	CALLMETHOD_DEF(debug_MIndex_utils, 0),
	CALLMETHOD_DEF(ByPos_MIndex_endIndex, 3),
	CALLMETHOD_DEF(SparseMIndex_endIndex, 4),
	CALLMETHOD_DEF(ByPos_MIndex_combine, 1),

/* match_pattern_at.c */
	CALLMETHOD_DEF(debug_match_pattern_at, 0),
	CALLMETHOD_DEF(XString_match_pattern_at, 8),
	CALLMETHOD_DEF(XStringSet_vmatch_pattern_at, 8),

/* match_pattern_boyermoore.c */
	CALLMETHOD_DEF(debug_match_pattern_boyermoore, 0),

/* match_pattern_shiftor.c */
	CALLMETHOD_DEF(debug_match_pattern_shiftor, 0),
	CALLMETHOD_DEF(bits_per_long, 0),

/* match_pattern_indels.c */
	CALLMETHOD_DEF(debug_match_pattern_indels, 0),

/* match_pattern.c */
	CALLMETHOD_DEF(debug_match_pattern, 0),
	CALLMETHOD_DEF(XString_match_pattern, 7),
	CALLMETHOD_DEF(XStringViews_match_pattern, 9),
	CALLMETHOD_DEF(XStringSet_vmatch_pattern, 7),

/* match_BOC.c */
	CALLMETHOD_DEF(debug_match_BOC, 0),
	CALLMETHOD_DEF(match_BOC_preprocess, 12),
	CALLMETHOD_DEF(match_BOC_exact, 16),

/* match_BOC2.c */
	CALLMETHOD_DEF(debug_match_BOC2, 0),
	CALLMETHOD_DEF(match_BOC2_preprocess, 9),
	CALLMETHOD_DEF(match_BOC2_exact, 13),

/* match_PWM.c */
	CALLMETHOD_DEF(PWM_score_starting_at, 4),
	CALLMETHOD_DEF(match_PWM, 5),

/* find_palindromes.c */
	CALLMETHOD_DEF(debug_find_palindromes, 0),
	CALLMETHOD_DEF(find_palindromes, 6),

/* PreprocessedTB_class.c */
	CALLMETHOD_DEF(debug_PreprocessedTB_class, 0),

/* match_pdict_Twobit.c */
	CALLMETHOD_DEF(debug_match_pdict_Twobit, 0),
	CALLMETHOD_DEF(build_Twobit, 3),

/* match_pdict_ACtree.c */
	CALLMETHOD_DEF(debug_match_pdict_ACtree, 0),
	CALLMETHOD_DEF(free_actree_nodes_buf, 0),
	CALLMETHOD_DEF(build_ACtree, 3),
	CALLMETHOD_DEF(ACtree_summary, 1),

/* BAB_class.c */
	CALLMETHOD_DEF(debug_BAB_class, 0),
	CALLMETHOD_DEF(IntegerBAB_new, 1),

/* match_pdict_ACtree2.c */
	CALLMETHOD_DEF(debug_match_pdict_ACtree2, 0),
	CALLMETHOD_DEF(ACtree2_nodebuf_max_nblock, 0),
	CALLMETHOD_DEF(ACtree2_nodeextbuf_max_nblock, 0),
	CALLMETHOD_DEF(ACtree2_print_nodes, 1),
	CALLMETHOD_DEF(ACtree2_summary, 1),
	CALLMETHOD_DEF(ACtree2_build, 5),

/* match_pdict.c */
	CALLMETHOD_DEF(debug_match_pdict, 0),
	CALLMETHOD_DEF(XString_match_pdict, 8),
	CALLMETHOD_DEF(XStringViews_match_pdict, 10),
	CALLMETHOD_DEF(XStringSet_vmatch_pdict, 10),

/* align_utils.c */
	CALLMETHOD_DEF(PairwiseAlignedXStringSet_nmatch, 4),
	CALLMETHOD_DEF(AlignedXStringSet_nchar, 1),
	CALLMETHOD_DEF(AlignedXStringSet_align_aligned, 2),
	CALLMETHOD_DEF(PairwiseAlignedFixedSubject_align_aligned, 3),
	CALLMETHOD_DEF(align_compareStrings, 6),

/* pmatchPattern.c */
	CALLMETHOD_DEF(lcprefix, 6),
	CALLMETHOD_DEF(lcsuffix, 6),

/* align_pairwiseAlignment.c */
	CALLMETHOD_DEF(XStringSet_align_pairwiseAlignment, 14),
	CALLMETHOD_DEF(XStringSet_align_distance, 12),

/* align_needwunsQS.c */
	CALLMETHOD_DEF(align_needwunsQS, 7),

/* strutils.c (belonged originally to old matchprobes package) */
	CALLMETHOD_DEF(MP_rna_revcomp, 1),
	CALLMETHOD_DEF(MP_dna_revcomp, 1),
	CALLMETHOD_DEF(MP_revstring, 1),
	CALLMETHOD_DEF(MP_complementSeq, 3),
	CALLMETHOD_DEF(MP_basecontent, 2),
	CALLMETHOD_DEF(MP_longestConsecutive, 2),

/* matchprobes.c (belonged originally to old matchprobes package) */
	CALLMETHOD_DEF(MP_matchprobes, 3),

	{NULL, NULL, 0}
};

void R_init_Biostrings(DllInfo *info)
{
	/* Lots of code around assumes that sizeof(Rbyte) == sizeof(char) */
	if (sizeof(Rbyte) != sizeof(char))
		error("sizeof(Rbyte) != sizeof(char)");
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);

/* RoSeq_utils.c */
	REGISTER_CCALLABLE(_new_STRSXP_from_RoSeqs);
	REGISTER_CCALLABLE(_new_RoSeqs_from_CharAEAE);
	REGISTER_CCALLABLE(_new_IRanges_from_RoSeqs);

/* XString_class.c */
	REGISTER_CCALLABLE(_DNAencode);
	REGISTER_CCALLABLE(_DNAdecode);
	REGISTER_CCALLABLE(_RNAencode);
	REGISTER_CCALLABLE(_RNAdecode);
	REGISTER_CCALLABLE(_get_XString_asRoSeq);

/* XStringSet_class.c */
	REGISTER_CCALLABLE(_get_XStringSet_baseClass);
	REGISTER_CCALLABLE(_get_XStringSet_length);
	REGISTER_CCALLABLE(_new_CachedXStringSet);
	REGISTER_CCALLABLE(_get_CachedXStringSet_elt_asRoSeq);
	REGISTER_CCALLABLE(_get_XStringSet_elt_asRoSeq);
	REGISTER_CCALLABLE(_new_XStringSet_from_RoSeqs);
	REGISTER_CCALLABLE(_set_XStringSet_names);
	REGISTER_CCALLABLE(_alloc_XStringSet);
	REGISTER_CCALLABLE(_write_RoSeq_to_CachedXStringSet_elt);
	REGISTER_CCALLABLE(_write_RoSeq_to_XStringSet_elt);

/* match_reporting.c */
	REGISTER_CCALLABLE(_init_match_reporting);
	REGISTER_CCALLABLE(_drop_reported_matches);
	REGISTER_CCALLABLE(_shift_match_on_reporting);
	REGISTER_CCALLABLE(_report_match);
	REGISTER_CCALLABLE(_reported_matches_asSEXP);

	return;
}

