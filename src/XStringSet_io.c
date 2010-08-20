/****************************************************************************
 *                       Read/write FASTA/FASTQ files                       *
 *                           Author: Herve Pages                            *
 ****************************************************************************/
#include "Biostrings.h"
#include "IRanges_interface.h"

static int debug = 0;

SEXP debug_XStringSet_io()
{
#ifdef DEBUG_BIOSTRINGS
	debug = !debug;
	Rprintf("Debug mode turned %s in file %s\n",
		debug ? "on" : "off", __FILE__);
#else
	Rprintf("Debug mode not available in file %s\n", __FILE__);
#endif
	return R_NilValue;
}

#define LINEBUF_SIZE 20001
static char errmsg_buf[200];



/****************************************************************************
 *  A. READING/WRITING FASTA FILES                                          *
 ****************************************************************************/

/*
 * Even if the standard FASTA markup is made of 1-letter sympols, we want to
 * be able to eventually support longer sympols.
 */
static const char *FASTA_comment_markup = ";", *FASTA_desc_markup = ">";

typedef struct fasta_loader {
	void (*load_desc)(struct fasta_loader *loader,
			  const cachedCharSeq *dataline);
	void (*load_empty_seq)(struct fasta_loader *loader);
	void (*load_seq_line)(struct fasta_loader *loader,
			      const cachedCharSeq *dataline);
	int nrec;
} FASTAloader;

/*
 * The FASTAINFO loader only loads the lengths (and optionally the names)
 * of the sequences.
 */

static CharAEAE descs_buf;
static IntAE seqlengths_buf;

static void FASTAINFO_load_desc(FASTAloader *loader,
		const cachedCharSeq *dataline)
{
	// This works only because dataline->seq is nul-terminated!
	append_string_to_CharAEAE(&descs_buf, dataline->seq);
	return;
}

static void FASTAINFO_load_empty_seq(FASTAloader *loader)
{
	IntAE_insert_at(&seqlengths_buf, seqlengths_buf.nelt, 0);
	return;
}

static void FASTAINFO_load_seq_line(FASTAloader *loader,
		const cachedCharSeq *dataline)
{
	seqlengths_buf.elts[seqlengths_buf.nelt - 1] += dataline->length;
	return;
}

static FASTAloader new_FASTAINFO_FASTAloader(SEXP use_descs)
{
	FASTAloader loader;

	seqlengths_buf = new_IntAE(0, 0, 0);
	if (LOGICAL(use_descs)[0]) {
		loader.load_desc = &FASTAINFO_load_desc;
		descs_buf = new_CharAEAE(0, 0);
	} else {
		loader.load_desc = NULL;
	}
	loader.load_empty_seq = &FASTAINFO_load_empty_seq;
	loader.load_seq_line = &FASTAINFO_load_seq_line;
	loader.nrec = 0;
	return loader;
}

/*
 * The FASTA_seqbuf loader.
 */

static cachedXVectorList FASTA_seqbuf;
static cachedCharSeq FASTA_seqbuf_elt;
static const int *FASTA_lkup;
static int FASTA_lkup_length;

static void FASTA_load_empty_seq(FASTAloader *loader)
{
	FASTA_seqbuf_elt = get_cachedXRawList_elt(&FASTA_seqbuf, loader->nrec);
	FASTA_seqbuf_elt.length = 0;
	return;
}

static void FASTA_load_seq_line(FASTAloader *loader,
		const cachedCharSeq *dataline)
{
	int i1;

	i1 = FASTA_seqbuf_elt.length;
	FASTA_seqbuf_elt.length += dataline->length;
	/* FASTA_seqbuf_elt.seq is a const char * so we need to cast it to
           char * before we can write to it */
	Ocopy_bytes_to_i1i2_with_lkup(i1, FASTA_seqbuf_elt.length - 1,
		(char *) FASTA_seqbuf_elt.seq, FASTA_seqbuf_elt.length,
		dataline->seq, dataline->length,
		FASTA_lkup, FASTA_lkup_length);
	return;
}

static FASTAloader new_FASTAloader(SEXP ans, SEXP lkup)
{
	FASTAloader loader;

	FASTA_seqbuf = cache_XVectorList(ans);
	if (lkup == R_NilValue) {
		FASTA_lkup = NULL;
	} else {
		FASTA_lkup = INTEGER(lkup);
		FASTA_lkup_length = LENGTH(lkup);
	}
	loader.load_desc = NULL;
	loader.load_empty_seq = &FASTA_load_empty_seq;
	loader.load_seq_line = &FASTA_load_seq_line;
	loader.nrec = 0;
	return loader;
}


/*
 * Ignore empty lines and lines starting with 'FASTA_comment_markup' like in
 * the original Pearson FASTA format.
 */
static const char *parse_FASTA_file(FILE *stream, int *recno,
		int nrec, int skip, FASTAloader *loader)
{
	int FASTA_comment_markup_length, FASTA_desc_markup_length,
	    lineno, load_record, new_record;
	char linebuf[LINEBUF_SIZE];
	cachedCharSeq dataline;

	FASTA_comment_markup_length = strlen(FASTA_comment_markup);
	FASTA_desc_markup_length = strlen(FASTA_desc_markup);
	lineno = 0;
	load_record = -1;
	while (fgets(linebuf, LINEBUF_SIZE, stream) != NULL) {
		lineno++;
		dataline.length = rtrimline(linebuf, -1);
		// dataline.length > LINEBUF_SIZE - 1 should never happen
		if (dataline.length >= LINEBUF_SIZE - 1) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "cannot read line %d, line is too long",
				 lineno);
			return errmsg_buf;
		}
		if (dataline.length == 0)
			continue; // we ignore empty lines
		if (strncmp(linebuf, FASTA_comment_markup,
				FASTA_comment_markup_length) == 0)
			continue; // we ignore comment lines
		dataline.seq = linebuf;
		new_record = strncmp(linebuf, FASTA_desc_markup,
				FASTA_desc_markup_length) == 0;
		if (new_record) {
			load_record = *recno >= skip;
			if (load_record && nrec >= 0 && *recno >= skip + nrec)
				return NULL;
			load_record = load_record && loader != NULL;
			if (load_record && loader->load_desc != NULL) {
				dataline.seq += FASTA_desc_markup_length;
				dataline.length -= FASTA_desc_markup_length;
				loader->load_desc(loader, &dataline);
			}
			if (load_record && loader->load_empty_seq != NULL)
				loader->load_empty_seq(loader);
			if (load_record)
				loader->nrec++;
			(*recno)++;
			continue;
		}
		if (load_record == -1) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "\"%s\" expected at beginning of line %d",
				 FASTA_desc_markup, lineno);
			return errmsg_buf;
		}
		if (load_record && loader->load_seq_line != NULL)
			loader->load_seq_line(loader, &dataline);
	}
	return NULL;
}

/* --- .Call ENTRY POINT --- */
SEXP fasta_info(SEXP efp_list, SEXP nrec, SEXP skip, SEXP use_descs)
{
	int nrec0, skip0, i, recno;
	FASTAloader loader;
	FILE *stream;
	SEXP ans, ans_names;
	const char *errmsg;

	nrec0 = INTEGER(nrec)[0];
	skip0 = INTEGER(skip)[0];
	loader = new_FASTAINFO_FASTAloader(use_descs);
	recno = 0;
	for (i = 0; i < LENGTH(efp_list); i++) {
		stream = R_ExternalPtrAddr(VECTOR_ELT(efp_list, i));
		errmsg = parse_FASTA_file(stream, &recno,
				nrec0, skip0, &loader);
		if (errmsg != NULL)
			error("reading FASTA file %s: %s",
			      STRING_ELT(GET_NAMES(efp_list), i), errmsg_buf);
	}
	PROTECT(ans = new_INTEGER_from_IntAE(&seqlengths_buf));
	if (LOGICAL(use_descs)[0]) {
		PROTECT(ans_names = new_CHARACTER_from_CharAEAE(&descs_buf));
		SET_NAMES(ans, ans_names);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP read_fasta_in_XStringSet(SEXP efp_list, SEXP nrec, SEXP skip,
		SEXP set_names, SEXP elementType, SEXP lkup)
{
	int nrec0, skip0, i, recno;
	SEXP ans_width, ans_names, ans;
	const char *element_type;
	char classname[40]; // longest string will be "DNAStringSet"
	FASTAloader loader;
	FILE *stream;

	nrec0 = INTEGER(nrec)[0];
	skip0 = INTEGER(skip)[0];
	PROTECT(ans_width = fasta_info(efp_list, nrec, skip, set_names));
	PROTECT(ans_names = GET_NAMES(ans_width));
	SET_NAMES(ans_width, R_NilValue);
	element_type = CHAR(STRING_ELT(elementType, 0));
	if (snprintf(classname, sizeof(classname), "%sSet", element_type)
	    >= sizeof(classname))
	{
		UNPROTECT(2);
		error("Biostrings internal error in "
		      "read_fasta_in_XStringSet(): 'elementType' too long");
	}
	PROTECT(ans = alloc_XRawList(classname, element_type, ans_width));
	_set_XStringSet_names(ans, ans_names);
	loader = new_FASTAloader(ans, lkup);
	recno = 0;
	for (i = 0; i < LENGTH(efp_list); i++) {
		stream = R_ExternalPtrAddr(VECTOR_ELT(efp_list, i));
		rewind(stream);
		parse_FASTA_file(stream, &recno, nrec0, skip0, &loader);
	}
	UNPROTECT(3);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP write_XStringSet_to_fasta(SEXP x, SEXP efp_list, SEXP width, SEXP lkup)
{
	cachedXStringSet X;
	int x_length, width0, lkup_length, i, j1, j2, dest_nbytes;
	FILE *stream;
	const int *lkup0;
	SEXP x_names, desc;
	cachedCharSeq X_elt;
	char linebuf[LINEBUF_SIZE];

	X = _cache_XStringSet(x);
	x_length = _get_cachedXStringSet_length(&X);
	stream = R_ExternalPtrAddr(VECTOR_ELT(efp_list, 0));
	width0 = INTEGER(width)[0];
	if (width0 >= LINEBUF_SIZE)
		error("'width' must be <= %d", LINEBUF_SIZE - 1);
	linebuf[width0] = 0;
	if (lkup == R_NilValue) {
		lkup0 = NULL;
		lkup_length = 0;
	} else {
		lkup0 = INTEGER(lkup);
		lkup_length = LENGTH(lkup);
	}
	x_names = get_XVectorList_names(x);
	for (i = 0; i < x_length; i++) {
		if (fputs(FASTA_desc_markup, stream) == EOF)
			error("write error");
		if (x_names != R_NilValue) {
			desc = STRING_ELT(x_names, i);
			if (desc == NA_STRING)
				error("'names(x)' contains NAs");
			if (fputs(CHAR(desc), stream) == EOF)
				error("write error");
		}
		if (fputs("\n", stream) == EOF)
			error("write error");
		X_elt = _get_cachedXStringSet_elt(&X, i);
		for (j1 = 0; j1 < X_elt.length; j1 += width0) {
			j2 = j1 + width0;
			if (j2 > X_elt.length)
				j2 = X_elt.length;
			dest_nbytes = j2 - j1;
			j2--;
			Ocopy_bytes_from_i1i2_with_lkup(j1, j2,
				linebuf, dest_nbytes,
				X_elt.seq, X_elt.length,
				lkup0, lkup_length);
			linebuf[dest_nbytes] = 0;
			if (fputs(linebuf, stream) == EOF)
				error("write error");
			if (fputs("\n", stream) == EOF)
				error("write error");
		}
	}
	return R_NilValue;
}



/****************************************************************************
 *  B. READING/WRITING FASTQ FILES                                          *
 ****************************************************************************/

static const char *FASTQ_line1_markup = "@", *FASTQ_line3_markup = "+";

typedef struct fastq_loader {
	void (*load_seqid)(struct fastq_loader *loader,
			   const cachedCharSeq *dataline);
	void (*load_seq)(struct fastq_loader *loader,
			 const cachedCharSeq *dataline);
	void (*load_qualid)(struct fastq_loader *loader,
			    const cachedCharSeq *dataline);
	void (*load_qual)(struct fastq_loader *loader,
			  const cachedCharSeq *dataline);
	int nrec;
} FASTQloader;

/*
 * The FASTQGEOM loader keeps only the common width of the sequences
 * (NA if the set of sequences is not rectangular).
 */

static int FASTQ_width;

static void FASTQGEOM_load_seq(FASTQloader *loader,
		const cachedCharSeq *dataline)
{
	if (loader->nrec == 0) {
		FASTQ_width = dataline->length;
		return;
	}
	if (FASTQ_width == NA_INTEGER || dataline->length == FASTQ_width)
		return;
	FASTQ_width = NA_INTEGER;
	return;
}

static FASTQloader new_FASTQGEOM_FASTQloader()
{
	FASTQloader loader;

	FASTQ_width = NA_INTEGER;
	loader.load_seqid = NULL;
	loader.load_seq = FASTQGEOM_load_seq;
	loader.load_qualid = NULL;
	loader.load_qual = NULL;
	loader.nrec = 0;
	return loader;
}


/*
 * The FASTQ loader.
 */

static cachedXVectorList FASTQ_seqbuf;
static const int *FASTQ_lkup;
static int FASTQ_lkup_length;

static void FASTQ_load_seq(FASTQloader *loader, const cachedCharSeq *dataline)
{
	cachedCharSeq FASTQ_seqbuf_elt;

	FASTQ_seqbuf_elt = get_cachedXRawList_elt(&FASTQ_seqbuf, loader->nrec);
	/* FASTQ_seqbuf_elt.seq is a const char * so we need to cast it to
           char * before we can write to it */
	Ocopy_bytes_to_i1i2_with_lkup(0, FASTQ_seqbuf_elt.length - 1,
		(char *) FASTQ_seqbuf_elt.seq, FASTQ_seqbuf_elt.length,
		dataline->seq, dataline->length,
		FASTQ_lkup, FASTQ_lkup_length);
	return;
}

static FASTQloader new_FASTQloader(SEXP ans, SEXP lkup)
{
	FASTQloader loader;

	FASTQ_seqbuf = cache_XVectorList(ans);
	if (lkup == R_NilValue) {
		FASTQ_lkup = NULL;
	} else {
		FASTQ_lkup = INTEGER(lkup);
		FASTQ_lkup_length = LENGTH(lkup);
	}
	loader.load_seqid = NULL;
	loader.load_seq = FASTQ_load_seq;
	loader.load_qualid = NULL;
	loader.load_qual = NULL;
	loader.nrec = 0;
	return loader;
}


/*
 * Ignore empty lines.
 */
static const char *parse_FASTQ_file(FILE *stream, int *recno,
		int nrec, int skip, FASTQloader *loader)
{
	int FASTQ_line1_markup_length, FASTQ_line3_markup_length,
	    lineno, lineinrecno, load_record;
	char linebuf[LINEBUF_SIZE];
	cachedCharSeq dataline;

	FASTQ_line1_markup_length = strlen(FASTQ_line1_markup);
	FASTQ_line3_markup_length = strlen(FASTQ_line3_markup);
	lineno = lineinrecno = 0;
	while (fgets(linebuf, LINEBUF_SIZE, stream) != NULL) {
		lineno++;
		dataline.length = rtrimline(linebuf, -1);
		// dataline.length > LINEBUF_SIZE - 1 should never happen
		if (dataline.length >= LINEBUF_SIZE - 1) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "cannot read line %d, line is too long",
				 lineno);
			return errmsg_buf;
		}
		if (dataline.length == 0)
			continue; // we ignore empty lines
		dataline.seq = linebuf;
		lineinrecno++;
		if (lineinrecno > 4)
			lineinrecno = 1;
		switch (lineinrecno) {
		    case 1:
			if (strncmp(linebuf, FASTQ_line1_markup,
					FASTQ_line1_markup_length) != 0) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "\"%s\" expected at beginning of line %d",
					 FASTQ_line1_markup, lineno);
				return errmsg_buf;
			}
			load_record = *recno >= skip;
			if (load_record && nrec >= 0 && *recno >= skip + nrec)
				return NULL;
			load_record = load_record && loader != NULL;
			if (load_record && nrec >= 0)
				load_record = *recno < skip + nrec;
			if (load_record && loader->load_seqid != NULL) {
				dataline.seq += FASTQ_line1_markup_length;
				dataline.length -= FASTQ_line1_markup_length;
				loader->load_seqid(loader, &dataline);
			}
		    break;
		    case 2:
			if (load_record && loader->load_seq != NULL)
				loader->load_seq(loader, &dataline);
		    break;
		    case 3:
			if (strncmp(linebuf, FASTQ_line3_markup,
					FASTQ_line3_markup_length) != 0) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "\"%s\" expected at beginning of line %d",
					 FASTQ_line3_markup, lineno);
				return errmsg_buf;
			}
			if (load_record && loader->load_qualid != NULL) {
				dataline.seq += FASTQ_line3_markup_length;
				dataline.length -= FASTQ_line3_markup_length;
				loader->load_qualid(loader, &dataline);
			}
		    break;
		    case 4:
			if (load_record && loader->load_qual != NULL)
				loader->load_qual(loader, &dataline);
			if (load_record)
				loader->nrec++;
			(*recno)++;
		    break;
		}
	}
	return NULL;
}

/* --- .Call ENTRY POINT --- */
SEXP fastq_geometry(SEXP efp_list, SEXP nrec, SEXP skip)
{
	int nrec0, skip0, i, recno;
	FASTQloader loader;
	FILE *stream;
	const char *errmsg;
	SEXP ans;

	nrec0 = INTEGER(nrec)[0];
	skip0 = INTEGER(skip)[0];
	loader = new_FASTQGEOM_FASTQloader();
	recno = 0;
	for (i = 0; i < LENGTH(efp_list); i++) {
		stream = R_ExternalPtrAddr(VECTOR_ELT(efp_list, i));
		errmsg = parse_FASTQ_file(stream, &recno,
				nrec0, skip0, &loader);
		if (errmsg != NULL)
			error("reading FASTQ file %s: %s",
			      STRING_ELT(GET_NAMES(efp_list), i), errmsg_buf);
	}
	PROTECT(ans = NEW_INTEGER(2));
	INTEGER(ans)[0] = loader.nrec;
	INTEGER(ans)[1] = FASTQ_width;
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT ---
 * 'set_names' is ignored.
 */
SEXP read_fastq_in_XStringSet(SEXP efp_list, SEXP nrec, SEXP skip,
		SEXP set_names, SEXP elementType, SEXP lkup)
{
	int nrec0, skip0, ans_length, i, recno;
	SEXP ans, ans_geom, ans_width;
	const char *element_type;
	char classname[40]; // longest string will be "DNAStringSet"
	FASTQloader loader;
	FILE *stream;

	nrec0 = INTEGER(nrec)[0];
	skip0 = INTEGER(skip)[0];
	PROTECT(ans_geom = fastq_geometry(efp_list, nrec, skip));
	ans_length = INTEGER(ans_geom)[0];
	PROTECT(ans_width = NEW_INTEGER(ans_length));
	if (ans_length != 0) {
		if (INTEGER(ans_geom)[1] == NA_INTEGER) {
			UNPROTECT(2);
			error("read_fastq_in_XStringSet(): FASTQ files with "
			      "variable sequence lengths are not supported yet");
		}
		for (recno = 0; recno < ans_length; recno++)
			INTEGER(ans_width)[recno] = INTEGER(ans_geom)[1];
	}
	element_type = CHAR(STRING_ELT(elementType, 0));
	if (snprintf(classname, sizeof(classname), "%sSet", element_type)
	    >= sizeof(classname))
	{
		UNPROTECT(2);
		error("Biostrings internal error in "
		      "read_fasta_in_XStringSet(): 'elementType' too long");
	}
	PROTECT(ans = alloc_XRawList(classname, element_type, ans_width));
	loader = new_FASTQloader(ans, lkup);
	recno = 0;
	for (i = 0; i < LENGTH(efp_list); i++) {
		stream = R_ExternalPtrAddr(VECTOR_ELT(efp_list, i));
		rewind(stream);
		parse_FASTQ_file(stream, &recno, nrec0, skip0, &loader);
	}
	UNPROTECT(3);
	return ans;
}

