/**
 * @file sam.c
 * @author K. S. Dorman
 *
 * This file is for reading sam files.  It also contains the code for Roshan's
 * model, which ultimately needs to move out.  A lot of variables are
 * hard-coded, and there is lots of testing code.  The safest way to use this
 * code is to look through main and then compile and run it once you have verified
 * that it is doing what you want.  It is not yet set up for ease of use.
 *
 * Ultimate, this code is meant to parse sam files, which have entries like:
@SQ	SN:BP48	LN:2001
@PG	ID:bwa	PN:bwa	VN:0.7.15-r1140	CL:bwa mem -t 8 -a targetsB.fsa Olin.fastq
M00259:19:000000000-AUAFV:1:1101:11435:2052	16	BP18	915	60	48M1I133M	*	0	0	TCATGGCCGTCTTCACTTTCGAGGATGAAATCACCTCCCCCGTGCCCCACTGCCAAACGTTACAATGCTATGAAGGATGCGGATTCTATCACCCCTAAGATTATTGATGACATCAAGAGTGTTGAAATCGTTGGGGGAAACGGTGGTCCTGGAACCATCAAGAAACTCACCATTGTCGAGGG	HHHHHGHHHHHGGHGGHHFHHHHHHGHHHHGFFFFFCGHECFFHHFGFGHHHHHHHHGHHFGGGGDFEDHGHHHHFFEHFHHGGGGGGGHHGHGHHHHFGFGEFHHHHHHFHFGHHHHFFFFEEFEFF?EE?GGGFFHGGFEF9AA0;1BFFGFB3G3BHGFGGF1;1B19D11>1BFA>>1	NM:i:12	MD:Z:41C7C5G1T9C11C2C2C0C23G21A48	AS:i:119	XS:i:0
 *
 * To compile:

gcc -O3 -Wall -pedantic -o sam sam.c sequence.c nuc.c qual.c error.c uthash.h fastq.c align.c io.c -lm -lncurses

 */

#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "sam.h"
#include "fastq.h"
#include "nuc.h"
#include "qual.h"
#include "uthash.h"
#include "io.h"
#include "array.h"

#define	N_FILES	4

double ll_align(sam_entry *se, unsigned char *ref);

char cigar_char[CIGAR_NCHAR] = {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'};

/**
 * Parse a cigar string from a file and fill a cigar object with the contents.
 * A cigar string is made up of a variable number of parts, matches, soft clips,
 * etc., each with a type and length.  The number of parts, called ashes, is
 * automatically determined and allocated appropriately.
 *
 * @param fp		file pointer
 * @param cig_in	pointer to cigar object to allocate and fill
 * @return		error status
 */
int read_cigar(FILE *fp, cigar **cig_in)
{
	int fxn_debug = ABSOLUTE_SILENCE;
	*cig_in = malloc(sizeof **cig_in);
	cigar *cig = *cig_in;
	long int fpos = ftell(fp);
	unsigned int len;
	char c;

	cig->n_ashes = 0;
	while ((c = fgetc(fp)) != '\t' && c != EOF) {
		ungetc(c, fp);
		if (fscanf(fp, "%u", &len) != 1)
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
				"Cigar string: expect number.\n");
		c = fgetc(fp);
		++cig->n_ashes;
	}
	debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "Size of cigar string: %u\n",
								cig->n_ashes);

	cig->ashes = malloc(cig->n_ashes * sizeof *cig->ashes);

	size_t i = 0;
	fseek(fp, fpos, SEEK_SET);
	while ((c = fgetc(fp)) != '\t' && c != EOF) {
		ungetc(c, fp);
		if (fscanf(fp, "%u", &cig->ashes[i].len) != 1)
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
				"Cigar string: exepect number.\n");
		c = fgetc(fp);
		debug_msg(fxn_debug >= DEBUG_II, fxn_debug,
			"Cigar string: %c%u\n", c, len);
		if (c == 'M')
			cig->ashes[i].type = CIGAR_MMATCH;
		else if (c == 'H')
			cig->ashes[i].type = CIGAR_HARD_CLIP;
		else if (c == 'S')
			cig->ashes[i].type = CIGAR_SOFT_CLIP;
		else if (c == 'I')
			cig->ashes[i].type = CIGAR_INSERTION;
		else if (c == 'D')
			cig->ashes[i].type = CIGAR_DELETION;
		else if (c == 'N')
			cig->ashes[i].type = CIGAR_SKIP;
		else if (c == 'P')
			cig->ashes[i].type = CIGAR_PAD;
		else if (c == 'X')
			cig->ashes[i].type = CIGAR_MISMATCH;
		else if (c == '=')
			cig->ashes[i].type = CIGAR_MATCH;
		else
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "Cigar "
				"string:  Expect [MHSIDNPX=]; read %c.\n", c);
		++i;
	}

	return NO_ERROR;
} /* read_cigar */

void write_cigar(FILE *fp, cigar *c)
{
	for (unsigned int i = 0; i < c->n_ashes; ++i)
		fprintf(fp, "%u%c", c->ashes[i].len,
					cigar_char[c->ashes[i].type]);
} /* write_cigar */

/**
 * Read and parse a BAM file.
 *
 * @param fp	file pointer
 * @param s_in	sam object to allocate and fill
 * @return	error status
 */
int read_bam(gzFile fp, sam **s_in)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	char c, c1, c2;
	unsigned int i = 0;
	size_t nchar = 0;
	size_t schar = 0, qchar = 0;
	size_t rchar = 0;
	int32_t i32 = 0, block_size;
	uint32_t ui32 = 0;
	z_off_t cpos, cpost;
	char magic[4] = "";
	sam *s;

	*s_in = malloc(sizeof **s_in);
	if (!*s_in)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "sam");
	s = *s_in;

	s->n_se = 0;
	s->n_ref = 0;

	if (gzread(fp, magic, 4) != 4) 
		return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "magic string");
	else
		fprintf(stderr, "Read magick: %.4s\n", magic);
	if (strncmp(magic, "BAM\1", 4))
		return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "not bam files");
	if (gzread(fp, &i32, sizeof(i32)) != sizeof(i32))/* read length of BAM header */
		return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "l_text");
	else
		fprintf(stderr, "Read l_text = %d\n", i32);
	if (gzseek(fp, i32, SEEK_CUR) < 0)/* skip BAM header */
		return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "header");
	if (gzread(fp, &i32, sizeof(i32)) != sizeof(i32))/* read no. references: n_ref */
		return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "n_ref");
	s->n_ref = i32;
	fprintf(stderr, "Number of references: %u\n", s->n_ref);
	cpos = gztell(fp);
	for (uint32_t i = 0; i < s->n_ref; ++i) {
		/* length of reference sequence name */
		if (gzread(fp, &i32, sizeof(i32)) != sizeof(i32))
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "l_name");
		fprintf(stderr, "Reference %u name is size %u\n", i, i32);
		rchar += i32 - 1;
		/* skip ref. name & length */
		if (gzseek(fp, i32 + sizeof(i32), SEEK_CUR) < 0)
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "name");
	}
	fprintf(stderr, "Number of name chars: %zu\n", rchar);
	if (gzseek(fp, cpos, SEEK_SET) < 0)
		return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "gzseek()");
	++rchar;
	s->ref_names = malloc(rchar * sizeof *s->ref_names);
	rchar = 0;
	for (uint32_t i = 0; i < s->n_ref; ++i) {
		if (gzread(fp, &i32, sizeof(i32)) != sizeof(i32))
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "l_name");
		fprintf(stderr, "Going to read reference name %u of size %u\n", i, i32);
		if (gzread(fp, &s->ref_names[rchar], i32) != i32)
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "name");
		if (gzseek(fp, sizeof(i32), SEEK_CUR) < 0)/* ignore l_ref */
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "name");
		fprintf(stderr, "Read reference %.*s\n", i32 - 1, &s->ref_names[rchar]);
		rchar += i32 - 1;
	}

	/* count number of reads */
	cpos = gztell(fp);
	while (gzread(fp, &block_size, sizeof(block_size)) == sizeof(block_size)) {
		cpost = gztell(fp);
		if (gzseek(fp, 2*sizeof(i32), SEEK_CUR) < 0)/* skip refID, pos */
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "block_size");
		if (gzread(fp, &ui32, sizeof(ui32)) != sizeof(ui32))
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "bin_mq_nl");
		nchar += ui32 & 0xFF;
		if (gzseek(fp, sizeof(ui32), SEEK_CUR) < 0)
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "skip flag_nc");
		if (gzread(fp, &i32, sizeof(i32)) != sizeof(i32))
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "bin_mq_nl");
		schar += i32;
		if (gzseek(fp, cpost, SEEK_SET) < 0)
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "gzseek()");
		if (gzseek(fp, block_size, SEEK_CUR) < 0)
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "skip block_size");
		++s->n_se;
		fprintf(stderr, "Read %zu at file position %zu\n", s->n_se, cpost);
	}
	fprintf(stderr, "Number of reads: %zu (name characters: %zu; sequence characters: %zu)\n", s->n_se, nchar, schar);
exit(0);
	/* WORKING: this is way too slow; perhaps we need htslib? */

} /* read_bam */

/**
 * Read and parse a sam file.
 *
 * @param fp	file pointer
 * @param s_in	sam object to allocate and fill
 * @return	error status
 */
int read_sam(FILE *fp, sam **s_in)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	char c, c1, c2;
	unsigned int i = 0;
	size_t nchar = 0;
	size_t schar = 0, qchar = 0;
	size_t rchar = 0;
	sam *s;

	*s_in = malloc(sizeof **s_in);
	if (!*s_in)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "sam");
	s = *s_in;

	s->n_se = 0;
	s->n_ref = 0;

	while ((c = fgetc(fp)) == '@' && c != EOF) {
		c1 = fgetc(fp);
		c2 = fgetc(fp);
		if (c1 == 'S' && c2 == 'Q') {	/* ref sequence */
			++s->n_ref;
			while ((c = fgetc(fp)) != ':' && c != EOF);
								/* SN: */
			while ((c = fgetc(fp)) != '\t' && c != EOF)
				++rchar;			/* ref */
			++rchar;				/* null */
			while ((c = fgetc(fp)) != '\n' && c != EOF);
								/* LN: */
		} else {
			while ((c = fgetc(fp)) != '\n' && c != EOF);
		}
	}
	if (c == EOF)
		return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "premature eof");
	if (!s->n_ref)
		return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "sam file must "
				"include references, i.e. @SQ header lines");
	ungetc(c, fp);
	while (!feof(fp)) {
		++s->n_se;
		if (!(s->n_se % 1000))
			fprintf(stderr, ".");
		while ((c = fgetc(fp)) != '\t' && c != EOF)	/* name */
			++nchar;
		if (c == EOF) {				/* incomplete record */
			--s->n_se;
			break;
		}
		++nchar;	/* null character */
		while ((c = fgetc(fp)) != '\t' && c != EOF);	/* flag */
		while ((c = fgetc(fp)) != '\t' && c != EOF);	/* ref */
		while ((c = fgetc(fp)) != '\t' && c != EOF);	/* pos */
		while ((c = fgetc(fp)) != '\t' && c != EOF);	/* qmap */
		while ((c = fgetc(fp)) != '\t' && c != EOF);	/* cigar */
		while ((c = fgetc(fp)) != '\t' && c != EOF);	/* next */
		while ((c = fgetc(fp)) != '\t' && c != EOF);	/* */
		while ((c = fgetc(fp)) != '\t' && c != EOF);	/* tlen */
		while ((c = fgetc(fp)) != '\t' && c != EOF)	/* sequence */
			++schar;
		while ((c = fgetc(fp)) != '\n' && c != EOF);
	}
	fprintf(stderr, "\n");
	debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "Number of entries: %zu\n",
								s->n_se);
	debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "Ref chars: %zu\n", rchar);
	debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "Name chars: %zu\n", nchar);
	debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "Sequence chars: %zu\n",
								schar);
	rewind(fp);

	s->ref_names = malloc(rchar * sizeof *s->ref_names);
	rchar = 0;
	s->n_ref = 0;
	unsigned int max_ref = 0;
	while ((c = fgetc(fp)) == '@' && c != EOF) {
		c1 = fgetc(fp);
		c2 = fgetc(fp);
		if (c1 == 'S' && c2 == 'Q') {	/* ref sequence */
			while ((c = fgetc(fp)) != ':' && c != EOF);
			if (fscanf(fp, "%s", &s->ref_names[rchar]) != 1)
				return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
							"invalid @SQ format");
			debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
				"Read reference: %s\n", &s->ref_names[rchar]);
			if (max_ref < strlen(&s->ref_names[rchar]))
				max_ref = strlen(&s->ref_names[rchar]);
			++s->n_ref;
			rchar += strlen(&s->ref_names[rchar]) + 1;
			while ((c = fgetc(fp)) != '\n' && c != EOF);
		} else {
			while ((c = fgetc(fp)) != '\n' && c != EOF);
		}
	}

	data_t *sdata = sequence_alloc(schar, nuc_sequence_opt(XY_ENCODING));
	data_t *qdata = sequence_alloc(schar, &_qual_sequence_opt);
	char *cdata = malloc(nchar * sizeof *cdata);
	char *ref_name = malloc((max_ref + 1) * sizeof *ref_name);
	sam_entry *se = malloc(s->n_se * sizeof *se);
	sequence *seqs = malloc(2 * s->n_se * sizeof *seqs);

	s->rchars = malloc(s->n_ref * sizeof *s->rchars);
	rchar = 0;
	for (unsigned int i = 0; i < s->n_ref; ++i) {
		s->rchars[i] = rchar;
		rchar += strlen(&s->ref_names[rchar]) + 1;
	}

	nchar = schar = qchar = 0;
	s->n_se = 0;
	s->n_mapping = 0;
	size_t j = 0;
	while (!feof(fp)) {

		ungetc(c, fp);
		/* read name */
		se[s->n_se].name = &cdata[nchar];
		if (fscanf(fp, "%s", se[s->n_se].name) != 1)
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
						"failure to read name");
		nchar += strlen(se[s->n_se].name) + 1;

		debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
			"Read name: %s\n", se[s->n_se].name);

		/* alignment flag */
		if (fscanf(fp, "%hu", &se[s->n_se].flag) != 1)
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
						"failure to read flag");
		debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
				"Read flag: %u\n", se[s->n_se].flag);

		/* reference name */
		if (fscanf(fp, "%s", ref_name) != 1)
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
					"failure to read reference");
		debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
				"Read reference: %s\n", ref_name);
		if (!(se[s->n_se].flag >> 2 & 1L)) {	/* mapped */
			++s->n_mapping;
			for (unsigned int i = 0; i < s->n_ref; ++i) {
				if (!strcmp(ref_name, s->ref_names + s->rchars[i])) {
					se[s->n_se].ref = i;
					break;
				}
			}
		}

		/* position */
		if (fscanf(fp, "%zu\t", &se[s->n_se].pos) != 1)
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
				"failure to read position");;
		debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
			"Read position: %zu\n", se[s->n_se].pos);

		while ((c = fgetc(fp)) != '\t' && c != EOF);	/* qmap */

		/* cigar */
		if (!(se[s->n_se].flag & (1 << 2))) {
			if (read_cigar(fp, &se[s->n_se].cig))
				return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
						"Invalid cigar string.\n");
		} else {
			while ((c = fgetc(fp)) != '\t' && c != EOF);	/* empty cigar */
		}

		while ((c = fgetc(fp)) != '\t' && c != EOF);	/* rnext */
		while ((c = fgetc(fp)) != '\t' && c != EOF);	/* pnext */
		while ((c = fgetc(fp)) != '\t' && c != EOF);	/* tlen */

		/* sequence */
		se[s->n_se].read = &seqs[j++];
		se[s->n_se].read->seq = &sdata[schar];
		i = 0;
		while ((c = fgetc(fp)) != '\t' && c != EOF)	/* sequence */
			write_char(se[s->n_se].read,
					nuc_sequence_opt(XY_ENCODING), i++,
						(data_t) char_to_xy(c));
		se[s->n_se].read->len = i;
		schar += i;
		debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
			"Read sequence %zu of length %u: ", s->n_se, i);
		debug_call(fxn_debug >= DEBUG_I, fxn_debug, fwrite_nuc_sequence(
				stderr, se[s->n_se].read, XY_ENCODING));
		debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, "\n");

		se[s->n_se].qual = &seqs[j++];
		se[s->n_se].qual->seq = &qdata[qchar];
		i = 0;
		while ((c = fgetc(fp)) != '\t' && c != EOF)	/* qual */
			write_char(se[s->n_se].qual, &_qual_sequence_opt, i++,
/* [KSD,TEMPORARY,WORKING] trying this censoring out */
#ifdef CENSOR_QUALITY_SCORES
					(data_t) char_to_censored_qual(c,
						MIN_ILLUMINA_QUALITY_SCORE,
						MAX_ILLUMINA_QUALITY_SCORE));
#else
					(data_t) char_to_qual(c));
#endif
		se[s->n_se].qual->len = i;
		qchar += i;
		debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
			"Read quality sequene %zu of length %u: ", s->n_se, i);
		debug_call(fxn_debug >= DEBUG_I, fxn_debug,
			fwrite_qual_sequence(stderr, se[s->n_se].qual));
		debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, "\n");

		while ((c = fgetc(fp)) != '\n' && c != EOF);	/* extra columns */
		c = fgetc(fp);

		++s->n_se;
	}
	s->se = se;
	return NO_ERROR;
}/* read_sam */

void write_sam(FILE *fp, sam *s)
{

	for (size_t i = 0; i < s->n_se; ++i)
		write_sam_entry(fp, s, &s->se[i]);
} /* write_sam */

void write_sam_entry(FILE *fp, sam *s, sam_entry *se)
{
	fprintf(fp, "%s\t", se->name);
	fprintf(fp, "%hu\t", se->flag);
	if (!(se->flag & (1 << 2)))
		fprintf(fp, "%s\t", s->ref_names + s->rchars[se->ref]);
	else
		fprintf(fp, "*\t");
	fprintf(fp, "%zu\t", se->pos);
	/* [TODO] no QMAP! */
	fprintf(fp, "0\t");
	if (!(se->flag & (1 << 2))) {
		write_cigar(fp, se->cig);
		fprintf(fp, "\t");
	} else {
		fprintf(fp, "*\t");
	}
	/* [TODO] no rnext, pnext, tlen */
	fprintf(fp, "*\t0\t0\t");
	fwrite_nuc_segment(fp, se->read, XY_ENCODING, 0, se->read->len);
	fprintf(fp, "\t");
	fwrite_qual_sequence(fp, se->qual);

	/* [TODO] no extra columns */
	fprintf(fp, "\n");
} /* write_sam_entry */

/**
 * Once a sam file is read, it can be convenient to hash the entries.  This code
 * can hash in several different ways, passed in via \par sh->type, to faciliate
 * access to aligned reads.
 *
 * 	HASH_REFERENCE: hash aligned reads to the reference sequence they are
 *	aligned to.  Unaligned reads are not hashed.
 *	HASH_READ: hash reads by the read sequence.
 *	HASH_NAME: hash reads by read name.
 *
 * @param s	pointer to sam object
 * @param sh	pointer to sam hash (see uthash.h)
 * @return	error status
 */
int fill_hash(sam *s, sam_hash *sh)
{
	sam_hash *entry = NULL;

	/* finish setting up the per reference indices */
	if (!sh || sh->type & HASH_REFERENCE) {
		size_t *stuff;
		stuff = malloc(s->n_mapping * sizeof **s->ref_list);	/* [TODO] sam::n_mapping contains far more than \sum_i sam::n_per_ref[i] */
		if (!stuff)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"sam::ref_list");
		for (size_t i = 0; i < s->n_ref; ++i) {
			s->ref_list[i] = stuff;
			stuff += s->n_per_ref[i];
			s->n_per_ref[i] = 0;
		}
	}
	for (size_t i = 0; i < s->n_se; ++i) {

		/* hash_sam() has filtered out some reads */
		if (s->se[i].exclude)
			continue;

		/* record this alignment under its ref */
		if (!sh || sh->type & HASH_REFERENCE)
			if (!(s->se[i].flag >> 2 & 1L))	/* mapped */
				s->ref_list[s->se[i].ref][
					s->n_per_ref[s->se[i].ref]++] = i;
		if (sh && sh->type & HASH_NAME)
			HASH_FIND(hh, sh, s->se[i].name,
				strlen(s->se[i].name) * sizeof *s->se[i].name,
									entry);
		else if (sh && sh->type & HASH_READ)
			HASH_FIND(hh, sh, s->se[i].read,
					seqbytes(s->se[i].read,
					nuc_sequence_opt(XY_ENCODING)), entry);
		else if (sh && !(sh->type & HASH_REFERENCE))
			/* no recognized hash type */
			return mmessage(ERROR_MSG, INTERNAL_ERROR,
						"Invalid hash type.\n");
		else	/* HASH_REFERENCE only */
			continue;
		if (!entry)
			return mmessage(ERROR_MSG, INTERNAL_ERROR,
					"Corrupt hash for entry %u.\n", i);

		if (!entry->indices) {
			entry->indices = malloc(entry->count
					* sizeof *entry->indices);
			if (!entry->indices)
				return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"sam_hash::indices");
			entry->indices[0] = i;
			entry->count = 1;
		} else {
			entry->indices[entry->count++] = i;
		}
	}
	return NO_ERROR;
} /* fill_hash */

/**
 * Combine pairs of reads from different files by name.  For example, for
 * Roshan's code, this function is used to combine homeologous pairs from the A
 * and B alignments into single hash.  This function can also filter reads by
 * various criteria, such as length or expected number of errors.
 *
 * @param mh		new hash
 * @param nfiles	number of files to hash
 * @param sds		sam list
 * @param rindex	desired reference index in same file
 * @return		number of reads hashed
 */
size_t hash_merge(merge_hash **mh, unsigned int nfiles, sam **sds, size_t *rindex)
{
	sam_entry *se;
	merge_hash *entry;
	size_t hash_size = 0;

	for (unsigned int j = 0; j < nfiles; ++j) {

		for (size_t i = 0; i < sds[j]->n_per_ref[rindex[j]]; ++i) {
			se = &sds[j]->se[sds[j]->ref_list[rindex[j]][i]];

			if (se->exclude)
				continue;

			HASH_FIND(hh, *mh, se->name, strlen(se->name)
						* sizeof *se->name, entry);

			if (!entry) {
				entry = malloc(sizeof *entry);
				entry->count = calloc(nfiles,
							sizeof *entry->count);
				entry->count[j] = 1;
				entry->indices = calloc(nfiles,
							sizeof *entry->indices);
				entry->nfiles = 1;
				entry->exclude = 0;
				HASH_ADD_KEYPTR(hh, *mh, se->name,
					strlen(se->name) * sizeof *se->name,
									entry);
				++hash_size;
			} else {
				++entry->count[j];
				if (entry->count[j] == 1)
					++entry->nfiles;
			}
		}
	}

	/* set up indices via a second readthrough */
	for (unsigned int j = 0; j < nfiles; ++j) {
		for (size_t i = 0; i < sds[j]->n_per_ref[rindex[j]]; ++i) {
			se = &sds[j]->se[sds[j]->ref_list[rindex[j]][i]];

			if (se->exclude)
				continue;

			HASH_FIND(hh, *mh, se->name, strlen(se->name)
						* sizeof *se->name, entry);

			if (!entry->indices[j]) {
				entry->indices[j] = malloc(entry->count[j]
						* sizeof *entry->indices[j]);
				entry->indices[j][0]
						= sds[j]->ref_list[rindex[j]][i];
				entry->count[j] = 1;
			} else {
				entry->indices[j][entry->count[j]++]
						= sds[j]->ref_list[rindex[j]][i];
			}
		}
	}
/*
	for (size_t i = 0; i < n2; ++i) {
		se = &s2->se[read_set2[i]];

		if (se->exclude)
			continue;

		HASH_FIND(hh, *mh, se->name,
				strlen(se->name) * sizeof *se->name, entry);
		if (!entry->indices[1]) {
			entry->indices[1] = malloc(entry->count[1]
						* sizeof *entry->indices[1]);
			entry->indices[1][0] = read_set2[i];
			entry->count[1] = 1;
		} else {
			entry->indices[1][entry->count[1]++] = read_set2[i];
		}
	}
*/
	return hash_size;
} /* hash_merge */


/**
 * For a merged hash, where multiple alignments per read are combined,
 * reset the alignments to cover the minimum covered region by soft-clipping.
 *
 * @param mh		pointer to the merged hash
 * @param nalign	number of alignments
 * @param sds		sam hashes, one per reference
 * @param start_pos	maxinimum extent of all reads per alignment
 *			(0-base reference index)
 */
int match_soft_clipping(merge_hash *mh, unsigned int nalign, sam **sds,
		size_t *start_pos)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	size_t n_reads = 0;

	debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "Maximum extents:\n");
	for (unsigned int j = 0; j < nalign; ++j)
		debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, "\t%zu-\n",
								start_pos[j]);

	for (merge_hash *me = mh; me != NULL; me = me->hh.next, ++n_reads) {
		unsigned int start_rpos = 0, end_rpos = UINT_MAX;

		if (me->exclude)
			continue;

		/* find minimum alignment extent w/o soft clipping */
		for (unsigned int j = 0; j < nalign; ++j) {
			sam_entry *se = &sds[j]->se[me->indices[j][0]];
			unsigned int rf_idx = se->pos - 1 - start_pos[j];
				/* guaranteed >= 0 */

			if (rf_idx > start_rpos)
				start_rpos = rf_idx;

			for (unsigned int i = 0; i < se->cig->n_ashes; ++i) {
				if (se->cig->ashes[i].type == CIGAR_DELETION
					|| se->cig->ashes[i].type == CIGAR_MATCH
					|| se->cig->ashes[i].type == CIGAR_MMATCH
					|| se->cig->ashes[i].type == CIGAR_MISMATCH
					|| se->cig->ashes[i].type == CIGAR_SKIP)
					rf_idx += se->cig->ashes[i].len;
			}

			if (rf_idx < end_rpos)
				end_rpos = rf_idx;

		}

		debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "Read minimum "
			"extent: %u-%u (0 index on extent: start_pos = %zu)\n",
					 start_rpos, end_rpos, start_pos[0]);

		if (end_rpos <= start_rpos) {
			me->exclude = 1;
			continue;
		}

		/* remake ashes to soft clip to minimum alignment extent */
		for (unsigned int j = 0; j < nalign; ++j) {
			sam_entry *se = &sds[j]->se[me->indices[j][0]];
			unsigned int rf_idx = se->pos - start_pos[j] - 1;
			unsigned int first_ash_nlen = 0, n_ashes = 0;
			unsigned int out = 0, diff_extent = 0;

			for (unsigned int i = 0; i < se->cig->n_ashes; ++i) {
				if (se->cig->ashes[i].type == CIGAR_DELETION
					|| se->cig->ashes[i].type == CIGAR_MATCH
					|| se->cig->ashes[i].type == CIGAR_MMATCH
					|| se->cig->ashes[i].type == CIGAR_MISMATCH
					|| se->cig->ashes[i].type == CIGAR_SKIP)
					rf_idx += se->cig->ashes[i].len;

				debug_msg(fxn_debug >= DEBUG_I, fxn_debug, 
					"Read %s (%zu), alignment %u: rf_idx=%u"
					" (%u - %u)\n", se->name, n_reads, j,
						rf_idx, start_rpos, end_rpos);

				/* --R--[--...--]--... OR --R--]--... */
				if (!n_ashes && rf_idx <= start_rpos) {
					diff_extent = 1;
					++n_ashes;
				} else if (rf_idx <= start_rpos) {
					continue;
				/* [--R--... OR [--R--]--... OR --[--R--... OR --[--R--]--... */
				} else if (!n_ashes && rf_idx > start_rpos
					&& rf_idx <= end_rpos) {
					if (start_rpos > se->pos - start_pos[j] - 1) {
						diff_extent = 1;
						n_ashes += 2;
					} else {
						++n_ashes;
					}
				/* --[--]--R--... OR [--]--R-- */
				} else if (!n_ashes && rf_idx > end_rpos) {
					diff_extent = 1;
					n_ashes += 2 + (start_rpos
						> se->pos - start_pos[j] - 1);
				/* ...--[--R--]--... OR ...--[--R--... */
				} else if (n_ashes && rf_idx <= end_rpos) {
					++n_ashes;
				/* ...--[--]--R--... */
				} else if (n_ashes && rf_idx > end_rpos) {
					diff_extent = 1;
					n_ashes += 2;
				}

				if (rf_idx > end_rpos)
					break;
			}

			debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug,
				"\n\tFile %u old cigar: ", j);
			for (unsigned int i = 0; i < se->cig->n_ashes; ++i)
				debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug,
						"%u%c", se->cig->ashes[i].len,
					cigar_char[se->cig->ashes[i].type]);
			debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, "\n");

			ash *new_ashes = se->cig->ashes;
			if (diff_extent) {
				debug_msg(fxn_debug >= DEBUG_I, fxn_debug, 
					"Readjusting read %zu from %u ashes to "
					"%u ashes\n", n_reads, se->cig->n_ashes,
									n_ashes);
				new_ashes = malloc(n_ashes * sizeof *se->cig);
			}
			rf_idx = se->pos - start_pos[j] - 1;
			n_ashes = 0;
			for (unsigned int i = 0; i < se->cig->n_ashes; ++i) {

				if (se->cig->ashes[i].type == CIGAR_DELETION
					|| se->cig->ashes[i].type == CIGAR_MATCH
					|| se->cig->ashes[i].type == CIGAR_MMATCH
					|| se->cig->ashes[i].type == CIGAR_MISMATCH
					|| se->cig->ashes[i].type == CIGAR_SKIP)
					rf_idx += se->cig->ashes[i].len;

				debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
					"%u%c: rf_idx = %u; n_ashes = %u, "
					"first_ash_nlen = %u (%u-%u)\n",
					se->cig->ashes[i].len,
					cigar_char[se->cig->ashes[i].type],
					rf_idx, n_ashes, first_ash_nlen,
							start_rpos, end_rpos);

				/* read nucleotides 5' of joint cover */
				if (!n_ashes && rf_idx < start_rpos) {

					/* add to 5' soft-clip if read nucs consumed */
					/* adjust alignment position if adding soft clip */
					if (se->cig->ashes[i].type == CIGAR_MATCH
						|| se->cig->ashes[i].type == CIGAR_MMATCH
						|| se->cig->ashes[i].type == CIGAR_MISMATCH) {
						first_ash_nlen += se->cig->ashes[i].len;
						se->pos += se->cig->ashes[i].len;
					} else if (se->cig->ashes[i].type == CIGAR_SOFT_CLIP
						|| se->cig->ashes[i].type == CIGAR_INSERTION) {
						first_ash_nlen += se->cig->ashes[i].len;
					} else if (se->cig->ashes[i].type == CIGAR_DELETION) {
						se->pos += se->cig->ashes[i].len;
					}

				/* next ash starts coincident with joint cover
				 * and read nucs have been consumed: combine all
				 * consumed read nucs into 5' soft clip
				 */
				} else if (!n_ashes && rf_idx == start_rpos && first_ash_nlen) {

					/* add to 5' soft-clip if read nucs consumed */
					if (se->cig->ashes[i].type == CIGAR_MATCH
						|| se->cig->ashes[i].type == CIGAR_MMATCH
						|| se->cig->ashes[i].type == CIGAR_MISMATCH) {
						first_ash_nlen += se->cig->ashes[i].len;
						se->pos += se->cig->ashes[i].len;
					} else if (se->cig->ashes[i].type == CIGAR_SOFT_CLIP
						|| se->cig->ashes[i].type == CIGAR_INSERTION) {
						first_ash_nlen += se->cig->ashes[i].len;
					} else if (se->cig->ashes[i].type == CIGAR_DELETION) {
						se->pos += se->cig->ashes[i].len;
					}
					new_ashes[n_ashes].type = CIGAR_SOFT_CLIP;
					new_ashes[n_ashes].len = first_ash_nlen;
					++n_ashes;

				debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
					"Adding first ash: %uS\n", first_ash_nlen);

				/* next ash starts coincident with joint cover
				 * and there has been no 5' consumption of read
				 * nucs: the only explanation is clipping, which
				 * we retain as is.
				 */
				} else if (!n_ashes && rf_idx == start_rpos) {

					new_ashes[n_ashes].type = se->cig->ashes[i].type;
					new_ashes[n_ashes].len = se->cig->ashes[i].len;
					++n_ashes;

					debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
						"Adding first ash: %u%c\n",
						 se->cig->ashes[i].len,
						 	cigar_char[
							se->cig->ashes[i].type]);

				/* next ash ends inside joint cover:
				 * soft clip includes any consumed read
				 * nucleotides, including part of current
				 * ash before joint cover, if any.
				 */
				} else if (!n_ashes && rf_idx > start_rpos
							&& rf_idx <= end_rpos) {

					/* add to 5' soft-clip if read nucs consumed */
					if (se->cig->ashes[i].type == CIGAR_MATCH
						|| se->cig->ashes[i].type == CIGAR_MMATCH
						|| se->cig->ashes[i].type == CIGAR_MISMATCH ) {
						first_ash_nlen += start_rpos
							+ se->cig->ashes[i].len - rf_idx;
						se->pos += start_rpos
							+ se->cig->ashes[i].len - rf_idx;
					} else if (se->cig->ashes[i].type == CIGAR_SOFT_CLIP
						|| se->cig->ashes[i].type == CIGAR_INSERTION) {
						first_ash_nlen += start_rpos
							+ se->cig->ashes[i].len - rf_idx;
					} else if (se->cig->ashes[i].type == CIGAR_DELETION) {
						se->pos += start_rpos
							+ se->cig->ashes[i].len - rf_idx;
					}
					if (first_ash_nlen) {
						new_ashes[n_ashes].type = CIGAR_SOFT_CLIP;
						new_ashes[n_ashes].len = first_ash_nlen;

						debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
							"Adding first ash: %uS (%u)\n",
								new_ashes[n_ashes].len,
									first_ash_nlen);

						++n_ashes;
					}

					new_ashes[n_ashes].type = se->cig->ashes[i].type;
					new_ashes[n_ashes].len = rf_idx - start_rpos;
					++n_ashes;

					debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
						"Adding first/second ash: %u%c\n",
						rf_idx - start_rpos, cigar_char[
							se->cig->ashes[i].type]);

				/* next ash ends beyond joint cover:
				 */
				} else if (!n_ashes && rf_idx > end_rpos) {

					/* add to 5' soft-clip if read nucs consumed */
					if (se->cig->ashes[i].type == CIGAR_MATCH
						|| se->cig->ashes[i].type == CIGAR_MMATCH
						|| se->cig->ashes[i].type == CIGAR_MISMATCH) {
						first_ash_nlen +=  start_rpos
							+ se->cig->ashes[i].len - rf_idx;
						se->pos += start_rpos
							+ se->cig->ashes[i].len - rf_idx;
					} else if (se->cig->ashes[i].type == CIGAR_INSERTION
						|| se->cig->ashes[i].type == CIGAR_SOFT_CLIP) {
						first_ash_nlen +=  start_rpos
							+ se->cig->ashes[i].len - rf_idx;
					} else if (se->cig->ashes[i].type == CIGAR_DELETION) {
						se->pos += start_rpos
							+ se->cig->ashes[i].len - rf_idx;
					}

					if (first_ash_nlen) {
						new_ashes[n_ashes].type = CIGAR_SOFT_CLIP;
						new_ashes[n_ashes].len = first_ash_nlen;

						debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
							"Adding first ash: %uS (%u)\n", 
							new_ashes[n_ashes].len, first_ash_nlen);

						++n_ashes;
					}

					new_ashes[n_ashes].type = se->cig->ashes[i].type;
					new_ashes[n_ashes].len = end_rpos - start_rpos;

					debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
						"Adding first/second ash: %u%c\n",
						end_rpos - start_rpos,
						cigar_char[se->cig->ashes[i].type]);

					++n_ashes;

					new_ashes[n_ashes].type = CIGAR_SOFT_CLIP;
					new_ashes[n_ashes].len = rf_idx - end_rpos;

					debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
						"Adding second/third ash: %u%c\n",
						rf_idx - end_rpos,
						cigar_char[CIGAR_SOFT_CLIP]);

					++n_ashes;
					out = 1;

				/* next ash is fully contained in joint cover:
				 * copy as is
				 */
				} else if (n_ashes && rf_idx <= end_rpos) {
					new_ashes[n_ashes].type = se->cig->ashes[i].type;
					new_ashes[n_ashes].len = se->cig->ashes[i].len;
					++n_ashes;

					debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
						"Adding fully contained later ash: %u%c\n",
						se->cig->ashes[i].len,
						cigar_char[se->cig->ashes[i].type]);

				/* next ash passes outside joint cover for the
				 * first time:
				 */
				} else if (n_ashes && rf_idx > end_rpos && !out) {
					new_ashes[n_ashes].type = se->cig->ashes[i].type;
					new_ashes[n_ashes].len = end_rpos + se->cig->ashes[i].len - rf_idx;

					debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
						"Adding later split ash: %u%c\n",
						new_ashes[n_ashes].len,
						cigar_char[se->cig->ashes[i].type]);

					++n_ashes;
					new_ashes[n_ashes].type = CIGAR_SOFT_CLIP;
					new_ashes[n_ashes].len = se->cig->ashes[i].len - new_ashes[n_ashes - 1].len;

					debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
						"Adding first version of last ash: %u%c\n",
						new_ashes[n_ashes].len,
						cigar_char[CIGAR_SOFT_CLIP]);

					++n_ashes;
					out = 1;

				/* next ash is still outside joint cover */
				} else if (n_ashes && rf_idx > end_rpos && out
					&& (se->cig->ashes[i].type == CIGAR_INSERTION
						|| se->cig->ashes[i].type == CIGAR_MATCH
						|| se->cig->ashes[i].type == CIGAR_MMATCH
						|| se->cig->ashes[i].type == CIGAR_MISMATCH
						|| se->cig->ashes[i].type == CIGAR_SOFT_CLIP)) {
					new_ashes[n_ashes-1].len += se->cig->ashes[i].len;

					debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
						"Adding to last ash: %u%c\n",
						new_ashes[n_ashes-1].len,
							cigar_char[CIGAR_SOFT_CLIP]);
				}
			}


			if (fxn_debug >= DEBUG_I) {
				debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug,
					"\tFile %u new cigar: ", j);
				for (unsigned int i = 0; i < n_ashes; ++i)
					debug_msg_cont(fxn_debug >= DEBUG_I,
						fxn_debug, "%u%c",
						new_ashes[i].len,
						cigar_char[new_ashes[i].type]);
				debug_msg_cont(fxn_debug >= DEBUG_I,
					fxn_debug, " (n_ashes = %u)\n", n_ashes);
			}


			se->cig->ashes = new_ashes;
			se->cig->n_ashes = n_ashes;
		}
	}

	return NO_ERROR;
} /* match_soft_clip */


/**
 * Hash a sam file.
 *
 * @param s		sam data
 * @param sh_in		new hash
 * @param hash_on	how to hash
 *
 * 	HASH_REFERENCE:	hash aligned reads to the reference sequence they are
 *			aligned to.  Unaligned reads are not hashed.
 *	HASH_READ: hash reads by the read sequence.
 *	HASH_NAME: hash reads by read name.
 *
 * @param rindex	index of desired reference
 *
 * ===FILTERS===
 * [TODO] Move these filters to read_sam, although read_sam does not do much
 * [TODO] interpretation.  Or drop them from memory, but that is not so easy.
 * @param drop_unmapped	filter: drop unmapped reads
 * @param drop_second	filter: discard secondary alignments
 * @param drop_sc	filter: discard soft-clips of minimum length drop_sc
 * @param drop_id	filter: discard indels of minimum length drop_id
 * @param min_len	filter: minimum length (0 does nothing)
 * @param max_len	filter: maximum length (0 does nothing)
 * @param max_exp_err	filter: maximum expected errrors (negative does nothing)
 * @return		error status
 */
int
hash_sam (
	sam *s,
	sam_hash **sh_in,
	int hash_on,
	size_t rindex,
	unsigned char drop_unmapped,
	unsigned char drop_second,
	unsigned char drop_sc,
	unsigned char drop_id,
	unsigned int min_length,
	unsigned int max_length,
	double max_exp_err)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_II;//DEBUG_I;//
	size_t n_unmapped = 0;
	size_t n_secondary = 0;
	size_t n_min_length = 0;
	size_t n_max_length = 0;
	size_t n_experr = 0;
	size_t n_soft_clip = 0;
	size_t n_hard_clip = 0;
	size_t n_indel = 0;
	sam_hash *entry;
	unsigned int len, len_sc, len_hc, len_id;

	s->hash_length = 0;

	/* prepare to record indices of mapped reads under each reference */
	if (hash_on & HASH_REFERENCE) {
		s->ref_list = malloc(s->n_ref * sizeof *s->ref_list);
		if (!s->ref_list)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"sam::ref_list");
		s->n_per_ref = malloc(s->n_ref * sizeof *s->n_per_ref);
		if (!s->n_per_ref)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"sam::n_per_ref");
		for (size_t i = 0; i < s->n_ref; ++i)
			s->n_per_ref[i] = 0;
	}

	for (size_t i = 0; i < s->n_se; ++i) {
		sam_entry *se = &s->se[i];
		se->exclude = 0;

		debug_msg(fxn_debug >= DEBUG_II, fxn_debug, "%zu: ", i);

		if (drop_unmapped && se->flag >> 2 & 1) {
			se->exclude = 1;
			debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, "remove %s b/c unmapped\n", se->name);
			++n_unmapped;
			continue;
		}

		/* discard secondary alignments */
		if (drop_second && se->flag >> 11 & 1) {
			se->exclude = 1;
			if (se->ref == rindex) {
				debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, "remove %s b/c secondary alignment\n", se->name);
				++n_secondary;
			} else {
				debug_msg_cont(fxn_debug >= DEBUG_II, fxn_debug, "remove %s b/c secondary alignment to wrong reference\n", se->name);
			}
			continue;
		}

		if (max_length || min_length || drop_sc || drop_id) {
			len = len_sc = len_hc = len_id = 0;
			for (unsigned j = 0; j < se->cig->n_ashes; ++j) {
				if (se->cig->ashes[j].type == CIGAR_SOFT_CLIP
					&& se->cig->ashes[j].len >= len_sc)
					len_sc = se->cig->ashes[j].len;
				if (se->cig->ashes[j].type == CIGAR_HARD_CLIP
					&& se->cig->ashes[j].len >= len_hc)
					len_hc = se->cig->ashes[j].len;
				if ((se->cig->ashes[j].type == CIGAR_INSERTION
					|| se->cig->ashes[j].type == CIGAR_DELETION)
					&& se->cig->ashes[j].len >= len_id)
					len_id = se->cig->ashes[j].len;
				if (se->cig->ashes[j].type == CIGAR_INSERTION
					|| se->cig->ashes[j].type == CIGAR_MATCH
					|| se->cig->ashes[j].type == CIGAR_MMATCH
					|| se->cig->ashes[j].type == CIGAR_MISMATCH
					|| se->cig->ashes[j].type == CIGAR_SOFT_CLIP)
					len += se->cig->ashes[j].len;
			}
			if (drop_sc && len_sc >= drop_sc) {
				se->exclude = 1;
				if (se->ref == rindex) {
					debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, "remove %s b/c soft clipped %uS\n", se->name, len_sc);
					++n_soft_clip;
				} else {
					debug_msg_cont(fxn_debug >= DEBUG_II, fxn_debug, "remove %s b/c soft clipped on wrong reference\n", se->name);
				}
				continue;

			/* assume if user doesn't want soft-clipping, they also don't want hard-clipping */
			} else if (drop_sc && len_hc >= drop_sc) {
				se->exclude = 1;
				if (se->ref == rindex) {
					debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, "remove %s b/c hard clipped\n", se->name);
					++n_hard_clip;
				} else {
					debug_msg_cont(fxn_debug >= DEBUG_II, fxn_debug, "remove %s b/c hard clipped on wrong reference\n", se->name);
				}
				continue;
			}
			if (drop_id && len_id >= drop_id) {
				se->exclude = 1;
				if (se->ref == rindex) {
					debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, "remove %s b/c of indels\n", se->name);
					++n_indel;
				} else {
					debug_msg_cont(fxn_debug >= DEBUG_II, fxn_debug, "remove %s b/c of indels on wrong reference (%u == %zu)\n", se->name, se->ref, rindex);
				}
				continue;
			}
			if (len > max_length) {
//mmessage(INFO_MSG, NO_ERROR, "Read %s length: %u\n", se->name, len);
				se->exclude = 1;
				if (se->ref == rindex) {
					debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, "remove %s b/c read too long\n", se->name);
					++n_max_length;
				} else {
					debug_msg_cont(fxn_debug >= DEBUG_II, fxn_debug, "remove %s b/c read too long aligned to wrong reference\n", se->name);
				}
				continue;
			}
			if (len < min_length) {
				se->exclude = 1;
				if (se->ref == rindex) {
					debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, "remove %s b/c read too short\n", se->name);
					++n_min_length;
				} else {
					debug_msg_cont(fxn_debug >= DEBUG_II, fxn_debug, "remove %s b/c read too short aligned to wrong reference\n", se->name);
				}
				continue;
			}
		}

		if (max_exp_err > 0) {
			double mee = 0;
			for (unsigned int j = 0; j < se->read->len; ++j)
				mee += qual_to_prob(se->qual, j);
			if (mee > max_exp_err) {
				se->exclude = 1;
				if (se->ref == rindex) {
					debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, "remove %s b/c too many expected errors\n", se->name);
					++n_experr;
				} else {
					debug_msg_cont(fxn_debug >= DEBUG_II, fxn_debug, "remove %s b/c too many expected errors aligned to wrong reference\n", se->name);
				}
				continue;
			}
		}

		/* count number mapped to each reference */
		if (hash_on & HASH_REFERENCE)
			if (!(se->flag >> 2 & 1L)) {
				++s->n_per_ref[se->ref];
				if (se->ref == rindex)
					debug_msg_cont(fxn_debug >= DEBUG_II, fxn_debug, "retain %s\n", se->name);
				else
					debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, "remove %s for not aligning to desired reference\n", se->name);
			}

		if (hash_on & HASH_NAME) {
			HASH_FIND(hh, *sh_in, se->name,
				strlen(se->name) * sizeof *se->name, entry);
			if (!entry) {
				entry = malloc(sizeof *entry);
				entry->idx = i;
				entry->count = 1;
				entry->indices = NULL;
				HASH_ADD_KEYPTR(hh, *sh_in, se->name,
					strlen(se->name) * sizeof *se->name,
									entry);
				++s->hash_length;
			} else {
				++entry->count;
			}
		} else if (hash_on & HASH_READ) {
			HASH_FIND(hh, *sh_in, se->read, seqbytes(se->read,
					nuc_sequence_opt(XY_ENCODING)), entry);
			if (!entry) {
				entry = malloc(sizeof *entry);
				entry->idx = i;
				entry->count = 1;
				entry->indices = NULL;
				HASH_ADD_KEYPTR(hh, *sh_in, se->read,
					seqbytes(se->read,
					nuc_sequence_opt(XY_ENCODING)), entry);
				++s->hash_length;
			} else {
				++entry->count;
			}
		} else if (!(hash_on & HASH_REFERENCE)) {
			return mmessage(ERROR_MSG, INTERNAL_ERROR,
						"invalid hash type\n");
		}
	}

	mmessage(INFO_MSG, NO_ERROR, "%zu total reads (reports below are for "
			"reads mapping to reference %u)\n", s->n_se, rindex);

	if (drop_unmapped)
		mmessage(INFO_MSG, NO_ERROR, "%zu unmapped reads filtered\n",
								n_unmapped);

	if (drop_second)
		mmessage(INFO_MSG, NO_ERROR, "%zu secondary alignments "
						"filtered\n", n_secondary);

	if (min_length || max_length)
		mmessage(INFO_MSG, NO_ERROR, "(%zu, %zu) reads filtered because"
			" length outside range (%u, %u)\n", n_min_length,
				n_max_length, min_length, max_length);

	if (drop_sc) {
		mmessage(INFO_MSG, NO_ERROR, "%zu reads filtered because soft "
			"clip length exceeds %u\n", n_soft_clip,
			(unsigned int) drop_sc);
		if (n_hard_clip)
			mmessage(INFO_MSG, NO_ERROR, "%zu reads filtered "
				"because hard clip length exceeds %u\n",
				n_hard_clip, (unsigned int) drop_sc);
	}

	if (drop_id)
		mmessage(INFO_MSG, NO_ERROR, "%zu reads filtered because indel "
			"length exceeds %u\n", n_indel, (unsigned int) drop_id);

	if (max_exp_err > 0)
		mmessage(INFO_MSG, NO_ERROR, "%zu reads filtered because "
			"expected number of errors exceeded %f\n", n_experr,
								max_exp_err);

	if (hash_on != HASH_REFERENCE)
		(*sh_in)->type = hash_on;
	fill_hash(s, *sh_in);

	return NO_ERROR;
} /* hash_sam */

/**
 * Output error data as:
 * read_position reference_nuc|reference_idx quality read_nuc
 *
 * @param fp	open writable file pointer
 * @param se	alignment entry from sam file (xy_t)
 * @param ref	reference sequence (iupac_t) (if NULL output index)
 * @param lprob	log probability of current alignment (NAN if no output)
 * @param errs	print only errors if true
 * @return	error status
 */
int output_error_data(FILE *fp, sam_entry *se, char_t const *ref, double lprob,
							int errs)
{
//	int fxn_debug = DEBUG_I;//ABSOLUTE_SILENCE;//DEBUG_III;//DEBUG_II;//
	size_t rf_index = se->pos - 1;
	unsigned int rd_index = 0;
	unsigned int num_hc = 0;

	if (errs && !ref) {
		mmessage(ERROR_MSG, INVALID_USER_INPUT, "Cannot output only "
					"errors without reference sequence.\n");
		return 1;
	}

	for (unsigned int i = 0; i < se->cig->n_ashes; ++i) {
		if (se->cig->ashes[i].type == CIGAR_HARD_CLIP) {
			num_hc += se->cig->ashes[i].len;
		} else if (se->cig->ashes[i].type == CIGAR_DELETION) {
			for (unsigned int j = 0; j < se->cig->ashes[i].len; ++j) {
				fprintf(fp, "%s", se->name);
				fprintf(fp, " %d", rd_index + num_hc);
				if (ref)
					fprintf(fp, " %c", iupac_to_char[
							ref[rf_index + j]]);
				fprintf(fp, " %zu", rf_index + j);
				fprintf(fp, " -1 -");
				if (!isnan(lprob))
					fprintf(fp, " %f", lprob);
				fprintf(fp, "\n");
			}
			rf_index += se->cig->ashes[i].len;
		} else if (se->cig->ashes[i].type == CIGAR_SOFT_CLIP) {
			rd_index += se->cig->ashes[i].len;
		} else if (se->cig->ashes[i].type == CIGAR_INSERTION) {
			for (unsigned int j = 0; j < se->cig->ashes[i].len; ++j) {
				fprintf(fp, "%s", se->name);
				fprintf(fp, " %d", rd_index + num_hc + j);
				if (ref)
					fprintf(fp, " -");
				fprintf(fp, " -1");
				fprintf(fp, " %d %c", 
					get_qual(se->qual, rd_index + j),
					xy_to_char[get_nuc(se->read, XY_ENCODING, rd_index + j)]);
				if (!isnan(lprob))
					fprintf(fp, " %f", lprob);
				fprintf(fp, "\n");
			}
			rd_index += se->cig->ashes[i].len;
		} else if (se->cig->ashes[i].type == CIGAR_MATCH
			|| se->cig->ashes[i].type == CIGAR_MMATCH
			|| se->cig->ashes[i].type == CIGAR_MISMATCH) {
			for (unsigned int j = 0; j < se->cig->ashes[i].len; ++j) {
				if (errs && iupac_to_char[ref[rf_index + j]]
					== xy_to_char[get_nuc(se->read,
						XY_ENCODING, rd_index + j)])
					continue;
				fprintf(fp, "%s", se->name);
				fprintf(fp, " %d", rd_index + num_hc + j);
				if (ref)
					fprintf(fp, " %c", iupac_to_char[
							ref[rf_index + j]]);
				fprintf(fp, " %zu", rf_index + j);
				fprintf(fp, " %d %c", 
					get_qual(se->qual, rd_index + j),
					xy_to_char[get_nuc(se->read,
					XY_ENCODING, rd_index + j)]);
				if (!isnan(lprob))
					fprintf(fp, " %f", lprob);
				fprintf(fp, "\n");
			}
			rd_index += se->cig->ashes[i].len;
			rf_index += se->cig->ashes[i].len;
		}
	}

	return NO_ERROR;
} /* output_error_data */

