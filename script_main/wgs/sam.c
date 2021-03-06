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
	
	cig->length_rf = 0;
	for (unsigned int l = 0; l < cig->n_ashes; ++l) {
		if (cig->ashes[l].type == CIGAR_DELETION ||
		    cig->ashes[l].type == CIGAR_MATCH ||
		    cig->ashes[l].type == CIGAR_MISMATCH ||
		    cig->ashes[l].type == CIGAR_MMATCH ||
		    cig->ashes[l].type == CIGAR_SKIP) {
			cig->length_rf += cig->ashes[l].len;
		}
	}
	
	return NO_ERROR;
} /* read_cigar */

/**
 * Read and parse a BAM file.
 *
 * @param fp	file pointer
 * @param s_in	sam object to allocate and fill
 * @return	error status
 */
#ifdef USE_GZIP
int read_bam(gzFile fp, sam **s_in)
{
//	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	size_t nchar = 0;
	size_t schar = 0;
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
	s->rname_ptr = malloc(rchar * sizeof(*s->rname_ptr));
	s->ref_names = malloc(s->n_ref * sizeof(*s->ref_names));
	rchar = 0;
	for (uint32_t i = 0; i < s->n_ref; ++i) {
		if (gzread(fp, &i32, sizeof(i32)) != sizeof(i32))
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "l_name");
		fprintf(stderr, "Going to read reference name %u of size %u\n", i, i32);
		if (gzread(fp, &s->rname_ptr[rchar], i32) != i32)
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "name");
		if (gzseek(fp, sizeof(i32), SEEK_CUR) < 0)/* ignore l_ref */
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "name");
		s->ref_names[i] = &s->rname_ptr[rchar];
		fprintf(stderr, "Read reference %.*s\n", i32 - 1, s->ref_names[i]);
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
#endif

/**
 * Read and parse a sam file.
 *
 * @param fp			file pointer
 * @param s_in			sam object to allocate and fill
 * @param append_strand_char	append mapping strand char to distinguish forward/reverse read
 * @param progress_monitor	progress monitor to stderr
 * @return	error status
 */
int read_sam(FILE *fp, sam **s_in, unsigned char append_strand_char, unsigned char progress_monitor)
{
	int fxn_debug = ABSOLUTE_SILENCE;//append_strand_char ? ABSOLUTE_SILENCE : DEBUG_I;//ABSOLUTE_SILENCE;//
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

	s->ref_list = NULL;
	s->n_per_ref = NULL;
	s->se = NULL;
	s->rname_ptr = NULL;
	s->ref_names = NULL;
	s->n_se = 0;
	s->n_ref = 0;

	while ((c = fgetc(fp)) == '@' && c != EOF) {
		c1 = fgetc(fp);
		c2 = fgetc(fp);
		if (c1 == 'S' && c2 == 'Q') {	/* ref sequence */
			++s->n_ref;
			while ((c = fgetc(fp)) != ':' && c != EOF);
			while ((c = fgetc(fp)) != '\t' && c != EOF)
				++rchar;
			++rchar;
			while ((c = fgetc(fp)) != '\n' && c != EOF);
		} else {
			while ((c = fgetc(fp)) != '\n' && c != EOF);
		}
	}
	if (c == EOF)
		return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "premature eof\n");
	if (!s->n_ref)
		return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "sam file must "
				"include references, i.e. @SQ header lines\n");
	ungetc(c, fp);
	while (!feof(fp)) {
		++s->n_se;
		//if (!(s->n_se % 1000))
		//	fprintf(stderr, ".");
		while ((c = fgetc(fp)) != '\t' && c != EOF)	/* name */
			++nchar;
		nchar += 1 + append_strand_char;	/* null character & +/- for strand */
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
	//fprintf(stderr, "\n");
	debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "Number of entries: %zu\n",
								s->n_se);
	debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "Name chars: %zu\n", nchar);
	debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "Sequence chars: %zu\n",
								schar);
	rewind(fp);

	s->rname_ptr = malloc(rchar * sizeof(*s->rname_ptr));
	s->ref_names = malloc(s->n_ref * sizeof(*s->ref_names));
	rchar = 0;
	s->n_ref = 0;
	unsigned int max_ref = 0;
	while ((c = fgetc(fp)) == '@' && c != EOF) {
		c1 = fgetc(fp);
		c2 = fgetc(fp);
		if (c1 == 'S' && c2 == 'Q') {	/* ref sequence */
			while ((c = fgetc(fp)) != ':' && c != EOF);
			if (fscanf(fp, "%s", &s->rname_ptr[rchar]) != 1)
				return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
							"invalid @SQ format");
			debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
				"Read reference: %s\n", &s->rname_ptr[rchar]);
			if (max_ref < strlen(&s->rname_ptr[rchar]))
				max_ref = strlen(&s->rname_ptr[rchar]);
			s->ref_names[s->n_ref] = &s->rname_ptr[rchar];
			rchar += strlen(s->ref_names[s->n_ref]) + 1;
			s->n_ref++;
			while ((c = fgetc(fp)) != '\n' && c != EOF);
		} else {
			while ((c = fgetc(fp)) != '\n' && c != EOF);
		}
	}

	data_t *sdata = sequence_alloc(schar, nuc_sequence_opt(XY_ENCODING));
	data_t *qdata = sequence_alloc(schar, &_qual_sequence_opt);
	char *cdata = malloc(nchar * sizeof(*cdata));
	char *ref_name = malloc((max_ref + 1) * sizeof(*ref_name));
	sam_entry *se = malloc(s->n_se * sizeof(*se));
	sequence *seqs = malloc(2 * s->n_se * sizeof(*seqs));
	nchar = schar = qchar = 0;
	s->n_se = 0;
	s->n_mapping = 0;
	size_t j = 0;
	while (!feof(fp)) {
		if (progress_monitor && !(s->n_se % 1000))
			fprintf(stderr, ".");
		ungetc(c, fp);
		/* read name */
		se[s->n_se].name = &cdata[nchar];
		if (fscanf(fp, "%s", se[s->n_se].name) != 1)
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
						"failure to read name");
		nchar += strlen(se[s->n_se].name) + 1 + append_strand_char;/* strand indicator */
		debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
			"Read name: %s\n", se[s->n_se].name);

		/* alignment flag */
		if (fscanf(fp, "%hu", &se[s->n_se].flag) != 1)
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
						"failure to read flag");
		debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
				"Read flag: %u\n", se[s->n_se].flag);

		if (append_strand_char && se[s->n_se].flag & 16UL)
			strcat(se[s->n_se].name, "-");
		else if (append_strand_char)
			strcat(se[s->n_se].name, "+");

		/* reference name */
		if (fscanf(fp, "%s", ref_name) != 1)
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
					"failure to read reference");
		debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
				"Read reference: %s\n", ref_name);
		if (!(se[s->n_se].flag >> 2 & 1U)) {	/* mapped */
			++s->n_mapping;
			rchar = 0;
			for (unsigned int i = 0; i < s->n_ref; ++i) {
				if (!strcmp(ref_name, s->ref_names[i])) {
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

		se[s->n_se].exclude = 0;	/* [KSD] avoid uninitialised value error */

		while ((c = fgetc(fp)) != '\n' && c != EOF);	/* extra columns */
		c = fgetc(fp);

		++s->n_se;
	}
	if (progress_monitor)
		fprintf(stderr, "\n");
	s->se = se;
	return NO_ERROR;
}/* read_sam */

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
 * @param n_ref	no. of total references
 * @return	error status
 */
//NOTICE:change se.ref to re.which_ref for the targeted location
int fill_hash(sam *s, sam_hash *sh, unsigned int n_ref)
{
	sam_hash *entry = NULL;

	/* finish setting up the per reference indices */
	if (!sh || sh->type & HASH_REFERENCE) {
		size_t *stuff = NULL;
		size_t nsize = 0;

		for (size_t i = 0; i < n_ref; ++i)
			nsize += s->n_per_ref[i];

		stuff = malloc(nsize * sizeof(**s->ref_list));

		if (!stuff)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"sam::ref_list");
		for (size_t i = 0; i < n_ref; ++i) {
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
				s->ref_list[s->se[i].which_ref][
					s->n_per_ref[s->se[i].which_ref]++] = i;
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
 * Combine reads appearing in different sam files by name.  For example, for
 * Roshan's code, this function is used to combine homeologous pairs from the A
 * and B alignments into single hash.
 *
 * @param mh		new hash
 * @param nfiles	number of files to hash
 * @param sds		sam list
 * @param rindex	optional, only hash reads previously hashed to this
 *			reference via previous HASH_REFERENCE call to hash_sam()
 * @return		number of reads hashed
 */
size_t hash_merge(merge_hash **mh, unsigned int nfiles, sam **sds, size_t *rindex)
{
	sam_entry *se;
	merge_hash *entry;
	size_t hash_size = 0;

	for (unsigned int j = 0; j < nfiles; ++j) {
		size_t n_reads = rindex ? sds[j]->n_per_ref[rindex[j]]
			: sds[j]->n_se;

		for (size_t i = 0; i < n_reads; ++i) {
			se = rindex
				? &sds[j]->se[sds[j]->ref_list[rindex[j]][i]]
				: &sds[j]->se[i];

			if (se->exclude)
				continue;
			
			HASH_FIND(hh, *mh, se->name, strlen(se->name)
						* sizeof(*se->name), entry);
			
			if (!entry) {
				entry = malloc(sizeof(*entry));
				entry->count = calloc(nfiles,
							sizeof(*entry->count));
				entry->count[j] = 1;
				entry->indices = calloc(nfiles,
							sizeof(*entry->indices));
				entry->nfiles = 1;
				entry->exclude = 0;
				HASH_ADD_KEYPTR(hh, *mh, se->name,
					strlen(se->name) * sizeof(*se->name),
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
						* sizeof(*se->name), entry);
			
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

/* reverse in plcae */
void reverse_in_place(sam_entry *se) {
	for (unsigned int i = 0; i < se->cig->n_ashes/2; ++i) {
		unsigned int tmp_type = se->cig->ashes[se->cig->n_ashes-i-1].type;
		se->cig->ashes[se->cig->n_ashes-i-1].type = se->cig->ashes[i].type;
		se->cig->ashes[i].type = tmp_type;
		unsigned int tmp_len = se->cig->ashes[se->cig->n_ashes-i-1].len;
		se->cig->ashes[se->cig->n_ashes-i-1].len = se->cig->ashes[i].len;
		se->cig->ashes[i].len = tmp_len;
	}
} /* reverse_in_place */


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
 * @param nref		optional number of references in external reference list
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
	size_t n_ref,
	unsigned char drop_unmapped,
	unsigned char drop_second,
	unsigned int drop_sc,
	unsigned int drop_id,
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

	/* if an external reference list has not been set up
	 * prepare to record indices of mapped reads under each reference
	 */
	if (hash_on & HASH_REFERENCE) {
		s->ref_list = malloc((n_ref || s->n_ref) * sizeof *s->ref_list);
		if (!s->ref_list)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"sam::ref_list");
		s->n_per_ref = malloc((n_ref || s->n_ref)
							* sizeof *s->n_per_ref);
		if (!s->n_per_ref)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"sam::n_per_ref");
		for (size_t i = 0; i < (n_ref || s->n_ref); ++i)
			s->n_per_ref[i] = 0;
	}

	for (size_t i = 0; i < s->n_se; ++i) {
		sam_entry *se = &s->se[i];
        
		/* KSD, TODO Repace with earlier check on previously set exclude. */
		/* discard targeted unmapped */
		if (se->exclude == 1)
			continue;
        
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
			debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, "remove %s secondary alignment\n", se->name);
			++n_secondary;
			continue;
		}

		if (max_length || min_length || drop_sc < UINT_MAX || drop_id < UINT_MAX) {
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
			if (drop_sc < UINT_MAX && len_sc >= drop_sc) {
				se->exclude = 1;
				debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, "remove %s b/c soft clipped %uS\n", se->name, len_sc);
				++n_soft_clip;
				continue;

			/* assume if user doesn't want soft-clipping, they also don't want hard-clipping */
			} else if (drop_sc < UINT_MAX && len_hc >= drop_sc) {
				se->exclude = 1;
				debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, "remove %s b/c hard clipped\n", se->name);
				++n_hard_clip;
				continue;
			}
			if (drop_id < UINT_MAX && len_id >= drop_id) {
				se->exclude = 1;
				debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, "remove %s b/c of indels\n", se->name);
				++n_indel;
				continue;
			}
			if (max_length && len > max_length) {
				se->exclude = 1;
				debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, "remove %s b/c read too long\n", se->name);
				++n_max_length;
				continue;
			}
			if (len < min_length) {
				se->exclude = 1;
				debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, "remove %s b/c read too short\n", se->name);
				++n_min_length;
				continue;
			}
		}

		if (max_exp_err > 0) {
			double mee = 0;
			for (unsigned int j = 0; j < se->read->len; ++j)
				mee += qual_to_prob(se->qual, j);
			if (mee > max_exp_err) {
				se->exclude = 1;
				debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, "remove %s alignment with too many expected errors\n", se->name);
				++n_experr;
				continue;
			}
		}

		/* count number mapped to each reference */
		/* for external reference indicated by argument n_ref > 0
		 * the sam_entry::which_ref must be already set in range
		 * {0,1,...,n_ref-1}
		 */
		if (hash_on & HASH_REFERENCE)
			if (!(se->flag >> 2 & 1L)) {
				++s->n_per_ref[se->which_ref];
				debug_msg_cont(fxn_debug >= DEBUG_II, fxn_debug, "retain %s\n", se->name);
			}

		if (hash_on & HASH_NAME) { // notice, consider reverse or forward strand
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

	mmessage(INFO_MSG, NO_ERROR, "%zu total reads\n", s->n_se);

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

	if (drop_sc < UINT_MAX) {
		mmessage(INFO_MSG, NO_ERROR, "%zu reads filtered because soft "
			"clip length exceeds %u\n", n_soft_clip,
			(unsigned int) drop_sc);
		if (n_hard_clip)
			mmessage(INFO_MSG, NO_ERROR, "%zu reads filtered "
				"because hard clip length exceeds %u\n",
				n_hard_clip, (unsigned int) drop_sc);
	}

	if (drop_id < UINT_MAX)
		mmessage(INFO_MSG, NO_ERROR, "%zu reads filtered because indel "
			"length exceeds %u\n", n_indel, drop_id);

	if (isfinite(max_exp_err))
		mmessage(INFO_MSG, NO_ERROR, "%zu reads filtered because "
			"expected number of errors exceeded %f\n", n_experr,
								max_exp_err);

	if (hash_on != HASH_REFERENCE)
		(*sh_in)->type = hash_on;

	fill_hash(s, *sh_in, n_ref);

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
 * @return	error status
 */
int output_error_data(FILE *fp, sam_entry *se, unsigned char *ref, double lprob)
{
//	int fxn_debug = DEBUG_I;//ABSOLUTE_SILENCE;//DEBUG_III;//DEBUG_II;//
	size_t rf_index = se->pos - 1;
	unsigned int rd_index = 0;
	unsigned int num_hc = 0;

	for (unsigned int i = 0; i < se->cig->n_ashes; ++i) {
		if (se->cig->ashes[i].type == CIGAR_HARD_CLIP) {
			num_hc += se->cig->ashes[i].len;
		} else if (se->cig->ashes[i].type == CIGAR_DELETION) {
			for (unsigned int j = 0; j < se->cig->ashes[i].len; ++j) {
				fprintf(fp, "%d ", se->index);
				fprintf(fp, "%d", rd_index + num_hc);
				if (ref)
					fprintf(fp, " %c", iupac_to_char[
							ref[rf_index + j]]);
				else
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
				fprintf(fp, "%d", rd_index + num_hc + j);
				if (ref)
					fprintf(fp, " -");
				else
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
				fprintf(fp, "%d", rd_index + num_hc + j);
				if (ref)
					fprintf(fp, " %c", 
						iupac_to_char[ref[rf_index + j]]);
				else
					fprintf(fp, " %zu", rf_index + j);
				fprintf(fp, " %d %c", 
					get_qual(se->qual, rd_index + j),
					xy_to_char[get_nuc(se->read, XY_ENCODING, rd_index + j)]);
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

void print_cigar(FILE *file, sam_entry *se)
{
	for (unsigned int i = 0; i < se->cig->n_ashes; ++i)
		fprintf(file, "%u%c", se->cig->ashes[i].len,
			cigar_char[se->cig->ashes[i].type]);
} /* print_cigar */
