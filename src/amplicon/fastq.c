#include <stdlib.h>
#include <limits.h>

#include "fastq.h"
#ifndef NO_ALIGNMENT
#include "align.h"
#endif
#include "math.h"
#include "error.h"

char const *file_type_name[NUM_FILE_TYPES] = {"fastq", "fasta"};

/**
 * Macro forward to the next newline or EOF.
 */
#define fforward(f, c, t) do {                                                 \
	(c) = fgetc((f));                                                      \
} while ((c) != EOF && (c) != (t));

/**
 * Macro forward to the next newline or EOF, and count characters consumed.
 */
#define fforward_cnt(f, c, t, a)                                               \
while (((c) = fgetc((f))) != EOF && (c) != (t)) {                              \
	(a)++;                                                                 \
}

/**
 * Macro forward to the next newline, record and count characters.
 */
#define fforward_copy(f, c, t, s, a)                                           \
while (((c) = fgetc((f))) != EOF && (c) != (t)) {                              \
	*((s) + (a)) = c;                                                      \
	(a)++;                                                                 \
}

/**
 * Convert quality score to probability.
 */
extern double error_prob(fastq_data *fqd, char_t c);

/**
 * Print all observed quality scores as probabilities.
 */
extern void fprint_error_probs(FILE *fp, fastq_data *fqd);

/**
 * Write one read of given length to file in R table format (space-separated
 * integers).
 */
extern void write_read_in_table(FILE *fp, char_t *read, unsigned int len);

/**
 * Return length of requested read.
 */
extern unsigned int read_length(fastq_data *fqd, unsigned int i);

/**
 * Return ASCII nucleotide for encoded character c, using encoding of fqd.
 */
extern inline char nuc(fastq_data *fqd, char_t c);


/**
 * Open fastq file and count number of reads.
 *
 * @param filename	fastq filename
 * @return		number of reads (0 for failure)
 */
unsigned int cnt_reads(char const * const filename)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	FILE *fp = fopen(filename, "r");

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "opening fastq file '%s'\n",
								 filename);

	if (fp == NULL) {
		mmessage(ERROR_MSG, FILE_OPEN_ERROR, filename);
		return 0;
	}

	unsigned int cnt = fcnt_reads(fp);

	fclose(fp);

	return cnt;
} /* cnt_reads */

unsigned int fcnt_reads(FILE *fp)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	int file_type = FASTQ_FILE;
	unsigned int nread = 0;
	char c;

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "entering\n");

	c = fgetc(fp);
	if (c != '@' && c != '>') {
		mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "not fast[aq] file\n");
		return 0;
	}

	if (c == '>')
		file_type = FASTA_FILE;

	do {

		if (nread) {
			c = fgetc(fp);
			if (c == EOF)
				return nread;
			if (c != '@' && c != '>') {
				mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
					"not fast[aq] file\n");
				return 0;
			}
		}

		/* fast forward through name */
		fforward(fp, c, '\n');

		if (c == EOF) {
			mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
				"premature end-of-file\n");
			return 0;
		}

		/* fast forward through read */
		fforward(fp, c, '\n');

		if (c == EOF) {
			mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
				"premature end-of-file\n");
			return 0;
		}

		if (file_type == FASTA_FILE)
			continue;

		/* fast forward through divider */
		fforward(fp, c, '\n');

		if (c == EOF) {
			mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
				"premature end-of-file\n");
			return 0;
		}

		/* fast forward through quality score */
		fforward(fp, c, '\n');

		nread++;
	} while (c != EOF);

	return nread;

} /* fcnt_reads */

/**
 * Index fast[aq] file.  Record byte index of the start of each read.
 *
 * @param fp	fast[aq] open file handle
 * @param fqd	fastq_data object pointer, with fastq_data.n_reads set
 * @return	error status
 */
int findex_reads(FILE *fp, fastq_data *fqd)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_II;//
	unsigned int i = 0;
	size_t nbytes = 0;
	char c;

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "entering\n");

	c = fgetc(fp);
	if (c != '@' && c != '>')
		return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
			"not fast[aq] file\n");

	fqd->file_type = FASTQ_FILE;
	if (c == '>')
		fqd->file_type = FASTA_FILE;

	fqd->index = malloc(fqd->n_reads * sizeof *fqd->index);

	if (!fqd->index)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"fastq_data.index");

	do {

		if (i) {
			c = fgetc(fp);
			if (c == EOF)
				return NO_ERROR;
			if (c != '@' && c != '>')
				return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
					"not fast[aq] file\n");
		}

		fqd->index[i++] = nbytes;

		debug_msg(DEBUG_II <= fxn_debug, fxn_debug, "Read %u at byte "
							"%u.\n", i, nbytes);

		++nbytes;

		/* fast forward through name */
		fforward_cnt(fp, c, '\n', nbytes);

		if (c == EOF)
			return FASTQ_EOF;
		++nbytes;	/* for newline */

		/* fast forward through read */
		fforward_cnt(fp, c, '\n', nbytes);

		if (c == EOF)
			return FASTQ_EOF;
		++nbytes;

		if (fqd->file_type == FASTA_FILE)
			continue;

		/* fast forward through divider */
		fforward_cnt(fp, c, '\n', nbytes);

		if (c == EOF)
			return FASTQ_EOF;
		++nbytes;

		/* fast forward through quality score */
		fforward_cnt(fp, c, '\n', nbytes);
		++nbytes;

	} while (c != EOF);

	return NO_ERROR;

} /* findex_reads */

/**
 * Open, read and store fastq file in newly allocated memory.
 *
 * @param filename	fastq filename
 * @param in_fqd	unallocated fastq object
 * @param fqo		fastq options pointer
 * @return		error code
 */
int read_fastq(const char *filename, fastq_data **in_fqd, fastq_options *fqo)
{
	int fxn_debug = ABSOLUTE_SILENCE;//SILENT;//DEBUG_III;//
	int err = NO_ERROR;
	FILE *fp = fopen(filename, "r");

	if (fxn_debug)
		fprintf(stderr, "%s:%d: opening fastq file '%s'\n", __func__,
			__LINE__, filename);

	if (fp == NULL)
		return mmessage(ERROR_MSG, FILE_OPEN_ERROR, filename);

	err = fread_fastq(fp, in_fqd, fqo);

	fclose(fp);

	return err;
} /* read_fastq */

int fread_fastq(FILE *fp, fastq_data **in_fqd, fastq_options *fqo)
{
	int fxn_debug = DEBUG_I;//ABSOLUTE_SILENCE;//SILENT;//DEBUG_III;//
	char_t *rptr, *qptr;
	unsigned int n_reads, n_read_bytes, n_last_read_bytes, n_bytes, *uiptr;
	unsigned int n_name_bytes;
	unsigned int n_total_name_bytes;
	char *nptr = NULL;
	unsigned int *nlptr = NULL;
	unsigned char c;
	unsigned char elen = 1;	/* equal length reads */
	fastq_data *fqd = *in_fqd;
	int err = NO_ERROR;

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "entering\n");

	c = fgetc(fp);
	if (c != '@' && c != '>')
		return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "not fast[aq] file\n");

	rewind(fp);

	if (fqd == NULL) {
		*in_fqd = malloc(sizeof **in_fqd);
		if (*in_fqd == NULL)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"fastq_data reads");
		fqd = *in_fqd;
		fqd->file_type = FASTQ_FILE;
		fqd->n_lengths = NULL;
		fqd->reads = NULL;
		fqd->quals = NULL;
		fqd->names = NULL;
		fqd->name_lengths = NULL;
		fqd->reference_seq = NULL;
		fqd->index = NULL;
	}

	if (c == '>') {
		fqd->file_type = FASTA_FILE;
		if (fqo->read_encoding == DEFAULT_ENCODING)
			fqo->read_encoding = IUPAC_ENCODING;
	}

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "file type = %s\n",
					file_type_name[fqd->file_type]);

	fqd->read_encoding = fqo->read_encoding == DEFAULT_ENCODING
		? XY_ENCODING : fqo->read_encoding;
	fqd->min_quality = MAX_ASCII_QUALITY_SCORE;
	fqd->max_quality = MIN_ASCII_QUALITY_SCORE;
	fqd->n_reads = 0;
	fqd->read_flag = NULL;
	fqd->site_flag = NULL;
	fqd->n_min_length = (unsigned int) -1;	/* maximum unsigned int */
	fqd->empty = 1;

	n_read_bytes = 0;	// bytes in current read
	n_bytes = 0;		// total bytes in reads
	n_name_bytes = 0;	// bytes in current name
	n_total_name_bytes = 0;	// total bytes in names
	n_last_read_bytes = 0;	// bytes in last read

	/* When this loop is complete, we will know the number of reads
	 * fqd->n_reads and the number of bytes n_bytes necessary to store them.
	 */
	do {
		err = read_read(fp, fqd, fqo, &n_read_bytes, NULL, NULL, NULL,
								&n_name_bytes);
		/* n_read_bytes: number of nucleotides that would be read
		 * n_name_bytes: number of bytes in name
		 */

		debug_msg(DEBUG_III <= fxn_debug, fxn_debug, "read read %u: "
			"err=%d (%s), length=%u, name length=%u (%u)\n",
			fqd->n_reads, err, fastq_error_message(err),
			n_read_bytes, n_name_bytes, n_bytes);

		/* valid read, increment read count and cumulative length */
		if (err == NO_ERROR || (err == FASTQ_EOF && n_name_bytes)) {

			n_bytes += n_read_bytes;
			if (n_read_bytes < fqd->n_min_length)
				fqd->n_min_length = n_read_bytes;
			fqd->n_reads++;
			n_total_name_bytes += n_name_bytes;
			debug_msg(DEBUG_III <= fxn_debug, fxn_debug, "read "
				"sequence %d of length %d (%u)\n", fqd->n_reads,
				n_read_bytes, n_bytes);

			/* reset last read length */
			n_last_read_bytes = 0;

		/* screen reads too short */
		} else if (err == FASTQ_READ_TOO_SHORT) {
			mmessage(WARNING_MSG, err, "%s in read %u: dropping "
						"read of length %u (%u)\n",
					fastq_error_message(err), fqd->n_reads,
							n_read_bytes, n_bytes);

			if (n_read_bytes > n_last_read_bytes)
				n_last_read_bytes = n_read_bytes;
			if (fqo->paired) {
				fqd->n_reads++;
				n_total_name_bytes += n_name_bytes;
			}

		/* screen reads too long */
		} else if (err == FASTQ_READ_TOO_LONG) {
			mmessage(WARNING_MSG, err, "%s in read %u: dropping "
						"read of length %u (%u)\n",
					fastq_error_message(err), fqd->n_reads,
							n_read_bytes, n_bytes);

			if (n_read_bytes > n_last_read_bytes)
				n_last_read_bytes = n_read_bytes;
			if (fqo->paired) {
				fqd->n_reads++;
				n_total_name_bytes += n_name_bytes;
			}

		/* tolerate ambiguous characters: drop these reads,
		 * but retain spot if paired reads
		 */
		} else if (err == FASTQ_AMBIGUOUS_READ_CHAR &&
			fqo->drop_ambiguous_reads) {

			mmessage(WARNING_MSG, err, "%s in read %u: dropping "
					"read (%u)\n", fastq_error_message(err),
							n_bytes, fqd->n_reads);
			if (n_read_bytes > n_last_read_bytes)
				n_last_read_bytes = n_read_bytes;
			if (fqo->paired) {
				fqd->n_reads++;
				n_total_name_bytes += n_name_bytes;
			}

		/* tolerate ambiguous characters: drop these nucleotides,
		 * or entire read if all ambiguous
		 */
		} else if (err == FASTQ_AMBIGUOUS_READ_CHAR &&
			fqo->drop_ambiguous_nucs) {

			mmessage(WARNING_MSG, err, "%s in read %u: dropping "
				"ambiguous bases, leaving %u bases (%u).\n",
				fastq_error_message(err), fqd->n_reads, 
							n_read_bytes, n_bytes);
			/* n_read_bytes includes ambiguous bases*/
			if (n_read_bytes) {
				++fqd->n_reads;
				n_bytes += n_read_bytes;
				n_total_name_bytes += n_name_bytes;
			} else if (fqo->paired) {
				fqd->n_reads++;
				n_total_name_bytes += n_name_bytes;
			}

			/* reset last read length */
			n_last_read_bytes = 0;

		/* tolerate ambiguous characters: change to iupac encoding */
		} else if (err == FASTQ_AMBIGUOUS_READ_CHAR
			&& fqo->read_encoding != XY_ENCODING) {

			if (fqd->read_encoding != IUPAC_ENCODING)
				mmessage(WARNING_MSG, err, "%s in read %u: "
					"using iupac encoding (%u)\n",
					fastq_error_message(err),
					fqd->n_reads + 1, n_bytes);
			fqd->read_encoding = IUPAC_ENCODING;
			fqd->n_reads++;
			n_bytes += n_read_bytes;
			n_total_name_bytes += n_name_bytes;

			/* reset last read length */
			n_last_read_bytes = 0;

		/* cannot tolerate ambiguous characters: abort */
		} else if (err == FASTQ_AMBIGUOUS_READ_CHAR) {
			mmessage(ERROR_MSG, err, "%s in read %u: cannot encode "
				"using XY encoding (%u)\n",
				fastq_error_message(err), fqd->n_reads + 1,
									n_bytes);
			if (n_read_bytes > n_last_read_bytes)
				n_last_read_bytes = n_read_bytes;
			break;

		/* tolerate invalid characters: drop read */
		} else if (err == FASTQ_INVALID_READ_CHAR
			|| err == FASTQ_INVALID_QUALITY_CHAR) {

			mmessage(WARNING_MSG, err, "%s: read %u will be "
				"discarded (%u)\n", fastq_error_message(err),
				fqd->n_reads + 1, n_bytes);
			if (n_read_bytes > n_last_read_bytes)
				n_last_read_bytes = n_read_bytes;

		/* other return codes indicate EOF or irrecoverable error */
		} else {
			break;
		}
	} while (err != FASTQ_EOF);

	if (err && err != FASTQ_EOF)
		return err;

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Found %u nucleotides in "
		"%u %s.\n", n_bytes, fqd->n_reads,
		fqd->file_type == FASTA_FILE ? "sequences" : "reads");

	/* add room to read in last read, even if will be discarded */
	/* also n_bytes includes dropped ambiguous nucleotides */
	n_bytes += n_last_read_bytes;
	fqd->reads = malloc(n_bytes * sizeof *fqd->reads);
	if (fqd->reads == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"fastq_data reads");

	fqd->quals = malloc(n_bytes * sizeof *fqd->quals);
	if (fqd->quals == NULL) {
		free(fqd->reads);
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"fastq_data quals");
	}

	fqd->n_lengths = malloc(fqd->n_reads * sizeof *fqd->n_lengths);
	if (fqd->n_lengths == NULL) {
		free(fqd->reads);
		free(fqd->quals);
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"fastq_data.n_lengths");
	}

	if (fqo->read_names) {
		/* extra bytes avoid seg fault if last read to be discarded */
		fqd->names = malloc((n_total_name_bytes
					+ (!n_read_bytes ? n_name_bytes : 0))
							* sizeof *fqd->names);
		fqd->name_lengths = malloc(fqd->n_reads * sizeof *fqd->name_lengths);
		if (!fqd->names || !fqd->name_lengths) {
			free(fqd->reads);
			free(fqd->quals);
			free(fqd->n_lengths);
			if (fqd->names)
				free(fqd->names);
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"fastq_data.names or fastq_data.name_lengths");
		}
		nptr = fqd->names;
		nlptr = fqd->name_lengths;
	}

	/* now we'll record the data */
	rewind(fp);

	rptr = fqd->reads;
	qptr = fqd->quals;
	uiptr = fqd->n_lengths;
	n_reads = 0;	/* double-check */
	fqd->n_max_length = 0;
	n_name_bytes = 0;
	n_bytes = 0;	/* debug */
	do {
		err = read_read(fp, fqd, fqo, &n_read_bytes, rptr, qptr, nptr, &n_name_bytes);

		/* process a read */
		if (err == NO_ERROR || (err == FASTQ_AMBIGUOUS_READ_CHAR
			&& fqo && fqo->drop_ambiguous_nucs)		/* drop nucs */
			|| (err == FASTQ_EOF && n_name_bytes)) {	/* last read */

			/* n_read_bytes excludes dropped ambiguous bases */
			if (n_read_bytes) {

				if (n_read_bytes > fqd->n_max_length)
					fqd->n_max_length = n_read_bytes;
				else if (n_read_bytes != fqd->n_min_length)
					elen = 0;	/* unequal length reads */
				n_reads++;
				if (fxn_debug > DEBUG_I) {
					unsigned char *str = display_sequence(rptr, n_read_bytes, fqd->read_encoding);
					fprintf(stderr, "%s:%d: reread sequence %d of length %d: %.*s (%d; %u)\n", __func__, __LINE__, n_reads, n_read_bytes, n_read_bytes, str, elen, n_bytes);
					free(str);
				}
				rptr += n_read_bytes;
				qptr += n_read_bytes;
				*uiptr = n_read_bytes;
				n_bytes += n_read_bytes;
				uiptr++;
				if (nptr) {
					*nlptr++ = n_name_bytes;
					nptr += n_name_bytes;
				}

			} else if (!n_read_bytes && fqo->paired) {

				/* placeholder only for proper pairing */
				n_reads++;
				*uiptr = n_read_bytes;
				uiptr++;
				if (nptr) {
					*nlptr++ = n_name_bytes;
					nptr += n_name_bytes;
				}
			}

		/* tolerable errors */
		} else if (err != FASTQ_INVALID_READ_CHAR
			&& err != FASTQ_AMBIGUOUS_READ_CHAR
			&& err != FASTQ_INVALID_QUALITY_CHAR
			&& err != FASTQ_READ_TOO_SHORT
			&& err != FASTQ_READ_TOO_LONG)
			break;
	} while (err != FASTQ_EOF);

	if (err && err != FASTQ_EOF)
		return err;

	if (n_reads != fqd->n_reads) {
		free(fqd->reads);
		free(fqd->quals);
		free(fqd->n_lengths);
		if (fqo->read_names) {
			free(fqd->names);
			free(fqd->name_lengths);
		}
		fprintf(stderr, "n_reads = %u; fastq::n_reads = %u\n", n_reads,
			fqd->n_reads);
		return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "fastq file");
	}

	if (fxn_debug >= DEBUG_III) {
		debug_msg(DEBUG_III <= fxn_debug, fxn_debug, "Reads 1-10:\n");
		rptr = fqd->reads;
		for (unsigned int i = 0; i < 10; ++i) {
			debug_msg(DEBUG_III <= fxn_debug, fxn_debug, "Read %2u:"
				" %.*s\n", i + 1, fqd->n_lengths[i],
				display_sequence(rptr, fqd->n_lengths[i],	/* [TODO, KSD] free memory */
				fqd->read_encoding));
			rptr += fqd->n_lengths[i];
		}
	}

	if (fqd->file_type == FASTQ_FILE) {
		debug_msg(SILENT <= fxn_debug, fxn_debug, "Minimum quality "
			"score: %c (%d)\n", fqd->min_quality,
			(int) fqd->min_quality);
		debug_msg(SILENT <= fxn_debug, fxn_debug, "Maximum quality "
			"score: %c (%d)\n", fqd->max_quality,
			(int) fqd->max_quality);
	}
	debug_msg(SILENT <= fxn_debug, fxn_debug, "Minimum %s length: %u\n",
			fqd->file_type==FASTA_FILE ? "sequence" : "read",
							fqd->n_min_length);
	debug_msg(SILENT <= fxn_debug, fxn_debug, "Maximum %s length: %u\n",
			fqd->file_type==FASTA_FILE ? "sequence" : "read",
							fqd->n_max_length);

	/* shift to 0-base qualities */
	qptr = fqd->quals;
	for (unsigned int i = 0; i < fqd->n_reads; ++i)
		for (unsigned int j = 0; j < fqd->n_lengths[i]; ++j) {
			*qptr = *qptr - fqd->min_quality;
			qptr++;
		}

	/* do not store lengths if not necessary */
	if (elen) {
		free(fqd->n_lengths);
		fqd->n_lengths = NULL;
	}

	/* fastq data object is now not empty! */
	fqd->empty = 0;

	return NO_ERROR;
} /* fread_fastq */

/**
 * Process a single read from a fastq file.
 *
 * @param fp	open fastq file handle
 * @param fqd	allocated fastq object
 * @param fqo	fastq options struct
 * @param len	pointer to memory to store length of current read
 * @param rptr	pointer to memory to store read base characters
 * @param qptr	pointer to memory to store quality characters
 * @param nptr	name pointer
 * @param nlen	name length
 *
 * @return	error code
 */
int read_read(FILE *fp, fastq_data *fqd, fastq_options *fqo, unsigned int *len,
	char_t *rptr, char_t *qptr, char *nptr, unsigned int *nlen)
{
	int fxn_debug = ABSOLUTE_SILENCE;//SILENT;//DEBUG_II;//DEBUG_I;//DEBUG_III;//
	int err = NO_ERROR;
	unsigned int qlen = 0;
	unsigned int *ambig_nuc = NULL, nambig_nuc = 0;
	int v_iupac, v_nuc;

	if (nlen)
		*nlen = 0;
	(*len) = 0;	/* signal incomplete read */

	/* check for @ at start of read record */
	char c = fgetc(fp);

	debug_msg(DEBUG_III <= fxn_debug, fxn_debug, "First char: '%c'\n", c);

	if (c != '@' && c != EOF && fqd->file_type == FASTQ_FILE)
		return FASTQ_FILE_FORMAT_ERROR;
	else if (c != '>' && c != EOF && fqd->file_type == FASTA_FILE)
		return FASTQ_FILE_FORMAT_ERROR;
	else if (c == EOF)
		return FASTQ_EOF;

	/* process read name */
	if (nptr)
		fforward_copy(fp, c, '\n', nptr, *nlen)
	else if (nlen)
		fforward_cnt(fp, c, '\n', *nlen)
	else
		fforward(fp, c, '\n');

	if (c == EOF)
		return FASTQ_PREMATURE_EOF;

	long int mark = ftell(fp);

	/* read or skip the nucleotide sequence */
	debug_msg(DEBUG_II <= fxn_debug, fxn_debug, "");
	do {
	while ((c = fgetc(fp)) != '\n' && c != EOF) {

		debug_msg_cont(DEBUG_II <= fxn_debug, fxn_debug, "%c", c);

		c = toupper(c);

		v_iupac = valid_iupac(&c);
		v_nuc = v_iupac ? valid_nucleotide(&c) : 0;

		/* record valid iupac */
		if (rptr && v_iupac && fqd->read_encoding == IUPAC_ENCODING) {
			*rptr = nuc_to_iupac[c - 'A'];
			rptr++;
			(*len)++;

		/* record valid nuc */
		} else if (rptr && v_nuc) {

			*rptr = (c >> 1) & 3;	/* unique encoding of A, C, G, T, but no others! */
			rptr++;
			(*len)++;

		/* invalid read: invalid character */
		} else if (!v_iupac) {

			err = FASTQ_INVALID_READ_CHAR;
			*len = 0;
			break;

		/* invalid read: cannot encode ambiguous nucleotide */
		} else if (!v_nuc && fqd->read_encoding == XY_ENCODING
			&& (!fqo || (fqo && !fqo->drop_ambiguous_nucs))) {

			err = FASTQ_AMBIGUOUS_READ_CHAR;
			*len = 0;
			break;

		/* ambiguous nucleotide that we'll drop */
		} else if (!v_nuc && fqd->read_encoding == XY_ENCODING) {

			err = FASTQ_AMBIGUOUS_READ_CHAR;
			if (qptr)	/* need to know full length */
				++(*len);

		/* ambiguous nucleotide we can encode */
		} else if (!v_nuc) {

			err = FASTQ_AMBIGUOUS_READ_CHAR;
			++(*len);

		/* valid nucleotide: count it */
		} else {

			++(*len);

		}

	}
	if (*len == 0)
		break;
	/* read through newline in middle of fasta record */
	if (fqd->file_type == FASTA_FILE && c == '\n') {
		/* peek at next char */
		c = fgetc(fp);
		if (c != EOF)
			ungetc(c, fp);
	}
	} while (fqd->file_type == FASTA_FILE && c != '>' && c != EOF);

	debug_msg_cont(DEBUG_II <= fxn_debug, fxn_debug, " (length = %u)\n", *len);

	/* invalid read: skip rest */
	if (err == FASTQ_INVALID_READ_CHAR) {

		mmessage(WARNING_MSG, FASTQ_INVALID_READ_CHAR, "%s '%c'\n",
			fastq_error_message(FASTQ_INVALID_READ_CHAR), c);
		fforward(fp, c, '\n');

	} else if (err == FASTQ_AMBIGUOUS_READ_CHAR
		&& fqo && fqo->drop_ambiguous_nucs) {
/*
		mmessage(WARNING_MSG, FASTQ_INVALID_READ_CHAR, "%s",
			fastq_error_message(FASTQ_AMBIGUOUS_READ_CHAR));
		fprintf(stderr, " will be dropped\n");
*/

		/* need to know which quality scores to skip */
		if (qptr) {

			fseek(fp, mark, SEEK_SET);
			ambig_nuc = calloc(*len, sizeof *ambig_nuc);
			if (!ambig_nuc)
				return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
					"fastq_data ambiguous nuc locations");
			qlen = 0;
			while ((c = fgetc(fp)) != '\n' && c != EOF) {
				v_nuc = valid_nucleotide(&c);
				if (!v_nuc)
					ambig_nuc[nambig_nuc++] = qlen;
				++qlen;
			}
			*len = *len - nambig_nuc;

		}

	/* ambiguous nucleotide: issue warning */
	} else if (err == FASTQ_AMBIGUOUS_READ_CHAR
		&& fqd->read_encoding == XY_ENCODING
		&& (!fqo || (fqo && !fqo->drop_ambiguous_nucs))) {

		mmessage(WARNING_MSG, FASTQ_AMBIGUOUS_READ_CHAR, "%s",
			fastq_error_message(FASTQ_AMBIGUOUS_READ_CHAR));

		fprintf(stderr, " '%c'", c);
		fforward(fp, c, '\n');
		fprintf(stderr, "\n");
	}

	if (fqo && fqo->max_length < UINT_MAX && *len > fqo->max_length)
		err = FASTQ_READ_TOO_LONG;
	if (fqo && fqo->min_length > 0 && *len < fqo->min_length)
		err = FASTQ_READ_TOO_SHORT;

	if (c == EOF) {
		if (fqd->file_type == FASTQ_FILE)
			return FASTQ_PREMATURE_EOF;
		return FASTQ_EOF;
	}

	if (fqd->file_type == FASTA_FILE)
		return NO_ERROR;

	/* fast forward through divider */
	fforward(fp, c, '\n');

	if (c == EOF)
		return FASTQ_PREMATURE_EOF;

	/* read or skip the quality score sequence */
	if (qptr && nambig_nuc) {
		nambig_nuc = 0;
		qlen = 0;
		while ((c = fgetc(fp)) != '\n' && c != EOF) {
			if (ambig_nuc[nambig_nuc] > qlen) {
				*qptr = (char_t) c;
				++qptr;
				++nambig_nuc;
			}
			qlen++;
		}
		qlen -= nambig_nuc;
	} else if (qptr) {
		while ((c = fgetc(fp)) != '\n' && c != EOF) {
			*qptr = (char_t) c;
			qptr++;
			qlen++;
			if (c < fqd->min_quality)
				fqd->min_quality = (char_t) c;
			if (c > fqd->max_quality)
				fqd->max_quality = (char_t) c;
		}
	} else {
		fforward_cnt(fp, c, '\n', qlen);
	}

	if (c == EOF) {
		/* incomplete quality score string */
		if (qlen < *len)
			return FASTQ_PREMATURE_EOF;
		return FASTQ_EOF;
	}

	return err;
} /* read_read */

int allocate_empty_fastq(fastq_data **in_fqd, fastq_options *fqo,
			unsigned int nreads, unsigned int read_length)
{
	int fxn_debug = ABSOLUTE_SILENCE;//SILENT;//DEBUG_III;//
	fastq_data *fqd = *in_fqd;

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "entering\n");

	if (fqd == NULL) {
		*in_fqd = malloc(sizeof **in_fqd);
		if (*in_fqd == NULL)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"fastq_data reads");
		fqd = *in_fqd;
		fqd->file_type = FASTQ_FILE;
		fqd->n_lengths = NULL;
		fqd->reads = NULL;
		fqd->quals = NULL;
		fqd->reference_seq = NULL;
		fqd->index = NULL;
	}

	fqd->empty = 1;
	fqd->read_encoding = fqo->read_encoding == DEFAULT_ENCODING
		? XY_ENCODING : fqo->read_encoding;
	fqd->min_quality = MIN_ASCII_QUALITY_SCORE;
	fqd->max_quality = MAX_ASCII_QUALITY_SCORE;
	fqd->n_reads = nreads;
	fqd->read_flag = NULL;
	fqd->site_flag = NULL;
	fqd->n_min_length = fqd->n_max_length = read_length;

	fqd->reads = malloc(nreads * read_length * sizeof *fqd->reads);
	if (fqd->reads == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"fastq_data reads");

	fqd->quals = malloc(nreads * read_length * sizeof *fqd->quals);
	if (fqd->quals == NULL) {
		free(fqd->reads);
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"fastq_data quals");
	}

	fqd->n_lengths = NULL;

	return NO_ERROR;
}/* allocate_empty_fastq */

/**
 * Allocate read flag vector.
 *
 * @param fqd	fastq data pointer with content
 * @return	error status
 */
int allocate_read_flag(fastq_data *fqd)
{
	if (!fqd->n_reads)
		return NO_ERROR;

	fqd->read_flag = malloc(fqd->n_reads * sizeof *fqd->read_flag);
	if (!fqd->read_flag)
		return MEMORY_ALLOCATION;
	else
		return NO_ERROR;
} /* allocate_read_flag */

/**
 * Allocate site flag vector.
 *
 * @param fqd	fastq data pointer with content
 * @return	error status
 */
int allocate_site_flag(fastq_data *fqd)
{
	if (!fqd->n_max_length)
		return NO_ERROR;
	
	fqd->site_flag = malloc(fqd->n_max_length * sizeof *fqd->site_flag);
	if (!fqd->site_flag)
		return MEMORY_ALLOCATION;
	else
		return NO_ERROR;
} /* allocate_site_flag */

void write_sequence_trimmed(FILE *fp, char_t const * const in_str,
	unsigned int len, int encoding, unsigned char *trim)
{
	for (unsigned int j = 0; j < len; ++j)
		if (!trim || !trim[j])
			fprintf(fp, "%c", encoding == IUPAC_ENCODING
				? iupac_to_char[(int) in_str[j]]
				: xy_to_char[(int) in_str[j]]);
} /* write_sequence_trimmed */

void write_sequence(FILE *fp, char_t const * const in_str, unsigned int len,
						int encoding)
{
	write_sequence_trimmed(fp, in_str, len, encoding, NULL);
} /* write_sequence */

unsigned char *display_sequence_trimmed(char_t const * const in_str,
	unsigned int len, int encoding, unsigned char *trim)
{
	unsigned int j, k = 0;
	unsigned char *str = malloc((len + 1) * sizeof *str);

	if (str == NULL) {
		mmessage(ERROR_MSG, MEMORY_ALLOCATION, "string to store read");
		return NULL;
	}

	for (j = 0; j < len; ++j)
		if (!trim || !trim[j])
			str[k++] = encoding == IUPAC_ENCODING
				? iupac_to_char[(int) in_str[j]]
				: xy_to_char[(int) in_str[j]];
	str[k] = '\0';

	return str;
} /* display_sequence_trimmed */

unsigned char *display_sequence(char_t const * const in_str,
	unsigned int len, int encoding)
{
	return display_sequence_trimmed(in_str, len, encoding, NULL);
} /* display_sequence */

void write_quals_trimmed(FILE *fp, char_t const * const in_str,
	unsigned int len, unsigned char min, unsigned char *trim)
{
	for (unsigned int j = 0; j < len; ++j)
		if (!trim || !trim[j])
			fprintf(fp, "%c", in_str[j] + min);
} /* write_quals_trimmed */

void write_quals(FILE *fp, char_t const * const in_str, unsigned int len, unsigned char min)
{
	write_quals_trimmed(fp, in_str, len, min, NULL);
} /* write_quals */

unsigned char *display_quals_trimmed(char_t const * const in_str,
	unsigned int len, char_t min, unsigned char *trim)
{
	unsigned int j, k = 0;
	unsigned char *str = malloc((len + 1) * sizeof *str);

	if (str == NULL) {
		mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"string to store qualities");
		return NULL;
	}

	for (j = 0; j < len; ++j)
		if (!trim || !trim[j])
			str[k++] = (unsigned char) (in_str[j] + min);
	str[k] = '\0';

	return str;
} /* display_quals_trimmed */

unsigned char *display_quals(char_t const * const in_str, unsigned int len, char_t min)
{
	return display_quals_trimmed(in_str, len, min, NULL);
} /* display_quals */

unsigned char *display_reverse_complement_trimmed(char_t const * const in_str,
		unsigned int len, int encoding, unsigned char *trim)
{
	unsigned int j, l, k = 0;
	unsigned char *str = malloc((len + 1) * sizeof *str);

	if (str == NULL) {
		mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"string to store reverse complement");
		return NULL;
	}

	for (j = len - 1, l = 0; l < len; --j, ++l)
		if (!trim || !trim[j])
			str[k++] = encoding == IUPAC_ENCODING
				? iupac_to_char[(int) iupac_to_rc[(int) in_str[j]]]
				: xy_to_char[(int) xy_to_rc[(int) in_str[j]]];
	str[k] = '\0';

	return str;
} /* display_reverse_complement_trimmed */

unsigned char *display_reverse_complement(char_t const * const in_str,
						unsigned int len, int encoding)
{
	return display_reverse_complement_trimmed(in_str, len, encoding, NULL);
} /* display_reverse_complement */

unsigned char *display_reverse_quals_trimmed(char_t const * const in_str,
	unsigned int len, char_t min, unsigned char *trim)
{
	unsigned int j, l, k = 0;
	unsigned char *str = malloc((len + 1) * sizeof *str);

	if (str == NULL) {
		mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"string to store qualities");
		return NULL;
	}

	for (j = len - 1, l = 0; l < len; --j, ++l)
		if (!trim || !trim[j])
			str[k++] = in_str[j] + min;
	str[k] = '\0';

	return str;
} /* display_reverse_quals_trimmed */

unsigned char *display_reverse_quals(char_t const * const in_str, unsigned int len, char_t min)
{
	return display_reverse_quals_trimmed(in_str, len, min, NULL);
} /* display_reverse_quals */

unsigned char *reverse_complement(char_t const *const in_str,
			unsigned int len, int encoding)
{
	unsigned char *str = malloc(len * sizeof *str);

	if (!str) {
		mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"string for reverse complement\n");
		return NULL;
	}
	for (unsigned int j = len - 1, l = 0; l < len; --j, ++l)
		str[l] = encoding == IUPAC_ENCODING
			? iupac_to_rc[(int) in_str[j]]
			: xy_to_rc[(int) in_str[j]];

	return str;
} /* reverse_complement */

void in_situ_reverse_complement(char_t * const str, unsigned int len, int encoding)
{
	char_t c;
	unsigned int halfway = len/2;
	for (unsigned int j = len - 1, l = 0; l < halfway; --j, ++l) {
		c = encoding == IUPAC_ENCODING
			? iupac_to_rc[(int) str[l]]
			: xy_to_rc[(int) str[l]];
		str[l] = encoding == IUPAC_ENCODING
			? iupac_to_rc[(int) str[j]]
			: xy_to_rc[(int) str[j]];
		str[j] = c;
	}
	if (len % 2)
		str[halfway] = encoding == IUPAC_ENCODING
			? iupac_to_rc[(int) str[halfway]]
			: xy_to_rc[(int) str[halfway]];
} /* in_situ_reverse_complement */

void in_situ_reverse(char_t * const str, unsigned int len)
{
	char_t c;
	unsigned int halfway = len / 2;
	for (unsigned int j = len - 1, l = 0; l < halfway; --j, ++l) {
		c = str[l];
		str[l] = str[j];
		str[j] = c;
	}
} /* in_situ_reverse */

int write_fastq(fastq_data *fqd, fastq_options *fqo)
{
	return write_fastq_marked_trimmed(fqd, fqo, NULL, 0, NULL, NULL);
} /* write_fastq */

/* abort
int write_fastq_selected_sorted_trimmed(fastq_data *fqd, fastq_options *fqo,
	unsigned int *id, unsigned int n_selected,
	int (*trim)(fastq_data *, char_t *, char_t *, unsigned int, unsigned int, void *), void *obj)
{
	unsigned int i, len;
	FILE *fp = fopen(fqo->outfile, fqo->append ? "a" : "w");

	if (!fp)
		return mmessage(ERROR_MSG, FILE_OPEN_ERROR, fqo->outfile);
	
	for (i = 0; i < n_selected; ++i) {
		len = read_length(fqd, i);
	}
} */

/**
 * Write fastq file with possible filtering or trimming of reads.
 * The mechanism for trimming is controlled by the caller via a callback
 * function, but the mechanism for filtering is by simple integer codes.
 * We can always make the filtering more complex via callback if we want.
 *
 * @param fqd		pointer to fastq_data object
 * @param fqo		pointer to fastq_options object
 * @param id		integer label for each read
 * @param selected_id	print reads with matching integer label
 * @param trim		trimming callback function
 * @param obj		void pointer to pass to callback function
 * @return		error status
 */
int write_fastq_marked_trimmed(fastq_data *fqd, fastq_options *fqo, unsigned int *id,
	unsigned int selected_id,
	int (*trim)(fastq_data *, char_t *, char_t *, unsigned int, unsigned int, void *), void *obj)
{
	unsigned int i, len;
	char_t *reads = fqd->reads;
	char_t *quals = fqd->quals;
	char *names = fqd->names;
	char_t *str = NULL;
	FILE *fp = fopen(fqo->outfile, fqo->append ? "a" : "w");

	if (!fp)
		return(mmessage(ERROR_MSG, FILE_OPEN_ERROR, fqo->outfile));

	for (i = 0; i < fqd->n_reads; ++i) {
		len = read_length(fqd, i);
		if (!id || id[i] == selected_id) {
			if (names)
				fprintf(fp, "%c%.*s", fqo->fasta ? '>' : '@',
						fqd->name_lengths[i], names);
			else
				fprintf(fp, "%c%u", fqo->fasta ? '>' : '@', i);
			if (fqo->casavize)
				fprintf(fp, ":1:1:1:1:0:0 %u:N:0:A",
						(unsigned int) fqo->casavize);
			fprintf(fp, "\n");
			if (trim && trim(fqd, reads, quals, i, len, obj))
				return(mmessage(ERROR_MSG, INTERNAL_ERROR,
							"trim function\n"));
			str = !fqo->reverse_complement
				? display_sequence_trimmed(reads, len,
					fqd->read_encoding, fqd->site_flag)
				: display_reverse_complement_trimmed(reads,
					len, fqd->read_encoding, fqd->site_flag);
			fprintf(fp, "%s\n", str);
			free(str);
			if (!fqo->fasta) {
				fprintf(fp, "+\n");
				str = !fqo->reverse_complement
					? display_quals_trimmed(quals, len,
						fqd->min_quality, fqd->site_flag)
					: display_reverse_quals_trimmed(quals,
						len, fqd->min_quality,
								fqd->site_flag);
				fprintf(fp, "%s\n", str);
				free(str);
			}
		}
		if (names)
			names += fqd->name_lengths[i];
		reads += len;
		if (!fqo->fasta)
			quals += len;
	}

	fclose(fp);

	return NO_ERROR;
} /* write_fastq_marked_trimmed */

/**
 * Write fastq data as R-style data table.
 *
 * @param fqd		fastq_data object pointer
 * @param filename	name of output file
 * @return		error status
 */
int write_table(fastq_data *fqd, char const *filename)
{
	return write_table_marked(fqd, filename, NULL, 0);
} /* write_table */

/**
 * Write fastq data as R-style data table.
 *
 * @param fqd		fastq_data object pointer
 * @param filename	name of output file
 * @param mark		integer mark for each read
 * @param val		value of mark to select read for printing
 * @return		error status
 */
int write_table_marked(fastq_data *fqd, char const *filename,
	unsigned int *mark, unsigned int val)
{
	unsigned int i, len;
	char_t *reads = fqd->reads;
	FILE *fp = fopen(filename, "w");

	if (!fp)
		return(mmessage(ERROR_MSG, FILE_OPEN_ERROR, filename));
	
	for (i = 0; i < fqd->n_reads; ++i) {
		if (!mark || mark[i] == val) {
			len = read_length(fqd, i);
			write_read_in_table(fp, reads, len);
			reads += len;
		}
	}
	fclose(fp);

	return NO_ERROR;
} /* write_table_marked */

/**
 * Free fastq object.
 *
 * @param fqd	fastq object pointer
 */
void free_fastq(fastq_data *fqd)
{
	if (fqd) {
		if (fqd->reads)
			free(fqd->reads);
		if (fqd->quals)
			free(fqd->quals);
		if (fqd->n_lengths)
			free(fqd->n_lengths);
		if (fqd->index)
			free(fqd->index);
		if (fqd->name_lengths)
			free(fqd->name_lengths);
		if (fqd->names)
			free(fqd->names);
		if (fqd->read_flag)
			free(fqd->read_flag);
		if (fqd->site_flag)
			free(fqd->site_flag);
		free(fqd);
	}
} /* free_fastq */

/**
 * Return human-friendly error message for fastq error code.
 *
 * @param err_no	error number
 * @return		error message string
 */
const char *fastq_error_message(int err_no)
{
	if (err_no == FASTQ_INVALID_READ_CHAR)
		return "illegal nucleotide";
	else if (err_no == FASTQ_AMBIGUOUS_READ_CHAR)
		return "ambiguous nucleotide";
	else if (err_no == FASTQ_INVALID_QUALITY_CHAR)
		return "illegal quality score";
	else if (err_no == FASTQ_INCOMPLETE_READ)
		return "incomplete read";
	else if (err_no == FASTQ_FILE_FORMAT_ERROR)
		return "invalid file format";
	else if (err_no == FASTQ_PREMATURE_EOF)
		return "premature end-of-file";
	else if (err_no == FASTQ_EOF)
		return "end-of-file";
	else if (err_no == FASTQ_READ_TOO_SHORT)
		return "read too short";
	else if (err_no == FASTQ_READ_TOO_LONG)
		return "read too long";
	else
		return "No error";
} /* fastq_error_message */

/**
 * Check if two reads are equal.  This function assumes there is biological
 * homology at read positions.
 *
 * This may be most useful outside fastq_data where the homology
 * relationship among reads is known.
 *
 * @param fqd	fastq object pointer
 * @param i	index of first read
 * @param j	index of second read
 * @return	<0, 0, >0
 */
int read_compare(fastq_data *fqd, unsigned int i, unsigned int j)
{
	/* same index */
	if (i == j)
		return 0;

	unsigned int len = fqd->n_max_length;
	size_t i_idx = len * i, j_idx = len * j;

	if (fqd->n_lengths) {

		/* assume reads of different lengths not equal */
		if (fqd->n_lengths[i] < fqd->n_lengths[j])
			return -1;
		else if (fqd->n_lengths[i] > fqd->n_lengths[j])
			return 1;

		/* they share a common length */
		len = fqd->n_lengths[i];

		/* find index of ith and jth read */
		i_idx = j_idx = 0;
		for (unsigned int l = 0; l < MIN(i, j); ++l) {
			i_idx += fqd->n_lengths[l];
			j_idx += fqd->n_lengths[l];
		}
		if (i < j)
			for (unsigned int l = i; l < j; ++l)
				j_idx += fqd->n_lengths[l];
		else
			for (unsigned int l = j; l < i; ++l)
				i_idx += fqd->n_lengths[l];
	}

	for (unsigned int l = 0; l < len; ++l) {
		if (fqd->reads[i_idx + l] < fqd->reads[j_idx + l])
			return -1;
		else if (fqd->reads[i_idx + l] > fqd->reads[j_idx + l])
			return 1;
	}
	return 0;
} /* read_compare */

#ifndef NO_ALIGNMENT
int pw_align_reads(fastq_data *fqd, char const * const rfile)
{
	int err = NO_ERROR;
	FILE *fp = fopen(rfile, "r");	/* fasta format */

	if (!fp)
		return mmessage(ERROR_MSG, FILE_OPEN_ERROR, rfile);

	char c = fgetc(fp);
	if (c != '>') {
		fclose(fp);
		return FILE_FORMAT_ERROR;
	}

	fforward(fp, c, '\n');

	unsigned int len = 0;
	while ((c = fgetc(fp)) != EOF) {
		if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
			len++;
		else if (c != '\n' && c != ' ' && c != '\t') {
			fclose(fp);
			return mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"Reference sequence contains ambiguous "
				"character: %c\n", c);
		}
	}

	fqd->reference_seq = malloc(len * sizeof *fqd->reference_seq);

	if (!fqd->reference_seq) {
		fclose(fp);
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"fastq_data.reference_seq");
	}

	rewind(fp);
	fforward(fp, c, '\n');

	for (unsigned int i = 0; i < len; ++i) {
		c = fgetc(fp);
		if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
			fqd->reference_seq[i] = (c >> 1) & 3;
		else
			--i;
	}

	fclose(fp);

	mmessage(INFO_MSG, NO_ERROR, "Reference sequence: ");
	for (unsigned int i = 0; i < len; ++i)
		fprintf(stderr, "%c", xy_to_char[(int) fqd->reference_seq[i]]);
	fprintf(stderr, "\n");

	char_t *rptr = fqd->reads;
	int score[NUM_NUCLEOTIDES][NUM_NUCLEOTIDES] = {{2, -3, -3, -2},
		{-3, 2, -2, -3}, {-3, -2, 2, -3}, {-2, -3, -3, 2}};
	//double const perr[] = {0.999401, 0.992814, 0.993413, 0.997006, 0.996407, 0.994910, 0.994311, 0.742627, 0.993713, 0.995808, 0.992216, 0.997305, 0.997006, 0.997904, 0.968593, 0.993114, 0.994611, 0.995808, 0.997305, 0.998204, 0.975150, 0.998503, 0.998204, 0.994611, 0.984731, 0.967365, 0.996707, 0.997904, 0.915868, 0.984431, 0.987126, 0.997605, 0.979042, 0.993114, 0.994311, 0.989820, 0.985629, 0.993114, 0.985329, 0.977246, 0.995808, 0.997305, 0.986527, 0.996108, 0.997006, 0.988623, 0.989532, 0.970659, 0.944346, 0.998204, 0.989820, 0.996707, 0.996707, 0.960778, 0.982934, 0.986527, 0.998204, 0.986527, 0.998503, 0.981737, 0.991916, 0.991018, 0.995210, 0.985329, 0.991916, 0.978789, 0.972156, 0.970060, 0.994963, 0.990719, 0.992814, 0.994611, 0.991916, 0.985917, 0.976700, 0.990419, 0.996707, 0.991018, 0.982635, 0.985329, 0.997305, 0.986228, 0.978144, 0.997006, 0.994311, 0.994012, 0.996108, 0.985928, 0.971856, 0.962874, 0.980838, 0.986228, 0.988323, 0.994311, 0.915194, 0.971512, 0.994311, 0.968862, 0.977545, 0.981437, 0.985329, 0.997006, 0.995210, 0.994012, 0.993450, 0.987183, 0.993413, 0.991916, 0.996379, 0.998503, 0.985629, 0.993513, 0.998503, 0.987725, 0.975449, 0.981138, 0.979641, 0.961770, 0.996108, 0.994910, 0.981437, 0.981437, 0.985329, 0.994311, 0.965015, 0.990075, 0.953293, 0.979341, 0.990120, 0.991617, 0.993114, 0.995210, 0.985928, 0.996707, 0.997904, 0.985329, 0.998204, 0.979341, 0.991916, 0.979641, 0.981437, 0.984132, 0.991617, 0.997006, 0.959880, 0.991916, 0.996707, 0.989521, 0.997006, 0.961386, 0.993413, 0.972156, 0.995509, 0.973653, 0.990120, 0.996707, 0.976647, 0.988323, 0.997605, 0.988922, 0.985629, 0.965269, 0.998503, 0.977545, 0.971557, 0.973952, 0.981737, 0.992814, 0.986527, 0.981737, 0.995509, 0.994311, 0.981737, 0.997006, 0.973054, 0.965569, 0.994311, 0.978144, 0.972455, 0.994311, 0.990719, 0.988623, 0.997305, 0.964072, 0.988623, 0.986527, 0.991018, 0.995210, 0.996407, 0.981437, 0.971788, 0.959581, 0.982335, 0.992515, 0.993713, 0.991617, 0.993114, 0.988623, 0.986826, 0.994611, 0.987126, 0.967504, 0.997305, 0.936483, 0.994910, 0.996707, 0.982335, 0.996407, 0.997305, 0.947006, 0.985940, 0.994910, 0.963473, 0.950299, 0.953593, 0.994311, 0.972156, 0.995210, 0.989222, 0.920475, 0.997305, 0.941018, 0.988024, 0.971557, 0.960479, 0.989222, 0.994910, 0.988323, 0.977246, 0.996707, 0.967365, 0.912167, 0.948165, 0.995509, 0.979940, 0.985050, 0.944311, 0.973353, 0.947305, 0.990719, 0.987126, 0.970958, 0.975150, 0.997006, 0.992814, 0.932335, 0.948802, 0.933832, 0.937504, 0.991317, 0.982335, 0.991617, 0.967365, 0.956287, 0.996108, 0.960180, 0.985329, 0.994311, 0.971257, 0.994611, 0.994311, 0.937126, 0.994311, 0.995210, 0.997904, 0.997305, 0.904790, 0.923653, 0.926048, 0.901982, 0.998503, 0.986826, 0.964970, 0.997605};

	for (unsigned int i = 0; i < fqd->n_reads; ++i) {
		unsigned int alen;
		unsigned char **aln = nw_alignment(fqd->reference_seq, rptr,
			len, read_length(fqd, i), score, -1, -1, 1,
			NULL, NULL, NULL, NULL, &err, &alen, NULL, 0,
							fqd->read_encoding);
		fprintf(stderr, "Read %u alignment length %u\n", i, alen);
		size_t ngap1 = 0, ngap2 = 0;
		for (size_t j = 0; j < alen; ++j) {
			if (aln[0][j] == '-') ngap1++;
			fprintf(stderr, "%c", aln[0][j] == '-'
				? '-' : xy_to_char[(int) aln[0][j]]);
		}
		fprintf(stderr, "\n");
		for (size_t j = 0; j < alen; ++j) {
			if (aln[1][j] == '-') ngap2++;
			fprintf(stderr, "%c", aln[1][j] == '-'
				? '-' : xy_to_char[(int) aln[1][j]]);
		}
		fprintf(stderr, "\nGaps: %lu %lu\n", ngap1, ngap2);
		if (ngap1 != ngap2) exit(0);
		rptr += read_length(fqd, i);
	}

	return NO_ERROR;
} /* pw_align_reads */
#endif

double read_distance(fastq_data *fqd, unsigned int i, unsigned int j)
{
	if (i == j)
		return 0;

	char_t *qptr1 = NULL, *qptr2 = NULL;
	char_t *rptr1 = NULL, *rptr2 = NULL;
	char_t *rptr = fqd->reads;
	char_t *qptr = fqd->quals;
	unsigned int len;

	for (unsigned int n = 0; n < fqd->n_reads; ++n) {
		if (i == n) {
			rptr1 = rptr;
			qptr1 = qptr;
		}
		if (j == n) {
			rptr2 = rptr;
			qptr2 = qptr;
		}
		if ((i < j && j == n) || (j < i && i == n))
			break;

		len = read_length(fqd, n);
		rptr += len;
		qptr += len;
	}

	len = fqd->n_lengths ? MIN(fqd->n_lengths[i], fqd->n_lengths[j])
		: fqd->n_max_length;
	return read_distance_ptr(fqd, len, rptr1, rptr2, qptr1, qptr2);
} /* read_distance */

double read_distance_ptr(fastq_data *fqd, unsigned int len, char_t *rptr1,
	char_t *rptr2, char_t *qptr1, char_t *qptr2)
{
	if (qptr1 == qptr2) {
fprintf(stderr, "here!\n");
		return 0;
	}
	double dis = 0;
	for (unsigned int n = 0; n < len; ++n) {
/*
		fprintf(stderr, "%c %c %f %f", fqd->read_encoding == IUPAC_ENCODING
                        ? iupac_to_char[(int) *rptr1]
                        : xy_to_char[(int) *rptr1],
			fqd->read_encoding == IUPAC_ENCODING
                        ? iupac_to_char[(int) *rptr2]
                        : xy_to_char[(int) *rptr2],
			error_prob(fqd, *qptr1), error_prob(fqd, *qptr2));
*/
		double ldis = 0;
		if (*rptr1 == *rptr2) {
			ldis += error_prob(fqd, *qptr1) * (1 - error_prob(fqd, *qptr2));
			ldis += error_prob(fqd, *qptr2) * (1 - error_prob(fqd, *qptr1));
			ldis += 2*error_prob(fqd, *qptr2) * error_prob(fqd, *qptr1)/3;
		} else {
			ldis += error_prob(fqd, *qptr1)/3 * (1 - error_prob(fqd, *qptr2));
			ldis += error_prob(fqd, *qptr2)/3 * (1 - error_prob(fqd, *qptr1));
			ldis += 2*error_prob(fqd, *qptr2) * error_prob(fqd, *qptr1) / 9;
		}
		dis += ldis;
//fprintf(stderr, ": %f (%f)\n", ldis, dis);
		qptr1++;
		qptr2++;
		rptr1++;
		rptr2++;
	}
	return dis;
} /* read_distance_ptr */

char const *read_name(fastq_data *fqd, unsigned int j)
{
	char *names = fqd->names;

	for (unsigned int i = 0; i < j; ++i)
		names += fqd->name_lengths[i];
	
	return names;
}/* read_name */

char const *next_read_name(fastq_data *fqd, char const **names, unsigned int j)
{
	char const *name = *names;

	*names += fqd->name_lengths[j];
	return name;
} /* next_read_name */

char *next_read_name_rw(fastq_data *fqd, char **names, unsigned int j)
{
	char *name = *names;

	*names += fqd->name_lengths[j];

	return name;
}/* next_read_name_rw */

int make_fastq_options(fastq_options **opt)
{
	fastq_options *op;
	*opt = malloc(sizeof **opt);
	if (!*opt)
		return(mmessage(ERROR_MSG, MEMORY_ALLOCATION, "fastq_data"));
	op = *opt;

	op->outfile = NULL;
	op->read_encoding = DEFAULT_ENCODING;
	op->reverse_complement = 0;
	op->fasta = 0;
	op->drop_ambiguous_reads = 0;
	op->drop_ambiguous_nucs = 0;
	op->read_names = 0;
	op->max_length = UINT_MAX;
	op->min_length = 0;
	op->paired = 0;
	op->casavize = 0;

	return NO_ERROR;
} /* make_fastq_options */

void free_fastq_options(fastq_options *opt)
{
	if (opt)
		free(opt);
} /* free_fastq_options */
