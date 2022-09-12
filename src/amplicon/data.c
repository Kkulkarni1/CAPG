/**
 * @file data.c
 * @author Karin S. Dorman
 *
 * Manipulate data structure for ampliclust.
 *
 * Note about formatting.  Line widths are at 80 characters, not because we live
 * in the 60's but to help force good coding and to reduce complexity.  Function
 * predeclarations may break this rule so that the entire prototype can be
 * found with a simple grep on the source code.
 */

#include <stdlib.h>
#include <string.h>

#include "data.h"
#include "options.h"
#include "model.h"
#include "simulate.h"
#include "util.h"
#include "math.h"
#include "error.h"
#include "hash.h"

int build_hash(data *dat);

/**
 * Setup data object.
 */
int make_data(data **dat, options *opt)
{
	UNUSED(opt);
	data *dp;
	*dat = malloc(sizeof **dat);

	if (*dat == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data object");

	dp = *dat;

	dp->fdata = NULL;
	dp->dmat = NULL;
	dp->qmat = NULL;

	dp->read_idx = NULL;
	dp->sample_size = 0;
	dp->max_read_length = 0;
	dp->min_read_length = 0;
	dp->lengths = NULL;

	dp->offset_ori = NULL;
	dp->offset = NULL;
	dp->max_offset = 0;
	dp->max_offset_ori = 0;
	dp->max_read_position = 0;
	dp->coverage = NULL;

	dp->seq_count = NULL;
	dp->hash_length = 0;

	return NO_ERROR;
} /* make_data */

/**
 * Load data into data object.  Once fastq file is read, we load
 * structural information into the data object.  Nothing we do now,
 * however, uses the data content (nucleotides and quality scores).
 * Use sync_data() to finalize data-dependent information.
 *
 * @param dat	pointer to data object
 * @param opt	pointer to options object
 * @param tbd	read nucleotides to be determined (simulation)
 * @return	error status
 */
int load_data(data *dat, options *opt, int tbd)
{
	int fxn_debug = ABSOLUTE_SILENCE;
	int err = NO_ERROR;

	if (!dat->fdata)
		return mmessage(ERROR_MSG, INVALID_USER_INPUT, "There is no "
			"data. See command-line flags: -s -f -m\n");

	/* currently, the sample is the complete data */
	dat->sample_size = dat->fdata->n_reads;

	if (opt->multi_stage_method && opt->sample_size < dat->fdata->n_reads)
		dat->sample_size = opt->sample_size;

	dat->max_read_length = dat->fdata->n_max_length;
	dat->min_read_length = dat->fdata->n_min_length;
	dat->n_quality = dat->fdata->max_quality - dat->fdata->min_quality + 1;
	dat->lengths = malloc(dat->sample_size * sizeof *dat->lengths);

	if (!dat->lengths)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data.lengths");

	/* [BUG?] first sample is first data::sample_size reads, assumption
	 * continues throughout this function with consequences
	 */
	if (dat->fdata->n_lengths)
		memcpy(dat->lengths, dat->fdata->n_lengths, dat->sample_size
			* sizeof *dat->lengths);
	else
		for (size_t i = 0; i < dat->sample_size; ++i)
			dat->lengths[i] = dat->max_read_length;

	/* allocate the index array of reads */
	dat->read_idx = malloc(dat->fdata->n_reads * sizeof *dat->read_idx);

	if (!dat->read_idx)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "dat.read_idx");

	for (size_t i = 0; i < dat->fdata->n_reads; ++i)
		dat->read_idx[i] = i;

	/* allocate short-cut pointers to nucleotide sequences */
	dat->dmat = malloc(dat->sample_size * sizeof *dat->dmat);

	if (!dat->dmat)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "dat.dmat");

	unsigned char *rptr = dat->fdata->reads;
	for (size_t i = 0; i < dat->sample_size; ++i) {
		dat->dmat[i] = rptr;
		rptr += dat->lengths[i];
	}

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Allocated %dx(%d) sequence "
		"matrix\n", dat->sample_size, dat->max_read_length);

	/* allocate short-cut pointers to quality sequences  */
	dat->qmat = malloc(dat->sample_size * sizeof *dat->qmat);

	if (!dat->qmat)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "dat.qmat");

	unsigned char *qptr = dat->fdata->quals;
	for (size_t i = 0; i < dat->sample_size; i++) {
		dat->qmat[i] = qptr;
		qptr += dat->lengths[i];
	}

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Allocated %dx(%d) quality "
		"matrix\n", dat->sample_size, dat->max_read_length);

	dat->max_read_position = dat->max_read_length;
	if (opt->offset_file) {
		dat->offset = malloc(dat->fdata->n_reads * sizeof *dat->offset);
		if (!dat->offset)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"data.offsets");
		if ((err = read_uints(opt->offset_file, dat->offset,
			dat->fdata->n_reads)))
			return err;
		for (size_t i = 0; i < dat->fdata->n_reads; ++i) {
			if (dat->offset[i] > dat->max_offset)
				dat->max_offset = dat->offset[i];
			if (dat->lengths[i] + dat->offset[i]
				> dat->max_read_position)
				dat->max_read_position = dat->lengths[i]
					+ dat->offset[i];
		}

	}
	dat->max_offset_ori = dat->max_offset;
	if (opt->multi_stage_method && dat->offset) {
		dat->offset_ori = malloc(dat->fdata->n_reads
						* sizeof *dat->offset_ori);
		if (!dat->offset_ori)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"data.offsets_ori");
		memcpy(dat->offset_ori, dat->offset,
				dat->fdata->n_reads * sizeof *dat->offset_ori);
	}

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "dat->max_offset:%i\n",
							dat->max_offset);

	dat->coverage = malloc(dat->max_read_position * sizeof *dat->coverage);

	if (!dat->coverage)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data.coverage");

	/* [TODO] If read lengths are constant and no offset, there is a much
	 * quicker way to set data::coverage.
	 */
	for (size_t j = 0; j < dat->max_read_position; ++j)
		dat->coverage[j] = 0;

	for (size_t i = 0; i < dat->sample_size; ++i)
		for (size_t j = 0; j < dat->lengths[i]; ++j)
			dat->coverage[j + (dat->offset ? dat->offset[i] : 0)]++;

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "coverage");

	/* check ancestor length matches read length, and if so, 
	 * convert to XY_ENCODING (comes in as char)
	 */
	if (opt->simo->ancestor) {
		if (strlen((char *)opt->simo->ancestor) != dat->max_read_length)
			return mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"Ancestor is not same length as data.\n");
		for (unsigned int j = 0; j < dat->max_read_length; ++j)
			opt->simo->ancestor[j] = (opt->simo->ancestor[j] >> 1) & 3L;
	}

	if (!dat->fdata->empty && !tbd)
		err = sync_data(dat);

	return err;
} /* load_data */


/**
 * Build hash of reads.
 *
 * @param dat	pointer to data object
 * @return	error status
 */
int build_hash(data *dat)
{
	int err = NO_ERROR;
	/* build hash table */
	dat->hash_length = 0;
	for (size_t i = 0; i < dat->sample_size; ++i) {
//fprintf(stderr, "%zu: %s\n", i, display_sequence(dat->dmat[i], dat->lengths[i], dat->fdata->read_encoding));
		dat->hash_length += add_sequence(&dat->seq_count, dat->dmat[i],
							dat->lengths[i], i);
	}

	/* store index of reads for all unique sequences */
	for (size_t i = 0; i < dat->sample_size; ++i)
		if ((err = add_seq_idx(dat->seq_count, dat->dmat[i],
							dat->lengths[i], i)))
			return err;
	
	/* sort hash table by count */
	sort_by_count(&dat->seq_count);

	return err;
} /* build_hash */


/**
 * [KSD: stays here because uses data_t data format, and fastq.[ch] unaware]
 * [BUG: this assumes xy_t type, but could be iupac_t; need to fix data_t issue]
 */
void fprint_fasta(FILE *fp, data_t *data, size_t n, size_t p, char const * const prefix) {
	for (size_t i = 0; i < n; ++i) {
		fprintf(fp, ">%s%lu\n", prefix, i);
		for (size_t j = 0; j < p; ++j)
			fprintf(fp, "%c", xy_to_char[(int)data[i*p + j]]);
		fprintf(fp, "\n");
	}
} /* fprint_fasta */

/**
 * [KSD: see fprint_fasta()]
 */
#ifdef USE_CURSES
void wprint_fasta(WINDOW *wp, data_t *data, size_t n, size_t p, char const * const prefix) {
	for (size_t i = 0; i < n; ++i) {
		wprintw(wp, ">%s%u\n", prefix, i);
		for (size_t j = 0; j < p; ++j)
			wprintw(wp, "%c", xy_to_char[(int)data[i*p + j]]);
		wprintw(wp, "\n");
	}
} /* wprint_fasta */
#endif

/**
 * [KSD: stays here because uses data_t data format, which fastq.[ch] do not know]
 */
void fprint_alignment2(FILE *fp, data_t **data, size_t n, size_t p) {
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < p; ++j)
			fprintf(fp, "%c", xy_to_char[(int) data[i][j]]);
		fprintf(fp, "\n");
	}
} /* fprint_alignment2 */

/**
 * [KSD: see fprint_alignment2()]
 */
#ifdef USE_CURSES
void wprint_alignment2(WINDOW *wp, data_t **data, size_t n, size_t p) {
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < p; ++j)
			wprintw(wp, "%c", xy_to_char[(int) data[i][j]]);
		wprintw(wp, "\n");
	}
} /* wprint_alignment2 */
#endif

/**
 * Check if two reads are equal where they overlap.
 *
 * @param dat	data object pointer
 * @param i	index of first read
 * @param j	index of second read
 * @return	<0, 0, >0
 */
int amplicon_read_compare(data *dat, size_t i, size_t j)
{
	/* same index */
	if (i == j)
		return 0;

	unsigned int len = MIN(dat->lengths[i], dat->lengths[j]);

	for (size_t l = 0; l < len; ++l) {

		if (dat->dmat[i][l] < dat->dmat[j][l])
			return -1;
		else if (dat->dmat[i][l] > dat->dmat[j][l])
			return 1;
	}

	return 0;

} /* amplicon_read_compare */

/**
 * Free data object.
 */
void free_data(data *dat)
{
	if (dat) {
		if (dat->fdata)
			free_fastq(dat->fdata);
		if (dat->dmat)
			free(dat->dmat);
		if (dat->qmat)
			free(dat->qmat);
		if (dat->read_idx)
			free(dat->read_idx);
		if (dat->lengths)
			free(dat->lengths);
		if (dat->offset)
			free (dat->offset);
		if (dat->offset_ori)
			free(dat->offset_ori);
		if (dat->coverage)
			free (dat->coverage);
		if (dat->seq_count)
			delete_all(dat->seq_count);
		free(dat);
	}
} /* free_data */

/**
 * Update the data object at each stage, preparing to run AmpliclustK on the
 * current subsample of the data.
 * Based on sample_idx
 * including dmat, qmat,lengths, offset, sample_size, coverage
 * max_read_length, max_read_position, max_offset, hash table
 * 
 * @param dat		pointer to data object
 * @param opt		pointer to options object
 * @param sample_size	size of sample we are preparing
 * @param sample_idx    the index of sample we are preparing 
 * @return		error status
 */
int update_data(data *dat, options *opt, size_t sample_size,size_t *sample_idx)
{

	int fxn_debug = ABSOLUTE_SILENCE;	//DEBUG_I;	//
	int err = NO_ERROR;
	size_t pre_sample_size = dat->sample_size;
	dat->sample_size = sample_size;

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "sample_size '%i'\n",
							 dat->sample_size);

	/* if sample size bigger than planned, we need realloc 
	 * dmat, qmat, lengths, offset, coverage */
	if (dat->sample_size > pre_sample_size ) {
		/* dmat, qmat */
		data_t **dmat = realloc(dat->dmat,
			dat->sample_size * sizeof *dat->dmat);
		data_t **qmat = realloc(dat->qmat,
			dat->sample_size * sizeof *dat->qmat);

		if (!dmat || !qmat )
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"realloc.data.dqmat");

		dat->dmat = dmat;
		dat->qmat = qmat;
		/*
		free(dat->dmat);
		free(dat->qmat);
		dat->dmat = NULL;
		dat->qmat = NULL;
		dat->dmat = malloc (dat->sample_size * sizeof *dat->dmat);
		dat->qmat = malloc (dat->sample_size * sizeof *dat->qmat);
		*/

		/* lengths, coverage */
		unsigned int *lengths = realloc(dat->lengths,
			sample_size * sizeof *dat->lengths);
		size_t *coverage = realloc(dat->coverage,
			dat->max_read_position * sizeof *dat->coverage);

		if (!lengths || !coverage)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"realloc.data.length.coverage");

		dat->lengths = lengths;
		dat->coverage = coverage;

		/* offset */
		if (opt->offset_file) {  /* dat->offset = NULL if there is no input of offset */
			unsigned int *offset = realloc(dat->offset,
				sample_size * sizeof *dat->offset);
			if (!offset)
				return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
					"realloc.data.offset");
			dat->offset = offset;
		}
	}

	 /* reallocate short-cut pointers to nucleotide and quality sequences */
	 /* [KSD, TODO] Aren't the first opt->sample_size pointers correct still? */
	 /*  opt->sample_size is constant while dat->sample_size is changing for different tasks */
	 
	 /* dmat, qmat */
	 unsigned char *rptr = dat->fdata->reads;
	 unsigned char *qptr = dat->fdata->quals;
	 for (size_t i = 0; i < dat->fdata->n_reads; ++i) {
		 for (size_t j = 0; j < dat->sample_size; ++j) {
			 if (sample_idx[j] == i) {
				 dat->dmat[j] = rptr;
				 dat->qmat[j] = qptr;
				 break;
			 }
		 }
		 rptr += read_length(dat->fdata, i);
		 qptr += read_length(dat->fdata, i);
	 }

	 debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Allocated %dx(%d) sequence"
		" matrix\n", dat->sample_size, dat->fdata->n_max_length);

	 unsigned int max_offset = 0;
	 unsigned int max_read_length = 0;
	 unsigned int max_read_position = 0;
	 for (size_t j = 0; j < dat->sample_size; ++j) {
		dat->lengths[j] = read_length(dat->fdata, sample_idx[j]);

		if (dat->lengths[j] > max_read_length)
			max_read_length = dat->lengths[j];

		if (dat->offset_ori) {
			dat->offset[j] = dat->offset_ori[sample_idx[j]];
			if (dat->offset[j] > max_offset)
				max_offset = dat->offset[j];
		}
		if (dat->lengths[j] + ( dat->offset ? dat->offset[j] : 0 )
			> max_read_position)
			max_read_position = dat->lengths[j]
				+ ( dat->offset ? dat->offset[j] : 0 );

		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "j:%i\n",
							dat->lengths[j]);
	}

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "max_offset: %i;"
		"max_read_length: %i;max_read_position: %i\n", max_offset,
		max_read_length, max_read_position);

	for (size_t j = 0; j < max_read_position; ++j)
		dat->coverage[j] = 0;

	for (size_t i = 0; i < dat->sample_size; ++i)
		for (size_t j = 0; j < dat->lengths[i]; ++j)
			++dat->coverage[j + (dat->offset ? dat->offset[i] : 0)];

	dat->max_offset = max_offset;
	dat->max_read_length = max_read_length;
	dat->max_read_position = max_read_position;

	/* rebuild a new hash table for the new data set */
	delete_all(dat->seq_count);
	dat->seq_count = NULL;
	dat->hash_length = 0;
	for (size_t i = 0; i < dat->sample_size; ++i)
		dat->hash_length += add_sequence(&dat->seq_count, dat->dmat[i],
							dat->lengths[i], i);

	/* store index of reads for all unique sequences */
	for (size_t i = 0; i < dat->sample_size; ++i)
		if((err = add_seq_idx(dat->seq_count, dat->dmat[i],
							dat->lengths[i], i)))
			return err;
	
	/* sort hash table by count */
	sort_by_count(&dat->seq_count);

	return err;
}/* update_data */


/**
 * Once data is loaded, build hash.
 *
 * @param dat	data pointer
 * @return	error status
 */
int sync_data(data *dat)
{
	return build_hash(dat);
} /* sync_data */
