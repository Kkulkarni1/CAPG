/**
 * @file fqmorph.c
 * @author Karin S. Dorman
 *
 * Manipulate fastq files (like fastx-toolkit)
 *
 * Note about formatting.  Line widths at 80 characters, not because
 * we live in the 60's but to help force good coding and to reduce
 * complexity.
 */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <unistd.h>

#include "fqmorph.h"
#include "fastq.h"
#include "order.h"
#include "cmdline.h"
#include "error.h"
#include "util.h"

int make_options(options **opt);
int make_data(data **in_dat);
int parse_options(options *opt, int argc, const char **argv);
void free_options(options *opt);

int sample(char const * const fastq_file, fastq_data **fd, unsigned int nsel);
int stratified_sample(options *opt, fastq_data **fqd);
int subset(char const * const filename, fastq_data **fqd, char const * const partition_file, int isel);
unsigned int *selected_group(char const * const partition_file, unsigned int nreads, int isel, unsigned int *nsel);
int read_selected_fastq(FILE *fp, fastq_data *fqd, unsigned int *selected, unsigned int nselected);
int trim_read_by_ee(fastq_data *fqd, char_t *read, char_t *qual, unsigned int idx, unsigned int len, void *vopt);
int trim_read_by_qs(fastq_data *fqd, char_t *read, char_t * qual, unsigned int idx, unsigned int len, void *vopt);
int trim_read_by_len(fastq_data *fqd, char_t *read, char_t *qual, unsigned int idx, unsigned int len, void *vopt);


int main(int argc, const char **argv)
{
	int err = EXIT_SUCCESS;	/* error code */
	options *opt = NULL;	/* run options */
	data *dat = NULL;	/* other data */
	fastq_data *fqd = NULL;	/* data object */
	fastq_options *fqo = NULL;	/* fastq options object */

	if ((err = make_options(&opt)))
		goto CLEAR_AND_EXIT;

	if ((err = make_fastq_options(&fqo)))
		goto CLEAR_AND_EXIT;
	
	opt->fqo = fqo;
	/* fqmorph defaults to XY_ENCODING */
	fqo->read_encoding = XY_ENCODING;

	if ((err = parse_options(opt, argc, argv)))
		goto CLEAR_AND_EXIT;

	/* copy fastq options into fastq options object */
	fqo->read_names = opt->read_names;
	fqo->casavize = opt->casavize;

	if ((err = make_data(&dat)))
		goto CLEAR_AND_EXIT;

	/* sample some reads */
	if (opt->nsample && !opt->partition_file
		&& (err = sample(opt->fastq_file, &fqd, opt->nsample))) {
			goto CLEAR_AND_EXIT;
	
	/* stratified sample reads from partition */
	} else if (opt->nsample && opt->partition_file
		&& (err = stratified_sample(opt, &fqd))) {
			goto CLEAR_AND_EXIT;
	
	/* select one group from partition */
	} else if (!opt->nsample && opt->s_class >= 0 && opt->partition_file
		&& (err = subset(opt->fastq_file, &fqd, opt->partition_file,
							opt->s_class))) {
			goto CLEAR_AND_EXIT;

	/* read all the data */
	} else if (!opt->nsample && opt->s_class < 0
		&& (!opt->partition_file || opt->output_format == TABLE_FORMAT)
		&& (err = read_fastq(opt->fastq_file, &fqd, fqo))) {
		goto CLEAR_AND_EXIT;

	}
	mmessage(INFO_MSG, NO_ERROR, "Number of reads: %u\n", fqd->n_reads);

	if ((opt->leading_ee > 0 || opt->trailing_ee > 0
		|| opt->leading_qs > 0 || opt->trailing_qs > 0
		|| opt->leading_nuc || opt->trailing_nuc || opt->max_length)
		&& (err = allocate_site_flag(fqd)))
		goto CLEAR_AND_EXIT;
		
	if ((opt->nsample || opt->s_class >= 0) && opt->partition_file)
		opt->partition_file = NULL;

	mmessage(INFO_MSG, NO_ERROR, "Data read.\n");

	if (opt->report_lengths)
		for (unsigned int i = 0; i < fqd->n_reads; ++i)
			fprintf(stdout, "%u\n", read_length(fqd, i));
	if (opt->report_names) {
		char *nptr = fqd->names;
		if (!opt->read_names) {
			mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"You must use --names option with --std names.\n");
			goto CLEAR_AND_EXIT;
		}
		for (unsigned int i = 0; i < fqd->n_reads; ++i)
			fprintf(stdout, "%.*s\n", fqd->name_lengths[i],
				next_read_name(fqd, (char const **) &nptr, i));
	}

	if (opt->cut_end) {
		if (opt->cut_start == 0) {
			err = mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"use 1-index with -c\n");
			goto CLEAR_AND_EXIT;
		}
		opt->cut_start--;
		opt->cut_end--;
		size_t n_bytes = 0;
		for (unsigned int i = 0; i < fqd->n_reads; ++i) {
			if (opt->cut_start > (fqd->n_lengths
				? fqd->n_lengths[i] : fqd->n_max_length)) {
				mmessage(ERROR_MSG, INVALID_USER_INPUT,
					"-c start after end of read %u\n", i);
				goto CLEAR_AND_EXIT;
			}
			if (opt->cut_end > (fqd->n_lengths
				? fqd->n_lengths[i] : fqd->n_max_length)) {
				mmessage(ERROR_MSG, INVALID_USER_INPUT,
					"-c end after end of read %u\n", i);
				goto CLEAR_AND_EXIT;
			}
			n_bytes += (fqd->n_lengths
				? fqd->n_lengths[i] : fqd->n_max_length)
				- (opt->cut_end - opt->cut_start + 1);
		}
		unsigned char *reads = malloc(n_bytes * sizeof *reads);
		unsigned char *quals = malloc(n_bytes * sizeof *reads);
		if (!reads || !quals) {
			mmessage(ERROR_MSG, MEMORY_ALLOCATION, "reads/quals");
			goto CLEAR_AND_EXIT;
		}
		unsigned char *rptr = reads, *qptr = quals;
		unsigned int idx = 0;
		for (unsigned int i = 0; i < fqd->n_reads; ++i) {
			fprintf(stderr, "Read %u idx=%u\n", i, idx);
			if (opt->cut_start) {
				memcpy(rptr, &fqd->reads[idx], opt->cut_start);
				memcpy(qptr, &fqd->quals[idx], opt->cut_start);
			}
			idx += opt->cut_end + 1;
			rptr += opt->cut_start;
			qptr += opt->cut_start;
			memcpy(rptr, &fqd->reads[idx], (fqd->n_lengths
				? fqd->n_lengths[i] : fqd->n_max_length)
				- (opt->cut_end + 1));
			memcpy(qptr, &fqd->quals[idx], (fqd->n_lengths
				? fqd->n_lengths[i] : fqd->n_max_length)
				- (opt->cut_end + 1));
			idx += (fqd->n_lengths ? fqd->n_lengths[i]
				: fqd->n_max_length)
				- (opt->cut_end + 1);
			rptr += (fqd->n_lengths ? fqd->n_lengths[i]
				: fqd->n_max_length)
				- (opt->cut_end + 1);
			qptr += (fqd->n_lengths ? fqd->n_lengths[i]
				: fqd->n_max_length)
				- (opt->cut_end + 1);
			if (fqd->n_lengths)
				fqd->n_lengths[i] = fqd->n_lengths[i] -
					(opt->cut_end - opt->cut_start + 1);
		}
		fqd->n_max_length = fqd->n_max_length - (opt->cut_end
			- opt->cut_start + 1);
		free(fqd->reads);
		fqd->reads = reads;
		free(fqd->quals);
		fqd->quals = quals;
	}

	if ((opt->max_ee > 0 || opt->name_match)
		&& (err = allocate_read_flag(fqd)))
		goto CLEAR_AND_EXIT;

	if (opt->expected_errors || opt->max_ee > 0) {
		unsigned char *qptr = fqd->quals;
		char *nptr = fqd->names;
		for (unsigned int i = 0; i < fqd->n_reads; ++i) {
			double eecnt = 0;
			for (unsigned int j = 0; j < read_length(fqd, i); ++j) {
				eecnt += error_prob(fqd, *qptr);
				++qptr;
			}
			if (!opt->expected_errors) {
				if (eecnt > opt->max_ee)
					fqd->read_flag[i] = 0;
				else
					fqd->read_flag[i] = 1;
				continue;
			}
			if (opt->read_names)
				fprintf(stdout, "Read %.*s: %f\n",
					fqd->name_lengths[i],
					next_read_name(fqd,
					(char const **) &nptr, i), eecnt);
			else
				fprintf(stdout, "Read %u: %f\n", i, eecnt);
		}
		if (opt->expected_errors)
			goto CLEAR_AND_EXIT;
	}

	/* filter out reads based on name */
	if (opt->name_match) {

		/* default: include all reads */
		for (unsigned int i = 0; i < fqd->n_reads; ++i)
			fqd->read_flag[i] = 1;

		for (unsigned int k = 0; k < opt->n_name_match; ++k) {
			char *nptr = fqd->names;
			unsigned int nread = 0;
			size_t slen = strlen(opt->name_match[k]);
			for (unsigned int i = 0; i < fqd->n_reads; ++i) {

				/* already removed */
				if (!fqd->read_flag[i]) {
					next_read_name(fqd, (char const **) &nptr, i);
					continue;
				}

				unsigned char nlen = fqd->name_lengths[i];

				/* too short to match */
				if (nlen < slen) {
					next_read_name(fqd, (char const **) &nptr, i);
					continue;
				}

				/* scan name */
				int nmatch = 0;	/* assume no match */
				for (unsigned int j = 0; j < nlen; ++j) {
					/* potential match starts here */
					if (opt->name_match[k][0] == nptr[j]) {
						nmatch = 1;
						for (unsigned int l = 0; l < slen; ++l) {
							/* potential match ends here */
							if (opt->name_match[k][l]
									!= nptr[j + l]) {
								nmatch = 0;
								break;
							}
						}
						/* match complete */
						if (nmatch) {
							++nread;
							fqd->read_flag[i] = 0;
							break;
						}
					}
				}
				next_read_name(fqd, (char const **) &nptr, i);
			}
			mmessage(INFO_MSG, NO_ERROR, "Filtering out %u reads of"
				" %u with name matching '%s'\n", nread,
					fqd->n_reads, opt->name_match[k]);
		}
	}

	if (opt->distn_file) {
		if (fqd->read_encoding != XY_ENCODING) {
			mmessage(ERROR_MSG, INTERNAL_ERROR, "Ambiguous "
						"nucleotides not allowed");
			goto CLEAR_AND_EXIT;
		}
		FILE *fp = fopen(opt->distn_file, "w");
		if (!fp) {
			mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt->distn_file);
			goto CLEAR_AND_EXIT;
		}
		unsigned int nqual = fqd->max_quality - fqd->min_quality + 1;
		double *freq = malloc(NUM_NUCLEOTIDES * nqual * sizeof *freq);
		if (!freq) {
			mmessage(ERROR_MSG, MEMORY_ALLOCATION, "freq");
			goto CLEAR_AND_EXIT;
		}
		for (unsigned int i = 0; i < nqual; ++i)
			for (unsigned int b = 0; b < NUM_NUCLEOTIDES; ++b)
				freq[i * NUM_NUCLEOTIDES + b] = 0;
		unsigned char *quals = fqd->quals;
		unsigned char *reads = fqd->reads;
		for (unsigned int i = 0; i < fqd->n_reads; ++i) {
			unsigned int len = read_length(fqd, i);
			for (unsigned int j = 0; j < len; ++j)
				freq[NUM_NUCLEOTIDES * quals[j] + reads[j]]++;
			reads += len;
			quals += len;
		}
		for (unsigned int i = 0; i < nqual; ++i) {
			double sum = 0;
			for (unsigned int b = 0; b < NUM_NUCLEOTIDES; ++b)
				sum += freq[i * NUM_NUCLEOTIDES + b];
			fprintf(fp, "%u", fqd->min_quality + i);
			for (unsigned int b = 0; b < NUM_NUCLEOTIDES; ++b) {
				freq[i * NUM_NUCLEOTIDES + b] /= sum;
				fprintf(fp, " %f", freq[i * NUM_NUCLEOTIDES + b]);
			}
			fprintf(fp, "\n");
		}
		free(freq);
		fclose(fp);
	}

	if (opt->nni_k) {
		if (opt->nni_k > 1) {
			mmessage(ERROR_MSG, INTERNAL_ERROR, "NNI requires k=1");
			goto CLEAR_AND_EXIT;
		}
		unsigned char *rptr1 = fqd->reads;
		unsigned char *qptr1 = fqd->quals;
		for (unsigned int i = 0; i < fqd->n_reads; ++i) {
			double dis = INFINITY, d;
			unsigned char *rptr2 = fqd->reads;
			unsigned char *qptr2 = fqd->quals;
			for (unsigned int j = 0; j < fqd->n_reads; ++j) {
				if (j != i) {
					unsigned int len = fqd->n_lengths ?
						MIN(fqd->n_lengths[i],
						fqd->n_lengths[j])
						: fqd->n_max_length;
					d = read_distance_ptr(fqd, len, rptr1,
						rptr2, qptr1, qptr2);
					//fprintf(stderr, "%u, %u: %.3f (%.3f)\n", i, j, d, dis);
					if (d < dis)
						dis = d;
				}
				rptr2 += fqd->n_lengths ? fqd->n_lengths[j]
					: fqd->n_max_length;
				qptr2 += fqd->n_lengths ? fqd->n_lengths[j]
					: fqd->n_max_length;
			}
			fprintf(stderr, "%e\n", dis);
			qptr1 += fqd->n_lengths ? fqd->n_lengths[i]
				: fqd->n_max_length;
			rptr1 += fqd->n_lengths ? fqd->n_lengths[i]
				: fqd->n_max_length;
		}
	}

	/* read partition file */
	if (opt->partition_file) {
		dat->cluster_id = malloc(fqd->n_reads * sizeof *dat->cluster_id);
		if (!dat->cluster_id) {
			mmessage(ERROR_MSG, MEMORY_ALLOCATION, "cluster_id");
			goto CLEAR_AND_EXIT;
		}
		err = read_uints(opt->partition_file, dat->cluster_id, fqd->n_reads);
	}

	if (opt->strip_slash || opt->strip_dash || opt->casavize) {
		char *nptr = fqd->names, *name = NULL;
		int j;
		for (unsigned int i = 0; i < fqd->n_reads; ++i) {
			name = next_read_name_rw(fqd, (char **) &nptr, i);
			for (j = fqd->name_lengths[i] - 1; j >= 0; --j)
				if ((opt->strip_slash && name[j] == '/')
					|| (opt->strip_dash && name[j] == '-')
					|| (opt->casavize
					&& (name[j] == ' ' || name[j] == ':' || name[j] == '/')))
					break;
			if (j >= 0)
				name[j] = '\0';
		}
	}

	/* get ready to output */
	fqo->reverse_complement = opt->reverse_complement;
	fqo->fasta = opt->output_format == FASTA_FORMAT;
	fqo->qs_fasta = opt->qs_to_fasta;
	fqo->append = opt->append;
	if (opt->split_file && !err) {

		unsigned max_k = 0;
		for (unsigned int i = 0; i < fqd->n_reads; ++i)
			if (dat->cluster_id[i] > max_k)
				max_k = dat->cluster_id[i];

		for (unsigned int k = 0; k < max_k; ++k) {
			char *str = malloc((strlen(opt->outfile) + 3 + k/10)
				* sizeof *str);
			if (!str) {
				mmessage(ERROR_MSG, MEMORY_ALLOCATION, "str");
				goto CLEAR_AND_EXIT;
			}
			sprintf(str, "%s.%u", opt->outfile, k);
			fqo->outfile = str;
			if (opt->leading_ee > 0 || opt->trailing_ee > 0)
				err = write_fastq_marked_trimmed(fqd, fqo,
					dat->cluster_id, k, trim_read_by_ee,
								(void *)opt);
			else if (opt->leading_qs > 0 || opt->trailing_qs > 0)
				err = write_fastq_marked_trimmed(fqd, fqo,
					dat->cluster_id, k, trim_read_by_qs,
								(void *)opt);
			else if (opt->leading_nuc || opt->trailing_nuc
				|| opt->max_length)
				err = write_fastq_marked_trimmed(fqd, fqo,
					dat->cluster_id, k, trim_read_by_len,
								(void *)opt);
			else
				err = write_fastq_marked_trimmed(fqd, fqo,
					dat->cluster_id, k, NULL, NULL);
		}
	/* request to output true cluster assignments as first column */
	} else if (opt->output_clusters && !err) {

		/* we'll open the out file ourselves */
		FILE *fp = fopen(opt->outfile, "w");
		if (!fp) {
			mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt->outfile);
			goto CLEAR_AND_EXIT;
		}

		/* violate the privacy of the fqd struct and access reads */
		unsigned char *reads = fqd->reads;
		for (unsigned int i = 0; i < fqd->n_reads; ++i) {
			unsigned int len = read_length(fqd, i);

			/* cluster assignment in first column */
			fprintf(fp, "%u ", dat->cluster_id[i]);

			/* write out read as usual */
			write_read_in_table(fp, reads, len);

			/* manually forward read pointer */
			reads += len;
		}

	} else if (!err && opt->outfile) {
		mmessage(INFO_MSG, NO_ERROR, "Will output in '%s'\n", opt->outfile);
		fqo->outfile = opt->outfile;
		if (opt->output_format != TABLE_FORMAT) {
			if (opt->leading_ee > 0 || opt->trailing_ee > 0)
				err = write_fastq_marked_trimmed(fqd, fqo,
					fqd->read_flag, 1, trim_read_by_ee,
								(void *)opt);
			else if (opt->leading_qs > 0 || opt->trailing_qs > 0)
				err = write_fastq_marked_trimmed(fqd, fqo,
					fqd->read_flag, 1, trim_read_by_qs,
								(void *)opt);
			else if (opt->leading_nuc || opt->trailing_nuc
				|| opt->max_length)
				err = write_fastq_marked_trimmed(fqd, fqo,
					fqd->read_flag, 1, trim_read_by_len,
								(void *)opt);
			else
				err = write_fastq_marked_trimmed(fqd, fqo,
					fqd->read_flag, 1, NULL, NULL);
		} else {
			err = write_table_marked(fqd, opt->outfile,
							fqd->read_flag, 1);
		}
	}

CLEAR_AND_EXIT:

	free_fastq(fqd);
	free_options(opt);
	free_fastq_options(fqo);

	if (err)
		mmessage(ERROR_MSG, err, "Error upon output: %d\n", err);

	return(err);
} /* main */

/**
 * Read a read partition file.  For each read in the fasta/fastq file, in the
 * order of the fasta/fastq file, assign it an integer cluster id.
 *
 * @param partition_file	name of partition file
 * @param nreads		size of partition
 * @param isel			selected group id
 * @param nsel			size of selected group
 * @return			ids of each entry
 */
unsigned int *selected_group(char const * const partition_file,
		unsigned int nreads, int isel, unsigned int *nsel)
{
	FILE *fp = fopen(partition_file, "r");
	int id;

	if (!fp) {
		mmessage(ERROR_MSG, FILE_OPEN_ERROR, partition_file);
		return NULL;
	}
	mmessage(INFO_MSG, NO_ERROR, "Opened partition file '%s'\n",
							partition_file);

	unsigned int *idx = malloc(nreads * sizeof *idx);
	if (!idx) {
		fclose(fp);
		mmessage(ERROR_MSG, MEMORY_ALLOCATION, "idx");
		return NULL;
	}

	*nsel = 0;
	for (unsigned int i = 0; i < nreads; ++i) {
		if (fscanf(fp, "%d", &id) != 1) {
			char c = fgetc(fp);
			if (c == 'N' && fgetc(fp) == 'A')
				continue;
			fclose(fp);
			mmessage(ERROR_MSG, FILE_FORMAT_ERROR, partition_file);
			return NULL;
		}
		if (id == isel)
			idx[(*nsel)++] = i;
	}
	fclose(fp);

	mmessage(INFO_MSG, NO_ERROR, "Found %u of %u matching reads.\n", *nsel, nreads);

	return idx;
} /* selected_group */

/**
 * Read a read partition file.  For each read in the fasta/fastq file, in that
 * order, read in its assigned integer cluster id.  The ids are
 * assumed to be NA or 0, 1, ..., where the maximum id is one less than the
 * number of distinct ids, and therefore the number of clusters.  As noted,
 * observations without cluster assignments are allowed, recorded as -1 in the
 * returned index and assumed to be NA in the file.
 *
 * @param partition_file	name of partition file
 * @param nreads		size of partition
 * @param max_id		number of groups (technically max. group id)
 * @return			ids of each entry
 */
int *read_partition_file(char const * const partition_file, unsigned int nreads,
	unsigned int *max_id)
{
	FILE *fp = fopen(partition_file, "r");

	if (!fp) {
		mmessage(ERROR_MSG, FILE_OPEN_ERROR, partition_file);
		return NULL;
	}
	fprintf(stderr, "Opened partition file '%s'\n", partition_file);

	int *idx = malloc(nreads * sizeof *idx);
	if (!idx) {
		fclose(fp);
		mmessage(ERROR_MSG, MEMORY_ALLOCATION, "idx");
		return NULL;
	}

	*max_id = 0;
	char c;
	while ((c = fgetc(fp)) == ' ');	/* fast forward through space */
	for (unsigned int i = 0; i < nreads; ++i) {
		if (c == 'N' && !feof(fp)) {
			idx[i] = -1;
			c = fgetc(fp);	/* 'A' */
			while ((c = fgetc(fp)) == ' ' && !feof(fp));
			continue;
		} else {
			ungetc(c, fp);
		}
		if (fscanf(fp, "%d", &idx[i]) != 1) {
			fclose(fp);
			mmessage(ERROR_MSG, FILE_FORMAT_ERROR, partition_file);
			return NULL;
		}
		while ((c = fgetc(fp)) == ' ' && !feof(fp));
		if ((int) *max_id < idx[i])
			*max_id = idx[i];
	}
	fclose(fp);

	return idx;
} /* read_partition_file */

int subset(char const * const filename, fastq_data **fqd, char const * const 
	partition_file, int isel)
{
//	int fxn_debug = DEBUG_I;//ABSOLUTE_SILENCE;//
	int err = NO_ERROR;
	FILE *fp = fopen(filename, "r");
	fastq_data *fqp;

	if (!fp)
		return mmessage(ERROR_MSG, FILE_OPEN_ERROR, filename);
	
	*fqd = malloc(sizeof **fqd);

	if (!*fqd)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "fastq_data");

	fqp = *fqd;

	fqp->n_reads = fcnt_reads(fp);
	mmessage(INFO_MSG, NO_ERROR, "Found %u reads in '%s'.\n", fqp->n_reads,
								filename);

	if (!fqp->n_reads)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "Failed to read "
			"input file '%s'\n", filename);
	
	//return mmessage(ERROR_MSG, INTERNAL_ERROR, "Coding not completed!");

	unsigned int nsel = 0;
	unsigned int *idx = selected_group(partition_file, fqp->n_reads, isel, &nsel);
	
	rewind(fp);
	findex_reads(fp, fqp);

	rewind(fp);
	err = read_selected_fastq(fp, fqp, idx, nsel);
	free(idx);

	return err;
} /* subset */

int sample(char const * const filename, fastq_data **fqd, unsigned int nsel)
{
	int fxn_debug = DEBUG_I;	//ABSOLUTE_SILENCE;	//
	int err = NO_ERROR;
	unsigned int i;
	fastq_data *fqp;

	FILE *fp = fopen(filename, "r");

	if (!fp)
		return mmessage(ERROR_MSG, FILE_OPEN_ERROR, filename);

	*fqd = malloc(sizeof **fqd);

	if (!*fqd)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "fastq_data");

	fqp = *fqd;

	fqp->n_reads = fcnt_reads(fp);

	if (!fqp->n_reads)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "Failed to read "
			"input file '%s'\n", filename);

	if (nsel > fqp->n_reads)
		return mmessage(ERROR_MSG, INVALID_USER_INPUT, "Cannot "
			"sample %u reads from %u total reads.\n",
			nsel, fqp->n_reads);

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Number of reads: %u\n",
								fqp->n_reads);

	rewind(fp);
	findex_reads(fp, fqp);

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Finished indexing reads\n");

	unsigned int *selected = malloc(nsel * sizeof *selected);

	if (!selected)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "selected");

	for (i = 0; i < nsel; ++i) {

		/* [KSD,TODO] This is not efficient method. */
		char repeated;
		do {
			selected[i] = rand() / (RAND_MAX + 1.) * fqp->n_reads;
			repeated = 0;
			for (unsigned int l = 0; l < i; ++l)
				if (selected[i] == selected[l]) {
					repeated = 1;
					break;
				}

			debug_msg(DEBUG_II <= fxn_debug, fxn_debug,
					"selected[%u] = %u\n", i, selected[i]);

		} while (repeated);
	}

	/* sort in increasing order */
	qsort((void *) selected, nsel, sizeof(*selected), uint_compare);

	for (i = 0; i < nsel; ++i)
		fprintf(stderr, " %u", selected[i]);
	fprintf(stderr, "\n");

	rewind(fp);
	findex_reads(fp, fqp);

	rewind(fp);
	err =  read_selected_fastq(fp, fqp, selected, nsel);
	free(selected);

	return err;
} /* sample */

int stratified_sample(options *opt, fastq_data **fqd)
{
	unsigned int nsel = opt->nsample;
	FILE *fp = fopen(opt->fastq_file, "r");

	if (!fp)
		return mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt->fastq_file);

	*fqd = malloc(sizeof **fqd);

	if (!*fqd)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "fastq_data");

	fastq_data *fqp = *fqd;

	fqp->n_reads = fcnt_reads(fp);
	fprintf(stderr, "Found %u reads in '%s'.\n", fqp->n_reads, opt->fastq_file);

	if (!fqp->n_reads)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "Failed to read "
			"input file '%s'\n", opt->fastq_file);

	unsigned int max_id = 0;

	int *idx = read_partition_file(opt->partition_file, fqp->n_reads, &max_id);

	if (!idx)
		return 1;
	
	mmessage(INFO_MSG, NO_ERROR, "Found %u clusters.\n", ++max_id);

	unsigned int *count = calloc(max_id, sizeof *count);

	if (!count) {
		free(idx);
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "count");
	}

	for (unsigned int i = 0; i < fqp->n_reads; ++i)
		if (idx[i] >= 0)
			++count[idx[i]];

	/* cumulative counts */
	unsigned int tsel = 0;
	for (unsigned int k = 0; k < max_id; ++k)
		tsel += count[k] > nsel ? nsel : count[k];

	mmessage(INFO_MSG, NO_ERROR, "Will sample %u reads.\n", tsel);

	unsigned int *selected = malloc(tsel * sizeof *selected);

	if (!selected) {
		free(count);
		free(idx);
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "selected");
	}

	for (unsigned int i = 0; i < tsel; ++i)
		selected[i] = fqp->n_reads;

	unsigned cnsel = 0;
	for (unsigned int k = 0; k < max_id; ++k) {
		unsigned int m = 0, id = 0, t = 0, cid = 0, i = 0;
		unsigned lnsel = count[k] > nsel ? nsel : count[k];
		while (m < lnsel) {
			double u = (double) rand() / RAND_MAX;

			if ( (count[k] - t) * u >= lnsel - m ) {
				++t;
			} else {
				id = t++;	/* chosen id */
				/* find it among the subset with this class */
				while (idx[i++] != (int) k || cid++ != id);
				selected[cnsel + m++] = i - 1;
				fprintf(stderr, "%uth sampled is %uth of %u in class %u: %u of %u\n", m, id, count[k], k, i, fqp->n_reads);
			}
        	}
		cnsel += lnsel;
	}

	free(count);

	/* sort in increasing order */
	qsort((void *) selected, tsel, sizeof(*selected), uint_compare);

	for (unsigned int i = 0; i < tsel; ++i)
		fprintf(stderr, " %u", selected[i]);
	fprintf(stderr, "\n");

	if (opt->partition_outfile) {
		FILE *fout = fopen(opt->partition_outfile, "w");
		if (!fout)
			mmessage(ERROR_MSG, FILE_OPEN_ERROR,
					opt->partition_outfile);
		for (unsigned int i = 0; i < tsel; ++i)
			fprintf(fout, " %d", idx[selected[i]]);
		fclose(fout);
	}

	free(idx);

	rewind(fp);
	findex_reads(fp, fqp);

	rewind(fp);
	int err =  read_selected_fastq(fp, fqp, selected, tsel);

	free(selected);

	fprintf(stderr, "Returning %d...\n", err);

	return err;
} /* stratified_sample */

int read_selected_fastq(FILE *fp, fastq_data *fqp, unsigned int *selected,
					unsigned int nsel)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	int err = NO_ERROR;
	size_t n_bytes = 0;
	unsigned int j = 0;
	unsigned int len;
	long int pos;

	for (unsigned int i = 0; i < nsel; ++i) {
		//fprintf(stderr, "Seeking %uth read %u (of %u)\n", i, selected[i], fqp->n_reads);
		while (j++ < fqp->index[selected[i]])
			fgetc(fp);
		pos = ftell(fp);
		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "current position: "
								"%ld\n", pos);
		err = read_read(fp, fqp, NULL, &len, NULL, NULL, NULL, NULL);
		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "read_read returns: "
				"%s (%d)\n", fastq_error_message(err), err);
		fseek(fp, pos + 1, SEEK_SET);	/* rewinds back through read */
		if (err == FASTQ_FILE_FORMAT_ERROR
			|| err == FASTQ_INVALID_READ_CHAR
			|| err == FASTQ_INVALID_QUALITY_CHAR)
			return FILE_FORMAT_ERROR;
		else if (err == FASTQ_EOF)
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
				"premature end-of-file\n");
		else if (err == FASTQ_AMBIGUOUS_READ_CHAR)
			return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
				"use \"--filter ambiguous\" command-line option"
				" to remove ambiguous nucleotides first\n");
		n_bytes += len;
	}

	debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "combined length: %u\n",
								 n_bytes);

	fqp->reads = malloc(n_bytes * sizeof *fqp->reads);
	if (!fqp->reads) {
		free(selected);
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "reads");
	}

	fqp->quals = NULL;
	if (fqp->file_type == FASTQ_FILE) {
		fqp->quals = malloc(n_bytes * sizeof *fqp->quals);
		if (!fqp->quals) {
			free(selected);
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "quals");
		}
	}

	fqp->n_lengths = malloc(nsel * sizeof *fqp->n_lengths);

	if (!fqp->n_lengths) {
		free(selected);
		free(fqp->reads);
		if (fqp->quals)
			free(fqp->quals);
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "n_lengths");
	}

	rewind(fp);
	j = 0;
	unsigned char *rptr = fqp->reads;
	unsigned char *qptr = fqp->quals;
	for (unsigned int i = 0; i < nsel; ++i) {
		//fprintf(stderr, "Reading %uth read %u (of %u)\n", i, selected[i], fqp->n_reads);
		while (j++ < fqp->index[selected[i]])
			fgetc(fp);
		pos = ftell(fp);
		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "current position: "
								"%ld\n", pos);
		err = read_read(fp, fqp, NULL, &fqp->n_lengths[i], rptr, qptr,
								NULL, NULL);
		debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "read_read returns: "
								"%d\n", err);
		fseek(fp, pos + 1, SEEK_SET);	/* rewinds back through read */
		rptr += fqp->n_lengths[i];
		if (qptr)
			qptr += fqp->n_lengths[i];
	}

	fqp->n_reads = nsel;
	fprintf(stderr, "Have selected and read %u reads\n", fqp->n_reads);

	return NO_ERROR;
} /* read_selected_fastq */

int trim_read_by_ee(fastq_data *fqd, char_t *read, char_t *qual,
		unsigned int idx, unsigned int len, void *vopt)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	UNUSED(read);
	UNUSED(idx);
	options *opt = (options *) vopt;
	char_t *qual_orig = qual;

	if (!fqd->site_flag)
		return INTERNAL_ERROR;

	/* trim nothing */
	memset(fqd->site_flag, 0, fqd->n_max_length * sizeof *fqd->site_flag);

	/* trim 5' */
	double eecnt = 0;
	for (unsigned int i = 1; i <= len; ++i) {
		eecnt += error_prob(fqd, *qual++);
debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Position %u (%c): %f\n", i, (char) *(qual - 1) + 33, eecnt / i);
		if (eecnt / i > opt->leading_ee) {
			fqd->site_flag[i - 1] = 1;
		} else {
debug_msg(i > 1 && DEBUG_I <= fxn_debug, fxn_debug, "Trimming %u nucs from 5' end.\n", i - 1);
			break;
		}
	}

	/* trim 3' */
	eecnt = 0;
	qual = qual_orig + len;
	for (unsigned int i = 1; i <= len; ++i) {
		eecnt += error_prob(fqd, *--qual);
debug_msg(DEBUG_I <= fxn_debug, fxn_debug, "Position %u (%c): %f\n", len - i + 1, (char) *(qual - 1) + 33, eecnt / i);
		if (eecnt / i > opt->trailing_ee)
			fqd->site_flag[len - i] = 1;
		else {
debug_msg(i > 1 && DEBUG_I <= fxn_debug, fxn_debug, "Trimming %u nucs from 3' end.\n", i - 1);
			break;
		}
	}

	return NO_ERROR;
} /* trim_read_by_ee */

int trim_read_by_qs(fastq_data *fqd, char_t *read, char_t *qual,
			unsigned int idx, unsigned int len, void *vopt)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	UNUSED(read);
	UNUSED(idx);
	options *opt = (options *) vopt;
	char_t *qual_orig = qual;

	if (!fqd->site_flag)
		return INTERNAL_ERROR;
	
	memset(fqd->site_flag, 0, fqd->n_max_length * sizeof *fqd->site_flag);

	for (unsigned int i = 0; i < len; ++i)
		if (*qual++ + fqd->min_quality - MIN_ASCII_QUALITY_SCORE
							< opt->leading_qs)
			fqd->site_flag[i] = 1;
		else {
debug_msg(i && DEBUG_I <= fxn_debug, fxn_debug, "Trimming %u nucs from 5' end\n", i);
			break;
		}

	qual = qual_orig + len;
	for (unsigned int i = 1; i <= len; ++i)
		if (*--qual + fqd->min_quality - MIN_ASCII_QUALITY_SCORE
							< opt->trailing_qs)
			fqd->site_flag[len - i] = 1;
		else {
debug_msg(i && DEBUG_I <= fxn_debug, fxn_debug, "Trimming %u nucs from 3' end\n", i - 1);
			break;
		}
	
	return NO_ERROR;
} /* trim_read_by_qs */

int trim_read_by_len(fastq_data *fqd, char_t *read, char_t *qual,
	unsigned int idx, unsigned int len, void *vopt)
{
	//int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_1;//
	UNUSED(read);
	UNUSED(qual);
	UNUSED(idx);
	options *opt = (options *) vopt;
	
	if (!fqd->site_flag)
		return INTERNAL_ERROR;

	memset(fqd->site_flag, 0, fqd->n_max_length * sizeof *fqd->site_flag);

	for (unsigned int i = 0; i < opt->leading_nuc; ++i)
		fqd->site_flag[i] = 1;
	
	if (opt->max_length && len - opt->leading_nuc - opt->trailing_nuc
							> opt->max_length) {
		for (unsigned int i = opt->leading_nuc + opt->max_length;
			i < len; ++i)
			fqd->site_flag[i] = 1;
	} else {

		for (unsigned int i = 1; i <= opt->trailing_nuc; ++i)
			fqd->site_flag[len - i] = 1;
	}

	return NO_ERROR;
} /* trim_read_by_len */

int make_data(data **in_dat)
{
	data *dat;

	*in_dat = malloc(sizeof **in_dat);
	if (!*in_dat)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data object");
	dat = *in_dat;
	dat->cluster_id = NULL;

	return NO_ERROR;
} /* make_data */

int make_options(options **opt) {
	options *op;

	*opt = malloc(sizeof **opt);
	if (!*opt)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "options object");

	op = *opt;

	op->fastq_file = NULL;
	op->outfile = NULL;
	op->fqo = NULL;
	op->distn_file = NULL;
	op->partition_file = NULL;
	op->partition_outfile = NULL;
	op->read_names = 0;
	op->s_class = -1;
	op->append = 0;
	op->split_file = 0;
	op->output_clusters = 0;
	op->output_format = FASTQ_FORMAT;
	op->reverse_complement = 0;
	op->cut_start = op->cut_end = 0;
	op->nsample = 0;
	op->expected_errors = 0;
	op->nni_k = 0;
	op->report_lengths = 0;
	op->report_names = 0;
	op->qs_to_fasta = 0;
	op->max_ee = 0;
	op->leading_ee = NAN;
	op->trailing_ee = NAN;
	op->leading_qs = 0;
	op->trailing_qs = 0;
	op->leading_nuc = 0;
	op->trailing_nuc = 0;
	op->max_length = 0;
	op->name_match = NULL;
	op->n_name_match = 0;

	op->strip_slash = 0;
	op->strip_dash = 0;
	op->casavize = 0;

	return NO_ERROR;
} /* make_options */

int parse_options(options *opt, int argc, const char **argv) {
	int i;
	size_t j;
	int err = NO_ERROR;
	char a;

	for (i = 1; i < argc; i++) {
		if (strlen(argv[i]) < 2)
			usage_error(argv, i, (void *)opt);
		j = 1;
		a = argv[i][j];
		while (a == '-' && ++j < strlen(argv[i]))
			a = argv[i][j];
		switch(a) {
			case 'a':
				opt->outfile = argv[++i];
				opt->append = 1;
				break;
			case 'c':
				opt->cut_start = read_uint(argc, argv, ++i,
					(void *)opt);
				opt->cut_end = read_uint(argc, argv, ++i,
					(void *)opt);
				break;
			case 'e':
				if (!strncmp(&argv[i][j], "enc", 3)) {
					if (i + 1 == argc)
						goto CMDLINE_ERROR;
					if (!strncmp(argv[++i], "iu", 2))
						opt->fqo->read_encoding
							= IUPAC_ENCODING;
				} else {
					opt->expected_errors = 1;
				}
				break;
			case 'r':
				opt->reverse_complement = 1;
				if (errno)
					goto CMDLINE_ERROR;
				break;
			case 'd':
				if (i == argc - 1)
					goto CMDLINE_ERROR;
				opt->distn_file = argv[++i];
				break;
			case 'f':
				if (!strcmp(&argv[i][j], "filter")) {
					if (i + 1 == argc)
						goto CMDLINE_ERROR;
					if (!strncmp(argv[i + 1], "amb", 3)) {
						++i;
						/* drop ambig nucs */
						if (i + 1 < argc && !strncmp(
							argv[i + 1], "nuc", 3)) {
							opt->fqo->drop_ambiguous_nucs = 1;
							++i;
						/* drop reads w/ ambig nucs */
						} else if (i + 1 < argc &&
							!strncmp(argv[i + 1], "rea", 3)) {
							opt->fqo->drop_ambiguous_reads = 1;
							++i;
						/* default is drop reads */
						} else {
							opt->fqo->drop_ambiguous_reads = 1;
						}
					} else if (!strncmp(argv[i + 1], "len", 3)) {
						++i;
						if (i + 1 < argc && !strncmp(
							argv[i + 1], "max", 3)) {
							if (i + 2 == argc)
								goto CMDLINE_ERROR;
							opt->fqo->max_length = strtoul(argv[i + 2], NULL, 0);
							i += 2;
						} else if (i + 1 < argc && !strncmp(
							argv[i + 1], "min", 3)) {
							if (i + 2 == argc)
								goto CMDLINE_ERROR;
							opt->fqo->min_length = strtoul(argv[i + 2], NULL, 0);
							i += 2;
						} else {
							goto CMDLINE_ERROR;
						}
					} else if (!strncmp(argv[i + 1], "exp", 3)) {
						opt->max_ee =
							read_cmdline_double(
							argc, argv, i + 2,
								(void *)opt);
						mmessage(INFO_MSG, NO_ERROR, 
							"Filtering on fewer "
							"than %f expected "
							"errors.\n",
							opt->max_ee);
						i += 2;
					} else if (!strncmp(argv[i + 1], "nam", 3)) {
						if (i + 2 >= argc)
							goto CMDLINE_ERROR;
						++opt->n_name_match;
						opt->read_names = 1;
						char const **nnm = realloc(
							opt->name_match,
							opt->n_name_match
							* sizeof *opt->name_match);
						if (!nnm)
							goto CMDLINE_ERROR;
						opt->name_match = nnm;
						opt->name_match[opt->n_name_match - 1]
							= argv[i + 2];
						i += 2;
					} else {
						goto CMDLINE_ERROR;
					}
				} else {
					opt->output_format = FASTA_FORMAT;
				}
				break;
			case 'i':
				if (i == argc - 1)
					goto CMDLINE_ERROR;
				opt->fastq_file = argv[++i];
				break;
			case 'o':
				if (i == argc - 1)
					goto CMDLINE_ERROR;
				opt->outfile = argv[++i];
				break;
			case 'p':
				if (i == argc - 1)
					goto CMDLINE_ERROR;
				opt->partition_file = argv[++i];
				opt->split_file = 1;
				break;
			case 'q':
				opt->qs_to_fasta = 1;
				mmessage(INFO_MSG, NO_ERROR, "Outputting "
					"quality scores.\n");
				break;
			case 'n':
				if (!strncmp(&argv[i][j], "na", 1)) {
					opt->read_names = 1;
					mmessage(INFO_MSG, NO_ERROR,
						"Will retain read names.\n");
					if (i + 1 == argc)
						continue;
					if (!strncmp(argv[i + 1], "strip-s",
									7)) {
						opt->strip_slash = 1;
						mmessage(INFO_MSG, NO_ERROR,
							"Will strip /# from "
							"names.\n");
						++i;
					} else if (!strncmp(argv[i + 1],
								"strip-d", 7)) {
						opt->strip_dash = 1;
						mmessage(INFO_MSG, NO_ERROR,
							"Will strip -* from "
							"names.\n");
						++i;
					} else if (!strncmp(argv[i + 1], "cas",
									3)) {
						if (i + 1 == argc ||
							argv[i + 2][0] == '-') {
							err = INVALID_CMD_ARGUMENT;
							goto CMDLINE_ERROR;
						}
						opt->casavize = atoi(argv[i + 2]);
						mmessage(INFO_MSG, NO_ERROR,
							"Will \"casavize\" "
							"the names.\n");
						i += 2;
					}
				} else {
					if (i == argc - 1)
						goto CMDLINE_ERROR;
					opt->nni_k = read_uint(argc, argv, ++i,
						(void *)opt);
				}
				break;
			case 's':
				if (!strncmp(&argv[i][j], "se", 2)) {
					opt->seed = read_uint(argc, argv, ++i,
						(void *)opt);
					srand(opt->seed);
				} else if (!strncmp(&argv[i][j], "su", 2)) {
					if (access(argv[i + 1], F_OK) != -1) {
						opt->partition_file = argv[++i];
						opt->s_class = read_uint(argc,
							argv, ++i, (void *)opt);
						mmessage(INFO_MSG, NO_ERROR, 
							"Class file: '%s'; "
							"Index: %u\n",
							opt->partition_file,
							opt->s_class);
					} else {
						err = INVALID_CMD_ARGUMENT;
						goto CMDLINE_ERROR;
					}
				} else if (!strncmp(&argv[i][j], "sa", 2)) {
					if (access(argv[i + 1], F_OK) != -1) {
						opt->partition_file = argv[++i];
						opt->nsample = read_uint(argc,
							argv, ++i, (void *)opt);
						if (argc > i
							&& argv[i][0] != '-')
							opt->partition_outfile
								= argv[++i];
						mmessage(INFO_MSG, NO_ERROR, 
							"Class file: '%s'; "
							"Sample: %u; Out class "
							"file: '%s'\n",
							opt->partition_file,
							opt->nsample,
							opt->partition_outfile ?
							opt->partition_outfile :
							NULL);
					} else {
						opt->nsample = read_uint(argc,
							argv, ++i, (void *)opt);
						mmessage(INFO_MSG, NO_ERROR,
							"Subsample %u\n",
							opt->nsample);
					}
				} else if (!strncmp(&argv[i][j], "st", 2)) {
					++i;
					if (i == argc)
						goto CMDLINE_ERROR;
					if (!strncmp(argv[i], "l", 1)) {
						opt->report_lengths = 1;
					} else if (!strncmp(argv[i], "n", 1)) {
						opt->report_names = 1;
					} else {
						err = INVALID_CMD_OPTION;
						goto CMDLINE_ERROR;
					}
				} else {
					err = INVALID_CMD_OPTION;
					goto CMDLINE_ERROR;
				}
				break;
			case 't':
				if (!strncmp(&argv[i][j], "ta", 2)) {
					opt->output_format = TABLE_FORMAT;
					/* user has provided an argument */
					if (i + 1 < argc && access(argv[i + 1], F_OK) != -1) {
						opt->partition_file = argv[++i];
						opt->output_clusters = 1;
					}
				} else if (!strncmp(&argv[i][j], "tr", 2)) {
					if (i + 1 == argc) {
						err = INVALID_CMD_ARGUMENT;
						goto CMDLINE_ERROR;
					}
					if (!strncmp(argv[i + 1], "exp", 3)) {
						++i;
						if (i + 1 < argc && !strncmp(
							argv[i + 1], "l", 1)) {
							if (i + 2 == argc)
								goto CMDLINE_ERROR;
							opt->leading_ee = 
								read_cmdline_double(
								argc, argv, i + 2,
									(void *)opt);
							mmessage(INFO_MSG, NO_ERROR,
								"Trimming "
								"maximum 5' "
								"window with "
								"average %f "
								"errors per "
								"site.\n",
								opt->leading_ee);
							i += 2;
						} else if (i + 1 < argc && !strncmp(
							argv[i + 1], "t", 1)) {
							if (i + 2 == argc)
								goto CMDLINE_ERROR;
							opt->trailing_ee = 
								read_cmdline_double(
								argc, argv, i + 2,
									(void *)opt);
							mmessage(INFO_MSG, NO_ERROR,
								"Trimming "
								"maximum 3' "
								"window with "
								"average %f "
								"errors per "
								"site.\n",
								opt->trailing_ee);
							i += 2;
						}
					} else if (!strncmp(argv[i + 1], "q", 1)) {
						++i;
						if (i + 1 < argc && !strncmp(
							argv[i + 1], "l", 1)) {
							if (i + 2 == argc)
								goto CMDLINE_ERROR;
							opt->leading_qs = atoi(argv[i+2]);
							mmessage(INFO_MSG, NO_ERROR,
								"Trimming "
								"5' nucleotides "
								"with quality "
								"score less "
								"than %d.\n",
								(int)opt->leading_qs);
							i += 2;
						} else if (i + 1 < argc && !strncmp(
							argv[i + 1], "t", 1)) {
							if (i + 2 == argc)
								goto CMDLINE_ERROR;
							opt->trailing_qs = atoi(argv[i+2]);
							mmessage(INFO_MSG, NO_ERROR,
								"Trimming "
								"3' nucleotides "
								"with quality "
								"score less "
								"than %d.\n",
								(int)opt->trailing_qs);
							i += 2;
						}
					} else if (!strncmp(argv[i + 1], "l", 1)) {
						++i;
						if (i == argc)
							goto CMDLINE_ERROR;
						if (!strncmp(argv[i + 1], "t", 1)) {
							if (i + 2 == argc)
								goto CMDLINE_ERROR;
							opt->trailing_nuc = atoi(argv[i+2]);
							mmessage(INFO_MSG, NO_ERROR,
								"Trimming %u 3'"
								"nucleotides.\n",
								opt->trailing_nuc);
							i += 2;
						} else if (!strncmp(argv[i + 1], "m", 1)) {
							if (i + 2 == argc)
								goto CMDLINE_ERROR;
							opt->max_length = atoi(argv[i+2]);
							mmessage(INFO_MSG, NO_ERROR,
								"Trimming to max "
								"length %u.\n",
								opt->max_length);
							i += 2;
						} else {
							++i;
							if (i == argc)
								goto CMDLINE_ERROR;
							opt->leading_nuc = atoi(argv[i]);
							mmessage(INFO_MSG, NO_ERROR,
								"Trimming %u 5'"
								"nucleotides.\n",
								opt->leading_nuc);
							j = 0;
							size_t len = strlen(argv[i]);
							while (j < len && argv[i][j] != ',')
								++j;
							if (j < len) {
								opt->trailing_nuc = atoi(&argv[i][j+1]);
								mmessage(INFO_MSG, NO_ERROR,
									"Trimming %u 3'"
									"nucleotides.\n",
									opt->trailing_nuc);
							}
						}
					}
				} else {
					err = INVALID_CMD_OPTION;
					goto CMDLINE_ERROR;
				}
				break;
			case 'h':
				fprint_usage(stderr, argv[0], opt);
				free_options(opt);
				exit(EXIT_SUCCESS);
			default:
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
		}
	}

	if (opt->fqo->drop_ambiguous_reads && opt->nsample)
		return mmessage(ERROR_MSG, INVALID_CMDLINE, "Cannot use "
			"\"--filter ambiguous\" and \"--select\" on the "
			"same command-line.\n");

	return err;

CMDLINE_ERROR:
	if (err == NO_ERROR) {
		err = INVALID_CMD_ARGUMENT;
		i--;
	}
	usage_error(argv, i, (void *)opt);
	return err;
} /* parse_options */

void free_options(options *opt) {
	if (opt) free(opt);
} /* free_options */

const char * fqmorph_error_message(int err_no) {
	if (err_no == FQMORPH_AMBIGUOUS_NUCLEOTIDE)
		return "ambiguous nucleotide";
	else if (err_no)
		return "unknown error";
	else
		return "";
} /* fqmorph_error_message */

void fprint_usage(FILE *fp, const char *cmdname, void *obj) {
	options *opt = (options *) obj;
	size_t start = strlen(cmdname) - 1;

	while (cmdname[start] != '/' && start) start--;
	if (cmdname[start] == '/') start++;

	for (size_t i = start; i < strlen(cmdname); ++i)
		fputc(toupper(cmdname[i]), fp);
	fprintf(fp, "(%d)\n", 1);
	fprintf(fp, "\nNAME\n\t%s - morph fastq files\n",
		&cmdname[start]);
	fprintf(fp, "\nSYNOPSIS\n\t%s [-see <see> -sel <sel> -r -f [-q] -c <c1> <c2>] -i <i> -[a|o] <o>...\n",
		&cmdname[start]);
	fprintf(fp, "\nDESCRIPTION\n\t%s maniplates fastq files (-i) in "
		"various ways, including to\n\treverse complement (-r), ouput "
		" as fasta (-f), cut sites (-c),\n\trandom select a subset of size <sel>, etc.\n",
		&cmdname[start]);
	fprintf(fp, "\nOPTIONS\n");
	fprintf(fp, "\t-a <ostr>\n\t\tAppend to this output file. [-a|-o required]\n");
	fprintf(fp, "\t-c <n> <m>\n\t\tCut sites <n> to <m> inclusive (1-index) out of each read.\n");
	fprintf(fp, "\t-d <dstr>\n\t\tBase per quality score distribution written to file <dstr>.\n");
	fprintf(fp, "\t-e\n\t\tCompute expected error per read.\n");
	fprintf(fp, "\t-encoding iupac|xy\n\t\tChoose read encoding, iupac for ambiguous nucleotides.\n");
	fprintf(fp, "\t-f\n\t\tOutput fasta format. [Default: %u]\n", opt->output_format == FASTA_FORMAT);
	fprintf(fp, "\t-q\n\t\tOutput quality scores in fasta format. [Default: %u]\n", opt->qs_to_fasta);
	fprintf(fp, "\t--filter <fstr>\n\t\tApply filter.\n");
	fprintf(fp, "\t\t\tambiguous reads: remove reads with ambiguous characters, such as 'N'\n");
	fprintf(fp, "\t\t\tambiguous nucleotides: remove nucleotides with ambiguous characters, such as 'N'\n");
	fprintf(fp, "\t\t\tlen[gth] max <n>: remove reads longer than given length\n");
	fprintf(fp, "\t\t\tlen[gth] min <n>: remove reads shorter than given length\n");
	fprintf(fp, "\t\t\texp[ected-errors] <f>: remove reads with more than <f> expected errors\n");
	fprintf(fp, "\t\t\tname <nstr>: remove reads with names matching <nstr>\n");
	fprintf(fp, "\t--trim <tstr>\n\t\tRead trimming.\n");
	fprintf(fp, "\t\t\texp[ected-errors] l[eading] <f>: trim maximum 5' window with more than <f> expected errors/site\n");
	fprintf(fp, "\t\t\texp[ected-errors] t[railing] <f>: trim maximum 3' window with more than <f> expected errors/site\n");
	fprintf(fp, "\t\t\tl[ength] <i1>[,<i2>]: trim leading <i1> nucleotides and trailing <i2> nucleotides\n");
	fprintf(fp, "\t\t\tl[ength] t[railing] <i1>: trim trailing <i1> nucleotides\n");
	fprintf(fp, "\t\t\tl[ength] m[ax] <i1>: trim reads to maximum length <i1> nucleotides\n");
	fprintf(fp, "\t-r\n\t\tOutput reverse complement reads. [Default: %u]\n", opt->reverse_complement);
	fprintf(fp, "\t-i <istr>\n\t\tSet the fastq input file.\n");
	fprintf(fp, "\t--st[dout] lengths|names\n\t\tReport lengths or names of all reads, one per line of stdout.\n");
	fprintf(fp, "\t-o <ostr>\n\t\tWrite to this output file. [-a|-o required]\n");
	fprintf(fp, "\t-p <pstr>\n\t\tUse <pstr> partition file to split fastq file\n");
	fprintf(fp, "\t--names\n\t\tRetain fastq read names and optionally manipulate them:\n");
	fprintf(fp, "\t\t\tstrip-slash: remove trailing /# from read names\n");
	fprintf(fp, "\t\t\tstrip-dash: remove trailing -.* from read names\n");
	fprintf(fp, "\t\t\tcasavize <n>: convert to casava 1.8 read ids with read number <n>\n");
	fprintf(fp, "\t-n <nint>\n\t\tOutput distance to <nint>-nearest neighbor\n");
	fprintf(fp, "\t--table [<tstr>]\n\t\tOutput table format, optionally placing cluster assignments from <tstr> in column 0. [Default: %u]\n", opt->output_format == TABLE_FORMAT);
	fprintf(fp, "\t-sa|--sample <suint>\n\t\tSelect random sample of this size. [OPTIONAL]\n");
	fprintf(fp, "\t-sa|--sample <sfile> <sint>\n\t\tSelect <sint> sample from each category in <sfile> which is merely partition as produced by ampliclust. [OPTIONAL]\n");
	fprintf(fp, "\t-se|--seed <sulong>\n\t\tSet random number seed. [OPTIONAL]\n");
	fprintf(fp, "\t-su|--subset <sfile> <sint>\n\t\tSelect reads with category <sint> in <sfile> which is merely partition as produced by ampliclust.\n");
	fprintf(fp, "\t-h\n\t\tThis help.\n");
	fprintf(fp, "\n");
	for (size_t i = start; i < strlen(cmdname); ++i)
		fputc(toupper(cmdname[i]), fp);
	fprintf(fp, "(%d)\n", 1);
} /* fprint_usage */
