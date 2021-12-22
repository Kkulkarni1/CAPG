/**
 * @file trim.c
 * @author Karin S. Dorman
 *
 * Trim the adapter, identify the barcode and remove.
 *
 * Note about formatting.  Line widths at 80 characters, not because
 * we live in the 60's but to help force good coding and to reduce
 * complexity.
 */

#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#include "trim.h"
#include "fastq.h"
#include "cmdline.h"
#include "error.h"
#include "io.h"

/* this object file will contain these inline function bodies */
extern int nuccmp(char_t nuc, data *dat, size_t i, size_t j);
extern size_t err_idx(iupac_t nuc, data *dat, size_t i, size_t j);
extern size_t nuc_idx(data *dat, size_t i, size_t j);
extern char nuc_char(data *dat, size_t i, size_t j);
extern void first_sequence(char_t * const cseq, char_t * const seq, size_t len);
extern int next_sequence(char_t * const cseq, char_t * const seq, size_t len);





int make_options(options **opt);
int parse_options(options *opt, int argc, const char **argv);
int sync_options(options *opt);
void free_options(options *opt);
int make_data(data **data, options *opt);
int finish_make_data(data *dat, options *opt);
void free_data(data *dat, options *opt);
int make_model(model **mod, data *dat, options *opt);
int initialize_model(model *mod, data *dat, options *opt);
void free_model(model *mod);

int trim_5primer(fastq_data *fqd, char_t *read, char_t *qual, unsigned int idx, unsigned int len, void *vopt);
int trim_to_barcode(fastq_data *fqd, char_t *read, char_t *qual, unsigned int idx, unsigned int len, void *vopt);


int main(int argc, const char **argv)
{
	int err = NO_ERROR;	/* error code */
	options *opt = NULL;	/* run options */
	fastq_options *fqo = NULL;/* fastq options */
	data *dat = NULL;	/* data object */
	model *mod = NULL;	/* model object */

	if ((err = make_options(&opt)))
		goto CLEAR_AND_EXIT;

	if ((err = parse_options(opt, argc, argv)))
		goto CLEAR_AND_EXIT;

	if ((err = make_data(&dat, opt)))
		goto CLEAR_AND_EXIT;

	if ((err = make_fastq_options(&fqo)))
		goto CLEAR_AND_EXIT;

	//fqo->read_encoding = XY_ENCODING;
	fqo->read_encoding = IUPAC_ENCODING;
	fqo->read_names = 1;

	if ((err = read_fastq(opt->fastq_file, &dat->fqd, fqo)))
		goto CLEAR_AND_EXIT;

	if ((err = allocate_read_flag(dat->fqd)))
		goto CLEAR_AND_EXIT;

	if ((err = allocate_site_flag(dat->fqd)))
		goto CLEAR_AND_EXIT;

	if ((err = finish_make_data(dat, opt)))
		goto CLEAR_AND_EXIT;

	if ((err = make_model(&mod, dat, opt)))
		goto CLEAR_AND_EXIT;


	/* repeated initializations of EM */
	mod->best_ll = -INFINITY;
	for (unsigned int i = 0; i < opt->n_init; ++i) {
		initialize_model(mod, dat, opt);

		em(dat, mod, opt);

		/* better solution! */
		if (mod->best_ll < mod->ll) {
			/* assign observations to clusters */
			mod->ll = e_step(dat, mod, opt);
			assign_index(dat, mod, opt, OPTIMAL);
			param_update(mod, opt, BEST);

			/* output information about better solution */
			fprintf(stderr, "pi: ");
			fprint_doubles(stderr, mod->best_pi, opt->n_indices + 1, 3, 1);
			fprintf(stderr, "delta: ");
			fprint_doubles(stderr, mod->best_delta, opt->primer_len, 3, 1);
			fprintf(stderr, "gamma:\n");
			print_gamma(stderr, mod->best_gamma, NUM_NUCLEOTIDES);
		}

		/* friendly output */
		mmessage(INFO_MSG, NO_ERROR, "");
		fprintf(stderr, ": %.3f (%.3f)\n", mod->ll, mod->best_ll);

	}

	fprint_doubles(stdout, mod->best_delta, opt->primer_len, 3, 1);
	fprint_doubles(stdout, mod->best_gamma, NUM_NUCLEOTIDES*NUM_NUCLEOTIDES, 3, 1);
	fprint_doubles(stdout, mod->best_beta, NUM_NUCLEOTIDES, 3, 1);
	fprint_doubles(stdout, mod->best_eta, NUM_NUCLEOTIDES, 3, 1);
	fprint_doubles(stdout, mod->best_pi, opt->n_indices + 1, 3, 1);

	if (0) {	/* for homework solution */
	for (unsigned int i = 0; i < dat->fqd->n_reads; ++i) {
		double max = 0, sum = 0;
		size_t max_state = 0;
		for (size_t k = 0; k < opt->n_states; ++k) {
			sum += mod->eik[i*opt->n_states + k];
			if (max < mod->eik[i*opt->n_states + k]) {
				max = mod->eik[i*opt->n_states + k];
				max_state = k;
			}
		}
		if (max_state == opt->n_states - 1)
			fprintf(stdout, "-1");
		else
			fprintf(stdout, "%zu", max_state);// % opt->n_primers);
/*
		for (size_t k = 0; k < opt->n_states; ++k)
			fprintf(stdout, " %f", mod->eik[i*opt->n_states +k]/sum);
*/
		fprintf(stdout, "\n");
	}
	}

/*
*/
	if (1) {
	/* string 160 nucs after primer */
	char const * const end_seq = "TCTTATACTTTTTCTTGTATTATTGTTGGGTCTTGTACAATTAATCCCTACAGATTCATTGAGATGTACTATTATTGTTTTGACATTGT";
	char_t *end = NULL;
	unsigned int *primer_len = malloc(dat->fqd->n_reads * sizeof *primer_len);
	unsigned int *posn_idx = malloc(dat->fqd->n_reads * sizeof *posn_idx);
	if (0) {
		end = malloc((strlen(end_seq) + 1) * sizeof *end);
		memcpy(end, (char_t *)end_seq, strlen(end_seq));
//		strncpy(end, end_seq, strlen(end_seq) + 1);
		for (size_t i = 0; i < strlen(end_seq); ++i)
			end[i] = nuc_to_iupac[end[i] - 'A'];
	}
	double err_prob = 0.02;
	double *pik = malloc((opt->n_indices + 1) * sizeof *pik);
	char *nptr = dat->fqd->names;
	for (unsigned int i = 0; i < dat->fqd->n_reads; ++i) {
		dat->fqd->read_flag[i] = 0;
		primer_len[i] = 0;
		//fprintf(stderr, "%u", i);
		for (unsigned int k = 0; k <= opt->n_indices; ++k)
			pik[k] = 0;
		double max = 0, sum = 0;
		size_t max_state = 0;
		for (size_t k = 0; k < opt->n_states; ++k) {
			if (k == opt->n_states - 1)
				pik[opt->n_indices] += mod->eik[i * opt->n_states + k];
			else
				pik[k / opt->n_primers] += mod->eik[i * opt->n_states + k];
			sum += mod->eik[i*opt->n_states + k];
			if (max < mod->eik[i*opt->n_states + k]) {
				max = mod->eik[i*opt->n_states + k];
				max_state = k;
			}
		}
		size_t primer_idx = max_state % opt->n_primers;		/* primer */
		unsigned int position_idx = max_state / opt->n_primers;	/* index */

		/* temporary: print out trimmed read */
if (0) {
		if (max_state == opt->n_states - 1)
			fprintf(stdout, "%s\n", display_sequence(dat->reads[i], read_length(dat->fqd, i), IUPAC_ENCODING));
		else if (read_length(dat->fqd, i) > 52 + position_idx)
			fprintf(stdout, "%s\n", display_sequence(&dat->reads[i][52 + position_idx],  read_length(dat->fqd, i) - position_idx - 52, IUPAC_ENCODING));
		else
			fprintf(stdout, "\n");
		continue;
}

		//fprintf(stderr, " %2u %2u %2u :", max_state, k, l);
		double max_position = 0;
		unsigned int aggregate_position_idx = 0;
		for (unsigned int k = 0; k <= opt->n_indices; ++k)
			if (max_position < pik[k]) {
				max_position = pik[k];
				aggregate_position_idx = k;
			}
if (0) {	/* BCBio444 */
		if (read_length(dat->fqd, i) < 52 || max_state == opt->n_states - 1)
			fprintf(stdout, "-1\n");
		else
			fprintf(stdout, "%u\n", aggregate_position_idx);
		continue;
}
		double ll = 0;
		if (max_state == opt->n_states - 1) {
			ll = log(mod->pi[opt->n_indices]);
			for (unsigned int j = 0; j < read_length(dat->fqd, i); ++j)
				ll += log(mod->eta[nuc_idx(dat, i, j)]);
/*
*/
			char *name = next_read_name_rw(dat->fqd, &nptr, i);
			for (unsigned int j = 0; j < dat->fqd->name_lengths[i]; ++j)
				if (name[j] == ' ') {
					name[j] = '\0';
					break;
				}
			fprintf(stdout, "%6u %-18.*s %3u %2u %-52s %-9s %-10s", i, strlen(name) < 18 ? (int)strlen(name) : 18, name,
				read_length(dat->fqd, i), aggregate_position_idx,
				"NA", "NA", "NA");
			for (size_t j = 0; j < opt->n_states; ++j)
				fprintf(stdout, " %4.2f", mod->eik[i*opt->n_states + j]);
			fprintf(stdout, " %9.3f %8s %8s %8s %7s\n", ll, "NA", "NA", "NA", "NA");
		} else {
			if (read_length(dat->fqd, i) < opt->primer_len) {
				next_read_name(dat->fqd, (char const **) &nptr, i);
				continue;
			}

			first_sequence(opt->cprimer, opt->primer, opt->primer_len);
			ll = log(mod->pi[position_idx]) - opt->n_primers*log(opt->n_primers);
			double llp = log(mod->pi[position_idx]) - opt->n_primers*log(opt->n_primers);
			double llp5 = 0, llp3 = 0;
			if (!isfinite(llp))
				exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "llp: %f = %f - %f\n", llp,
					log(mod->pi[position_idx]), opt->n_primers*log(opt->n_primers)));
			do {
//fprintf(stderr, "cprimer = %s\n", display_sequence(opt->cprimer, opt->primer_len, IUPAC_ENCODING));
				int n = -position_idx;
				int past_barcode = 0;
				for (unsigned int j = 0; j < read_length(dat->fqd, i); ++j) {
					if (n >= 0 && (unsigned int)n < opt->primer_len) {
//fprintf(stderr, "Primer position %u matching site %u cp=%d p=%d\n", n, j, opt->cprimer[n], opt->primer[n]);
						if (popcnt[(int)opt->primer[n]] < NUM_NUCLEOTIDES) {
							if (nuccmp(opt->cprimer[n], dat, i, j)) {
								ll += log(1 - mod->delta[j]) + log(mod->gamma[err_idx(opt->cprimer[n], dat, i, j)]);
								llp += log(1 - mod->delta[j]) + log(mod->gamma[err_idx(opt->cprimer[n], dat, i, j)]);
								if (past_barcode)
									llp3 += log(1 - mod->delta[j]) + log(mod->gamma[err_idx(opt->cprimer[n], dat, i, j)]);
								else
									llp5 += log(1 - mod->delta[j]) + log(mod->gamma[err_idx(opt->cprimer[n], dat, i, j)]);
								if (!isfinite(llp))
									exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "llp[%u %u %u]: %f += %f + %f (%f)\n", n, j, opt->primer_len, llp, log(1 - mod->delta[j]), log(mod->gamma[err_idx(opt->cprimer[n], dat, i, j)]), mod->delta[j]));
							} else {
								ll += log(mod->delta[j]);
								llp += log(mod->delta[j]);
								if (past_barcode)
									llp3 += log(mod->delta[j]);
								else
									llp5 += log(mod->delta[j]);
								if (!isfinite(llp))
									exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "llp[%u]: %f += %f\n", j, llp, log(mod->delta[j])));
							}
						} else {	/* barcode */
							ll += log(mod->beta[nuc_idx(dat, i, j)]);
							past_barcode = 1;
						}
					} else if (n < 0) {
						ll += log(mod->beta[nuc_idx(dat, i, j)]);
					} else
						ll += log(mod->eta[nuc_idx(dat, i, j)]);
					n++;
				}
			} while (next_sequence(opt->cprimer, opt->primer, opt->primer_len));
			first_sequence(opt->cprimer, opt->primer, opt->primer_len);
			for (size_t j = 0; j < primer_idx; ++j)
				next_sequence(opt->cprimer, opt->primer, opt->primer_len);

			double ll_end = 0;
			unsigned int cnt = 0;
			if (end) {
				for (unsigned int n = 0, j = 160 + opt->primer_len + aggregate_position_idx; n < opt->primer_len && j < read_length(dat->fqd, i); ++j, ++n) {
					if (nuccmp(end[n], dat, i, j))
						ll_end += log(err_prob) + log(mod->gamma[err_idx(end[n], dat, i, j)]);
					else
						ll_end += log(1 - err_prob);
					cnt++;
				}
			}
			char *name = next_read_name_rw(dat->fqd, &nptr, i);
			for (unsigned int j = 0; j < dat->fqd->name_lengths[i]; ++j)
				if (name[j] == ' ') {
					name[j] = '\0';
					break;
				}
			fprintf(stdout, "%6u %-18.*s %3u %2u(%2u) %-52s %9s %s", i, strlen(name) < 18 ? (int)strlen(name) : 18, name,
				read_length(dat->fqd, i), aggregate_position_idx, position_idx,
//				display_sequence(opt->cprimer, opt->primer_len, fqo->read_encoding),
				display_sequence(&dat->reads[i][position_idx], read_length(dat->fqd, i) < opt->primer_len ? read_length(dat->fqd, i) : opt->primer_len, fqo->read_encoding),
				display_sequence(&dat->reads[i][18 + position_idx], 9, fqo->read_encoding),
				display_sequence(&dat->reads[i][230 + 52 + position_idx], 10, fqo->read_encoding));
			for (size_t j = 0; j < opt->n_states; ++j)
				fprintf(stdout, " %4.2f", mod->eik[i*opt->n_states + j]);
			fprintf(stdout, " %9.3f %8.3f %8.3f %8.3f", ll, llp, llp5, llp3);
			if (cnt)
				fprintf(stdout, " %7.3f\n", ll_end/cnt);
			else
				fprintf(stdout, " %7s\n", "NA");
			if (llp3 > -100) {	/* WARNING: hard-coded */
				dat->fqd->read_flag[i] = 1;
				primer_len[i] = position_idx + opt->primer_len;
				posn_idx[i] = position_idx;
			}
		}
	}
	free(end);
	if (opt->trimmed_outfile || opt->trim_outfile) {
		fqo->outfile = opt->trimmed_outfile;
		if ((err = write_fastq_marked_trimmed(dat->fqd, fqo,
			dat->fqd->read_flag, 1, trim_5primer, (void *)primer_len)))
			goto CLEAR_AND_EXIT;
		fqo->outfile = opt->trim_outfile;
		if ((err = write_fastq_marked_trimmed(dat->fqd, fqo,
			dat->fqd->read_flag, 1, trim_to_barcode, (void *)posn_idx)))
			goto CLEAR_AND_EXIT;
	}
	}

CLEAR_AND_EXIT:

	free_data(dat, opt);
	free_options(opt);
	free_model(mod);

	return(EXIT_FAILURE);
} /* main */

/**
 * Trim 5' primer from read.
 *
 * @param fqd	fastq data structure
 * @param read	pointer to current read
 * @param qual	pointer to current quality scores
 * @param idx	index of current read
 * @param len	length of current read
 * @param vopt	optional information passed as a void pointer
 * @return	error status
 */
int trim_5primer(fastq_data *fqd, char_t *read, char_t *qual,
			unsigned int idx, unsigned int len, void *vopt)
{
	if (!fqd->site_flag)
		return INTERNAL_ERROR;

	UNUSED(read);
	UNUSED(qual);
	UNUSED(len);
	unsigned int *primer_len = (unsigned int *) vopt;

	/* trim nothing */
	memset(fqd->site_flag, 0, fqd->n_max_length * sizeof *fqd->site_flag);

	for (unsigned i = 0; i < primer_len[idx]; ++i)
		fqd->site_flag[i] = 1;

	return NO_ERROR;
} /* trim_5primer */

int trim_to_barcode(fastq_data *fqd, char_t *read, char_t *qual,
		unsigned int idx, unsigned int len, void *vopt)
{
	if (!fqd->site_flag)
		return INTERNAL_ERROR;
	
	UNUSED(read);
	UNUSED(qual);
	unsigned int *posn_idx = (unsigned int *) vopt;

	memset(fqd->site_flag, 0, fqd->n_max_length * sizeof *fqd->site_flag);

	unsigned int bc_pos = posn_idx[idx] + 18;	/* WARNING: hard-coded */
	for (unsigned i = 0; i < bc_pos; ++i)
		fqd->site_flag[i] = 1;
	for (unsigned i = bc_pos + 9; i < len; ++i)
		fqd->site_flag[i] = 1;

	return NO_ERROR;
} /* trim_to_barcode */

/**
 * Make options object.
 *
 * @param opt	pointer to options object
 * @return	error status
 */
int make_options(options **opt)
{
	options *op;

	*opt = malloc(sizeof **opt);
	if (!*opt)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "options object");

	op = *opt;

	op->fastq_file = NULL;
	op->in_primer = NULL;
	op->primer = op->cprimer = NULL;

	op->n_indices = 4;	/* 0-3 random bases inserted */

	op->n_init = 1;
	op->n_iter = 100;
	op->epsilon = 1e-2;
	op->tolerance = 1e-12;
	op->info = QUIET;
	op->estimate = ESTIMATE_ALL;

	op->trimmed_outfile = NULL;
	op->trim_outfile = NULL;

	return NO_ERROR;
} /* make_options */

/**
 * Parse command-line options.
 *
 * @param opt	pointer to options object
 * @param argc	number of command-line arguments
 * @param argv	the command lines arguments
 * @return	error status
 */
int parse_options(options *opt, int argc, const char **argv)
{
	int i, j;
	int err = NO_ERROR;
	char a;

	for (i = 1; i < argc; i++) {
		if (strlen(argv[i]) < 2)
			usage_error(argv, i, (void *)opt);
		j = 1;
		a = argv[i][j];
		while (a == '-' && ++j < (int) strlen(argv[i]))
			a = argv[i][j];
		switch(a) {
			case 'x':
				opt->n_indices = read_uint(argc, argv, ++i,
					(void *)opt);
				if (!opt->n_indices)
					opt->n_indices = 1;
				if (errno)
					goto CMDLINE_ERROR;
				break;
			case 'i':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				opt->n_init = read_uint(argc, argv, ++i,
					(void *)opt);
				if (errno)
					goto CMDLINE_ERROR;
				break;
			case 'f':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				opt->fastq_file = argv[++i];
				break;
			case 'o':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				if (argv[i][j+1] == 'f')
					opt->trimmed_outfile = argv[++i];
				else if (argv[i][j+1] == 't')
					opt->trim_outfile = argv[++i];
				else
					goto CMDLINE_ERROR;
				break;
			case 'p':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				opt->in_primer = argv[++i];
				break;
			case 's':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				opt->seed = read_uint(argc, argv, ++i,
					(void *)opt);
				srand(opt->seed);
				break;
			case 'c':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				if (!strcmp(argv[++i], "PI"))
					opt->estimate ^= ESTIMATE_PI;
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

	err = sync_options(opt);
	return err;

CMDLINE_ERROR:
	if (err == NO_ERROR) {
		err = INVALID_CMD_ARGUMENT;
		i--;
	}
	usage_error(argv, i, (void *)opt);
	return err;
} /* parse_options */

/**
 * Finish preparing the options object.
 *
 * @param opt	options pointer
 */
int sync_options(options *opt)
{
	unsigned int i;
	char const * const str = "GCCTTGCCACACGCTCAGNNNNNNNNNGTTGTAAYTTCTAGRTCCCCTCCTG";

	opt->primer_len = opt->in_primer ? strlen(opt->in_primer) : strlen(str);

	opt->primer = malloc(opt->primer_len * sizeof *opt->primer);
	opt->cprimer = malloc(opt->primer_len * sizeof *opt->cprimer);

	if (!opt->primer || !opt->cprimer)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"options::primer");

	for (i = 0; i < opt->primer_len; ++i)
		opt->primer[i] = (char_t) (opt->in_primer ? opt->in_primer[i] : str[i]);

//	strncpy(opt->primer, opt->in_primer ? opt->in_primer : str, opt->primer_len);
	memcpy(opt->cprimer, opt->primer, opt->primer_len);
//	strncpy(opt->cprimer, opt->primer, opt->primer_len);


	/* \ref options.primer contains raw chars that we need to encode. */

	for (i = 0; i < opt->primer_len; ++i)
		opt->primer[i] = nuc_to_iupac[opt->primer[i] - 'A'];

	/* count the number of distinct primers */
	opt->n_primers = 1;
	for (i = 0; i < opt->primer_len; ++i)
		/* non-barcode */
		if (popcnt[(int)opt->primer[i]] < NUM_NUCLEOTIDES)
			opt->n_primers *= popcnt[(int)opt->primer[i]];

	/* the number of possible hidden states: primer x position */
	opt->n_states = opt->n_indices * opt->n_primers + 1;

	return NO_ERROR;
} /* sync_options */

/**
 * Free options object.
 *
 * @param opt	 pointer to allocate options object
 */
void free_options(options *opt) {
	if (opt) free(opt);
} /* free_options */

/**
 * Create a data object.
 *
 * @param dat	data object pointer
 * @param opt	filled-in options object
 * @return	error status
 */
int make_data(data **dat, options *opt)
{
	UNUSED(opt);
	data *dp;

	*dat = malloc(sizeof **dat);

	if (!*dat)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data object");

	dp = *dat;

	dp->fqd = NULL;

	return NO_ERROR;
} /* make_data */

/**
 * Finish making data object.
 *
 * @param dat	pointer to previously allocated data object
 * @param opt	pointer to filled-in options object
 * @return	error status
 */
int finish_make_data(data *dat, options *opt)
{
	//int fxn_debug = ABSOLUTE_SILENCE;
	UNUSED(opt);
	unsigned int i;
	char_t *rptr = dat->fqd->reads;

	dat->index = malloc(dat->fqd->n_reads * sizeof *dat->index);
	dat->best_index = malloc(dat->fqd->n_reads * sizeof *dat->index);
	dat->optimal_index = malloc(dat->fqd->n_reads * sizeof *dat->index);

	if (!dat->index || !dat->best_index || !dat->optimal_index)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "index arrays");

	dat->reads = malloc(dat->fqd->n_reads * sizeof *dat->reads);

	if (!dat->reads)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "read pointers");

	for (i = 0; i < dat->fqd->n_reads; ++i) {
		dat->reads[i] = rptr;
		rptr += read_length(dat->fqd, i);
	}

	return NO_ERROR;
} /* finish_make_data */

/**
 * Free data object.
 *
 * @param dat	data object
 * @param opt	options object
 */
void free_data(data *dat, options *opt)
{
	UNUSED(opt);
	if (!dat) return;
	if (dat->fqd) free_fastq(dat->fqd);
	free(dat);
} /* free_data */

/**
 * Create model object.
 *
 * @param mod	pointer to empty model object
 * @param dat	pointer to data object
 * @param opt	pointer to options object
 * @return	error status
 */
int make_model(model **mod, data *dat, options *opt)
{
	model *mp;

	*mod = malloc(sizeof **mod);
	if (!*mod)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "model object");

	mp = *mod;

	mp->pi = malloc((opt->n_indices + 1) * sizeof *mp->pi);
	mp->npi = malloc((opt->n_indices + 1) * sizeof *mp->npi);
	mp->best_pi = malloc((opt->n_indices + 1) * sizeof *mp->best_pi);

	if (!mp->pi || !mp->npi || !mp->best_pi)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "model::pi");

	mp->gamma = malloc(NUM_NUCLEOTIDES * NUM_NUCLEOTIDES * sizeof *mp->gamma);
	mp->ngamma = malloc(NUM_NUCLEOTIDES * NUM_NUCLEOTIDES * sizeof *mp->ngamma);
	mp->best_gamma = malloc(NUM_NUCLEOTIDES * NUM_NUCLEOTIDES
		* sizeof *mp->best_gamma);

	if (!mp->gamma || !mp->ngamma || !mp->best_gamma)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "model::gamma");

	mp->delta = malloc((opt->primer_len + opt->n_indices - 1) * sizeof *mp->delta);
	mp->ndelta = malloc((opt->primer_len + opt->n_indices - 1) * sizeof *mp->ndelta);
	mp->best_delta = malloc((opt->primer_len + opt->n_indices - 1) * sizeof *mp->delta);

	if (!mp->delta || !mp->ndelta || !mp->best_delta)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "model::delta");

	mp->eta = malloc(NUM_NUCLEOTIDES * sizeof *mp->eta);
	mp->neta = malloc(NUM_NUCLEOTIDES * sizeof *mp->neta);
	mp->best_eta = malloc(NUM_NUCLEOTIDES * sizeof *mp->best_eta);

	if (!mp->eta || !mp->neta || !mp->best_eta)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "model::eta");

	mp->beta = malloc(NUM_NUCLEOTIDES * sizeof *mp->beta);
	mp->nbeta = malloc(NUM_NUCLEOTIDES * sizeof *mp->nbeta);
	mp->best_beta = malloc(NUM_NUCLEOTIDES * sizeof *mp->best_beta);

	if (!mp->beta || !mp->nbeta || !mp->best_beta)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "model::beta");
/*
	mp->n_quality = dat->fqd->n_quality;
	mp->lambda0 = malloc(dat->fqd->n_length * mp->n_quality
		* sizeof *mp->nlambda0);
	mp->nlambda0 = malloc(dat->fqd->n_length * mp->n_quality
		* sizeof *mp->nlambda0);
	mp->best_lambda0 = malloc(dat->fqd->n_length * mp->n_quality
		* sizeof *mp->best_lambda0);
	mp->lambda1 = malloc(dat->fqd->n_length * mp->n_quality
		* sizeof *mp->lambda1);
	mp->nlambda1 = malloc(dat->fqd->n_length * mp->n_quality
		* sizeof *mp->nlambda1);
	mp->best_lambda1 = malloc(dat->fqd->n_length * mp->n_quality
		* sizeof *mp->best_lambda1);

	if (!mp->lambda0 || !mp->nlambda0 || !mp->best_lambda0 || !mp->lambda1
		|| !mp->nlambda1 || !mp->best_lambda1)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "mpel::lambda");
*/
	mp->eik = malloc(dat->fqd->n_reads * opt->n_states * sizeof *mp->eik);

	if (mp->eik == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "model::eik");

	return NO_ERROR;
} /* make_model */

/**
 * Initialize the model.
 *
 * @param mod	pointer to allocated model object
 * @param dat	pointer to data object
 * @param opt	pointer to options object
 * @return	error status
 */
int initialize_model(model *mod, data *dat, options *opt)
{

	/* mini E-step: randomly assign hidden states */
	for (unsigned int i = 0; i < dat->fqd->n_reads; ++i)
		for (size_t k = 0; k < opt->n_states; ++k)
			mod->eik[i * opt->n_states + k] = 0.;

	/* randomize hidden state assignments */
	for (unsigned int i = 0; i < dat->fqd->n_reads; ++i) {
		size_t j = (size_t) (rand() / (RAND_MAX + 1.) * opt->n_states);
		mod->eik[i * opt->n_states + j] = 1.;
	}

	if (!(opt->estimate & ESTIMATE_PI))
		for (unsigned int i = 0; i <= opt->n_indices; ++i)
			mod->npi[i] = 1. / (opt->n_indices + 1);

	m_step(dat, mod, opt);

	param_update(mod, opt, DEFAULT);

	fprint_doubles(stderr, mod->pi, opt->n_indices + 1, 3, 1);
	fprint_doubles(stderr, mod->delta, opt->n_indices + opt->primer_len - 1, 3, 1);
	print_gamma(stderr, mod->gamma, NUM_NUCLEOTIDES);

	return NO_ERROR;
} /* initialize_model */

/**
 * Print gamma parameters to file.
 *
 * @param fp	open file pointer
 * @param gamma	transition matrix parameter
 * @param n	number of stats
 */
void print_gamma(FILE *fp, double *gamma, size_t n)
{
	size_t i, j;

	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			if (gamma[i*n + j] < 1e-2)
				fprintf(fp, " %8.2e", gamma[i*n + j]);
			else
				fprintf(fp, " %8.3f", gamma[i*n + j]);
		}
		fprintf(fp, "\n");
	}
} /* print_gamma */

/**
 * Primt lambda parameter to file.
 *
 * @param fp		open file pointer
 * @param lambda	parameters
 * @param n		
 * @param l		
 */
void print_lambda(FILE *fp, double *lambda, size_t n, size_t l)
{
	size_t i, j;

	for (i = 0; i < n; ++i) {
		fprintf(fp, "%3zu", i);
		for (j = 0; j < l; ++j) {
			if (lambda[i*l + j] < 1e-2)
				fprintf(fp, " %8.2e", lambda[i*l + j]);
			else
				fprintf(fp, " %8.3f", lambda[i*l + j]);
		}
		fprintf(fp, "\n");
	}
} /* print_lambda */

/**
 * Free model object.
 *
 * @param mod	pointer to allocate model object
 */
void free_model(model *mod)
{
	if (!mod) return;
	if (mod->pi) free(mod->pi);
	if (mod->npi) free(mod->npi);
	if (mod->best_pi) free(mod->best_pi);
	if (mod->delta) free(mod->delta);
	if (mod->ndelta) free(mod->ndelta);
	if (mod->best_delta) free(mod->best_delta);
	if (mod->gamma) free(mod->gamma);
	if (mod->ngamma) free(mod->ngamma);
	if (mod->best_gamma) free(mod->best_gamma);
	if (mod->eta) free(mod->eta);
	if (mod->neta) free(mod->neta);
	if (mod->best_eta) free(mod->best_eta);
	free(mod);
} /* free_model */

/**
 * Display trim-related error messages.
 *
 * @param err_no	number of error to print
 * @return		character string
 */
char * trim_error_message(int err_no)
{
	if (err_no == TRIM_AMBIGUOUS_NUCLEOTIDE)
		return "ambiguous nucleotide";
	else if (err_no)
		return "unknown error";
	else
		return "";
} /* trim_error_message */

/**
 * Print command-line usage.
 *
 * @param fp		Open file pointer
 * @param cmdname	Name of executable command
 * @param obj		Void pointer to additional information
 */
void fprint_usage(FILE *fp, const char *cmdname, void *obj)
{
	options *opt = (options *) obj;
	size_t start = strlen(cmdname) - 1;

	while (cmdname[start] != '/' && start) start--;
	if (cmdname[start] == '/') start++;

	for (size_t i = start; i < strlen(cmdname); ++i)
		fputc(toupper(cmdname[i]), fp);
	fprintf(fp, "(%d)\n", 1);
	fprintf(fp, "\nNAME\n\t%s - trim amplicon sequences\n",
		&cmdname[start]);
	fprintf(fp, "\nSYNOPSIS\n\t%s [-s <sulong> -i <iuint>] -x <xuint>"
		"-f <fstr>...\n", &cmdname[start]);
	fprintf(fp, "\nDESCRIPTION\n\t%s identifies known primer sequence in "
		"NGS reads from fastq file <fstr> among <xuint> possible "
		"positions, randomly initialize <iuint> times after setting "
		"random number seed <sulong>.\n", &cmdname[start]);
	fprintf(fp, "\nOPTIONS\n");
	fprintf(fp, "\t-x <xuint>\n\t\tMaximum number of random nucleotides "
		"possible at 5' end of read (%d).\n", opt->n_indices);
	fprintf(fp, "\t-c <cstr>\n\t\tTurn off estimation of variable "
		"(choices: PI) [DOES NOT WORK]\n");
	fprintf(fp, "\t-i <iuint>\n\t\tSet desired number of initializations "
		"[OPTIONAL; DEFAULT: %u].\n", opt->n_init);
	fprintf(fp, "\t-f <fstr>\n\t\tSet the fastq input file.\n");
	fprintf(fp, "\t-of <ostr>\n\t\tSet the fastq output file for trimmed reads.\n");
	fprintf(fp, "\t-ot <ostr>\n\t\tSet the fastq output file for trimmings.\n");
	fprintf(fp, "\t-p <pstr>\n\t\tSet the primer.\n");
	fprintf(fp, "\t-s <sulong>\n\t\tSet random number seed. [OPTIONAL]\n");
	fprintf(fp, "\t-h\n\t\tThis help.\n");
	fprintf(fp, "\n");
	for (size_t i = start; i < strlen(cmdname); ++i)
		fputc(toupper(cmdname[i]), fp);
	fprintf(fp, "(%d)\n", 1);
} /* fprint_usage */
