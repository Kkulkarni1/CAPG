/**
 * @file capg.c
 * @author K. S. Dorman
 *
 * This file contains the code for genotyping allotetraploids.
 *
 * To compile:
 make capg
 *
 * PLAN (searchable keyword) for whole-genome-sequencing genotyping:
 * - no need to filter by reference
 * - index reads by reference AND position in reference
 * - a good first test: have it genotype all targets for a single genotype
 */


/* subgenome A                         E = read alignment extent
 * \                                   T = target region,        A = aligned region
 *  \ read 1
 *   \\____           __/              E    T    A                           A     E T
 *    \__________________/             |____|____|_________..._______________|_____|_|
 *     __________________    =======> _____________________...__________________________
 *    /      _____     _ \            |          | |                         |   |     |
 *   /                  \ \           T          A E                         A   T     E
 *                       \
 *   subgenome B
 */

#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <unistd.h>
#include <zlib.h>

#include "capg.h"
#include "io.h"
#include "cmdline.h"
#include "array.h"
#include "order.h"
#include "pick_reads.h"
#include "nmath.h"

/**
 * State of a read nucleotide, including pointer to containing alignment,
 * probability of error, position in alignment, true nucleotide, observed
 * read nucleotide, and quality score.
 * http://www.catb.org/esr/structure-packing/
 */


double ll_align(sam_entry *se, unsigned int i, unsigned char *ref, mlogit_stuff *vptr, unsigned char *show, size_t start_rf, int debug);
int update_vcf(options *opt, int *fail[N_FILES], size_t *hpos[N_FILES]);
uint64_t get_haplotype_id(sam_entry *se, size_t *haplotype, unsigned int n_segregating);
int test_equal_homolog_coverage(merge_hash *mh, ref_info *rfi, double **ll, char_t ref_base[N_FILES], char *covers, xy_t *obs_nuc, qual_t *obs_q, unsigned int *obs_rpos, int g_max[N_FILES], xy_t nuc1, xy_t nuc2, double pvals[N_FILES]);
unsigned int ns_key_len = 3 * sizeof(data_t) + sizeof(unsigned int);


/**
 * Test equal coverage of homologous chromosomes.
 *
 * @param mh		merged hash of aligned reads
 * @param ref_info	reference information
 * @param ll		log likelihood of reads aligned to each subgenome
 * @param ref_base	reference bases
 * @param covers	indicate if retained reads cover the site
 * @param obs_nuc	observed nucleotides in reads covering the site
 * @param obs_q		observed quality scores of nucleotides covering the site
 * @param obs_rpos	other information for error model
 * @param g_max		the genotype call
 * @param nuc1		the major allele
 * @param nuc2		the minor allele
 * @param debug_level	inherited debugging level
 * @param pvals		up to two pvals
 * @return		error status
 */
int test_equal_homolog_coverage(merge_hash *mh, ref_info *rfi, double **all,
	char_t ref_base[N_FILES], char *covers, xy_t *obs_nuc, qual_t *obs_q,
	unsigned int *obs_rpos, int g_max[N_FILES], xy_t nuc1, xy_t nuc2,
	double pvals[N_FILES])
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_II;//
	mlogit_stuff mls = {NULL, 0};
	double epsilon = 1e-6;
	double gamma[N_FILES] = {g_max[0]/2., g_max[1]/2.};
	double eta = 0.5;	/* assumes N_FILES == 2 */
	double new_gamma[N_FILES];
	double new_eta, gamma1, gamma2, eta_h1;
	double lpi[2*N_FILES];
	double ll1;
	double ll = -INFINITY, pll;
	double pp[2*N_FILES], cll[2*N_FILES];
	double log_half = log(0.5);
	size_t n_cover, n_read;
	unsigned int iter = 0, max_iter = 100;

	/* estimate log likelihood under H1 */
	do {
		/* initialize log likelihood */
		pll = ll;
		ll = 0;

		/* compute current mixing proportions */
		lpi[0] = lpi[1] = log(eta);
		if (g_max[0] == 1) {
			lpi[0] += log(gamma[0]);
			lpi[1] += log(1 - gamma[0]);
		}
		lpi[2] = lpi[3] = log(1 - eta);
		if (g_max[1] == 1) {
			lpi[2] += log(gamma[1]);
			lpi[3] += log(1 - gamma[1]);
		}

		new_eta = 0;
		new_gamma[0] = 0;
		new_gamma[1] = 0;	/* N_FILES == 2 */
		n_cover = 0;
		n_read = 0;

		double total1 = 0, total2 = 0;

		for (merge_hash *me = mh; me != NULL; me = me->hh.next) {

			/* skip excluded reads */
			if (me->exclude)
				continue;

			/* skip read not covering site */
			if (!covers[n_read]) {
				++n_read;
				continue;
			}

			mls.pos = obs_rpos[n_cover];

			double max = 0;

			for (int i = 0; i < 2*N_FILES; ++i) {

				xy_t tmp = obs_nuc[n_cover];
				if (i>=2 && rfi->strand_B)
					tmp = xy_to_rc[tmp];
				cll[i] = lpi[i] + all[i/2][n_read]
					- sub_prob_given_q_with_encoding(
						ref_base[i/2],				/* homozygous ref base */
						tmp, IUPAC_ENCODING, XY_ENCODING, obs_q[n_cover], 1, (void *)&mls)
					+ sub_prob_given_q_with_encoding(
						 !g_max[i/2] ? xy_to_iupac[nuc1]	/* new genotype */
						: g_max[i/2] == 2 ? xy_to_iupac[nuc2]
						: !(i%2) ? xy_to_iupac[nuc1] : xy_to_iupac[nuc2],
						obs_nuc[n_cover], IUPAC_ENCODING, XY_ENCODING, obs_q[n_cover], 1, (void *)&mls);
debug_msg(fxn_debug >= DEBUG_III, fxn_debug, "Read %zu: genotype %c -> %c emitting %c/%c has ll %f\n", n_read, iupac_to_char[ref_base[i/2]], xy_to_char[!g_max[i/2] ? nuc1 : g_max[i/2] == 2 ? nuc2 : !(i%2) ? nuc1: nuc2], xy_to_char[tmp], xy_to_char[obs_nuc[n_cover]], cll[i]);
				if (cll[i] > max)
					max = cll[i];
			}

			double sum = 0;
			for (int i = 0; i < 2*N_FILES; ++i) {
				pp[i] = exp(cll[i] - max);
				sum += pp[i];
			}
			for (int i = 0; i < 2*N_FILES; ++i)
				pp[i] /= sum;
			new_eta += pp[0] + pp[1];
			if (g_max[0] == 1) {
				new_gamma[0] += pp[0];
				total1 += pp[0] + pp[1];
			}
			if (g_max[1] == 1) {
				new_gamma[1] += pp[2];
				total2 += pp[2] + pp[3];
			}

			ll += log(sum) + max;

			++n_cover;
			++n_read;
		}
		eta = new_eta / n_cover;
		if (g_max[0] == 1)
			gamma[0] = new_gamma[0] / total1;
		if (g_max[1] == 1)
			gamma[1] = new_gamma[1] / total2;
		debug_msg(fxn_debug >= DEBUG_II, fxn_debug, "eta = %f, gamma1 = %f, gamma2 = %f, ll = %f, rel = %e\n", eta, gamma[0], gamma[1], ll, (pll - ll) / ll);
	} while (iter++ < max_iter && (ll - pll) > -ll * epsilon);

	if (g_max[0] == 1)
		gamma1 = gamma[0];
	else if (!g_max[0])
		gamma1 = 1;
	else
		gamma1 = 0;
	if (g_max[1] == 1)
		gamma2 = gamma[1];
	else if (!g_max[1])
		gamma2 = 1;
	else gamma2 = 0;
	eta_h1 = eta;

	ll1 = ll;
	ll = -INFINITY;

	for (int j = 0; j < N_FILES; ++j) {
		if (g_max[j] != 1)
			continue;

		do {
			/* initialize log likelihood */
			pll = ll;
			ll = 0;
	
			/* compute current mixing proportions */
			lpi[0] = lpi[1] = log(eta);
			if (g_max[0] == 1) {
				if (!j) {
					lpi[0] += log_half;
					lpi[1] += log_half;
				} else {
					lpi[0] += log(gamma[0]);
					lpi[1] += log(1 - gamma[0]);
					new_gamma[0] = 0;
				}
			}
			lpi[2] = lpi[3] = log(1 - eta);
			if (g_max[1] == 1) {
				if (j) {
					lpi[2] += log_half;
					lpi[3] += log_half;
				} else {
					lpi[2] += log(gamma[1]);
					lpi[3] += log(1 - gamma[1]);
					new_gamma[1] = 0;
				}
			}
	
			new_eta = 0;
			n_cover = 0;
			n_read = 0;
	
			double total1 = 0, total2 = 0;
	
			for (merge_hash *me = mh; me != NULL; me = me->hh.next) {
	
				/* skip excluded reads */
				if (me->exclude)
					continue;
	
				/* skip read not covering site */
				if (!covers[n_read]) {
					++n_read;
					continue;
				}
	
				mls.pos = obs_rpos[n_cover];
	
				double max = 0;
	
				for (int i = 0; i < 2*N_FILES; ++i) {
	
					xy_t tmp = obs_nuc[n_cover];
					if (i>=2 && rfi->strand_B)
						tmp = xy_to_rc[tmp];
					cll[i] = lpi[i] + all[i/2][n_read]
						- sub_prob_given_q_with_encoding(
							ref_base[i/2],				/* homozygous ref base */
							tmp, IUPAC_ENCODING, XY_ENCODING, obs_q[n_cover], 1, (void *)&mls)
						+ sub_prob_given_q_with_encoding(
							 !g_max[i/2] ? xy_to_iupac[nuc1]	/* new genotype */
							: g_max[i/2] == 2 ? xy_to_iupac[nuc2]
							: !(i%2) ? xy_to_iupac[nuc1] : xy_to_iupac[nuc2],
							obs_nuc[n_cover], IUPAC_ENCODING, XY_ENCODING, obs_q[n_cover], 1, (void *)&mls);
		
debug_msg(fxn_debug >= DEBUG_III, fxn_debug, "Read %zu: genotype %c -> %c emitting %c/%c has ll %f\n", n_read, iupac_to_char[ref_base[i/2]], xy_to_char[!g_max[i/2] ? nuc1 : g_max[i/2] == 2 ? nuc2 : !(i%2) ? nuc1: nuc2], xy_to_char[tmp], xy_to_char[obs_nuc[n_cover]], cll[i]);
					if (cll[i] > max)
						max = cll[i];
				}
	
				double sum = 0;
				for (int i = 0; i < 2*N_FILES; ++i) {
					pp[i] = exp(cll[i] - max);
					sum += pp[i];
				}
				for (int i = 0; i < 2*N_FILES; ++i)
					pp[i] /= sum;
				new_eta += pp[0] + pp[1];
				if (j && g_max[0] == 1) {
					new_gamma[0] += pp[0];
					total1 += pp[0] + pp[1];
				}
				if (!j && g_max[1] == 1) {
					new_gamma[1] += pp[2];
					total2 += pp[2] + pp[3];
				}
	
				ll += log(sum) + max;
	
				++n_cover;
				++n_read;
			}
			eta = new_eta / n_cover;
			if (j && g_max[0] == 1)
				gamma[0] = new_gamma[0] / total1;
			if (!j && g_max[1] == 1)
				gamma[1] = new_gamma[1] / total2;
			debug_msg(fxn_debug >= DEBUG_II, fxn_debug, "eta = %f, gamma1 = %f, gamma2 = %f, ll = %f, rel = %e\n", eta, (g_max[0] && j) ? gamma[0] : 0.5, (g_max[1] && !j) ? gamma[1] : 0.5, ll, (pll - ll) / ll);
		} while (iter++ < max_iter && (ll - pll) > -ll * epsilon);
	
		double lrt = 2 * (ll1 - ll);
		pvals[j] = pchisq(lrt, 1, 0, 0);
		mmessage(INFO_MSG, NO_ERROR,
			"Equal coverage test (subgenome %u): eta = %f, gamma1 = %f, gamma2 = %f vs. eta = %f, gamma1 = %f, gamma2 = %f; lrt = %f (%f %f); pval = %e\n", j,
			eta_h1, gamma1, gamma2, eta, j && g_max[0] == 1 ? gamma[0] : g_max[0] / 2., !j && g_max[1] == 1 ? gamma[1] : g_max[1] / 2., 2 * (ll1 - ll), ll1, ll, pvals[j]);
	}

	return 0;
} /* test_equal_homolog_coverage */

/**
 * Update vcf files after post-hoc filters.  The sites that have failed the
 * filter(s) are indicated in fail, currently only segregating sites.
 * To identify the offending sites, the 0-based position is in hpos.
 *
 * @param opt	options object pointer
 * @param fail	failed segregating sites for each subgenome
 * @param hpos	position of segregating sites
 * @error	error status
 */
int update_vcf(options *opt, int *fail[N_FILES],
				size_t *hpos[N_FILES])
{
	const char *tmpfile_template = "tmp_vcfXXXXXX";

	for (int i = 0; i < N_FILES; ++i) {
		unsigned int nsegregating = 0;
		int fpd;
		char *tmpfile = NULL;
		FILE *fpr = NULL;
		FILE *fpt = NULL;

		if (!opt->vcf_files[i])
			continue;

		if (!fail[i])
			continue;

		fpr = fopen(opt->vcf_files[i], "r");

		if (!fpr)
			return mmessage(ERROR_MSG, FILE_OPEN_ERROR,
							opt->vcf_files[i]);

		tmpfile = malloc((strlen(tmpfile_template) + 1)
							* sizeof *tmpfile);

		if (!tmpfile)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
								"tmpfile");

		strcpy(tmpfile, tmpfile_template);

		fpd = mkstemp(tmpfile);
		fpt = fdopen(fpd, "w");

		if (!fpt)
			return mmessage(ERROR_MSG, FILE_OPEN_ERROR, tmpfile);

		char c;
		while (!feof(fpr)) {
			unsigned int pos;
			int already_filtered = 0;
			const char *pass = "PASS";

			c = fgetc(fpr);

			if (feof(fpr))
				break;

			/* copy header */
			while (!feof(fpr) && c == '#') {
				fputc(c, fpt);	/* leading # */
				while (!feof(fpr) && (c = fgetc(fpr)) != '\n')
					fputc(c, fpt);
				fputc(c, fpt);	/* newline */
				c = fgetc(fpr);
			}

			/* skip first tab-separated columns */
			fputc(c, fpt);	/* first char or tab */
			while (!feof(fpr) && (c = fgetc(fpr)) != '\t')
				fputc(c, fpt);

			if (fscanf(fpr, "%u", &pos) != 1)
				return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
							opt->vcf_files[i]);

			fprintf(fpt, "\t%d", pos);
			c = fgetc(fpr);	/* tab */

			/* skip next four tab-separated columns */
			for (int i = 0; i < 4; ++i) {
				fputc(c, fpt);	/* first char or tab */
				while (!feof(fpr) && (c = fgetc(fpr)) != '\t')
					fputc(c, fpt);
			}
			fputc(c, fpt);	/* tab */

			c = fgetc(fpr);	/* read first char in FILTER */
			if (c == 'P') {	/* maybe PASS */
				unsigned int cnt = 0;

				/* skip rest of PASS */
				do {
					c = fgetc(fpr);
					++cnt;
				} while (!feof(fpr) && cnt < strlen("PASS") && pass[cnt] == c);

				/* nope, something else: output it */
				if (c != '\t') {
					already_filtered = 1;
					fprintf(fpt, "PASS%c", c);
					while (!feof(fpr) && (c = fgetc(fpr)) != '\t')
						fputc(c, fpt);
				}
			} else if (c != '.') {	/* other filters */
				already_filtered = 1;
				fputc(c, fpt);
				while (!feof(fpr) && (c = fgetc(fpr)) != '\t')
					fputc(c, fpt);
			} else {
				c = fgetc(fpr);	/* tab */
			}

			/* one of the segregating sites */
			if (pos - 1 == hpos[i][nsegregating]) {

				/* this one failed: currently just 1 test */
				if (fail[i][nsegregating]) {
					if (already_filtered)
						fputc(';', fpt);

					fprintf(fpt, "sc5");	/* add failed filter */
				
				/* otherwise it passed */
				} else if (!already_filtered) {
					fprintf(fpt, "PASS");
				}
				++nsegregating;

			/* non-segregating site, not filtered by us */
			} else if (!already_filtered) {
				fprintf(fpt, "PASS");
			}

			fputc(c, fpt);		/* tab after FILTER */
			while (!feof(fpr) && (c = fgetc(fpr)) != '\n')
				fputc(c, fpt);	/* rest of line */
			if (!feof(fpr))
				fputc(c, fpt);	/* newline */
		}

		fclose(fpr);

		char *command = malloc((strlen("mv") + strlen(tmpfile_template)
			+ strlen(opt->vcf_files[i]) + 3) * sizeof *command);

		if (!command) {
			free(tmpfile);
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "command");
		}

		sprintf(command, "cp %s %s", tmpfile, opt->vcf_files[i]);

		fclose(fpt);
		system(command);

		sprintf(command, "rm %s", tmpfile);
		system(command);

		free(tmpfile);
		free(command);
	}

	return NO_ERROR;
} /* update_vcf */

/**
 * Read error parameter file from mlogit fit.
 * 64 C 1226 C 0.997768688734120945
 *
 * @param param_file	name of parameter file
 * @return		hash of read probabilities
 */
nuc_state *read_param_file(char const *param_file)
{
	FILE *fp = fopen(param_file, "r");

	if (!fp) {
		fprintf(stderr, "Could not open parameter file '%s'.\n",
			param_file);
		return NULL;
	}

	/* count number of lines excluding header */
	unsigned int n_states = 0;
	char c = fgetc(fp), cprev = 0;
	int header = (c > 57);
	while (!feof(fp)) {
		if (c == '\n')
			++n_states;
		cprev = c;
		c = fgetc(fp);
	}
	if (cprev != '\n')
		++n_states;
	if (header)
		--n_states;

	/* allocate data */
	nuc_state *ns_hash = NULL;
	nuc_state *itm = NULL;
	nuc_state *ns = malloc(n_states * sizeof *ns);
	if (!ns) {
		fprintf(stderr, "Could not allocate memory.\n");
		return NULL;
	}

	rewind(fp);

	/* skip first line */
	if (header)
		while(!feof(fp) && (c = fgetc(fp)) != '\n');

	int i = 0;
	do {
		if (fscanf(fp, "%hhu %c %u %c %lf", &ns[i].quality, &ns[i].true_nuc,
			   &ns[i].pos, &ns[i].read_nuc, &ns[i].prob) != 5) {
			if (feof(fp))
				break;
			fprintf(stderr, "Error reading file '%s', entry %d.\n",
				param_file, i);
			free(ns);
			return NULL;
		}
		--ns[i].pos;	/* 0 index */
		//		ns[i].quality -= MIN_ASCII_QUALITY_SCORE;
		ns[i].true_nuc = (ns[i].true_nuc >> 1) & 3;
		ns[i].read_nuc = (ns[i].read_nuc >> 1) & 3;
		if (ns[i].prob == 0)
			ns[i].prob = DBL_MIN;
		//fprintf(stderr, "Read state (%c %c %u %u) to hash %p.\n", xy_to_char[ns[i].true_nuc], xy_to_char[ns[i].read_nuc], ns[i].quality, ns[i].pos, (void *)ns_hash);
		HASH_FIND(hh, ns_hash, &ns[i].pos, ns_key_len, itm);
		if (!itm) {
			//fprintf(stderr, "Adding state (%c %c %u %u) to hash %p.\n", xy_to_char[ns[i].true_nuc], xy_to_char[ns[i].read_nuc], ns[i].quality, ns[i].pos, (void *)ns_hash);
			HASH_ADD(hh, ns_hash, pos, ns_key_len, &ns[i]);
			++i;
		}
	} while (!feof(fp));

	nuc_state *nns = realloc(ns, i * sizeof *ns);
	if (!nns) {
		fprintf(stderr, "ERROR: Could not reallocate shrunken hash.\n");
		return NULL;
	}
	ns = nns;

	fclose(fp);

	return ns_hash;
} /* read_param_file */

///**
// * Find medium from vector x[start:end_inclusive].
// *
// * @param x		vector of numeric values
// * @param start		first index to include in dataset
// * @param end_inclusive	last index to include in dataset
// * @return		median
// */
//static double median(double *x, int start, int end_inclusive)
//{
//	int size = end_inclusive - start + 1;
//	if (size <= 0) {
//		printf("Array slice cannot be empty\n");
//		exit(1);
//	}
//	int m = start + size / 2;
//	if (size % 2) return x[m];
//	return (x[m - 1] + x[m]) / 2.0;
//}/* median */

///**
// * Find mimimum, lower quartile, median, upper quartile and maximum of
// * a set x_len of numeric values.
// *
// * @param[in]	x	vector of numeric values
// * @param[out]	result	length-5 vector of result statistics
// * @param[in]	x_len	length of vector x
// * @return		error status
// */
//static int fivenum(double *x, double *result, int x_len)
//{
//	int i, m, lower_end;
//
//	for (i = 0; i < x_len; i++) {
//		if (x[i] != x[i]) {
//			printf("Unable to deal with arrays containing NaN\n\n");
//			return 1;
//		}
//	}
//	qsort(x, x_len, sizeof(double), double_compare);
//	result[0] = x[0];
//	result[2] = median(x, 0, x_len - 1);
//	result[4] = x[x_len - 1];
//	m = x_len / 2;
//	lower_end = (x_len % 2) ? m : m - 1;
//	result[1] = median(x, 0, lower_end);
//	result[3] = median(x, m, x_len - 1);
//	return NO_ERROR;
//}/* fivenum */

/**
 * Probability of substituting true nucleotide s for read nucleotide r
 * given covariate quality score q.
 *
 * typedef double (*sub_prob_fxn)(data_t, data_t, data_t, void *);
 *
 * @param[in]	s	true nucleotide
 * @param[in]	r	read nucleotide
 * @param[in]	q	quality score
 * @param[in]	vptr	additional data (mlogit hash)
 * @return		probability
 */
double mlogit_sub_prob(data_t s, data_t r, data_t q, void *vptr)
{
	mlogit_stuff *ms = (mlogit_stuff *) vptr;
	nuc_state *ns = ms->ns;
	nuc_state c_ns, *itm;

	c_ns.true_nuc = s;
	c_ns.read_nuc = r;
	c_ns.quality = q;
	c_ns.pos = ms->pos;
	HASH_FIND(hh, ns, &c_ns.pos, ns_key_len, itm);
	if (!itm) {
		//fprintf(stderr, "ERROR: Could not find nucleotide state (%c %c %u %u).\n", xy_to_char[c_ns.true_nuc], xy_to_char[c_ns.read_nuc], c_ns.quality, c_ns.pos);
		return 1;
	}
	return itm->prob;
}/* mlogit_sub_prob */

/**
 * Log probability of substitution (see mlogit_sub_prob).
 */
double mlogit_sub_lprob(data_t s, data_t r, data_t q, void *vptr)
{
	mlogit_stuff *ms = (mlogit_stuff *) vptr;
	nuc_state *ns = ms->ns;
	nuc_state c_ns, *itm;

	c_ns.true_nuc = s;
	c_ns.read_nuc = r;
	c_ns.quality = q;
	c_ns.pos = ms->pos;
	HASH_FIND(hh, ns, &c_ns.pos, ns_key_len, itm);
	if (!itm) {
		//		fprintf(stderr, "ERROR: Could not find nucleotide state (%c %c %u %u).\n", xy_to_char[c_ns.true_nuc], xy_to_char[c_ns.read_nuc], c_ns.quality, c_ns.pos);
		//exit(0);
		return 0;
	}
	//	fprintf(stderr, "itm (%c %c %u %u) -> prob: %e\n", xy_to_char[c_ns.true_nuc], xy_to_char[c_ns.read_nuc], c_ns.quality, c_ns.pos, itm->prob);
	return log(itm->prob);
}/* mlogit_sub_lprob */

/* roshan's research */
int main(int argc, const char *argv[])
{
	int debug_level = QUIET;//ABSOLUTE_SILENCE;//MINIMAL;//DEBUG_I;//
	int err = NO_ERROR;
	int karin_version = 1;
	fastq_data *fds[N_FILES] = {NULL, NULL};//, NULL, NULL};
	unsigned int fs_index[N_FILES] = {0, 0};//, 0, 0}; // index of alignment to accept from subgenomic alignment file: currently assume there is one alignment per target, so these are not changed
	sam *sds[N_FILES] = {NULL, NULL};
	fastq_options fop = {.read_encoding = IUPAC_ENCODING, .read_names = 1};
	sam_hash *by_name[N_FILES] = {NULL, NULL};//, NULL, NULL};
	FILE *fp = NULL;
	mlogit_stuff mls = {NULL, 0};
	options_rf opt_rf;
	ref_info *rf_info = NULL;
	options opt;

	/* process the command line */
	default_options(&opt);
	if (parse_options_capg(&opt, argc, argv))
		exit(mmessage(ERROR_MSG, INVALID_CMDLINE, ""));

	/* default reference options */
	default_options_rf(&opt_rf);
	opt_rf.sam_file = opt.sam_file;

	/* set names for extracted references */
	size_t rfile_len = strlen(opt.extracted_rf) + strlen(".fsa") + (int)(log10(N_FILES) + 1) + 1;
	for (int i = 0; i < N_FILES; ++i) {
		opt_rf.extracted_rf[i] = malloc(rfile_len * sizeof(opt_rf.extracted_rf[i]));
		sprintf(opt_rf.extracted_rf[i], "%s%d.fsa", opt.extracted_rf, i);
	}

	/* chromosome names */
	char *csome_names[N_FILES];
	for (unsigned int j = 0; j < N_FILES; ++j) {
		size_t i = 0;

		while (i < strlen(opt.ref_names[j]) && opt.ref_names[j][i++] != opt_rf.delim_ref);
		csome_names[j] = malloc(i-- * sizeof(*csome_names[j]));
		strncpy(csome_names[j], opt.ref_names[j], i);
		csome_names[j][i] = '\0';
		mmessage(INFO_MSG, NO_ERROR, "Chromosome %u is '%s'\n", j, csome_names[j]);
	}

	/* read sam file with target alignments and retain information about
	 * selected target
	 */
	make_targets_info(&opt_rf, &rf_info, opt.ref_names);
	
	/* read sam files of read alignments */
	for (unsigned int j = 0; j < N_FILES; ++j) {
		fp = fopen(opt.sbam_files[j], "r");
		if (!fp)
			exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR,
					  opt.sbam_files[j]));
		read_sam(fp, &sds[j], 1, 1);
		fclose(fp);
	}

	/* pick the reads aligned to homoeologous target, setting
	 * sam_entry::which_ref = 0 for selected reads.
	 */
	pickreads(rf_info, sds, (char const **)csome_names);

	/* index surviving reads to single target reference, while
	 * dropping additional reads according to user-supplied filters
	 */
	for (unsigned int j = 0; j < N_FILES; ++j) {
		
		hash_sam(sds[j], &by_name[j], HASH_REFERENCE, 1, //my_refs[j], rf_info->ref_sam->n_se,
			opt.drop_unmapped, opt.drop_secondary,
			opt.drop_soft_clipped, opt.drop_indel, opt.min_length,
			opt.max_length, opt.max_eerr);
		
		mmessage(INFO_MSG, NO_ERROR, "Number of %u alignments: %zu\n",
			 j, sds[j]->n_per_ref[0]);

		/* TODO,KSD This is slow, looping through EVERY read alignment. Why are we doing it?
		 * YD: Because forward and reverse reads may share the same name and would hash as multiple alignments of same read without modifying the read name to indicate strand.
		 */
		/*
		char strand;
		unsigned char found = 0;

		for (size_t i = 0; i < sds[j]->n_se; ++i) {
			sam_entry *se = &sds[j]->se[i];
			
			// skip unmapped
			if ((se->flag & (1 << 2)))
				continue;

			// and those already filtered, including not mapping to target
			if (se->exclude == 1)
				continue;

			if (!strcmp(opt.ref_names[j], se->ref_name)) {
				se->name_s = NULL;
				if ((se->flag & 16) == 0) {
					strand = '+';
					// flip the strand if A is aligned to reverse complement of B
					if (j && rf_info->strand_B == 1)
						strand = '-';
				} else {
					strand = '-';
					if (j && rf_info->strand_B == 1)
						strand = '+';
				}
	
				size_t length = strlen(se->name) + 1 + 1;
				se->name_s = malloc(length);
				sprintf(se->name_s, "%s%c", se->name, strand);
					
				my_refs[j] = se->which_ref;
				found = 1;
			}
		}
		if (!found)
			exit(mmessage(ERROR_MSG, INVALID_USER_INPUT, "no "
					  "reference '%s' in fasta file '%s'",
					  opt.ref_names[j], opt.sbam_files[j]));
		*/
	}

	/* merge reads for given reference */
	merge_hash *mh = NULL;

	size_t ref_indices[N_FILES] = {0,0};
	size_t n_read = hash_merge(&mh, N_FILES, sds, ref_indices);
	debug_msg(debug_level > ABSOLUTE_SILENCE, debug_level, "Number of reads"
		" aligned to target in AT LEAST one subgenome: %zu\n", n_read);
	
	if (opt_rf.fastq_file)
		output_selected_reads(opt_rf.fastq_file, sds, mh);
	
	size_t n_included_reads = 0;
	input *in = NULL;
	if (opt.ampliclust_command)
		make_input(&in);
	
	/* exclude reads not aligned to both genomes;
	 * identify extent of alignments on reference genome
	 */
	size_t start_reference[N_FILES];/* 0-based start of reference segment in each subgenome, inclusive */
	size_t end_reference[N_FILES];	/* 0-based end of reference segment in each subgenome, exclusive */
	for (unsigned int i = 0; i < N_FILES; ++i) {
		start_reference[i] = SIZE_MAX;
		end_reference[i] = 0;
	}

	n_read = 0;
	for (merge_hash *me = mh; me != NULL; me = me->hh.next) {
		
		if (me->nfiles != N_FILES) {
			me->exclude = 1;
			mmessage(INFO_MSG, NO_ERROR, "Read %s does not "
				 "align to all genomes (skipping).\n",
				me->indices[0]	// TODO: hack
				? sds[0]->se[me->indices[0][0]].name
				: sds[1]->se[me->indices[1][0]].name);
			continue;
		}
		
		/* force one alignment per sub-genome */
		for (unsigned int j = 0; j < N_FILES; ++j) {
			
			if (me->count[j] > 1) {	/* can happen if not proper pair */
				me->exclude = 1;
				mmessage(INFO_MSG, NO_ERROR,
					"Read %s aligns %dx in genome %u.\n",
					 sds[j]->se[me->indices[j][0]].name, me->count[j], j);
				break;
			}
		}
		if (me->exclude)
			continue;

		/* use this read to find extent of coverage */
		for (unsigned int j = 0; j < N_FILES; ++j) {
			sam_entry *se = &sds[j]->se[me->indices[j][0]];
			size_t rf_pos = se->pos - 1;

			if (start_reference[j] > rf_pos)
				start_reference[j] = rf_pos;
			rf_pos += se->cig->length_rf;
			if (end_reference[j] < rf_pos)
				end_reference[j] = rf_pos;
		}

		++n_read;
	}

	if (!n_read)
		exit(mmessage(ERROR_MSG, INTERNAL_ERROR,
				  "No read aligns to selected genome %s.\n"));

	for (unsigned int i = 0; i < N_FILES; ++i)
		mmessage(INFO_MSG, NO_ERROR, "Read extent on subgenome %u: [%u, %u)\n",
			 i, start_reference[i], end_reference[i]);

	/* extend reference regions to at least cover the user-selected
	 * target region
	 */
	if (start_reference[0] > rf_info->start_A)
		start_reference[0] = rf_info->start_A;
	if (end_reference[0] < rf_info->end_A)
		end_reference[0] = rf_info->end_A;
	if (start_reference[1] > rf_info->start_B)
		start_reference[1] = rf_info->start_B;
	if (end_reference[1] < rf_info->end_B)
		end_reference[1] = rf_info->end_B;

	for (unsigned int i = 0; i < N_FILES; ++i)
		mmessage(INFO_MSG, NO_ERROR, "Selected region on subgenome %u: [%u, %u)\n",
			 i, start_reference[i], end_reference[i]);

	/* TODO: handle >2 subgenomes in next two steps */

	/* record information about the reference alignment */
	match_pair(rf_info);

	/* match soft-clips so likelihood computed from same data in all subgenomes */
	match_soft_clipping(mh, N_FILES, sds, rf_info->strand_B);
	
	double **pp = malloc(N_FILES * sizeof(*pp));
	double **ll = malloc(N_FILES * sizeof(*ll));

	for (unsigned int j = 0; j < N_FILES; ++j) {

		pp[j] = malloc(n_read * sizeof(**pp));	/* overestimated */
		ll[j] = malloc(n_read * sizeof(**ll));	/* overestimated */
	}
	
	if (opt.param_file) {
		debug_msg(debug_level > ABSOLUTE_SILENCE, debug_level,
			  "Reading parameter file '%s'.\n", opt.param_file);
		mls.ns = read_param_file(opt.param_file);
		sub_prob = mlogit_sub_prob;
		sub_lprob = mlogit_sub_lprob;
	}
	
	/* open fastq file for output */
	if (opt.ampliclust_command || opt.write_fastq_and_quit) {
		fp = fopen(opt.ac_fastq_file, "w");
		if (opt.ampliclust_command) {
			MAKE_1ARRAY(in->not_input, n_read);	/* overestimate */
			in->n_hash_excluded = 0;
		}
	}
	
	double *mll = NULL;
	if (opt.min_log_likelihood > 0)
		MAKE_1ARRAY(mll, n_read);
	
	/* extract the desired region of reference sequences from whole genome
	 * reference FASTA files: output as FASTA files; recall this region
	 * contains the target region, but may be longer
	 */
	for (int j = 0; j < N_FILES; ++j) {

		mmessage(INFO_MSG, NO_ERROR, "samtools index region %zu-%zu in reference %s\n", start_reference[j] + 1, end_reference[j], csome_names[j]);

		extract_ref(opt_rf.samtools_command, csome_names[j], //opt.ref_names[j],
			start_reference[j], end_reference[j],
			opt.fsa_files[j], opt_rf.extracted_rf[j], &opt_rf);

		/* read in reference genomes */
		if ((err = read_fastq(opt_rf.extracted_rf[j], &fds[j], &fop)))
			exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "Reading '%s' "
					  "failed with error '%s' (%d).\n",
					  opt_rf.extracted_rf[j], fastq_error_message(err),
					  err));
	}

	/* compute posterior probability read from each subgenome;
	 * screen reads on minimum log likelihood
	 */
	n_read = 0;
	size_t coverA = 0;	/* number of reads assigned to subgenome A */
	size_t n_addl_excluded = 0;
	double A_expected_coverage = 0;
	unsigned char show = opt.display_alignment;

	//int tmp_debug = 0;
	for (merge_hash *me = mh; me != NULL; me = me->hh.next) {
		
		if (me->exclude) {
			if (opt.ampliclust_command)
				in->not_input[in->n_hash_excluded++] = n_read;
			continue;
		}
		
		double sum = 0, max = -INFINITY;
		//show = n_read == 955;

		for (unsigned int j = 0; j < N_FILES; ++j) {
			sam_entry *se = &sds[j]->se[me->indices[j][0]];
			// se->pos should be adjusted since the reference is a selected region, but do not need to consider the RC of B, since read is aligned to forward of B
			pp[j][n_read] = ll_align(se, n_read,
				&fds[j]->reads[fs_index[j]], &mls, &show,
				se->pos - start_reference[j] - 1,	/* 0-based, relative start of alignment in previously extracted reference */
				debug_level); //, n_read == 955 ? DEBUG_III : debug_level);
//						fprintf(stderr, "%s = %lf  ", se->name, pp[j][n_read]);
			ll[j][n_read] = pp[j][n_read];
			//if (isnan(ll[j][n_read])) tmp_debug = DEBUG_II;
			if (max < pp[j][n_read])
				max = pp[j][n_read];
		}
//		fprintf(stderr, "\n");
		
		if (opt.min_log_likelihood > 0) {
			mll[n_read] = max;
		} else if (isfinite(opt.min_log_likelihood)
			   && max < opt.min_log_likelihood) {
			me->exclude = 1;
			mmessage(INFO_MSG, NO_ERROR, "Read %s excluded for max "
				 "alignment log likelihood %f (below threshold "
				 "%f).\n", sds[0]->se[me->indices[0][0]].name,
				 		max, opt.min_log_likelihood);
			++n_addl_excluded;
			if (opt.ampliclust_command)
				in->not_input[in->n_hash_excluded++] = n_read;
			show = opt.display_alignment;
			continue;
		}
		
		debug_msg(show || debug_level > QUIET, debug_level, "Log likelihoods:");
		for (unsigned int j = 0; j < N_FILES; ++j) {
			debug_msg_cont(show || debug_level > QUIET, debug_level,
					   " %f", pp[j][n_read]);
			pp[j][n_read] = exp(pp[j][n_read] - max);
			sum += pp[j][n_read];
		}

		double max_p = -INFINITY;
		unsigned int max_p_index = 0;

		debug_msg_cont(show || debug_level > QUIET, debug_level, "\n");
		debug_msg(show || debug_level > QUIET, debug_level, "Probabilities:");
		for (unsigned int j = 0; j < N_FILES; ++j) {
			pp[j][n_read] = pp[j][n_read] / sum;
			debug_msg_cont(show || debug_level > QUIET, debug_level,
					   " %f", pp[j][n_read]);
			if (max_p < pp[j][n_read]) {
				max_p = pp[j][n_read];
				max_p_index = j;
			}
			if (opt.error_file)
				output_error_data(opt.error_file,
						  &sds[j]->se[me->indices[j][0]],
						  &fds[j]->reads[fs_index[j]],
						  log(pp[j][n_read]));
		}
		if (max_p_index > 1)
			debug_msg(debug_level > QUIET, debug_level, "***!!!***");
		debug_msg_cont(show || debug_level > QUIET, debug_level, "\n");
		
		debug_msg(debug_level > QUIET, debug_level, "Log Probabilities:");
		for (unsigned int j = 0; j < N_FILES; ++j)
			debug_msg_cont(debug_level > QUIET, debug_level, " %f",
					   log(pp[j][n_read]));
		debug_msg_cont(debug_level > QUIET, debug_level, "\n");
		
		/* write fastq of selected reads for amplici (now A and B many contain different reads)*/
		if (opt.ampliclust_command || opt.write_fastq_and_quit) {
			sam_entry *se = &sds[0]->se[me->indices[0][0]];
			fprintf(fp, "@%s ppA=%.6e\n", se->name, pp[0][n_read]);
			fwrite_nuc_segment(fp, se->read, XY_ENCODING, 0,
					   se->read->len);
			fprintf(fp, "\n+\n");
			fwrite_qual_sequence(fp, se->qual);
			fprintf(fp, "\n");
			
			/* record error free probability */
			++n_included_reads;
		}
		
		if (pp[0][n_read] > 0.5)
			++coverA;
		A_expected_coverage += pp[0][n_read];
		++n_read;
		show = opt.display_alignment;

	}
	if (opt.min_log_likelihood > 0) {
		qsort(mll, n_read, sizeof(double), double_compare);
		debug_msg(1, debug_level, "mll: ");
		fprint_doubles(stderr, mll, n_read, 6, 1);
		fprintf(stderr, "%f %f %f %f %f %f %f\n", mll[(int)(n_read*0.025)], mll[(int)(n_read*0.05)], mll[(int)(n_read*0.10)], mll[(int)(n_read*0.5)], mll[(int)(n_read*0.9)], mll[(int)(n_read*.95)], mll[(int)(n_read*.975)]);
	} else if (isfinite(opt.min_log_likelihood)) {
		mmessage(INFO_MSG, NO_ERROR, "%zu reads fail to achieve minimum"
			 " log likelihood %f.\n", n_addl_excluded,
			 opt.min_log_likelihood);
	}
	
	for (unsigned int j = 0; j < N_FILES; ++j) {
		double *new_pp = realloc(pp[j], n_read * sizeof **pp);
		double *new_ll = realloc(ll[j], n_read * sizeof **ll);

		if (!new_pp || !new_ll)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
					"(reallocate) posterior probability (%zu).\n",
					n_read);
		pp[j] = new_pp;
		ll[j] = new_ll;
	}
	
	if (opt.error_file)
		fclose(opt.error_file);
	
	debug_msg(debug_level > QUIET, debug_level, "Read coverage probability:"
		  " %f %f\n", (double) coverA / n_read,
		  (double) (n_read - coverA) / n_read);
	mmessage(INFO_MSG, NO_ERROR, "Expected coverage: %f %f (%zu)\n",
		 A_expected_coverage, n_read - A_expected_coverage, n_read);
	
	double B_expected_coverage = n_read - A_expected_coverage;
	double min_expected_coverage = MIN(A_expected_coverage,
					   B_expected_coverage);
	/* finish write selected reads as fastq and end */
	if (opt.ampliclust_command || opt.write_fastq_and_quit) {
		
		fclose(fp);
		
		if (opt.write_fastq_and_quit)
			return mmessage(INFO_MSG, NO_ERROR, "Finished writing "
					"fastq file '%s'.\n", opt.ac_fastq_file);
		
		in->n_observation = n_read;
		
		unsigned int *new_not_input = realloc(in->not_input,
							  n_read * sizeof *in->not_input);
		if (!new_not_input)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
					"(reallocate) not_input (%zu).\n", n_read);
		in->not_input = new_not_input;
		in->proptest = opt.proptest_screen;
		
		/* SHOULD WE CALL R SCRIPT HERE TO DO THE TEST AND TAKE THE RESULT TO FILTER MORE READS:
		 NOW JUST TAKE THE RESULTS FROM THE COMMAND LINE
		 */
		
		/* call ampliclust: it will give 3 files: err.out, res.out, res.fa. res.out will be used to next step */
		unsigned int cmd_len = strlen(opt.ampliclust_command) + strlen(opt.ac_fastq_file) + strlen(" -f -lb  -ll -Inf -o ") + strlen(opt.ac_outfile) + (int)(log10(opt.ac_low_bound) + 1) + 8;
		mmessage(INFO_MSG, NO_ERROR, "Length of command: %u\n", cmd_len);
		char *command = malloc(cmd_len * sizeof *command);
		sprintf(command, "%s -f %s -lb %f -ll -Inf -o %s",
			opt.ampliclust_command, opt.ac_fastq_file, opt.ac_low_bound, opt.ac_outfile);
		
		mmessage(INFO_MSG, NO_ERROR, "Running ampliclust: '%s'\n",
			 command);
		system(command);
		free(command);
		
		command = malloc((strlen(opt.ac_outfile) + strlen(".out") + 1) * sizeof *command);
		sprintf(command, "%s.out", opt.ac_outfile);
		FILE *ampli_results = fopen(command, "r");
		if (!ampli_results)
			return mmessage(ERROR_MSG, FILE_OPEN_ERROR, command);
		free(command);
		
		in->n_observation = n_read;
		
		read_ampliclust_results(ampli_results, in);
		
		mmessage(INFO_MSG, NO_ERROR, "Excluded reads:");
		unsigned int n_hash = 0, n_haplotype = 0;
		size_t n_excluded = 0;
		
		for (merge_hash *me = mh; me != NULL; me = me->hh.next) {
			if (me->exclude)
				continue;
			if (in->assignment[n_hash] >= (unsigned int) (4U - in->proptest)) {
				fprintf(stderr, " %d", n_hash);
				++n_excluded;
				if (in->assignment[n_hash] > n_haplotype)
					n_haplotype = in->assignment[n_hash];
			}
			++n_hash;
		}
		
		fprintf(stderr, " (%zu total from %u haplotypes)\n",
			n_excluded, n_haplotype - 3 + in->proptest);
		
		if (n_excluded) {
			size_t n_idx1 = 0, n_idx2 = 0;
			double **new_pp = malloc(N_FILES * sizeof **pp);
			if (!new_pp)
				return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"posterior probability.\n");
			n_read = n_read - n_excluded;
			for (int j = 0; j < N_FILES; ++j)
				new_pp[j] = malloc(n_read * sizeof *new_pp);
			A_expected_coverage = 0;
			if (opt.min_log_likelihood > 0)
				fprintf(stderr, "Maximum log likelihoods:");
			for (merge_hash *me = mh; me != NULL; me = me->hh.next) {
				if (me->exclude)
					continue;
				if (in->assignment[n_idx1] >= (unsigned int) (4U - in->proptest)) {
					me->exclude = 1;
				} else {
					for (int j = 0; j < N_FILES; ++j)
						new_pp[j][n_idx2] = pp[j][n_idx1];
					A_expected_coverage += new_pp[0][n_idx2];
					if (opt.min_log_likelihood > 0)
						fprintf(stderr, " %f", mll[n_idx1]);
					++n_idx2;
				}
				++n_idx1;
			}
			if (opt.min_log_likelihood > 0)
				fprintf(stderr, "\n");
			for (int j = 0; j < N_FILES; ++j) {
				free(pp[j]);
				pp[j] = new_pp[j];
			}
			free(new_pp);
			B_expected_coverage = n_read - A_expected_coverage;
			min_expected_coverage = MIN(A_expected_coverage,
							B_expected_coverage);
			mmessage(INFO_MSG, NO_ERROR, "Updated expected coverage"
				 ": %f %f (%zu)\n", A_expected_coverage,
				 B_expected_coverage, n_read);
			if (A_expected_coverage < opt.min_expected_coverage
				|| B_expected_coverage < opt.min_expected_coverage)
				return mmessage(ERROR_MSG, INTERNAL_ERROR,
						"Coverage is too low.\n");
			
		}
	}
//	if (min_expected_coverage < opt.min_expected_coverage)
//		return mmessage(ERROR_MSG, INTERNAL_ERROR,
//				"Sugenomic coverage is too low: A=%.2f B=%.2f.\n",
//				A_expected_coverage, B_expected_coverage);
	if (opt.min_log_likelihood > 0)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "Stopping after "
				"writing maximum log likelihoods.\n");
	
	xy_t *obs_nuc = malloc(n_read * sizeof *obs_nuc);
	qual_t *obs_q = malloc(n_read * sizeof *obs_q);
	char *genome_src = malloc(n_read * sizeof *genome_src);
	char *covers = malloc(n_read * sizeof *covers);
	unsigned int *obs_rpos = malloc(n_read * sizeof *obs_rpos);
	int *rd_idxA = malloc(n_read * sizeof *rd_idxA);
	size_t num_nuc[NUM_NUCLEOTIDES];
	double ebaseA[NUM_NUCLEOTIDES];
	double ebaseB[NUM_NUCLEOTIDES];
	xy_t nuc1, nuc2 = 0, nuc3 = 0;
	size_t num_nuc1, num_nuc2, num_nuc3;
	
	/* open and write header of vcf files */
	FILE *vcf_fp[N_FILES];
	for (unsigned int i = 0; i < N_FILES; ++i) {
		vcf_fp[i] = NULL;
		if (opt.vcf_files[i])
			vcf_fp[i] = fopen(opt.vcf_files[i], "w");
		if (!vcf_fp[i])
			mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt.vcf_files[i]);
		print_vcf_header(vcf_fp[i], opt.vcf_opt, opt.fsa_files[i],
							 opt.sample_name);
   	}
	
	/* store information for post-hoc tests of allelic coverage */
	/* [BUG,KSD] no memory allocation checks */
	size_t region_len = 0;
	size_t *haplotype_posns[N_FILES];	/* positions of het calls */
	unsigned int n_segregatingA = 0, n_segregatingB = 0;
	double *hapA_prop = NULL, *hapB_prop = NULL, *hapA_covg = NULL, *hapB_covg = NULL;
	xy_t *hapA_dom_nuc = NULL, *hapB_dom_nuc = NULL;

	if (opt.posthoc_coverage_test) {
		region_len = end_reference[1] - start_reference[1];
		if (end_reference[0] - start_reference[0] > region_len)
			region_len = end_reference[0] - start_reference[0];

		haplotype_posns[0] = calloc(region_len, sizeof *haplotype_posns[0]);
		haplotype_posns[1] = calloc(region_len, sizeof *haplotype_posns[1]);
	
		/* proportion of reads matching dominant haplotype in each subgenome */
		hapA_prop = malloc(region_len * sizeof *hapA_prop);
		hapB_prop = malloc(region_len * sizeof *hapB_prop);
	
		/* coverage of dominant haplotype in each subgenome */
		hapA_covg = malloc(region_len * sizeof *hapA_covg);
		hapB_covg = malloc(region_len * sizeof *hapB_covg);
	
		/* dominant nucleotides at each position of dominant haplotype */
		hapA_dom_nuc = calloc(region_len, sizeof *hapA_dom_nuc);
		hapB_dom_nuc = calloc(region_len, sizeof *hapB_dom_nuc);
	}

	size_t start_target[N_FILES];
	size_t end_target[N_FILES];

	start_target[0] = rf_info->start_A;	/* 0-based, inclusive */
	end_target[0] = rf_info->end_A;		/* 0-based, exclusive */
	start_target[1] = rf_info->start_B;
	end_target[1] = rf_info->end_B;
	
	//debug_level = DEBUG_II;
	fprintf(stderr, "Start genotyping\n");
	FILE *final_out = NULL;
	if (opt.output_file) {
		final_out = fopen(opt.output_file, "w");
		fprintf(final_out, "ChromA  ChromB	PositionA	PositionB	Genotype_call	Call_A genome	Call_Bgenome	Major allele	Minor allele	PP(0,0)	PP(0,1)	PP(0,2)	PP(1,0)	PP(1,1)	PP(1,2)	PP(2,0)	PP(2,1)	PP(2,2)	PA(0)	PA(1)	PA(2)	PB(0)	PB(1)  PB(2)	CovA	CovB\n");
	}
	//the reference index of A AND B should be adjusted according to the cigar string, this only genotype on the site that are not -/A or A/-
	/* finally: march along reference positions and genotype */
	for (size_t posA = start_target[0]; posA < end_target[0]; ++posA) {
	//for (size_t posA = 77536234; posA < 77536235; ++posA) {
		//ref_entry *re = rf_info;
		size_t target_a = posA; 			/* 0-based, absolute position and ... */
		size_t site = target_a - start_target[0];	/* relative position in aligned region of subgenome A */
		
		/* skip insertions in genome A; insertions in genome B are ignored by loop over posA */
		if (rf_info->map_A_to_B[site] == -1) {
			debug_msg(debug_level > ABSOLUTE_SILENCE, debug_level,
				  "Site %zu is not a homologous position, skip\n", site);
			continue;
		}

		size_t target_b = start_target[1] + rf_info->map_A_to_B[site]; /* 0-based, absolute position within aligned region of genome B */

		debug_msg(debug_level > ABSOLUTE_SILENCE, debug_level,
			  "Site A %zu (target_a = %zu; target_b = %zu, "
			  "A_align_B = %d; start_target[0] = %u)\n", site, target_a, target_b, rf_info->map_A_to_B[site], start_target[0]);

		/* count each of alleles across all reads */
		for (int b = 0; b < NUM_NUCLEOTIDES; ++b) {
			num_nuc[b] = 0;
			ebaseA[b] = 0;
			ebaseB[b] = 0;
		}
		size_t n_cover_A = 0;
		char_t ref_base[N_FILES] = {IUPAC_A, IUPAC_A};	/* [KSD,TODO] Assumes N_FILES == 2 */
		size_t ref_pos[N_FILES] = {target_a, target_b};	/* [KSD,TODO] Assumes N_FILES == 2 */

		/* extract read position aligned to genome A */
		int no_alt_allele = 0;
		n_read = 0;
		for (merge_hash *me = mh; me != NULL; me = me->hh.next) {
			
			if (me->exclude)
				continue;
			
			/* A alignment */
			sam_entry *se = &sds[0]->se[me->indices[0][0]];
			unsigned int rd_idx = 0;
			size_t rf_idx = se->pos - 1;

			rd_idxA[n_read] = -1;	/* default: assume not covering */
//debug_msg(debug_level > ABSOLUTE_SILENCE, debug_level, "Read %s (%u), rf_idx=%zu, n_ashes=%u\n", se->name, n_read, rf_idx, se->cig->n_ashes);
			
			if (rf_idx > target_a) {
				++n_read;
				continue;
			}
			
			for (unsigned int j = 0; j < se->cig->n_ashes; ++j) {
				
				debug_msg(debug_level > DEBUG_I, debug_level, "Read %s (%u) ash %u%c (%u), rd_idx=%d, target_a=%d, rf_idx=%d\n",
					se->name, n_read, se->cig->ashes[j].len, cigar_char[se->cig->ashes[j].type], j, rd_idx, target_a, rf_idx);
				/* reference consumed */
				if (se->cig->ashes[j].type == CIGAR_DELETION
					|| se->cig->ashes[j].type == CIGAR_SKIP) {
					
					/* read deletes or skips desired site */
					if (rf_idx + se->cig->ashes[j].len > target_a)
						break;
					
					rf_idx += se->cig->ashes[j].len;
					continue;
					
					/* read consumed */
				} else if (se->cig->ashes[j].type == CIGAR_SOFT_CLIP
					   || se->cig->ashes[j].type == CIGAR_INSERTION) {
					
					rd_idx += se->cig->ashes[j].len;
					continue;
					
					/* neither consumed: HARD_CLIP, PAD */
				} else if (se->cig->ashes[j].type != CIGAR_MATCH
					   && se->cig->ashes[j].type != CIGAR_MMATCH
					   && se->cig->ashes[j].type != CIGAR_MISMATCH) {
					continue;
				}
				
				/* read and reference consumed */
				/* desired site within this ash */
				if (rf_idx + se->cig->ashes[j].len > target_a) {
					rd_idxA[n_read] = rd_idx + target_a - rf_idx;
					ref_base[0] = fds[0]->reads[fs_index[0] + target_a - start_reference[0]];
					debug_msg(debug_level >= DEBUG_I, debug_level, "Read %s (%u) covers target %c at position %u in %u%c from (%zu-%zu) with rd_idx = %lu\n",
						se->name, n_read, iupac_to_char[ref_base[0]], target_a, se->cig->ashes[j].len, cigar_char[se->cig->ashes[j].type], rf_idx, rf_idx + se->cig->ashes[j].len, rd_idx + target_a - rf_idx);
					n_cover_A++;
					break;
				}
				
				rf_idx += se->cig->ashes[j].len;
				rd_idx += se->cig->ashes[j].len;
			}
			++n_read;
		}

		debug_msg(debug_level > ABSOLUTE_SILENCE, debug_level,
			  "Site A %zu (target_a = %zu; target_b = %zu, "
			  "start_target[0] = %u): A coverage = %zu (%zu)\n", site,
			  ref_pos[0], ref_pos[1], start_target[0], n_cover_A, n_read);
		
		/* no reads cover this site after accounting for soft-clipping */
		if (!n_cover_A)
			continue;
		
		/* extract set of nucleotides aligned to this position in both
		 * alignments
		 */
		size_t n_cover = 0;
		n_read = 0;
		for (merge_hash *me = mh; me != NULL; me = me->hh.next) {
			
			if (me->exclude)
				continue;

			covers[n_read] = 0;	/* assume: not covering */

			if (rd_idxA[n_read] == -1) {
				++n_read;
				continue;
			}
			
			/* B alignment */
			sam_entry *se = &sds[1]->se[me->indices[1][0]];
			unsigned int rd_idx = 0;
			size_t rf_idx = //se->pos - 1;
				rf_info->strand_B 
					? se->pos - 1 + se->cig->length_rf - 1 // 0-based
					: se->pos - 1;

			debug_msg(debug_level > DEBUG_I, debug_level, "Read "
				"%s (%u)%s, pos=%u, rf_len=%u, target_b=%zu, "
				"rf_idx=%zu, cigar=", 
				me->indices[0]	// TODO: hack
				? sds[0]->se[me->indices[0][0]].name
				: sds[1]->se[me->indices[1][0]].name,
				n_read, rf_info->strand_B ? " [rc B]" : "",
				se->pos, se->cig->length_rf,
				target_b, rf_idx);
			debug_call(debug_level > ABSOLUTE_SILENCE, debug_level, print_cigar(stderr, se));
			debug_msg_cont(debug_level > ABSOLUTE_SILENCE, debug_level, "\n");
/**/

			
			/* aligned to forward strand of B subgenome */
			if (!rf_info->strand_B) {
				if (rf_idx > target_b) {
					n_read++;
					continue;
				}
				for (unsigned int j = 0; j < se->cig->n_ashes; ++j) {
					
					debug_msg(debug_level > DEBUG_I, debug_level,
						"Read %s (%u) cigar %u%c (%u), "
						"rd_idxA[%u] = %d, rd_idx=%d, "
						"target_b=%d, rf_idx=%d\n", se->name,
						n_read, se->cig->ashes[j].len,
						cigar_char[se->cig->ashes[j].type],
						j, n_read, rd_idxA[n_read], rd_idx,
						target_b, rf_idx);

					/* reference consumed */
					if (se->cig->ashes[j].type == CIGAR_DELETION
						|| se->cig->ashes[j].type == CIGAR_SKIP) {
						
						/* read deletes or skips desired site */
						if (rf_idx + se->cig->ashes[j].len > target_b)
							break;
						
						rf_idx += se->cig->ashes[j].len;
						continue;
						
					/* read consumed */
					} else if (se->cig->ashes[j].type == CIGAR_SOFT_CLIP
						   || se->cig->ashes[j].type == CIGAR_INSERTION) {
						
						rd_idx += se->cig->ashes[j].len;
						continue;
						
					/* neither consumed: HARD_CLIP, PAD */
					} else if (se->cig->ashes[j].type != CIGAR_MATCH
						   && se->cig->ashes[j].type != CIGAR_MMATCH
						   && se->cig->ashes[j].type != CIGAR_MISMATCH) {
						continue;
					}
					
					/* read and reference consumed */
					/* dataset is set of nucleotides that align to both A
					 * and B at same position
					 */
					if (rf_idx + se->cig->ashes[j].len > target_b
						&& rd_idxA[n_read] == (int) (rd_idx + target_b - rf_idx)) {
						covers[n_read] = 1;
						obs_rpos[n_cover] = (int) (rd_idx + target_b - rf_idx);
						obs_nuc[n_cover] = get_nuc(se->read, XY_ENCODING, obs_rpos[n_cover]);
						ebaseA[obs_nuc[n_cover]] += pp[0][n_read];
						ebaseB[obs_nuc[n_cover]] += 1 - pp[0][n_read];
						++num_nuc[obs_nuc[n_cover]];
						obs_q[n_cover] = get_qual(se->qual, obs_rpos[n_cover]);
						ref_base[1] = fds[1]->reads[fs_index[1] + target_b - start_reference[1]];
						debug_msg(debug_level >= DEBUG_I, debug_level, "[SUCCESS] Read %s (%u) covers target %c at position %u in %u%c from (%zu-%zu) with nucleotide %c at rd_idx = %lu and pp[A]=%f\n",
							se->name, n_read, iupac_to_char[ref_base[1]], target_b, se->cig->ashes[j].len, cigar_char[se->cig->ashes[j].type], rf_idx, rf_idx + se->cig->ashes[j].len, xy_to_char[obs_nuc[n_cover]], obs_rpos[n_cover], pp[0][n_read]);
						++n_cover;
						break;
					} else if (rf_idx + se->cig->ashes[j].len > target_b) {
						debug_msg(debug_level > DEBUG_I, debug_level, "[FAILURE] Read %s (%u) covers target position %u in %u%c from (%zu-%zu) but rd_idx = %lu\n",
							se->name, n_read, target_b, se->cig->ashes[j].len, cigar_char[se->cig->ashes[j].type], rf_idx, rf_idx + se->cig->ashes[j].len, rd_idx + target_b - rf_idx);
					}
					rf_idx += se->cig->ashes[j].len;
					rd_idx += se->cig->ashes[j].len;
				}
			/* B is reverse complemented in alignment to A: reverse complement the sam alignment */
			} else {
				if (rf_idx < target_b) {
					n_read++;
					continue;
				}
				
				for (unsigned int j = se->cig->n_ashes; j-- > 0; ) {
					debug_msg(debug_level > DEBUG_I, debug_level, "Read %s (%u) cigar %u%c (%u), rd_idxA[%u] = %d, rd_idx=%d, target_b=%d, rf_idx=%d\n",
						se->name, n_read, se->cig->ashes[j].len, cigar_char[se->cig->ashes[j].type], j, n_read, rd_idxA[n_read], rd_idx, target_b, rf_idx);
					/* reference consumed */
					if (se->cig->ashes[j].type == CIGAR_DELETION
						|| se->cig->ashes[j].type == CIGAR_SKIP) {
						
						/* read deletes or skips desired site */
						if (rf_idx - se->cig->ashes[j].len < target_b)
							break;
						
						rf_idx -= se->cig->ashes[j].len;
						continue;
						
					/* read consumed */
					} else if (se->cig->ashes[j].type == CIGAR_SOFT_CLIP
						   || se->cig->ashes[j].type == CIGAR_INSERTION) {
						
						rd_idx += se->cig->ashes[j].len;
						continue;
						
					/* neither consumed: HARD_CLIP, PAD */
					} else if (se->cig->ashes[j].type != CIGAR_MATCH
						   && se->cig->ashes[j].type != CIGAR_MMATCH
						   && se->cig->ashes[j].type != CIGAR_MISMATCH) {
						continue;
					}
					
					/* read and reference consumed */
					/* dataset is set of nucleotides that align to both A
					 * and B at same homoeologous position
					 */
					if (rf_idx - se->cig->ashes[j].len < target_b
						&& rd_idxA[n_read] == (int) (rd_idx + rf_idx - target_b)) {
						covers[n_read] = 1;
						obs_rpos[n_cover] = (int) (rd_idx + rf_idx - target_b);
						obs_nuc[n_cover] = get_nuc(se->read, XY_ENCODING, se->read->len - obs_rpos[n_cover] - 1);
						/* REVERSE complement the nuc! */
						obs_nuc[n_cover] = xy_to_rc[obs_nuc[n_cover]];
						ebaseA[obs_nuc[n_cover]] += pp[0][n_read];
						ebaseB[obs_nuc[n_cover]] += 1 - pp[0][n_read];
						++num_nuc[obs_nuc[n_cover]];
						obs_q[n_cover] = get_qual(se->qual, se->read->len - obs_rpos[n_cover] - 1);
						ref_base[1] = fds[1]->reads[fs_index[1] + target_b - start_reference[1]]; // need to reverse complement B
						debug_msg(debug_level >= DEBUG_I, debug_level, "[SUCCESS] Read %s (%u) covers target %c at position %u in %u%c from (%zu-%zu) with nucleotide %c at rd_idx = %lu and pp[A]=%f\n",
							se->name, n_read, iupac_to_char[ref_base[1]], target_b, se->cig->ashes[j].len, cigar_char[se->cig->ashes[j].type], rf_idx, rf_idx + se->cig->ashes[j].len, xy_to_char[obs_nuc[n_cover]], obs_rpos[n_cover], pp[0][n_read]);
						++n_cover;
						break;
					} else if (rf_idx - se->cig->ashes[j].len < target_b) {
						debug_msg(debug_level > DEBUG_I, debug_level, "[FAILURE] Read %s (%u) covers target position %u in %u%c from (%zu-%zu) but rd_idx = %lu\n",
							se->name, n_read, target_b, se->cig->ashes[j].len, cigar_char[se->cig->ashes[j].type], rf_idx, rf_idx + se->cig->ashes[j].len, rd_idx + rf_idx - target_b);
					}
					rf_idx -= se->cig->ashes[j].len;
					rd_idx += se->cig->ashes[j].len;
				}
			}
			++n_read;
		}

		//debug_msg_cont(debug_level > QUIET, debug_level, "\n");
		debug_msg(debug_level > ABSOLUTE_SILENCE, debug_level,
			  "Site %zu (target_a = %zu; target_b = %zu, "
			  "start_target[0] = %u): A + B coverage = %zu (%zu)\n",
			  site, ref_pos[0], ref_pos[1], start_target[0], n_cover,
			  						n_read);

		/* no reads cover this site IN BOTH GENOME */
		if (!n_cover)
			continue;
		
		mmessage(INFO_MSG, NO_ERROR, "Expected counts genome A: ");
		fprint_doubles(stderr, ebaseA, NUM_NUCLEOTIDES, 3, 1);
		mmessage(INFO_MSG, NO_ERROR, "Expected counts genome B: ");
		fprint_doubles(stderr, ebaseB, NUM_NUCLEOTIDES, 3, 1);
		
		double penalty_A = 0, penalty_B = 0;
		double sumA = 0, sumB = 0;
		double mina = 0, minb = 0;

		xy_t modeA = XY_A, modeB = XY_A;
		for (unsigned int i = 0; i < NUM_NUCLEOTIDES; ++i) {
			if (ebaseA[i] > mina) {
				mina = ebaseA[i];
				modeA = i;
			}
			if (ebaseB[i] > minb) {
				minb = ebaseB[i];
				modeB = i;
			}
			sumA += ebaseA[i];
			sumB += ebaseB[i];
		}
		double ecoverage[N_FILES] = {sumA, sumB};	/* [KSD,TODO] assumes N_FILES == 2 */
		/* [KSD,TODO,BUG?] If current site is only subreference
		 * mismatch in a distance and this is an allelic SNP, then CAPG
		 * will lose coverage of the nonmatching allele, leading to 
		 * screening away of evidence. Compute expected coverage WITHOUT
		 * current site? Or is this just part of the loss of 
		 * non-identifiable cases?
		 */

		if (sumA == 0 || sumB == 0) {	/* [KSD,TODO] allow user to impose a minimum bound here */
			debug_msg(debug_level > ABSOLUTE_SILENCE, debug_level,
				"Expected coverage of site %zu in at least "
				"one genome is 0: skip site.\n", site);
			continue;
		}
		
		/* identify most and second most common allele */
		nuc1 = nuc2 = XY_A;
		num_nuc1 = num_nuc[nuc1];
		num_nuc2 = num_nuc3 = 0;
		for (int b = 1; b < NUM_NUCLEOTIDES; ++b)
			if (num_nuc[b] > num_nuc1) {
				nuc3 = nuc2;
				num_nuc3 = num_nuc2;
				nuc2 = nuc1;
				num_nuc2 = num_nuc1;
				nuc1 = b;
				num_nuc1 = num_nuc[b];
			} else if (num_nuc[b] > num_nuc2) {
				nuc3 = nuc2;
				num_nuc3 = num_nuc2;
				nuc2 = b;
				num_nuc2 = num_nuc[b];
			} else if (num_nuc[b] > num_nuc3) {
				nuc3 = b;
				num_nuc3 = num_nuc[b];
			}
		if (nuc2 == nuc1) {
			no_alt_allele = 1;
			nuc2 = XY_C;
		} else if (!num_nuc[nuc2]) {
			no_alt_allele = 1;
		}
		
		/* check for possible 3rd allele: we do not handle this yet */
		if (num_nuc3 > opt.biallelic_screen
			* min_expected_coverage / 2) {
			debug_msg(debug_level > QUIET, debug_level, "Evidence "
				  "of third nucleotide, will not genotype this"
				  "site (%c=%zu, %c=%zu, %c=%zu).  To change "
				  "the behavior, see command-line option "
				  "--biallelic\n", xy_to_char[nuc1], num_nuc1,
				  xy_to_char[nuc2], num_nuc2,
				  xy_to_char[nuc3], num_nuc3);
			for (int i = 0; i < N_FILES; ++i) {
				int one_alt = 0;
				
				if (!vcf_fp[i])
					continue;
				if (i == 0)
					fprintf(vcf_fp[i], "%s", rf_info->name_A);
				else
					fprintf(vcf_fp[i], "%s", rf_info->name_B);
				fprintf(vcf_fp[i], "\t%lu\t.\t%c\t", ref_pos[i] + 1,
					iupac_to_char[ref_base[i]]);

				if (xy_to_iupac[nuc1] != ref_base[i]) {
					fprintf(vcf_fp[i], "%c",
						xy_to_char[nuc1]);
					one_alt = 1;
				}
				if (!no_alt_allele && nuc2 != nuc1
					&& xy_to_iupac[nuc2] != ref_base[i]) {
					if (one_alt)
						fputc(',', vcf_fp[i]);
					fprintf(vcf_fp[i], "%c",
						xy_to_char[nuc2]);
					one_alt = 1;
				}
				if (!one_alt)
					fputc('.', vcf_fp[i]);
				fprintf(vcf_fp[i], "\t.\tal2\t.\t.\t.\n");
			}
			continue;
		}

		fprintf(stderr, "Observed nucleotides (%zu): ", n_cover);
		for (unsigned int i = 0; i < n_cover; ++i)
			fprintf(stderr, "%c", xy_to_char[obs_nuc[i]]);
		fprintf(stderr, "\n");
		fprintf(stderr, "Observed   qualities (%zu): ", n_cover);
		for (unsigned int i = 0; i < n_cover; ++i)
			fprintf(stderr, "%c", (char)(obs_q[i]
						+ MIN_ASCII_QUALITY_SCORE));
		fprintf(stderr, "\n");
		
		/* estimate posterior probability of each genotype */
		/*mls.pos = posA - start_target[0];*/
		
		double lprob[9], gprob[9];
		double mprob = 0;
		int g1_max = 0, g2_max = 0;
		
		if (karin_version) {/* to replace old version */
			
			double tmp1, tmp2, tmp1a, tmp1b, tmp2a, tmp2b;
			double max = -INFINITY, den = 0;
			/* consider each genotype at the current locus */
			for (int g1 = 0; g1 <= 2; ++g1) {
				for (int g2 = 0; g2 <= 2; ++g2) {
					n_read = 0;
					n_cover = 0;
					lprob[g1 * 3 + g2] = 0;
					for (merge_hash *me = mh; me != NULL; me = me->hh.next) {
						if (me->exclude)
							continue;
						if (covers[n_read]) {
							mls.pos = obs_rpos[n_cover];

							/* log likelihood of A alignment */
							tmp1a = sub_prob_given_q_with_encoding(
											 ref_base[0],			/* homozygous ref base */
											 obs_nuc[n_cover], IUPAC_ENCODING, XY_ENCODING, obs_q[n_cover], 1, (void *)&mls);
							tmp1b = sub_prob_given_q_with_encoding(
											 !g1 ? xy_to_iupac[nuc1]	/* new genotype */
											 : g1 == 2 ? xy_to_iupac[nuc2]
											 : xy_to_iupac[nuc1] | xy_to_iupac[nuc2],
											 obs_nuc[n_cover], IUPAC_ENCODING, XY_ENCODING, obs_q[n_cover], 1, (void *)&mls);
							tmp1 = ll[0][n_read] - tmp1a + tmp1b;

							/* log likelihood of B alignment */
							/* [BROKEN, KSD 2/14/22] [FIXED, KSD 2/20/22] obs_nuc
							 * was reverse complemented to match subgenome A, but
							 * ref_base[1] is still on reverse-complemented
							 * subgenome B
							 */
							xy_t tmp = obs_nuc[n_cover];
							if (rf_info->strand_B) //RC if B is reversed
								tmp = xy_to_rc[obs_nuc[n_cover]];
							tmp2a = sub_prob_given_q_with_encoding(
											 ref_base[1],			/* homozygous ref base */
											 tmp, IUPAC_ENCODING, XY_ENCODING, obs_q[n_cover], 1, (void *)&mls);
							tmp2b = sub_prob_given_q_with_encoding(
											 !g2 ? xy_to_iupac[nuc1]	/* new genotype */
											 : g2 == 2 ? xy_to_iupac[nuc2]
											 : xy_to_iupac[nuc1] | xy_to_iupac[nuc2],
											 obs_nuc[n_cover], IUPAC_ENCODING, XY_ENCODING, obs_q[n_cover], 1, (void *)&mls);
							tmp2 = ll[1][n_read] - tmp2a + tmp2b;
							
							/* combine assuming uniform prior */
							lprob[g1 * 3 + g2] += log(exp(tmp1) + exp(tmp2));
							debug_msg(debug_level >= DEBUG_II, debug_level, "[Ref=%c, Alt=%c] Read %zu (%c,%d): Subgenome A %c%c -> %c%c M=(%d,%d): %f (%f-%f); Subgenome B %c%c -> %c%c M=(%d,%d): %f (%f-%f) (%f)\n", 
								xy_to_char[nuc1], xy_to_char[nuc2], n_cover, xy_to_char[obs_nuc[n_cover]], obs_q[n_cover],
								iupac_to_char[ref_base[0]], iupac_to_char[ref_base[0]],
								g1<=1 ? xy_to_char[nuc1] : xy_to_char[nuc2], g1>=1 ? xy_to_char[nuc2] : xy_to_char[nuc1], g1, g2, tmp1, tmp1b, tmp1a,
								iupac_to_char[ref_base[1]], iupac_to_char[ref_base[1]],
								g2<=1 ? xy_to_char[nuc1] : xy_to_char[nuc2], g2>=1 ? xy_to_char[nuc2] : xy_to_char[nuc1], g1, g2, tmp2, tmp2b, tmp2a, lprob[g1 * 3 + g2]);
							++n_cover;
						}
						++n_read;
					}
					if (max < lprob[g1 * 3 + g2])
						max = lprob[g1 * 3 + g2];
				}
			}
			for (int g1 = 0; g1 <= 2; ++g1)
				for (int g2 = 0; g2 <= 2; ++g2) {
					gprob[g1 * 3 + g2] = exp(lprob[g1 * 3 + g2] - max);
					den += gprob[g1 * 3 + g2];
				}
			for (int g1 = 0; g1 <= 2; ++g1)
				for (int g2 = 0; g2 <= 2; ++g2) {
					gprob[g1 * 3 + g2] /= den;
					fprintf(stderr, "M = (%u, %u) %c%c at site (%zu, %zu): %e\n", g1, g2, xy_to_char[nuc1], xy_to_char[nuc2], ref_pos[0] + 1, ref_pos[1] + 1, gprob[g1*3 + g2]);
					if (gprob[g1*3 + g2] > mprob) {
						mprob = gprob[g1*3 + g2];
						g1_max = g1;
						g2_max = g2;
					}
				}
			
		} else {
		for (int g1 = 0; g1 <= 2; ++g1) {
			for (int g2 = 0; g2 <= 2; ++g2) {
				double sum = 0;
				
				/* repeat simulation */
				for (unsigned int b = 0; b < opt.n_sample; ++b) {
					double max = -INFINITY, den = 0;
					
					/* simulate alignments (or genome source) */
					n_read = 0; //here n_read should be correct, should be total read defined before
					n_cover = 0;
					for (merge_hash *me = mh; me != NULL; me = me->hh.next) {
						if (me->exclude)
							continue;
						
						if (covers[n_read]) {
							genome_src[n_cover] = 'A';
							if (rand() / (RAND_MAX + 1.) > pp[0][n_read])
								genome_src[n_cover] = 'B';
							//fprintf(stderr, "Assigning %c to genome %c (%f)\n", xy_to_char[obs_nuc[n_read]], genome_src[n_cover], pp[0][n_read]);
							++n_cover;
						}
						++n_read;
					}
					//size_t dbg_site = 246;
					/* compute Pr(M|R,A) */
					for (int g1p = 0; g1p <= 2; ++g1p) {
						for (int g2p = 0; g2p <= 2; ++g2p) {
							
							/* compute likelihood of this genotype */
							n_read = 0;
							n_cover = 0;
							lprob[g1p*3 + g2p] = 0; //change it to a binomial prior, nuc2 and nuc1, p = 0.5
							for (merge_hash *me = mh; me != NULL; me = me->hh.next) {
								if (me->exclude)
									continue;
								if (covers[n_read]) {
									/*
									 if (posA == dbg_site) {
									 sam_entry *se = &sds[0]->se[me->indices[0][0]];
									 fprintf(stderr, "Read %zu (%s) is %c w/ q %u at position %u from genome %c with %c:",
									 n_cover, se->name, xy_to_char[obs_nuc[n_cover]], obs_q[n_cover], obs_rpos[n_cover], genome_src[n_cover],
									 genome_src[n_cover] == 'A'
									 ? iupac_to_char[!g1p ? xy_to_iupac[nuc1] : g1p == 2 ? xy_to_iupac[nuc2] : xy_to_iupac[nuc1] | xy_to_iupac[nuc2]]
									 : iupac_to_char[!g2p ? xy_to_iupac[nuc1] : g2p == 2 ? xy_to_iupac[nuc2] : xy_to_iupac[nuc1] | xy_to_iupac[nuc2]]);
									 }
									 */
									mls.pos = obs_rpos[n_cover];
									/* assume source is A genome */
									if (genome_src[n_cover] == 'A') {
										lprob[g1p*3 + g2p] += sub_prob_given_q_with_encoding(
																	 !g1p ? xy_to_iupac[nuc1]			/* A genome is nuc1nuc1 */
																	 : g1p == 2 ? xy_to_iupac[nuc2]		/* A genome is nuc2nuc2 */
																	 : xy_to_iupac[nuc1] | xy_to_iupac[nuc2],	/* A genome is nuc1nuc2 */
																	 obs_nuc[n_cover], IUPAC_ENCODING, XY_ENCODING, obs_q[n_cover], 1, (void *)&mls);
										/*
										 if (posA == dbg_site)
										 fprintf(stderr, " %e", sub_prob_given_q_with_encoding(
										 !g1p ? xy_to_iupac[nuc1]
										 : g1p == 2 ? xy_to_iupac[nuc2]
										 : xy_to_iupac[nuc1] | xy_to_iupac[nuc2],
										 obs_nuc[n_cover], IUPAC_ENCODING, XY_ENCODING, obs_q[n_cover], 1, (void *)&mls));
										 */
										if (g1p == 1 && penalty_A != 0){
											lprob[g1p*3 + g2p] += penalty_A;}
										/* assume source is B genome */
									} else {
										lprob[g1p*3 + g2p] += sub_prob_given_q_with_encoding(
																	 !g2p ? xy_to_iupac[nuc1]			/* B genome is nuc1nuc1 */
																	 : g2p == 2 ? xy_to_iupac[nuc2]		/* B genome is nuc2nuc2 */
																	 : xy_to_iupac[nuc1] | xy_to_iupac[nuc2],	/* B genome is nuc1nuc2 */
																	 obs_nuc[n_cover], IUPAC_ENCODING, XY_ENCODING, obs_q[n_cover], 1, (void *)&mls);
										/*
										 if (posA == dbg_site)
										 fprintf(stderr, " %e", sub_prob_given_q_with_encoding(
										 !g2p ? xy_to_iupac[nuc1]
										 : g2p == 2 ? xy_to_iupac[nuc2]
										 : xy_to_iupac[nuc1] | xy_to_iupac[nuc2],
										 obs_nuc[n_cover], IUPAC_ENCODING, XY_ENCODING, obs_q[n_cover], 1, (void *)&mls));
										 */
										if (g2p == 1 && penalty_B != 0){
											lprob[g1p*3 + g2p] += penalty_B;}
									}
									/*
									 if (posA == dbg_site)
									 fprintf(stderr, " (%f)\n", lprob[g1p*3 + g2p]);
									 */
									++n_cover;
								}
								++n_read;
							}
							if (max < lprob[g1p*3 + g2p])
								max = lprob[g1p*3 + g2p];
							/*
							 if (posA == dbg_site)
							 fprintf(stderr, "%u: (%u, %u) %f\n", b, g1p, g2p, lprob[g1p*3 + g2p]);
							 */
						}
					}
					
					/* scaled normalization */
					for (int g1p = 0; g1p <= 2; ++g1p)
						for (int g2p = 0; g2p <= 2; ++g2p)
							den += exp(lprob[g1p*3 + g2p] - max);
					sum += exp(lprob[g1*3 + g2] - max) / den;
					/*
					 if (posA == dbg_site)
					 fprintf(stderr, "%u (%u, %u): %e (%e; q=%u)\n", b, g1, g2, sum, max, obs_q[n_cover]);
					 */
				}
				gprob[g1*3 + g2] = sum / opt.n_sample;
				fprintf(stderr, "M = (%u, %u) %c%c at site (%zu, %zu): %e\n", g1, g2, xy_to_char[nuc1], xy_to_char[nuc2], ref_pos[0] + 1, ref_pos[1] + 1, sum / opt.n_sample);
				if (gprob[g1*3 + g2] > mprob) {
					mprob = gprob[g1*3 + g2];
					g1_max = g1;
					g2_max = g2;
				}
			}
		}
		}
		fprintf(stderr, "Genotype (%4zu, %4zu, %3zu, %3zu): %c%c/%c%c (%f) [",
			target_a + 1, target_b + 1, target_a - start_target[0], target_b - start_target[1],
			g1_max < 2 ? xy_to_char[nuc1] : xy_to_char[nuc2],
			g1_max ? xy_to_char[nuc2] : xy_to_char[nuc1],
			g2_max < 2 ? xy_to_char[nuc1] : xy_to_char[nuc2],
			g2_max ? xy_to_char[nuc2] : xy_to_char[nuc1], mprob);
		if (final_out) {
			fprintf(final_out, "%s  %s  %4zu	%4zu",
					opt.ref_names[0], opt.ref_names[1], target_a + 1, target_b + 1);
			fprintf(final_out, "	%c%c/%c%c   %c%c	%c%c   %c   %c",
					g1_max < 2 ? xy_to_char[nuc1] : xy_to_char[nuc2],
					g1_max ? xy_to_char[nuc2] : xy_to_char[nuc1],
					g2_max < 2 ? xy_to_char[nuc1] : xy_to_char[nuc2],
					g2_max ? xy_to_char[nuc2] : xy_to_char[nuc1],
					g1_max < 2 ? xy_to_char[nuc1] : xy_to_char[nuc2],
					g1_max ? xy_to_char[nuc2] : xy_to_char[nuc1],
					g2_max < 2 ? xy_to_char[nuc1] : xy_to_char[nuc2],
					g2_max ? xy_to_char[nuc2] : xy_to_char[nuc1],
					nuc1, nuc2);
		}

		double prob_heterozygoteA = 0, prob_heterozygoteB = 0;
		double prob_heterozygote[N_FILES] = {0, 0};	/* [KSD,TODO] Assumes N_FILES == 2 */

		for (int g1 = 0; g1 <= 2; ++g1)
			for (int g2 = 0; g2 <= 2; ++g2) {
				if (final_out)
					fprintf(final_out, "   %f", gprob[g1*3 + g2]);
				fprintf(stderr, " %f", gprob[g1*3 + g2]);
				if (g1 == 1)
					prob_heterozygoteA += gprob[g1*3 + g2];
				if (g2 == 1)
					prob_heterozygoteB += gprob[g1*3 + g2];
			}
		if (final_out) {
			fprintf(final_out, "   %f   %f   %f   %f   %f   %f	%f  %f\n",
					gprob[0] + gprob[1] + gprob[2],
					gprob[3] + gprob[4] + gprob[5],
					gprob[6] + gprob[7] + gprob[8],
					gprob[0] + gprob[3] + gprob[6],
					gprob[1] + gprob[4] + gprob[7],
					gprob[2] + gprob[5] + gprob[8],
					sumA, sumB);
		}
		prob_heterozygote[0] = prob_heterozygoteA;
		prob_heterozygote[1] = prob_heterozygoteB;
		fprintf(stderr, "]%s\n", g1_max == 1 || g2_max == 1 ? "***"
			: abs(g1_max - g2_max) > 1 ? "+++" : "");
		fprintf(stderr, "prob_heterozygoteA %lf, prob_heterozygoteB %lf\n",
			prob_heterozygoteA, prob_heterozygoteB);
		int g_max[N_FILES] = {g1_max, g2_max};
		double ect_pvals[N_FILES] = {1, 1};
		if (opt.equal_homolog_coverage_test
			&& (g_max[0] == 1 || g_max[1] == 1))
				test_equal_homolog_coverage(mh, rf_info, ll,
					ref_base, covers, obs_nuc, obs_q,
					obs_rpos, g_max, nuc1, nuc2,
					ect_pvals);
		/* write out results to vcf files */
		
		for (int i = 0; i < N_FILES; ++i) {
			int one_alt = 0;
			
			if (!vcf_fp[i])
				continue;
			
			double pe = i
				? (1 - gprob[g_max[i]] - gprob[3 + g_max[i]] - gprob[6 + g_max[i]])
				: (1 - gprob[3 * g_max[i] + 2] - gprob[3 * g_max[i] + 1] - gprob[3 * g_max[i]]);
			
			if (i == 0)
				fprintf(vcf_fp[i], "%s", rf_info->name_A);
			else
				fprintf(vcf_fp[i], "%s", rf_info->name_B);
			fprintf(vcf_fp[i], "\t%lu\t.\t%c\t",
				ref_pos[i] + 1, iupac_to_char[ref_base[i]]);
			if (xy_to_iupac[nuc1] != ref_base[i]) {
				fprintf(vcf_fp[i], "%c", xy_to_char[nuc1]);
				one_alt = 1;
			}
			if (!no_alt_allele && nuc2 != nuc1
				&& xy_to_iupac[nuc2] != ref_base[i]) {
				if (one_alt)
					fputc(',', vcf_fp[i]);
				fprintf(vcf_fp[i], "%c", xy_to_char[nuc2]);
				one_alt = 1;
			}
			if (!one_alt)
				fputc('.', vcf_fp[i]);
			/* NOTE: Currently we do not output ALT alleles
			 * from the other subgenome.
			 */
			if (opt.posthoc_coverage_test && g_max[i] == 1 &&
				prob_heterozygote[i] < opt.min_genotype_post_prob) {
				int mgq = opt.min_genotype_post_prob < 1
					? MIN(99, (int) (-10 * log10(
					1 - opt.min_genotype_post_prob))) : 99;
				fprintf(vcf_fp[i], "\t.\tgq%d\t.\tGT:DP:GQ", mgq);
			} else {
				fprintf(vcf_fp[i], "\t.\tPASS\t.\tGT:DP:GQ");
			}
			if (opt.equal_homolog_coverage_test && g_max[i] == 1)
				fprintf(vcf_fp[i], ":ET");
			if (opt.vcf_opt->output_gl)
				fprintf(vcf_fp[i], ":GL");
			fputc('\t', vcf_fp[i]);

			/* reference allele is dominant allele */
			if (xy_to_iupac[nuc1] == ref_base[i]) {
				if (g_max[i] == 2)
					fprintf(vcf_fp[i], "1/1");
				else if (g_max[i] == 1)
					fprintf(vcf_fp[i], "0/1");
				else
					fprintf(vcf_fp[i], "0/0");
				
				/* reference allele is subdominant allele */
			} else if (xy_to_iupac[nuc2] == ref_base[i]) {
				if (g_max[i] == 2)
					fprintf(vcf_fp[i], "0/0");
				else if (g_max[i] == 1)
					fprintf(vcf_fp[i], "0/1");
				else
					fprintf(vcf_fp[i], "1/1");
				
				/* reference allele is neither of 2 dominant alleles */
			} else {
				if (g_max[i] == 2)
					fprintf(vcf_fp[i], "2/2");
				else if (g_max[i] == 1)
					fprintf(vcf_fp[i], "1/2");
				else
					fprintf(vcf_fp[i], "1/1");
			}
			fprintf(vcf_fp[i], ":%d:%d",
				(int) (ecoverage[i] + 0.5),
				pe > 0 ? MIN(99, (int) (-10 * log10(pe))) : 99);
			
			if (opt.equal_homolog_coverage_test && g_max[i] == 1)
				fprintf(vcf_fp[i], ":%.1f",
						fabs(log10(ect_pvals[i])));
			
			if (!opt.vcf_opt->output_gl) {
				fputc('\n', vcf_fp[i]);
				continue;
			}

			/* use profile log likelihood */
			if (xy_to_iupac[nuc1] == ref_base[i]) {
				if (!i)
					fprintf(vcf_fp[i], ":%.2f,%.2f,%.2f\n",
						lprob[0 + g_max[1]] / log(10),
						lprob[3 + g_max[1]] / log(10),
						lprob[6 + g_max[1]] / log(10));
				else
					fprintf(vcf_fp[i], ":%.2f,%.2f,%.2f\n",
						lprob[g_max[0]] / log(10),
						lprob[g_max[0] + 1] / log(10),
						lprob[g_max[0] + 2] / log(10));
			} else if (xy_to_iupac[nuc2] == ref_base[i]) {
				if (!i)
					fprintf(vcf_fp[i], ":%.2f,%.2f,%.2f\n",
						lprob[6 + g_max[1]] / log(10),
						lprob[3 + g_max[1]] / log(10),
						lprob[0 + g_max[1]] / log(10));
				else
					fprintf(vcf_fp[i], ":%.2f,%.2f,%.2f\n",
						lprob[g_max[0] + 2] / log(10),
						lprob[g_max[0] + 1] / log(10),
						lprob[g_max[0] + 0] / log(10));
			} else {
				if (!i)
					fprintf(vcf_fp[i], ":%.2f,%.2f,%.2f\n",
						lprob[0 + g_max[1]] / log(10),
						lprob[3 + g_max[1]] / log(10),
						lprob[6 + g_max[1]] / log(10));
				else
					fprintf(vcf_fp[i], ":%.2f,%.2f,%.2f\n",
						lprob[g_max[0]] / log(10),
						lprob[g_max[0] + 1] / log(10),
						lprob[g_max[0] + 2] / log(10));
			}
		}
		/* if heterozygous call with confidence record haplotype */
		if (opt.posthoc_coverage_test) {
			if (prob_heterozygote[0] >= opt.min_genotype_post_prob) {
				haplotype_posns[0][n_segregatingA] = target_a;
				hapA_prop[n_segregatingA] = ebaseA[modeA] / ecoverage[0];
				hapA_covg[n_segregatingA] = ecoverage[0];
				hapA_dom_nuc[n_segregatingA++] = modeA;
			}
			if (prob_heterozygote[1] >=  opt.min_genotype_post_prob) {
				haplotype_posns[1][n_segregatingB] = target_b;
				hapB_prop[n_segregatingB] = ebaseB[modeB] / ecoverage[1];
				hapB_covg[n_segregatingB] = ecoverage[1];
				hapB_dom_nuc[n_segregatingB++] = modeB;
			}
		}
	}
	if (final_out)
		fclose(final_out);

	for (unsigned int i = 0; i < N_FILES; ++i)
		if (vcf_fp[i])
			fclose(vcf_fp[i]);

	/* KSD,TODO Delete all this! */
	if (opt.posthoc_coverage_test) {
		
		unsigned int buffer_size = 32, buffer_block = 32, n_buffer = 1;
		unsigned int n_haplotype[2] = {0, 0};
		uint64_t *haplotype_id[2] = {
		malloc(buffer_size * sizeof *haplotype_id[0]),
		malloc(buffer_size * sizeof *haplotype_id[1])};
		unsigned int *haplotype_cnt[2] = {
		calloc(buffer_size, sizeof *haplotype_cnt[0]),
		calloc(buffer_size, sizeof *haplotype_cnt[1])};
		/* [KSD, TODO] check allocations */
		
		/* count haplotype occurrences in both genomes */
		n_read = 0;
		for (merge_hash *me = mh; me != NULL; me = me->hh.next) {
			if (me->exclude)
				continue;

//			fprintf(stderr, "Processing read %zu", n_read);
			/* find subgenomic assignment or continue if ambiguously aligned */
			unsigned int sgenome = 0;
			for (unsigned int j = 1; j < N_FILES; ++j)
				if (pp[j][n_read] >= opt.min_alignment_post_prob) {
					sgenome = j;
					break;
				}
			if (!sgenome && pp[sgenome][n_read]
						< opt.min_alignment_post_prob) {
				++n_read;
				continue;
			}
//			fprintf(stderr, " assigned to subgenome %s.", sgenome?"B":"A");

			/* get haplotype: packed 2-bit nucleotides */
			sam_entry *se = &sds[sgenome]->se[me->indices[sgenome][0]];
			uint64_t id = get_haplotype_id(se,
				sgenome ? haplotype_posns[1] : haplotype_posns[0],
				sgenome ? n_segregatingB : n_segregatingA);
			//fprintf(stderr, "id = %lu\n", id);

			/* increase count for this haplotype or ... */
			unsigned int exists = 0, hpos = 0;
			for (unsigned int j = 0; j < n_haplotype[sgenome]; ++j)
				if (haplotype_id[sgenome][j] == id) {
					++haplotype_cnt[sgenome][j];
					hpos = j;
					exists = 1;
					break;
				}
			/* add new haplotype */
			if (!exists) {
				if (n_haplotype[sgenome] == buffer_size) {
					buffer_size = ++n_buffer * buffer_block;
					uint64_t *new_id = realloc(haplotype_id[0], buffer_size * sizeof *haplotype_id[0]);
					if (!new_id)
						return mmessage(ERROR_MSG, INTERNAL_ERROR, "Ran out of memory!\n");
					haplotype_id[0] = new_id;
					new_id = realloc(haplotype_id[1], buffer_size * sizeof *haplotype_id[1]);
					if (!new_id)
						return mmessage(ERROR_MSG, INTERNAL_ERROR, "Ran out of memory!\n");
					haplotype_id[1] = new_id;
					unsigned int *new_cnt = realloc(haplotype_cnt[0], buffer_size * sizeof *haplotype_cnt[0]);
					if (!new_cnt)
						return mmessage(ERROR_MSG, INTERNAL_ERROR, "Ran out of memory!\n");
					memset(new_cnt + (n_buffer - 1) * buffer_block, 0, buffer_block * sizeof *new_cnt);
					haplotype_cnt[0] = new_cnt;
					new_cnt = realloc(haplotype_cnt[1], buffer_size * sizeof *haplotype_cnt[1]);
					if (!new_cnt)
						return mmessage(ERROR_MSG, INTERNAL_ERROR, "Ran out of memory!\n");
					memset(new_cnt + (n_buffer - 1) * buffer_block, 0, buffer_block * sizeof *new_cnt);
					haplotype_cnt[1] = new_cnt;
				}
				hpos = n_haplotype[sgenome]++;
				haplotype_id[sgenome][hpos] = id;
				haplotype_cnt[sgenome][hpos] = 1;
			}
//			fprintf(stderr, "Found haplotype %lu (%u)\n", id, haplotype_cnt[sgenome][hpos]);
			++n_read;
		}

		/* for each subgenome, test hypothesis of equal coverage of two
		 * most common haplotypes, assuming there are two haplotypes
		 */
		
		int *fail_test[2] = {NULL, NULL};
		for (unsigned int i = 0; i < 2; ++i) {
			unsigned int max = 0, smax = 0, total = 0;
			unsigned int idx = 0;

			if (i ? n_segregatingB : n_segregatingA)
				fail_test[i] = calloc(i ? n_segregatingB
					: n_segregatingA, sizeof *fail_test[i]);

			for (unsigned int j = 0; j < n_haplotype[i]; ++j) {
				if (haplotype_cnt[i][j] > max) {
					smax = max;
					max = haplotype_cnt[i][j];
					idx = j;
				} else if (haplotype_cnt[i][j] > smax) {
					smax = haplotype_cnt[i][j];
				}
				total += haplotype_cnt[i][j];
			}
			uint64_t id = haplotype_id[i][idx];
			unsigned int n = max + smax;
			/* Wilson score interval with continuity correction */
			double phat = (double) max / n;
			double z = 1.959963984540054;
			double z2 = z*z;
			double sq = sqrt(z2 - 1./n + 4*n*phat*(1-phat)
								+ 4*phat - 2);
			double den = 2*(n+z2);
			double wminus = (2*n*phat + z2 - (z*sq + 1))/den;
			double wplus = (2*n*phat + z2 + (z*sq + 1))/den;
			if (wminus < 0) wminus = 0;
			if (wplus > 1) wplus = 1;
			unsigned int covers_p5 = 0;
			if (wminus <= 0.5 && wplus >= 0.5)
				covers_p5 = 1;

			double x2 = 2*n*(((double) max / n - 0.5) * ((double) max / n - 0.5)
				 + ((double) smax / n - 0.5) * ((double) smax / n - 0.5));
			double pval = pchisq(x2, 1, 0, 0);
			if (i ? n_segregatingB : n_segregatingA) {

				debug_msg(1, 1, "Subgenome %s haplotype ",
								i?"B":"A");
				for (unsigned int j = 0; j < (i ? n_segregatingB
							: n_segregatingA); ++j)
					debug_msg_cont(1, 1, "%c",
						xy_to_char[3U & (id >> (2*j))]);
				debug_msg_cont(1, 1, " has coverage of %u/%u "
					" (%u : %u of total %u) high-confidence"
					" reads (X2 = %f; p-value %e)\n",
					haplotype_cnt[i][idx],
					n, max, smax, total, x2, pval);
				debug_msg(1, 1, "Subgenome %s coverage "
					"proportion confidence interval: %f "
					"%f (%u)\n", i?"B":"A", wminus, wplus,
								covers_p5);

				debug_msg(1, 1, "Subgenome %s modal alleles ",
								i?"B":"A");
				for (unsigned int j = 0; j < (i ? n_segregatingB
							: n_segregatingA); ++j)
					debug_msg_cont(1, 1, "%c", xy_to_char[
						i ? hapB_dom_nuc[j]
							: hapA_dom_nuc[j]]);
				debug_msg_cont(1, 1, " coverage:");
				for (unsigned int j = 0; j < (i ? n_segregatingB
						: n_segregatingA); ++j) {
					x2 = 2 * (i ? hapB_covg[j] : hapA_covg[j])
						* ( ((i ? hapB_prop[j] : hapA_prop[j]) - 0.5) * ((i ? hapB_prop[j] : hapA_prop[j]) - 0.5)
						+ ( -(i ? hapB_prop[j] : hapA_prop[j]) + 0.5) * ( -(i ? hapB_prop[j] : hapA_prop[j]) + 0.5));
					pval = pchisq(x2, 1, 0, 0);
					if (pval < 0.05)
						fail_test[i][j] = 1;
					debug_msg_cont(1, 1, " %f of %f (%e)",
						i ?  hapB_prop[j] : hapA_prop[j],
						i ? hapB_covg[j] : hapA_covg[j],
									pval);
				}
				debug_msg_cont(1, 1, "\n");
			}
		}
		update_vcf(&opt, fail_test, haplotype_posns);

		if (fail_test[0])
			free(fail_test[0]);
		if (fail_test[1])
			free(fail_test[1]);
		free(haplotype_id[0]);
		free(haplotype_id[1]);
		free(haplotype_cnt[0]);
		free(haplotype_cnt[1]);
	}

	for (unsigned int j = 0; j < N_FILES; ++j)
		if (pp[j])
			free(pp[j]);
	if (pp)
		free(pp);
	if (obs_nuc)
		free(obs_nuc);
	if (obs_q)
		free(obs_q);
	if (genome_src)
		free(genome_src);
	if (covers)
		free(covers);
	if (rd_idxA)
		free(rd_idxA);
	if (opt.ampliclust_command)
		free_input(in);
} /* main */

uint64_t get_haplotype_id(sam_entry *se, size_t *haplotype,
						unsigned int n_segregating)
{
	size_t rf_index = se->pos - 1;
	unsigned int id = 0;
	unsigned int n_ash = 0;
	unsigned int rd_index = 0;

	for (unsigned int i = 0; i < n_segregating; ++i) {
		size_t next_rpos = haplotype[i];
		//fprintf(stderr, "\nLooking for position %zu (%u) ", next_rpos, n_ash);
		while (n_ash < se->cig->n_ashes) {
			if (se->cig->ashes[n_ash].type == CIGAR_INSERTION
				|| se->cig->ashes[n_ash].type == CIGAR_SOFT_CLIP) {
				rd_index += se->cig->ashes[n_ash].len;
				++n_ash;
				continue;
			} else if (se->cig->ashes[n_ash].type == CIGAR_DELETION) {
				rf_index += se->cig->ashes[n_ash].len;
				++n_ash;
				continue;
			} else if (se->cig->ashes[n_ash].type != CIGAR_MATCH
				&& se->cig->ashes[n_ash].type != CIGAR_MISMATCH
				&& se->cig->ashes[n_ash].type != CIGAR_MMATCH) {
				++n_ash;
				continue;
			}
			//fprintf(stderr, "[%zu, %zu)", rf_index, rf_index + se->cig->ashes[n_ash].len);
			if (rf_index + se->cig->ashes[n_ash].len > next_rpos
				&& rf_index <= next_rpos) {
				//fprintf(stderr, "%c", xy_to_char[get_nuc(se->read, XY_ENCODING, rd_index + next_rpos - rf_index)]);
				id |= get_nuc(se->read, XY_ENCODING,
					rd_index + next_rpos - rf_index) << (i*2);
				break;	/* while */
			}
			rd_index += se->cig->ashes[n_ash].len;
			rf_index += se->cig->ashes[n_ash].len;
			++n_ash;
		}
	}
	//fprintf(stderr, "\n");

	return id;
} /* get_haplotype_id */


/**
 * Log likelihood of alignment.  Compute log likelihood of alignment
 * assuming quality scores are literal and all substitutions equally
 * likely.
 *
 * @param se		alignment entry from sam file (xy_t)
 * @param rd_id		index of read
 * @param ref		reference sequence (iupac_t)
 * @param vptr		mlogit stuff
 * @param in_show	show alignments
 * @param start_rf	starting position of alignment relative to extracted 
 *			target regions from references
 * @param in_debug	debugging
 * @return		log likelihood
 */
double ll_align(sam_entry *se, unsigned int rd_id, unsigned char *ref,
		mlogit_stuff *mls, unsigned char *in_show, size_t start_rf,
		int in_debug)
{

	size_t rf_index = start_rf;		/* starting reference position */
	unsigned int rd_index = 0;		/* starting position in read */

/* for debugging: output selected alignments */
/*
	size_t rf_pos1 = 296, rf_pos2 = 402;
	char nuc1 = 'G', nuc2 = 'C';
	unsigned int flag1 = 0, flag2 = 0;

	for (unsigned int i = 0; i < se->cig->n_ashes; ++i) {

		if (se->cig->ashes[i].type == CIGAR_DELETION) {
			rf_index += se->cig->ashes[i].len;
			continue;
		} else if (se->cig->ashes[i].type == CIGAR_SOFT_CLIP) {
			rd_index += se->cig->ashes[i].len;
			continue;
		} else if (se->cig->ashes[i].type == CIGAR_INSERTION) {
			rd_index += se->cig->ashes[i].len;
			continue;
		} else if (se->cig->ashes[i].type == CIGAR_HARD_CLIP) {
			continue;
		} else if (se->cig->ashes[i].type != CIGAR_MATCH
			   && se->cig->ashes[i].type != CIGAR_MISMATCH
			   && se->cig->ashes[i].type != CIGAR_MMATCH) {
			continue;
		}
		if (rf_index + se->cig->ashes[i].len > rf_pos1
			&& rf_index < rf_pos1
			&& xy_to_char[get_nuc(se->read, XY_ENCODING,
				rd_index + rf_pos1 - rf_index)] == nuc1) {
			flag1 = 1;
		}
		if (rf_index + se->cig->ashes[i].len > rf_pos2
			&& rf_index < rf_pos2
			&& xy_to_char[get_nuc(se->read, XY_ENCODING,
				rd_index + rf_pos2 - rf_index)] == nuc2) {
			flag2 = 1;
		}
		rf_index += se->cig->ashes[i].len;
		rd_index += se->cig->ashes[i].len;
	}
	fprintf(stderr, "flag1=%u, flag2=%u\n", flag1, flag2);
	if (flag1)// && flag2)
		*in_show = 1;
	rf_index = se->pos - 1;
	rd_index = 0;
 */

/* end debugging code */

	unsigned char show = *in_show;

	int fxn_debug = in_debug;//ABSOLUTE_SILENCE;//DEBUG_II;//DEBUG_III;//DEBUG_II;//
	/* control display during debugging */
	int display_reverse_complement = 1;	/* show rc if so aligned */
	int display_dot = 1;			/* display dot for match */
	int guide_posts = (fxn_debug || show) && 1;
	/* display mark every 10 nucs */
	int display_after = (fxn_debug || show) && 1;
	/* do not draw as process */
	int display_qual = (fxn_debug || show) && 1;
	/* display quality scores */


	double ll = 0;				/* initial alignment log likelihood */
	size_t align_len = 0;			/* alignment length */
	int reverse_complement = display_reverse_complement && se->flag >> 4 & 1;
	iupac_t *align_display = NULL;
	qual_t *qual_display = NULL;
	size_t align_index = 0;

	for (unsigned int i = 0; i < se->cig->n_ashes; ++i)
		if (se->cig->ashes[i].type == CIGAR_DELETION
			|| se->cig->ashes[i].type == CIGAR_INSERTION
			|| se->cig->ashes[i].type == CIGAR_MATCH
			|| se->cig->ashes[i].type == CIGAR_MMATCH
			|| se->cig->ashes[i].type == CIGAR_MISMATCH)
			align_len += se->cig->ashes[i].len;

	if (reverse_complement && (fxn_debug >= DEBUG_II || show)
		&& !display_qual && !display_after)
		align_display = malloc(align_len * sizeof *align_display);
	if (display_after)
		align_display = malloc(align_len * sizeof *align_display);
	if (display_qual)
		qual_display = malloc(align_len * sizeof *qual_display);

	debug_msg(fxn_debug >= DEBUG_II || show, fxn_debug, "Read = %u, "
		  "Name = %s, Flag = %u, Pos = %u, Align. len = %u, Cigar = ",
		  rd_id, se->name, se->flag, se->pos, align_len);
	for (unsigned int i = 0; i < se->cig->n_ashes; ++i)
		debug_call(fxn_debug >= DEBUG_II || show, fxn_debug,
			fprintf(stderr, "%u%c", se->cig->ashes[i].len,
				cigar_char[se->cig->ashes[i].type]));
	debug_msg_cont(fxn_debug >= DEBUG_II || show, fxn_debug, "\n");
	debug_msg(fxn_debug >= DEBUG_II || show, fxn_debug, "Ref : ");

	unsigned int out = 0;
	for (unsigned int i = 0; i < se->cig->n_ashes; ++i) {

		/* march through alignment; optimally output reference */
		if (se->cig->ashes[i].type == CIGAR_DELETION) {
			if (fxn_debug >= DEBUG_II || show || (fxn_debug && fxn_debug <= global_debug_level))
				for (size_t j = 0; j < se->cig->ashes[i].len; ++j) {
					if (!reverse_complement && !display_after)
						fprintf(stderr, "%c",
							iupac_to_char[ref[rf_index + j]]);
					else if (!display_qual && (fxn_debug || show))
						align_display[align_index++] = ref[rf_index + j];
					if (!reverse_complement && display_qual) {
						align_display[align_index] = ref[rf_index + j];
						qual_display[align_index++] = 0;
					}
					if (!reverse_complement && !display_after && guide_posts && !(++out % 10))
						fprintf(stderr, "|");
				}
			rf_index += se->cig->ashes[i].len;
			continue;
		} else if (se->cig->ashes[i].type == CIGAR_SOFT_CLIP) {
			if (fxn_debug >= DEBUG_II || show || (fxn_debug && fxn_debug <= global_debug_level))
				for (size_t j = 0; j < se->cig->ashes[i].len; ++j) {
					if (!reverse_complement && !display_after)
						fputc('.', stderr);
					if (!reverse_complement && !display_after && guide_posts && !(++out % 10))
						fprintf(stderr, "|");
				}
			rd_index += se->cig->ashes[i].len;
			continue;
		} else if (se->cig->ashes[i].type == CIGAR_INSERTION) {
			if (fxn_debug >= DEBUG_II || show || (fxn_debug && fxn_debug <= global_debug_level))
				for (size_t j = 0; j < se->cig->ashes[i].len; ++j) {
					if (!reverse_complement && !display_after)
						fputc('-', stderr);
					else if (!display_qual)
						align_display[align_index++] = 0;
					if (display_qual) {
						align_display[align_index] = 0;
						qual_display[align_index++] =
						(char) get_qual(se->qual,
								rd_index + j) +
						MIN_ASCII_QUALITY_SCORE;
					}
					if (!reverse_complement && !display_after && guide_posts && !(++out % 10))
						fprintf(stderr, "|");
				}
			rd_index += se->cig->ashes[i].len;
			continue;
		} else if (se->cig->ashes[i].type == CIGAR_HARD_CLIP) {
			continue;
		} else if (se->cig->ashes[i].type != CIGAR_MATCH
			   && se->cig->ashes[i].type != CIGAR_MISMATCH
			   && se->cig->ashes[i].type != CIGAR_MMATCH) {
			mmessage(WARNING_MSG, NO_ERROR,
				 "Unhandled cigar string: %c.\n",
				 cigar_char[se->cig->ashes[i].type]);
			continue;
		}

		for (size_t j = 0; j < se->cig->ashes[i].len; ++j) {
			mls->pos = rd_index + j;
			double llt = sub_prob_given_q_with_encoding(ref[rf_index + j],
									get_nuc(se->read, XY_ENCODING, rd_index + j),
									IUPAC_ENCODING, XY_ENCODING,
									get_qual(se->qual, rd_index + j), 1, (void *) mls);
			ll += llt;
			debug_msg(fxn_debug >= DEBUG_III, fxn_debug, "%u (%u): %c -> %c (%c): %f (%f)\n", rf_index + j, j, iupac_to_char[ref[rf_index + j]], xy_to_char[get_nuc(se->read, XY_ENCODING, rd_index + j)], (char)get_qual(se->qual, rd_index + j) + MIN_ASCII_QUALITY_SCORE, llt, ll);
			if (!reverse_complement && !display_after
				&& (fxn_debug || show))
				debug_msg_cont(fxn_debug >= DEBUG_I || show, fxn_debug,
						   "%c", iupac_to_char[ref[rf_index + j]]);
			else if (!display_qual && (fxn_debug || show))
				align_display[align_index++] = ref[rf_index + j];
			if (display_qual) {
				align_display[align_index] = ref[rf_index + j];
				qual_display[align_index++] = (char) get_qual(
										  se->qual, rd_index + j)
				+ MIN_ASCII_QUALITY_SCORE;
			}
			if (!reverse_complement && !display_after && guide_posts && !(++out % 10))
				debug_call(fxn_debug >= DEBUG_I || show, fxn_debug, fprintf(stderr, "|"));
		}
		rf_index += se->cig->ashes[i].len;
		rd_index += se->cig->ashes[i].len;
	}

	if ((fxn_debug >= DEBUG_II || show || (fxn_debug && fxn_debug <= global_debug_level)) && reverse_complement) {
		for (size_t j = align_len; j > 0; --j) {
			fputc(!align_display[j - 1] ? '-' :
				  iupac_to_char[iupac_to_rc[
							align_display[j - 1]]], stderr);
			if (!((align_len - j + 1) % 10) && guide_posts)
				fputc('|', stderr);
		}
		fprintf(stderr, "\n");
		if (qual_display) {
			debug_msg(fxn_debug >= DEBUG_II || show, fxn_debug, "Qual: ");
			for (size_t j = align_len; j > 0; --j) {
				fputc(!qual_display[j - 1]
					  ? ' ' : qual_display[j - 1], stderr);
				if (!((align_len - j + 1) % 10) && guide_posts)
					fputc('|', stderr);
			}
			fprintf(stderr, "\n");
		}
	}
	if (display_after && !reverse_complement && (fxn_debug >= DEBUG_II || show || (fxn_debug && fxn_debug <= global_debug_level))) {
		for (size_t j = 0; j < align_len; ++j) {
			fputc(!align_display[j] ? '-' :
				  iupac_to_char[align_display[j]], stderr);
			if (!((j + 1) % 10) && guide_posts)
				fputc('|', stderr);
		}
		fprintf(stderr, "\n");
	}
	if (qual_display && !reverse_complement && (fxn_debug >= DEBUG_II || show || (fxn_debug && fxn_debug <= global_debug_level))) {
		debug_msg(fxn_debug >= DEBUG_II || show, fxn_debug, "Qual: ");
		for (size_t j = 0; j < align_len; ++j) {
			fputc(!qual_display[j] ? ' ' : qual_display[j], stderr);
			if (!((j + 1) % 10) && guide_posts)
				fputc('|', stderr);
		}
		fprintf(stderr, "\n");
	}

	/* optionally output read to show alignment */
	if (fxn_debug >= DEBUG_II || show || (fxn_debug && fxn_debug <= global_debug_level)) {
		align_index = 0;
		debug_msg(fxn_debug >= DEBUG_II || show, fxn_debug, "Read: ");
		rf_index = start_rf;
		rd_index = 0;
		out = 0;
		for (unsigned int i = 0; i < se->cig->n_ashes; ++i) {
//			flag1 = flag2 = 0;
			if (se->cig->ashes[i].type == CIGAR_HARD_CLIP) {
				//				rd_index += se->cig->ashes[i].len;
			} else if (se->cig->ashes[i].type == CIGAR_DELETION) {
				for (size_t j = 0; j < se->cig->ashes[i].len;
					 ++j) {
					if (!reverse_complement && !display_after)
						fputc('-', stderr);
					else
						align_display[align_index++] = 0;
					if (!reverse_complement && !display_after && guide_posts && !(++out % 10))
						fprintf(stderr, "|");
				}
				rf_index += se->cig->ashes[i].len;
			} else if (se->cig->ashes[i].type == CIGAR_SOFT_CLIP) {
				if (fxn_debug >= DEBUG_II || show) {
					if (!reverse_complement && !display_after)
						fwrite_nuc_segment(stderr,
								   se->read, XY_ENCODING,
								   rd_index, rd_index
								   + se->cig->ashes[i].len);
					out += se->cig->ashes[i].len;
					if (!reverse_complement && !display_after && guide_posts && !(out % 10))
						fprintf(stderr, "|");
				}
				rd_index += se->cig->ashes[i].len;
			} else if (se->cig->ashes[i].type == CIGAR_INSERTION) {
				if (fxn_debug >= DEBUG_II || show) {
					if (!reverse_complement && !display_after)
						fwrite_nuc_segment(stderr, se->read,
							   XY_ENCODING, rd_index, rd_index
							   + se->cig->ashes[i].len);
					else
						for (size_t j = 0; j < se->cig->ashes[i].len; ++j)
							align_display[align_index++] = xy_to_iupac[get_nuc(se->read, XY_ENCODING, rd_index + j)];
					out += se->cig->ashes[i].len;
					if (!reverse_complement && !display_after && guide_posts && !(out % 10))
						fprintf(stderr, "|");
				}
				rd_index += se->cig->ashes[i].len;
			} else if (se->cig->ashes[i].type == CIGAR_MATCH
				   || se->cig->ashes[i].type == CIGAR_MMATCH
				   || se->cig->ashes[i].type == CIGAR_MISMATCH) {
/*
if (show)
	fprintf(stderr, " rf_index=%zu -> %zu > %zu", rf_index, rf_index + se->cig->ashes[i].len, rf_pos1);
if (rf_index + se->cig->ashes[i].len > rf_pos1) {
	fprintf(stderr, "Site %zu: %c\n", rf_pos1, xy_to_char[get_nuc(se->read, XY_ENCODING, rd_index + rf_pos1 - rf_index)]);
	flag1 = 1;
}
if (rf_index + se->cig->ashes[i].len > rf_pos2) {
	fprintf(stderr, "Site %zu: %c\n", rf_pos2, xy_to_char[get_nuc(se->read, XY_ENCODING, rd_index + rf_pos2 - rf_index)]);
	flag2 = 1;
}
*/
				for (size_t j = 0; j < se->cig->ashes[i].len;
					 ++j) {
					data_t nuc = get_nuc(se->read,
								 XY_ENCODING, rd_index + j);
//if (flag1 || flag2)
//if (show) fprintf(stderr, " %zu=%c", rd_index + j, xy_to_char[nuc]);
					if (!reverse_complement && !display_after
						&& (fxn_debug || show))
						fprintf(stderr, "%c", display_dot && iupac_to_xy[
												 ref[rf_index + j]] == nuc
							? '.' : xy_to_char[nuc]);
					else if (display_after)
						align_display[align_index++] = display_dot && ref[rf_index + j] == xy_to_iupac[nuc] ? 15 : xy_to_iupac[nuc];
					if (!reverse_complement && !display_after && guide_posts && !(++out % 10))
						fprintf(stderr, "|");
				}
				rd_index += se->cig->ashes[i].len;
				rf_index += se->cig->ashes[i].len;
			}
		}
		if (reverse_complement && (fxn_debug || show || (fxn_debug && fxn_debug <= global_debug_level))) {
			for (size_t j = align_len; j > 0; --j) {
				fprintf(stderr, "%c", !align_display[j - 1] ? '-' : align_display[j - 1] == 15 ? '.' : iupac_to_char[iupac_to_rc[align_display[j - 1]]]);
				if (!((align_len - j + 1) % 10) && guide_posts)
					fputc('|', stderr);
			}
			fprintf(stderr, "\n");
		}
		if (!reverse_complement && display_after) {
			for (size_t j = 0; j < align_len; ++j) {
				fprintf(stderr, "%c", !align_display[j] ? '-' : align_display[j] == 15 ? '.' : iupac_to_char[align_display[j]]);
				if (!((j + 1) % 10) && guide_posts)
					fputc('|', stderr);
			}
			fprintf(stderr, "\n");
		}
	}

	if (align_display)
		free(align_display);
	if (qual_display)
		free(qual_display);

	return ll;
} /* ll_align */

/**
 * Read the output file produced by ampliclust.  The input fastq file to ampliCI
 * excludes all reads already screened out by capg.c, including those marked
 * as excluded in the merge hash.  This code is recording indices in the merge
 * hash that should be excluded because of closeness to low abundance
 * haplotypes.
 *
 * This definition and that of make_input() and free_input() may belong
 * in a different file, but they do not belong in sam.c, which is focused
 * on parsing and storing data from sam files.
 *
 * @param fp	open handle of ampliclust output file
 * @param in	store data in this structure
 * @return error status
 */
int read_ampliclust_results(FILE *fp, input *in)
{
	char c;
	unsigned int i, k, count = 0;

	MAKE_1ARRAY(in->assignment, in->n_observation);

	/* get read assignments (indices of assigned haplotype) and
	 * haplotype estimated true abundances
	 */
	while (!feof(fp)) {
		c = fgetc(fp);
		if (c == 'K') {
			fgetc(fp);	/* discard ':' */
			fscanf(fp, "%u", &in->K);
		} else if (c == 'a') {
			c = fgetc(fp);
			if (c == 's') {	/* assignments! */
				while (c != ':')
					c = fgetc(fp);
				for (i = 0; i < in->n_observation; ++i)
					fscanf(fp, "%u", &in->assignment[i]);
			}
		}
	}

	fclose(fp);
	return NO_ERROR;

	/* record indices of haplotypes to exclude by low estimated abundance */
	/* Suppose coverage of subgenome A is w times that of subgenome B.
	 * Homeoelogous SNPs only: w/(1+w), 1/(1+w)
	 * Homeoelogous SNPs, subgenomic SNP segregating in one subgenome only: w/(1+w), 0.5/(1+w), 0.5/(1+w)
	 * Homeoelogous SNPs, subgenomic SNPs segregating in both subgenomes: 0.5w/(1+w), 0,5w/(1+w), 0.5/(1+w), 0.5/(1+w)
	 * in->proptest is number of hypotheses rejected, tested last to first
	 */

	CMAKE_1ARRAY(in->exclude, in->K);
	for (i = 4 - in->proptest; i < in->K; ++i)
		in->exclude[count++] = i;

	/* record merge hash indices of reads to exclude */
	for (i = 0; i < in->n_observation; ++i) {
		for (k = 0; k < count; ++k)
			if (in->assignment[i] == in->exclude[k]) {
				in->n_excluded_id++;
				break;
			}
	}
	mmessage(INFO_MSG, NO_ERROR, "Excluding %zu reads from %u haplotypes "
		"(%zu kept).\n", in->n_excluded_id, count, 4 - in->proptest,
		in->n_observation - in->n_excluded_id);

	MAKE_1ARRAY(in->exclude_id, in->n_excluded_id);

	in->n_excluded_id = 0;
	for (i = 0; i < in->n_observation; ++i)
		for (k = 0; k < count; ++k)
			if (in->assignment[i] == in->exclude[k]) {
				in->exclude_id[in->n_excluded_id++] = i;
				break;
			}

	return NO_ERROR;
} /* read_ampliclust_results */

/**
 * Make input data for results given by ampliclust.
 *
 * @param in	pointer to input structure
 * @return error status
 */
int make_input(input **in) {

	*in = malloc(sizeof **in);

	if (*in == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "input object");

	(*in)->assignment = NULL;
	(*in)->exclude = NULL;
	(*in)->exclude_id = NULL;
	(*in)->not_input = NULL;
	(*in)->n_excluded = 0;
	(*in)->n_excluded_id = 0;
	(*in)->n_hash_excluded = 0;
	(*in)->K = 0;
	(*in)->proptest = 0;
	(*in)->n_prop_exclude = 0;
	(*in)->n_observation = 0;

	return NO_ERROR;
}/* make_input */

/**
 * Free input.
 *
 * @param in input structure
 */
void free_input(input *in) {

	if (!in)
		return;

	free(in->assignment);
	free(in->exclude);
	free(in->exclude_id);
	free(in->not_input);

	free(in);
}/* free_input */

int default_options(options *opt)
{
	opt->display_alignment = 0;
	opt->drop_unmapped = 1;
	opt->drop_secondary = 1;
	opt->drop_soft_clipped = UINT_MAX;
	opt->drop_indel = UINT_MAX;
	opt->proptest_screen = 0;
	opt->weight_penalty = 1;
	opt->min_length = 0;
	opt->max_length = 0;
	opt->min_log_likelihood = -INFINITY;
	opt->n_sample = 100;
	opt->max_eerr = INFINITY;
	opt->param_file = NULL;
	opt->error_file = NULL;
	opt->biallelic_screen = 0.5;
//	opt->ampliclust_file = NULL;
	opt->max_quality_score = INT_MAX;
	opt->use_bam = 0;
	for (int i = 0; i < N_FILES; ++i) {
		opt->sbam_files[i] = NULL;
		opt->fsa_files[i] = NULL;
		opt->ref_names[i] = NULL;
		opt->vcf_files[i] = NULL;
	}
	opt->output_file = NULL;
	opt->sample_name = NULL;
	opt->extracted_rf = "extracted";
	opt->sam_file = NULL;
	opt->ampliclust_command = NULL;
	opt->ac_fastq_file = "union.fastq";
	opt->ac_outfile = "amplici";
	opt->ac_low_bound = 1.5;
	opt->coverage_file = "coverage.txt";
	opt->write_fastq_and_quit = 0;
	opt->min_expected_coverage = 5.;
	opt->posthoc_coverage_test = 0;
	opt->equal_homolog_coverage_test = 0;
	opt->min_alignment_post_prob = 0.99;
	opt->min_genotype_post_prob = 0.99;
	make_default_vcf_options(&opt->vcf_opt);
	
	return NO_ERROR;
} /* default_options */

int parse_options_capg(options *opt, int argc, const char **argv)
{
	int i, j;
	int err = NO_ERROR;
	char a;

	for (i = 1; i < argc; ++i) {
		if (strlen(argv[i]) < 2)
			usage_error(argv, i, (void *)opt);
		j = 1;
		a = argv[i][j];
		while (a == '-' && ++j < (int) strlen(argv[i]))
			a = argv[i][j];
		switch(a) {

		/* cases within switch not indented */
//		case 'a':
//			if (!strncmp(&argv[i][j], "ampliclust_o", 12)
//				|| !strncmp(&argv[i][j], "ampliclust-o", 12)) {
//				opt->ac_outfile = argv[++i];
//				mmessage(INFO_MSG, NO_ERROR, "Ampliclust "
//					"output file base name: '%s'\n",
//						opt->ac_outfile);
//			} else if (!strncmp(&argv[i][j], "ampliclust_f", 12)
//				|| !strncmp(&argv[i][j], "ampliclust-f", 12)) {
//				opt->ac_fastq_file = argv[++i];
//				mmessage(INFO_MSG, NO_ERROR, "Ampliclust "
//					"fastq filename: '%s'\n",
//						opt->ac_fastq_file);
//			} else if (!strncmp(&argv[i][j], "ampliclust_l", 12)
//				|| !strncmp(&argv[i][j], "ampliclust-l", 12)) {
//				opt->ac_low_bound = read_cmdline_double(argc,
//					argv, ++i, opt);
//				mmessage(INFO_MSG, NO_ERROR, "Ampliclust low"
//					 "er bound: %f\n", opt->ac_low_bound);
//			} else {
//				opt->ampliclust_command = argv[++i];
//				mmessage(INFO_MSG, NO_ERROR, "Ampliclust "
//							"command: '%s'\n",
//						opt->ampliclust_command);
//			}
//			break;
		case 'b':
			if (!strncmp(&argv[i][j], "bam", 3)) {
				if (i + N_FILES >= argc) {
					err = mmessage(ERROR_MSG,
							   INVALID_CMD_ARGUMENT, "Too few "
							   "arguments to --bam_files "
							   "command-line option.\n");
					goto CMDLINE_ERROR;
				}
				opt->use_bam = 1;
				mmessage(INFO_MSG, NO_ERROR, "BAM files:");
				for (j = 0; j < N_FILES; ++j) {
					opt->sbam_files[j] = argv[++i];
					fprintf(stderr, " %s",
						opt->sbam_files[j]);
				}
				fprintf(stderr, "\n");
			} else if (!strncmp(&argv[i][j], "bia", 3)) {
				opt->biallelic_screen = read_cmdline_double(
										argc, argv, ++i, opt);
				mmessage(INFO_MSG, NO_ERROR, "Dropping sites "
					 "with third allele above %f of minimum "
					 "estimated subgenomic coverage.\n",
					 opt->biallelic_screen);
			} else {
				goto CMDLINE_ERROR;
			}
			break;
		case 'd':
			opt->display_alignment = !opt->display_alignment;
			mmessage(INFO_MSG, NO_ERROR, "Display alignments on "
				 "stderr? %s\n", opt->display_alignment
				 ? "yes" : "no");
			break;
		case 'e':
			if (!strncmp(&argv[i][j], "ex", 2)) {
				opt->max_eerr = read_cmdline_double(
							argc, argv, ++i, opt);
				mmessage(INFO_MSG, NO_ERROR, "Dropping reads "
					 "with more than %f expected errors.\n",
								 opt->max_eerr);
			} else if (!strncmp(&argv[i][j], "eq", 2)) {
				opt->equal_homolog_coverage_test = 1;
				mmessage(INFO_MSG, NO_ERROR, "Will run "
						"equal coverage test.\n");
				/* legacy option: you want? try -po
				opt->posthoc_coverage_test = 1;
				mmessage(INFO_MSG, NO_ERROR, "Performing "
					"post-hoc confidence interval of equal "
					"homologous chromosome coverage\n");*/
			} else if (!strncmp(&argv[i][j], "error_d", 7)
				   || !strncmp(&argv[i][j], "error-d", 7)) {
				opt->error_file = fopen(argv[++i], "w");
				if (!opt->error_file) {
					err = mmessage(ERROR_MSG,
							   INVALID_CMD_ARGUMENT,
							   "Could not open file '%s'\n",
							   argv[i]);
					goto CMDLINE_ERROR;
				}
				mmessage(INFO_MSG, NO_ERROR, "Error data file: "
					 "'%s'\n", argv[i]);
			} else if (!strncmp(&argv[i][j], "er", 2)) {
				opt->param_file = argv[++i];
				if (access(opt->param_file, F_OK) == -1) {
					err = mmessage(ERROR_MSG,
							   INVALID_CMD_ARGUMENT,
							   "Could not open file '%s'.\n",
							   opt->param_file);
					goto CMDLINE_ERROR;
				}
				mmessage(INFO_MSG, NO_ERROR, "Error "
					 "probabilities file: '%s'\n",
					 opt->param_file);
			}
			break;
		case 'f':
			if (i + N_FILES >= argc) {
				err = mmessage(ERROR_MSG, INVALID_CMD_ARGUMENT,
						   "Too few arguments to --fsa_files "
						   "command-line option.\n");
				goto CMDLINE_ERROR;
			}
			mmessage(INFO_MSG, NO_ERROR, "Fasta files:");
			for (j = 0; j < N_FILES; ++j) {
				opt->fsa_files[j] = argv[++i];
				fprintf(stderr, " %s",
					opt->fsa_files[j]);
			}
			fprintf(stderr, "\n");
			break;
		case 'g':
			if (!strncmp(&argv[i][j], "geno", 4)) {
				mmessage(INFO_MSG, NO_ERROR, "Sam file of aligned targets:");
				opt->sam_file = argv[++i];
				fprintf(stderr, " %s",
					opt->sam_file);
				fprintf(stderr, "\n");
			} else if (!strncmp(&argv[i][j], "gl", 2)) {
				opt->vcf_opt->output_gl
				= !opt->vcf_opt->output_gl;
				mmessage(INFO_MSG, NO_ERROR,
					 "Outputing GL to vcf files: %s\n",
					 opt->vcf_opt->output_gl ? "yes" : "no");
			}
			break;
		case 'h':
			fprint_usage(stderr, argv[0], opt);
			exit(EXIT_SUCCESS);
		case 'i':
			opt->drop_indel = read_uint(argc, argv, ++i, opt);
			mmessage(INFO_MSG, NO_ERROR, "Dropping reads with indel"
				 " longer than %u in either alignment.\n",
				 opt->drop_indel);
			break;
		case 'j':
			mmessage(INFO_MSG, NO_ERROR, "Name for extracted reference file:");
			opt->extracted_rf = argv[++i];
			fprintf(stderr, " %s", opt->extracted_rf);
			fprintf(stderr, "\n");
			break;
		case 'l':
			opt->min_log_likelihood = read_cmdline_double(argc,
								argv, ++i, opt);
			mmessage(INFO_MSG, NO_ERROR, "Dropping reads with log "
						"likelihood less than %f.\n",
							opt->min_log_likelihood);
			break;

		case 'm':
			if (!strncmp(&argv[i][j], "min_s", 5)
				|| !strncmp(&argv[i][j], "min-s", 5)) {
				opt->min_expected_coverage
					= read_cmdline_double(argc, argv, ++i,
									opt);
			} else if (!strncmp(&argv[i][j], "mi", 2)) {
				opt->min_length = read_uint(argc, argv, ++i, opt);
				mmessage(INFO_MSG, NO_ERROR, "Minimum read "
					 "length: %u\n", opt->min_length);
			} else if (!strncmp(&argv[i][j], "ma", 2)) {
				opt->max_length = read_uint(argc, argv, ++i, opt);
				mmessage(INFO_MSG, NO_ERROR, "Maximum read "
					 "length: %u\n", opt->max_length);
			}
			break;
		case 'n':
			if (i + 1 >= argc) {
				err = mmessage(ERROR_MSG, INVALID_CMD_ARGUMENT,
						"Option --name needs an argument.\n");
				goto CMDLINE_ERROR;
			}
			opt->sample_name = argv[++i];
			break;
		case 'r':
			if (i + N_FILES >= argc) {
				err = mmessage(ERROR_MSG, INVALID_CMD_ARGUMENT,
						   "Too few arguments to --ref_names "
						   "command-line option.\n");
				goto CMDLINE_ERROR;
			}
			mmessage(INFO_MSG, NO_ERROR, "Reference names:");
			for (j = 0; j < N_FILES; ++j) {
				opt->ref_names[j] = argv[++i];
				fprintf(stderr, " %s",
					opt->ref_names[j]);
			}
			fprintf(stderr, "\n");
			break;
		case 's':
			if (!strncmp(&argv[i][j], "se", 2)) {
				opt->drop_secondary = !opt->drop_secondary;
				mmessage(INFO_MSG, NO_ERROR, "%s secondary "
					 "alignments\n", opt->drop_secondary
					 ? "Dropping" : "Keeping");
			} else if (!strncmp(&argv[i][j], "so", 2)) {
				opt->drop_soft_clipped
					= read_uint(argc, argv, ++i, opt);
				mmessage(INFO_MSG, NO_ERROR, "Dropping reads "
					 "with soft clip >= %u in either "
					 "alignment.\n", opt->drop_soft_clipped);
			} else if (!strncmp(&argv[i][j], "samp", 4)) {
				opt->n_sample = read_uint(argc, argv, ++i, opt);
				mmessage(INFO_MSG, NO_ERROR, "%u Monte Carlo "
					 "samples.\n", opt->n_sample);
			} else if (!strncmp(&argv[i][j], "sam", 3)) {
				if (i + N_FILES >= argc) {
					err = mmessage(ERROR_MSG,
							   INVALID_CMD_ARGUMENT, "Too few "
							   "arguments to --sam_files "
							   "command-line option.\n");
					goto CMDLINE_ERROR;
				}
				mmessage(INFO_MSG, NO_ERROR, "Sam files:");
				for (j = 0; j < N_FILES; ++j) {
					opt->sbam_files[j] = argv[++i];
					fprintf(stderr, " %s",
						opt->sbam_files[j]);
				}
				fprintf(stderr, "\n");
			}
			break;

		case 'u':
			opt->drop_unmapped = !opt->drop_unmapped;
			mmessage(INFO_MSG, NO_ERROR, "%s unmapped reads\n",
				 opt->drop_unmapped ? "Dropping" : "Keeping");
			break;
				
		case 'v':
			if (i + N_FILES >= argc)
				goto CMDLINE_ERROR;
			for (j = 0; j < N_FILES; ++j)
				opt->vcf_files[j] = argv[++i];
			if (!opt->vcf_opt && (err = make_default_vcf_options(&opt->vcf_opt)))
				goto CMDLINE_ERROR;
			break;
		case 'o':
			opt->output_file = argv[++i];
			mmessage(INFO_MSG, NO_ERROR,
				 "Final output file: %s\n",
				 opt->output_file);
			break;
		case 'w':
			opt->write_fastq_and_quit = 1;
			mmessage(INFO_MSG, NO_ERROR, "Will write selected "
				"reads in fastq file and quit.\n");
			break;

		case 'p':
			if (!strncmp(&argv[i][j], "po", 2)) {	/* hidden legacy */
				opt->equal_homolog_coverage_test = 1;
				mmessage(INFO_MSG, NO_ERROR, "Will run "
						"equal coverage test.\n");
				opt->posthoc_coverage_test = 1;
				mmessage(INFO_MSG, NO_ERROR, "Performing "
					"post-hoc confidence interval of equal "
					"homologous chromosome coverage\n");
				break;
			}
			if (i + 1 >= argc)
				goto CMDLINE_ERROR;
			opt->proptest_screen = (unsigned char)
						read_uint(argc, argv, ++i, opt);
			mmessage(INFO_MSG, NO_ERROR, "Dropping reads assigned "
				"to all but %u most-abundant haplotypes\n",
				 4 - opt->proptest_screen);
			break;

		case 'c':
			opt->weight_penalty = read_cmdline_double(argc, argv,
								++i, opt);
			mmessage(INFO_MSG, NO_ERROR, "Add tunning parameter for penalty: %f\n",
							opt->weight_penalty);
			break;

		default:
			err = INVALID_CMD_OPTION;
			goto CMDLINE_ERROR;
		}
	}
	opt->vcf_opt->equal_homolog_coverage_test = opt->equal_homolog_coverage_test;
	opt->vcf_opt->posthoc_coverage_test = opt->posthoc_coverage_test;
	return err;

CMDLINE_ERROR:
	if (err == NO_ERROR) {
		err = INVALID_CMD_ARGUMENT;
		i--;
	}
	usage_error(argv, i, (void *)opt);
	return err;
} /* parse_options */

void fprint_usage(FILE *fp, const char *cmdname, void *obj) {
	options *opt = (options *) obj;
	size_t start = strlen(cmdname) - 1;

	while (cmdname[start] != '/' && start) start--;
	if (cmdname[start] == '/') start++;

	for (size_t i = start; i < strlen(cmdname); ++i)
		fputc(toupper(cmdname[i]), fp);
	fprintf(fp, "(%d)\n", 1);
	fprintf(fp, "\nNAME\n\t%s - genotype tetraploids\n", &cmdname[start]);
	fprintf(fp, "\nSYNOPSIS\n\t%s --sam_files <fsam1> <fsam2> --fsa_files "
		"<fsa1> <fsa2> --ref_names <sref1> <sref2> --g <refsam> --j <reffsa>\n\t\t"
		"[--vcf_files --min-subgenomic-coverage <dbl>]\n\t\t"
		"[--min <int> --max <int> --expected-errors <dbl> --indel <int> --loglike <dbl>"
		" --secondary --soft-clipped <int> --coverage <dbl>]\n\t\t"
//		"[--p <int> --amp <exe> [--ampliclust-f <ffastq> --ampliclust-o <str> --ampliclust-l <dbl>]]\n\t\t"
		"[ --o <fout> --error_file|--error_data <ferr>] ...]\n",
		&cmdname[start]);
	fprintf(fp, "\nDESCRIPTION\n\t%s genotypes tetraploids' targeted genome regions using"
		" screened reads in <fsam1> and <fsam2> aligned to the whole genome references"
		" from fasta files <fsa1> <fsa2>\n",
		&cmdname[start]);
	fprintf(fp, "\nNOTICE\n\t%s requires the reference names (appeared in <refsam>) to contain"
		" start (0 based) and end position (1 based) of the targeted genome regions "
		"relative to the whole genome. Using ':' to seperate original genome name "
		"(appeared in <fsam1> and <fsam2>) and region index, '-' to seperate start "
		"and end position (e.g. chr1:0-11, which starts at 1st position in 'chr1' and "
		"end at 11th position, length is 11 bases)\n",
		&cmdname[start]);
	fprintf(fp, "\nOPTIONS\n");
	fprintf(fp, "\t+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	fprintf(fp, "\tInput: all required\n");
	fprintf(fp, "\t+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	fprintf(fp, "\t--fsa_files <fsa1> <fsa2>\n\t\tSpecify fasta files "
		"containing subgenomic reference sequences [use samtools faidx to index the files]"
		" (Default: none).\n");
	fprintf(fp, "\t--sam_files <fsam1> <fsam2>\n\t\tSpecify sam files "
				"containing alignments (Default: none)\n");
//	fprintf(fp, "\t--bam_files <fbam1> <fbam2>\n\t\tSpecify bam files "
//				"containing alignments (Default: none)\n");
	fprintf(fp, "\t--ref_names <sref1> <sref2>\n\t\tSpecify names of "
		"subgenomic references for target region; must exist in sam "
						" files (Default: none)\n");
	fprintf(fp, "\t--geno <refsam>\n\t\tSpecify name of sam file of aligning "
		"<sref2> to <sref1> (Default: none)\n");
	fprintf(fp, "\t--j <reffsa>\n\t\tSpecify prefix of targeted fsa files to be extracted "
		" by samtools [Set samtools in system PATH] (Default: extracted)\n");
	fprintf(fp, "\t+++++++++++++++++\n");
	fprintf(fp, "\tOutput:  all optional except for --vcf_files\n");
	fprintf(fp, "\t+++++++++++++++++\n");
//	fprintf(fp, "\t--censor <int>\n\t\tCurrent code censors quality scores "
//		"at maximum 41 [This option is not used!]\n");
	fprintf(fp, "\t--display_alignment\n\t\tDisplay alignments in stderr "
		"output (Default: %s).\n", opt->display_alignment ? "yes" : "no");
	fprintf(fp, "\t--vcf_files FILE1 FILE2\n"
		"\t\tGenotyping output in one vcf file per subgenome (Default: none).\n");
	fprintf(fp, "\t--gl\n"
		"\t\tToggle GL output to vcf files (Default: %s).\n", opt->vcf_opt->output_gl ? "yes" : "no");
	fprintf(fp, "\t--name STRING\n"
		"\t\tName of accession/individual/genotype; used in vcf header (Default: none).\n");
	fprintf(fp, "\t--o <fout>\n\t\tSpecify file with "
		"the final output (Default: none).\n");
	fprintf(fp, "\t++++++++++++++++++++++++++++++\n");
	fprintf(fp, "\tError recalibration:  optional\n");
	fprintf(fp, "\t++++++++++++++++++++++++++++++\n");
	fprintf(fp, "\t--error_file <ferr>\n\t\tSpecify file with "
		"estimates of error rates (Default: none).\n");
	fprintf(fp, "\t--error_data <ferr>\n\t\tSpecify file to output "
		"observed errors (Default: none).\n");
//	fprintf(fp, "\t++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
//	fprintf(fp, "\tScreening paralogs and other contaminants:  optional\n");
//	fprintf(fp, "\t++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
//	fprintf(fp, "\t--ampliclust <sampliclust>\n\t\tDrop reads aligning to "
//		"low abundance haplotypes identified by amplicon clusterer "
//				"executable <sampliclust> (Default: no).\n");
//	fprintf(fp, "\t\tUses auxiliary files \"%s\" \"%s.fa\", and \"%s.out\" "
//			"(see --ampliclust_fastq)\n", opt->ac_fastq_file,
//					opt->ac_outfile, opt->ac_outfile);
//	fprintf(fp, "\t--ampliclust_fastq <sampliclust_f>\n\t\tName of fastq "
//		"file where reads written for amplicon clusterer (Default: %s)."
//						"\n", opt->ac_fastq_file);
//	fprintf(fp, "\t--write-fastq\n\t\tJust write fastq file (Default: "
//		"%s).\n", opt->write_fastq_and_quit ? "yes" : "no");
//	fprintf(fp, "\t\t\tSee --ampliclust_fastq to name the file.\n");
//	fprintf(fp, "\t--ampliclust_output <sampliclust_output>\n\t\tName of "
//			"amplicon clusterer output files (Default: %s.fa and "
//				"%s.out).\n", opt->ac_outfile, opt->ac_outfile);
//	fprintf(fp, "\t--ampliclust_low_bound <fampliclust_lb>\n\t\t"
//		"amplicon clusterer lower bound (-lb option to amplici)."
//					" (Default: %f)\n", opt->ac_low_bound);
	fprintf(fp, "\t+++++++++++++++++++++++++++++++++++++++++++\n");
	fprintf(fp, "\tScreening reads, coverage checks:  optional\n");
	fprintf(fp, "\t+++++++++++++++++++++++++++++++++++++++++++\n");
	fprintf(fp, "\t--eq <pint>\n\t\tEqual coverage test. (Default: %s)\n", opt->equal_homolog_coverage_test ? "yes" : "no");
	fprintf(fp, "\t--expected_errors <dbl>\n\t\tDiscard reads with more "
		"than <dbl> expected errors (Default: %f).\n", opt->max_eerr);
	fprintf(fp, "\t--indel <i>\n\t\tDrop reads with alignments containing "
		"more than <i> indels (Default: %u)\n", opt->drop_indel);
	fprintf(fp, "\t--loglik <l>\n\t\tDrop reads with log likelihood less "
		"than <l> (Default: %f)\n", opt->min_log_likelihood);
	fprintf(fp, "\t--biallelic FLOAT\n"
		"\t\tSkip site if third allele >100*FLOAT%% of minimum subgenomic coverage (Default: %.1f).\n", opt->biallelic_screen);
	fprintf(fp, "\t--min <dbl>\n\t\tDrop reads shorter than <dbl> (Default:"
		" %u)\n", opt->min_length);
	fprintf(fp, "\t--max <dbl>\n\t\tDrop reads longer than <dbl> (Default: "
		"%u)\n", opt->max_length);
	fprintf(fp, "\t--secondary\n\t\tDrop secondary alignments (Default: "
		"%s)\n", opt->drop_secondary ? "yes" : "no");
	fprintf(fp, "\t--coverage <c>\n\t\tTuning parameter for penalty (Default: %.1f).\n", opt->weight_penalty);
	/* hide legacy CI
	fprintf(fp, "\t--eq\n\t\tPost-hoc test of equal "
		"coverage of homologous chromosomes. (Default: %s)\n",
		opt->posthoc_coverage_test ? "no" : "yes");*/
	fprintf(fp, "\t--soft-clipped <s>\n\t\tDrop reads where either "
		"alignment is clipped by <s> or more nucleotides (Default: %du\n",
		opt->drop_soft_clipped);
	fprintf(fp, "\t--unmapped\n\t\tDrop reads unmapped in either "
		"alignment (Default: %s)\n", opt->drop_unmapped ? "yes" : "no");
	fprintf(fp, "\t+++++++++++++++++++++++++++++++\n");
	fprintf(fp, "\tEstimation/Inference:  optional\n");
	fprintf(fp, "\t+++++++++++++++++++++++++++++++\n");
//	fprintf(fp, "\t--sample <nsamp>\n\t\tNumber of Monte Carlo samples "
//		"(Default: %u)\n", opt->n_sample);
	fprintf(fp, "\t--min_subgenomic_coverage <c>\n\t\tAbort if subgenomic "
				"coverage drops below <c> (Default: %.1f).\n",
						opt->min_expected_coverage);
} /* fprint_usage */
