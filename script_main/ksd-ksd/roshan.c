/**
 * @file roshan.c
 * @author K. S. Dorman
 *
 * This file contains the code for genotyping allotetraploid
 * amplicon data.
 *
 * To compile:
 make roshan
 *
 * The number of subgenomes is hard-coded to two (allotetraploid) in places.
 *
 * Working on MERGING WGS and amplicon
 */

#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <unistd.h>
#include <zlib.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "roshan.h"
#include "sam.h"
#include "vcf.h"
#include "fastq.h"
#include "nuc.h"
#include "qual.h"
#include "uthash.h"
#include "io.h"
#include "cmdline.h"
#include "array.h"
#include "order.h"
#include "nmath.h"
#include "align.h"
#include "pick_reads.h"

/**
 * State of a read nucleotide, including pointer to containing alignment,
 * probability of error, position in alignment, true nucleotide, observed
 * read nucleotide, and quality score.
 * http://www.catb.org/esr/structure-packing/
 */
typedef struct nuc_state {
	UT_hash_handle hh;		/* 56 bytes w/ pointer alignment */
	double prob;			/* 8 bytes */
	unsigned int pos;		/* 4 */
	data_t true_nuc, read_nuc;	/* 1, 1 bytes */
	data_t quality;			/* 1 */
	/* total: 71 bytes round up to
	 * 72 with 1 byte slop
	 */
} nuc_state;

typedef struct mlogit_stuff {
	nuc_state *ns;
	int pos;
} mlogit_stuff;

unsigned int ns_key_len = 3 * sizeof(data_t) + sizeof(unsigned int);

int default_options(options *opt);
int parse_options(options *opt, int argc, const char **argv);

nuc_state *read_param_file(char const *param_file);
double mlogit_sub_prob(data_t s, data_t r, data_t q, void *vptr);
double mlogit_sub_lprob(data_t s, data_t r, data_t q, void *vptr);

double ll_align(sam_entry *se, unsigned int i, unsigned char *ref, mlogit_stuff *vptr, unsigned char *show);

double max_multinom_ll(unsigned int p, unsigned int *cnts);
double four_haplotype(double w, void *param);
double three_haplotype(double w, void *param);
int discard_haplotypes(options *opt, input *input);
unsigned int exclude_cluster_reads(merge_hash *mh, input *in, unsigned int id);
int genotype_by_clustering(options *opt, input *in, merge_hash *mh, sam *sds[N_SUBGENOMES], char_t const * const rseqs[N_SUBGENOMES], unsigned int const rlens[N_SUBGENOMES]);
uint64_t get_haplotype_id(sam_entry *se, size_t *haplotype, unsigned int n_segregating);
int get_reference_aligned_haplotypes(input *in, fastq_data *fd, char_t const *subgenomic_profiles[2*N_SUBGENOMES], char_t const * const rseqs[N_SUBGENOMES], unsigned int const rlens[N_SUBGENOMES]);
int load_chromosomes(options *opt, input *in);
int align_haplotypes(fastq_data **fd, options *opt);
int get_aligned_haplotypes(input *in, options *opt, fastq_data **fd, char_t const *haplotypes[2*N_SUBGENOMES]);
int test_equal_homolog_coverage(merge_hash *mh, double **ll, char_t ref_base[N_SUBGENOMES], char *covers, xy_t *obs_nuc, qual_t *obs_q, unsigned int *obs_rpos, int g_max[N_SUBGENOMES], xy_t nuc1, xy_t nuc2, int debug_level, double pvals[N_SUBGENOMES]);

int read_amplici_results(FILE *fp, input *in, options *opt);
int make_input(input **in);
void free_input(input *in);

int update_vcf(options *opt, int *fail[N_SUBGENOMES], size_t *hpos[N_SUBGENOMES]);

/* roshan's research */
int main(int argc, const char *argv[])
{
	int karin_version = 1;
	int debug_level = QUIET;//ABSOLUTE_SILENCE;//MINIMAL;//DEBUG_I;//
	int err = NO_ERROR;
		/* command line options */
	options opt;
		/* command line options related to reference */
	options_rf opt_rf;
		/* information about subsetted/targetted reference */
	ref_info *rf_info = NULL;
		/* reference sequences in FASTA format */
	fastq_data *fds[N_SUBGENOMES] = {NULL, NULL};
		/* indices of selected references */
	unsigned int fs_index[N_SUBGENOMES] = {0, 0};
		/* reference sequences */
	char_t const *ref_seqs[N_SUBGENOMES] = {NULL, NULL};
		/* reference lengths: by design are equal */
	unsigned int ref_lengths[N_SUBGENOMES];
		/* sam data */
	sam *sds[N_SUBGENOMES];
	fastq_options fop = {.read_encoding = IUPAC_ENCODING, .read_names = 1};
		/* sam hashed by reference name */
	sam_hash *by_name[N_SUBGENOMES] = {NULL, NULL};
		/* indices of selected references in sam files */
	size_t my_refs[N_SUBGENOMES] = {0, 0};	
	FILE *fp = NULL;
	gzFile gzfp;
	mlogit_stuff mls = {NULL, 0};
	double **pp = NULL;	/* posterior alignment probability */
	double **ll = NULL;	/* log likelihood of alignment */

	/* parse command line */
	default_options(&opt);
	if ((err = parse_options(&opt, argc, argv)))
		exit(err);

	if (!opt.amplicon) {
		default_rf_options(&opt_rf);
		opt_rf.sam_file = opt.ref_alignment;

		/* choose names for fasta files to store selected reference */
		if (!opt.subref_fsas[0]) {
			size_t sr_fsa_len = strlen(opt.subref_fsa_base)
				+ strlen(".fsa") + log10(N_SUBGENOMES + 1) + 1;
			for (int i = 0; i < N_SUBGENOMES; ++i) {
				opt.subref_fsas[i] = malloc(sr_fsa_len
						* sizeof *opt.subref_fsas[i]);
				sprintf(opt.subref_fsas[i], "%s%d.fsa",
							opt.subref_fsa_base, i);
			}
		}

		/* store information about targeted reference regions */
		make_targets_info(opt_rf, &rf_info);

	} else {	/* MERGING: remove? */

		/* read in sam files and extract reads aligned to selected reference */
		for (unsigned int j = 0; j < N_SUBGENOMES; ++j) {
	
			/* read in reference genomes */
			if ((err = read_fastq(opt.fsa_files[j], &fds[j], &fop)))
				exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "Reading '%s' "
					      "failed with error '%s' (%d).\n",
					      opt.fsa_files[j], fastq_error_message(err),
					      err));
	
			/* find the desired reference in each sam file */
			size_t index = 0;
			for (unsigned int i = 0; i < fds[j]->n_reads; index +=
				fds[j]->name_lengths[i++]) {
			     	unsigned int len = read_length(fds[j], i);
				if (!strncmp(&fds[j]->names[index], opt.ref_names[j],
					     strlen(opt.ref_names[j]))) {
					    ref_lengths[j] = len;
					break;
				}
				fs_index[j] += len;
			}
			ref_seqs[j] = &fds[j]->reads[fs_index[j]];
		}
	}

	/* read sam/bam files with read alignments */
	for (unsigned int j = 0; j < N_SUBGENOMES; ++j) {
		/* read sam/bam files */
		if (opt.use_bam) {
			gzfp = gzopen(opt.sbam_files[j], "r");
			if (!gzfp)
				exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR,
					      opt.sbam_files[j]));
			read_bam(gzfp, &sds[j]);
			gzclose(gzfp);
		} else {
			fp = fopen(opt.sbam_files[j], "r");
			if (!fp)
				exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR,
					      opt.sbam_files[j]));
			read_sam(fp, &sds[j]);
			fclose(fp);
		}
		if (opt.amplicon) {	/* MERGING: remove? */
			/* find selected references index in sam files */
			size_t rchar = 0;
			unsigned char found = 0;
	
			for (size_t i = 0; i < sds[j]->n_ref; ++i) {
				if (!strcmp(opt.ref_names[j], &sds[j]->ref_names[rchar])) {
					my_refs[j] = i;
					found = 1;
					break;
				}
				rchar += strlen(&sds[j]->ref_names[rchar]) + 1;
			}
			if (!found)
				exit(mmessage(ERROR_MSG, INVALID_USER_INPUT, "no "
					      "reference '%s' in fasta file '%s'",
					      opt.ref_names[j], opt.sbam_files[j]));
	
			/* hash sam file to reference */
			hash_sam(sds[j], &by_name[j], HASH_REFERENCE, my_refs[j],
				 opt.drop_unmapped, opt.drop_secondary,
				 opt.drop_soft_clipped, opt.drop_indel,
				 opt.min_length, opt.max_length, opt.max_eerr);
	
			mmessage(INFO_MSG, NO_ERROR, "Number of %u alignments: %zu\n",
				 j, sds[j]->n_per_ref[my_refs[j]]);
		}
	}

	/* MERGING: apply universally */
	if (!opt.amplicon) {

		/* pick reads aligned to targetted regions */
		pickreads(rf_info, &opt_rf, sds);	/* MERGING: sets which_ref == -1 until excluded: why can't we use exclude? */
		mmessage(INFO_MSG, NO_ERROR, "assign the reads to targeted regions done\n");

//		char strand;
//	        for (unsigned int j = 0; j < N_SUBGENOMES; ++j) {
//			/* find selected references index in sam files */
//			unsigned char found = 0;
//
//			for (size_t i = 0; i < sds[j]->n_se; ++i) {	/* MERGING: why did we change this loop? */
//				sam_entry *se = &sds[j]->se[i];
//	
//				/* skip unmapped */
//				if ((se->flag & (1 << 2)))
//					continue;
//				if (se->which_ref == -1)	/* MERGING: cannot we use excluded here? */
//					continue;
//	
//				if (!strcmp(opt.ref_names[j], se->ref_name)) {
//					se->name_s = NULL;
//	
//					/* strand for hashing on strand and name */
//					if ((se->flag & 16) == 0) {
//						strand = "+";
//						// flip the strand if A is aligned to reverse complement of B
//						if (j == 1 && rf_info->info[se->which_ref].strand_B == 1)	/* TODO: assumes N_SUBGENOMES == 2 */
//							strand = "-";
//					} else {
//						strand = "-";
//						if (j == 1 && rf_info->info[se->which_ref].strand_B == 1)	/* TODO: assumes N_SUBGENOMES == 2 */
//							strand = "+";
//					}
//	
//					size_t length = strlen(se->name) + 2;
//					se->name_s = malloc(length);
//					sprintf(se->name_s, "%s%s", se->name, strand);	/* TODO: why encode strand information in a big character string w/ repetitive information */
//	
//					my_refs[j] = se->which_ref; // this my_refs index should be adjusted	/* TODO: repeated again and again...why? */
//					found = 1;
//				}
//			}
//			printf("\n");
//			if (!found)
//				exit(mmessage(ERROR_MSG, INVALID_USER_INPUT, "no "
//					      "reference '%s' in fasta file '%s'",
//					      opt.ref_names[j], opt.sbam_files[j]));
//	
//			/* hash sam file to reference (use n_se since some references are repeated in the targted sam file) */
//			hash_sam(sds[j], &by_name[j], HASH_REFERENCE, my_refs[j], rf_info->ref_sam->n_se,
//				 opt.drop_unmapped, opt.drop_secondary,
//				 opt.drop_soft_clipped, opt.drop_indel,
//				 opt.min_length, opt.max_length, opt.max_eerr);
//	
//			mmessage(INFO_MSG, NO_ERROR, "Number of %u alignments: %zu\n",
//				 j, sds[j]->n_per_ref[my_refs[j]]);
//	
//		}
	}


	/* merge reads for given reference */
	merge_hash *mh = NULL;

	size_t n_read = hash_merge(&mh, N_SUBGENOMES, sds, my_refs);
	debug_msg(debug_level > ABSOLUTE_SILENCE, debug_level,
			  "Number of aligned reads: %zu\n", n_read);
	
	if (!opt.amplicon && opt_rf.fastq_file)	/* MERGING: why? */
		output_selected_reads(opt_rf.fastq_file, sds, mh);

	/* prepare for AmpliCI-based clustering, filtering, and genotyping */
	size_t n_included_reads = 0;
	input *in = NULL;
	if (opt.amplici_command)
		make_input(&in);

	/* count number of alignments and obtain their extent on reference */
	size_t n_align = 0;
	size_t start_pos[N_SUBGENOMES];	/* 0-based starting position */
	size_t end_pos[N_SUBGENOMES];	/* 0-based ending position */
	for (unsigned int i = 0; i < N_SUBGENOMES; ++i) {
		start_pos[i] = SIZE_MAX;
		end_pos[i] = 0;
	}

	/* exclude reads not aligned to both subgenomes or aligned multiply
	 * to one subgenome; identify extent of alignments on reference genome
	 */
	n_read = 0;
	for (merge_hash *me = mh; me != NULL; me = me->hh.next) {
		sam_entry *se;
		size_t rf_pos;

		if (me->nfiles != N_SUBGENOMES) {
			me->exclude = 1;
			mmessage(INFO_MSG, NO_ERROR, "Read %s does not "
				 "align to all genomes (skipping).\n",
				 me->indices[0]	/* TODO: hack */
				 ? sds[0]->se[me->indices[0][0]].name
				 : sds[1]->se[me->indices[1][0]].name);
			continue;
		}

		/* force one alignment per sub-genome */
		for (unsigned int j = 0; j < N_SUBGENOMES; ++j) {

			if (me->count[j] > 1)
				exit(mmessage(ERROR_MSG, INTERNAL_ERROR,
					      "Read %u aligns twice in genome %u.\n",
					      j, sds[j]->se[me->indices[j][0]].name));

			n_align += me->count[j];
			se = &sds[j]->se[me->indices[j][0]];
			if (start_pos[j] > se->pos - 1)
				start_pos[j] = se->pos - 1;
			rf_pos = se->pos - 1;
			for (unsigned int i = 0; i < se->cig->n_ashes; ++i) {
				if (se->cig->ashes[i].type == CIGAR_DELETION
				    || se->cig->ashes[i].type == CIGAR_MATCH
				    || se->cig->ashes[i].type == CIGAR_MMATCH
				    || se->cig->ashes[i].type == CIGAR_MISMATCH
				    || se->cig->ashes[i].type == CIGAR_SKIP)
					rf_pos += se->cig->ashes[i].len;
			}
			if (end_pos[j] < rf_pos)
				end_pos[j] = rf_pos;

		}

		++n_read;
	}

	if (n_read == 0)
		exit(mmessage(ERROR_MSG, INTERNAL_ERROR,
			"No read aligns to selected genomic regions.\n"));


	for (unsigned int i = 0; i < N_SUBGENOMES; ++i) {
// MERGING:		--end_pos[i];	/* 0-based index */
		mmessage(INFO_MSG, NO_ERROR, "File %u extent: %u - %u\n",
			 i, start_pos[i], end_pos[i]);
	}

	/* [TODO, BUG] Fails to account for indels between subgenome references;
	 * [TODO, BUG] current solution is to pre-align the references,
	 * [TODO, BUG] inserting N for gaps.
	 *
	 * match soft clips so two subgenomic alignments are trimmed to cover
	 * the same region; additional reads may be dropped if there is no
	 * alignment overlap
	 */
//	if(!rf_info->info[my_refs[1]].strand_B)
//		match_soft_clipping(mh, N_SUBGENOMES, sds, start_pos, end_pos,
//					rf_info->info[my_refs[1]].strand_B);
	match_soft_clipping(mh, N_SUBGENOMES, sds, start_pos);

	/* optionally read error rates callibrated by R's mlogit */ 
	if (opt.param_file) {
		debug_msg(debug_level > ABSOLUTE_SILENCE, debug_level,
			  "Reading parameter file '%s'.\n", opt.param_file);
		mls.ns = read_param_file(opt.param_file);
		sub_prob = mlogit_sub_prob;
		sub_lprob = mlogit_sub_lprob;
	}

	/* compute posterior probability of alignment to each subgenome */
	pp = malloc(N_SUBGENOMES * sizeof *pp);
	ll = malloc(N_SUBGENOMES * sizeof *ll);

	for (unsigned int j = 0; j < N_SUBGENOMES; ++j) {
		debug_msg(debug_level > ABSOLUTE_SILENCE, debug_level, "Genome "
			  "%u extent: %zu - %zu\n", j, start_pos[j], end_pos[j]);

		if (j && end_pos[j] - start_pos[j]
		    != end_pos[j-1] - start_pos[j-1])
			mmessage(WARNING_MSG, NO_ERROR, "WARNING:  Genome %u "
				 "and %u alignment regions differ in length.\n",
				 j, j - 1);
		pp[j] = malloc(n_read * sizeof **pp);	/* overestimated */
		ll[j] = malloc(n_read * sizeof **ll);	/* overestimated */
	}

	/* [KSD, TODO, IDEA] Do not normalize this here.  Instead leave as log likelihood that then gets adjusted by calculations at the focal locus. */

	/* open fastq file for output of so-far selected reads */
	if (opt.amplici_command || opt.write_fastq_and_quit) {
		fp = fopen(opt.ac_fastq_file, "w");
		if (opt.amplici_command)
			in->n_hash_excluded = 0;
	}

	/* optional screen: maximum log likelihood of subgenomic alignments */
	double *mll = NULL;
	if (opt.min_log_likelihood > 0)
		MAKE_1ARRAY(mll, n_read);

	/* compute posterior probability read from each subgenome;
	 * screen reads on minimum log likelihood or post. prob.
	 */
	n_read = 0;
	size_t coverA = 0;	/* number of reads assigned to subgenome A */
	size_t n_ll_exclude = 0, n_pp_exclude = 0;
	double A_expected_coverage = 0;
	unsigned char show = opt.display_alignment;

	for (merge_hash *me = mh; me != NULL; me = me->hh.next) {

		if (me->exclude)
			continue;

		double sum = 0, max = -INFINITY;

		for (unsigned int j = 0; j < N_SUBGENOMES; ++j) {

			ll[j][n_included_reads] = ll_align(
				 &sds[j]->se[me->indices[j][0]], n_read,
				 &fds[j]->reads[fs_index[j]], &mls, &show);
			if (max < ll[j][n_included_reads])
				max = ll[j][n_included_reads];
		}

		if (opt.min_log_likelihood > 0) {
			mll[n_read] = max;
		} else if (isfinite(opt.min_log_likelihood)
					&& max < opt.min_log_likelihood) {
			me->exclude = 1;
			++n_ll_exclude;
			show = opt.display_alignment;
			++n_read;
			continue;
		}

		debug_msg(show || debug_level > QUIET, debug_level,
							"Log likelihoods:");
		for (unsigned int j = 0; j < N_SUBGENOMES; ++j) {
			sum += exp(ll[j][n_included_reads] - max);
			debug_msg_cont(show || debug_level > QUIET, debug_level,
					       " %f", ll[j][n_included_reads]);
		}
		debug_msg_cont(show || debug_level > QUIET, debug_level, "\n");
		debug_msg(show || debug_level > QUIET, debug_level,
							"Probabilities:");
		double max_p = -INFINITY;
		unsigned int max_p_index = 0;
		for (unsigned int j = 0; j < N_SUBGENOMES; ++j) {
			pp[j][n_included_reads] = exp(ll[j][n_included_reads]
								- max) / sum;
			debug_msg_cont(show || debug_level > QUIET, debug_level,
				       " %f", pp[j][n_included_reads]);
			if (max_p < pp[j][n_included_reads]) {
				max_p = pp[j][n_included_reads];
				max_p_index = j;
			}
			if (opt.error_file)
				output_error_data(opt.error_file,
					  &sds[j]->se[me->indices[j][0]],
						  &fds[j]->reads[fs_index[j]],
					  log(pp[j][n_included_reads]), 0);
		}
		if (max_p_index > 1)
			debug_msg(debug_level > QUIET, debug_level, "***!!!***");
		debug_msg_cont(show || debug_level > QUIET, debug_level, "\n");

		debug_msg(debug_level > QUIET, debug_level, "Log Probabilities:");
		for (unsigned int j = 0; j < N_SUBGENOMES; ++j)
			debug_msg_cont(debug_level > QUIET, debug_level, " %f",
				       log(pp[j][n_included_reads]));

		if (max_p < opt.min_posterior_alignment_prob) {
			me->exclude = 1;
			++n_pp_exclude;
			show = opt.display_alignment;
			++n_read;
			debug_msg_cont(debug_level > QUIET, debug_level,
				" posterior assignment probability < %f ... "
				"REMOVING!!\n", opt.min_posterior_alignment_prob);
			continue;
		}
		debug_msg_cont(debug_level > QUIET, debug_level, "\n");

		/* write fastq of selected reads for amplici */
		if (opt.amplici_command || opt.write_fastq_and_quit) {
			sam_entry *se = &sds[0]->se[me->indices[0][0]];
			fprintf(fp, "@%s ppA=%.6e\n", se->name,
						pp[0][n_included_reads]);
			fwrite_nuc_segment(fp, se->read, XY_ENCODING, 0,
							se->read->len);
			fprintf(fp, "\n+\n");
			fwrite_qual_sequence(fp, se->qual);
			fprintf(fp, "\n");

		}

		if (pp[0][n_included_reads] > 0.5)
			++coverA;
		A_expected_coverage += pp[0][n_included_reads];
		++n_included_reads;
		++n_read;
		show = opt.display_alignment;

	}

	/* display to stderr: sorted maximum subgenomic log likelihoods */
	if (opt.min_log_likelihood > 0) {
		qsort(mll, n_read, sizeof(double), double_compare);
		debug_msg(1, debug_level, "mll: ");
		fprint_doubles(stderr, mll, n_read, 6, 1);
		fprintf(stderr, "%f %f %f %f %f %f %f\n", mll[(int)(n_read*0.025)], mll[(int)(n_read*0.05)], mll[(int)(n_read*0.10)], mll[(int)(n_read*0.5)], mll[(int)(n_read*0.9)], mll[(int)(n_read*.95)], mll[(int)(n_read*.975)]);
	}

	if (isfinite(opt.min_log_likelihood) && n_ll_exclude)
		mmessage(INFO_MSG, NO_ERROR, "%zu reads excluded by failing to "
			"reach minimum log likelihood %f.\n", n_ll_exclude,
							opt.min_log_likelihood);
	if (opt.min_posterior_alignment_prob > 0 && n_pp_exclude)
		mmessage(INFO_MSG, NO_ERROR, "%zu reads excluded by failing to "
			"reach minimum posterior alignment probability %f.\n",
				n_pp_exclude, opt.min_posterior_alignment_prob);

	/* discard posterior probabilities of discarded reads */
	for (unsigned int j = 0; j < N_SUBGENOMES; ++j) {
		double *new_pp = realloc(pp[j], n_included_reads * sizeof **pp);
		double *new_ll = realloc(ll[j], n_included_reads * sizeof **ll);
		if (!new_pp || !new_ll)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"(reallocate) posterior probability (%zu).\n",
								n_included_reads);
		pp[j] = new_pp;
		ll[j] = new_ll;
	}

	/* command-line option --error_data outputs error data for mlogit fit */
	if (opt.error_file)
		fclose(opt.error_file);

	debug_msg(debug_level > QUIET, debug_level, "Read coverage probability:"
				" %f %f\n", (double) coverA / n_included_reads,
		(double) (n_included_reads - coverA) / n_included_reads);
	mmessage(INFO_MSG, NO_ERROR, "Expected coverage: %f %f (%zu)\n",
		A_expected_coverage, n_included_reads - A_expected_coverage,
								n_included_reads);

	double B_expected_coverage = n_included_reads - A_expected_coverage;
	double min_expected_coverage = MIN(A_expected_coverage,
							B_expected_coverage);

	/* write reads as fastq and end if that was all we wanted;
	 * otherwise, use AmpliCI to detect paralogs and screen reads
	 * genotyping-by-clustering is in here
	 */
	if (opt.amplici_command || opt.write_fastq_and_quit) {

		fclose(fp);

		if (opt.write_fastq_and_quit) {
			mmessage(INFO_MSG, NO_ERROR, "Finished writing fastq "
					"file '%s'.\n", opt.ac_fastq_file);
			goto EXIT_NOW;
		}

		in->n_observation = n_included_reads;

		/* perform haplotype abundance screen to id likely paralogs:
		 * we expect 2-4 distinct haplotypes (not 1 because amplicons
		 * are selected to have homoeologous SNPs), but while we do not
		 * assume equal subgenomic coverage, we do assume equal allelic
		 * coverage
		 */
		if (opt.proptest_screen < 0)
			if ((err = discard_haplotypes(&opt, in)))
				return err;

		/* or take them from the command-line (deprecated legacy) */
		if (!in->K) {
			in->proptest = opt.proptest_screen;

			/* call amplici: it will give 3 files: *.out, *.fa */
			unsigned int cmd_len = strlen(opt.amplici_command)
					+ strlen(opt.ac_fastq_file)
					+ strlen(" -f  -lb  -ll  -o ")
					+ strlen(opt.ac_outfile) + 8
					+ (int)(log10(opt.ac_low_bound) + 2)
					+ (int)(log10(opt.ac_ll_threshold) + 1);
			char *command = malloc(cmd_len * sizeof *command);
			sprintf(command, "%s -f %s -lb %.1f -ll %.0f -o %s",
				opt.amplici_command, opt.ac_fastq_file,
				opt.ac_low_bound, opt.ac_ll_threshold,
							opt.ac_outfile);
	
			mmessage(INFO_MSG, NO_ERROR, "Running amplici to screen"
						" paralogs: '%s'\n", command);
			system(command);
			free(command);
	
			cmd_len = strlen(opt.ac_outfile) + strlen(".out") + 1;
			command = malloc(cmd_len * sizeof *command);
			sprintf(command, "%s.out", opt.ac_outfile);
			fp = fopen(command, "r");
			if (!fp) {
				mmessage(ERROR_MSG, FILE_OPEN_ERROR, command);
				free(command);
				goto EXIT_NOW;
			}
			free(command);
	
			if ((err = read_amplici_results(fp, in, &opt)))
				return err;
		}

		/* Now we will further drop haplotypes if they do not have
		 * expected subgenomic coverage; so far we have only checked
		 * the assumption of equal allelic coverage without verifying
		 * subgenomic assignments of reads from each haplotype.  We
		 * only get subgenomic assignments by using reference sequences.
		 *
		 * Let p1, p2, p3, and p4 be the relative abundance of top
		 * four most abundant haplotypes.  Let p1A, p2A, p3A, and p4A
		 * be the proportion of reads from each haplotype cluster
		 * assigned to A subgenome.  We expect:
		 * H04a: p1=p2=p3=p4 (equal subgenomic & allelic coverage)
		 *	p1A, p2A high; p3A, p4A low or vice versa or
		 *      p1A, p3A high; p2A, p4A low or vice versa and so on...
		 * H04: p1=p2, p3=p4 (equal allelic coverage)
		 *	p1A, p2A high; p3A, p4A low or vice versa
		 * H03a: p1, p2=p3 (equal allelic coverage, true SNP)
		 *	p1A high, p2A, p3A low or vice versa
		 * H03b: p1=p2, p3 (equal allelic coverage, true SNP)
		 *	p1A, p2A high, p3A low or vice versa
		 * H02: p1, p2 (homoeologous SNP)
		 *	p1A high, p2A low or vice versa
		 */
		/* [TODO,KSD] move to function */

		/* cluster indices of top sorted haplotypes */
		for (unsigned int k = 0; k < MIN(4, in->K); ++k)
			in->included_clusters[k] = in->sort_id[k];

		/* exclude all remaining clusters as presumed paralogs */
		for (unsigned int k = 4 - in->proptest; k < in->K; ++k)
			in->excluded_clusters[in->sort_id[k]] = 1;

		mmessage(INFO_MSG, NO_ERROR, "Cluster sizes: ");
		fprint_uints(stderr, in->cluster_sizes, in->K, 2, 1);
		mmessage(INFO_MSG, NO_ERROR, "Excluded clusters: ");
		fprint_uints(stderr, in->excluded_clusters, in->K, 1, 1);
		mmessage(INFO_MSG, NO_ERROR, "Keeping clusters: ");
		fprint_uints(stderr, in->included_clusters, MIN(in->K, 
						4 - in->proptest), 1, 1);
		mmessage(INFO_MSG, NO_ERROR, "Excluded reads:");

		/* initialize p1A, p2A, p3A, and p4A */
		unsigned int n_excluded = 0;	/* number to be excluded */
		in->coverage_A[0] = in->coverage_A[1] = in->coverage_A[2]
			= in->coverage_A[3] = 0;/* proportion A coverage of 
						 * top 4 clusters
						 */

		/* exclude reads of likely paralog; sum expected A coverage
		 * on retained haplotypes
		 */
		n_read = 0;
		for (merge_hash *me = mh; me != NULL; me = me->hh.next) {
			if (me->exclude)
				continue;
			if (in->assignment[n_read] == NA_ASSIGNMENT ||
				in->excluded_clusters[in->assignment[n_read]]) {
				fprintf(stderr, " %zu", n_read);
				++n_excluded;
			} else {
				/* looks slightly dangerous as no guarantee of 
				 * 4 haplotypes, but it short-circuits
				 */

				if (in->included_clusters[0]
					== in->assignment[n_read])
					in->coverage_A[0] += pp[0][n_read];
				else if (in->included_clusters[1]
					== in->assignment[n_read])
					in->coverage_A[1] += pp[0][n_read];
				else if (in->included_clusters[2]
					== in->assignment[n_read])
					in->coverage_A[2] += pp[0][n_read];
				else
					in->coverage_A[3] += pp[0][n_read];
			}
			++n_read;
		}

		fprintf(stderr, " (%u total from %u haplotypes)\n",
							n_excluded, MIN(in->K,
						in->K - 4 + in->proptest));

		for (unsigned int k = 0; k < MIN(4, in->K); ++k)
			in->proportion_A[k] = in->coverage_A[k]
				/ in->cluster_sizes[in->included_clusters[k]];

		mmessage(INFO_MSG, NO_ERROR, "Fraction expected A in top "
			"four largest clusters: %.3f (%u) %.3f (%u) %.3f (%u) "
			"%.3f (%u)\n",
			in->proportion_A[0], in->included_clusters[0],
			in->proportion_A[1], in->included_clusters[1],
			in->proportion_A[2], in->included_clusters[2],
			in->proportion_A[3], in->included_clusters[3]);

		/**
		 * check for problems with subgenomic coverage 
		 * At the end, input object will be updated with
		 * - input::reject_state: accepted hypothesis
		 * - input::excluded_clusters: cluster_id -> 0|1
		 * - input::included_clusters: cluster_ids of kept clusters
		 * - n_excluded: total number of excluded reads
		 * but the offending reads will NOT YET be excluded.
		 */

		/* Accept H04 (4 haplotypes), reject H04a (equal subgenomic coverage)
		 * Truth is: p1 = p2; p3 = p4 and coverage restrictions
		 *	|p12A - 0.5| >= 0.5 - z	: impure subgenome 1
		 *	|p34A - 0.5| >= 0.5 - z	: impure subgenome 2
		 *	sign(p12A - 0.5)	: map to same subgenome
		 *		= sign(p34A - 0.5)
		 * where z is \par options.max_misalignment_proportion \in [0, 0.5)
		 */
		double midz = 0.5 - opt.max_misalignment_proportion;
		if (in->reject_state == ACCEPT_FOUR &&
			in->pval_equal_subgenomic_coverage < opt.proptest_alpha) {
			
			double A12mid = -0.5 + (in->coverage_A[0] + in->coverage_A[1]) / 
				(in->cluster_sizes[in->included_clusters[0]]
				+ in->cluster_sizes[in->included_clusters[1]]);
			double A34mid = -0.5 + (in->coverage_A[2] + in->coverage_A[3]) / 
				(in->cluster_sizes[in->included_clusters[2]]
				+ in->cluster_sizes[in->included_clusters[3]]);

			/* violates coverage restrictions */
			if (fabs(A12mid) < midz				/* ~50% A coverage on subgenome 1 */
				|| fabs(A34mid) < midz			/* ~50% B coverage on subgenome 2 */
				|| SIGN(A12mid) == SIGN(A34mid)) {	/* >50% A coverage on both subgenomes */

				mmessage(INFO_MSG, NO_ERROR, "H04: impure "
					"subgenomic coverage (Proportion A: "
					"%.2f %.2f).  Discarding fourth most "
					"abundant haplotype (%u) and its "
					"reads:", A12mid + 0.5, A34mid + 0.5,
							in->sort_id[3]);

				/* tentatively accept H03b, will be tested */
				in->reject_state = ACCEPT_THREE_B;
				++in->proptest;
				in->excluded_clusters[in->sort_id[3]] = 1;
				n_excluded += exclude_cluster_reads(mh, in,
							in->sort_id[3]);
				fprintf(stderr, " (%u total from %u "
						"haplotypes)\n", n_excluded,
						in->K - 4 + in->proptest);
			} else {
				mmessage(INFO_MSG, NO_ERROR, "Accepting H04: "
								"p1=p2,p3=p4.\n");
			}

		/* Accept H04 (4 haplotypes), accept H04a (equal subgenomic coverage)
		 * Truth is: p1 = p2 = p3 = p4, but we don't know subgenomic identity.
		 */
		} else if (in->reject_state == ACCEPT_FOUR) {
			/* split 1: 1100 */
			double A12mid = -0.5 + (in->coverage_A[0] + in->coverage_A[1]) / 
				(in->cluster_sizes[in->included_clusters[0]]
				+ in->cluster_sizes[in->included_clusters[1]]);
			double A34mid = -0.5 + (in->coverage_A[2] + in->coverage_A[3]) / 
				(in->cluster_sizes[in->included_clusters[2]]
				+ in->cluster_sizes[in->included_clusters[3]]);
			/* split 2: 1010 */
			double A13mid = -0.5 + (in->coverage_A[0] + in->coverage_A[2]) / 
				(in->cluster_sizes[in->included_clusters[0]]
				+ in->cluster_sizes[in->included_clusters[2]]);
			double A24mid = -0.5 + (in->coverage_A[1] + in->coverage_A[3]) / 
				(in->cluster_sizes[in->included_clusters[1]]
				+ in->cluster_sizes[in->included_clusters[3]]);
			/* split 3: 1001 */
			double A14mid = -0.5 + (in->coverage_A[0] + in->coverage_A[3]) / 
				(in->cluster_sizes[in->included_clusters[0]]
				+ in->cluster_sizes[in->included_clusters[3]]);
			double A23mid = -0.5 + (in->coverage_A[1] + in->coverage_A[2]) / 
				(in->cluster_sizes[in->included_clusters[1]]
				+ in->cluster_sizes[in->included_clusters[2]]);

			/* neither possible split into subgenomes pure enough */
			if ((fabs(A12mid) < midz
					|| fabs(A34mid) < midz
					|| SIGN(A12mid) == SIGN(A34mid))
				&& (fabs(A13mid) < midz
					|| fabs(A24mid) < midz
					|| SIGN(A13mid) == SIGN(A24mid))
				&& (fabs(A14mid) < midz
					|| fabs(A23mid) < midz
					|| SIGN(A14mid) == SIGN(A23mid))) {
				mmessage(INFO_MSG, NO_ERROR, "H04a: impure "
					"subgenomic coverage (Proportion A: "
					"%.2f, %.2f; %.2f, %.2f; %.2f, %.2f)."
							"  Rejecting H04a.\n",
						A12mid + 0.5, A34mid + 0.5,
						A13mid + 0.5, A24mid + 0.5,
						A14mid + 0.5, A23mid + 0.5);

				/* look for another valid hypothesis */

				double A1mid = in->proportion_A[0] - 0.5;
				double A2mid = in->proportion_A[1] - 0.5;
				double A3mid = in->proportion_A[2] - 0.5;
				double A4mid = in->proportion_A[3] - 0.5;

				/* 3 pure clusters: set ACCEPT_THREE_A
				 * so no test of ACCEPT_THREE_B triggered
				 */
				if (fabs(A1mid) >= midz		/* 1, 2, */
					&& fabs(A2mid) >= midz	/* and 3 */
					&& fabs(A3mid) >= midz
					&& fabs(A4mid) < midz
					&& (SIGN(A1mid) != SIGN(A2mid)
						|| SIGN(A1mid) != SIGN(A3mid)
						|| SIGN(A2mid) != SIGN(A3mid))) {

					in->reject_state = ACCEPT_THREE_A;
					++in->proptest;
					in->excluded_clusters[in->sort_id[3]] = 1;
					n_excluded += exclude_cluster_reads(mh,
						in, in->sort_id[3]);
					fprintf(stderr, " (%u total from %u "
						"haplotypes)\n", n_excluded,
						in->K - 4 + in->proptest);
					mmessage(INFO_MSG, NO_ERROR, "Accepting"
							" H03a: p1, p2=p3.\n");
				} else if (fabs(A1mid) >= midz	/* 1, 2 */
					&& fabs(A2mid) >= midz	/* and 4 */
					&& fabs(A4mid) >= midz
					&& fabs(A3mid) < midz
					&& (SIGN(A1mid) != SIGN(A2mid)
						|| SIGN(A1mid) != SIGN(A4mid)
						|| SIGN(A2mid) != SIGN(A4mid))) {

					in->reject_state = ACCEPT_THREE_A;
					++in->proptest;
					in->excluded_clusters[in->sort_id[2]] = 1;
					n_excluded += exclude_cluster_reads(mh,
						in, in->sort_id[2]);
					fprintf(stderr, " (%u total from %u "
						"haplotypes)\n", n_excluded,
						in->K - 4 + in->proptest);
					mmessage(INFO_MSG, NO_ERROR, "Accepting"
							" H03a: p1, p2=p3.\n");
				} else if (fabs(A1mid) >= midz	/* 1, 3 */
					&& fabs(A3mid) >= midz	/* and 4 */
					&& fabs(A4mid) >= midz
					&& fabs(A2mid) < midz
					&& (SIGN(A1mid) != SIGN(A3mid)
						|| SIGN(A1mid) != SIGN(A4mid)
						|| SIGN(A3mid) != SIGN(A4mid))) {

					in->reject_state = ACCEPT_THREE_A;
					++in->proptest;
					in->excluded_clusters[in->sort_id[1]] = 1;
					n_excluded += exclude_cluster_reads(mh,
						in, in->sort_id[1]);
					fprintf(stderr, " (%u total from %u "
						"haplotypes)\n", n_excluded,
						in->K - 4 + in->proptest);
					mmessage(INFO_MSG, NO_ERROR, "Accepting"
							" H03a: p1, p2=p3.\n");
				} else if (fabs(A2mid) >= midz	/* 2, 3 */
					&& fabs(A3mid) >= midz	/* and 4 */
					&& fabs(A4mid) >= midz
					&& fabs(A1mid) < midz
					&& (SIGN(A2mid) != SIGN(A3mid)
						|| SIGN(A2mid) != SIGN(A4mid)
						|| SIGN(A3mid) != SIGN(A4mid))) {

					in->reject_state = ACCEPT_THREE_A;
					++in->proptest;
					in->excluded_clusters[in->sort_id[0]] = 1;
					n_excluded += exclude_cluster_reads(mh,
						in, in->sort_id[0]);
					fprintf(stderr, " (%u total from %u "
						"haplotypes)\n", n_excluded,
						in->K - 4 + in->proptest);
					mmessage(INFO_MSG, NO_ERROR, "Accepting"
							" H03a: p1, p2=p3.\n");

				/* two pure clusters */
				} else if (fabs(A1mid) >= midz	/* 1, 2 */
					&& fabs(A2mid) >= midz
					&& SIGN(A1mid) != SIGN(A2mid)) {

					in->reject_state = ACCEPT_TWO;
					in->proptest += 2;
					in->excluded_clusters[in->sort_id[2]] = 1;
					in->excluded_clusters[in->sort_id[3]] = 1;
					n_excluded += exclude_cluster_reads(mh,
						in, in->sort_id[2]);
					n_excluded += exclude_cluster_reads(mh,
						in, in->sort_id[3]);
					fprintf(stderr, " (%u total from %u "
						"haplotypes)\n", n_excluded,
						in->K - 4 + in->proptest);
					mmessage(INFO_MSG, NO_ERROR, "Accepting"
							" H02: p1, p2.\n");
				} else if (fabs(A1mid) >= midz	/* 1, 3 */
					&& fabs(A3mid) >= midz
					&& SIGN(A1mid) != SIGN(A3mid)) {

					in->reject_state = ACCEPT_TWO;
					in->proptest += 2;
					in->excluded_clusters[in->sort_id[1]] = 1;
					in->excluded_clusters[in->sort_id[3]] = 1;
					n_excluded += exclude_cluster_reads(mh,
						in, in->sort_id[1]);
					n_excluded += exclude_cluster_reads(mh,
						in, in->sort_id[3]);
					fprintf(stderr, " (%u total from %u "
						"haplotypes)\n", n_excluded,
						in->K - 4 + in->proptest);
					mmessage(INFO_MSG, NO_ERROR, "Accepting"
							" H02: p1, p2.\n");
				} else if (fabs(A1mid) >= midz	/* 1, 4 */
					&& fabs(A4mid) >= midz
					&& SIGN(A1mid) != SIGN(A4mid)) {

					in->reject_state = ACCEPT_TWO;
					in->proptest += 2;
					in->excluded_clusters[in->sort_id[1]] = 1;
					in->excluded_clusters[in->sort_id[2]] = 1;
					n_excluded += exclude_cluster_reads(mh,
						in, in->sort_id[1]);
					n_excluded += exclude_cluster_reads(mh,
						in, in->sort_id[2]);
					fprintf(stderr, " (%u total from %u "
						"haplotypes)\n", n_excluded,
						in->K - 4 + in->proptest);
					mmessage(INFO_MSG, NO_ERROR, "Accepting"
							" H02: p1, p2.\n");
				} else if (fabs(A2mid) >= midz	/* 2, 3 */
					&& fabs(A3mid) >= midz
					&& SIGN(A2mid) != SIGN(A3mid)) {

					in->reject_state = ACCEPT_TWO;
					in->proptest += 2;
					in->excluded_clusters[in->sort_id[0]] = 1;
					in->excluded_clusters[in->sort_id[3]] = 1;
					n_excluded += exclude_cluster_reads(mh,
						in, in->sort_id[0]);
					n_excluded += exclude_cluster_reads(mh,
						in, in->sort_id[3]);
					fprintf(stderr, " (%u total from %u "
						"haplotypes)\n", n_excluded,
						in->K - 4 + in->proptest);
					mmessage(INFO_MSG, NO_ERROR, "Accepting"
							" H02: p1, p2.\n");
				} else if (fabs(A2mid) >= midz	/* 2, 4 */
					&& fabs(A4mid) >= midz
					&& SIGN(A2mid) != SIGN(A4mid)) {

					in->reject_state = ACCEPT_TWO;
					in->proptest += 2;
					in->excluded_clusters[in->sort_id[0]] = 1;
					in->excluded_clusters[in->sort_id[2]] = 1;
					n_excluded += exclude_cluster_reads(mh,
						in, in->sort_id[0]);
					n_excluded += exclude_cluster_reads(mh,
						in, in->sort_id[2]);
					fprintf(stderr, " (%u total from %u "
						"haplotypes)\n", n_excluded,
						in->K - 4 + in->proptest);
					mmessage(INFO_MSG, NO_ERROR, "Accepting"
							" H02: p1, p2.\n");
				} else if (fabs(A3mid) >= midz	/* 3, 4 */
					&& fabs(A4mid) >= midz
					&& SIGN(A3mid) != SIGN(A4mid)) {

					in->reject_state = ACCEPT_TWO;
					in->proptest += 2;
					in->excluded_clusters[in->sort_id[0]] = 1;
					in->excluded_clusters[in->sort_id[1]] = 1;
					n_excluded += exclude_cluster_reads(mh,
						in, in->sort_id[0]);
					n_excluded += exclude_cluster_reads(mh,
						in, in->sort_id[1]);
					fprintf(stderr, " (%u total from %u "
						"haplotypes)\n", n_excluded,
						in->K - 4 + in->proptest);
					mmessage(INFO_MSG, NO_ERROR, "Accepting"
							" H02: p1, p2.\n");

				/* no good: <= 1 pure cluster */
				} else {

					return mmessage(ERROR_MSG,
						INTERNAL_ERROR, "Found fewer "
						"than two subgenomes with clean"
						" enough coverage: Refusing to "
						"genotype!\n");
				}

			} else {
				mmessage(INFO_MSG, NO_ERROR, "Accepting H04a: "
								"p1=p2=p3=p4.\n");
			}

		/* H03a: p1, p2 = p3 */
		} else if (in->reject_state == ACCEPT_THREE_A) {
			double A23mid = -0.5 + (in->coverage_A[1] + in->coverage_A[2]) /
				(in->cluster_sizes[in->included_clusters[1]]
				+ in->cluster_sizes[in->included_clusters[2]]);
			double A1mid = in->proportion_A[0] - 0.5;

			/* coverage purity test fails for H03a */
			if (fabs(A1mid) < midz
				|| fabs(A23mid) < midz
				|| SIGN(A1mid) == SIGN(A23mid)) {

				double A2mid = in->proportion_A[1] - 0.5;
				double A3mid = in->proportion_A[2] - 0.5;
				mmessage(INFO_MSG, NO_ERROR, "H03a: impure "
					"subgenomic coverage (Proportion A: "
					"%.2f %.2f).  Rejecting H03a.\n", 
					in->proportion_A[0], A23mid + 0.5);

				/* check for purity on two clusters */
				if (fabs(A1mid) >= midz	/* 1, 2 */
					&& fabs(A2mid) >= midz
					&& SIGN(A1mid) != SIGN(A2mid)) {

					mmessage(INFO_MSG, NO_ERROR,
						"Discarding third largest "
						"cluster (%u) and reads:",
							in->sort_id[2]);
					in->reject_state = ACCEPT_TWO;
					++in->proptest;
					in->excluded_clusters[in->sort_id[2]] = 1;
					n_excluded += exclude_cluster_reads(mh,
						in, in->sort_id[2]);
					fprintf(stderr, " (%u total from %u "
						"haplotypes)\n", n_excluded,
						in->K - 4 + in->proptest);
					mmessage(INFO_MSG, NO_ERROR, "Accepting"
							" H02: p1, p2.\n");

				} else if (fabs(A1mid) >= midz	/* 1, 3 */
					&& fabs(A3mid) >= midz
					&& SIGN(A1mid) != SIGN(A3mid)) {

					mmessage(INFO_MSG, NO_ERROR,
						"Discarding second largest "
						"cluster (%u) and reads:",
							in->sort_id[1]);
					in->reject_state = ACCEPT_TWO;
					++in->proptest;
					in->excluded_clusters[in->sort_id[1]] = 1;
					n_excluded += exclude_cluster_reads(mh,
						in, in->sort_id[1]);
					unsigned int k = in->included_clusters[1];
					in->included_clusters[1] = in->included_clusters[2];
					in->included_clusters[2] = k;
					fprintf(stderr, " (%u total from %u "
						"haplotypes)\n", n_excluded,
						in->K - 4 + in->proptest);
					mmessage(INFO_MSG, NO_ERROR, "Accepting"
							" H02: p1, p2.\n");

				/* fewer than two pure clusters */
				} else {

					return mmessage(ERROR_MSG,
						INTERNAL_ERROR, "Found fewer "
						"than two subgenomes with clean"
						" enough coverage: Refusing to "
						"genotype!\n");
				}
			} else {
				mmessage(INFO_MSG, NO_ERROR, "Accepting H03a: "
								"p1, p2=p3.\n");
			}
		} else if (in->reject_state == ACCEPT_TWO) {
			double A1mid = in->proportion_A[0] - 0.5;
			double A2mid = in->proportion_A[1] - 0.5;

			/* coverage purity test fails for H02 */
			if (fabs(A1mid) < midz
				|| fabs(A2mid) < midz
				|| SIGN(A1mid) == SIGN(A2mid)) {

				return mmessage(ERROR_MSG, INTERNAL_ERROR,
					"Found fewer than two clusters with "
					"clean subgenomic coverage:  Refusing "
					"to genotype!\n");
			} else {
				mmessage(INFO_MSG, NO_ERROR, "Accepting H02: "
								"p1, p2.\n");
			}
		}


		/* H03b: p1 = p2, p3 */
		if (in->reject_state == ACCEPT_THREE_B) {
			double p12A = (in->coverage_A[0] + in->coverage_A[1]) /
				(in->cluster_sizes[in->included_clusters[0]]
					+ in->cluster_sizes[in->included_clusters[1]]);
			double A12mid = p12A - 0.5;
			double A3mid = in->proportion_A[2] - 0.5;

			/* coverage purity test fails for H03b */
			if (fabs(A12mid) < midz || fabs(A3mid) < midz
				|| SIGN(A12mid) == SIGN(A3mid)) {

				double A1mid = in->proportion_A[0] - 0.5;
				double A2mid = in->proportion_A[1] - 0.5;
				mmessage(INFO_MSG, NO_ERROR, "H03b: bad "
					"subgenomic coverage (Proportion A: "
					"%.2f %.2f).  Rejecting H03b.\n",
						A12mid + 0.5, A3mid + 0.5);

				/* cluster 1 and 2 are pure and cover A and B */
				if (fabs(A1mid) >= midz && fabs(A2mid) >= midz
					&& SIGN(A1mid) != SIGN(A2mid)) {

					mmessage(INFO_MSG, NO_ERROR,
						"Discarding third largest "
						"cluster (%u) and reads: ",
						in->sort_id[2]);
					in->reject_state = ACCEPT_TWO;
					++in->proptest;
					in->excluded_clusters[in->sort_id[2]] = 1;
					n_excluded += exclude_cluster_reads(mh,
						in, in->sort_id[2]);
					fprintf(stderr, " (%u total from %u haplotypes)\n",
						n_excluded, in->K - 4 + in->proptest);
				} else {

					return mmessage(ERROR_MSG,
						INTERNAL_ERROR, "Found fewer "
						"than two subgenomes with clean"
						" enough coverage: Refusing to "
						"genotype!\n");
				}
			} else {
				mmessage(INFO_MSG, NO_ERROR, "Accepting H03b: "
								"p1=p2, p3.\n");
			}
		}

		mmessage(INFO_MSG, NO_ERROR, "Keeping final clusters:");
		fprint_uints(stderr, in->included_clusters, MIN(in->K, 
						4 - in->proptest), 1, 1);

		/* no more to do if we are genotyping by clustering */
		if (opt.genotype_by_clustering) {
			genotype_by_clustering(&opt, in, mh, sds, ref_seqs, ref_lengths);
			goto EXIT_NOW;
		}

		/* n_read is number of reads prior to exclusions for paralogs */

		/* commit to exclude n_excluded paralog reads */
		if (n_excluded) {
			
			size_t n_idx1 = 0;	/* 0 .. n_read */
			size_t n_idx2 = 0;	/* 0 .. n_read - n_excluded */
			double **new_pp = malloc(N_SUBGENOMES * sizeof **pp);
			double **new_ll = malloc(N_SUBGENOMES * sizeof **ll);

			if (!new_pp || !new_ll)
				return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
					"posterior probability.\n");
			n_read -= n_excluded;
			for (int j = 0; j < N_SUBGENOMES; ++j) {
				new_pp[j] = malloc(n_read * sizeof *new_pp);
				new_ll[j] = malloc(n_read * sizeof *new_ll);
			}
			A_expected_coverage = 0;
			if (opt.min_log_likelihood > 0)
				fprintf(stderr, "Maximum log likelihoods:");
			for (merge_hash *me = mh; me != NULL; me = me->hh.next) {
				if (me->exclude)
					continue;
				if (in->assignment[n_idx1] == NA_ASSIGNMENT
					|| in->excluded_clusters[in->assignment[n_idx1]]) {
					me->exclude = 1;
				} else {
					for (int j = 0; j < N_SUBGENOMES; ++j) {
						new_pp[j][n_idx2] = pp[j][n_idx1];
						new_ll[j][n_idx2] = ll[j][n_idx1];
					}
					/*sam_entry *se = &sds[0]->se[me->indices[0][0]];
					mmessage(INFO_MSG, NO_ERROR, "@%s ppA=%.6e\n", se->name, new_pp[0][n_idx2]);*/
					A_expected_coverage += new_pp[0][n_idx2];
					if (opt.min_log_likelihood > 0)
						fprintf(stderr, " %f", mll[n_idx1]);
					++n_idx2;
				}
				++n_idx1;
			}
			if (opt.min_log_likelihood > 0)
				fprintf(stderr, "\n");
			for (int j = 0; j < N_SUBGENOMES; ++j) {
				free(pp[j]);
				free(ll[j]);
				pp[j] = new_pp[j];
				ll[j] = new_ll[j];
			}
			free(new_pp);
			free(new_ll);
			B_expected_coverage = n_read - A_expected_coverage;
			min_expected_coverage = MIN(A_expected_coverage,
							B_expected_coverage);
			mmessage(INFO_MSG, NO_ERROR, "Updated expected coverage"
					": %f %f (%zu)\n", A_expected_coverage,
						B_expected_coverage, n_read);
		}

		if (opt.fastq_after_filter || opt.sam_after_filter[0]) {
			FILE *ffp[N_SUBGENOMES] = {fopen(
				opt.fastq_after_filter ?  opt.fastq_after_filter
					: opt.sam_after_filter[0], "w")};

			for (unsigned int i = 0; i < N_SUBGENOMES; ++i)
				ffp[i] = fopen(opt.sam_after_filter[i], "w");

			if (!ffp[0] || (opt.sam_after_filter[1] && !ffp[1]))
				return mmessage(ERROR_MSG, INTERNAL_ERROR,
					"Could not open file '%s' for "
					"writing.\n", opt.fastq_after_filter
					|| !ffp[0] ? opt.sam_after_filter[0]
						: opt.sam_after_filter[1]);

			n_included_reads = 0;

			for (merge_hash *me = mh; me != NULL; me = me->hh.next) {

				if (me->exclude)
					continue;

				if (opt.fastq_after_filter) {
					sam_entry *se
						= &sds[0]->se[me->indices[0][0]];
					fprintf(ffp[0], "@%s ppA=%.6e\n",
								se->name,
						pp[0][n_included_reads]);
					fwrite_nuc_segment(ffp[0], se->read,
						XY_ENCODING, 0, se->read->len);
					fprintf(ffp[0], "\n+\n");
					fwrite_qual_sequence(ffp[0], se->qual);
					fprintf(ffp[0], "\n");
				} else {
					for (unsigned int i = 0; i
							< N_SUBGENOMES; ++i) {
						sam_entry *se = &sds[i]->se[
							me->indices[i][0]];
						write_sam_entry(ffp[i], sds[i], se);
					}
				}
				++n_included_reads;
			}

			fclose(ffp[0]);
			for (unsigned int i = 1; i < N_SUBGENOMES; ++i)
				fclose(ffp[i]);
		}
	}

	if (min_expected_coverage < opt.min_expected_coverage)
		return mmessage(ERROR_MSG, INTERNAL_ERROR,
			"Sugenomic coverage is too low: A=%.2f B=%.2f.\n",
			A_expected_coverage, B_expected_coverage);

	if (opt.min_log_likelihood > 0)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "Stopping after "
					"writing maximum log likelihoods.\n");

	/* store pileup information */
	xy_t *obs_nuc = malloc(n_read * sizeof *obs_nuc);	/* bases */
	qual_t *obs_q = malloc(n_read * sizeof *obs_q);		/* qualities */
	char *genome_src = malloc(n_read * sizeof *genome_src);	/* subgenome */
	char *covers = malloc(n_read * sizeof *covers);		/* gap? */

	/* we need to make sure the same read base aligns to the homoeologous
	 * positions; rd_idxA is the read position aligned to current reference
	 * position as per subgenome A alignment; obs_rpos is the shared read
	 * position aligned to the current reference position in the subset of
	 * reads where both alignments agree on the aligned read base
	 */
	unsigned int *rd_idxA = malloc(n_read * sizeof *rd_idxA);
	unsigned int *obs_rpos = malloc(n_read * sizeof *obs_rpos);

	/* store nucleotide counts for determining ref and alt alleles */
	size_t num_nuc[NUM_NUCLEOTIDES];	/* total number of each nuc */
	xy_t nuc1, nuc2 = 0, nuc3 = 0;		/* top-most abundant 3 nucs */
	size_t num_nuc1, num_nuc2, num_nuc3;	/* and their counts */

	/* record expected nucleotide counts in each subgenome; only expected
	 * because we do not know the true alignment
	 */
	double ebaseA[NUM_NUCLEOTIDES];		/* expected in subgenome A */
	double ebaseB[NUM_NUCLEOTIDES];		/* expected in subgenome B */

	/* Use the commented code to focus on a particular subset of sites for debugging */
	/*
	 start_pos[0] = 70;
	 end_pos[0] = 80;
	 start_pos[1] = 70;
	 end_pos[1] = 80;
	 */

	size_t region_len = end_pos[1] - start_pos[1] + 1;
	if (end_pos[0] - start_pos[0] > region_len)
		region_len = end_pos[0] - start_pos[0] + 1;

	/* store information for post-hoc tests of allelic coverage */
	/* [BUG,KSD] no memory allocation checks */
	size_t *haplotype_posns[N_SUBGENOMES];	/* positions of het calls */
	haplotype_posns[0] = calloc(region_len, sizeof *haplotype_posns[0]);
	haplotype_posns[1] = calloc(region_len, sizeof *haplotype_posns[1]);

	/* proportion of reads matching dominant haplotype in each subgenome */
	double *hapA_prop = malloc(region_len * sizeof *hapA_prop);
	double *hapB_prop = malloc(region_len * sizeof *hapB_prop);

	/* coverage of dominant haplotype in each subgenome */
	double *hapA_covg = malloc(region_len * sizeof *hapA_covg);
	double *hapB_covg = malloc(region_len * sizeof *hapB_covg);

	/* dominant nucleotides at each position of dominant haplotype */
	xy_t *hapA_dom_nuc = calloc(region_len, sizeof *hapA_dom_nuc);
	xy_t *hapB_dom_nuc = calloc(region_len, sizeof *hapB_dom_nuc);

	/* number of segregating sites in each subgenome */
	unsigned int n_segregatingA = 0, n_segregatingB = 0;

	/* open and write header of vcf files */
	FILE *vcf_fp[N_SUBGENOMES];
	for (unsigned int i = 0; i < N_SUBGENOMES; ++i) {
		vcf_fp[i] = NULL;
		if (opt.vcf_files[i])
			vcf_fp[i] = fopen(opt.vcf_files[i], "w");
		else
			continue;
		if (!vcf_fp[i]) {
			mmessage(ERROR_MSG, FILE_OPEN_ERROR,
							opt.vcf_files[i]);
			goto EXIT_AFTER_GENOTYPING_STARTS;
		}
		print_vcf_header(vcf_fp[i], opt.vcf_opt, opt.fsa_files[i],
							opt.sample_name);
	}

	/* local debug level */
	debug_level = DEBUG_I;

	/* commit to allotetraploid */

	/* finally: march along reference positions and genotype */
	for (size_t posA = start_pos[0], posB = start_pos[1];
		posA <= end_pos[0] || posB <= end_pos[1]; ++posA, ++posB) {

		double prob_heterozygote[N_SUBGENOMES] = {0, 0};
		size_t ref_pos[N_SUBGENOMES] = {posA, posB};
		char_t ref_base[N_SUBGENOMES] = {IUPAC_A, IUPAC_A};
		int no_alt_allele = 0;

#ifdef ALLOW_DEBUGGING
if (opt.debugging_site >= 0 && opt.debugging_site < (int)posA)
	exit(0);
#endif

		/* count each of alleles across all reads */
		for (int b = 0; b < NUM_NUCLEOTIDES; ++b) {
			num_nuc[b] = 0;
			ebaseA[b] = 0;
			ebaseB[b] = 0;
		}

		/* extract read position aligned to genome A */
		size_t target_a = posA;
		size_t target_b = posB;
		size_t site = target_a - start_pos[0];
		iupac_t ref_allele[N_SUBGENOMES] = {
			fds[0]->reads[fs_index[0] + posA], 
			fds[1]->reads[fs_index[1] + posB]};

		debug_msg(debug_level > ABSOLUTE_SILENCE, debug_level,
			  "Site %zu (target_a = %zu; target_b = %zu, "
			  "start_pos[0] = %u)\n", site, target_a, target_b, start_pos[0]);

		n_read = 0;
		size_t n_cover = 0;
		for (merge_hash *me = mh; me != NULL; me = me->hh.next) {

			if (me->exclude)
				continue;

			/* A alignment */
			sam_entry *se = &sds[0]->se[me->indices[0][0]];
			unsigned int rd_idx = 0;
			size_t rf_idx = se->pos - 1;

			if (rf_idx > target_a) {
				++n_read;
				continue;
			}

			for (unsigned int j = 0; j < se->cig->n_ashes; ++j) {

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
					ref_base[0] = fds[0]->reads[fs_index[0] + rd_idx + target_a - rf_idx];
					/*
					 debug_msg(debug_level > ABSOLUTE_SILENCE,
					 debug_level, "Read %u (A): rd_idx = "
					 "%u\n", n_read, rd_idxA[n_read]);
					 */
					/*
					 debug_msg_cont(debug_level > QUIET, debug_level,
					 "%c w/ quality %u",
					 xy_to_char[obs_nuc[n_cover]],
					 obs_q[n_cover]);
					 */
					n_cover++;
					break;
				}

				rf_idx += se->cig->ashes[j].len;
				rd_idx += se->cig->ashes[j].len;
			}
			++n_read;
		}

		debug_msg(debug_level > ABSOLUTE_SILENCE, debug_level,
			  "Site %zu (target_a = %zu; target_b = %zu, "
			  "start_pos[0] = %u): A coverage = %zu (%zu)\n", site,
			  target_a, target_b, start_pos[0], n_cover, n_read);

		/* no reads cover this site after accounting for soft-clipping */
		if (!n_cover)
			continue;

		/* extract set of nucleotides aligned to this position in
		 * BOTH alignments
		 */
		n_cover = 0;
		n_read = 0;

		for (merge_hash *me = mh; me != NULL; me = me->hh.next) {

			if (me->exclude)
				continue;

			/* B alignment */
			sam_entry *se = &sds[1]->se[me->indices[1][0]];
			unsigned int rd_idx = 0;
			size_t rf_idx = se->pos - 1;
			covers[n_read] = 0;

			if (rf_idx > target_b) {
				++n_read;
				continue;
			}

			for (unsigned int j = 0; j < se->cig->n_ashes; ++j) {

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
				/*
				if (rf_idx + se->cig->ashes[j].len > target_b)
					debug_msg(debug_level > ABSOLUTE_SILENCE,
						debug_level, "Read %u (B): rd_idx = "
						"%u\n", n_read,
						rd_idx + target_b - rf_idx);
				 */
				if (rf_idx + se->cig->ashes[j].len > target_b
				    && rd_idxA[n_read] == rd_idx + target_b - rf_idx) {

					covers[n_read] = 1;
					obs_rpos[n_cover] = rd_idx + target_b - rf_idx;
					obs_nuc[n_cover] = get_nuc(se->read,
								   XY_ENCODING, obs_rpos[n_cover]);
					ref_base[1] = fds[1]->reads[fs_index[1] + rd_idx + target_b - rf_idx];
					ebaseA[obs_nuc[n_cover]] += pp[0][n_read];
					ebaseB[obs_nuc[n_cover]] += 1 - pp[0][n_read];
					++num_nuc[obs_nuc[n_cover]];
					obs_q[n_cover] = get_qual(se->qual,
								  obs_rpos[n_cover]);
					/*
					 fprintf(stderr, "Read %zu, Nuc %c (%u) at %u from A wp %e, from B wp %e (%e; %e)\n", n_read, xy_to_char[obs_nuc[n_cover]], obs_q[n_cover], obs_rpos[n_cover], pp[0][n_read], 1 - pp[0][n_read], ebaseA[obs_nuc[n_cover]], ebaseB[obs_nuc[n_cover]]);
					if (!covers[n_read]) {
						++num_nuc[obs_nuc[n_cover]];
						obs_q[n_cover] = get_qual(se->qual,
						obs_rpos[n_cover]);
						covers[n_read] = 1;
					}
					debug_msg_cont(debug_level > QUIET, debug_level,
						"%c", xy_to_char[get_nuc(se->read,
						XY_ENCODING, obs_rpos[n_cover]);
					 */

					//				if (!me->count[0]) {	// pick up observations not aligned to genome A
					//					obs_nuc[n_cover] = xy_to_rc[get_nuc(se->read, XY_ENCODING, rd_idx + target_b - rf_idx)];
					//					obs_q[n_cover] = get_qual(se->qual, rd_idx + target_b - rf_idx);
					//					++num_nuc[obs_nuc[n_cover]];
					//					covers[n_read] = 1;
					//					++n_cover;
					//					++n_read;
					//				}
					++n_cover;
					break;
				}

				rf_idx += se->cig->ashes[j].len;
				rd_idx += se->cig->ashes[j].len;
			}
			++n_read;
		}

		//debug_msg_cont(debug_level > QUIET, debug_level, "\n");

		debug_msg(debug_level > ABSOLUTE_SILENCE, debug_level,
			  "Site %zu (target_a = %zu; target_b = %zu, "
			  "start_pos[0] = %u): B coverage = %zu (%zu)\n", site,
			  target_a, target_b, start_pos[0], n_cover, n_read);
		
		debug_msg(debug_level > QUIET, debug_level,
						"Expected counts genome A:");
		debug_call(debug_level > QUIET, debug_level,
			fprint_doubles(stderr, ebaseA, NUM_NUCLEOTIDES, 3, 1));
		debug_msg(debug_level > QUIET, debug_level,
						"Expected counts genome B:");
		debug_call(debug_level > QUIET, debug_level,
			fprint_doubles(stderr, ebaseB, NUM_NUCLEOTIDES, 3, 1));

		/* mode nucleotide in each subgenome */
		double ecoverage[N_SUBGENOMES] = {0, 0};	/* subgenomic exp. coverage */
		xy_t modeA = XY_A, modeB = XY_A;
		double mode_heightA = 0, mode_heightB = 0;
		for (unsigned int i = 0; i < NUM_NUCLEOTIDES; ++i) {
			ecoverage[0] += ebaseA[i];
			ecoverage[1] += ebaseB[i];
			if (ebaseA[i] > mode_heightA) {
				mode_heightA = ebaseA[i];
				modeA = i;
			}
			if (ebaseB[i] > mode_heightB) {
				mode_heightB = ebaseB[i];
				modeB = i;
			}
		}
		
		debug_msg(debug_level > ABSOLUTE_SILENCE, debug_level,
			"A %f; B %f; ratios %f; %f (%f)\n", ecoverage[0], 
			ecoverage[1], ecoverage[0]/A_expected_coverage,
			ecoverage[1]/B_expected_coverage, opt.coverage_screen);

		/* screen subgenomic coverage relative to subgenomic read
		 * assignments from alignment; low number indicates likely
		 * indel difference between subgenomes imposed by our external
		 * alignment of the references; we do not genotype such sites
		 * as it is either regular diploid genotyping or there are
		 * issues of indel variants which requires further thought
		 */
		if (ecoverage[0] / A_expected_coverage < opt.coverage_screen ||
			ecoverage[1] / B_expected_coverage < opt.coverage_screen) {

			debug_msg(debug_level > QUIET, debug_level, "Coverage "
				"rate is too low, will not genotype site!\n");
			for (int i = 0; i < N_SUBGENOMES; ++i) {
				if (!vcf_fp[i])
					continue;
				fprintf(vcf_fp[i], "%s\t%lu\t.\t%c\t.\t.\t"
					"c%2.0f\t.\t.\t.\n", opt.ref_names[i],
					ref_pos[i] + 1,
					iupac_to_char[ref_allele[i]],
					100*opt.coverage_screen);
			}
			continue;
		}

		/* identify three most abundant alleles */
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

		/* if only one A allele in reads, choose C as alternate just
		 * for genotyping: this is like proposing a straw man, which is
		 * fine when using the default error model where nucleotide
		 * identity is not important
		 */
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
				"site (%c=%zu, %c=%zu, %c=%zu).  To change the "
				"behavior, see command-line option "
				"--biallelic\n", xy_to_char[nuc1], num_nuc1,
						xy_to_char[nuc2], num_nuc2,
						xy_to_char[nuc3], num_nuc3);
			for (int i = 0; i < N_SUBGENOMES; ++i) {
				int one_alt = 0;

				if (!vcf_fp[i])
					continue;

				fprintf(vcf_fp[i], "%s\t%lu\t.\t%c\t",
					opt.ref_names[i], ref_pos[i] + 1,
					iupac_to_char[ref_allele[i]]);

				if (xy_to_iupac[nuc1] != ref_allele[i]) {
					fprintf(vcf_fp[i], "%c",
						xy_to_char[nuc1]);
					one_alt = 1;
				}
				if (!no_alt_allele && nuc2 != nuc1
					&& xy_to_iupac[nuc2] != ref_allele[i]) {
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

		/* finally, we prepare to genotype */
		debug_msg(debug_level > QUIET, debug_level,
			  "Observed nucleotides (%u): ", n_cover);
		for (unsigned int i = 0; i < n_cover; ++i)
			debug_msg_cont(debug_level > QUIET, debug_level, "%c",
				       xy_to_char[obs_nuc[i]]);
		debug_msg_cont(debug_level > QUIET, debug_level, "\n");
		debug_msg(debug_level > QUIET, debug_level,
			  "Observed   qualities (%u): ", n_cover);
		for (unsigned int i = 0; i < n_cover; ++i)
			debug_msg_cont(debug_level > QUIET, debug_level,
				       "%c", (char)(obs_q[i] + MIN_ASCII_QUALITY_SCORE));
		debug_msg_cont(debug_level > QUIET, debug_level, "\n");

		/* estimate posterior probability of each genotype */

		double lprob[9], gprob[9] = {0,0,0,0,0,0,0,0,0};
		double avg_lprob[9] = {0,0,0,0,0,0,0,0,0};
		int g_max[N_SUBGENOMES] = {0, 0};
		double mprob = 0;

		if (karin_version && !opt.karin_old) {/* to replace old version */

		double tmp1, tmp2;
		double max = -INFINITY, den = 0;

		/* consider each genotype at the current locus */
		for (int g1 = 0; g1 <= 2; ++g1) {
			for (int g2 = 0; g2 <= 2; ++g2) {

				n_read = 0;
				n_cover = 0;
				lprob[g1 * 3 + g2] = 0;

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

					/* log likelihood of subgenomic A
					 * alignment is pre-computed
					 * log-likelihood minus the old
					 * plus the new contribution
					 * from current site
					 */
					tmp1 = ll[0][n_read]
						- sub_prob_given_q_with_encoding(
							ref_base[0],			/* homozygous ref base */
							obs_nuc[n_cover], IUPAC_ENCODING, XY_ENCODING, obs_q[n_cover], 1, (void *)&mls)
						+ sub_prob_given_q_with_encoding(
							 !g1 ? xy_to_iupac[nuc1]	/* new genotype */
							: g1 == 2 ? xy_to_iupac[nuc2]
							: xy_to_iupac[nuc1] | xy_to_iupac[nuc2],
							obs_nuc[n_cover], IUPAC_ENCODING, XY_ENCODING, obs_q[n_cover], 1, (void *)&mls);
#ifdef ALLOW_DEBUGGING
	if ((int) posA == opt.debugging_site)
		fprintf(stderr, "g1=%d, g2=%d (%c%c), rA=%c (%e) -> %c (%d): %f %f Delta = %f\n",
			g1, g2, xy_to_char[nuc1], xy_to_char[nuc2], iupac_to_char[ref_base[0]],
			pp[0][n_read], xy_to_char[obs_nuc[n_cover]], obs_q[n_cover],
			ll[0][n_read], ll[1][n_read],
			sub_prob_given_q_with_encoding(
				!g1 ? xy_to_iupac[nuc1] : g1 == 2 ? xy_to_iupac[nuc2]
				: xy_to_iupac[nuc1] | xy_to_iupac[nuc2],
				obs_nuc[n_cover], IUPAC_ENCODING, XY_ENCODING, obs_q[n_cover], 1, (void *)&mls)
			- sub_prob_given_q_with_encoding(
				ref_base[0], obs_nuc[n_cover], IUPAC_ENCODING, XY_ENCODING, obs_q[n_cover], 1, (void *)&mls));
#endif

					/* log likelihood of B alignment */
					tmp2 = ll[1][n_read]
						- sub_prob_given_q_with_encoding(
							ref_base[1],			/* homozygous ref base */
							obs_nuc[n_cover], IUPAC_ENCODING, XY_ENCODING, obs_q[n_cover], 1, (void *)&mls)
						+ sub_prob_given_q_with_encoding(
							 !g2 ? xy_to_iupac[nuc1]	/* new genotype */
							: g2 == 2 ? xy_to_iupac[nuc2]
							: xy_to_iupac[nuc1] | xy_to_iupac[nuc2],
							obs_nuc[n_cover], IUPAC_ENCODING, XY_ENCODING, obs_q[n_cover], 1, (void *)&mls);

					/* combine assuming uniform prior */
					lprob[g1 * 3 + g2] += log(exp(tmp1) + exp(tmp2));

#ifdef ALLOW_DEBUGGING
	if ((int) posA == opt.debugging_site)
		fprintf(stderr, "g1=%d, g2=%d (%c%c), rB=%c (%e) -> %c (%d): %f %f Delta = %f : %f\n",
			g1, g2, xy_to_char[nuc1], xy_to_char[nuc2], iupac_to_char[ref_base[1]],
			pp[1][n_read], xy_to_char[obs_nuc[n_cover]], obs_q[n_cover],
			ll[0][n_read], ll[1][n_read],
			sub_prob_given_q_with_encoding(
				!g2 ? xy_to_iupac[nuc1]
				: g2 == 2 ? xy_to_iupac[nuc2]
				: xy_to_iupac[nuc1] | xy_to_iupac[nuc2],
				obs_nuc[n_cover], IUPAC_ENCODING, XY_ENCODING, obs_q[n_cover], 1, (void *)&mls) - sub_prob_given_q_with_encoding(
				ref_base[1], obs_nuc[n_cover], IUPAC_ENCODING, XY_ENCODING, obs_q[n_cover], 1, (void *)&mls), lprob[g1 * 3 + g2]);
#endif


					++n_cover;
					++n_read;
				}
				if (max < lprob[g1 * 3 + g2])
					max = lprob[g1 * 3 + g2];
			}
		}

		/* normalize */
		for (int g1 = 0; g1 <= 2; ++g1)
			for (int g2 = 0; g2 <= 2; ++g2) {
				gprob[g1 * 3 + g2] = exp(lprob[g1 * 3 + g2] - max);
				den += gprob[g1 * 3 + g2];
			}

		for (int g1 = 0; g1 <= 2; ++g1)
			for (int g2 = 0; g2 <= 2; ++g2) {
				gprob[g1 * 3 + g2] /= den;
 				fprintf(stderr, "M = (%u, %u) %c%c at site "
					"(%zu, %zu): %e\n", g1, g2,
					xy_to_char[nuc1], xy_to_char[nuc2],
					posA + 1, posB + 1, gprob[g1*3 + g2]);
				if (gprob[g1 * 3 + g2] > mprob) {
					mprob = gprob[g1 * 3 + g2];
					g_max[0] = g1;
					g_max[1] = g2;
				}
			}
		} else {	/* end new genotyping (karin_version) */

		for (unsigned int b = 0; b < opt.n_sample; ++b) {
			double max = -INFINITY, den = 0;

			/* repeat simulation */

			/* simulate alignments (or genome source) */
			n_read = 0;
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
			/* compute Pr(M|R,A) */
			for (int g1 = 0; g1 <= 2; ++g1) {
				for (int g2 = 0; g2 <= 2; ++g2) {

					/* compute likelihood of this genotype */
					n_read = 0;
					n_cover = 0;
					lprob[g1*3 + g2] = 0;	/* uniform prior on genotypes */
					for (merge_hash *me = mh; me != NULL; me = me->hh.next) {
						if (me->exclude)
							continue;
						if (covers[n_read]) {
							mls.pos = obs_rpos[n_cover];
							/* assume source is A genome */
							if (genome_src[n_cover] == 'A') {
								lprob[g1*3 + g2] += sub_prob_given_q_with_encoding(
									     !g1 ? xy_to_iupac[nuc1]			/* A genome is nuc1nuc1 */
									     : g1 == 2 ? xy_to_iupac[nuc2]		/* A genome is nuc2nuc2 */
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
							/* assume source is B genome */
							} else {
								lprob[g1*3 + g2] += sub_prob_given_q_with_encoding(
									!g2 ? xy_to_iupac[nuc1]				/* B genome is nuc1nuc1 */
									: g2 == 2 ? xy_to_iupac[nuc2]			/* B genome is nuc2nuc2 */
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
							}
							++n_cover;
						}
						++n_read;
					}
					if (max < lprob[g1*3 + g2])
						max = lprob[g1*3 + g2];
#ifdef ALLOW_DEBUGGING
					if ((int) posA == opt.debugging_site)
						fprintf(stderr, "%u (%u, %u) %f\n", b, g1, g2, lprob[g1*3 + g2]);
#endif

				}
			}

			/* normalize */
			for (int g1 = 0; g1 <= 2; ++g1)
				for (int g2 = 0; g2 <= 2; ++g2) {
					avg_lprob[g1*3 + g2] += lprob[g1*3 + g2];
					lprob[g1*3 + g2] = exp(lprob[g1*3 + g2] - max);
					den += lprob[g1*3 + g2];
				}

			for (int g1 = 0; g1 <= 2; ++g1)
				for (int g2 = 0; g2 <= 2; ++g2) {
					lprob[g1*3 + g2] /= den;
					gprob[g1*3 + g2] += lprob[g1*3 + g2];
				}

		}

		/* average across Monte Carlo samples */
		for (int g1 = 0; g1 <= 2; ++g1)
			for (int g2 = 0; g2 <= 2; ++g2) {
				gprob[g1*3 + g2] /= opt.n_sample;
				avg_lprob[g1*3 + g2] /= opt.n_sample;
 				fprintf(stderr, "M = (%u, %u) %c%c at site (%zu, %zu): %e\n", g1, g2, xy_to_char[nuc1], xy_to_char[nuc2], posA + 1, posB + 1, gprob[g1*3 + g2]);
				if (gprob[g1*3 + g2] > mprob) {
					mprob = gprob[g1*3 + g2];
					g_max[0] = g1;
					g_max[1] = g2;
				}
			}
		} /* end original genotype */

		/* stderr output */
		fprintf(stderr, "Genotype (%4zu, %4zu, %3zu, %3zu): %c%c/%c%c (%f) [",
			posA + 1, posB + 1, posA - start_pos[0], posB - start_pos[1],
			g_max[0] < 2 ? xy_to_char[nuc1] : xy_to_char[nuc2],
			g_max[0] ? xy_to_char[nuc2] : xy_to_char[nuc1],
			g_max[1] < 2 ? xy_to_char[nuc1] : xy_to_char[nuc2],
			g_max[1] ? xy_to_char[nuc2] : xy_to_char[nuc1], mprob);

		/* compute heterozygote probability */
		for (int g1 = 0; g1 <= 2; ++g1)
			for (int g2 = 0; g2 <= 2; ++g2) {
				fprintf(stderr, " %f", gprob[g1*3 + g2]);
				if (g1 == 1)
					prob_heterozygote[0] += gprob[g1*3 + g2];
				if (g2 == 1)
					prob_heterozygote[1] += gprob[g1*3 + g2];
			}

		/* extra flare to identify allelic and homoeologous SNPs */
		fprintf(stderr, "]%s\n", g_max[0] == 1 || g_max[1] == 1 ? "***"
			: abs(g_max[0] - g_max[1]) > 1 ? "+++" : "");

		double ect_pvals[N_SUBGENOMES] = {1, 1};
		if (opt.equal_homolog_coverage_test
			&& (g_max[0] == 1 || g_max[1] == 1))
				test_equal_homolog_coverage(mh, ll, ref_base,
					covers, obs_nuc, obs_q, obs_rpos, g_max,
					nuc1, nuc2, debug_level, ect_pvals);

		/* write out results to vcf files */
		for (int i = 0; i < N_SUBGENOMES; ++i) {
			int one_alt = 0;

			if (!vcf_fp[i])
				continue;

			double pe = i
				? (1 - gprob[g_max[i]] - gprob[3 + g_max[i]] - gprob[6 + g_max[i]])
				: (1 - gprob[3 * g_max[i] + 2] - gprob[3 * g_max[i] + 1] - gprob[3 * g_max[i]]);

			fprintf(vcf_fp[i], "%s\t%lu\t.\t%c\t",
				opt.ref_names[i], ref_pos[i] + 1,
				iupac_to_char[ref_allele[i]]);

			if (xy_to_iupac[nuc1] != ref_allele[i]) {
				fprintf(vcf_fp[i], "%c", xy_to_char[nuc1]);
				one_alt = 1;
			}
			if (!no_alt_allele && nuc2 != nuc1
				&& xy_to_iupac[nuc2] != ref_allele[i]) {
				if (one_alt)
					fputc(',', vcf_fp[i]);
				fprintf(vcf_fp[i], "%c", xy_to_char[nuc2]);
				one_alt = 1;
			}
			if (!one_alt)
				fputc('.', vcf_fp[i]);
			/* NOTE: Currently we are not output ALT alleles
			 * from the other subgenome.
			 */
			if (opt.posthoc_coverage_test && g_max[i] == 1 &&
				prob_heterozygote[i] < opt.phc_min_genotype_pp) {
				int mgq = opt.phc_min_genotype_pp < 1
					? MIN(99, (int) (-10 * log10(
					1 - opt.phc_min_genotype_pp))) : 99;
				fprintf(vcf_fp[i], "\t.\tgq%d\t.\tGT:DP:GQ",
								mgq);
			} else {
				fprintf(vcf_fp[i], "\t.\tPASS\t.\tGT:DP:GQ");
			}
			if (opt.equal_homolog_coverage_test && g_max[i] == 1)
				fprintf(vcf_fp[i], ":ET");
			if (opt.vcf_opt->output_gl)
				fprintf(vcf_fp[i], ":GL");

			fputc('\t', vcf_fp[i]);

			/* reference allele is dominant allele */
			if (xy_to_iupac[nuc1] == ref_allele[i]) {
				if (g_max[i] == 2)
					fprintf(vcf_fp[i], "1/1");
				else if (g_max[i] == 1)
					fprintf(vcf_fp[i], "0/1");
				else
					fprintf(vcf_fp[i], "0/0");

			/* reference allele is subdominant allele */
			} else if (xy_to_iupac[nuc2] == ref_allele[i]) {
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
			if (xy_to_iupac[nuc1] == ref_allele[i]) {
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
			} else if (xy_to_iupac[nuc2] == ref_allele[i]) {
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
		if (prob_heterozygote[0] >= opt.phc_min_genotype_pp) {
			haplotype_posns[0][n_segregatingA] = posA;
			hapA_prop[n_segregatingA] = ebaseA[modeA] / ecoverage[0];
			hapA_covg[n_segregatingA] = ecoverage[0];
			hapA_dom_nuc[n_segregatingA++] = modeA;
		}
		if (prob_heterozygote[1] >=  opt.phc_min_genotype_pp) {
			haplotype_posns[1][n_segregatingB] = posB;
			hapB_prop[n_segregatingB] = ebaseB[modeB] / ecoverage[1];
			hapB_covg[n_segregatingB] = ecoverage[1];
			hapB_dom_nuc[n_segregatingB++] = modeB;
		}
	}

/*
	fprintf(stderr, "Genotype A has %u segregating sites:", n_segregatingA);
	for (unsigned int i = 0; i < n_segregatingA; ++i)
		fprintf(stderr,  " %zu", haplotype_posns[0][i]);
	fprintf(stderr, "\nGenotype B has %u segregating sites:", n_segregatingB);
	for (unsigned int i = 0; i < n_segregatingB; ++i)
		fprintf(stderr,  " %zu", haplotype_posns[1][i]);
	fprintf(stderr, "\n");
 */
	for (unsigned int i = 0; i < N_SUBGENOMES; ++i)
		if (vcf_fp[i])
			fclose(vcf_fp[i]);


	/* post hoc coverage tests: under equal allelic coverage assumption,
	 * each haplotype should about 50% subgenomic coverage; we test
	 * for 50:50 coverage of the two dominant haplotypes to remove the
	 * degrading impact of errors that leach coverage
	 *
	 * however, the above tests do not require differences in the two
	 * dominant haplotypes at every segregating site, so some sites may
	 * have alleles with much higher than 50% coverage, if the second
	 * highest haplotype has the same allele  We also check each
	 * segregating site for 50\%, in this case ignoring leach of errors.
	 *
	 * [TODO,KSD] Move to a function, simplify logic.  Maybe check if
	 * most and second-most abundant are 50:50 per site.
	 */
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

			//fprintf(stderr, "Processing read %zu", n_read);
			/* find subgenomic assignment or continue if ambiguously aligned */
			unsigned int sgenome = 0;
			for (unsigned int j = 1; j < N_SUBGENOMES; ++j)
				if (pp[j][n_read] >= opt.phc_min_alignment_pp) {
					sgenome = j;
					break;
				}
			if (!sgenome && pp[sgenome][n_read]
						< opt.phc_min_alignment_pp) {
				++n_read;
				continue;
			}
			//fprintf(stderr, " assigned to subgenome %s.", sgenome?"B":"A");

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
			//fprintf(stderr, "Found haplotype %lu (%u)\n", id, haplotype_cnt[sgenome][hpos]);
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


EXIT_AFTER_GENOTYPING_STARTS:

	if (obs_nuc)
		free(obs_nuc);
	if (obs_q)
		free(obs_q);
	if (obs_rpos)
		free(obs_rpos);
	if (genome_src)
		free(genome_src);
	if (covers)
		free(covers);
	if (rd_idxA)
		free(rd_idxA);
	if (haplotype_posns[0])
		free(haplotype_posns[0]);
	if (haplotype_posns[1])
		free(haplotype_posns[1]);
	if (hapA_prop)
		free(hapA_prop);
	if (hapB_prop)
		free(hapB_prop);
	if (hapA_covg)
		free(hapA_covg);
	if (hapB_covg)
		free(hapB_covg);
	if (hapA_dom_nuc)
		free(hapA_dom_nuc);
	if (hapB_dom_nuc)
		free(hapB_dom_nuc);

EXIT_NOW:
	for (unsigned int j = 0; j < N_SUBGENOMES; ++j)
		if (pp[j])
			free(pp[j]);
	if (pp)
		free(pp);
	if (opt.amplici_command)
		free_input(in);
	
	return EXIT_SUCCESS;
} /* main */


/**
 * Test equal coverage of homologous chromosomes.
 *
 * @param mh		merged hash of aligned reads
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
int test_equal_homolog_coverage(merge_hash *mh, double **all,
	char_t ref_base[N_SUBGENOMES], char *covers, xy_t *obs_nuc,
	qual_t *obs_q, unsigned int *obs_rpos, int g_max[N_SUBGENOMES], 
	xy_t nuc1, xy_t nuc2, int debug_level, double pvals[N_SUBGENOMES])
{
	mlogit_stuff mls = {NULL, 0};
	double epsilon = 1e-6;
	double gamma[N_SUBGENOMES] = {g_max[0]/2., g_max[1]/2.};
	double eta = 0.5;	/* assumes N_SUBGENOMES == 2 */
	double new_gamma[N_SUBGENOMES];
	double new_eta;
	double lpi[2*N_SUBGENOMES];
	double ll1;
	double ll = -INFINITY, pll;
	double pp[2*N_SUBGENOMES], cll[2*N_SUBGENOMES];
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
		new_gamma[1] = 0;	/* N_SUBGENOMES == 2 */
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

			for (int i = 0; i < 2*N_SUBGENOMES; ++i) {

				cll[i] = lpi[i] + all[i/2][n_read]
					- sub_prob_given_q_with_encoding(
						ref_base[i/2],				/* homozygous ref base */
						obs_nuc[n_cover], IUPAC_ENCODING, XY_ENCODING, obs_q[n_cover], 1, (void *)&mls)
					+ sub_prob_given_q_with_encoding(
						 !g_max[i/2] ? xy_to_iupac[nuc1]	/* new genotype */
						: g_max[i/2] == 2 ? xy_to_iupac[nuc2]
						: !(i%2) ? xy_to_iupac[nuc1] : xy_to_iupac[nuc2],
						obs_nuc[n_cover], IUPAC_ENCODING, XY_ENCODING, obs_q[n_cover], 1, (void *)&mls);
	
				if (cll[i] > max)
					max = cll[i];
			}

			double sum = 0;
			for (int i = 0; i < 2*N_SUBGENOMES; ++i) {
				pp[i] = exp(cll[i] - max);
				sum += pp[i];
			}
			for (int i = 0; i < 2*N_SUBGENOMES; ++i)
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
		//fprintf(stderr, "eta = %f, gamma1 = %f, gamma2 = %f, ll = %f, rel = %e\n", eta, gamma[0], gamma[1], ll, (pll - ll) / ll);
	} while (iter++ < max_iter && (ll - pll) > -ll * epsilon);

	ll1 = ll;

	for (int j = 0; j < N_SUBGENOMES; ++j) {
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

			for (int i = 0; i < 2*N_SUBGENOMES; ++i) {

				cll[i] = lpi[i] + all[i/2][n_read]
					- sub_prob_given_q_with_encoding(
						ref_base[i/2],				/* homozygous ref base */
						obs_nuc[n_cover], IUPAC_ENCODING, XY_ENCODING, obs_q[n_cover], 1, (void *)&mls)
					+ sub_prob_given_q_with_encoding(
						 !g_max[i/2] ? xy_to_iupac[nuc1]	/* new genotype */
						: g_max[i/2] == 2 ? xy_to_iupac[nuc2]
						: !(i%2) ? xy_to_iupac[nuc1] : xy_to_iupac[nuc2],
						obs_nuc[n_cover], IUPAC_ENCODING, XY_ENCODING, obs_q[n_cover], 1, (void *)&mls);
	
				if (cll[i] > max)
					max = cll[i];
			}

			double sum = 0;
			for (int i = 0; i < 2*N_SUBGENOMES; ++i) {
				pp[i] = exp(cll[i] - max);
				sum += pp[i];
			}
			for (int i = 0; i < 2*N_SUBGENOMES; ++i)
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
		//fprintf(stderr, "eta = %f, gamma1 = %f, gamma2 = %f, ll = %f, rel = %e\n", eta, (g_max[0] && j) ? gamma[0] : 0.5, (g_max[1] && !j) ? gamma[1] : 0.5, ll, (pll - ll) / ll);
	} while (iter++ < max_iter && (ll - pll) > -ll * epsilon);

	double lrt = 2 * (ll1 - ll);
	pvals[j] = pchisq(lrt, 1, 0, 0);
	debug_msg(debug_level > QUIET, debug_level,
		"Equal coverage test: eta = %f; gamma1 = %f; gamma2 = %f; lrt = %f (%f %f); pval = %e\n",
		eta, j && g_max[0] == 1 ? gamma[0] : g_max[0] / 2., !j && g_max[1] == 1 ? gamma[1] : g_max[1] / 2., 2 * (ll1 - ll), ll1, ll, pvals[j]);
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
int update_vcf(options *opt, int *fail[N_SUBGENOMES],
				size_t *hpos[N_SUBGENOMES])
{
	const char *tmpfile_template = "tmp_vcfXXXXXX";

	for (int i = 0; i < N_SUBGENOMES; ++i) {
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
 * Exclude additional reads from the merge hash that have been assigned to 
 * cluster id.  Actually nothing is done, except print out the read index
 * and count the number of newly excluded reads.
 *
 * @param mh	merge hash pointer
 * @param in	input object pointer
 * @param id	id of cluster whose reads are to be removed
 * @return	number of newly excluded reads
 */
unsigned int exclude_cluster_reads(merge_hash *mh, input *in, unsigned int id)
{
	unsigned int n_read = 0;
	unsigned int n_excluded = 0;

	for (merge_hash *me = mh; me != NULL; me = me->hh.next) {
		if (me->exclude)
			continue;
		if (in->assignment[n_read] == id) {
			fprintf(stderr, " %d", n_read);
			++n_excluded;
		}
		++n_read;
	}

	return n_excluded;
} /* exclude_cluster_reads */

/**
 * Hash read on segregating sites.
 *
 * @param se		entry for current read
 * @param haplotype	segregating positions
 * @param n_segregating	number of segregating sites
 * @return		hashed haplotype
 */
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
 * @param se	alignment entry from sam file (xy_t)
 * @param rd_id	index of read
 * @param ref	reference sequence (iupac_t)
 * @param vptr	mlogit stuff
 * @param show	show alignments
 * @return	log likelihood
 */
double ll_align(sam_entry *se, unsigned int rd_id, unsigned char *ref,
		mlogit_stuff *mls, unsigned char *in_show)
{

	size_t rf_index = se->pos - 1;		/* starting reference position */
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

	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_III;//DEBUG_II;//
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
	if (fxn_debug >= DEBUG_II || show)
		for (unsigned int i = 0; i < se->cig->n_ashes; ++i)
			fprintf(stderr, "%u%c", se->cig->ashes[i].len,
				cigar_char[se->cig->ashes[i].type]);
	debug_msg_cont(fxn_debug >= DEBUG_II || show, fxn_debug, "\n");
	debug_msg(fxn_debug >= DEBUG_II || show, fxn_debug, "Ref : ");

	unsigned int out = 0;
	for (unsigned int i = 0; i < se->cig->n_ashes; ++i) {

		/* march through alignment; optimally output reference */
		if (se->cig->ashes[i].type == CIGAR_DELETION) {
			if (fxn_debug >= DEBUG_II || show)
				for (size_t j = 0; j < se->cig->ashes[i].len; ++j) {
					if (!reverse_complement && !display_after
					    && (fxn_debug || show))
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
			if (fxn_debug >= DEBUG_II || show)
				for (size_t j = 0; j < se->cig->ashes[i].len; ++j) {
					if (!reverse_complement && !display_after)
						fputc('.', stderr);
					if (!reverse_complement && !display_after && guide_posts && !(++out % 10))
						fprintf(stderr, "|");
				}
			rd_index += se->cig->ashes[i].len;
			continue;
		} else if (se->cig->ashes[i].type == CIGAR_INSERTION) {
			if (fxn_debug >= DEBUG_II || show)
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
				fprintf(stderr, "|");
		}
		rf_index += se->cig->ashes[i].len;
		rd_index += se->cig->ashes[i].len;
	}

	if ((fxn_debug >= DEBUG_II || show) && reverse_complement) {
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
	if (display_after && !reverse_complement) {
		for (size_t j = 0; j < align_len; ++j) {
			fputc(!align_display[j] ? '-' :
			      iupac_to_char[align_display[j]], stderr);
			if (!((j + 1) % 10) && guide_posts)
				fputc('|', stderr);
		}
		fprintf(stderr, "\n");
	}
	if (qual_display && !reverse_complement) {
		debug_msg(fxn_debug >= DEBUG_II || show, fxn_debug, "Qual: ");
		for (size_t j = 0; j < align_len; ++j) {
			fputc(!qual_display[j] ? ' ' : qual_display[j], stderr);
			if (!((j + 1) % 10) && guide_posts)
				fputc('|', stderr);
		}
		fprintf(stderr, "\n");
	}

	/* optionally output read to show alignment */
	if (fxn_debug >= DEBUG_II || show) {
		align_index = 0;
		debug_msg(fxn_debug >= DEBUG_II || show, fxn_debug, "Read: ");
		rf_index = se->pos - 1;
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
		if (reverse_complement && (fxn_debug || show)) {
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
 * Count the number (of the top 4) haplotypes to discard.
 *
 * @param opt	options pointer
 * @param input	input object
 * @return	error status
 */
int discard_haplotypes(options *opt, input *input)
{
	int err = NO_ERROR;
	FILE *fp = NULL;
	unsigned int abund[2*N_SUBGENOMES] = {0, 0, 0, 0};	/* estimated abundance of haplotypes */
	double ll1, xmin, ll0 = 0, lrt, pval;
	unsigned int cmd_len = 0;
	char *command = NULL;

	if (!opt->amplici_complete) {
		cmd_len = strlen(opt->amplici_command)
			+ strlen(opt->ac_fastq_file)
			+ strlen(" -f  -lb 10.5 -m 4") + 1;	/* BUG: assumes \par opt.ac_low_bound < 100 */
		command = malloc(cmd_len * sizeof *command);

		sprintf(command, "%s -f %s -lb %.1f -m 4",
			opt->amplici_command, opt->ac_fastq_file,
						opt->ac_low_bound);
		mmessage(INFO_MSG, NO_ERROR, "Running amplici for top four "
				"most observed haplotypes: '%s'\n", command);
	
		fp = popen(command, "r");
		if (!fp)
			return mmessage(ERROR_MSG, PIPE_OPEN_ERROR, command);

		for (int j = 0; j < 2*N_SUBGENOMES; ++j)
			if (fscanf(fp, "%u", &abund[j]) != 1)
				return mmessage(ERROR_MSG, PIPE_READ_ERROR,
					"failed to read integer %d from command"
							"'%s'\n", j, command);
		free(command);
	}
			
	if (opt->amplici_complete || abund[0] < opt->proptest_min_freq) {
		cmd_len = strlen(opt->amplici_command)
			+ strlen(opt->ac_fastq_file) + strlen(opt->ac_outfile)
				+ strlen(" -f  -lb  -ll  -o ") + 1
				+ (int)(log10(opt->ac_low_bound) + 3)
				+ (int)(isfinite(opt->ac_ll_threshold) ? 2
					: log10(fabs(opt->ac_ll_threshold))) + 2;
		command = malloc(cmd_len * sizeof *command);
		sprintf(command, "%s -f %s -lb %.1f -ll %.0f -o %s",		/* BUG: assumes no scatter reads */
			opt->amplici_command, opt->ac_fastq_file,
			opt->ac_low_bound, opt->ac_ll_threshold,
							opt->ac_outfile);

		if (!opt->amplici_complete) {
			mmessage(INFO_MSG, NO_ERROR, "Most abundant sequence is"
				" observed only %u times.\n", abund[0]);
			mmessage(INFO_MSG, NO_ERROR, "Rerunning amplici for "
				"top four most populated clusters: '%s'\n",
								command);
		} else {
			mmessage(INFO_MSG, NO_ERROR, "Running amplici: '%s'\n",
								command);
		}

		system(command);
		free(command);

		cmd_len = strlen(opt->ac_outfile) + strlen(".out") + 1;
		command = malloc(cmd_len * sizeof *command);
		sprintf(command, "%s.out", opt->ac_outfile);

		fp = fopen(command, "r");
		if (!fp)
			return mmessage(ERROR_MSG, FILE_OPEN_ERROR, command);

		if ((err = read_amplici_results(fp, input, opt)))
			return err;

		for (unsigned int j = 0; j < MIN(2*N_SUBGENOMES, input->K); ++j)
			abund[j] = input->cluster_sizes[input->sort_id[j]];

	}

	mmessage(INFO_MSG, NO_ERROR, "Most abundant haplotypes/clusters: ");
	fprint_uints(stderr, abund, 2*N_SUBGENOMES, 3, 1);

	/* start keeping all 4 haplotypes */
	opt->proptest_screen = 0;
	
	/* entering strictly allotetraploid area */

	/* H40: 4 haplotypes with decreasing abundance a, a, b, b */
	if (input->K > 3) {
		input->reject_state = ACCEPT_FOUR;
		ll1 = max_multinom_ll(4, abund);
		input->mle_w = brent(0.01, 100, four_haplotype, (void *) abund);
		ll0 = - four_haplotype(input->mle_w, (void *) abund);
		lrt = -2 * (ll0 - ll1);
		pval = pchisq(lrt, 2, 0, 0);
		mmessage(INFO_MSG, NO_ERROR, "H04 (p1=p2, p3=p4): double SNP "
			"p-value %e (lrt=%f, %f %f) %s\n", pval, lrt,
			0.5 * input->mle_w / (1 + input->mle_w),
			0.5 / (1 + input->mle_w),
			pval < opt->proptest_alpha ? "REJECT" : "ACCEPT");
	} else {
		pval = 0;
	}
	if (pval >= opt->proptest_alpha) {
		ll1 = ll0;
		ll0 = - four_haplotype(1, (void *) abund);
		lrt = -2 * (ll0 - ll1);
		input->pval_equal_subgenomic_coverage = pchisq(lrt, 1, 0, 0);
		mmessage(INFO_MSG, NO_ERROR, "H04a (p1=p2=p3=p4): double SNP, "
			"equal subgenomic coverage p-value %e (lrt=%f) %s\n",
				input->pval_equal_subgenomic_coverage, lrt,
			input->pval_equal_subgenomic_coverage
				< opt->proptest_alpha ? "REJECT" : "ACCEPT");

	/* reject H40: there are two tests for 3 haplotypes */
	} else if (pval < opt->proptest_alpha) {
		input->reject_state = ACCEPT_THREE_A;
		++opt->proptest_screen;	/* remove 1 of 4 so far */

		/* H30a: 3 haplotypes of decreasing abundance a, b, b */
		if (input->K > 2) {
			ll1 = max_multinom_ll(3, abund);
			xmin = brent(0.01, 100, three_haplotype, (void *) abund);
			ll0 = - three_haplotype(xmin, (void *) abund);
			lrt = -2 * (ll0 - ll1);
			pval = pchisq(lrt, 1, 0, 0);
			mmessage(INFO_MSG, NO_ERROR, "H03a (p1, p2=p3): single "
				"SNP p-value %e (lrt=%f, %f %f) %s\n", pval,
				lrt, xmin / (1 + xmin), 0.5/(1+xmin),
				pval < opt->proptest_alpha ? "REJECT"
								: "ACCEPT");
		}
		if (pval < opt->proptest_alpha) {
			input->reject_state = ACCEPT_THREE_B;
			/* H30b: 3 haplotypes of decreasing abundance a, a, b */
			/* swap */
			if (input->K > 2) {
				ll0 = abund[0];
				abund[0] = abund[2];
				abund[2] = ll0;
				xmin = brent(0.01, 100, three_haplotype,
								(void *) abund);
				ll0 = - three_haplotype(xmin, (void *) abund);
				lrt = -2 * (ll0 - ll1);
				pval = pchisq(lrt, 1, 0, 0);
				mmessage(INFO_MSG, NO_ERROR, "H03b (p1=p2, p3):"
					"single SNP p-value %e (lrt=%f, %f %f) "
					"%s\n", pval, lrt, 1/(1 + xmin),
					0.5*xmin/(1+xmin), pval
					< opt->proptest_alpha ? "REJECT"
								: "ACCEPT");
			}
			if (pval < opt->proptest_alpha) {
				input->reject_state = ACCEPT_TWO;
				++opt->proptest_screen;
			} else {
				input->mle_w = xmin;
			}
		} else {
			input->mle_w = xmin;
		}
	}

	input->proptest = opt->proptest_screen;

	return err;
} /* discard_haplotypes */


/**
 * Four haplotype model where there are true heterozygous SNPs in both
 * subgenomes, excluding multinomial coefficient.
 *
 * @param w	sampling bias (optimizing over this)
 * @param param	void pointer to parameters
 * @return	negative log likelihood without multinomial coefficient
 */
double four_haplotype(double w, void *param)
{
	unsigned int *abun = (unsigned int *) param;
	unsigned sum = 0;
	double ll = 0;

	for (int j = 0; j < 4; ++j) {
		sum += abun[j];
		if (j < 2)
			ll += abun[j] * log(w);
	}
	ll += sum * log(0.5 / (1 + w));

//	fprintf(stderr, "four_haplotype(%f: %f, %f); %f\n", w, 0.5 * w / (1 + w), 0.5 / (1 + w), ll);

	return -ll;
}/* four_haplotype */

/**
 * Three haplotype model where there are true allelic SNPs in one subgenome.
 *
 * @param w	parameter to maximize, sampling bias
 * @param param	other parameters
 * @return	negative log likelihood without multinomial coefficient
 */
double three_haplotype(double w, void *param)
{
	unsigned int *abun = (unsigned int *) param;
	unsigned int sum = 0;
	double ll = 0;

	for (int j = 0; j < 3; ++j) {
		sum += abun[j];
		if (j < 1)
			ll += abun[j] * log(w);
		else
			ll += abun[j] * log(0.5);
	}
	ll -= sum * log(1 + w);

//	fprintf(stderr, "three_haplotype(%f: %u*%f, (%u,%u)*%f): %f\n", w, abun[0], w/(1+w), abun[1], abun[2], 0.5/(1+w), ll);

	return -ll;
} /* three_haplotype */

/**
 * Maximum log likelihood for multinomial distribution.
 *
 * @param p	number of categories
 * @param cnts	counts for each category
 * @return	maximum log likelihood
 */
double max_multinom_ll(unsigned int p, unsigned int *cnts)
{
	unsigned int sum = 0;	/* 1, 2, 3, ..., n */
	double ll = 0;

	for (unsigned int i = 0; i < p; ++i) {
		ll += cnts[i] * log(cnts[i]);
		sum += cnts[i];
/*		for (unsigned int j = 1; j <= cnts[i]; ++j)
			ll += log(sum++) - log(j);
 */
	}
	ll -= sum * log(sum);

	return ll;
} /* max_multinom_ll */

/**
 * Read the output file produced by amplici.  The input fastq file to ampliCI
 * excludes all reads already screened out by roshan.c, including those marked
 * as excluded in the merge hash.  This code is recording indices in the merge
 * hash that should be excluded because of closeness to low abundance
 * haplotypes.
 *
 * This definition and that of make_input() and free_input() may belong
 * in a different file, but they do not belong in sam.c, which is focused
 * on parsing and storing data from sam files.
 *
 * @param fp	open handle of amplici output file
 * @param in	store data in this structure
 * @param opt	options object pointer
 * @return error status
 */
int read_amplici_results(FILE *fp, input *in, options *opt)
{
	char c;
	unsigned int i;

	MAKE_1ARRAY(in->assignment, in->n_observation);

	/* get read assignments (indices of assigned haplotype) and
	 * haplotype estimated true abundances
	 */
	while (!feof(fp)) {
		c = fgetc(fp);
		if (c == 'K') {		/* number of clusters K */
			fgetc(fp);	/* discard ':' */
			fscanf(fp, "%u", &in->K);
			in->cluster_sizes = malloc(in->K
						* sizeof *in->cluster_sizes);
			if (!in->cluster_sizes)
				return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"input:cluster_sizes");
			in->excluded_clusters = calloc(in->K,
						sizeof *in->excluded_clusters);
			if (!in->excluded_clusters)
				return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"input:excluded_clusters");
			in->sort_id = malloc(in->K * sizeof *in->sort_id);
			if (!in->sort_id)
				return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"input:sort_id");
			for (unsigned int k = 0; k < in->K; ++k)
				in->sort_id[k] = k;
			if (opt->genotype_by_clustering) {
				in->cluster_to_csome = calloc(in->K,
						sizeof *in->cluster_to_csome);
				if (!in->cluster_to_csome)
					return mmessage(ERROR_MSG,
						MEMORY_ALLOCATION,
						"input::cluster_to_csome");
			}

		} else if (c == 'a') {
			c = fgetc(fp);
			if (c == 's') {	/* assignments! */
				while ((c = fgetc(fp)) != ':');
				for (i = 0; i < in->n_observation; ++i)
					if (fscanf(fp, "%u", &in->assignment[i])
									!= 1) {
						while ((c = fgetc(fp)) == ' ');
						if (c == 'N') {	/* NA */
							in->assignment[i]
								= UINT_MAX;
							c = fgetc(fp);	/* A */
						} else {
							break;
						}
					}
				if (i < in->n_observation)
					return mmessage(ERROR_MSG,
						FILE_FORMAT_ERROR, "Invalid "
						"assignments line in amplici "
								"results.\n");
			}
		} else if (c == 's') {	/* sizes */
			c = fgetc(fp);
			if (c != 'i')
				continue;
			while ((c = fgetc(fp)) != ':');
			for (i = 0; i < in->K; ++i)
				fscanf(fp, "%u", &in->cluster_sizes[i]);	/* BUG, TODO, KSD: NO CHECK */
			break;	/* assumes particular order */
		}
	}

	fclose(fp);

	/* sort haplotype clusters by size */
	reverse_order_uint_quicksort_with_index(in->cluster_sizes,
						in->sort_id, in->K);
	for (unsigned int k = 0; k < in->K; ++k)
		fprintf(stderr, "Cluster_id %u of size %u goes to sort_id %zu\n", k, in->cluster_sizes[k], in->sort_id[k]);

	return NO_ERROR;

} /* read_amplici_results */


/********** genotype by clustering **********/


/**
 * Genotype by clustering.
 *
 * After finding the read-supported haplotypes (up to four), we align the
 * haplotypes to the corresponding subgenomic reference and report up to two
 * alternate alleles for the up to two haplotypes assigned to each subgenome.
 * To assess confidence, we take the nucleotides from reads assigned to each
 * subgenome (via their assignment to the 2 haplotypes assigned to that
 * subgenome) and compute the likelihood of these reads given the three
 * possible genotypes.  If there is only one haplotype assigned to the
 * subgenome, we look for the most supported alternative allele in the reads
 * assigned to the subgenome (via the single haplotype assigned to that
 * subgenome).  If there is no alternative allele, we ...
 *
 * @param opt	options object
 * @param in	input object
 * @param mh	merge hash
 * @param sds	sam objects
 * @param rseqs	reference sequences
 * @return	error status
 */
int genotype_by_clustering(options *opt, input *in, merge_hash *mh,
	sam *sds[N_SUBGENOMES], char_t const * const rseqs[N_SUBGENOMES],
					unsigned int const rlens[N_SUBGENOMES])
{
	int err = NO_ERROR;

	/* entering allotetraploid specific portion */
	if (N_SUBGENOMES != 2)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "Genotyping by "
			"clustering now only written for allotetraploids.\n");

	/* open vcf files for writing and output header */
	FILE *fp[N_SUBGENOMES] = {NULL, NULL};
	for (int i = 0; i < N_SUBGENOMES; ++i) {

		if (!opt->vcf_files[i])
			continue;
		fp[i] = fopen(opt->vcf_files[i], "w");

		if (!fp[i]) {
			err = mmessage(ERROR_MSG, FILE_OPEN_ERROR,
							opt->vcf_files[i]);
			goto GENOTYPE_BY_CLUSTERING_RETURN;
		}

		print_vcf_header(fp[i], opt->vcf_opt, opt->fsa_files[i],
							opt->sample_name);
	}

	/* align selected haplotypes and load one for each chromosome */
	char_t const *haplotypes[2*N_SUBGENOMES] = {NULL, NULL, NULL, NULL};
	fastq_data *fd = NULL;

	if ((err = get_aligned_haplotypes(in, opt, &fd, haplotypes)))
		goto GENOTYPE_BY_CLUSTERING_RETURN;

	/* create profile for each subgenome: ambiguous nucleotide
	 * representing pair of nucleotides observed at a locus
	 */
	char_t *subgenomic_profiles[N_SUBGENOMES] = {NULL, NULL};

	for (int i = 0; i < N_SUBGENOMES; ++i) {
		subgenomic_profiles[i] = malloc(fd->n_max_length
					* sizeof *subgenomic_profiles[i]);
		if (!subgenomic_profiles[i]) {
			err = mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"subgenomic_profile");
			goto GENOTYPE_BY_CLUSTERING_RETURN;
		}
		for (unsigned int j = 0; j < fd->n_max_length; ++j)
			subgenomic_profiles[i][j]
				= haplotypes[i*2][j] | haplotypes[i*2 + 1][j];
	}

	/* align subgenomic reference to haplotype profiles */
	if ((err = get_reference_aligned_haplotypes(in, fd, (const char_t **)
			subgenomic_profiles, (const char_t **) rseqs, rlens)))
		goto GENOTYPE_BY_CLUSTERING_RETURN;

	/* subgenomic coverage estimated as size of clusters assigned to 
	 * subgenome
	 */
	unsigned int subgenomic_coverage[N_SUBGENOMES] = {0, 0};
	subgenomic_coverage[0]
		= in->cluster_sizes[in->included_clusters[in->csome_to_inc[0]]]
			+ (in->csome_to_inc[0] != in->csome_to_inc[1]
			? in->cluster_sizes[in->included_clusters[
						in->csome_to_inc[1]]] : 0);
	subgenomic_coverage[1]
		= in->cluster_sizes[in->included_clusters[in->csome_to_inc[2]]]
			+ (in->csome_to_inc[2] != in->csome_to_inc[3]
			? in->cluster_sizes[in->included_clusters[
						in->csome_to_inc[3]]] : 0);

		/* number of gaps in each subgenomic alignment */
	unsigned int n_gaps[2*N_SUBGENOMES] = {0,0,0,0};
		/* current position within subgenomic profile */
	unsigned int pos_ref[N_SUBGENOMES] = {0,0};
		/* current position in profile alignment to reference */
	unsigned int ref_prof_align_pos[N_SUBGENOMES] = {0,0};
		/* current site is indel variant: don't genotype */
	int has_gap;

	/* genotype each position of haplotype msa */
	for (unsigned int j = 0; j < fd->n_max_length; ++j) {
		unsigned int n_read;
		unsigned int no_alternate_allele = 0;
		char_t nucA1 = haplotypes[0][j];
		char_t nucA2 = haplotypes[1][j];
		char_t nucB1 = haplotypes[2][j];
		char_t nucB2 = haplotypes[3][j];
		char_t hap_nucs[2][2] = {{nucA1, nucA2}, {nucB1, nucB2}};
		char_t hap_nuc_array[4] = {nucA1, nucA2, nucB1, nucB2};
		char_t nucA2b = nucA2, nucB2b = nucB2;
		char_t nuc1 = nucA1, nuc2 = nucB2;
		char_t alt_nucs[2];
		int ref_is_gap[2] = {0,0};	/* insertion in haplotypes */

		/* keep track of position in each reference and
		 * reference/profile alignment
		 */
		for (int i = 0; i < N_SUBGENOMES; ++i) {

			/* gaps in the profile; deleted reference */
			while (in->ref_profile_alignments[i][1][
					ref_prof_align_pos[i]] == '-') {

				/* nothing to genotype, but output position */
				if (fp[i]) {
					fprintf(fp[i], "%s\t%u\t.\t%c\t*\t",
						opt->ref_names[i],
						ref_is_gap[i] && pos_ref[i] ? pos_ref[i] : pos_ref[i] + 1,
						iupac_to_char[rseqs[i][pos_ref[i]]]);
					fprintf(fp[i], ".\t.\t.\t.\t.\n");
				}
				++pos_ref[i];	/* advance reference position */
				++ref_prof_align_pos[i];
			}

			/* gaps in reference; profile insertion */
			if (in->ref_profile_alignments[i][0][ref_prof_align_pos[i]] == '-') {
				/* nothing to do: will output previous ref base and position */
				/* next reference base if no reference yet output */
				ref_is_gap[i] = 1;
			}
			++ref_prof_align_pos[i];
		}

#ifdef ALLOW_DEBUGGING
if (opt->debugging_site == (int) j)
fprintf(stderr, "Site (%u, %c) / (%u, %c): %c %c %c %c (%c %c %c %c)\n",
	pos_ref[0] + 1, iupac_to_char[rseqs[0][pos_ref[0]]],
	pos_ref[1] + 1, iupac_to_char[rseqs[1][pos_ref[1]]],
	iupac_to_char[nucA1], iupac_to_char[nucA2],
	iupac_to_char[nucB1], iupac_to_char[nucB2],
	iupac_to_char[haplotypes[0][j]], iupac_to_char[haplotypes[1][j]], iupac_to_char[haplotypes[2][j]], iupac_to_char[haplotypes[3][j]]);
#endif

		/* count gaps in aligned haplotypes */
		has_gap = 0;
		for (int i = 0; i < 4; ++i)
			if (hap_nuc_array[i] == IUPAC_X) {
				++n_gaps[i];
				has_gap = 1;
			}

		/* if there is a gap report genotype, but no confidence */
		if (has_gap) {
			mmessage(INFO_MSG, NO_ERROR, "Site %u has gap, not "
						"genotyping.\n", j + 1);
			for (int i = 0; i < N_SUBGENOMES; ++i) {
				int alt_output = 0;

				if (!fp[i])
					continue;
	
				fprintf(fp[i], "%s\t%u\t.\t%c\t",
						opt->ref_names[i],
						ref_is_gap[i] && pos_ref[i] ? pos_ref[i] : pos_ref[i] + 1,
						iupac_to_char[rseqs[i][
							ref_is_gap[i]
							&& pos_ref[i]
							? (pos_ref[i] - 1)
							: pos_ref[i]]]);

				/* most abundant haplotype is gap */
				if (!hap_nucs[i][0]) {
					fputc('*', fp[i]);
					alt_output = 1;

				/* reference is gap, haplotype insertion */
				} else if (ref_is_gap[i]) {
					fprintf(fp[i], "%c", iupac_to_char[
							hap_nucs[i][0]]);
					alt_output = 1;

				/* haplotype differs from reference */
				} else if (hap_nucs[i][0]
						!= rseqs[i][pos_ref[i]]) {
					fprintf(fp[i], "%c", iupac_to_char[
							hap_nucs[i][0]]);
					alt_output = 1;
				}

				/* 2nd allele is missing */
				if (hap_nucs[i][1] != hap_nucs[i][0]
							&& !hap_nucs[i][1]) {
					if (alt_output)
						fputc(',', fp[i]);
					fputc('*', fp[i]);

				/* 2nd allele differs from first */
				} else if (ref_is_gap[i] && hap_nucs[i][1]
							!= hap_nucs[i][0]) {
					if (alt_output)
						fputc(',', fp[i]);
					fprintf(fp[i], "%c",
						iupac_to_char[hap_nucs[i][1]]);
				} else if (!ref_is_gap[i] && hap_nucs[i][1]
					!= hap_nucs[i][0] && hap_nucs[i][1]
						!= rseqs[i][pos_ref[i]]) {
					if (alt_output)
						fputc(',', fp[i]);
					fprintf(fp[i], "%c",
						iupac_to_char[hap_nucs[i][1]]);

				/* no information */
				} else if (!alt_output) {
					fprintf(fp[i], ".");
				}

				fprintf(fp[i], "\t.\t%s\t.\tGT\t",
					in->n_indels[i] ? "sgi,hmi": "hmi");

				if (ref_is_gap[i]) {
					fprintf(fp[i], "1/%d\n", hap_nucs[i][0]
						== hap_nucs[i][1] ? 1 : 2);
				} else if (hap_nucs[i][0]
						== rseqs[i][pos_ref[i]]) {
					fprintf(fp[i], "0/%d\n", hap_nucs[i][1]
						== hap_nucs[i][0] ? 0 : 1);
				} else if (hap_nucs[i][1]
						== rseqs[i][pos_ref[i]]) {
					fprintf(fp[i], "0/%d\n", hap_nucs[i][1]
						== hap_nucs[i][0] ? 0 : 1);
				} else {
					fprintf(fp[i], "1/1\n");
				}
				if (!ref_is_gap[i])
					++pos_ref[i];
			}
#ifdef ALLOW_DEBUGGING
if (opt->debugging_site == (int) j)
	fprintf(stderr, "Site has gap, continuing...\n");
#endif

			continue;
		}

		/* find two possible alleles */
		if (nucA2 != nucA1)
			nuc2 = nucA2;
		else if (nucB1 != nucA1)
			nuc2 = nucB1;

		/* no alternative allele in genotype, but we need one to
		 * assess confidence
		 */
		unsigned int n_cover = 0;
		if (nuc1 == nuc2) {
			unsigned int num_nuc[NUM_NUCLEOTIDES];
			unsigned int subgenomic_cover[N_SUBGENOMES];
			unsigned int num_nuc1 = 0, num_nuc2 = 0;

			for (unsigned int i = 0; i < NUM_NUCLEOTIDES; ++i)
				num_nuc[i] = 0;
			for (unsigned int i = 0; i < N_SUBGENOMES; ++i)
				subgenomic_cover[i] = 0;

			n_read = 0;
			for (merge_hash *me = mh; me != NULL; me = me->hh.next) {

				/* skip excluded reads */
				if (me->exclude)
					continue;

				/* ... and those with NA assignment or in
				 * excluded (paralog) clusters
				 */
				if (in->assignment[n_read] == NA_ASSIGNMENT ||
					in->excluded_clusters[in->assignment[n_read]]) {
					++n_read;
					continue;
				}

				unsigned int assign = in->cluster_to_csome[
							in->assignment[n_read]];
				/* take read from subgenome A alignment */
	                        sequence *read = sds[0]->se[me->indices[0][0]].read;
				if (j >= read->len + n_gaps[assign]) {
					++n_read;
					continue;
				}
				++subgenomic_cover[assign / N_SUBGENOMES];
#ifdef ALLOW_DEBUGGING
if (opt->debugging_site == (int) j)
fprintf(stderr, "Read %s %c (%u %u)\n", sds[0]->se[me->indices[0][0]].name, xy_to_char[get_nuc(read, XY_ENCODING, j - n_gaps[assign])], subgenomic_cover[0], subgenomic_cover[1]);
#endif
				++num_nuc[get_nuc(read, XY_ENCODING, j - n_gaps[assign])];
				++n_cover;
				++n_read;
			}

			/* should not happen */
			if (!n_cover) {
				mmessage(WARNING_MSG, NO_ERROR, "No read "
					"coverage at site %u.\n", j + 1);
				break;
			}

                	/* screen subgenomic coverage relative to subgenomic
			 * coverage assessed by cluster sizes; see comments
			 * around similar test for regular genotyping
			 */
			if ((double) subgenomic_cover[0] / subgenomic_coverage[0]
				< opt->coverage_screen || (double) subgenomic_cover[1]
				/ subgenomic_coverage[1] < opt->coverage_screen) {

	                        mmessage(INFO_MSG, NO_ERROR, "Coverage rate is "
					"too low, will not genotype site!\n");
				for (int i = 0; i < N_SUBGENOMES; ++i) {
					if (!fp[i])
						continue;
					fprintf(fp[i], "%s\t%u\t.\t%c\t."
						"\t.\tc%2.0f\t.\t.\t.\n",
						opt->ref_names[i],
						ref_is_gap[i] && pos_ref[i] ? pos_ref[i] : pos_ref[i] + 1,
						iupac_to_char[rseqs[i][
							ref_is_gap[i]
							&& pos_ref[i]
							? (pos_ref[i] - 1)
							: pos_ref[i]]],
						100*opt->coverage_screen);
					if (!ref_is_gap[i])
						++pos_ref[i];
				}
				continue;
			}


#ifdef ALLOW_DEBUGGING
if (opt->debugging_site == (int) j)
	fprintf(stderr, "nuc1==nuc2: %u %u %u %u\n", num_nuc[0], num_nuc[1], num_nuc[2], num_nuc[3]);
#endif
			/* find top three most abundant nucleotides */
			for (int b = 0; b < NUM_NUCLEOTIDES; ++b) {
				if (num_nuc[b] > num_nuc1) {
					nuc2 = nuc1;
					num_nuc2 = num_nuc1;
					nuc1 = xy_to_iupac[b];
					num_nuc1 = num_nuc[b];
				} else if (num_nuc[b] > num_nuc2) {
					nuc2 = xy_to_iupac[b];
					num_nuc2 = num_nuc[b];
				}
			}

			/* if no second allele observed, pick any other
			 * as a strawman; will not report in vcf
			 */
			if (nuc2 == nuc1) {
				no_alternate_allele = 1;
				nuc2 = (nuc1 << 1) & 15L;
				if (!nuc2)
					nuc2 = IUPAC_A;
			}

			/* should not happen */
			if (nuc1 != nucA1 && nuc2 != nucA1) {
				err = mmessage(ERROR_MSG, INTERNAL_ERROR,
					"Found more abundant nucleotides %c%c "
					"than the call (%c, %c, %c, %c) at site"
					" %u!\n",
					iupac_to_char[nuc1], iupac_to_char[nuc2],
					iupac_to_char[nucA1], iupac_to_char[nucA2],
					iupac_to_char[nucB1], iupac_to_char[nucB1],
									j + 1);
				goto GENOTYPE_BY_CLUSTERING_RETURN;
			}

			if (nuc2 == nucA1) {
				data_t tmp = nuc2;
				nuc2 = nuc1;
				nuc1 = tmp;
			}

			/* could check for biallelic SNP, but since this 
			 * code only triggered for non-SNPs, don't do it
			 */

		} else {
			n_cover = 1;
		}

		/* nothing to genotype with; warning already issued */
		if (!n_cover)
			break;

		/* nucA2b and nucB2b are stand-in alleles for genotyping;
		 * they are not cluster-supported if haplotypes are
		 * homozygous at this site
		 */
		if (nucA1 == nucA2)
			nucA2b = nuc2;
		else
			nucA2b = nucA2;
		if (nucB1 == nucB2) {
			if (nucB1 == nuc1)
				nucB2b = nuc2;
			else
				nucB2b = nuc1;
		} else {
			nucB2b = nucB2;
		}
		alt_nucs[0] = nucA2b;
		alt_nucs[1] = nucB2b;

#ifdef ALLOW_DEBUGGING
if (opt->debugging_site == (int) j)
fprintf(stderr, "Site (%u, %c) / (%u, %c): %c %c (%c) %c %c (%c)\n",
	pos_ref[0] + 1, iupac_to_char[rseqs[0][pos_ref[0]]],
	pos_ref[1] + 1, iupac_to_char[rseqs[1][pos_ref[1]]],
	iupac_to_char[nucA1], iupac_to_char[nucA2], iupac_to_char[nucA2b],
	iupac_to_char[nucB1], iupac_to_char[nucB2], iupac_to_char[nucB2b]);
#endif

		/* compute the log likelihood & posterior probability of each 
		 * of 3 possible genotypes in each subgenome
		 */
		int max_g[2] = {0, 0};	/* max. post. gen. */
		double lprob[2][3] = {{0,0,0},{0,0,0}};
		double gprob[2][3] = {{0,0,0}, {0,0,0}};
		double max[2] = {-INFINITY, -INFINITY};	/* max ll */
		double den[2] = {0, 0};	/* normalize */

		for (int g = 0; g <= 2; ++g) {

			/* compute likelihood of this genotype,
			 * assuming read cluster assignments determine
			 * source nucleotide from haplotype
			 */
			n_read = 0;
			for (merge_hash *me = mh; me != NULL; me = me->hh.next) {

				/* skip excluded reads and those without
				 * assignments or assigned to paralogs
				 */
				if (me->exclude)
					continue;

				if (in->assignment[n_read] == NA_ASSIGNMENT
						|| in->excluded_clusters[
						in->assignment[n_read]]) {
					++n_read;
					continue;
				}

				/* nuc & qual from 1st A alignment */
				unsigned int assign = in->cluster_to_csome[
							in->assignment[n_read]];
				sam_entry *se = &sds[0]->se[me->indices[0][0]];
				data_t read_nuc = get_nuc(se->read, XY_ENCODING,
							j - n_gaps[assign]);
				data_t qual = get_qual(se->qual,
							j - n_gaps[assign]);

				/* subgenomic assignment from cluster */
				data_t ref_nuc = g == 1 ? (nucA1 | nucA2b)
					: !g ? nucA1 : nucA2b;
				int a_source = 1;

				if (in->assignment[n_read]
						== in->included_clusters[
							in->csome_to_inc[2]]
					|| in->assignment[n_read]
						== in->included_clusters[
							in->csome_to_inc[3]]) {
					ref_nuc = g == 1 ? (nucB1 | nucB2b)
						: !g ? nucB1 : nucB2b;
					a_source = 0;
				}

				double tmp = sub_prob_given_q_with_encoding(
					ref_nuc, read_nuc,
					IUPAC_ENCODING, XY_ENCODING,
					     qual, 1, (void *) NULL);
#ifdef ALLOW_DEBUGGING
if (opt->debugging_site == (int) j)
fprintf(stderr, "Site %u of subgenome %c, %c (g=%d): read %u, %c %d\n", 
	a_source ? pos_ref[0] : pos_ref[1],
	a_source ? 'A' : 'B',
	iupac_to_char[ref_nuc], g, n_read, xy_to_char[read_nuc], qual);
#endif

				if (a_source)
					lprob[0][g] += tmp;
				else 
					lprob[1][g] += tmp;
				++n_read;
			}
			if (max[0] < lprob[0][g]) {
				max[0] = lprob[0][g];
				max_g[0] = g;
			}
			if (max[1] < lprob[1][g]) {
				max[1] = lprob[1][g];
				max_g[1] = g;
			}
		}

		/* normalize */
		for (int g = 0; g <= 2; ++g) {
			gprob[0][g] = exp(lprob[0][g] - max[0]);
			den[0] += gprob[0][g];
			gprob[1][g] = exp(lprob[1][g] - max[1]);
			den[1] += gprob[1][g];
		}
		for (int g = 0; g <= 2; ++g) {
			gprob[0][g] /= den[0];
			gprob[1][g] /= den[1];
		}

		/* vcf output */
		for (int i = 0; i < N_SUBGENOMES; ++i) {
			int alt_output = 0;

			if (!fp[i])
				continue;

			fprintf(fp[i], "%s\t%u\t.\t%c\t",
					opt->ref_names[i], 
					ref_is_gap[i] && pos_ref[i] ? pos_ref[i] : pos_ref[i] + 1,
					iupac_to_char[rseqs[i][
						ref_is_gap[i]
						&& pos_ref[i]
						? (pos_ref[i] - 1)
						: pos_ref[i]]]);

			/* reference is gap; haplotype 1 too */
			if (ref_is_gap[i] && !hap_nucs[i][0]) {
				fputc('*', fp[i]);
				alt_output = 1;
			} else if (ref_is_gap[i]) {
				fprintf(fp[i], "%c",
						iupac_to_char[hap_nucs[i][0]]);
				alt_output = 1;
			} else if (!ref_is_gap[i] && hap_nucs[i][0]
						!= rseqs[i][pos_ref[i]]) {
				fprintf(fp[i], "%c",
						iupac_to_char[hap_nucs[i][0]]);
				alt_output = 1;
			}

			/* ref is gap, hap 1 is not gap, hap 2 is gap */
			if (ref_is_gap[i] && !hap_nucs[i][1]
							&& hap_nucs[i][0]) {
				if (alt_output)
					fputc(',', fp[i]);
				fputc('*', fp[i]);
			/* ref is gap, hap 2 is unseen allele */
			} else if (ref_is_gap[i] && hap_nucs[i][1]
					&& hap_nucs[i][1] != hap_nucs[i][0]) {
				if (alt_output)
					fputc(',', fp[i]);
				fprintf(fp[i], "%c",
						iupac_to_char[hap_nucs[i][1]]);
			/* ref is not gap, hap 2 is only gap */
			} else if (!ref_is_gap[i] && !hap_nucs[i][1]
				&& hap_nucs[i][1] != hap_nucs[i][0]) {
				if (alt_output)
					fputc(',', fp[i]);
				fputc('*', fp[i]);
			/* ref is not gap, hap 2 is unseen allele */
			} else if (!ref_is_gap[i]
					&& hap_nucs[i][1] != hap_nucs[i][0]
					&& hap_nucs[i][1] != rseqs[i][pos_ref[i]]) {
				if (alt_output)
					fputc(',', fp[i]);
				fprintf(fp[i], "%c",
						iupac_to_char[hap_nucs[i][1]]);
			/* ref is gap, hap1 == hap2, but another allele observed */
			} else if (ref_is_gap[i] && !no_alternate_allele) {
				if (alt_output)
					fputc(',', fp[i]);
				fprintf(fp[i], "%c",
						iupac_to_char[alt_nucs[i]]);
			/* ref is not gap, hap1 == hap2, another allele observed */
			} else if (!no_alternate_allele
				&& alt_nucs[i] != rseqs[i][pos_ref[i]]) {
				if (alt_output)
					fputc(',', fp[i]);
				fprintf(fp[i], "%c",
						iupac_to_char[alt_nucs[i]]);
			/* no variant */
			} else if (!alt_output) {
				fprintf(fp[i], ".");
			}
			fprintf(fp[i], "\t.\t%s\t.\tGT:DP:GQ",
				in->n_indels[i] ? "sgi" : "PASS");

			if (opt->vcf_opt->output_gl)
				fprintf(fp[i], ":GL");
			fputc('\t', fp[i]);

			if (!ref_is_gap[i] && hap_nucs[i][0] == rseqs[i][pos_ref[i]])
				fprintf(fp[i], "0/%d",
					hap_nucs[i][1] == hap_nucs[i][0] ? 0 : 1);
			else if (!ref_is_gap[i] && hap_nucs[i][1]
							== rseqs[i][pos_ref[i]])
				fprintf(fp[i], "0/%d",
					hap_nucs[i][0] == hap_nucs[i][1] ? 0 : 1);
			else
				fprintf(fp[i], "1/%d",
					hap_nucs[i][1] == hap_nucs[i][0] ? 1 : 2);
			fprintf(fp[i], ":%d:%d", subgenomic_coverage[i],
				gprob[i][max_g[i]] < 1 ? MIN(99, (int)(-10
				* log10(1 - gprob[i][max_g[i]]))) : 99);
			if (opt->vcf_opt->output_gl)
				fprintf(fp[i], ":%.2f,%.2f,%.2f",
					lprob[i][0]/log(10), lprob[i][1]/log(10),
							 lprob[i][2]/log(10));
			fputc('\n', fp[i]);
			if (!ref_is_gap[i])
				++pos_ref[i];
		}

		for (int g1 = 0; g1 <= 2; ++g1)
			for (int g2 = 0; g2 <= 2; ++g2)
 				fprintf(stderr, "M = (%u, %u) %c%c/%c%c at site"
					"%u: %e\n", g1, g2,
					g1 == 2 ? iupac_to_char[nucA2b]
							: iupac_to_char[nucA1],
					g1 ? iupac_to_char[nucA2b]
							: iupac_to_char[nucA1],
					g2 == 2 ? iupac_to_char[nucB2b]
							: iupac_to_char[nucB1],
					g2 ? iupac_to_char[nucB2b]
							: iupac_to_char[nucB1],
					j + 1, gprob[0][g1] * gprob[1][g2]);

		fprintf(stderr, "Genotype (%4u): %c%c/%c%c (%f) [", j + 1,
			max_g[0] == 2 ? iupac_to_char[nucA2b]
						: iupac_to_char[nucA1],
			max_g[0] ? iupac_to_char[nucA2b]
						: iupac_to_char[nucA1],
			max_g[1] == 2 ? iupac_to_char[nucB2b]
						: iupac_to_char[nucB1],
			max_g[1] ? iupac_to_char[nucB2b]
						: iupac_to_char[nucB1],
			gprob[0][max_g[0]] * gprob[1][max_g[1]]);
		for (int g1 = 0; g1 <= 2; ++g1)
			for (int g2 = 0; g2 <= 2; ++g2)
				fprintf(stderr, " %f", gprob[0][g1]
							* gprob[1][g2]);
		fprintf(stderr, "]%s\n", max_g[0] == 1 || max_g[1] == 1 ? "***"
			: abs(max_g[0] - max_g[1]) > 1 ? "+++" : "");

#ifdef ALLOW_DEBUGGING
		if (opt->debugging_site >= 0 && opt->debugging_site < (int)j)
			exit(0);
#endif
	}

GENOTYPE_BY_CLUSTERING_RETURN:
	for (int i = 0; i < N_SUBGENOMES; ++i)
		if (fp[i])
			fclose(fp[i]);
	
	return err;
}/* genotype_by_clustering */

/**
 * Align each subgenic reference to the profile alignment of haplotypes
 * assigned to that subgenome.  In allotetraploid, there are at most
 * two haplotypes assigned to each subgenome, so this is equivalent to
 * aligning to IUPAC-encoded sequences.  We assume there is a match
 * if the reference matches each nucleotide in profile.  Each non-
 * matching allele costs a mismatch penalty.  Deletions in the profile
 * cost nothing; they have already been penalized.
 *
 * @param in			pointer to input object
 * @param fd			alignment of selected haplotypes
 * @param subgenomic_profiles	reference to profile alignments
 * @param rseqs			reference sequences
 * @param rlens			lengths of reference sequences
 * @return			error status
 */
int get_reference_aligned_haplotypes(input *in, fastq_data *fd,
	char_t const *subgenomic_profiles[2*N_SUBGENOMES],
		char_t const * const rseqs[N_SUBGENOMES],
		unsigned int const rlens[N_SUBGENOMES])
{
	int err = NO_ERROR;
	unsigned int alen;
	double ascore;

	for (int i = 0; i < N_SUBGENOMES; ++i) {
		in->ref_profile_alignments[i] = nw_seq_profile_align(
			rseqs[i], subgenomic_profiles[i],
			rlens[i], fd->n_max_length, 1, -2, -2, -1, 0,
						&err, &alen, &ascore);
		if (err)
			break;
	}

	return err;
} /* get_reference_aligned_haplotypes */


/**
 * Align and reload chosen haplotypes.  Also counts indels within each
 * subgenomic alignment.
 *
 * @param in		pointer to input object
 * @param opt		pointer to options object
 * @param fd		address of pointer to fasta object
 * @param haplotypes	pointer to start of each aligned haplotype sequence
 * @return		error status
 */
int get_aligned_haplotypes(input *in, options *opt, fastq_data **fd,
	char_t const *haplotypes[2*N_SUBGENOMES])
{
	int err = NO_ERROR;

	/* open AmpliCI's FASTA file with inferred haplotypes [TODO, LIBRARY] */
	char *fsa_file = malloc((strlen(opt->ac_outfile) + 4)
							* sizeof *fsa_file);
	fastq_options fop = {.read_encoding = XY_ENCODING, .read_names = 1};

	if (!fsa_file)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"AmpliCI FASTA outfile.\n");

	sprintf(fsa_file, "%s.fa", opt->ac_outfile);

	if ((err = read_fastq(fsa_file, fd, &fop)))
		return err;

	free(fsa_file);

	/* identify which of the 2, 3, or 4 retained top haplotypes belong to
	 * each subgenome; store ids in input::csome_to_inc, maybe duplicating
	 * if there is only one cluster assigned to a subgenome (BEWARE!)
	 */
	if ((err = load_chromosomes(opt, in)))
		return err;

for (int i = 0; i < 4; ++i) {
	fprintf(stderr, "Csome %d maps to include id %u.\n", i, in->csome_to_inc[i]);
	fprintf(stderr, "Inc id %d to cluster id %u\n", i, in->included_clusters[i]);
	fprintf(stderr, "Proportion A %d: %f\n", i, in->proportion_A[i]);
}

	/* write fasta file with retained haplotypes */
	unsigned int new_cluster_id[4];	/* map include_id (0, 1, ..., keep) to alignment file index */
	unsigned int n_found = 0;
	unsigned int n_kept = 2*N_SUBGENOMES - in->proptest;

	/* find new index of kept clusters in alignment file */
	for (unsigned int k = 0; k < in->K; ++k) {
		if (in->excluded_clusters[k])
			continue;
		for (unsigned int j = 0; j < n_kept; ++j)
			if (in->included_clusters[j] == k) {
				new_cluster_id[j] = n_found++;
				break;
			}
	}

	unsigned int *selected = calloc((*fd)->n_reads, sizeof *selected);
	selected[in->included_clusters[in->csome_to_inc[0]]] = 1;
	selected[in->included_clusters[in->csome_to_inc[1]]] = 1;
	selected[in->included_clusters[in->csome_to_inc[2]]] = 1;
	selected[in->included_clusters[in->csome_to_inc[3]]] = 1;
	fop.fasta = 1;
	fop.outfile = opt->mafft_command ? opt->mafft_infile
		: opt->clustalo_infile;

	write_fastq_marked_trimmed(*fd, &fop, selected, 1, NULL, NULL);

	free(selected);

	/* align selected haplotypes and reload into fd object */
	if ((err = align_haplotypes(fd, opt)))
		return err;

	/* get pointers to aligned haplotypes for each csome */
	for (unsigned int k = 0; k < 2*N_SUBGENOMES; ++k)
		haplotypes[k] = new_cluster_id[in->csome_to_inc[k]]
					* (*fd)->n_max_length + (*fd)->reads;
	
	for (unsigned int j = 0; j < (*fd)->n_max_length; ++j) {
		char_t nucA1 = haplotypes[0][j];
		char_t nucA2 = haplotypes[1][j];
		char_t nucB1 = haplotypes[2][j];
		char_t nucB2 = haplotypes[3][j];
#ifdef ALLOW_DEBUGGING
		if (opt->debugging_site == (int) j)
			fprintf(stderr, "Debugging site %u: %c%c/%c%c\n", opt->debugging_site, iupac_to_char[nucA1], iupac_to_char[nucA2], iupac_to_char[nucB1], iupac_to_char[nucB2]);
#endif

		if (nucA1 != nucA2 && (!nucA1 || !nucA2))
			++in->n_indels[0];
		if (nucB1 != nucB2 && (!nucB1 || !nucB2))
			++in->n_indels[1];
	}

	return err;
} /* get_aligned_haplotypes */


/**
 * Based on results of coverage tests, get indices of the allotetraploid
 * chromosomes.  See roshan.h for mappings involved.
 *
 * @param in	pointer to input object
 * @return 	error status
 */
int load_chromosomes(options *opt, input *in)
{
	if (in->reject_state == ACCEPT_FOUR) {
		if (in->proportion_A[0] > 0.5
			&& in->proportion_A[1] > 0.5) {
			in->csome_to_inc[0] = 0;
			in->csome_to_inc[1] = 1;
			in->csome_to_inc[2] = 2;
			in->csome_to_inc[3] = 3;
		} else if (in->proportion_A[0] > 0.5
			&& in->proportion_A[2] > 0.5) {
			in->csome_to_inc[0] = 0;
			in->csome_to_inc[1] = 2;
			in->csome_to_inc[2] = 1;
			in->csome_to_inc[3] = 3;
		} else if (in->proportion_A[0] > 0.5
			&& in->proportion_A[3] > 0.5) {
			in->csome_to_inc[0] = 0;
			in->csome_to_inc[1] = 3;
			in->csome_to_inc[2] = 1;
			in->csome_to_inc[3] = 2;
		} else if (in->proportion_A[1] > 0.5
			&& in->proportion_A[2] > 0.5) {
			in->csome_to_inc[0] = 1;
			in->csome_to_inc[1] = 2;
			in->csome_to_inc[2] = 0;
			in->csome_to_inc[3] = 3;
		} else if (in->proportion_A[1] > 0.5
			&& in->proportion_A[3] > 0.5) {
			in->csome_to_inc[0] = 1;
			in->csome_to_inc[1] = 3;
			in->csome_to_inc[2] = 0;
			in->csome_to_inc[3] = 2;
		} else if (in->proportion_A[2] > 0.5
			&& in->proportion_A[3] > 0.5) {
			in->csome_to_inc[0] = 2;
			in->csome_to_inc[1] = 3;
			in->csome_to_inc[2] = 0;
			in->csome_to_inc[3] = 1;
		}
	} else if (in->reject_state == ACCEPT_THREE_A) {
		if (in->proportion_A[0] > 0.5) {	/* A, B, B */
			in->csome_to_inc[0] = in->csome_to_inc[1] = 0;
			in->csome_to_inc[2] = 1;
			in->csome_to_inc[3] = 2;
		} else {				/* B, A, A */
			in->csome_to_inc[0] = 1;
			in->csome_to_inc[1] = 2;
			in->csome_to_inc[2] = in->csome_to_inc[3] = 0;
		}
	} else if (in->reject_state == ACCEPT_THREE_B) {
		if (in->proportion_A[0] > 0.5) {	/* A, A, B */
			in->csome_to_inc[0] = 0;
			in->csome_to_inc[1] = 1;
			in->csome_to_inc[2] = in->csome_to_inc[3] = 2;
		} else {				/* B, B, A */
			in->csome_to_inc[0] = in->csome_to_inc[1] = 2;
			in->csome_to_inc[2] = 0;
			in->csome_to_inc[3] = 1;
		}
	} else if (in->reject_state == ACCEPT_TWO) {
		if (in->proportion_A[0] > 0.5) {
			in->csome_to_inc[0] = in->csome_to_inc[1] = 0;
			in->csome_to_inc[2] = in->csome_to_inc[3] = 1;
		} else {
			in->csome_to_inc[0] = in->csome_to_inc[1] = 1;
			in->csome_to_inc[2] = in->csome_to_inc[3] = 0;
		}
	} else {
		return mmessage(ERROR_MSG, INTERNAL_ERROR,
					"Invalid reject state.\n");
	}

	/* when genotyping by clustering, we need to map 
	 * cluster id to chromosome to know alignment of haplotypes
	 */
	if (opt->genotype_by_clustering)
		for (unsigned int k = 0; k < 2*N_SUBGENOMES; ++k)
			in->cluster_to_csome[in->included_clusters[
						in->csome_to_inc[k]]] = k;


	return NO_ERROR;
} /* load_chromosomes */


/**
 * Align the selected sequences and reload the alignment.
 *
 * We prefer mafft over clustal omega as in our experience it routinely
 * finds better alignments (smaller pairwise distances).  We like muscle
 * too, but there is no option to retain sequence input order, and it
 * is a minor hassle to reorder that we currently avoid.
 *
 * @param fd		pointer to unaligned data
 * @param opt		pointer to options with information about aligner
 * @return		error status
 */
int align_haplotypes(fastq_data **fd, options *opt)
{
	int err = NO_ERROR;
	char *cmd = NULL;
	fastq_options fop = {.read_encoding = XY_ENCODING, .read_names = 1,
		.fasta = 1};

	if (opt->mafft_command) {
		fop.outfile = opt->mafft_infile;
		cmd = malloc((strlen(opt->mafft_command)
			+ strlen(opt->mafft_infile) + strlen(opt->mafft_outfile)
			+ strlen(" -i  -o  --output-order=input-order --force "))
			* sizeof *cmd);
		if (!cmd)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "cmd");
		sprintf(cmd, "%s --retree 2 --inputorder %s > %s",
					opt->mafft_command, opt->mafft_infile,
							opt->mafft_outfile);
	} else if (opt->clustalo_command) {
		cmd = malloc((strlen(opt->clustalo_command)
			+ strlen(opt->clustalo_infile)
			+ strlen(opt->clustalo_outfile)
			+ strlen(" -i  -o  --output-order=input-order --force "))
			* sizeof *cmd);
		if (!cmd)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "cmd");
		sprintf(cmd, "%s -i %s -o %s --output-order=input-order --force",
				opt->clustalo_command, opt->clustalo_infile,
							opt->clustalo_outfile);
	}

	mmessage(INFO_MSG, NO_ERROR, "Running aligner: \"%s\"\n", cmd);
	system(cmd);
	free_fastq(*fd);
	*fd = NULL;
	fop.fasta = 0;
	fop.read_encoding = IUPAC_ENCODING;
	err = read_fastq(opt->mafft_command ? opt->mafft_outfile
				: opt->clustalo_outfile, fd, &fop);
	free(cmd);

	return err;
} /* align_haplotypes */


/********** multinomial logistic regression error model **********/

/**
 * Read error parameter file from mlogit fit.
 *
 * Format per line:
 * 64 C 1226 C 0.997768688734120945
 * quality_score true_base reference_position read_base prob_no_error
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


/**
 * Probability of substituting read nucleotide r for true nucleotide s
 * given covariate quality score q under multinomial logistic regression
 * fit.
 *
 * Satisfies function pointer typedef defined in qual.h:
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


/**
 * Make input data for results given by amplici.
 *
 * @param in	pointer to input structure
 * @return error status
 */
int make_input(input **in) {

	*in = malloc(sizeof **in);

	if (*in == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "input object");

	(*in)->assignment = NULL;
	(*in)->excluded_clusters = NULL;
	(*in)->cluster_to_csome = NULL;
	(*in)->cluster_sizes = NULL;
	(*in)->sort_id = NULL;
	(*in)->n_excluded = 0;
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

	
	if (in->cluster_sizes)
		free(in->cluster_sizes);
	if (in->assignment)
		free(in->assignment);
	if (in->sort_id)
		free(in->sort_id);
/*
	if (in->exclude)
		free(in->exclude);
	if (in->exclude_id)
		free(in->exclude_id);
*/
	if (in->cluster_to_csome)
		free(in->cluster_to_csome);
	if (in->excluded_clusters)
		free(in->excluded_clusters);

	free(in);
}/* free_input */

/**
 * Set default options.
 *
 * @param opt	options pointer
 * @return	error status
 */
int default_options(options *opt)
{
	opt->amplicon = 1;
	opt->display_alignment = 0;
	opt->drop_unmapped = 1;
	opt->drop_secondary = 1;
	opt->drop_soft_clipped = INT_MAX;
	opt->drop_indel = INT_MAX;
	opt->proptest_screen = -1;
	opt->proptest_alpha = 0.05;
	opt->proptest_min_freq = 0;
	opt->coverage_screen = 0.9;	/* default changed from 0.5 on 5/27/20 */
	opt->biallelic_screen = 0.5;
	opt->max_misalignment_proportion = 0.3;
	opt->min_length = 0;
	opt->max_length = INT_MAX;
	opt->min_log_likelihood = -INFINITY;
	opt->min_posterior_alignment_prob = 0;
	opt->n_sample = 100;
	opt->max_eerr = INFINITY;
	opt->param_file = NULL;
	opt->error_file = NULL;
	opt->max_quality_score = INT_MAX;
	opt->use_bam = 0;
	for (int i = 0; i < N_SUBGENOMES; ++i) {
		opt->sbam_files[i] = NULL;
		opt->fsa_files[i] = NULL;
		opt->ref_names[i] = NULL;
		opt->vcf_files[i] = NULL;
		opt->subref_fsas[i] = NULL;
	}
	opt->sample_name = NULL;
	opt->genotype_by_clustering = 0;
	opt->amplici_command = NULL;
	opt->amplici_complete = 1;
	opt->ac_fastq_file = "amplici.fastq";
	opt->ac_outfile = "amplici";
	opt->ac_low_bound = 1.5;
	opt->ac_ll_threshold = -INFINITY;
	opt->clustalo_command = NULL;	//"/bin/clustalo";
	opt->clustalo_infile = "selected_haplotypes.fa";
	opt->clustalo_outfile = "selected_haplotypes.co.fa";
	opt->mafft_command = NULL;
	opt->mafft_infile = "selected_haplotypes.fa";
	opt->mafft_outfile = "selected_haplotypes.mft.fa";
	opt->coverage_file = "coverage.txt";
	opt->write_fastq_and_quit = 0;
	opt->fastq_after_filter = NULL;
	opt->min_expected_coverage = 5.;

	opt->subref_fsa_base = "subsetted_refs";

	opt->posthoc_coverage_test = 0;
	opt->equal_homolog_coverage_test = 0;
	opt->phc_min_alignment_pp = 0.99;
	opt->phc_min_genotype_pp = 0.99;

	opt->debugging_site = -1;
	opt->karin_old = 0;

	opt->vcf_opt = NULL;
	make_default_vcf_options(&opt->vcf_opt);

	return NO_ERROR;
} /* default_options */

/**
 * Parse command line.
 *
 * @param opt	options pointer
 * @param argc	number of words on the command line
 * @param argv	the words on the command line
 * @return	error status
 */
int parse_options(options *opt, int argc, const char **argv)
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
		case 'a':
			if (i + 1 >= argc)
				goto CMDLINE_ERROR;
			if (!strncmp(&argv[i][j], "amplici_o", 9)
				|| !strncmp(&argv[i][j], "amplici-o", 9)
				|| !strncmp(&argv[i][j], "ampliclust_o", 12)
				|| !strncmp(&argv[i][j], "ampliclust-o", 12)) {
				opt->ac_outfile = argv[++i];
				mmessage(INFO_MSG, NO_ERROR, "Amplici "
					"output file base name: '%s'\n",
						opt->ac_outfile);
			} else if (!strncmp(&argv[i][j], "amplici_f", 9)
				|| !strncmp(&argv[i][j], "amplici-f", 9)
				|| !strncmp(&argv[i][j], "ampliclust_f", 12)
				|| !strncmp(&argv[i][j], "ampliclust-f", 12)) {
				opt->ac_fastq_file = argv[++i];
				mmessage(INFO_MSG, NO_ERROR, "Amplici "
					"fastq filename: '%s'\n",
						opt->ac_fastq_file);
			} else if (!strncmp(&argv[i][j], "amplici_l", 9)
				|| !strncmp(&argv[i][j], "amplici-l", 9)
				|| !strncmp(&argv[i][j], "ampliclust_l", 12)
				|| !strncmp(&argv[i][j], "ampliclust-l", 12)) {
				opt->ac_low_bound = read_cmdline_double(argc,
					argv, ++i, opt);
				mmessage(INFO_MSG, NO_ERROR, "Amplici low"
					 "er bound: %f\n", opt->ac_low_bound);
			} else if (!strncmp(&argv[i][j], "amplici_m", 9)
				|| !strncmp(&argv[i][j], "amplici-m", 9)
				|| !strncmp(&argv[i][j], "ampliclust_m", 12)
				|| !strncmp(&argv[i][j], "ampliclust-m", 12)) {
				opt->amplici_complete = 0;
				opt->proptest_min_freq = read_uint(argc, argv,
								++i, opt);
				if (opt->proptest_min_freq == 0)
					opt->amplici_complete = 1;
				if (opt->amplici_complete)
					mmessage(INFO_MSG, NO_ERROR,
						"Subgenomic coverage test based"
						" on AmpliCI-estimated cluster "
						"abundances.\n");
				else
					mmessage(INFO_MSG, NO_ERROR,
						"Subgenomic coverage test based"
						" on AmpliCI hash counts if "
						"maximum count exceeds %u.\n",
						opt->proptest_min_freq);
			} else if (!strncmp(&argv[i][j], "amplici_t", 9)
				|| !strncmp(&argv[i][j], "amplici-t", 9)
				|| !strncmp(&argv[i][j], "ampliclust_t", 12)
				|| !strncmp(&argv[i][j], "ampliclust-t", 12)) {
				opt->ac_ll_threshold = read_cmdline_double(
							argc, argv, ++i, opt);
				mmessage(INFO_MSG, NO_ERROR, "Amplici lower "
					"bound on log likelihood: %f\n",
					opt->ac_ll_threshold);
			} else if (!strncmp(&argv[i][j], "ali", 3)) {
				opt->mafft_infile = opt->clustalo_infile
					= argv[++i];
				if (i + 1 >= argc)
					goto CMDLINE_ERROR;
				opt->mafft_outfile = opt->clustalo_outfile
					= argv[++i];
			} else {
				opt->amplici_command = argv[++i];
				mmessage(INFO_MSG, NO_ERROR, "Amplici "
							"command: '%s'\n",
						opt->amplici_command);
			}
			break;
		case 'b':
			if (!strncmp(&argv[i][j], "bam", 3)) {
				if (i + N_SUBGENOMES >= argc) {
					err = mmessage(ERROR_MSG,
					       INVALID_CMD_ARGUMENT, "Too few "
					       "arguments to --bam_files "
					       "command-line option.\n");
					goto CMDLINE_ERROR;
				}
				opt->use_bam = 1;
				mmessage(INFO_MSG, NO_ERROR, "BAM files:");
				for (j = 0; j < N_SUBGENOMES; ++j) {
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
			if (!strncmp(&argv[i][j], "dis", 3)) {
				opt->display_alignment = !opt->display_alignment;
				mmessage(INFO_MSG, NO_ERROR, "Display "
					"alignments on stderr? %s\n",
					opt->display_alignment ? "yes" : "no");
			} else if (!strncmp(&argv[i][j], "deb", 3)) {
				if (i + 1 >= argc)
					goto CMDLINE_ERROR;
#ifdef ALLOW_DEBUGGING
				opt->debugging_site = read_int(argc, argv, ++i,
									opt);
				--opt->debugging_site;
				mmessage(INFO_MSG, NO_ERROR, "Debugging "
					"genotype call at site %u.\n",
							opt->debugging_site);
#else
				++i;
				mmessage(INFO_MSG, NO_ERROR, "Debugging "
					"disabled.  Ignoring use of "
					"--debugging_site option.\n");
#endif
			} else if (!strncmp(&argv[i][j], "drop", 4)) {
				if (i + 1 >= argc)
					goto CMDLINE_ERROR;
				opt->proptest_screen = read_int(argc, argv, ++i,
									opt);
				mmessage(INFO_MSG, NO_ERROR, "Dropping reads "
					"assigned to all but %d most-abundant "
					"haplotypes\n", 4 - opt->proptest_screen);
			} else {
				goto CMDLINE_ERROR;
			}
			break;
		case 'e':
			if (!strncmp(&argv[i][j], "ex", 2)) {
				opt->max_eerr =
				read_cmdline_double(argc, argv, ++i, opt);
				mmessage(INFO_MSG, NO_ERROR, "Dropping reads "
					 "with more than %f expected errors.\n",
					 opt->max_eerr);
			} else if (!strncmp(&argv[i][j], "eq", 2)) {
				opt->posthoc_coverage_test = 1;
				mmessage(INFO_MSG, NO_ERROR, "Performing "
					"post-hoc confidence interval of equal "
					"homologous chromosome coverage\n");
				if (i + 1 < argc && argv[i][0] != '-') {
					opt->phc_min_alignment_pp =
						read_cmdline_double(argc, argv,
							++i, opt);
					mmessage(INFO_MSG, NO_ERROR, "Minimum "
						"alignment posterior "
						"probability: %.3f\n",
						opt->phc_min_alignment_pp);
				}
				if (i + 1 < argc && argv[i][0] != '-') {
					opt->phc_min_genotype_pp =
						read_cmdline_double(argc, argv,
							++i, opt);
					mmessage(INFO_MSG, NO_ERROR, "Minimum "
						"genotype posterior "
						"probability: %.3f\n",
						opt->phc_min_genotype_pp);
				}
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
			if (i + N_SUBGENOMES >= argc) {
				err = mmessage(ERROR_MSG, INVALID_CMD_ARGUMENT,
				       "Too few arguments to --fasta_files "
					       "command-line option.\n");
				goto CMDLINE_ERROR;
			}
			mmessage(INFO_MSG, NO_ERROR, "Fasta files:");
			for (j = 0; j < N_SUBGENOMES; ++j) {
				opt->fsa_files[j] = argv[++i];
				fprintf(stderr, " %s",
					opt->fsa_files[j]);
			}
			fprintf(stderr, "\n");
			break;
		case 'g':
			if (!strncmp(&argv[i][j], "gen", 3)) {
				opt->genotype_by_clustering = 1;
				mmessage(INFO_MSG, NO_ERROR, "Will do genotype "
							"by clustering.\n");
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
			if (i + 1 >= argc)
				goto CMDLINE_ERROR;
			opt->drop_indel = read_int(argc, argv, ++i, opt);
			mmessage(INFO_MSG, NO_ERROR, "Dropping reads with indel"
				 " longer than %d in either alignment.\n",
				 opt->drop_indel);
			break;

		case 'k':
			opt->karin_old = 1;
			break;

		case 'l':
			if (i + 1 >= argc)
				goto CMDLINE_ERROR;
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
			} else if (!strncmp(&argv[i][j], "mis", 3)) {
				opt->max_misalignment_proportion
					= read_cmdline_double(argc, argv, ++i,
									 opt);
				if (opt->max_misalignment_proportion >= 0.5)
					opt->max_misalignment_proportion = 0.49;
				mmessage(INFO_MSG, NO_ERROR, "Maximum "
					"misalignment rate: %d\n",
					 opt->max_misalignment_proportion);
			} else if (!strncmp(&argv[i][j], "min_p", 5)
				|| !strncmp(&argv[i][j], "min-p", 5)) {
				opt->min_posterior_alignment_prob
					= read_cmdline_double(argc, argv, ++i,
									opt);
				mmessage(INFO_MSG, NO_ERROR, "Minimum posterior"
					" probability of alignment: %f.\n",
					opt->min_posterior_alignment_prob);
			} else if (!strncmp(&argv[i][j], "min", 3)) {
				opt->min_length = read_int(argc, argv, ++i, opt);
				mmessage(INFO_MSG, NO_ERROR, "Minimum read "
					 "length: %d\n", opt->min_length);
			} else if (!strncmp(&argv[i][j], "max", 3)) {
				opt->max_length = read_int(argc, argv, ++i, opt);
				mmessage(INFO_MSG, NO_ERROR, "Maximum read "
					 "length: %d\n", opt->max_length);
			} else if (!strncmp(&argv[i][j], "maf", 3)) {
				opt->mafft_command = argv[++i];
				mmessage(INFO_MSG, NO_ERROR, "Using mafft "
							"command '%s'.\n",
							opt->mafft_command);
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
			if (!strncmp(&argv[i][j], "ref_a", 5)) {
				if (i + 1 >= argc) {
					err = mmessage(ERROR_MSG,
						INVALID_CMD_ARGUMENT, "Option "
						"--ref_align needs argument.\n");
					goto CMDLINE_ERROR;
				}
				opt->ref_alignment = argv[++i];
				mmessage(INFO_MSG, NO_ERROR, "Reference "
					"alignment: %s", opt->ref_alignment);
			} else if (!strncmp(&argv[i][j], "ref_f", 5)) {
				if (i + N_SUBGENOMES >= argc) {
					err = mmessage(ERROR_MSG,
						INVALID_CMD_ARGUMENT, "Too few "
						"arguments to --ref_fsa "
						"command-line option.\n");
					goto CMDLINE_ERROR;
				}
				mmessage(INFO_MSG, NO_ERROR, "Reference fasta:");
				for (j = 0; j < N_SUBGENOMES; ++j) {
					opt->ref_fsas[j] = argv[++i];
					fprintf(stderr, " %s",
						opt->ref_fsas[j]);
				}
				fprintf(stderr, "\n");
			} else if (!strncmp(&argv[i][j], "ref_n", 5)) {
				if (i + N_SUBGENOMES >= argc) {
					err = mmessage(ERROR_MSG,
						INVALID_CMD_ARGUMENT, "Too few "
						"arguments to --ref_names "
						"command-line option.\n");
					goto CMDLINE_ERROR;
				}
				mmessage(INFO_MSG, NO_ERROR, "Reference names:");
				for (j = 0; j < N_SUBGENOMES; ++j) {
					opt->ref_names[j] = argv[++i];
					fprintf(stderr, " %s",
						opt->ref_names[j]);
				}
				fprintf(stderr, "\n");
			}
			break;
		case 's':
			if (!strncmp(&argv[i][j], "se", 2)) {
				opt->drop_secondary = !opt->drop_secondary;
				mmessage(INFO_MSG, NO_ERROR, "%s secondary "
					 "alignments\n", opt->drop_secondary
					 ? "Dropping" : "Keeping");
			} else if (!strncmp(&argv[i][j], "so", 2)) {
				opt->drop_soft_clipped
				= read_int(argc, argv, ++i, opt);
				mmessage(INFO_MSG, NO_ERROR, "Dropping reads "
					 "with soft clip >= %d in either "
					 "alignment.\n", opt->drop_soft_clipped);
			} else if (!strncmp(&argv[i][j], "samp", 4)) {
				opt->n_sample = read_uint(argc, argv, ++i, opt);
				mmessage(INFO_MSG, NO_ERROR, "%u Monte Carlo "
					 "samples.\n", opt->n_sample);
			} else if (!strncmp(&argv[i][j], "sam", 3)) {
				if (i + N_SUBGENOMES >= argc) {
					err = mmessage(ERROR_MSG,
						       INVALID_CMD_ARGUMENT, "Too few "
						       "arguments to --sam_files "
						       "command-line option.\n");
					goto CMDLINE_ERROR;
				}
				mmessage(INFO_MSG, NO_ERROR, "Sam files:");
				for (j = 0; j < N_SUBGENOMES; ++j) {
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

		case 'w':
			if (!strncmp(&argv[i][j], "write-f", 7)) {
				if (i + 1 < argc && argv[i+1][0] != '-') {
					opt->fastq_after_filter = argv[++i];
					mmessage(INFO_MSG, NO_ERROR, "Will write "
						"filtered reads in fastq file '%s'\n",
						opt->fastq_after_filter);
				}
			} else if (!strncmp(&argv[i][j], "write-s", 7)) {
				mmessage(INFO_MSG, NO_ERROR, "Will write "
					"filtered reads to sam files: ");
				for (j = 0; j < N_SUBGENOMES; ++j) {
					if (i + 1 >= argc) {
						err = mmessage(ERROR_MSG,
							INVALID_CMD_ARGUMENT,
							"Too few arguments to "
							"--write-sam option\n");
						goto CMDLINE_ERROR;
					}
					opt->sam_after_filter[j] = argv[++i];
					mmessage(INFO_MSG, NO_ERROR, " %s",
						opt->sam_after_filter[j]);
				}
				mmessage(INFO_MSG, NO_ERROR, "\n");
			} else {
				opt->write_fastq_and_quit = 1;
				mmessage(INFO_MSG, NO_ERROR, "Will write "
					"selected reads in fastq file and "
					"quit.\n");
			}
			break;

		case 'p':
			if (!strncmp(&argv[i][j], "po", 2)) {
				opt->equal_homolog_coverage_test = 1;
				mmessage(INFO_MSG, NO_ERROR, "Will run "
						"equal coverage test.\n");
				break;
			}
			if (i + 1 >= argc)
				goto CMDLINE_ERROR;
			opt->proptest_screen = read_int(argc, argv, ++i, opt);
			mmessage(INFO_MSG, NO_ERROR, "Dropping reads assigned "
				"to all but %d most-abundant haplotypes\n",
				 4 - opt->proptest_screen);
			break;

		case 'v':
			if (i + N_SUBGENOMES >= argc)
				goto CMDLINE_ERROR;
			for (j = 0; j < N_SUBGENOMES; ++j)
				opt->vcf_files[j] = argv[++i];
			if (!opt->vcf_opt && (err
				= make_default_vcf_options(&opt->vcf_opt)))
				goto CMDLINE_ERROR;
			break;
		case 'c':
			if (i + 1 >= argc)
				goto CMDLINE_ERROR;
			if (!strncmp(&argv[i][j], "cov", 3)) {
				opt->coverage_screen = read_cmdline_double(
							argc, argv, ++i, opt);
				mmessage(INFO_MSG, NO_ERROR, "Dropping sites "
					"with coverage drop above: %f\n",
							opt->coverage_screen);
			} else if (!strncmp(&argv[i][j], "clu", 3)) {
				opt->clustalo_command = argv[++i];
				mmessage(INFO_MSG, NO_ERROR, "Using clustalo "
							"command '%s'.\n",
							opt->clustalo_command);
			} else {
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
			}
			break;

		default:
			err = INVALID_CMD_OPTION;
			goto CMDLINE_ERROR;
		}
	}

	struct stat st;
	if (opt->amplici_command) 
		if (stat(opt->amplici_command, &st) || !(st.st_mode & S_IXUSR))
			err = mmessage(ERROR_MSG, INVALID_CMDLINE, "AmpliCI "
				"command '%s' is not executable.\n",
						opt->amplici_command);

	if (opt->genotype_by_clustering && !opt->amplici_command)
		err = mmessage(ERROR_MSG, INVALID_CMDLINE, "To perform "
			"genotype-by-clustering, you must provide the AmpliCI "
			"command with option --amplici\n");

	if (opt->genotype_by_clustering && !opt->clustalo_command
							&& !opt->mafft_command)
		err = mmessage(ERROR_MSG, INVALID_CMDLINE, "To perform "
			"genotype-by-clustering, you must provide an aligner: "
			"see --clustalo or --mafft.\n");

	if (opt->fastq_after_filter && opt->genotype_by_clustering)
		err = mmessage(ERROR_MSG, INVALID_CMDLINE, "Cannot write "
			"filtered reads when genotyping-by-clustering.\n");

	if (opt->fastq_after_filter && !opt->amplici_command)
		err = mmessage(ERROR_MSG, INVALID_CMDLINE, "Cannot write "
			"filtered reads without running AmpliCI.\n");

	if (!opt->mafft_command && opt->clustalo_command &&
		(stat(opt->clustalo_command, &st) || !(st.st_mode & S_IXUSR)))
		err = mmessage(ERROR_MSG, INVALID_CMDLINE, "Clustalo "
			"command '%s' is not executable.\n",
					opt->clustalo_command);

	if (opt->mafft_command && (stat(opt->mafft_command, &st)
						|| !(st.st_mode & S_IXUSR)))
		err = mmessage(ERROR_MSG, INVALID_CMDLINE, "Clustalo command "
			"'%s' is not executable.\n", opt->mafft_command);

	opt->vcf_opt->genotype_by_clustering = opt->genotype_by_clustering;
	opt->vcf_opt->coverage_screen = opt->coverage_screen;
	opt->vcf_opt->posthoc_coverage_test = opt->posthoc_coverage_test;
	opt->vcf_opt->phc_min_genotype_pp = opt->phc_min_genotype_pp;
	opt->vcf_opt->equal_homolog_coverage_test = opt->equal_homolog_coverage_test;

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
 * Print program usage, called from command-line utilities.
 *
 * @param fp		print to this file handle
 * @param cmdname	name of the executable
 * @param obj		void point to options object
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
	fprintf(fp, "\nNAME\n\t%s - genotype tetraploids\n", &cmdname[start]);
	fprintf(fp, "\nSYNOPSIS\n"
		"\t%s --sam_files SAM1 SAM2 --fasta_files FSA1 FSA2 --ref_names REF1 REF2\n"
		"\t\t[[--genotype_by_clustering [--alignment FILE1 FILE2]]\n"
		"\t\t[--sample INT --min-subgenomic-coverage FLOAT]\n"
		"\t\t[--min INT --max INT --expected-errors FLOAT --indel INT --loglik FLOAT\n"
		"\t\t --min-posterior FLOAT --secondary --soft-clipped INT]\n"
		"\t\t[--coverage FLOAT --biallelic FLOAT --equal_coverage_test [FLOAT1 FLOAT2]]\n"
		"\t\t[--drop INT --amplici EXE [--amplici-f FILE --amplici-o STRING --amplici-l FLOAT]]\n"
		"\t\t[--error_file|--error_data FILE] ...]\n", &cmdname[start]);
/***************************ruler*for*80*character*line*************************************/
	fprintf(fp, "\nDESCRIPTION\n"
		"\t%s genotypes allotetraploids using reads in SAM1 and SAM2 aligned to\n"
		"\tREF1 and REF2 references from fasta files FSA1 FSA2.\n", &cmdname[start]);
	fprintf(fp,
		"\tSAM files typically contain reads from a single individual, genotype, or\n"
		"\taccession aligned to multiple amplified targets, but %s genotypes\n"
		"\tone individual at one amplicon.\n", &cmdname[start]);
	fprintf(fp, "\nOPTIONS\n");
//	fprintf(fp, "\t+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	fprintf(fp, "\nInput (required):\n");
//	fprintf(fp, "\t+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	fprintf(fp, "\t--fasta_files FILE1 FILE2\n"
		"\t\tSubgenomic reference fasta files (Default: none).\n");
	fprintf(fp, "\t\tDEPRECATED: see --ref_fasta_files\n");
	fprintf(fp, "\t--ref_fasta_files FILE1 FILE2\n"
		"\t\tSubgenomic reference fasta files (Default: none).\n");
	fprintf(fp, "\t--sam_files FILE1 FILE2\n"
		"\t\tSAM files with reads aligned to each subgenome (Default: none)\n");
/* HIDDEN OPTION
 *	fprintf(fp, "\t--bam_files FILE1 FILE2\n"
 *		"\t\tBAM files with alignments to each subgenome (Default: none)\n");
 */
	fprintf(fp, "\t--ref_names STRING1 STRING2\n"
		"\t\tNames of subgenomic reference target regions (Default: none)\n");
	fprintf(fp, "\t--ref_alignment FILE\n"
		"\t\tSAM file containing alignment of references (Default: none)\n");

//	fprintf(fp, "\t+++++++++++++++++\n");
	fprintf(fp, "\nOutput (optional):\n");
//	fprintf(fp, "\t+++++++++++++++++\n");
/* UNIMPLEMENTED OPTION
 *	fprintf(fp, "\t--censor INT\n"
 *		"\t\tCode censors quality scores at maximum 41 [Option not used!]\n");
 */
	fprintf(fp, "\t--display_alignment\n"
		"\t\tDisplay alignments in stderr output (Default: %s).\n",
					opt->display_alignment ? "yes" : "no");
	fprintf(fp, "\t--vcf_files FILE1 FILE2\n"
		"\t\tGenotyping output in one vcf file per subgenome (Default: none).\n");
	fprintf(fp, "\t--subref_fasta_files FILE|FILE1 FILE2\n"
		"\t\tSubsetted reference regions output to these FASTA files (Default: %s[12].fsa)\n", opt->subref_fsa_base);
	fprintf(fp, "\t--gl\n"
		"\t\tToggle GL output to vcf files (Default: %s).\n", opt->vcf_opt->output_gl ? "yes" : "no");
	fprintf(fp, "\t--name STRING\n"
		"\t\tName of accession/individual/genotype; used in vcf header (Default: sample).\n");

//	fprintf(fp, "\t+++++++++++++++++++++++++++++++\n");
	fprintf(fp, "\nEstimation/Inference (optional):\n");
//	fprintf(fp, "\t+++++++++++++++++++++++++++++++\n");
	fprintf(fp, "\t--genotype_by_clustering\n"
		"\t\tGenotype by clustering (Default: no).\n");
	fprintf(fp, "\t\tRequires command-line arguments --amplici and --clustalo or --mafft.\n");
	fprintf(fp, "\t--clustalo FILE\n"
		"\t\tThe clustal omega executable (Default: %s).\n", opt->clustalo_command);
	fprintf(fp, "\t--mafft FILE\n"
		"\t\tThe mafft executable (Default: %s).\n", opt->mafft_command);
	fprintf(fp, "\t--alignment FILE1 FILE2\n"
		"\t\tAlignment input FILE1 and output FILE2 (Default: %s %s)\n",
					 opt->clustalo_infile, opt->clustalo_outfile);
	fprintf(fp, "\t--misalignment_rate FLOAT\n"
		"\t\tMaximum allowed subgenomic misalignment rate in [0, 0.5) (Default: %.2f).\n", opt->max_misalignment_proportion);
	fprintf(fp, "\t\tTolerated proportion of reads from one subgenome aligning to the other.\n");
//	fprintf(fp, "\t--sample INT\n"
//		"\t\tNumber of Monte Carlo samples (Default: %u)\n", opt->n_sample);

//	fprintf(fp, "\t++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	fprintf(fp, "\nScreening paralogs and other contaminants (optional):\n");
//	fprintf(fp, "\t++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	fprintf(fp, "\t-p, --drop INT\n"
		"\t\tDrop reads aligning to paralogs; -1 to automate (Default: %d)\n", opt->proptest_screen);
	fprintf(fp, "\t\tDrop INT of four most abundant haplotypes if specified.\n");
	fprintf(fp, "\t\tRequires command-line argument --amplici.\n");
	fprintf(fp, "\t--amplici EXE\n"
		"\t\tThe amplicon denoiser software (Default: none).\n");
	fprintf(fp, "\t\tWrites auxiliary files \"%s\" \"%s.fa\", and \"%s.out\".\n",
			opt->ac_fastq_file, opt->ac_outfile, opt->ac_outfile);
	fprintf(fp, "\t\tSee https://github.com/DormanLab/AmpliCI for more information.\n");
	fprintf(fp, "\t--amplici_fastq FILE\n"
		"\t\tSelected reads output to this FASTQ file for denoising (Default: %s).\n", opt->ac_fastq_file);
	fprintf(fp, "\t--write-fastq [FILE]\n"
		"\t\tWrite fastq file for AmpliCI and quit (Default: %s).\n", opt->write_fastq_and_quit ? "yes" : "no");
	fprintf(fp, "\t\tSee --amplici_fastq to name the file.\n");
	fprintf(fp, "\t\tWith optional argument write named fastq file after AmpliCI paralog filtering (Default: %s).\n", opt->fastq_after_filter ? opt->fastq_after_filter : "none");
	fprintf(fp, "\t--write-sam [FILE1 FILE2]\n"
		"\t\tWrite SAM files after AmpliCI paralog filtering (Default: %s, %s).\n", opt->sam_after_filter[0] ? opt->sam_after_filter[0] : "none", opt->sam_after_filter[1] ? opt->sam_after_filter[1] : "none");
	fprintf(fp, "\t--amplici_output BASENAME\n"
		"\t\tAmplicon denoiser output files (Default: %s.fa and %s.out).\n", opt->ac_outfile, opt->ac_outfile);
	fprintf(fp, "\t--amplici_low_bound FLOAT\n"
		"\t\tAmplicon denoiser lower abundance bound (Default: %f).\n", opt->ac_low_bound);
	fprintf(fp, "\t--amplici_threshold FLOAT\n"
		"\t\tAmplicon denoiser log likelihood threshold screens outlier reads (Default: %f).\n", opt->ac_ll_threshold);
	fprintf(fp, "\t--amplici_min_freq INT\n"
		"\t\tKeep haplotypes using unique reads frequency if most common > INT times (Default: %u).\n", opt->proptest_min_freq);
	fprintf(fp, "\t\tIf 0, use inferred cluster sizes to select kept haplotypes.\n");

//	fprintf(fp, "\t+++++++++++++++++++++++++++++++++++++++++++\n");
	fprintf(fp, "\nScreening reads, coverage checks (optional):\n");
//	fprintf(fp, "\t+++++++++++++++++++++++++++++++++++++++++++\n");
	fprintf(fp, "\t--expected_errors FLOAT\n"
		"\t\tDiscard reads with more than FLOAT expected errors (Default: %f).\n", opt->max_eerr);
	fprintf(fp, "\t--indel INT\n"
		"\t\tDrop reads whose alignments have more than INT indels (Default: %d)\n", opt->drop_indel);
	fprintf(fp, "\t--loglik FLOAT\n"
		"\t\tDrop reads with log likelihood < FLOAT (Default: %f).\n", opt->min_log_likelihood);
	fprintf(fp, "\t--min_posterior FLOAT\n"
		"\t\tDrop reads with maximum posterior alignment probability < FLOAT (Default: %f).\n", opt->min_posterior_alignment_prob);
	fprintf(fp, "\t--min INT\n"
		"\t\tDrop reads shorter than INT (Default: %d)\n", opt->min_length);
	fprintf(fp, "\t--max INT\n"
		"\t\tDrop reads longer than INT (Default: %d)\n", opt->max_length);
	fprintf(fp, "\t--secondary\n"
		"\t\tDrop secondary alignments (Default: %s)\n", opt->drop_secondary ? "yes" : "no");
	fprintf(fp, "\t--soft-clipped INT\n"
		"\t\tDrop reads where either alignment is clipped by > INT nucleotides (Default: %d)\n", opt->drop_soft_clipped);
	fprintf(fp, "\t--unmapped\n"
		"\t\tDrop reads unmapped in either alignment (Default: %s)\n", opt->drop_unmapped ? "yes" : "no");
	fprintf(fp, "\t--min_subgenomic_coverage FLOAT\n"
		"\t\tDo not genotype if expected subgenomic coverage < FLOAT (Default: %.1f).\n", opt->min_expected_coverage);
	fprintf(fp, "\t--coverage FLOAT\n"
		"\t\tSkip site if local subgenomic coverage < 100*FLOAT%% of expected subgenomic coverage (Default: %.1f).\n", opt->coverage_screen);
	fprintf(fp, "\t\tUse to skip subgenomic indel variants; not applicable if --genotype_by_clustering.\n");
	fprintf(fp, "\t--biallelic FLOAT\n"
		"\t\tSkip site if third allele >100*FLOAT%% of smaller subgenomic coverage (Default: %.1f).\n", opt->biallelic_screen);
	fprintf(fp, "\t\tUse --genotype_by_clustering to handle >2 alleles.\n");
	fprintf(fp, "\t--equal_coverage_test [FLOAT1 FLOAT2]\n"
		"\t\tPost-hoc test of equal allelic coverage (Default: %s, %.2f, %.2f).\n",
		opt->posthoc_coverage_test ? "yes" : "no", opt->phc_min_alignment_pp, opt->phc_min_genotype_pp);
	fprintf(fp, "\t\tOptional arguments are:\n"
		"\t\tFLOAT1: minimum posterior probability to assign read to subgenome\n"
		"\t\tFLOAT2: minimum posterior probability of heterozygosity to perform test\n");
	fprintf(fp, "\t\tNot applicable if --genotype_by_clustering.\n");

//	fprintf(fp, "\t++++++++++++++++++++++++++++++\n");
	fprintf(fp, "\nError recalibration (optional):\n");
//	fprintf(fp, "\t++++++++++++++++++++++++++++++\n");
	fprintf(fp, "\t--error_file FILE\n"
		"\t\tFile with error rates estimates (Default: none).\n");

	fprintf(fp, "\t--error_data FILE\n"
		"\t\tWrite observed errors to this file (Default: none).\n");

	fprintf(fp, "\nDebugging (optional):\n");
	fprintf(fp, "\t--site INT\n"
		"\t\tExtra debugging output for this site (Default: none).\n");

	fprintf(fp, "\n\nMore information: https://github.com/DormanLab/...\n");
} /* fprint_usage */
