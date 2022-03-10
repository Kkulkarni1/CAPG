//
//  pick_reads.c
//  pick reads aligned to each targeted region
//
//  Created by Yudi Zhang on 3/24/20.
//  Copyright © 2020 Yudi Zhang. All rights reserved.
//

#include <stdlib.h>
#include <string.h>

#include "fastq.h"
#include "nuc.h"
#include "qual.h"
#include "uthash.h"
#include "io.h"
#include "cmdline.h"
#include "error.h"
#include "pick_reads.h"
#include "order.h"

int match_indels(merge_hash *mh, sam **sds, ref_info *rfi);

int soft_clip_alignment(sam_entry *se, unsigned int fiveprime_sc, unsigned int threeprime_sc, int rc);
int match_indel(merge_hash *me, sam **sds, ref_info *rfi, unsigned int start_rd_idx, unsigned int end_rd_idx, unsigned int sg_idx);


void default_options_rf(options_rf *opt)
{
	opt->sam_file = NULL;
	opt->filter_unmapped = 1;
	opt->delim_ref = ':';
	opt->delim_len = '-';
	opt->samtools_command = "samtools";
	for (int i = 0; i < N_FILES; ++i) {
		opt->rsam_files[i] = NULL;
		opt->fsa_files[i] = NULL;
		opt->extracted_rf[i] = NULL;
	}
	opt->fastq_file = NULL;
	
} /* default_options_rf */

/**
 * Read sam file containing alignment(s) of homoeologous region(s) of 
 * subgenomic references. Record names of aligned references, start
 * and end locations of regions along genome (0-based), start and end
 * locations of aligned region (1-based). These are stored in ref
 * entry objects. Original sam_entry objects are also kept. Both in
 * ref_info object.
 *
 *
 * @param opt		reference options object
 * @param ref_in	reference information object, to be created
 * @param target_names	target names
 * @return		error status
 */
int make_targets_info(options_rf *opt, ref_info **ref_in,
	char const **target_names)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	sam *sd_ref  = NULL;
	ref_info *rfi = NULL;
	
	FILE *fp = fopen(opt->sam_file, "r");

	if (!fp)
		exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt->sam_file));

	read_sam(fp, &sd_ref, 0, 0);	/* assumes XY_ENCODING */
	fclose(fp);

	*ref_in = malloc(sizeof(**ref_in));
	if (!*ref_in)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "ref_in");
	rfi = *ref_in;

	for (size_t j = 0; j < N_FILES; ++j) {
		rfi->alignment_start[j] = SIZE_MAX;
		rfi->alignment_end[j] = 0;
		rfi->read_to_ref[j] = NULL;
	}
	rfi->read_len = 0;

	/* get the length of each reference name */
	int found = 0;

	//store the name[include chrosome] and starting and ending positions in the targeted A, B genome
	// i.e. split aradu.V14167.gnm2.chr02:696301-697301
	for (size_t i = 0; i < sd_ref->n_se; ++i) {
		sam_entry *se = &sd_ref->se[i];
		char const *rname = sd_ref->ref_names[se->ref];

		/* filter unmapped */
		if (opt->filter_unmapped && se->flag >> 2 & 1)
			continue;

		if (strcmp(se->name, target_names[1]))
			continue;

		if (found)
			exit(mmessage(ERROR_MSG, INVALID_USER_INPUT, "Target "
				"'%s' has two primary alignments in '%s'\n",
				target_names[1], opt->sam_file));
		found = 1;

		rfi->rf_idx = i;
		rfi->name_A = NULL;
		rfi->name_B = NULL;
		
		/* strand of B, in order to find the reads */
		if ((se->flag & 16) == 0)
			rfi->strand_B = 0;
		else
			rfi->strand_B = 1;

		/* parse subgenome B csome name, start and end */
		unsigned int idx_B = 0;
		while (idx_B < strlen(se->name) && se->name[idx_B] != opt->delim_ref)
			++idx_B;
		if (!idx_B || idx_B == strlen(se->name) || idx_B == strlen(se->name) - 1)
			exit(mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"Names are wrong format in SAM file of "
				"subgenomic alignments!"));
		rfi->name_B = malloc(idx_B + 1);
		strncpy(rfi->name_B, se->name, idx_B);
		rfi->name_B[idx_B] = '\0';
		if (opt->legacy_region_specification)			/* command-line provides 0-based, inclusive */
			rfi->start_B = atoi(&se->name[++idx_B]);	/* 0-based, inclusive */
		else							/* command-line provides 1-based, incclusive */
			rfi->start_B = atoi(&se->name[++idx_B]) - 1;	/* 0-based, inclusive */
		while (idx_B < strlen(se->name) && se->name[idx_B] != opt->delim_len)
			++idx_B;
		if (++idx_B >= strlen(se->name))
			exit(mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"Names are wrong format in SAM file of "
				"subgenomic alignments!"));
		rfi->end_B = atoi(&se->name[idx_B]);		/* 0-based, exclusive or 1-based inclusive */

		debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
			  "name: %s start %zu end %zu\n", rfi->name_B,
						rfi->start_B, rfi->end_B);

		/* parse subgenome A csome name, start and end */
		unsigned int idx_A = 0;
		while (idx_A < strlen(rname) && rname[idx_A] != opt->delim_ref)
			++idx_A;
		if (!idx_A || idx_A == strlen(se->name) || idx_A == strlen(se->name) - 1)
			exit(mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"Names are wrong format in SAM file of "
				"subgenomic alignments!"));
		rfi->name_A = malloc(idx_A + 1);
		strncpy(rfi->name_A, rname, idx_A);
		rfi->name_A[idx_A] = '\0';
		if (opt->legacy_region_specification)			/* command-line provides 0-based, inclusive */
			rfi->start_A = atoi(&rname[++idx_A]);		/* 0-based, inclusive */
		else							/* command-line provides 1-based, inclusive */
			rfi->start_A = atoi(&rname[++idx_A]) - 1;	/* 0-based, inclusive */
		while (idx_A < strlen(rname) && rname[idx_A] != opt->delim_len)
			++idx_A;
		if (++idx_A >= strlen(rname))
			exit(mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"Names are wrong format in SAM file of "
				"subgenomic alignments!"));
		rfi->end_A = atoi(&rname[idx_A]);	/* 0-based, exclusive */

		debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
			  "name: %s start %zu end %zu\n", rfi->name_A,
						rfi->start_A, rfi->end_B);

	}

	if (!found)
		exit(mmessage(ERROR_MSG, INVALID_USER_INPUT, "Target '%s' not "
			"found aligned in '%s' SAM file.\n", target_names[1],
			opt->sam_file));

	rfi->ref_sam = sd_ref;

	return NO_ERROR;
}/* make_targets_info */

/**
 * Exclude all reads not aligned to target regions and count the ones left in
 * sam::n_per_ref.
 *
 * The reads are initially aligned to the whole genome to reduce paralogous
 * read alignment, so the reference name in the sam file will just be a
 * chromosome name. Thus, to find reads aligned to the target, they must both
 * match the target chromosome and traverse the target region.
 *
 * @param ref_info	information about homoeologous aligned reference regions
 * @param opt		options about homeologous reference regions
 * @param sds		sam objects of reads aligned to each subgenome
 * @param csome_names	names of chromosomes containing target region
 * @return		error status
 */
int pickreads(ref_info *rfi, sam **sds, char const **csome_names)
{
	//int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	unsigned int j, m;
	size_t n_unmapped = 0, n_not_target = 0;


	/* find reads aligning to desired target */
	for (j = 0; j < N_FILES; ++j) {
		//sds[j]->n_per_ref = calloc(1, sizeof(*sds[j]->n_per_ref));
		//sds[j]->ref_list = malloc(sizeof(*sds[j]->ref_list));

		//if (!sds[j]->n_per_ref)
		//	return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "n_per_ref");

		for (m = 0; m < sds[j]->n_se; ++m) {
			sam_entry *se = &sds[j]->se[m];
			char const *rname = NULL;
			cigar *cig = se->cig;

			/* skip unmapped reads */
			if (se->flag & (1 << 2)) {
				++n_unmapped;
				se->exclude = 1;
				continue;
			}

			rname = sds[j]->ref_names[se->ref];
//fprintf(stderr, "se->ref = %u\n", se->ref);
//fprintf(stderr, "rname = %s\n", sds[j]->ref_names[se->ref]);

			/* skip reads not mapping to target */
			if (strcmp(rname, csome_names[j])) {
				++n_not_target;
				se->exclude = 1;
				/* there could be many such reads! */
				/*mmessage(INFO_MSG, NO_ERROR, "Read %s (%u) "
					"excluded because it does not align to "
					"target.\n", se->name, m);*/
				continue;
			}

			/* this read maps to this reference csome: allocate and assign a reference name for the target (NOT MEMORY EFFICIENT) */
			if (!strcmp(rname, csome_names[j])) {
				size_t start_pos[N_FILES] = {rfi->start_A, rfi->start_B};	/* 0 based, inclusive (from chrom_name:start-end in --ref_names and --geno file) */
				size_t end_pos[N_FILES] = {rfi->end_A, rfi->end_B};	/* 0-based, exclusive */
				size_t rf_index_s = se->pos - 1;			/* 0-based, inclusive */
				size_t rf_index_e = rf_index_s + cig->length_rf;	/* 0-based, exclusive */

//				printf("%zu %zu || %zu %zu\n", rf_index_s, rf_index_e, start_pos[j], end_pos[j]);

				/* and it maps within the target region */
				if ((rf_index_s >= start_pos[j] && rf_index_e <= end_pos[j]) ||	/* read contained within target */
					(rf_index_s <= start_pos[j] && rf_index_e > start_pos[j]) ||	/* read crosses 5' end of target */
					(rf_index_s < end_pos[j] && rf_index_e >= end_pos[j])) {	/* read crosses 3' end of target */

					//++sds[j]->n_per_ref[0];
					/* psych out hash_sam so correct read pairs hashed together */
					if (j && rfi->strand_B) {
						size_t lidx = strlen(se->name)-1;
						se->name[lidx] = se->name[lidx] == '+' ? '-' : '+';
					}
					se->which_ref = 0;
				}  else {
					++n_not_target;
					se->exclude = 1;
				}
			}
					/* there could be many such reads! */
					/*mmessage(INFO_MSG, NO_ERROR, "Read %s (%u) "
						"excluded because it does not align to "
						"target.\n", se->name, m);*/
				/*
					size_t length = strlen(rname) + 3
						+ (int)(log10(end_pos[j]) + 1) + (int)(log10(start_pos[j]) + 1);

					se->ref_name = malloc(length);
					if (!se->ref_name)
						return mmessage(ERROR_MSG, 
							MEMORY_ALLOCATION, "ref_name");
					sprintf(se->ref_name, "%s%s%zu%s%zu",
						rname, opt->delim_ref,
						start_pos[j], opt->delim_len,
						end_pos[j]);
					se->which_ref = rfi->ref_idx;
//					debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "REF_ID: %d REF_NAME: %s \n", se->which_ref, se->ref_name);
				*/
		}
	}

	mmessage(INFO_MSG, NO_ERROR, "Number of unmapped reads: %zu\n", n_unmapped);
	mmessage(INFO_MSG, NO_ERROR, "Number of alignments not mapping to target: %zu\n", n_not_target);

	return NO_ERROR;
}/* pickreads */

/**
 * Extract selected regions from fasta file. Samtools uses region specification
 * chr:from-to, where from and to are 1-based inclusive.
 *
 * @param samtools_command	samtools executable
 * @param ref_name		name of target region
 * @param ref_start		start of region, 0-based, inclusive
 * @param ref_end		end of region, 0-based, exclusive
 * @param ref_file		fasta file to extract regions from
 * @param ext_rf		fasta file to write with chosen region
 * @return			error status
 */
int extract_ref(char const *samtools_command, char const *ref_name, size_t ref_start,
		size_t ref_end, char const *ref_file, char const *ext_rf, options_rf *opt)
{
	FILE *fp;

	// index the whole reference genome file
	if (!ref_file || !(fp = fopen(ref_file, "r")))
		return(mmessage(ERROR_MSG, FILE_NOT_FOUND, ref_file));
	fclose(fp);
	
	size_t cmd_len = strlen(samtools_command) + strlen(ref_file)
		+ strlen(ref_name) + 2
		+ (int)(log10(ref_start + 1) + 1) + (int)(log10(ref_end) + 1)
		+ strlen(" faidx   -o ") + strlen(ext_rf) + 1;

	mmessage(INFO_MSG, NO_ERROR, "Length of command: %u\n", cmd_len);

	char *command = malloc(cmd_len * sizeof *command);
	if (!command)
		return(mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"samtools command"));
    
	sprintf(command, "%s faidx %s %s%c%zu%c%zu -o %s", samtools_command,
		ref_file, ref_name, opt->delim_ref, ref_start + 1,
		opt->delim_len, ref_end, ext_rf);
	
	mmessage(INFO_MSG, NO_ERROR, "Running samtools: '%s'\n", command);

	system(command);

	free(command);

	return NO_ERROR;
	
}/* extract_ref */

/**
 * Read command line options have to do with reference information.
 *
 * @param opt	options object for reference information
 * @param argc	number of arguments on command-line
 * @param argv	the arguments on the command-line
 * @return	error status
 */
int parse_rf_options(options_rf *opt, int argc, char *argv[])
{
	int argv_idx = 1;
	
	while (argv_idx < argc) {
		int j = 0;
		int len = strlen(argv[argv_idx]);
		while (j < len && argv[argv_idx][j] == '-')
			++j;
		char c = argv[argv_idx][j];
		switch (c) {
			case 'u':
				opt->filter_unmapped = 0;
				break;
			case 'l':
				opt->delim_len = argv[++argv_idx][0];
				break;
			case 'r':
				opt->delim_ref = argv[++argv_idx][0];
				break;
			case 'f':
				opt->sam_file = argv[++argv_idx];
				break;
			case 'b':
				for (j = 0; j < N_FILES; ++j) {
					opt->rsam_files[j] = argv[++argv_idx];
					fprintf(stderr, " %s", opt->rsam_files[j]);
				}
				fprintf(stderr, "\n");
				break;
			case 'd':
				for (j = 0; j < N_FILES; ++j) {
					opt->fsa_files[j] = argv[++argv_idx];
					fprintf(stderr, " %s", opt->fsa_files[j]);
				}
				fprintf(stderr, "\n");
				break;
			case 's':
				opt->samtools_command = argv[++argv_idx];
				mmessage(INFO_MSG, NO_ERROR, "Samtools "
					 "command: '%s'\n",
					 opt->samtools_command);
				break;
			case 'h':
			default:
				if (c != 'h')
					mmessage(ERROR_MSG, INVALID_CMD_OPTION,
						 argv[argv_idx]);
				usage_error((const char **)argv, argv_idx, opt);
				return 1;
		}
		++argv_idx;
	}
	return 0;
} /* parse_rf_options */

/**
 * Output reads aligned to chosen target, including reads not aligned to both.
 *
 * @param f	output FASTA file
 * @param sds	sam records of reads
 * @param mh	merged sam records
 */
void output_selected_reads(char const *f, sam **sds, merge_hash *mh)
{
	FILE *fpp = fopen(f, "w");

	if (!fpp)
		exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR, f));

	for (merge_hash *me = mh; me != NULL; me = me->hh.next) {
		sam_entry *se;

		if (me->nfiles != N_FILES) {
			if (me->indices[0])
				se = &sds[0]->se[me->indices[0][0]];
			else
				se = &sds[1]->se[me->indices[1][0]];
		} else {
			se = &sds[0]->se[me->indices[0][0]];
		}
		
		fprintf(fpp, "@%s\n", se->name);
		fwrite_nuc_segment(fpp, se->read, XY_ENCODING, 0,
				   se->read->len);
		fprintf(fpp, "\n+\n");
		fwrite_qual_sequence(fpp, se->qual);
		fprintf(fpp, "\n");
		continue;
	}
	fclose(fpp);
} /* output_selected_reads */

/**
 * Compute mapping from 0-based read index to 0-based reference index relative
 * to the target region (not alignment region) by alignment in SAM file. Maps
 * to -1 if maps to reference nucleotide not within aligned target. Resulting
 * map stored in ref_info::read_to_ref[N_FILES].
 *
 * @param rfi	reference information
 * @param sds	read alignments
 * @param me	merged hash of read alignments
 * @return	error status
 */
int index_read_to_ref(ref_info *rfi, sam *sds[N_FILES], merge_hash *me)
{
	sam_entry *se = &sds[0]->se[me->indices[0][0]];

	/* (re)allocate space for mapping */
	if (rfi->read_len != se->read->len) {
		for (unsigned int j = 0; j < N_FILES; ++j) {
			if (rfi->read_to_ref[j])
				free(rfi->read_to_ref[j]);
			rfi->read_to_ref[j] = malloc(se->read->len * sizeof(*rfi->read_to_ref[j]));
			if (!rfi->read_to_ref[j])
				exit(mmessage(ERROR_MSG, MEMORY_ALLOCATION, "ref_info::read_to_ref"));
		}
		rfi->read_len = se->read->len;
	}

	/* map read index to target index */
	for (unsigned int j = 0; j < N_FILES; ++j) {
		if (j)
			se = &sds[j]->se[me->indices[j][0]];
		//print_cigar(stderr, se);
		//fprintf(stderr, "Alignment %u (%zu):", j, se->read->len);

		size_t rf_idx = se->pos - 1;
		size_t rd_idx = 0;

		for (unsigned int i = 0; i < se->cig->n_ashes; ++i) {
			if (se->cig->ashes[i].type == CIGAR_DELETION) {
				rf_idx += se->cig->ashes[i].len;
				continue;
			} else if (se->cig->ashes[i].type == CIGAR_SOFT_CLIP) {
				for (unsigned int k = 0; k < se->cig->ashes[i].len; ++k) {
					rfi->read_to_ref[j][rd_idx + k] = ALIGN_SOFT_CLIP;	/* nowhere, but don't count */
					//fprintf(stderr, " %zu=%d", rd_idx + k, rfi->read_to_ref[j][rd_idx + k]);
				}
				rd_idx += se->cig->ashes[i].len;
				continue;
			} else if (se->cig->ashes[i].type == CIGAR_INSERTION) {
				for (unsigned int k = 0; k < se->cig->ashes[i].len; ++k) {
					rfi->read_to_ref[j][rd_idx + k] = ALIGN_INSERTION;	/* nowhere */
					//fprintf(stderr, " %zu=%d", rd_idx + k, rfi->read_to_ref[j][rd_idx + k]);
				}
				rd_idx += se->cig->ashes[i].len;
				continue;
			} else if (se->cig->ashes[i].type == CIGAR_HARD_CLIP) {
				continue;
			} else if (se->cig->ashes[i].type != CIGAR_MATCH
				   && se->cig->ashes[i].type != CIGAR_MISMATCH
				   && se->cig->ashes[i].type != CIGAR_MMATCH) {
				continue;
			}

			for (unsigned int k = 0; k < se->cig->ashes[i].len; ++k) {
				/* position maps to target region */
				if (!j && rf_idx + k >= (j ? rfi->start_B : rfi->start_A) && rf_idx + k < (j ? rfi->end_B : rfi->end_A))
					rfi->read_to_ref[j][rd_idx + k] = rf_idx + k - (j ? rfi->start_B : rfi->start_A);
				else if (!j)	/* outside target */
					rfi->read_to_ref[j][rd_idx + k] = ALIGN_OUTSIDE;
				else if (j && rf_idx + k >= (j ? rfi->start_B : rfi->start_B) && rf_idx + k < (j ? rfi->end_B : rfi->end_A))
					rfi->read_to_ref[j][rd_idx + k] = rf_idx + k - (j ? rfi->start_B : rfi->start_A);
				else
					rfi->read_to_ref[j][rd_idx + k] = ALIGN_OUTSIDE;
				//fprintf(stderr, " %zu=%d", rd_idx + k, rfi->read_to_ref[j][rd_idx + k]);
			}
			rd_idx += se->cig->ashes[i].len;
			rf_idx += se->cig->ashes[i].len;
		}
		//fprintf(stderr, "\n");
	}

	return NO_ERROR;
} /* index_read_to_ref */

/**
 * Match indels between read pairs to maximize homoeologous alignments.
 *
 * [TODO,BUG,KSD] This whole function assumes N_FILES == 2
 *
 * @param mh	merge hash
 * @param sds	sam files of read alignments
 * @param rfi	reference information
 * @return	error status
 */
int match_indels(merge_hash *mh, sam **sds, ref_info *rfi)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	mlogit_stuff mls = {NULL, 0};

	for (merge_hash *me = mh; me != NULL; me = me->hh.next) {

		if (me->exclude)
			continue;

		index_read_to_ref(rfi, sds, me);

		sam_entry *se = &sds[0]->se[me->indices[0][0]];	/* sgA is default alignment */
		unsigned int rd_idx = 0;
		unsigned int last_rd_idx = 0;	/* last read index involved in homoeologous alignment */
		unsigned int last_rf_idx, last_other_rf_idx;		/* last sgA and sgB index involved in sgA (default) alignment */
		unsigned int last_alt_rf_idx, last_alt_other_rf_idx;	/* last sgA and sgB index involved in sgB (alternate) alignment */
		unsigned char in_discrepancy = 0;
		double lla = 0, llb = 0;

		for (unsigned int i = 0; i < se->cig->n_ashes; ++i) {
			if (se->cig->ashes[i].type == CIGAR_SKIP) {
				return(mmessage(ERROR_MSG, INTERNAL_ERROR, "cannot handle cigar skips"));
			} else if (se->cig->ashes[i].type == CIGAR_SOFT_CLIP) {
				/* these have already been matched across read alignments */
				rd_idx += se->cig->ashes[i].len;
				continue;
			} else if (se->cig->ashes[i].type == CIGAR_DELETION) {
				if (in_discrepancy) {
					/* what to do? */
					for (unsigned int j = 0; j < se->cig->ashes[i].len; ++j) {
						mls.pos = rd_idx + j;
						/* force a mismatch to penalize deletions worse than substitution */
						lla += sub_prob_given_q_with_encoding(IUPAC_A,
							//get_nuc(se->read, XY_ENCODING, rd_idx + j),
							XY_C,	/* force mismatch */
							IUPAC_ENCODING, XY_ENCODING,
							MAX_ILLUMINA_QUALITY_SCORE, 1, (void *) &mls);
							//get_qual(se->qual, rd_idx + j), 1, (void *) &mls);
					}
				}
				continue;
			} else if (se->cig->ashes[i].type != CIGAR_MMATCH
				&& se->cig->ashes[i].type != CIGAR_MATCH
				&& se->cig->ashes[i].type != CIGAR_MISMATCH
				&& se->cig->ashes[i].type != CIGAR_INSERTION) {
				continue;
			}

			for (unsigned int j = 0; j < se->cig->ashes[i].len; ++j) {
				int rf_idx = rfi->read_to_ref[0][rd_idx + j];		/* sgA index in A alignment */
				int other_rf_idx = rfi->read_to_ref[1][rd_idx + j];	/* sgB index in B alignment */
				int alt_rf_idx = other_rf_idx >= 0
					? rfi->map_B_to_A[other_rf_idx]	/* sgA index implied by B alignment */
					: other_rf_idx;
				int alt_other_rf_idx = rf_idx >= 0
					? rfi->map_A_to_B[rf_idx]	/* sgB index implied by A alignment */
					: rf_idx;

				/* exit discrepancy */
				if (rf_idx >= 0 && rf_idx == alt_rf_idx) {
					if (other_rf_idx != alt_other_rf_idx)
						exit(mmessage(ERROR_MSG, INTERNAL_ERROR,
							"Homoeologous alignment in sgA, but not sgB!\n"));
					if (in_discrepancy && rd_idx + j > last_rd_idx + 1) {
						if (lla > llb) {
							match_indel(me, sds, rfi, last_rd_idx, rd_idx, 0);
						} else if (llb > lla) {
							match_indel(me, sds, rfi, last_rd_idx, rd_idx, 1);
						} else {
							double r = (double) rand() / RAND_MAX;
	
							if (r < 0.5)
								match_indel(me, sds, rfi, last_rd_idx, rd_idx, 0);
							else
								match_indel(me, sds, rfi, last_rd_idx, rd_idx, 1);
						}
						lla = llb = 0;
					}
					in_discrepancy = 1;
					last_rd_idx = rd_idx;	/* last read index involved in homoeologous alignment */
					last_rf_idx = rf_idx;				/* last sgA index aligned in sgA alignment */
					last_other_rf_idx = other_rf_idx;		/* last sgB index aligned in sgB alignment */
					last_alt_rf_idx = alt_rf_idx;			/* last sgA index aligned in sgB alignment */
					last_alt_other_rf_idx = alt_other_rf_idx;	/* last sgB index aligned in sgA alignment */
					

				/* continue in discrepancy */
				} else if (in_discrepancy && (rf_idx < 0 || rf_idx != alt_rf_idx)) {

					iupac_t rn;

					mls.pos = rd_idx + j;

					/* log likelihood of sgA alignment */
					/* read to sgA */
					if (rf_idx >= 0) {
						/* there has been deletion */
						if ((unsigned int) rf_idx > last_rf_idx + 1) {
							lla += (rf_idx - last_rf_idx - 1) * sub_prob_given_q_with_encoding(
								IUPAC_A, XY_C,
								IUPAC_ENCODING, XY_ENCODING,
								MAX_ILLUMINA_QUALITY_SCORE, 1, (void *) &mls);
						}
						rn = rfi->ref[0][(size_t) rf_idx + rfi->start_A - rfi->alignment_start[0]];
						lla += sub_prob_given_q_with_encoding(rn,
							get_nuc(se->read, XY_ENCODING, rd_idx + j),
							IUPAC_ENCODING, XY_ENCODING,
							get_qual(se->qual, rd_idx + j), 1, (void *) &mls);
						last_rf_idx = rf_idx;
					} else {	/* CIGAR_INSERTION */
						/* think I decided not to penalize insertion
						lla += sub_prob_given_q_with_encoding(IUPAC_A,
							XY_C,
							IUPAC_ENCODING, XY_ENCODING,
							MAX_ILLUMINA_QUALITY_SCORE, 1, (void *) &mls);
						*/
					}
					/* read to sgB via sgA alignment */
					if (alt_other_rf_idx >= 0) {
						if ((unsigned int) alt_other_rf_idx > last_alt_other_rf_idx + 1) {
							lla += (alt_other_rf_idx - last_other_rf_idx - 1) * sub_prob_given_q_with_encoding(
								IUPAC_A, XY_C,
								IUPAC_ENCODING, XY_ENCODING,
								MAX_ILLUMINA_QUALITY_SCORE, 1, (void *) &mls);
						}
						rn = rfi->ref[1][(size_t) alt_other_rf_idx + rfi->start_B - rfi->alignment_start[1]];
						lla += sub_prob_given_q_with_encoding(rn,
							get_nuc(se->read, XY_ENCODING, rd_idx + j),
							IUPAC_ENCODING, XY_ENCODING,
							get_qual(se->qual, rd_idx + j), 1, (void *) &mls);
						last_alt_other_rf_idx = alt_other_rf_idx;
					} else {
						/* insertion */
					}

					/* log likelihood of sgB alignment */
					/* read to sgA via sgB alignment */
					if (alt_rf_idx >= 0) {
						if ((unsigned int) alt_rf_idx > last_alt_rf_idx + 1) {
							llb += (alt_rf_idx - last_alt_rf_idx - 1) * sub_prob_given_q_with_encoding(
								IUPAC_A, XY_C,
								IUPAC_ENCODING, XY_ENCODING,
								MAX_ILLUMINA_QUALITY_SCORE, 1, (void *) &mls);
						}
						rn = rfi->ref[0][(size_t) alt_rf_idx + rfi->start_A - rfi->alignment_start[0]];
						/* log likelihood of this alternative */
						llb += sub_prob_given_q_with_encoding(rn,
							get_nuc(se->read, XY_ENCODING, rd_idx + j),
							IUPAC_ENCODING, XY_ENCODING,
							get_qual(se->qual, rd_idx + j), 1, (void *) &mls);
						last_alt_rf_idx = alt_rf_idx;
					}
					/* read to sgB via sgB alignment */
					if (other_rf_idx >= 0) {
						if ((unsigned int) other_rf_idx > last_other_rf_idx + 1) {
							llb += (other_rf_idx - last_other_rf_idx - 1) * sub_prob_given_q_with_encoding(
								IUPAC_A, XY_C,
								IUPAC_ENCODING, XY_ENCODING,
								MAX_ILLUMINA_QUALITY_SCORE, 1, (void *) &mls);
						}
						/* read to sgB */
						rn = rfi->ref[1][(size_t) other_rf_idx + rfi->start_B - rfi->alignment_start[1]];
						/* log likelihood of this alternative */
						llb += sub_prob_given_q_with_encoding(rn,
							get_nuc(se->read, XY_ENCODING, rd_idx + j),
							IUPAC_ENCODING, XY_ENCODING,
							get_qual(se->qual, rd_idx + j), 1, (void *) &mls);
						last_other_rf_idx = other_rf_idx;
					}
				}
			}
			rd_idx += se->cig->ashes[i].len;
		}
	}

	return NO_ERROR;
} /* match_indels */

/**
 * Match indels within a discrepant region between two homoeologously
 * aligned positions.
 *
 * [TODO,KSD,BUG] Assumes N_FILES == 2.
 *
 * @param me		merge_hash with paired reads
 * @param sds		sam entries
 * @param rfi		reference information
 * @param start_rd_idx	index of read nucleotide of 5' homoeologous alignment
 * @param end_rd_idx	index of read nucleotide of 3' homoeologous alignment
 * @param sg_idx	which subgenome has the "better" alignment
 * @return		error status
 */
int match_indel(merge_hash *me, sam **sds, ref_info *rfi, unsigned int start_rd_idx, unsigned int end_rd_idx, unsigned int sg_idx)
{

	/* this alignment */
	sam_entry *se = &sds[!sg_idx]->se[me->indices[!sg_idx][0]];
	unsigned char in_region = 0;

	/* other subgenome reference index aligned to homoeologous site right before discrepant region */
	int other_rf_idx = rfi->read_to_ref[sg_idx][start_rd_idx];
	unsigned int last_other_rf_idx = (unsigned int) other_rf_idx;	/* safe cast */
	int prev_ash = ALIGN_MATCH;	/* previous homoeologous site is a match */
	unsigned int n_old_ashes = 0, n_new_ashes = 0;
	unsigned int rd_idx = 0;

	/* count number of old and new ashes in the discrepant region */
	for (unsigned int j = 0; j < se->cig->n_ashes; ++j) {
		if (se->cig->ashes[j].type == CIGAR_SOFT_CLIP) {
			rd_idx += se->cig->ashes[j].len;
			if (in_region) /* cannot be */
				exit(mmessage(ERROR_MSG, INTERNAL_ERROR,
					"Cannot enter a soft clip while in a discrepant region.\n"));
			continue;
		} else if (se->cig->ashes[j].type == CIGAR_DELETION) {
			if (in_region) ++n_old_ashes;
			continue;
		} else if (se->cig->ashes[j].type != CIGAR_MMATCH
			&& se->cig->ashes[j].type != CIGAR_MATCH
			&& se->cig->ashes[j].type != CIGAR_MISMATCH
			&& se->cig->ashes[j].type != CIGAR_INSERTION) {
			continue;
		}

		/* handle insertion or match */

		if (in_region)
			++n_old_ashes; /* presume part of discrepant region */

		/* walk through insertion or match */
		for (unsigned int k = 0; k < se->cig->ashes[j].len; ++k) {

			/* stepping out of region */
			if (in_region && rd_idx + k >= end_rd_idx) {
				if (!k)	/* actually ends with this new ash */
					--n_old_ashes;
				in_region = 0;
				break;

			/* continuing in region */
			} else if (in_region) {
				other_rf_idx = rfi->read_to_ref[!sg_idx][rd_idx + k];		/* other reference index aligned to this read nucleotide */
				int alt_rf_idx = other_rf_idx >= 0
					? sg_idx ? rfi->map_A_to_B[other_rf_idx]
						: rfi->map_B_to_A[other_rf_idx]	/* A ref index or ALIGN_INSERTION or ALIGN_DELETION */
					: ALIGN_INSERTION;

				/* there must have been a deletion in other subgenome alignment */
				if (other_rf_idx >= 0 && (unsigned int) other_rf_idx > last_other_rf_idx + 1) {
					++n_new_ashes;
					prev_ash = ALIGN_DELETION;
				}

				/* entering new match region */
				if (alt_rf_idx >= 0 && prev_ash != ALIGN_MATCH) {
					++n_new_ashes;
					prev_ash = ALIGN_MATCH;

				/* entering new ash */
				} else if (alt_rf_idx != prev_ash) {
					++n_new_ashes;
					prev_ash = alt_rf_idx;
				}

				/* record last consumed index of sgB reference */
				if (other_rf_idx >= 0)
					last_other_rf_idx = other_rf_idx;

			/* entering region: must enter homoeologous INDEL */
			} else if (!in_region && rd_idx + k > start_rd_idx) {
				other_rf_idx = rfi->read_to_ref[!sg_idx][rd_idx + k];
				int alt_rf_idx = other_rf_idx >= 0
					? sg_idx ? rfi->map_A_to_B[other_rf_idx]
						: rfi->map_B_to_A[other_rf_idx]	/* ALIGN_INSERTION or ALIGN_DELETION */
					: ALIGN_INSERTION;

				/* there must have been a deletion in sgB alignment */
				if (other_rf_idx >= 0 && (unsigned int) other_rf_idx > last_other_rf_idx + 1) {
					++n_new_ashes;
					prev_ash = ALIGN_DELETION;
				}

				++n_new_ashes;
				prev_ash = alt_rf_idx;
				in_region = 1;

				/* record last consumed index of sgB reference */
				if (other_rf_idx >= 0)
					last_other_rf_idx = other_rf_idx;
			}
		}
		rd_idx += se->cig->ashes[j].len;
		if (rd_idx >= end_rd_idx)
			break;
	}

	/* allocate new ashes */
	ash *ashes = se->cig->ashes;
	if (n_new_ashes != n_old_ashes) {
		ash *ashes = malloc((se->cig->n_ashes + n_new_ashes - n_old_ashes) * sizeof(*ashes));

		if (!ashes)
			return(mmessage(ERROR_MSG, MEMORY_ALLOCATION, "ashes"));
	}

	/* build the new cigar for this subgenome alignment based on other selected subgenome alignment */
	rd_idx = 0;	/* this is read nucleotide with homoeologous alignment right before the discrepant region */
	last_other_rf_idx = rfi->read_to_ref[sg_idx][start_rd_idx];
	prev_ash = ALIGN_MATCH;
	unsigned int n_ash = 0;
	for (unsigned int j = 0; j < se->cig->n_ashes; ++j) {
		if (se->cig->ashes[j].type == CIGAR_SOFT_CLIP) {
			ashes[n_ash].type = CIGAR_SOFT_CLIP;
			ashes[n_ash++].len = se->cig->ashes[j].len;
			rd_idx += se->cig->ashes[j].len;
			continue;
		} else if (se->cig->ashes[j].type == CIGAR_DELETION) {
			ashes[n_ash].type = CIGAR_DELETION;
			ashes[n_ash++].len = se->cig->ashes[j].len;
			continue;
		} else if (se->cig->ashes[j].type != CIGAR_MMATCH
			&& se->cig->ashes[j].type != CIGAR_MATCH
			&& se->cig->ashes[j].type != CIGAR_MISMATCH
			&& se->cig->ashes[j].type != CIGAR_INSERTION) {
			continue;
		}

		/* insertion or match */
		for (unsigned int k = 0; k < se->cig->ashes[j].len; ++k) {

			/* stepping out of region */
			if (in_region && rd_idx + k >= end_rd_idx) {
				/* take rest of old ash */
				if (prev_ash != se->cig->ashes[j].type) {
					ashes[n_ash].type = se->cig->ashes[j].type;
					ashes[n_ash++].len = se->cig->ashes[j].len - k;
				}
				break;

			/* continuing in region */
			} else if (in_region) {
				other_rf_idx = rfi->read_to_ref[!sg_idx][rd_idx + k];
				int alt_rf_idx = other_rf_idx >= 0
					? sg_idx
						? rfi->map_A_to_B[other_rf_idx]	/* A ref index or ALIGN_INSERTION or ALIGN_DELETION */
						: rfi->map_B_to_A[other_rf_idx]	/* A ref index or ALIGN_INSERTION or ALIGN_DELETION */
					: ALIGN_INSERTION;	/* between homoeologous alignments, this is the only option for a non-aligned read nucleotide */

				/* there must have been a deletion in other subgenome alignment */
				if (other_rf_idx >= 0 && (unsigned int) other_rf_idx > last_other_rf_idx + 1) {
					ashes[n_ash].type = CIGAR_DELETION;
					ashes[n_ash++].len = other_rf_idx - last_other_rf_idx;
					prev_ash = ALIGN_DELETION;
				}

				/* entering new match region */
				if (alt_rf_idx >= 0 && prev_ash != ALIGN_MATCH) {
					ashes[n_ash].type = ALIGN_MATCH;
					ashes[n_ash].len = 1;
					prev_ash = ALIGN_MATCH;

				/* entering new ash */
				} else if (alt_rf_idx != prev_ash) {
					ashes[++n_ash].type = alt_rf_idx == ALIGN_INSERTION ? CIGAR_INSERTION: CIGAR_DELETION;
					ashes[n_ash].len = 1;
					prev_ash = alt_rf_idx;

				/* continuing in same ash */
				} else {
					++ashes[n_ash].len;
				}

				/* record last consumed index of sgB reference */
				if (other_rf_idx >= 0)
					last_other_rf_idx = other_rf_idx;

			/* entering region: must enter homoeologous INDEL, because
			 * if it was a match, we'd just have another homoeologous
			 * match
			 */
			} else if (!in_region && rd_idx + k > start_rd_idx) {
				other_rf_idx = rfi->read_to_ref[!sg_idx][rd_idx + k];
				int alt_rf_idx = other_rf_idx >= 0
					? sg_idx
						? rfi->map_A_to_B[other_rf_idx]
						: rfi->map_B_to_A[other_rf_idx]	/* ALIGN_INSERTION or ALIGN_DELETION */
					: ALIGN_INSERTION;

				/* there must have been a deletion in other subgenome alignment */
				if (other_rf_idx >= 0 && (unsigned int) other_rf_idx > last_other_rf_idx + 1) {
					ashes[n_ash].type = CIGAR_DELETION;
					ashes[n_ash++].len = other_rf_idx - last_other_rf_idx;
					prev_ash = ALIGN_DELETION;
				}

				/* entering new match region */
				if (alt_rf_idx >= 0 && prev_ash != ALIGN_MATCH) {
					ashes[n_ash].type = ALIGN_MATCH;
					ashes[n_ash].len = 1;
					prev_ash = ALIGN_MATCH;

				/* entering new ash */
				} else if (alt_rf_idx != prev_ash) {
					/* insertion in other subgenome with respect to this subgenome indicates a deletion in read wrt this subgenome alignment */
					ashes[n_ash].type = alt_rf_idx == ALIGN_INSERTION ? CIGAR_DELETION : CIGAR_INSERTION;
					ashes[n_ash++].len = 1;
				}

				prev_ash = alt_rf_idx;
				in_region = 1;

				/* record last consumed index of sgB reference */
				if (other_rf_idx >= 0)
					last_other_rf_idx = other_rf_idx;
			}
		}
		rd_idx += se->cig->ashes[j].len;
		if (rd_idx >= end_rd_idx)
			break;
	}
	if (n_new_ashes > n_old_ashes) {
		free(se->cig->ashes);
		se->cig->n_ashes += n_new_ashes - n_old_ashes;
	} else if (n_old_ashes > n_new_ashes) {
		free(se->cig->ashes);
		se->cig->n_ashes -= n_old_ashes - n_new_ashes;
	}
	se->cig->ashes = ashes;

	return NO_ERROR;
} /* match_indel */

/**
 * For a merged hash, reset alignments to have same level of soft-clipping of the
 * read in all alignments. This function may further exclude some reads if they
 * their aligned portion is completely soft-clipped.
 *
 *
 * @param mh		pointer to the merged hash
 * @param nalign	number of subgenomes (alignments per read)
 * @param sds		sam hashes, one per subgenome
 * @param b_rc		indicate if B subgenome reference is reverse
 *			complemented relative to A subgenome reference
 * @return		number of reads excluded by soft-clipping
 */
int match_soft_clipping(merge_hash *mh, unsigned int nalign, sam **sds,
	unsigned int b_rc)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	size_t n_reads_excluded = 0;

	debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "match_soft_clipping()\n");

	for (merge_hash *me = mh; me != NULL; me = me->hh.next) {
		unsigned int fiveprime_sc = 0, threeprime_sc = 0;

		if (me->exclude)
			continue;

		/* compute the maximum soft-clipping at each end */
		for (unsigned int j = 0; j < nalign; ++j) {
			sam_entry *se = &sds[j]->se[me->indices[j][0]];
			unsigned int lidx = se->cig->n_ashes - 1;

			if (j & b_rc) {
				if (se->cig->ashes[lidx].type == CIGAR_SOFT_CLIP
					&& fiveprime_sc < se->cig->ashes[lidx].len)
					fiveprime_sc = se->cig->ashes[lidx].len;
				if (se->cig->ashes[0].type == CIGAR_SOFT_CLIP
					&& threeprime_sc < se->cig->ashes[0].len)
					threeprime_sc = se->cig->ashes[0].len;
			} else {
				if (se->cig->ashes[0].type == CIGAR_SOFT_CLIP
					&& fiveprime_sc < se->cig->ashes[0].len)
					fiveprime_sc = se->cig->ashes[0].len;
				if (se->cig->ashes[lidx].type == CIGAR_SOFT_CLIP
					&& threeprime_sc < se->cig->ashes[lidx].len)
					threeprime_sc = se->cig->ashes[lidx].len;
			}
		}

		debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "Read %s 5' soft clip %u; 3' soft clip %u\n", sds[0]->se[me->indices[0][0]].name, fiveprime_sc, threeprime_sc);

		for (unsigned int j = 0; j < nalign; ++j) {
			sam_entry *se = &sds[j]->se[me->indices[j][0]];

			 if (soft_clip_alignment(se, fiveprime_sc,
			 			threeprime_sc, j && b_rc)) {
				++n_reads_excluded;
				mmessage(INFO_MSG, NO_ERROR, "Read %s excluded "
					"by soft-clipping entire alignment in "
					"subgenome %u.\n", se->name, j);
				me->exclude = 1;
				break;
			}
		}
	}

	return(n_reads_excluded);
} /* match_soft_clipping */

/**
 * Soft-clip alignment to match 5' and 3' soft clipping, which should always 
 * equal or exceed current soft clip.
 *
 * @param se		current alignment
 * @param fiveprime_sc	desired 5' soft clip
 * @param threeprime_sc	desired 3' soft clip
 * @param rc		current alignment is to subgenomic reference that is
 *			reverse complemented relative to first subgenomic
 *			reference
 * @return		will alignment be soft-clipped out of existence?
 */
int soft_clip_alignment(sam_entry *se, unsigned int fiveprime_sc,
					unsigned int threeprime_sc, int rc)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	unsigned int lidx = se->cig->n_ashes - 1;
	ash *new_ashes = NULL;
	unsigned int prime_sc = 0;
	unsigned int len = 0;
	int direction = 0;

	for (unsigned int k = 0; k < 2; ++k) {	/* 5' or 3' end */
		if (rc) {
			len = k 	/* 3' */
				? se->cig->ashes[0].type == CIGAR_SOFT_CLIP ? se->cig->ashes[0].len : 0
				: se->cig->ashes[lidx].type == CIGAR_SOFT_CLIP ? se->cig->ashes[lidx].len : 0;
			prime_sc = k ? threeprime_sc : fiveprime_sc;
			direction = k ? 0 : 1;
		} else {
			len = k 	/* 5' */
				? se->cig->ashes[lidx].type == CIGAR_SOFT_CLIP ? se->cig->ashes[lidx].len : 0
				: se->cig->ashes[0].type == CIGAR_SOFT_CLIP ? se->cig->ashes[0].len : 0;
			prime_sc = k ? threeprime_sc : fiveprime_sc;
			direction = k ? 1 : 0;
		}

		if (!direction && len < prime_sc) {	/* need to extend soft clip at left end */
			unsigned int idx = 0;
			int n_ashes = se->cig->ashes[0].type == CIGAR_SOFT_CLIP ? 0 : 1;

debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "Read %s with cigar ", se->name);
debug_call(fxn_debug >= DEBUG_I, fxn_debug, print_cigar(stderr, se));
debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, " to be soft-clip extended from %u to %u on left side\n", len, prime_sc);

			for (idx = 1 - n_ashes; idx < se->cig->n_ashes; ++idx) {
				/* soft-clip ashes that consume
				 * read until the desired soft
				 * clip length is achieved; adjust
				 * mapping position if ash consumes
				 * reference
				 */
				if (se->cig->ashes[idx].type == CIGAR_MATCH
					|| se->cig->ashes[idx].type == CIGAR_MMATCH
					|| se->cig->ashes[idx].type == CIGAR_MISMATCH
					|| se->cig->ashes[idx].type == CIGAR_INSERTION) {
					if (len + se->cig->ashes[idx].len > prime_sc) {
						se->cig->ashes[idx].len -= prime_sc - len;
						if (se->cig->ashes[idx].type != CIGAR_INSERTION) {
							se->pos += prime_sc - len;
							se->cig->length_rf -= prime_sc - len;
						}
						break;
					} else {	/* lose an entire ash */
						--n_ashes;
						len += se->cig->ashes[idx].len;
						if (se->cig->ashes[idx].type != CIGAR_INSERTION) {
							se->pos += se->cig->ashes[idx].len;
							se->cig->length_rf -= se->cig->ashes[idx].len;
						}
					}
				/* soft-clipping something consuming reference, but not read */
				} else if (se->cig->ashes[idx].type == CIGAR_DELETION
					|| se->cig->ashes[idx].type == CIGAR_SKIP) {
					--n_ashes;
					se->pos += se->cig->ashes[idx].len;
				} else if (se->cig->ashes[idx].type == CIGAR_SOFT_CLIP) {
					return 1;
				}
			}
			if (n_ashes)
				new_ashes = malloc((se->cig->n_ashes + n_ashes) * sizeof(*new_ashes));
			else
				new_ashes = se->cig->ashes;
			for (unsigned int i = idx; i < se->cig->n_ashes; ++i) { 
				new_ashes[i - idx + 1].type = se->cig->ashes[i].type;
				new_ashes[i - idx + 1].len = se->cig->ashes[i].len;
			}
			new_ashes[0].type = CIGAR_SOFT_CLIP;
			new_ashes[0].len = prime_sc;
			if (n_ashes) {
				free(se->cig->ashes);
				se->cig->ashes = new_ashes;
				se->cig->n_ashes += n_ashes;
				lidx += n_ashes;
			}

debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "New cigar ");
debug_call(fxn_debug >= DEBUG_I, fxn_debug, print_cigar(stderr, se));
debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, "\n");
		} else if (direction && len < prime_sc) {
			unsigned int idx = lidx;
			int n_ashes = se->cig->ashes[lidx].type == CIGAR_SOFT_CLIP ? 0 : 1;

debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "Read %s with cigar ", se->name);
debug_call(fxn_debug >= DEBUG_I, fxn_debug, print_cigar(stderr, se));
debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, " to be soft-clip extended from %u to %u on right side\n", len, prime_sc);

			for (idx = se->cig->n_ashes + n_ashes - 1; idx-- > 0; ) {
				if (se->cig->ashes[idx].type == CIGAR_MATCH
					|| se->cig->ashes[idx].type == CIGAR_MMATCH
					|| se->cig->ashes[idx].type == CIGAR_MISMATCH
					|| se->cig->ashes[idx].type == CIGAR_INSERTION) {
					if (len + se->cig->ashes[idx].len > prime_sc) {
						se->cig->ashes[idx].len -= prime_sc - len;
						if (se->cig->ashes[idx].type != CIGAR_INSERTION)
							se->cig->length_rf -= prime_sc - len;
						break;
					} else {
						--n_ashes;
						len += se->cig->ashes[idx].len;
						if (se->cig->ashes[idx].type != CIGAR_INSERTION)
							se->cig->length_rf -= se->cig->ashes[idx].len;
					}
				/* soft-clipping something consuming reference, but not read */
				} else if (se->cig->ashes[idx].type == CIGAR_DELETION
					|| se->cig->ashes[idx].type == CIGAR_SKIP) {
					--n_ashes;
				} else if (se->cig->ashes[idx].type == CIGAR_SOFT_CLIP) {
					return 1;
				}
			}
			if (n_ashes)
				new_ashes = malloc((se->cig->n_ashes + n_ashes) * sizeof(*new_ashes));
			else
				new_ashes = se->cig->ashes;
			for (unsigned int i = 0; i <= idx; ++i) {
				new_ashes[i].type = se->cig->ashes[i].type;
				new_ashes[i].len = se->cig->ashes[i].len;
			}
			se->cig->n_ashes += n_ashes;
			new_ashes[se->cig->n_ashes - 1].type = CIGAR_SOFT_CLIP;
			new_ashes[se->cig->n_ashes - 1].len = prime_sc;
			if (n_ashes) {
				free(se->cig->ashes);
				se->cig->ashes = new_ashes;
				lidx += n_ashes;
			}
debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "New cigar ");
debug_call(fxn_debug >= DEBUG_I, fxn_debug, print_cigar(stderr, se));
debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, "\n");
		}
	}

	return 0;
} /* soft_clip_alignment */


/**
 * For a merged hash, where multiple alignments per read are combined,
 * reset the alignments to cover the minimum covered region by soft-clipping
 * as necessary. If a reference is reverse complemented relative to the first
 * reference in reference alignments, reverse the cigar string temporarily
 * while soft-clipping. An alternative, more complicated solution is to
 * extend the alignment.
 *
 *
 * @param mh		pointer to the merged hash
 * @param nalign	number of alignments
 * @param sds		sam hashes, one per reference
 * @param start_pos	minimum start position of all reads per alignment
 *			(0-base reference index)
 * @param end_pos	maximum end position of all reads per alignment
 *			(1-base reference index) [it is start_pos + length of ref region]
 * @param B_strand	indicate if B subgenome reference is reverse
 *			complemented relative to A subgenome reference
 * [TODO,KSD] B_strand should be array unsigned int[N_FILES]
 * [question: should we consider the homolegous region of A and B here?]
 * [answer: yes, and it was a bug to ignore it]
 */
int match_extent(merge_hash *mh, unsigned int nalign, sam **sds,
		size_t *start_pos, size_t *end_pos, unsigned int B_strand)
{
	int fxn_debug = DEBUG_I;//ABSOLUTE_SILENCE;//
	size_t n_reads = 0;

	debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "Maximum extents:\n");
	for (unsigned int j = 0; j < nalign; ++j)
		debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, "\t%zu-%zu\n",
						start_pos[j], end_pos[j]);

	for (merge_hash *me = mh; me != NULL; me = me->hh.next, ++n_reads) {
		unsigned int start_rpos = 0, end_rpos = UINT_MAX;

		if (me->exclude)
			continue;

		/* find minimum alignment extent w/o soft clipping (need to consider reverse complemented of B) */
		// If B is RC, then reset read algned to genome B: se->pos to be the one relative to the reversed B
		unsigned int rf_idxb = 0;	/* KSD,BUG,TODO Should be array unsigned int[N_FILES] */
		unsigned int tmp_len = 0;


		/* find minimum extent of alignments on reference 0 */
		for (unsigned int j = 0; j < nalign; ++j) {
			sam_entry *se = &sds[j]->se[me->indices[j][0]];
			unsigned int rf_idx = 0;
            
			/* B reference is reverse complemented: reverse ashes and
			 * starting alignment position is ending alignment position
			 */
			if (j && B_strand) {
                
				if (se->cig->ashes[0].type == CIGAR_SOFT_CLIP)
					tmp_len = se->cig->ashes[0].len;
                /* KSD,TODO Better way is to reverse in place: for (i = 0; i < se2->cig->n_ashes/2; ++i) {unsigned int tmp_type = se2->cig->ashes[se2->cig->n_ashes-i-1].type; se2->cig->ashes[se2->cig->n_ashes-i-1].type = se2->cig->ashes[i].type; se2->cig->ashes[i].type = tmp_type;}
                 * Also, copy and reverse should be done by inline functions, just in case the ashes struct changes definition, although if you do in place, there is less risk since the newly defined parts of ash will not be affected.
                 */
				// Good to know, reverse in place is more efficient
				reverse_in_place(se);
                
				/* starting index is ending index: precompute */
				/*           1-based, exclusive (1-based, inclusive + length = length + 1) */
				rf_idxb = end_pos[j] + 1   - (se->pos + se->cig->length_rf); /* 0-based position from 5' end, inclusive: heavily verified 2/13/22 */

				/* subgenome reverse complemented */
				rf_idx = rf_idxb;	// 0-based 5' start
			} else {
				rf_idx = se->pos - 1 - start_pos[j]; // 0-based 5' start on 0th reference
					/* guaranteed >= 0 */
			}

			if (rf_idx > start_rpos)
				start_rpos = rf_idx;
			rf_idx += se->cig->length_rf;
			if (rf_idx < end_rpos)
				end_rpos = rf_idx;

		}

		debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "Read minimum "
			"extent: %u-%u (0 index on extent: start_pos = %zu)\n",
					 start_rpos, end_rpos, start_pos[0]);

		if (end_rpos <= start_rpos) {
			me->exclude = 1;
			mmessage(INFO_MSG, NO_ERROR, "Read %s excluded because "
				"it does not align to reference [this should "
				"not happen!].\n", sds[0]->se[me->indices[0][0]].name);
			continue;
		}

		/* remake ashes to soft clip to minimum alignment extent */
		for (unsigned int j = 0; j < nalign; ++j) {
			sam_entry *se = &sds[j]->se[me->indices[j][0]];
			unsigned int rf_idx = 0;
			unsigned int first_ash_nlen = 0, n_ashes = 0;
			unsigned int out = 0, diff_extent = 0;
            
			if (j && B_strand)	/* subgenome reverse complemented */
				rf_idx = rf_idxb;	/* Why precomputed? */
			else
				rf_idx = se->pos - 1 - start_pos[j];

			debug_msg(fxn_debug >= DEBUG_I, fxn_debug, 
				"Read %s (%zu), alignment %u: rf_idx=%u"
				" (%u - %u)\n", se->name, n_reads, j,
					rf_idx, start_rpos, end_rpos);

			/* count number of new ashes after matching soft clips */
			for (unsigned int i = 0; i < se->cig->n_ashes; ++i) {	/* recall: cigar reversed if B reverse complement relative to A */
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

				/* Drawings: R indicates 0-based start of next cigar (rf_idx),
				 * [,] indicate inclusive start, end of non-soft-clipped extent of current read pair
				 * - indicates extent of current cigar
				 */
				/* --R--[--... OR --R--... (R,[ coincident): will need to soft-clip */
				if (!n_ashes && rf_idx <= start_rpos) {	/* <= because rf_idx is at 0-based start of next ash, 1-based end of current ash */
					diff_extent = 1;
					++n_ashes;
				/* current ash dropped (included in new soft-clip) */
				} else if (rf_idx <= start_rpos) {
					continue;
				/* [--R--... (R,] coincident) OR [--R--]--... OR --[--R--... (R,] coincident) OR --[--R--]--... */
				} else if (!n_ashes && rf_idx > start_rpos
					&& rf_idx <= end_rpos) {
					/* --[-R--]-- or --[--R-- (R,] coincident): will need to add soft clip */
					if (start_rpos > se->pos - start_pos[j] - 1) {
						diff_extent = 1;
						n_ashes += 2;
					/* --[-R--... OR --[--R--... (R,[ coincident): ash kept */
					} else {
						++n_ashes;
					}
				/* --[--]--R--... OR [--]--R--: will need to add soft-clip, maybe 2 */
				} else if (!n_ashes && rf_idx > end_rpos) {
					diff_extent = 1;
							/* --RR[RRR]RR-- */
					n_ashes += 2 + (start_rpos
						> se->pos - start_pos[j] - 1);
				/* ...--[--R--]--... OR ...--[--R--... (R,] coincident)*/
				} else if (n_ashes && rf_idx <= end_rpos) {
					++n_ashes;
				/* ...--[--]--R--...: need to add 3' soft-clip, possible 5' soft-clip already added */
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
			debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, " [%zu, %zu)\n", se->pos, se->cig->length_rf);

			ash *new_ashes = se->cig->ashes;
			if (diff_extent) {
				debug_msg(fxn_debug >= DEBUG_I, fxn_debug, 
					"Readjusting read %zu from %u ashes to "
					"%u ashes\n", n_reads, se->cig->n_ashes,
									n_ashes);
				new_ashes = malloc(n_ashes * sizeof(*se->cig->ashes));	/* [KSD,BUG] used to be sizeof *se->cig, but would have just allocated too much space */
			}
			
			if (j && B_strand)	/* subgenome reverse complemented */
				rf_idx = rf_idxb;
			else
				rf_idx = se->pos - 1 - start_pos[j];

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
					"first_ash_nlen = %u (extent: %u-%u)\n",
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
						"Adding first ash1: %uS\n", first_ash_nlen);

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
						"Adding first ash2: %u%c\n",
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
							"Adding first ash3: %uS (%u)\n",
								new_ashes[n_ashes].len,
									first_ash_nlen);

						++n_ashes;
					}

					new_ashes[n_ashes].type = se->cig->ashes[i].type;
					new_ashes[n_ashes].len = rf_idx - start_rpos;
					++n_ashes;

					debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
						"Adding first/second ash1: %u%c\n",
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
							"Adding first ash4: %uS (%u)\n", 
							new_ashes[n_ashes].len, first_ash_nlen);

						++n_ashes;
					}

					new_ashes[n_ashes].type = se->cig->ashes[i].type;
					new_ashes[n_ashes].len = end_rpos - start_rpos;

					debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
						"Adding first/second ash2: %u%c\n",
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
					"\tread %s new cigar: ", se->name);	/* TODO,KSD just use name */
				for (unsigned int i = 0; i < n_ashes; ++i)
					debug_msg_cont(fxn_debug >= DEBUG_I,
						fxn_debug, "%u%c",
						new_ashes[i].len,
						cigar_char[new_ashes[i].type]);
				debug_msg_cont(fxn_debug >= DEBUG_I,
					fxn_debug, " (n_ashes = %u)\n", n_ashes);
			}

			// if B is RC, this need to be reverted again for computing the likelihood...
			se->cig->ashes = new_ashes;
			se->cig->n_ashes = n_ashes;
           
			if (B_strand && j) {	/* TODO also shouldn't you do outside the for (j...) loop anyway? */
				// revert the sigar string

				reverse_in_place(se);
				/* recompute se->pos if the first ash is S but the old is not */
				if (se->cig->ashes[0].type == CIGAR_SOFT_CLIP)
					se->pos += (se->cig->ashes[0].len - tmp_len);
//                fprintf(stderr, "%d %d|len %zu: \n", se->cig->ashes[0].len, tmp_len, se->pos);
			}

			/* recompute se->cig->length_rf */
			se->cig->length_rf  = 0;
			
			for (unsigned int l = 0; l < se->cig->n_ashes; ++l) {
//				fprintf(stderr, "len %d: , type %d   ", se->cig->ashes[j].len, se->cig->ashes[j].type);
				if (se->cig->ashes[l].type == CIGAR_DELETION ||
					se->cig->ashes[l].type == CIGAR_MATCH ||
					se->cig->ashes[l].type == CIGAR_MISMATCH ||
					se->cig->ashes[l].type == CIGAR_MMATCH ||
					se->cig->ashes[l].type == CIGAR_SKIP)
					se->cig->length_rf += se->cig->ashes[l].len;
			}
//				fprintf(stderr, "\n");
		}
	}

	return NO_ERROR;
} /* match_extent */

/**
 * Find homoloeologous positions for pair of aligned reference sequences at
 * targeted region, -1, -2 means not mapped; -1 for indel, -2 for end gaps
 *
 * NOTE: the read (NUC) given by mummer is not reversed complemented even when
 * the flag so indicates
 *
 * \A target
 *  \ (pos)       2I           /
 *   \____________  __________/
 *   /------------------ -----\
 *  / (soft clip)       1D     \ (soft clip)
 *                               B target
 * B_to_A (B not reverse complemented):
 * B index =01234567890123456789-0123456 (length 27)
 * A index ===345678901234--5678901234==
 * A_to_B (B not reverse complemented):
 * A index 012345678901234--567890123456 (length 27)
 * B index ===234567890123456789-01234==
 *
 * B_to_A (B reverse complemented):
 * B index =65432109876543210987-6543210
 * A index ===345678901634--5678901234==
 * A_to_B (B reverse complemented):
 * A index 012345678901234--567890123456
 * B index ===432109876543210987-65432==
 * @param rfi	ref_inf object to fill with information
 */
void match_pair(ref_info *rfi)
{
	sam_entry *se = &rfi->ref_sam->se[rfi->rf_idx];
	size_t length = rfi->end_A - rfi->start_A;
	size_t lengthb = rfi->end_B - rfi->start_B;
	size_t rd_idx = 0;	/* 0-based index in reference 1 (B) */

	// in this case, we ignore the insertion of B
	rfi->map_A_to_B = malloc(length * sizeof(*rfi->map_A_to_B));
	rfi->map_B_to_A = malloc(lengthb * sizeof(*rfi->map_B_to_A));
	
	for (size_t j = 0; j < length; ++j)
		rfi->map_A_to_B[j] = ALIGN_SOFT_CLIP;//-2;
	for (size_t j = 0; j < lengthb; ++j)
		rfi->map_B_to_A[j] = ALIGN_SOFT_CLIP;//-2;
	
	size_t rf_idx = se->pos - 1;	/* 0-based position in reference 0 (A) */
	
	for (unsigned int i = 0; i < se->cig->n_ashes; ++i) {
		if (se->cig->ashes[i].type == CIGAR_SOFT_CLIP
		   || se->cig->ashes[i].type == CIGAR_HARD_CLIP) { /* though HC does not consume read, we include it in the targeted genomes alignments since length info in the name is used */
			rd_idx += se->cig->ashes[i].len;
		} else if (se->cig->ashes[i].type == CIGAR_SKIP) {
		} else if (se->cig->ashes[i].type == CIGAR_DELETION) {
			for (unsigned int m = 0; m < se->cig->ashes[i].len; ++m)
				rfi->map_A_to_B[rf_idx + m] = ALIGN_DELETION;//-1;
			rf_idx += se->cig->ashes[i].len;
		} else if (se->cig->ashes[i].type == CIGAR_INSERTION) {
			for (unsigned int m = 0; m < se->cig->ashes[i].len; ++m) {
				if (rfi->strand_B)
					rfi->map_B_to_A[lengthb - (rd_idx + m) - 1] = ALIGN_INSERTION;//-1;
				else
					rfi->map_B_to_A[rd_idx + m] = ALIGN_INSERTION;//-1;
			}
			rd_idx += se->cig->ashes[i].len;
		} else if (se->cig->ashes[i].type == CIGAR_MATCH
			   || se->cig->ashes[i].type == CIGAR_MMATCH
			   || se->cig->ashes[i].type == CIGAR_MISMATCH) {
			for (unsigned int m = 0; m < se->cig->ashes[i].len; ++m) {
				rfi->map_A_to_B[rf_idx + m] = rd_idx + m;
				if (rfi->strand_B)
					rfi->map_B_to_A[lengthb - (rd_idx + m) - 1] = rf_idx + m;
				else
					rfi->map_B_to_A[rd_idx + m] = rf_idx + m;
			}
			rf_idx += se->cig->ashes[i].len;
			rd_idx += se->cig->ashes[i].len;
		}
	}


	if (rfi->strand_B) // if B reverse strand, then reverse mapping
		for (size_t j = 0; j < length; ++j)
			if (rfi->map_A_to_B[j] >= 0)
				rfi->map_A_to_B[j] = lengthb - rfi->map_A_to_B[j] - 1;

/*
	for (size_t j = 0; j < length; ++j)
		fprintf(stderr, "A %zu -> B %d\n", j, rfi->map_A_to_B[j]);
	for (size_t j = 0; j < lengthb; ++j)
		fprintf(stderr, "B %zu -> A %d\n", j, rfi->map_B_to_A[j]);
 */
} /* match_pair() */

//void fprint_usage(FILE *fp, char const * const cmd, void *vopt)
//{
//	options_rf *opt = (options_rf *) vopt;
//	
//	fprintf(fp, "Usage: %s -f <ref_sam_file> <sam_file>\n", cmd);
//	fprintf(fp, "\t-f <ref_sam_file>\t\tWrite fastq output in file '%s'.\n",
//		opt->sam_file);
//	fprintf(fp, "\t-u[nmapped]\t\tRemove unmapped reads.\n");
//	fprintf(fp, "\t-b <sam_file>\t\tRead sam file to parse.\n");
//	fprintf(fp, "\t-d <sam_file>\t\tRead reference fsa file to parse.\n");
//	fprintf(fp, "\t-s <samtools_command>\t\tSamtools location.\n");
//} /* fprint_usage */

//void ll_all_align(sam *sds, const char *ref_file) {
//	mlogit_stuff mls = {NULL, 0};
//	size_t i, j, n_read = 0;
//	fastq_data *fds;
//	fastq_options fop = {.read_encoding = IUPAC_ENCODING, .read_names = 1};
//	int err = NO_ERROR;
//	FILE *fp = NULL;
//	fp = fopen(ref_file, "r");
//	if (!fp)
//		exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR,
//			      ref_file));
//	if ((err = read_fastq(ref_file, &fds, &fop)))
//		exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "Reading '%s' "
//			      "failed with error '%s' (%d).\n",
//			      ref_file, fastq_error_message(err),
//			      err));
//
//	unsigned int fs_index[fds->n_reads];
//	for (unsigned int m = 0; m < fds->n_reads; ++m)
//		fs_index[m] += read_length(fds, m);
//
//	size_t rchar = 0;
//	unsigned char found = 0;
//	char *strand;
//	for (size_t i = 0; i < sds->n_se; ++i) {
//		sam_entry *se = &sds->se[i];
//		se->name_s = NULL;
//
//		/* strand for hashing on strand and name */
//		if ((se->flag & 16) == 0) {
//			strand = "+";
//		} else {
//			strand = "-";
//		}
//
//		size_t length = strlen(se->name) + strlen(strand) + 1;
//		se->name_s = malloc(length);
//		sprintf(se->name_s, "%s%s", se->name, strand);
//
//	}
//
//	for (j = 0; j < sds->n_se; ++j) {
//		sam_entry *se = &sds->se[j];
//		for (i = 0; i < sds->n_ref; ++i) {
//			if (se->ref == i) {
//				se->ll_aln = ll_align(se, n_read,
//						      &fds->reads[fs_index[j]], &mls, 0, se->pos - 1);
//				break;
//			}
//		}
//		fprintf(stderr, "%lf", se->ll_aln);
//		n_read++;
//	}
//	fclose(fp);
//}
