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

int soft_clip_alignment(sam_entry *se, unsigned int fiveprime_sc, unsigned int threeprime_sc, int rc);


void default_options_rf(options_rf *opt)
{
	opt->sam_file = NULL;
	opt->filter_unmapped = 1;
	opt->delim_ref = ":";
	opt->delim_len = "-";
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
 * @return		error status
 */
int make_targets_info(options_rf opt, ref_info **ref_in)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	sam *sd_ref  = NULL;
	FILE *fp = fopen(opt.sam_file, "r");

	if (!fp)
		exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt.sam_file));

	read_sam(fp, &sd_ref);	/* assumes XY_ENCODING */
	fclose(fp);

	ref_info *rf_info;
	
	*ref_in = malloc(sizeof **ref_in);
	if (!*ref_in)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "ref_in");
	rf_info = *ref_in;

	/* get the length of each reference name */
	size_t *rchar = NULL;
	rchar = malloc(sd_ref->n_ref * sizeof(*rchar));
	size_t rchar_in = 0;
	for (unsigned int i = 0; i < sd_ref->n_ref; ++i) {
		rchar[i] = rchar_in;
		rchar_in += strlen(&sd_ref->ref_names[rchar_in]) + 1;
	}
	//
	//	for (size_t i = 0; i < sd_ref->n_se; ++i) {
	//		sam_entry *se = &sd_ref->se[i];
	//		printf("%d\t", se->ref);
	//		printf("%s %s\n", se->name, &sd_ref->ref_names[rchar[se->ref]]);
	//	}
	
	//store the name[include chrosome] and starting and ending positions in the targeted A, B genome
	// i.e. split aradu.V14167.gnm2.chr02:696301-697301
	ref_entry *re = malloc(sd_ref->n_se * sizeof *re);
	
	for (size_t i = 0; i < sd_ref->n_se; ++i) {
		sam_entry *se = &sd_ref->se[i];

		re[i].name_A = NULL;
		re[i].name_B = NULL;

		/* filter unmapped */
		if (opt.filter_unmapped && se->flag >> 2 & 1)
			continue;
		
		/* strand of A, B, in order to find the reads */
//		re[i].strand_A = 0;
		if ((se->flag & 16) == 0)
			re[i].strand_B = 0;
		else
			re[i].strand_B = 1;

		/* parse subgenome B csome name, start and end */
		char temp_B[strlen(se->name) + 1];
		strcpy(temp_B, se->name);
		char *ptr_B = strtok(temp_B, opt.delim_ref);
		re[i].name_B = malloc(strlen(ptr_B) + 1);
		strcpy(re[i].name_B, ptr_B);
		ptr_B = strtok(NULL, opt.delim_ref);
		char *ptr_B_pos = strtok(ptr_B, opt.delim_len);
		re[i].start_B = atoi(ptr_B_pos);
		ptr_B_pos = strtok(NULL, opt.delim_len);
		re[i].end_B = atoi(ptr_B_pos);

		debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
			  "name: %s start %zu end %zu\n", re[i].name_B,
						re[i].start_B, re[i].end_B);

		/* parse subgenome A csome name, start and end */
		char temp_A[strlen(&sd_ref->ref_names[rchar[se->ref]]) + 1];
		strcpy(temp_A, &sd_ref->ref_names[rchar[se->ref]]);
		char *ptr_A = strtok(temp_A, opt.delim_ref);
		re[i].name_A = malloc(strlen(ptr_A) + 1);
		strcpy(re[i].name_A, ptr_A);
		ptr_A = strtok(NULL, opt.delim_ref);
		char *ptr_A_pos = strtok(ptr_A, opt.delim_len);
		re[i].start_A = atoi(ptr_A_pos);
		ptr_A_pos = strtok(NULL, opt.delim_len);
		re[i].end_A = atoi(ptr_A_pos);

		debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
			  "name: %s start %zu end %zu\n", re[i].name_A,
						re[i].start_A, re[i].end_A);

		/* extract alignment start and end positions in subgenome A */
//		re[i].real_sA = re[i].start_A + se->pos; // 1 based start_A
//		re[i].real_eA = re[i].real_sA + se->cig->length_rf - 1;

		/* length of subgenome B consumed in alignment */
		size_t length = 0;
		for (unsigned int m = 0; m < se->cig->n_ashes; ++m)
			if (se->cig->ashes[m].type == CIGAR_INSERTION
				|| se->cig->ashes[m].type == CIGAR_MATCH
				|| se->cig->ashes[m].type == CIGAR_MMATCH
				|| se->cig->ashes[m].type == CIGAR_MISMATCH)
				length += se->cig->ashes[m].len;

		/*   1    100               200  230
		 * A -----|------^---vv-----|-----
		 * B -----|    ----^-v------|-----------
		 *   450  500  510          597        1087
		 *
		 * start_A = 100, end_A = 200
		 * start_B = 500, end_B = 597
		 * real_sA = 100, real_eA = 200
		 * real_sB = 510, real_eB = 597
		 *
		 *        100               200
		 * A -----|------^---vv-----|-----
		 * B -----|    ----^-v------|----------- (reverse complemented)
		 *   647  597  587          500        10
		 *
		 * start_A = 100, end_A = 200
		 * start_B = 500, end_B = 597
		 * real_sA = 100, real_eA = 200
		 * real_sB = 500, real_eB = 587
		 */

//		if (re[i].strand_B) { // if B reversed
//			fprintf(stderr, "Genome B is reverse complemented\n");
//			re[i].real_eB = re[i].end_B;
//			if (se->cig->ashes[se->cig->n_ashes - 1].type
//							== CIGAR_SOFT_CLIP
//				|| se->cig->ashes[se->cig->n_ashes - 1].type
//							 == CIGAR_HARD_CLIP)
//				re[i].real_eB -= se->cig->ashes[se->cig->n_ashes - 1].len;
//			re[i].real_sB = re[i].real_eB - length + 1;
//		} else {
//			re[i].real_sB = re[i].start_B + 1; //1 based
//			if (se->cig->ashes[0].type == CIGAR_SOFT_CLIP
//				|| se->cig->ashes[0].type == CIGAR_HARD_CLIP)
//				re[i].real_sB += se->cig->ashes[0].len;
//			re[i].real_eB = re[i].real_sB + length - 1; // 1-based
//		}
	}
	rf_info->info = re;
	rf_info->ref_sam = sd_ref;
	free(rchar);
	return NO_ERROR;
}/* make_targets_info */

/**
 * Find reads aligned to homeologous reference regions. The reads are aligned to
 * whole subgenome A and separately to whole subgenome B. We are focused on 
 * subsets of homeologous target regions. Here, we seek those reads that align
 * to any of these homoeologous regions.
 *
 * @param ref_info	information about homoeologous aligned reference regions
 * @param opt		options about homeologous reference regions
 * @param sds		sam file objects of reads aligned to subgenomes
 * @return		error status
 */
int pickreads(ref_info *ref_info, options_rf *opt, sam **sds)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	unsigned int i, j, m;
	size_t *rchar[N_FILES];

	/* get start index of each reference name */
	for (j = 0; j < N_FILES; ++j) {
		size_t rchar_in = 0;

		rchar[j] = malloc(sds[j]->n_ref * sizeof **rchar);
		
		for (i = 0; i < sds[j]->n_ref; ++i) {
			rchar[j][i] = rchar_in;

			debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "%zu %s\n",
				rchar[j][i], &sds[j]->ref_names[rchar[j][i]]);
			
			rchar_in += strlen(&sds[j]->ref_names[rchar_in]) + 1;
		}
//		for (m = 0; m < sds[j]->n_se; ++m) {
//			sam_entry *se = &sds[j]->se[m];
//			if ((se->flag & (1 << 2)))
//				continue;
//			cigar *cig = se->cig;
//			printf("%d: %zu ", m, cig->length_rf);
//			printf("%s\t",  &sds[j]->ref_names[rchar[j][se->ref]]);
//		}
//		printf("\n");
	}

	/* find the reference region the read aligns to, if any */
	for (j = 0; j < N_FILES; ++j) {
		for (m = 0; m < sds[j]->n_se; ++m) {
			sam_entry *se = &sds[j]->se[m];
			cigar *cig = se->cig;

			/* [KSD, TODO] Restore which_ref to unsigned int.  Set se->exclude = 1 by default and unset here if you find a reference. Or if we may have already set exclude,
			 * use a local found_ref = 0, and if you find no matching ref here, set se->exclude = 1.
			 */
			unsigned int found_ref = 0;
//			se->which_ref = -1;
//			se->ref_name = NULL;

			/* skip unmapped reads */
			if (se->flag & (1 << 2))
				continue;

			/* [TODO,KSD] slow; write a direct map or use htslib
			 * to extract reads for each reference alignment; leave as is for now
			 */
			for (i = 0; i < ref_info->ref_sam->n_se; ++i) {
				sam_entry *rse = &ref_info->ref_sam->se[i];
				ref_entry *re = &ref_info->info[i];

				/* search only primary reference alignments */
				if (rse->flag >> 11 & 1)
					continue;

				char *ref_names[N_FILES] = {re->name_A, re->name_B};
//				printf("%s %s\n ", ref_names[j], &sds[j]->ref_names[rchar[j][se->ref]]);

				/* this read maps to this reference csome */
				if (!strcmp(&sds[j]->ref_names[rchar[j][se->ref]], ref_names[j])) {
					size_t start_pos[N_FILES] = {re->start_A, re->start_B};	// 0 based, inclusive (from chrom_name:start-end)
					size_t end_pos[N_FILES] = {re->end_A, re->end_B};	// 1-based, inclusive (0-based, exclusive)
					size_t rf_index_s = se->pos - 1;
					size_t rf_index_e = rf_index_s + cig->length_rf;

//					printf("%zu %zu || %zu %zu\n", rf_index_s, rf_index_e, start_pos[j], end_pos[j]);

					/* and it maps within the target region */
					if ((rf_index_s >= start_pos[j] && rf_index_e <= end_pos[j]) ||		/* read contained within target */
					    (rf_index_s <= start_pos[j] && rf_index_e > start_pos[j]) ||	/* read crosses 5' end of target */
					    (rf_index_s < end_pos[j] && rf_index_e >= end_pos[j])) {		/* read crosses 3' end of target */
						size_t length = strlen(ref_names[j]) + strlen(opt->delim_len) + strlen(opt->delim_ref) + (int)(log10(end_pos[j]) + 1) + 1;

						if (start_pos[j] != 0)
							length += (int)(log10(start_pos[j]) + 1);
						else
							length += 1;
						se->ref_name = malloc(length);
						sprintf(se->ref_name, "%s%s%zu%s%zu", ref_names[j], opt->delim_ref, start_pos[j], opt->delim_len, end_pos[j]);
						se->which_ref = i;
						found_ref = 1;
//						debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "REF_ID: %d REF_NAME: %s \n", se->which_ref, se->ref_name);
						break;
					}
				}
			}
			if (!found_ref) {
				se->exclude = 1;
				/*mmessage(INFO_MSG, NO_ERROR, "Read %s (%u) "
					"excluded because it does not align to "
					"one of references.\n", se->name_s, m);*/
			}
		}
	}

	for (j = 0; j < N_FILES; ++j)
		free(rchar[j]);

	return NO_ERROR;
}/* pickreads */

/**
 * Extract selected regions from fasta file. [Requires that someone ran samtools
 * faidx on the reference file previously.] Actually does not require this.
 *
 * TODO,KSD Check for prior run of faidx; provide helpful error message if not.
 *  I think samtools does not require index the reference files before extract the regions, if there is no index files, they will create one
 *
 * @param samtools_command	samtools executable
 * @param region		samtools region specification: chr:from-to
 * @param ref_file		fasta file to extract regions from
 * @param ext_rf		fasta file to write with chosen region
 */
void extract_ref(const char *samtools_command, char *region,
				const char *ref_file, const char *ext_rf)
{
	// index the whole reference genome file
	if (!ref_file)
		mmessage(ERROR_MSG, FILE_NOT_FOUND, "Reference file '%s' not found\n", ref_file);
    
	unsigned int cmd_len = strlen(samtools_command) + strlen(ref_file)
		+ strlen(" faidx   -o ") + strlen(region) + strlen(ext_rf) + 1;
	char *command = NULL;

	mmessage(INFO_MSG, NO_ERROR, "Length of command: %u\n", cmd_len);

	command = malloc(cmd_len * sizeof *command);	/* KSD,TODO,BUG Check for malloc() failure; this is a repeated TODO in many places. */
	if (!command)
		mmessage(ERROR_MSG, MEMORY_ALLOCATION, "samtools command");
    
	sprintf(command, "%s faidx %s %s -o %s",
		samtools_command, ref_file, region, ext_rf);
	
	mmessage(INFO_MSG, NO_ERROR, "Running samtools: '%s'\n", command);

	system(command);

	free(command);
	
}/* extract_ref */

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
				opt->delim_len = argv[++argv_idx];
				break;
			case 'r':
				opt->delim_ref = argv[++argv_idx];
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

// output the reads aligned to the target including the one not aligned to both
void output_selected_reads(const char *f, sam **sds, merge_hash *mh) {
	FILE *fpp = NULL;
	fpp = fopen(f, "w");
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

		//debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "Read %s 5' soft clip %u; 3' soft clip %u\n", sds[0]->se[me->indices[0][0]].name_s, fiveprime_sc, threeprime_sc);

		for (unsigned int j = 0; j < nalign; ++j) {
			sam_entry *se = &sds[j]->se[me->indices[j][0]];

			 if (soft_clip_alignment(se, fiveprime_sc,
			 			threeprime_sc, j && b_rc)) {
				++n_reads_excluded;
				mmessage(INFO_MSG, NO_ERROR, "Read %s excluded "
					"by soft-clipping entire alignment in "
					"subgenome %u.\n", se->name_s, j);
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

debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "Read %s with cigar ", se->name_s);
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

debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "Read %s with cigar ", se->name_s);
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
				"not happen!].\n", sds[0]->se[me->indices[0][0]].name_s);
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
			debug_msg_cont(fxn_debug >= DEBUG_I, fxn_debug, " [%zu, %zu)\n", se->pos, se->length_rf);

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
					"\tread %s new cigar: ", se->name_s);	/* TODO,KSD just use name */
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

void delta_position(ref_info *rf_info, size_t ref_id)
{
	ref_entry *re = &rf_info->info[ref_id];
	sam_entry *se = &rf_info->ref_sam->se[ref_id];
	unsigned int rd_idx = 0;

	size_t length = re->end_A - re->start_A;

	if (!length)
		return;

	re->delta_len = calloc(length, sizeof(*re->delta_len));

	/* first ash */
	if (se->cig->ashes[0].type == CIGAR_SOFT_CLIP ||
		se->cig->ashes[0].type == CIGAR_HARD_CLIP) {
		rd_idx += se->cig->ashes[0].len;
		re->delta_len[rd_idx] = se->pos - se->cig->ashes[0].len;
	} else if (se->cig->ashes[0].type == CIGAR_MATCH
		   || se->cig->ashes[0].type == CIGAR_MMATCH
		   || se->cig->ashes[0].type == CIGAR_MISMATCH) {
		rd_idx += se->cig->ashes[0].len;
	} else if (se->cig->ashes[0].type == CIGAR_INSERTION) {
		for (unsigned int j = 0; j <= se->cig->ashes[0].len; ++j)	/* include first nucleotide in next ash */
			re->delta_len[j] -= j + 1;
		rd_idx += se->cig->ashes[0].len;
	} else if (se->cig->ashes[0].type == CIGAR_DELETION) {
		re->delta_len[0] = se->cig->ashes[0].len;
	}

	for (unsigned int i = 1; i < se->cig->n_ashes; ++i) {
		if (se->cig->ashes[i].type == CIGAR_SOFT_CLIP ||	/* 3' clip */
			se->cig->ashes[i].type == CIGAR_HARD_CLIP) {
			for (unsigned int j = 1; j < se->cig->ashes[i].len; ++j)
				re->delta_len[rd_idx + j] = re->delta_len[rd_idx];
			rd_idx += se->cig->ashes[i].len;
		} else if (se->cig->ashes[i].type == CIGAR_MATCH
			   || se->cig->ashes[i].type == CIGAR_MMATCH
			   || se->cig->ashes[i].type == CIGAR_MISMATCH) {
			for (unsigned int j = 1; j < se->cig->ashes[i].len; ++j)
				re->delta_len[rd_idx + j] = re->delta_len[rd_idx];
			rd_idx += se->cig->ashes[i].len;
		} else if (se->cig->ashes[i].type == CIGAR_INSERTION) {
			for (unsigned int j = 0; j < se->cig->ashes[i].len; ++j)
				re->delta_len[rd_idx + j] = re->delta_len[rd_idx] - j - 1;
			rd_idx += se->cig->ashes[i].len;
		}
		if (i + 1 < se->cig->n_ashes) {	/* pass to first nucleotide of next ash */
			if (se->cig->ashes[i].type == CIGAR_DELETION)
				re->delta_len[rd_idx] += se->cig->ashes[i].len;
			else
				re->delta_len[rd_idx] = re->delta_len[rd_idx - 1];
		}
	}
} /* delta_position */

// find homology position pair for the selected reference, -1 means not mapped (or deletion in B)
// notice the read (NUC) given by mummer is not reversed complemented even if the flag shows it is
void match_pair(ref_info *rf_info, size_t ref_id)
{
	sam *sd_ref = rf_info->ref_sam;
	sam_entry *se = &sd_ref->se[ref_id];
	ref_entry *re = &rf_info->info[ref_id];
	size_t length = 0;
	unsigned int rd_idx = 0;	/* index in reference 1 */

	// get the length of paired aligned reference  (length of reference)
//	for (unsigned int i = 0; i < se->cig->n_ashes; ++i) {
//		if (se->cig->ashes[i].type == CIGAR_SOFT_CLIP
//		    || se->cig->ashes[i].type == CIGAR_HARD_CLIP
//		    || se->cig->ashes[i].type == CIGAR_INSERTION
//		    || se->cig->ashes[i].type == CIGAR_MATCH
//		    || se->cig->ashes[i].type == CIGAR_MMATCH
//		    || se->cig->ashes[i].type == CIGAR_MISMATCH)
//			length += se->cig->ashes[i].len;
//	}
	// in this case, we ignore the insertion of B
	length = re->end_A - re->start_A; // [start_A, end_A)
	re->idx_map = malloc(length * sizeof(*re->idx_map));
	
	for (size_t j = 0; j < length; ++j)
		re->idx_map[j] = -1;
	
	int rf_idx = se->pos - 1; // 0-based
	
	for (unsigned int i = 0; i < se->cig->n_ashes; ++i) {
		if (se->cig->ashes[i].type == CIGAR_SOFT_CLIP
		   || se->cig->ashes[i].type == CIGAR_HARD_CLIP) { /* though HC does not consume read, we include it in the targeted genomes alignments since length info in the name is used */
			rd_idx += se->cig->ashes[i].len;
			continue;
		} else if (se->cig->ashes[i].type == CIGAR_DELETION
			    || se->cig->ashes[i].type == CIGAR_SKIP) {
			rf_idx += se->cig->ashes[i].len;
		} else if (se->cig->ashes[i].type == CIGAR_INSERTION) {
			rd_idx += se->cig->ashes[i].len;
		} else if (se->cig->ashes[i].type == CIGAR_MATCH
			   || se->cig->ashes[i].type == CIGAR_MMATCH
			   || se->cig->ashes[i].type == CIGAR_MISMATCH) {
			for (unsigned int m = rf_idx; m < rf_idx + se->cig->ashes[i].len; ++m)
				re->idx_map[m] = rd_idx + m - rf_idx;
			rf_idx += se->cig->ashes[i].len;
			rd_idx += se->cig->ashes[i].len;
		}
	}
	size_t lengthb = re->end_B - re->start_B;
//	if (se->cig->ashes[se->cig->n_ashes - 1].type == CIGAR_SOFT_CLIP
//	    || se->cig->ashes[se->cig->n_ashes - 1].type == CIGAR_HARD_CLIP)
//		lengthb -= se->cig->ashes[se->cig->n_ashes - 1].len;
	if(re->strand_B) // if B reverse strand, then reverse mapping
		for (size_t j = 0; j < length; ++j)
			if(re->idx_map[j] != -1)
				re->idx_map[j] = lengthb - re->idx_map[j] - 1;
//	fprintf(stderr, "2 references alignment map\n");
//	for (size_t j = 0; j < length; ++j)
//		fprintf(stderr, "%d ", re->idx_map[j]);
//	fprintf(stderr, "\n");
}

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
