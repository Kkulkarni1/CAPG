//
//  pick_reads.c
//  pick reads aligned to each targeted region
//
//  Created by Yudi Zhang on 3/24/20.
//  Copyright © 2020 Yudi Zhang. All rights reserved.
//
// MERGING TMP: temporarily commented until I can figure this out
// MERGING TMP ADD: temporarily added until I figure this out

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


void default_rf_options(options_rf *opt)
{
	opt->sam_file = NULL;
	opt->filter_unmapped = 1;
	opt->delim_ref = ":";
	opt->delim_len = "-";
	opt->samtools_command = "samtools";
	for (int i = 0; i < N_SUBGENOMES; ++i) {
		opt->rsam_files[i] = NULL;
		opt->fsa_files[i] = NULL;
		opt->extracted_rf[i] = NULL;
	}
	opt->fastq_file = NULL;
	
} /* default_rf_options */

int make_targets_info(options_rf opt, ref_info **ref_in)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	sam *sd_ref  = NULL;
	ref_info *rf_info;

	/* reference alignment in sam format */
	FILE *fp = fopen(opt.sam_file, "r");
	
	if (!fp)
		exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt.sam_file));
	
	read_sam(fp, &sd_ref);	/* assumes XY_ENCODING [TODO] Is this a good idea? */
	fclose(fp);

	/* allocate reference information object */
	*ref_in = malloc(sizeof **ref_in);
	if (!*ref_in)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "ref_in");
	rf_info = *ref_in;

	/* get the length of each reference name */
	size_t *rchar = NULL;
	rchar = malloc(sd_ref->n_ref * sizeof *rchar);
	size_t rchar_in = 0;

	for (unsigned int i = 0; i < sd_ref->n_ref; ++i) {
		rchar[i] = rchar_in;
		rchar_in += strlen(&sd_ref->ref_names[rchar_in]) + 1;
	}

	/* store the name, including chrosome, start and end positions in genome
	 * i.e., split aradu.V14167.gnm2.chr02:696301-697301
	 */
	ref_entry *re = malloc(sd_ref->n_se * sizeof *re);
	
	for (size_t i = 0; i < sd_ref->n_se; ++i) {
		sam_entry *se = &sd_ref->se[i];

		re[i].name_A = NULL;
		re[i].name_B = NULL;

		/* filter unmapped */
		if (opt.filter_unmapped && se->flag >> 2 & 1)
			continue;
		
		/* strand of A, B, in order to find the reads */
		re[i].strand_A = 0;
		if ((se->flag & 16) == 0)	/* RC(B) aligned to A */
			re[i].strand_B = 0;
		else
			re[i].strand_B = 1;

		/* find name, start and end in subgenome B */
		char temp_B[strlen(se->name) + 1];
		strcpy(temp_B, se->name);
		char *ptr_B = strtok(temp_B, opt.delim_ref);
		re[i].name_B = malloc(strlen(ptr_B) + 1);
		strcpy(re[i].name_B, ptr_B);
		ptr_B = strtok(NULL, opt.delim_ref);	/* TODO: not thread safe */
		char *ptr_B_pos = strtok(ptr_B, opt.delim_len);
		re[i].start_B = atoi(ptr_B_pos);	/* TODO: left size_t, right int */
		ptr_B_pos = strtok(NULL, opt.delim_len);
		re[i].end_B = atoi(ptr_B_pos);
		debug_msg(fxn_debug >= DEBUG_I, fxn_debug,
			  "name: %s start %zu end %zu\n", re[i].name_B, re[i].start_B, re[i].end_B);

		/* find name, start and end in subgenome A */
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
			  "name: %s start %zu end %zu\n", re[i].name_A, re[i].start_A, re[i].end_A);
		// extract the real aligned regions (if it do is 0-based!)
		re[i].real_sA = re[i].start_A + se->pos; // 1 based start_A
// MERGING TMP		re[i].real_eA = re[i].real_sA + se->cig->length_rf - 1; //
		
		size_t length = 0;
		for (unsigned int m = 0; m < se->cig->n_ashes; ++m)
			if (se->cig->ashes[m].type == CIGAR_INSERTION
			    || se->cig->ashes[m].type == CIGAR_MATCH
			    || se->cig->ashes[m].type == CIGAR_MMATCH
			    || se->cig->ashes[m].type == CIGAR_MISMATCH)
				length += se->cig->ashes[m].len;
		if (re[i].strand_B) { // if B reversed
//			fprintf(stderr, "Genome B is reverse complemented\n");
			re[i].real_eB = re[i].end_B;
			if (se->cig->ashes[se->cig->n_ashes - 1].type == CIGAR_SOFT_CLIP || se->cig->ashes[se->cig->n_ashes - 1].type == CIGAR_HARD_CLIP)
				re[i].real_eB -= se->cig->ashes[se->cig->n_ashes - 1].len;
			re[i].real_sB = re[i].real_eB - length + 1;
		} else {
			re[i].real_sB = re[i].start_B + 1; //1 based
			if (se->cig->ashes[0].type == CIGAR_SOFT_CLIP || se->cig->ashes[0].type == CIGAR_HARD_CLIP)
				re[i].real_sB += se->cig->ashes[0].len;
			re[i].real_eB = re[i].real_sB + length - 1; // 1-based
		}
	}
	rf_info->info = re;
	rf_info->ref_sam = sd_ref;
	free(rchar);
	return NO_ERROR;
}

int pickreads(ref_info *ref_info, options_rf *opt, sam **sds) {
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//
	
	size_t **rchar = NULL;
	rchar = malloc(N_SUBGENOMES * sizeof(*rchar));
	unsigned int i, j, m;
	
	for (j = 0; j < N_SUBGENOMES; ++j){
		// get the reference names index
		rchar[j] = malloc((sds[j]->n_ref) * sizeof(**rchar));
		size_t rchar_in = 0;
		
		for (i = 0; i < sds[j]->n_ref; ++i) {
			rchar[j][i] = rchar_in;
			debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "%zu %s\n ", rchar[j][i], &sds[j]->ref_names[rchar[j][i]]);
			
// MERGING TMP			rchar_in += strlen(&sds[j]->ref_names[rchar_in]) + 1;
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
	
	for (j = 0; j < N_SUBGENOMES; ++j) {
		for (m = 0; m < sds[j]->n_se; ++m) {
			sam_entry *se = &sds[j]->se[m];
			se->which_ref = -1;
// MERGING TMP			se->ref_name = NULL;
			if ((se->flag & (1 << 2)))
				continue;
			cigar *cig = se->cig;
			for (i = 0; i < ref_info->ref_sam->n_se; ++i) {
				sam_entry *rse = &ref_info->ref_sam->se[i];
				if (rse->flag >> 11 & 1) //skip secondary
					continue;
				ref_entry *re = &ref_info->info[i];
				char *ref_names[N_SUBGENOMES] = {re->name_A, re->name_B};
				size_t start_pos[N_SUBGENOMES] = {re->start_A, re->start_B}; // 0 based
				size_t end_pos[N_SUBGENOMES] = {re->end_A, re->end_B}; // 1-based
//				printf("%s %s\n ", ref_names[j], &sds[j]->ref_names[rchar[j][se->ref]]);
				if (!strcmp(&sds[j]->ref_names[rchar[j][se->ref]], ref_names[j])) {
					size_t len_ref = end_pos[j] - start_pos[j];
					size_t rf_index_s = se->pos - 1;
// MERGING TMP					size_t rf_index_e = rf_index_s + cig->length_rf;
					size_t rf_index_e = 0;// MERGING TMP ADD
//					printf("%zu %zu || %zu %zu\n", rf_index_s, rf_index_e, start_pos[j], end_pos[j]);
// MERGING TMP					if(cig->length_rf < len_ref) {
					if (1) {// MERGING TMP ADD
						if ((rf_index_s >= start_pos[j] && rf_index_e <= end_pos[j]) ||
						    (rf_index_s < start_pos[j] && rf_index_e >= start_pos[j]) ||
						    (rf_index_s < end_pos[j] && rf_index_e >= end_pos[j])) {
							size_t length = strlen(ref_names[j]) + strlen(opt->delim_len) + strlen(opt->delim_ref) + (int)(log10(end_pos[j]) + 1) + 1;
							if (start_pos[j] != 0) {
								length += (int)(log10(start_pos[j]) + 1);
							} else
								length += 1;
// MERGING TMP							se->ref_name = malloc(length);
// MERGING TMP							sprintf(se->ref_name, "%s%s%zu%s%zu", ref_names[j], opt->delim_ref, start_pos[j], opt->delim_len, end_pos[j]);
							se->which_ref = i;
// MERGING TMP							debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "REF_ID: %d REF_NAME: %s \n", se->which_ref, se->ref_name);
							break;
						}
					} else {
						if(start_pos[j] >= rf_index_s && end_pos[j] <= rf_index_e) {
							size_t length = strlen(ref_names[j]) + strlen(opt->delim_len) + strlen(opt->delim_ref) + (int)(log10(start_pos[j]) + 1) + (int)(log10(end_pos[j]) + 1) + 1;
// MERGING TMP							se->ref_name = malloc(length);
// MERGING TMP							sprintf(se->ref_name, "%s%s%zu%s%zu", ref_names[j], opt->delim_ref, start_pos[j], opt->delim_len, end_pos[j]);
							se->which_ref = i;
// MERGING TMP							debug_msg(fxn_debug >= DEBUG_I, fxn_debug, "REF_ID: %d REF_NAME: %s \n", se->which_ref, se->ref_name);
							break;
						}
					}
				}
			}
		}
	}
	for (j = 0; j < N_SUBGENOMES; ++j)
		free(rchar[j]);
	free(rchar);
	return NO_ERROR;
}

// extract the region after using samtools, do samtools faidx reference.fasta first, then samtools faidx ref.fasta name:1-100
void extract_ref(const char *samtools_command, char *region, const char *ref_file, const char *ext_rf) {
	
	unsigned int cmd_len = strlen(samtools_command) + strlen(ref_file) + strlen(" faidx ") + strlen(region) + strlen(" -o ") + strlen(ext_rf) + 8;
	mmessage(INFO_MSG, NO_ERROR, "Length of command: %u\n", cmd_len);
	char *command = malloc(cmd_len * sizeof *command);
	sprintf(command, "%s faidx %s %s -o %s",
		samtools_command, ref_file, region, ext_rf);
	
	mmessage(INFO_MSG, NO_ERROR, "Running samtools: '%s'\n",
		 command);
	system(command);
	free(command);
	
}

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
				for (j = 0; j < N_SUBGENOMES; ++j) {
					opt->rsam_files[j] = argv[++argv_idx];
					fprintf(stderr, " %s", opt->rsam_files[j]);
				}
				fprintf(stderr, "\n");
				break;
			case 'd':
				for (j = 0; j < N_SUBGENOMES; ++j) {
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
		if (me->nfiles != N_SUBGENOMES) {
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
}

// find homology position pair for the selected reference, -1 means not mapped (or deletion in B)
// notice the read (NUC) given by mummer is not reversed complemented even if the flag shows it is
void match_pair(ref_info *rf_info, size_t my_refs) {
	sam *sd_ref = rf_info->ref_sam;
	sam_entry *se = &sd_ref->se[my_refs];
	ref_entry *re = &rf_info->info[my_refs];
	size_t length = 0;
	int rd_idx = 0;
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
	length = re->end_A - re->start_A; // ROSHAN'S index end is not actually end position, it + 1
	re->idx_map = NULL;
	re->idx_map = malloc(length * sizeof(*re->idx_map));
	
	for (size_t j = 0; j < length; ++j)
		re->idx_map[j] = -1;
	
	int rf_idx = se->pos - 1; // 0-based
	
	for (unsigned int i = 0; i < se->cig->n_ashes; ++i) {
		if (se->cig->ashes[i].type == CIGAR_SOFT_CLIP
		   || se->cig->ashes[i].type == CIGAR_HARD_CLIP) { // although HC does not consume read, in the targeted genomes alignments, we include it since the length info in the name is used
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
	fprintf(stderr, "2 references alignment map\n");
	for (size_t j = 0; j < length; ++j)
		fprintf(stderr, "%d ", re->idx_map[j]);
	fprintf(stderr, "\n");
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
