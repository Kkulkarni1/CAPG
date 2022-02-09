//
//  pick_reads.h
//  project
//
//  Created by Yudi Zhang on 3/24/20.
//  Copyright © 2020 Yudi Zhang. All rights reserved.
//

#ifndef pick_reads_h
#define pick_reads_h

#include <stdint.h>
#include "sam.h"
#include "capg.h"
#include "array.h"

typedef struct _ref_info ref_info;
typedef struct _ref_entry ref_entry;
typedef struct _ref_options options_rf;

struct _ref_options {
	char const *sam_file;			// reference sam file
	const char *rsam_files[N_FILES]; 	// read sam files
	const char *samtools_command;
	char *extracted_rf[N_FILES]; 		// targeted reference fsa files
	const char *fsa_files[N_FILES];	// reference fsa files
	char filter_unmapped;
	char *delim_ref;
	char *delim_len;
	const char *fastq_file;
};

/**
 * Store alignment of subgenomic references in homoeologous region.
 */
struct _ref_entry {
	char *name_A;		/*<! name/csome of subgenome A */
	char *name_B;		/*<! name/csome of subgenome B */
	size_t start_A; 	/*<! homoeologous region start in subgenome A */
	size_t start_B;		/*<! homoeologous region start in subgenome B */
	size_t end_A;		/*<! homoeologous region end in subgenome A */
	size_t end_B;		/*<! homoeologous region end in subgenome B */
	/* KSD,QUESTION,TODO The following appear to be set but not used. Delete? */
    // yeah, delete them
//	size_t real_sA; 	/*<! homoeologous region alignment start in subgenome A */
//	size_t real_eA;		/*<! homoeologous region alignment end in subgenome A */
//	size_t real_sB;		/*<! homoeologous region alignment start in subgenome B */
//	size_t real_eB;		/*<! homoeologous region alignment end in subgenome B */
	/* KSD,QUESTION,TODO strand_A is not used, right? strand_B should be [N_FILES] array */
    //[N_FILES] array for different ployploid?
//	unsigned int strand_A;  	/* 0:forward, 1:reverse */
	unsigned int strand_B;
	int *idx_map;	 	/*<! which base in B is aligned to each base in A */
};


struct _ref_info {
	ref_entry *info;
	sam *ref_sam;
};

void ll_all_align(sam *sds, const char *ref_file);
void default_options_rf(options_rf *opt);
int make_targets_info(options_rf opt, ref_info **ref_info);
int pickreads(ref_info *ref_info, options_rf *opt, sam **sds);
int parse_rf_options(options_rf *opt, int argc, char *argv[]);
void extract_ref(const char *samtools_command, char *region, const char *ref_file, const char *ext_rf);
void output_selected_reads(const char *f, sam **sds, merge_hash *mh);
void match_pair(ref_info *rf_info, size_t my_refs);

#endif /* pick_reads_h */
