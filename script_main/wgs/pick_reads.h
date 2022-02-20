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
//typedef struct _ref_entry ref_entry;
typedef struct _ref_options options_rf;

struct _ref_options {
	char const *sam_file;			// reference sam file
	char const *rsam_files[N_FILES]; 	// read sam files
	char const *samtools_command;
	char *extracted_rf[N_FILES]; 		// targeted reference fsa files
	char const *fsa_files[N_FILES];	// reference fsa files
	char filter_unmapped;
	char delim_ref;
	char delim_len;
	char const *fastq_file;
};

/**
 * Store alignment of subgenomic references in homoeologous region.
 */
//struct _ref_entry {
struct _ref_info {
	sam *ref_sam;
	char *name_A;		/*<! name/csome of subgenome A */
	char *name_B;		/*<! name/csome of subgenome B */
	int *map_A_to_B; 	/*<! which base in B is aligned to each base in A within target region */

	size_t rf_idx;		/*<! index of sam entry for selected target */
	size_t start_A; 	/*<! selected homoeologous region in subgenome A, 1-based, inclusive */
	size_t start_B;		/*<! selected homoeologous region in subgenome B, 1-based, inclusive */
	size_t end_A;		/*<! selected homoeologous region in subgenome A, 1-based, exclusive */
	size_t end_B;		/*<! selected homoeologous region in subgenome B, 1-based, exclusive */
	unsigned int strand_B;

	//ref_entry *info;
};

void ll_all_align(sam *sds, char const *ref_file);
void default_options_rf(options_rf *opt);
int make_targets_info(options_rf *opt, ref_info **ref_info, char const *ref_names[]);
int pickreads(ref_info *ref_info, sam **sds, char const **csome_names);
int parse_rf_options(options_rf *opt, int argc, char *argv[]);
int extract_ref(char const *samtools_command, char const *ref_name, size_t ref_start, size_t ref_end, char const *ref_file, char const *ext_rf, options_rf *opt);
void output_selected_reads(char const *f, sam **sds, merge_hash *mh);
void match_pair(ref_info *rf_info);
int match_soft_clipping(merge_hash *mh, unsigned int nalign, sam **sds, unsigned int b_rc);
int match_extent(merge_hash *mh, unsigned int nalign, sam **sds, size_t *start_pos, size_t *end_pos, unsigned int B_strand);

#endif /* pick_reads_h */
