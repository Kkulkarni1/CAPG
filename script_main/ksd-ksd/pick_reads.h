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
#include "roshan.h"
#include "array.h"

typedef struct _ref_info ref_info;
typedef struct _ref_entry ref_entry;
typedef struct _ref_options options_rf;

struct _ref_options {
	char const *sam_file;			/*<! sam file of ref alignment */
	const char *rsam_files[N_SUBGENOMES];	/*<! sam files of reads */
	const char *samtools_command;		/*<! samtools executable */
	char *extracted_rf[N_SUBGENOMES]; 	/*<! fsa files of targeted ref */
	const char *fsa_files[N_SUBGENOMES];	/*<! reference fsa files */
	char filter_unmapped;			/*<! filter unmapped reads? */
	char *delim_ref;			/*<! char to tokenize ref names */
	char *delim_len;			/*<! char to tokenize ref posns */
	const char *fastq_file;			/*<! ? */
};

struct _ref_entry {
	char *name_A;
	char *name_B;
	size_t start_A; 		/*<! start position on subgenome A */
	size_t start_B;			/*<! start position on subgenome B */
	size_t end_A;			/*<! end position on subgenome A */
	size_t end_B;			/*<! end position on subgenome B */
	size_t real_sA; 		/*<! start/end after alignment */
	size_t real_eA;
	size_t real_sB;
	size_t real_eB;
	unsigned int strand_A;  	/* 0:forward, 1:reverse */
	unsigned int strand_B;
	int *idx_map;		 	/* which base in B is aligned to which in A  */
};


struct _ref_info {
	ref_entry *info;
	sam *ref_sam;
};

void ll_all_align(sam *sds, const char *ref_file);
void default_rf_options(options_rf *opt);
int make_targets_info(options_rf opt, ref_info **ref_info);
int pickreads(ref_info *ref_info, options_rf *opt, sam **sds);
int parse_rf_options(options_rf *opt, int argc, char *argv[]);
void extract_ref(const char *samtools_command, char *region, const char *ref_file, const char *ext_rf);
void output_selected_reads(const char *f, sam **sds, merge_hash *mh);
void match_pair(ref_info *rf_info, size_t my_refs);

#endif /* pick_reads_h */
