#ifndef __SAM_H__
#define __SAM_H__

#include <stdint.h>
#include <zlib.h>

#include "sequence.h"
#include "uthash.h"

/**
 * Censor quality scores at maximum Illumina quality score.
 * Comment to stop censoring.
 */
//#define CENSOR_QUALITY_SCORES

/*
@SQ	SN:BP48	LN:2001
@PG	ID:bwa	PN:bwa	VN:0.7.15-r1140	CL:bwa mem -t 8 -a targetsB.fsa Olin.fastq
M00259:19:000000000-AUAFV:1:1101:11435:2052	16	BP18	915	60	48M1I133M	*	0	0	TCATGGCCGTCTTCACTTTCGAGGATGAAATCACCTCCCCCGTGCCCCACTGCCAAACGTTACAATGCTATGAAGGATGCGGATTCTATCACCCCTAAGATTATTGATGACATCAAGAGTGTTGAAATCGTTGGGGGAAACGGTGGTCCTGGAACCATCAAGAAACTCACCATTGTCGAGGG	HHHHHGHHHHHGGHGGHHFHHHHHHGHHHHGFFFFFCGHECFFHHFGFGHHHHHHHHGHHFGGGGDFEDHGHHHHFFEHFHHGGGGGGGHHGHGHHHHFGFGEFHHHHHHFHFGHHHHFFFFEEFEFF?EE?GGGFFHGGFEF9AA0;1BFFGFB3G3BHGFGGF1;1B19D11>1BFA>>1	NM:i:12	MD:Z:41C7C5G1T9C11C2C2C0C23G21A48	AS:i:119	XS:i:0
*/

typedef struct _sam sam;
typedef struct _sam_entry sam_entry;
typedef struct _cigar cigar;
typedef struct _ash ash;

enum {
	CIGAR_MMATCH,	/* (mis)match: ref, read consumed */
	CIGAR_INSERTION,/* insertion: only read consumed */
	CIGAR_DELETION,	/* deletion: only ref consumed */
	CIGAR_SKIP,	/* skip: only ref consumed */
	CIGAR_SOFT_CLIP,/* soft clip: only read consumed */
	CIGAR_HARD_CLIP,/* hard clip: neither consumed */
	CIGAR_PAD,	/* pad: neither consumed */
	CIGAR_MATCH,	/* match: read, ref consumed */
	CIGAR_MISMATCH,	/* mismatch: read, ref consumed */
	CIGAR_NCHAR
};

extern char cigar_char[CIGAR_NCHAR];

struct _ash {
	unsigned int type : 4;
	unsigned int len;
};

struct _cigar {
	size_t length_rf; 		/* total consumed reference length */
	unsigned int n_ashes;
	ash *ashes;
};

struct _sam_entry {			/* see VCF specification */
	char *name;			/* NAME */
	sequence *read;			/* sequence of read */
	sequence *qual;			/* qualities of read */
	cigar *cig;			/* CIGAR */
	size_t pos;			/* POS */
	double ll_aln;			/* alignment likelihood */
	unsigned int which_ref; 	/* which external reference */
	unsigned int index; 		/* index used to output data */
	uint32_t ref;			/* index of reference: see sam */
	uint16_t flag;			/* FLAG */
	unsigned char exclude;		/* exclude this SAM entry */
//	char *ref_name;			/* ref_name(w/ location): A:1-20 */
//	char *name_s; 			/* read name w/ +:forward, -:reverse */
};

struct _sam {
	sam_entry *se;		/* sam entries */
	char *rname_ptr;	/* contiguous memory for reference names */
	char **ref_names;	/* pointers to reference names */

	/* for indexing to reference, either internal or external */
	size_t *n_per_ref;	/* no. entries per [external] ref */
	size_t **ref_list;	/* entries per [external] ref */

	size_t n_se;		/* no. sam entries: alignments or reads */
	uint32_t n_ref;		/* no. of references */
	size_t hash_length;	/* length of hash */
	size_t n_mapping;	/* no. of mappings */
};

#define HASH_NAME	1
#define HASH_READ	2
#define HASH_REFERENCE	4

typedef struct _sam_hash sam_hash;
typedef struct _merge_hash merge_hash;

struct _sam_hash {
	size_t idx;
	unsigned int count;
	size_t *indices;
	int type;
	UT_hash_handle hh;
};

struct _merge_hash {
	UT_hash_handle hh;
	size_t *count;
	size_t **indices;
	unsigned int nfiles;
	unsigned char exclude;	/* for internal use in sam.* */
};

/* input */
int read_sam(FILE *fp, sam **s_in, unsigned char, unsigned char);
#ifdef USE_GZIP
int read_bam(gzFile gzfp, sam **s_in);
#endif

/* processing */
int hash_sam(sam *s, sam_hash **sh_in, int hash_on, size_t n_ref, unsigned char drop_unmapped, unsigned char drop_second, unsigned int drop_soft_clipped, unsigned int drop_indel, unsigned int min_length, unsigned int max_length, double max_exp_err);
size_t hash_merge(merge_hash **mh, unsigned int nfiles, sam **sds, size_t *rindex);
int fill_hash(sam *s, sam_hash *sh, unsigned int n_ref);
/* output */
int output_error_data(FILE *fp, sam_entry *se, unsigned char *ref, double lprob);
void reverse_in_place(sam_entry *se);
void print_cigar(FILE *file, sam_entry *se);


#endif
