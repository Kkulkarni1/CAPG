#ifndef __VCF_H__
#define __VCF_H__

#define VCF_VERSION 4.2
#define ROSHAN_VERSION 1.0

typedef struct _vcf_options vcf_options;

struct _vcf_options {
	int output_gl;			/*<! output genotype likelihoods */
	
	int genotype_by_clustering;
	double coverage_screen;
	int posthoc_coverage_test;
    int equal_homolog_coverage_test;
	double phc_min_genotype_pp;
};


void print_vcf_header(FILE *fp, vcf_options *vopt, const char *ref_file, const char *sample_id);
int make_default_vcf_options(vcf_options **vopt);
#endif
