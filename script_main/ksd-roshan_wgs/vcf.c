#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "vcf.h"
#include "math.h"
#include "error.h"

void print_vcf_header(FILE *fp, vcf_options *vopt, const char *ref_file, const char *sample_id)
{
    char timestamp[9];
    time_t now = time(NULL);

    fprintf(fp, "##fileformat=VCFv%.1f\n", VCF_VERSION);
    strftime(timestamp, 9, "%Y%m%d", localtime(&now));
    fprintf(fp, "##fileDate=%s\n", timestamp);
    fprintf(fp, "##source=roshanV%.1f\n", ROSHAN_VERSION);
    fprintf(fp, "##reference=%s\n", ref_file);
    fprintf(fp, "##phasing=no\n");
//    fprintf(fp, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n");
    if (!vopt->genotype_by_clustering)
        fprintf(fp, "##FILTER=<ID=al2,Description=\"Found evidence of more than 2 alleles\">\n");
    if (!vopt->genotype_by_clustering && vopt->coverage_screen > 0)
        fprintf(fp, "##FILTER=<ID=c%2.0f,Description=\"Subgenomic "
            "coverage at a site less than %2.0f%% of total "
            "subgenomic coverage\">\n", vopt->coverage_screen * 100,
                        vopt->coverage_screen * 100);
    if (!vopt->genotype_by_clustering && vopt->posthoc_coverage_test) {
        int mgq = vopt->phc_min_genotype_pp < 1
            ? MIN(99, (int) (-10 * log10(1 - vopt->phc_min_genotype_pp)))
            : 99;
        fprintf(fp, "##FILTER=<ID=sc5,Description=\"Equal "
            "genomic coverage test failed at level 0.05\">\n");
        fprintf(fp, "##FILTER=<ID=gq%d,Description=\"Genotype "
            "quality low, posterior probability below %.2f\">\n",
            mgq, vopt->phc_min_genotype_pp);
    }
    if (!vopt->genotype_by_clustering && vopt->equal_homolog_coverage_test)
        fprintf(fp, "##FORMAT=<ID=ET,Number=1,Type=Float,Description="
                "\"-log10(pvalue) of equal coverage test\">\n");
    if (vopt->genotype_by_clustering) {
        fprintf(fp, "##FILTER=<ID=sgi,Description=\"Subgenomic "
            "indels\">\n");
        fprintf(fp, "##FILTER=<ID=hmi,Description=\"Homoeologous "
            "indel at this site\">\n");
    }
    fprintf(fp, "##FORMAT=<ID=GT,Number=1,Type=String,Description="
                            "\"Genotype\">\n");
    fprintf(fp, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description="
                            "\"Read Depth\">\n");
    if (vopt->output_gl)
        fprintf(fp, "##FORMAT=<ID=GL,Number=G,Type=Float,Description="
                        "\"Genotype Likelihoods\">\n");
/*
    fprintf(fp, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description="
                "\"Phred-scaled Genotype Likelihoods\">\n");
*/
    fprintf(fp, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="
                        "\"Genotype Quality\">\n");
    fprintf(fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                    "%s\n", sample_id ? sample_id : "sample");
} /* print_vcf_header */

int make_default_vcf_options(vcf_options **vopt)
{
	*vopt = malloc(sizeof **vopt);
	
	if (!*vopt)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "vcf_options");
	
	(*vopt)->output_gl = 0;
	(*vopt)->genotype_by_clustering = 0;
	(*vopt)->posthoc_coverage_test = 0;
	(*vopt)->coverage_screen = 0;
	
	return NO_ERROR;
} /* make_default_vcf_options */
