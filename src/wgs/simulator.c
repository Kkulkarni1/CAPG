#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>
#include <limits.h>
#include <unistd.h>
#define MATHLIB_STANDALONE 1
#include <Rmath.h>

#include "simulator.h"
#include "cmdline.h"
#include "io.h"
#include "array.h"

#define PLOIDY 4
int main(int argc, const char *argv[]) {
	int err = NO_ERROR;
	simu_options opt;
	simu_dat *dat = NULL;

	make_simu_opt(&opt);
	if (parse_opt(&opt, argc, argv))
		exit(mmessage(ERROR_MSG, INVALID_CMDLINE, ""));
	
	fastq_data *fds = NULL;
	fastq_options fop = {.read_encoding = XY_ENCODING, .read_names = 1};
	// if use samtools
	//	if (opt.fsa_file) {
	//		char *region = NULL;
	//		unsigned long least, most;
	//		most = rand() % (opt.length + 1); // within 1-opt.length
	//		least = most - opt.len_N;
	//		char temp[strlen(opt.ref_name) + 1];
	//		strcpy(temp, opt.ref_name);
	//		char *ptr = strtok(temp, opt.delim_ref);
	//		size_t length = strlen(ptr) + strlen("-") + strlen(":") + (int)(floor(log10((least))) + 1) + (int)(floor(log10((most))) + 1) + 1;
	//		region = malloc(length);
	//		sprintf(region, "%s%s%zu%s%zu", ptr, ":", least, "-", most);
	//		mmessage(INFO_MSG, NO_ERROR, "sampled region: %s\n", region);
	//		extract_ref(opt.samtools_command, region, opt.fsa_file, opt.extracted_rf);
	//	}
	
	/* read in sampled reference genome */
	if ((err = read_fastq(opt.extracted_rf, &fds, &fop)))
		exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "Reading '%s' "
			      "failed with error '%s' (%d).\n",
			      opt.extracted_rf, fastq_error_message(err), err));
	opt.len_N = fds->n_max_length;

	/* sample another genome with homologous SNPs */
	if (make_data(&dat))
		exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "Make data error \n"));
	
	if (load_data(dat, &opt, fds))
		exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "Load data error \n"));
	
	return(err);
} /* main */

void fprint_fsa(FILE *fp, char_t **data, size_t n, size_t p, char const * const prefix) {
	for (size_t i = 0; i < n; ++i) {
		fprintf(fp, ">%s%lu\n", prefix, i);
		for (size_t j = 0; j < p; ++j)
			fprintf(fp, "%c", xy_to_char[(int)data[i][j]]);
		fprintf(fp, "\n");
	}
}

void write_sam(FILE *fp, char_t *B, size_t p, char *name_A, char *name_B) {
	fprintf(fp, "@SQ	SN:");
	fprintf(fp, "%s	LN:%lu\n", name_A, p);
	fprintf(fp, "%s\t0\t%s\t1\t255\t%luM\t*\t0\t0\t", name_B, name_A, p);
	for (size_t j = 0; j < p; ++j)
		fprintf(fp, "%c", xy_to_char[(int)B[j]]);
	fprintf(fp, " * \n");
}

void fprint_seq(FILE *fp, char_t *data, size_t p, char const * const prefix) {
	fprintf(fp, ">%s\n", prefix);
	for (size_t j = 0; j < p; ++j)
		fprintf(fp, "%c", xy_to_char[(int)data[j]]);
	fprintf(fp, "\n");
}

void make_simu_opt(simu_options *opt)
{
	opt->out_sam = NULL;
	opt->extracted_rf = NULL;
	opt->sub_ref_b = NULL;
	opt->fsa_file = NULL;
	opt->ref_name = NULL;
	opt->len_N = 0;
	opt->out_file = "sample";
	opt->seed = 0;
	opt->alpha = 0.05;
	opt->beta = 0.1;
	opt->heter_rate = 0.001;
	opt->homo_rate = 0.001;
	opt->substitution_rate = 0.3333333;
	opt->mismatch_prob = 0;
	opt->imperfect_ref = 0;
	opt->num_ind = 1;
	
	opt->ART_command = NULL;
	opt->error_file1 = NULL;
	opt->error_file2 = NULL;
	opt->fq_file = "sim";
	opt->coverage = 40;
	opt->length = 150;
	
	opt->bwa_command = NULL;
	opt->samAB = NULL;
} /* make_options */

int make_data(simu_dat **dat) {
	simu_dat *dp;

	*dat = malloc(sizeof(**dat));
	if (*dat == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data object");
	dp = *dat;
	dp->heter_loci = NULL;
	dp->homo_loci = NULL;
	dp->seq_A = NULL;
	dp->seq_B = NULL;
	dp->seq_A2 = NULL;
	dp->seq_B2 = NULL;
	dp->ind = NULL;
	return NO_ERROR;
} /* make_data */

int parse_opt(simu_options *opt, int argc, const char **argv)
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
			case 'f':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				mmessage(INFO_MSG, NO_ERROR, "Fasta file:");
				opt->fsa_file = argv[++i];
				fprintf(stderr, " %s", opt->fsa_file);
				fprintf(stderr, "\n");
				break;
			case 'h':
				fprint_usage(stderr, argv[0], opt);
				exit(EXIT_SUCCESS);
			case 'r':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				opt->out_sam = argv[++i];
				mmessage(INFO_MSG, NO_ERROR, "Reference "
					"alignment: %s\n", opt->out_sam);
				break;
			case 's':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				opt->seed = read_uint(argc, argv, ++i, (void *)opt);
				srand(opt->seed);
				mmessage(INFO_MSG, NO_ERROR, "Seed: %lu\n", opt->seed);
				break;
			case 'o':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				opt->out_file = argv[++i];
				mmessage(INFO_MSG, NO_ERROR, "Out fsa file "
					"base: %s\n", opt->out_file);
				break;
			case 'e':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				opt->extracted_rf = argv[++i];
				mmessage(INFO_MSG, NO_ERROR, "Fasta file:");
				fprintf(stderr, " %s", opt->extracted_rf);
				fprintf(stderr, "\n");
				break;
			case 'p':
				if (i + 1 == argc)
					goto CMDLINE_ERROR;
				opt->mismatch_prob = atof(argv[++i]);
				if (opt->mismatch_prob > 0)
					opt->imperfect_ref = 1;
				else
					opt->mismatch_prob = 0;
				mmessage(INFO_MSG, NO_ERROR, "Probability of "
					"mismatch in subgenome A: %f\n",
							opt->mismatch_prob);
				break;
			case 'j':
				opt->homo_rate = read_cmdline_double(argc, argv, ++i, (void *)opt);
				break;
			case 'g':
				opt->heter_rate = read_cmdline_double(argc, argv, ++i, (void *)opt);
				break;
			case 'a':
				opt->alpha = read_cmdline_double(argc, argv, ++i, (void *)opt);
				break;
			case 'b':
				opt->beta = read_cmdline_double(argc, argv, ++i, (void *)opt);
				break;
			case 'n':
				opt->num_ind = read_uint(argc, argv, ++i, (void *)opt);
				mmessage(INFO_MSG, NO_ERROR, "No. of individuals:");
				fprintf(stderr, " %d", opt->num_ind);
				fprintf(stderr, "\n");
				break;
			case 'l':
				opt->length = read_cmdline_double(argc, argv, ++i, (void *)opt);
				break;
			case 'i':
				opt->coverage = read_cmdline_double(argc, argv, ++i, (void *)opt);
				break;
			case 'c':
				opt->error_file1 = argv[++i];
				break;
			case 'd':
				opt->error_file2 = argv[++i];
				break;
			case 'm':
				opt->fq_file = argv[++i];
				break;
			case 't':
				opt->ART_command = argv[++i];
				mmessage(INFO_MSG, NO_ERROR, "ART command:");
				fprintf(stderr, " %s", opt->ART_command);
				fprintf(stderr, "\n");
				break;
			case 'w':
				opt->bwa_command = argv[++i];
				mmessage(INFO_MSG, NO_ERROR, "bwa command:");
				fprintf(stderr, " %s", opt->bwa_command);
				fprintf(stderr, "\n");
				break;
			case 'k':
				opt->samAB = argv[++i];
				break;
			default:
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
		}
	}
	
	return err;
	
CMDLINE_ERROR:
	if (err == NO_ERROR) {
		err = INVALID_CMD_ARGUMENT;
		i--;
	}
	usage_error(argv, i, (void *)opt);
	return err;
} /* parse_options */

int load_data(simu_dat *dat, simu_options *opt, fastq_data *fds)
{
	unsigned int i, j, m;
	unsigned int count = 0;
	char_t complement[3];
	double r;

	dat->seq_A = malloc(opt->len_N * sizeof(*dat->seq_A));
	dat->seq_B = malloc(opt->len_N * sizeof(*dat->seq_B));
	dat->ref_A = malloc(opt->len_N * sizeof(*dat->ref_A));
	dat->seq_A2 = malloc(opt->len_N * sizeof(*dat->seq_A2));
	dat->seq_B2 = malloc(opt->len_N * sizeof(*dat->seq_B2));
	dat->homo_loci = malloc(opt->len_N * sizeof(*dat->homo_loci));
	dat->heter_loci = malloc(opt->len_N * sizeof(*dat->heter_loci));
	dat->mm_loci = malloc(opt->len_N * sizeof(*dat->homo_loci));
	MAKE_2ARRAY(dat->ind, 4, opt->len_N);
	
	dat->seq_A = fds->reads;

	fprintf(stderr, "start simulation\n");

	/* simulate errors in subgenome A reference */
	for (j = 0; j < opt->len_N; ++j) {

		if (!opt->imperfect_ref) {
			dat->ref_A[j] = dat->seq_A[j];
			dat->mm_loci[j] = 0;
			continue;
		}

		/* simulate errors in subgenome reference A */
		count = 0;
		r = rand() / (RAND_MAX + 1.);
		if (r <= opt->mismatch_prob) {
			dat->mm_loci[j] = 1;
			for(i = 0; i < 4; ++i)
				if (i != dat->seq_A[j])
					complement[count++] = i;
			double rn = rand() / (RAND_MAX + 1.);
			if (rn <= opt->substitution_rate)
				dat->ref_A[j] = complement[0];
			else if (rn <= (opt->substitution_rate + opt->substitution_rate))
				dat->ref_A[j] = complement[1];
			else
				dat->ref_A[j] = complement[2];
		} else {
			dat->mm_loci[j] = 0;
			dat->ref_A[j] = dat->seq_A[j];
		}
	}
	fprintf(stderr, "reference A simulated, mismatch rate %lf, mutated:\n", opt->mismatch_prob);

	/* read subgenome B reference */
	if (fds->n_reads == 2) {
		unsigned int cnt_h = 0;

		dat->seq_B = fds->reads + opt->len_N;

		for (j = 0; j < opt->len_N; ++j) {
			if (dat->seq_A[j] != dat->seq_B[j]) {
				dat->homo_loci[j] = 1;
				++cnt_h;
			} else {
				dat->homo_loci[j] = 0;
			}
		}

		fprintf(stderr, "genome B read from '%s', observed homoeologous rate %lf, mutated:\n", opt->extracted_rf, (double) cnt_h / opt->len_N);

	/* simulate reference B */
	} else {

		for (j = 0; j < opt->len_N; ++j) {

			r = rand() / (RAND_MAX + 1.);
			count = 0;

			if (r <= opt->homo_rate) {
				dat->homo_loci[j] = 1;
				for(i = 0; i < 4; ++i)
					if (i != dat->seq_A[j])
						complement[count++] = i;
				// sample substitution nuc
				double rn = rand() / (RAND_MAX + 1.);
				//			fprintf(stderr, "%f ",rn);
				if (rn <= opt->substitution_rate)
					dat->seq_B[j] = complement[0];
				else if (rn <= (opt->substitution_rate + opt->substitution_rate))
					dat->seq_B[j] = complement[1];
				else
					dat->seq_B[j] = complement[2];
			} else {
				dat->homo_loci[j] = 0;
				dat->seq_B[j] = dat->seq_A[j];
			}
	
		}
		fprintf(stderr, "genome B simulated, homoeologous rate %lf, mutated:\n", opt->homo_rate);
	}

	/* truth to stderr */
	for (j = 0; j < opt->len_N; ++j) {
		if (dat->homo_loci[j] == 1)
			fprintf(stderr, "%d: %d|%d\n", j, dat->seq_A[j], dat->seq_B[j]);
		if (dat->mm_loci[j])
			fprintf(stderr, "mm: %d: %d|%d\n", j, dat->seq_A[j], dat->ref_A[j]);
	}

	// ref A nad B in the same file for HMM method to use
	FILE *fp = NULL;
	size_t len = strlen(opt->fsa_file) + 4 + 1;
	char *fsa_file = malloc(len);

	sprintf(fsa_file, "%s.fsa", opt->fsa_file);
	if (fsa_file) {
		fp = fopen(fsa_file, "w");
		if (!fp)
			exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR, fsa_file));
	}
	size_t name_len = strlen("Genome_A:0-") + (int)log10(opt->len_N) + 1 + 1;
	char *name_A = malloc(name_len);
	char *name_B = malloc(name_len);

	sprintf(name_A, "Genome_A:0-%u", opt->len_N);
	sprintf(name_B, "Genome_B:0-%u", opt->len_N);
	if (opt->imperfect_ref)
		fprint_seq(fp, dat->ref_A, opt->len_N, name_A);
	else
		fprint_seq(fp, dat->seq_A, opt->len_N, name_A);
	fprint_seq(fp, dat->seq_B, opt->len_N, name_B);
	fclose(fp);
	fp = NULL;
	if (opt->out_sam) {
		fp = fopen(opt->out_sam, "w");
		if (!fp)
			exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt->out_sam));
	}
	// sam file: alignment of A and B
	write_sam(fp, dat->seq_B, opt->len_N, name_A, name_B);
	fclose(fp);
	fp = NULL;
	len = len + 1;
	char *fsaA_file = malloc(len);
	char *fsaB_file = malloc(len);
	// ref A nad B in the seperate files for alignment of reads
	sprintf(fsaA_file, "%sA.fsa", opt->fsa_file);
	fp = fopen(fsaA_file, "w");
	if (opt->imperfect_ref)
		fprint_seq(fp, dat->ref_A, opt->len_N, "Genome_A");
	else
		fprint_seq(fp, dat->seq_A, opt->len_N, "Genome_A");
	fclose(fp);
	fp = NULL;
	sprintf(fsaB_file, "%sB.fsa", opt->fsa_file);
	fp = fopen(fsaB_file, "w");
	fprint_seq(fp, dat->seq_B, opt->len_N, "Genome_B");
	fclose(fp);
	fp = NULL;
	if (opt->bwa_command) {
		call_bwa(opt, fsaA_file, NULL, NULL, NULL);
		call_bwa(opt, fsaB_file, NULL, NULL, NULL);
	}

	if (!opt->num_ind)
		return NO_ERROR;

	// make heter loci
	for (j = 0; j < opt->len_N; ++j) {
		if (dat->homo_loci[j]) {
			dat->heter_loci[j] = 0;
			continue;
		}
		double r = rand() / (RAND_MAX + 1.);
		if (r <= (opt->heter_rate + opt->heter_rate)) {
			double rn = rand() / (RAND_MAX + 1.);
			//			fprintf(stderr, "%f ",rn);
			if (rn <= 0.5)
				dat->heter_loci[j] = 1;
			else
				dat->heter_loci[j] = 2;
		} else {
			dat->heter_loci[j] = 0;
		}
	}
	// make the rest genome A, B
	for (j = 0; j < opt->len_N; ++j) {
		char_t complement[3];
		unsigned int count = 0;
		if(dat->heter_loci[j] == 0) {
			dat->seq_A2[j] = dat->seq_A[j];
			dat->seq_B2[j] = dat->seq_B[j];
		} else if (dat->heter_loci[j] == 1) {
			dat->seq_B2[j] = dat->seq_B[j];
			double rn = rand() / (RAND_MAX + 1.);
			for(i = 0; i < PLOIDY; ++i)
				if (i != dat->seq_A[j])
					complement[count++] = i;
			if (rn <= opt->substitution_rate)
				dat->seq_A2[j] = complement[0];
			else if (rn <= (opt->substitution_rate + opt->substitution_rate))
				dat->seq_A2[j] = complement[1];
			else
				dat->seq_A2[j] = complement[2];
			// randomly select A1 or A2 to mutate
			double which_seq = rand() / (RAND_MAX + 1.);
			char_t tmp = dat->seq_A2[j];
			if (which_seq > 0.5) {
				dat->seq_A2[j] = dat->seq_A[j];
				dat->seq_A[j] = tmp;
			}
		} else {
			dat->seq_A2[j] = dat->seq_A[j];
			double rn = rand() / (RAND_MAX + 1.);
			for(i = 0; i < PLOIDY; ++i)
				if (i != dat->seq_B[j])
					complement[count++] = i;
			if (rn <= opt->substitution_rate)
				dat->seq_B2[j] = complement[0];
			else if (rn <= (opt->substitution_rate + opt->substitution_rate))
				dat->seq_B2[j] = complement[1];
			else
				dat->seq_B2[j] = complement[2];
			// randomly select A1 or A2 to mutate
			double which_seq = rand() / (RAND_MAX + 1.);
			char_t tmp = dat->seq_B2[j];
			if (which_seq <= 0.5) {
				dat->seq_B2[j] = dat->seq_B[j];
				dat->seq_B[j] = tmp;
			}
		}
	}
	fprintf(stderr, "\nallelic SNPs simulated, with rate %lf, mutated:\n", opt->heter_rate);
	for (j = 0; j < opt->len_N; ++j) {
		if(dat->heter_loci[j] == 1)
			fprintf(stderr, "A: %d: %d|%d \n", j, dat->seq_A[j], dat->seq_A2[j]);
		else if(dat->heter_loci[j] == 2)
			fprintf(stderr, "B: %d: %d|%d \n", j, dat->seq_B[j], dat->seq_B2[j]);
	}
	set_seed(opt->seed, opt->seed);
	opt->prop_allele = rbeta(opt->alpha, opt->beta);
	double p = opt->prop_allele * opt->prop_allele;
	double pq = 2 * opt->prop_allele * (1 - opt->prop_allele);
	fprintf(stderr, "\nhwe prob:%lf, p^2: %lf 2pq: %lf\n", opt->prop_allele, p, pq);
	
	for (m = 0; m < opt->num_ind; ++m) {
		fprintf(stderr, "Individual %d\n", m);
		for (j = 0; j < opt->len_N; ++j) {
			int hwe = 0; // indicate which genotype
			if (dat->homo_loci[j] == 0 && dat->heter_loci[j] == 0)
				for (i = 0; i < PLOIDY; ++i)
					dat->ind[i][j] = dat->seq_A[j];
			else if (dat->homo_loci[j] == 1 && dat->heter_loci[j] == 0) {
				for (i = 0; i < 2; ++i)
					dat->ind[i][j] = dat->seq_A[j];
				for (i = 2; i < 4; ++i)
					dat->ind[i][j] = dat->seq_B[j];
				fprintf(stderr, "++ %d: %c%c|%c%c\n", j,xy_to_char[dat->ind[0][j]],xy_to_char[dat->ind[1][j]], xy_to_char[dat->ind[2][j]], xy_to_char[dat->ind[3][j]]);
			}
			else if (dat->homo_loci[j] == 0 && dat->heter_loci[j] != 0){
				double rn = rand() / (RAND_MAX + 1.);
				if (rn <= p)
					hwe = 0;
				else if (rn <= p + pq)
					hwe = 1;
				else
					hwe = 2;
				if (dat->heter_loci[j] == 1) {
					for (i = 2; i < 4; ++i)
						dat->ind[i][j] = dat->seq_B[j];
					if (hwe == 0) {
						for (i = 0; i < 2; ++i)
							dat->ind[i][j] = dat->seq_A[j];
					} else if (hwe == 1) {
						fprintf(stderr, "** ");
						dat->ind[0][j] = dat->seq_A[j];
						dat->ind[1][j] = dat->seq_A2[j];
						fprintf(stderr, "%d: %c%c|%c%c\n", j, xy_to_char[dat->ind[0][j]],xy_to_char[dat->ind[1][j]], xy_to_char[dat->ind[2][j]], xy_to_char[dat->ind[3][j]]);
					} else {
						for (i = 0; i < 2; ++i)
							dat->ind[i][j] = dat->seq_A2[j];
					}
					if (dat->seq_B[j] != dat->ind[0][j] && hwe != 1) {
						fprintf(stderr, "++ %d: %c%c|%c%c\n", j, xy_to_char[dat->ind[0][j]],xy_to_char[dat->ind[1][j]], xy_to_char[dat->ind[2][j]], xy_to_char[dat->ind[3][j]]);
					}
				} else {
					for (i = 0; i < 2; ++i)
						dat->ind[i][j] = dat->seq_A[j];
					if (hwe == 0) {
						for (i = 2; i < 4; ++i)
							dat->ind[i][j] = dat->seq_B[j];
					} else if (hwe == 1) {
						fprintf(stderr, "** ");
						dat->ind[2][j] = dat->seq_B[j];
						dat->ind[3][j] = dat->seq_B2[j];
						fprintf(stderr, "%d: %c%c|%c%c\n", j,xy_to_char[dat->ind[0][j]],xy_to_char[dat->ind[1][j]], xy_to_char[dat->ind[2][j]], xy_to_char[dat->ind[3][j]]);
					} else {
						for (i = 2; i < 4; ++i)
							dat->ind[i][j] = dat->seq_B2[j];
					}
					if (dat->seq_A[j] != dat->ind[2][j] && hwe != 1) {
						fprintf(stderr, "++ %d: %c%c|%c%c\n", j, xy_to_char[dat->ind[0][j]],xy_to_char[dat->ind[1][j]], xy_to_char[dat->ind[2][j]], xy_to_char[dat->ind[3][j]]);
					}
				}
			}
		}
		int le = 1;
		if (m != 0)
			le = m;
		
		size_t length = strlen(opt->out_file) + (int)log10(le) + 5 + 1;
		char *file = malloc(length);
		sprintf(file, "%s%u.fsa", opt->out_file, m);
		if (file) {
			fp = fopen(file, "w");
			if (!fp)
				exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR, file));
		}
		fprint_fsa(fp, dat->ind, PLOIDY, opt->len_N, "H");
		fclose(fp);
		fp = NULL;
		if (opt->ART_command) {
			length = strlen(opt->fq_file) + (int)log10(le) + 1 + 1;
			char *fq_out = malloc(length);
			sprintf(fq_out, "%s%u", opt->fq_file, m);
			call_art(opt, fq_out, file);
			if (opt->bwa_command) {
				size_t l = strlen(opt->samAB) + (int)log10(le) + 6 + 1;
				char *sam_out = malloc(l);
				char *fq1 = malloc(length + 3 + 1);
				char *fq2 = malloc(length + 3 + 1);
				if (opt->error_file2) {
					sprintf(fq1, "%s1.fq", fq_out);
					sprintf(fq2, "%s2.fq", fq_out);
					sprintf(sam_out, "%s%uA.sam", opt->samAB, m);
					call_bwa(opt, fsaA_file, fq1, fq2, sam_out);
					sprintf(sam_out, "%s%uB.sam", opt->samAB, m);
					call_bwa(opt, fsaB_file, fq1, fq2, sam_out);
				} else {
					sprintf(fq1, "%s.fq", fq_out);
					sprintf(sam_out, "%s%uA.sam", opt->samAB, m);
					call_bwa(opt, fsaA_file, fq1, NULL, sam_out);
					sprintf(sam_out, "%s%uB.sam", opt->samAB, m);
					call_bwa(opt, fsaB_file, fq1, NULL, sam_out);
				}
			}
		}
	}
	return NO_ERROR;
} /* load_data */

void call_art(simu_options *opt, char *fq_out, char *fq_in) {
	if (opt->error_file2) {
		unsigned int cmd_len = strlen(opt->ART_command) + strlen(opt->error_file1) + strlen(opt->error_file2) + strlen(" -l -1 -2 -f -o -i -p -m 300 -s 10") + strlen(fq_out) + strlen(fq_in) + (int)(log10(opt->coverage) + 1) + (int)(log10(opt->length) + 1) + 8;
		mmessage(INFO_MSG, NO_ERROR, "Length of command: %u\n", cmd_len);
		char *command = malloc(cmd_len * sizeof *command);
		sprintf(command, "%s -1 %s -2 %s -l %d -i %s -f %d -o %s -p -m 300 -s 10",
			opt->ART_command, opt->error_file1, opt->error_file2, opt->length, fq_in, opt->coverage, fq_out);
		mmessage(INFO_MSG, NO_ERROR, "Running ART: '%s'\n",
			 command);
		system(command);
		free(command);
	} else {
		unsigned int cmd_len = strlen(opt->ART_command) + strlen(opt->error_file1) + strlen(" -l -1 -f -o -i ") + strlen(fq_out) + strlen(fq_in) + (int)(log10(opt->coverage) + 1) + (int)(log10(opt->length) + 1) + 8;
		mmessage(INFO_MSG, NO_ERROR, "Length of command: %u\n", cmd_len);
		char *command = malloc(cmd_len * sizeof *command);
		sprintf(command, "%s -1 %s -l %d -i %s -f %d -o %s",
			opt->ART_command, opt->error_file1, opt->length, fq_in, opt->coverage, fq_out);
		mmessage(INFO_MSG, NO_ERROR, "Running ART: '%s'\n",
			 command);
		system(command);
		free(command);
	}
	
}

void call_bwa(simu_options *opt, char *ref_in, char *reads1, char *reads2, char *sam) {
	if(!sam) {
		unsigned int cmd_len = strlen(opt->bwa_command) + strlen(" index ") + strlen(ref_in) + 8;
		mmessage(INFO_MSG, NO_ERROR, "Length of command: %u\n", cmd_len);
		char *command = malloc(cmd_len * sizeof *command);
		sprintf(command, "%s index %s",
			opt->bwa_command, ref_in);
		system(command);
		free(command);
	} else {
		if (opt->error_file2) {
			unsigned int cmd_len = strlen(opt->bwa_command) + strlen(" mem ") + strlen(ref_in) + strlen(reads1)
			+ strlen(reads2) + strlen(" > ") + strlen(" ") + strlen(sam) + 8;
			char *command = malloc(cmd_len * sizeof *command);
			sprintf(command, "%s mem %s %s %s > %s",
				opt->bwa_command, ref_in, reads1, reads2, sam);
			mmessage(INFO_MSG, NO_ERROR, "Running bwa: '%s'\n",
				 command);
			system(command);
			free(command);
		} else {
			unsigned int cmd_len = strlen(opt->bwa_command) + strlen(" mem ") + strlen(ref_in) + strlen(reads1)
			+ strlen(" > ") + strlen(" ") + strlen(sam) + 8;
			char *command = malloc(cmd_len * sizeof *command);
			sprintf(command, "%s mem %s %s > %s",
				opt->bwa_command, ref_in, reads1, sam);
			mmessage(INFO_MSG, NO_ERROR, "Running bwa: '%s'\n",
				 command);
			system(command);
			free(command);
		}
		
	}
	
}


void fprint_usage(FILE *fp, const char *cmdname, void *obj) {
	simu_options *opt = (simu_options *) obj;
	size_t start = strlen(cmdname) - 1;
	
	while (cmdname[start] != '/' && start) start--;
	if (cmdname[start] == '/') start++;
	
	for (size_t i = start; i < strlen(cmdname); ++i)
		fputc(toupper(cmdname[i]), fp);
	fprintf(fp, "(%d)\n", 1);
	fprintf(fp, "\nNAME\n\t%s - Generate simulated individuals\n", &cmdname[start]);
	fprintf(fp, "\nSYNOPSIS\n\t%s -e <fsa> -f <rf_fsa> -j <homo_rate> -g <heter_rate> -a <alpha> \n -b <beta> -s <seed> -o <output> -n <sample> -r <sam>\t\n", &cmdname[start]);
	fprintf(fp, "\nOPTIONS\n");
	fprintf(fp, "\t-e <fsa> \n\t\tSubgenome A (and optionally B) references in FASTA format.\n");
	fprintf(fp, "\t-r <sam> \n\t\tOutput name of SAM file with alignment of A and B reference.\n");
	fprintf(fp, "\t-f <rf_fsa> \n\t\tOutput basename of FASTA files for simulated reference genomes\n");
	fprintf(fp, "\t-o <outfile> \n\t\tOutput filenames of individual genomes (Default: %s)\n", opt->out_file);
	fprintf(fp, "\t-j <homo_rate>\n\t\tSpecify homoeologous rate\n");
	fprintf(fp, "\t-g <heter_rate>\n\t\tSpecify heterozygous rate\n");
	fprintf(fp, "\t-p <mm_rate>\n\t\tSpecify mismatch rate in reference A relative to true subgenome A\n");
	fprintf(fp, "\t-a <alpha>\n\t\tSpecify the shape parameter in beta distribution\n");
	fprintf(fp, "\t-b <beta>\n\t\tSpecify the rate parameter in beta distribution\n");
	fprintf(fp, "\t-s <seed>\n\t\tRandom number generator seed\n");
	fprintf(fp, "\t-n <sample>\n\t\tNo. of individuals simulated\n");
	fprintf(fp, "\t-l <length>\n\t\tSpecify length of simulated reads\n");
	fprintf(fp, "\t-i <coverage>\n\t\tSpecify coverage of the individual genome\n");
	fprintf(fp, "\t-c <error_file1>\n\t\tSpecify error file 1\n");
	fprintf(fp, "\t-d <error_file2>\n\t\tSpecify error file 2\n");
	fprintf(fp, "\t-m <fq_file>\n\t\tFastq file name (Default: %s)\n", opt->fq_file);
	fprintf(fp, "\t-t <ART>\n\t\tART command(only use it and bwa with -n = 1 and one coverage)\n");
	fprintf(fp, "\t-w <bwa>\n\t\tBWA command\n");
	fprintf(fp, "\t-k <rsam>\n\t\tRead aligned to reference name\n");
	fprintf(fp, "\t-h \n\t\tThis help\n");
} /* fprint_usage */

