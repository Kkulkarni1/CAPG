/**
 * @file run_sam.c
 * @author K. S. Dorman
 *
 * This file will be the interface for processing sam files.  Currently all
 * it can do is prepare BCB 568 homework.  Also see roshan.c for another
 * use of sam.c.
 *
 * To compile:
make sam
 */

#include <stdlib.h>
#include <string.h>

#include "sam.h"
#include "fastq.h"
#include "nuc.h"
#include "qual.h"
#include "uthash.h"
#include "io.h"
#include "cmdline.h"
#include "error.h"

void usage(char const * const msg);

typedef struct _options options;

struct _options {
	unsigned int min_length;
	unsigned int max_length;
	unsigned int start_posn;
	unsigned int end_posn;
	char const * ref_name;
	char const * fasta_file;
	char const * fastq_file;
	char const * sam_file;
	int print_lengths;
	int sample_id;
	unsigned int errors_only;
	unsigned int forward;
	unsigned int reverse;
	char filter_unmapped;
};

void default_options(options *opt)
{
	opt->min_length = 0;
	opt->max_length = 0;
	opt->start_posn = 0;
	opt->end_posn = 0;
	opt->ref_name = NULL;
	opt->fasta_file = NULL;
	opt->fastq_file = NULL;
	opt->sam_file = NULL;
	opt->print_lengths = 0;
	opt->sample_id = -1;
	opt->filter_unmapped = 1;
	opt->forward = 1;
	opt->reverse = 1;
	opt->errors_only = 0;
} /* default_options */

int parse_options(options *opt, int argc, char *argv[]);

int main(int argc, char *argv[])
{
	//int debug_level = ABSOLUTE_SILENCE;//
	int err = NO_ERROR;
	options opt;
	uint32_t rf_idx = -1;
	sam *sd = NULL;
	fastq_data *fqd = NULL;		/* fastq data object */
        fastq_options *fqo = NULL;	/* fastq options object */
	char_t const *ref_seq = NULL;	/* reference sequence */


	default_options(&opt);
	if ((err = parse_options(&opt, argc, argv))) {
		if (err > 0)
			exit(mmessage(ERROR_MSG, INVALID_CMDLINE, ""));
		else
			exit(0);
	}

	FILE *fp = fopen(opt.sam_file, "r");

	if (!fp)
		exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt.sam_file));

	read_sam(fp, &sd);	/* assumes XY_ENCODING */
	fclose(fp);
	fp = NULL;

	/* open fastq file for later writing */
	if (opt.fastq_file) {
		fp = fopen(opt.fastq_file, "w");
		if (!fp)
			exit(mmessage(ERROR_MSG, FILE_OPEN_ERROR,
							opt.fastq_file));
		mmessage(INFO_MSG, NO_ERROR, "Opened fastq file '%s' for "
						"writing.\n", opt.fastq_file);
	}

	/* read fasta file containing the reference sequences */
	if (opt.fasta_file) {
		if ((err = make_fastq_options(&fqo)))
               		goto CLEAR_AND_EXIT;
		fqo->read_names = 1;
		if ((err = read_fastq(opt.fasta_file, &fqd, fqo)))
	                goto CLEAR_AND_EXIT;
		mmessage(INFO_MSG, NO_ERROR, "Read fasta file '%s'.\n",
							opt.fasta_file);
	}

	/* find sam index of selected reference and sequence if fasta file */
	if (opt.ref_name) {
		int rchar = 0;
		int found = 0;
		for (unsigned int i = 0; i < sd->n_ref; ++i) {
			if (!strcmp(opt.ref_name, &sd->ref_names[rchar])) {
				rf_idx = i;
				found = 1;
				break;
			}
			rchar += strlen(&sd->ref_names[rchar]) + 1;
		}
		if (!found)
			exit(mmessage(ERROR_MSG, INVALID_USER_INPUT, "Bad "
				"reference name: '%s'\n", opt.ref_name));
		mmessage(INFO_MSG, NO_ERROR, "Reference '%s' found at index "
						"%u.\n", opt.ref_name, rf_idx);

		/* extract reference sequence from fasta file */
		if (opt.fasta_file) {
			found = 0;
			char *nam = fqd->names;
			char_t *reads = fqd->reads;
			size_t rnam_len = strlen(opt.ref_name);
			for (unsigned int i = 0; i < fqd->n_reads; ++i) {
				if (rnam_len == fqd->name_lengths[i] &&
					!strncmp(nam, opt.ref_name, rnam_len)) {
					mmessage(INFO_MSG, NO_ERROR, "Loaded "
						"reference '%s' from '%s'.\n",
						opt.ref_name, opt.fasta_file);
					ref_seq = reads;
					found = 1;
					break;
				}
				nam += fqd->name_lengths[i];
				reads += fqd->name_lengths[i];
			}
			if (!found)
				exit(mmessage(ERROR_MSG, INVALID_USER_INPUT,
					"Reference '%s' not found in fasta "
					"file '%s'.\n", opt.ref_name,
							opt.fasta_file));
		}
	}

	for (size_t i = 0; i < sd->n_se; ++i) {
		sam_entry *se = &sd->se[i];

		/* filter unmapped */
		if (opt.filter_unmapped && se->flag >> 2 & 1)
			continue;

		if (opt.print_lengths) {
			printf("%zu\n", se->read->len);
			continue;
		}

		/* filter too short or long */
		if ((opt.min_length && se->read->len < opt.min_length)
			|| (opt.max_length && se->read->len > opt.max_length))
			continue;

		/* filter mapped to wrong reference */
		if (opt.ref_name && se->ref != rf_idx)
			continue;

		/* filter mapped to wrong location */
		if (opt.ref_name && (opt.start_posn || opt.end_posn) &&
			(se->pos < opt.start_posn || se->pos >= opt.end_posn))
			continue;

		/* filter forward or reverse strands */
		if (opt.forward != opt.reverse) {
			if (opt.forward && (se->flag >> 4 & 1UL))
				continue;
			else if (opt.reverse && !(se->flag >> 4 & 1UL))
				continue;
		}

		/* output to fastq file */
		if (opt.fastq_file) {
			fprintf(fp, "@%s\n", se->name);
			fwrite_nuc_segment(fp, se->read, XY_ENCODING, 0,
								se->read->len);
			fprintf(fp, "\n+\n");
			fwrite_qual_sequence(fp, se->qual);
			fprintf(fp, "\n");
			continue;
		}

		/* output error information */
		output_error_data(stdout, se, ref_seq, NAN, opt.errors_only);
	}
CLEAR_AND_EXIT:
	return err;
} /* main */

/**
 * Parse command line options and arguments.
 *
 * @param opt	options object
 * @param argc	number of command-line arguments
 * @param argv	command-line arguments
 * @return	error status
 */
int parse_options(options *opt, int argc, char *argv[])
{
	int err = NO_ERROR;
	int argv_idx = 1;

	if (argc < 2) {
		fprint_usage(stderr, argv[0], (void *)opt);
		return 1;
	}
	while (argv_idx < argc) {
		int j = 0;
		int len = strlen(argv[argv_idx]);
		while (j < len && argv[argv_idx][j] == '-')
			++j;
		if (!j) {
			opt->sam_file = argv[argv_idx++];
			continue;
		}
		char c = argv[argv_idx][j];
		switch (c) {
		case 'e':
			opt->errors_only = 1;
			break;
		case 'f':
			if (!strncmp(&argv[argv_idx][j], "fastq", 5)) {
				opt->fastq_file = argv[++argv_idx];
				mmessage(INFO_MSG, NO_ERROR, "Writing data to "
					"fastq file '%s'.\n", opt->fastq_file);
			} else if (!strncmp(&argv[argv_idx][j], "fasta", 5)) {
				opt->fasta_file = argv[++argv_idx];
				mmessage(INFO_MSG, NO_ERROR, "Reading reference"
					" sequences from fasta file '%s'.\n",
								opt->fasta_file);
			} else if (!strncmp(&argv[argv_idx][j], "fo", 2)) {
				opt->reverse = 0;
				mmessage(INFO_MSG, NO_ERROR, "Discarding reads "
					"aligned to reverse strand.\n");
			}
			break;
		case 'l':
			if (argv_idx + 1 == argc
				|| argv[argv_idx + 1][0] == '-') {
				opt->print_lengths = 1;
				break;
			}
			opt->min_length = atoi(argv[++argv_idx]);
			opt->max_length = atoi(argv[++argv_idx]);
			break;
		case 'i':
			opt->sample_id = atoi(argv[++argv_idx]);
			break;
		case 'r':
			if (!strncmp(&argv[argv_idx][j], "rev", 3)) {
				opt->forward = 0;
				mmessage(INFO_MSG, NO_ERROR, "Discarding reads "
					"aligned to forward strand.\n");
				break;
			}
			opt->ref_name = argv[++argv_idx];
			if (argv_idx + 1 < argc
				&& argv[argv_idx + 1][0] != '-') {
				++argv_idx;
				if (sscanf(argv[argv_idx],
					"%u-%u", &opt->start_posn,
						&opt->end_posn) != 2) {
					err = INVALID_CMD_ARGUMENT;
					usage_error((const char **)argv,
						argv_idx, opt);
					return err;
				}
				if (opt->start_posn)
					--opt->start_posn;
				mmessage(INFO_MSG, NO_ERROR, "Selecting reads "
					"aligned to reference '%s', positions "
					"[%u, %u).\n", opt->ref_name,
					opt->start_posn, opt->end_posn);
			} else {
				mmessage(INFO_MSG, NO_ERROR, "Selecting reads "
					"aligned to reference '%s'\n",
								opt->ref_name);
			}
			break;
		case 'u':
			opt->filter_unmapped = 0;
			break;
		case 'h':
			fprint_usage(stderr, argv[0], (void *)opt);
			return -1;
		default:
			usage_error((const char **)argv, argv_idx, opt);
			return 1;
		}
		++argv_idx;
	}
	return err;
} /* parse_options */

/**
 * Print usage information.
 *
 * @param fp	open file handle
 * @param cmd	name of program
 * @param vopt	options object
 */
void fprint_usage(FILE *fp, char const * const cmd, void *vopt)
{
	options *opt = (options *) vopt;

	fprintf(fp, "NAME\n\t%s - process sam files\n\n", cmd);
	fprintf(fp, "SYNOPSIS\n\t%s [-l <min> <max> --ref <ref_name> "
		"--for|--rev --fasta <fasta_file> --unmapped"
		"--errors|--length|--fastq <fastq_file> -i id] <sam_file>\n",
									cmd);
	fprintf(fp, "\nDESCRIPTION\n");
	fprintf(fp, "\tRead and parse sam files.  By default, outputs one line "
		"per nucleotide read.  Request only errors by selecting "
		"reference and providing fasta file with references (error-only"
		" version does not yet handle multiple references).  Use "
		"--fastq to output data in fastq format.  Select reads by "
		"--length, --forward or --reverse strand, aligning to "
		"--reference, including optional region in --reference.\n");
	fprintf(fp, "\nOPTIONS\n");
	fprintf(fp, "\t--e[rrors]\n\t\tOutput data on errors only.\n");
	fprintf(fp, "\t--fasta <fasta_file>\n\t\tRead fasta input file "
			"'%s' with reference sequences.\n", opt->fasta_file);
	fprintf(fp, "\t--fastq <fastq_file>\n\t\tWrite fastq output in file "
						"'%s'.\n", opt->fastq_file);
	fprintf(fp, "\t--fo[rward]\n\t\tOnly extract reads mapped to forward "
			"strand (DEFAULT: %s).\n", opt->forward ? "yes" : "no");
	fprintf(fp, "\t--i[dentifier] <id>\n\t\tFirst output column is this id "
					"(DEFAULT: %d).\n", opt->sample_id);
	fprintf(fp, "\t--l[ength] <min> <max>\n\t\tRemove reads with length "
		"outside range (DEFAULT: %d-%d).\n", opt->min_length,
							opt->max_length);
	fprintf(fp, "\t--length\n\t\tOutput read lengths, one per line "
			"(DEFAULT: %s).\n", opt->print_lengths?"yes":"no");
	fprintf(fp, "\t--ref[erence] <ref_name> [<start>-<end>]\n\t\tRemove "
		"reads not aligning to reference '%s' in region %u to %u.\n",
				opt->ref_name, opt->start_posn, opt->end_posn);
	fprintf(fp, "\t--rev[erse]\n\t\tOnly extract reads mapped to reverse "
			"strand (DEFAULT: %s).\n", opt->reverse ? "yes" : "no");
	fprintf(fp, "\t--u[nmapped]\n\t\tRemove unmapped reads.\n");
	fprintf(fp, "\t<sam_file>\n\t\tRequired sam file to parse.\n");
} /* fprint_usage */
