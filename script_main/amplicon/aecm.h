/**
 * @file aecm.h
 * @author Karin S. Dorman
 *
 * Header file for aecm enums and extern functions.
 */

#ifndef __H_AECM__
#define __H_AECM__

#include "data.h"
#include "model_options.h"
#include "model.h"
#include "run_ampliclust.h"

/**
 * EM cycle.
 */
enum {
	FIRST,			/*<! first E step */
	SECOND,			/*<! second E step */
	PRINT_TRUE_READ_LLS,	/*<! print read lls, true parameters */
	PRINT_ESTD_READ_LLS	/*<! print read lls, true parameters */
};

/**
 * AECM errors.
 */
enum {
	AECM_NO_ERROR,			/*<! no AECM errors */
	AECM_EXCEED_MAX_ITERATIONS,	/*<! over maximum allowed iterations */
	AECM_ASCENT_VIOLATION,		/*<! AECM log likelihood decline */
};

/**
 * Which copy of parameters.
 */
enum {
	DEFAULT_COPY,
	NEW_COPY,
	BEST_COPY,
	OPTIMAL_COPY
};

int aecm(data *dat, model *mod, model_options *opt);
double e_step(data *dat, model *mod, model_options *opt, int second);
void m_other(data *dat, model *mod, model_options *opt);
int param_update(model *mod, data *dat, model_options *opt, int which);
int param_copy(model *to_model, model *from_mod, data *dat, model_options *opt);
const char *aecm_error_message(int err_no);
int m_JC69(unsigned char * hap, unsigned char * anc, double *dist, unsigned int K, unsigned int len);
double e_JC69(unsigned char * hap, unsigned char * anc, double *dist, unsigned int K, unsigned int len);
void haplotype_delta_hamming_matrix(model *mod, data *dat, model_options *opt, int, int);

#endif
