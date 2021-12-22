/**
 * @file em.c
 * @author Karin S. Dorman
 *
 * [TODO] Model nucleotide distribution of ambiguous primer nucleotides?
 */
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "error.h"
#include "trim.h"

int em(data *dat, model *mod, options *opt) {
	double rdelta;
	mod->iter = 0;
	mod->ll = -INFINITY;
	do {
		/* E step */
		mod->nll = e_step(dat, mod, opt);
/* debugging code
		if (mod->nll < mod->ll)
			return mmessage(INFO_MSG, NO_ERROR, "%5d: %15.3f (>=" 
				"%5.3f)****\n", mod->iter, mod->nll, mod->ll);
*/

		/* M step */
		m_step(dat, mod, opt);


		rdelta = (mod->ll - mod->nll)/mod->ll;
/* debugging code
		if (mod->nll < mod->ll)
			return mmessage(INFO_MSG, NO_ERROR, "%5d: %15.3f "
				"(%5.3e)****\n", mod->iter, mod->nll, rdelta);
*/
		if (opt->info > SILENT)
			mmessage(INFO_MSG, NO_ERROR, "%5d: %15.3f (%5.3e)\n",
				mod->iter, mod->nll, rdelta);

		/* store latest parameter updates in model::n* variables */
		param_update(mod, opt, DEFAULT);

		mod->iter++;

		/* terminate iterations */
		if (mod->iter >= opt->n_iter)
			return mmessage(WARNING_MSG, EM_EXCEED_MAX_ITERATIONS,
				"%s: %d\n",
				em_error_message(EM_EXCEED_MAX_ITERATIONS),
				opt->n_iter);

		if (isfinite(mod->ll) && rdelta < opt->epsilon)
			return NO_ERROR;
	} while (1);
} /* em */

int param_update(model *mod, options *opt, int best) {
	if (best) {
		mod->best_ll = mod->ll;
		memcpy(mod->best_pi, mod->pi, (opt->n_indices + 1)
					* sizeof *mod->pi);
		memcpy(mod->best_delta, mod->delta, (opt->primer_len
				+ opt->n_indices - 1) * sizeof *mod->delta);
		memcpy(mod->best_eta, mod->eta, NUM_NUCLEOTIDES
						* sizeof *mod->eta);
		memcpy(mod->best_beta, mod->beta, NUM_NUCLEOTIDES
					* sizeof *mod->beta);
		memcpy(mod->best_gamma, mod->gamma, NUM_NUCLEOTIDES
					* NUM_NUCLEOTIDES * sizeof *mod->gamma);
	} else {
		mod->ll = mod->nll;
		memcpy(mod->pi, mod->npi, (opt->n_indices + 1)
					* sizeof *mod->pi);
		memcpy(mod->delta, mod->ndelta, (opt->primer_len
				+ opt->n_indices - 1) * sizeof *mod->delta);
		memcpy(mod->eta, mod->neta, NUM_NUCLEOTIDES * sizeof *mod->eta);
		memcpy(mod->beta, mod->nbeta, NUM_NUCLEOTIDES
					* sizeof *mod->beta);
		memcpy(mod->gamma, mod->ngamma, NUM_NUCLEOTIDES
					* NUM_NUCLEOTIDES * sizeof *mod->gamma);
	}
	return 0;
} /* param_update */

double e_step(data *dat, model *mod, options *opt) {
	int fxn_debug = SILENT;//mod->iter > 0;//ABSOLUTE_SILENCE;//QUIET;//
	unsigned int i, j, k;
	int l;
	size_t m;
	double sum, tmp, ll = 0., max;

	/* [TODO] clean up the following code, which is too deeply nested */

	/* for each read */
	for (i = 0; i < dat->fqd->n_reads; ++i) {
		max = -INFINITY;	/* to hold max_k e_{ik} */

		m = i * opt->n_states;

		/* consider possible locations of primer (or absence of) */
		for (k = 0; k <= opt->n_indices; ++k) {

			/* probability of primer at this position */
			mod->eik[m] = log(mod->pi[k]);

			if (!isfinite(mod->eik[m])) {
				fprintf(stderr, "eik[%u][%u] = %f because pi[%u] = %e\n", i, k, mod->eik[m], k, mod->pi[k]);
				exit(0);
			}

			/* probability of no primer */
			if (k == opt->n_indices) {

				/* every position in read is iid background */
				for (j = 0; j < read_length(dat->fqd, i); ++j)
					mod->eik[m] += log(mod->eta[nuc_idx(dat, i, j)]);
/* DEBUGGING */
				if (!isfinite(mod->eik[m])) {
					fprintf(stderr, "eik[%u][%u] = %f because pi[%u] = %e\n", i, k, mod->eik[m], k, mod->pi[k]);
					exit(0);
				}

				break;

			}

			first_sequence(opt->cprimer, opt->primer, opt->primer_len);
			/* [TODO] assumes all primers are equally likely: estimate them? */
			mod->eik[m] -= log(opt->n_primers);

			do {
				l = -k;

				/* for every position in read */
				for (j = 0; j < read_length(dat->fqd, i); ++j) {

					if (fxn_debug > SILENT)
						fprintf(stderr, "%s:%d: read %u, index %u, site %u", __func__, __LINE__, i, k, j);

					/* within the primer: model errors */
					if (l >= 0 && (unsigned int) l < opt->primer_len) {

						/* barcode */
						if (popcnt[(int)opt->primer[l]] == NUM_NUCLEOTIDES)
							mod->eik[m] += log(mod->beta[nuc_idx(dat, i, j)]);

						/* no read error */
						else if (!nuccmp(opt->cprimer[l], dat, i, j))
							mod->eik[m] += log(mod->delta[j]);

						/* read error */
						else
							mod->eik[m] += log(1 - mod->delta[j]) 
								+ log(mod->gamma[err_idx(opt->cprimer[l], dat, i, j)]);

/* save time as this no longer triggers
						if (!isfinite(mod->eik[m]))
							exit(mmessage(INFO_MSG, NO_ERROR,
								"eik[%u][%u] = %f at site %u (primer site %u), iteration %u (p=%c r=%c %e 1-d=%e g=%e)\n",
								i, k, mod->eik[m], j, l, mod->iter, iupac_to_char[(int)opt->cprimer[l]],
								xy_to_char[(int)dat->reads[i][j]], -log(opt->n_primers), 1-mod->delta[j],
								mod->gamma[err_idx(opt->cprimer[l], dat, i, j)]));
*/
					} else if (l < 0) {
						mod->eik[m] += log(mod->beta[nuc_idx(dat, i, j)]);
					} else {
						mod->eik[m] += log(mod->eta[nuc_idx(dat, i, j)]);
/* save time as this no longer triggers
						if (!isfinite(mod->eik[m]))
							exit(mmessage(INFO_MSG, NO_ERROR,
								"eik[%u][%u] = %f at site %u (primer site %u), iteration %u (r=%c %e)\n",
								i, k, mod->eik[m], j, l, mod->iter,
								nuc_char(dat, i, j), mod->eta[nuc_idx(dat, i, j)]));
*/
					}
					if (fxn_debug > SILENT)
						fprintf(stderr, " = %e\n", mod->eik[m]);
					l++;
				}
				if (max < mod->eik[m])
					max = mod->eik[m];
				m++;
			} while (next_sequence(opt->cprimer, opt->primer, opt->primer_len));

		}
		if (max < mod->eik[m])
			max = mod->eik[m];


		/* normalize */
		sum = 0.;
		m = i * opt->n_states;
		for (k = 0; k < opt->n_states; ++k) {
//fprintf(stderr, "Read %u, state %u: %f", i, k, mod->eik[m]);
			mod->eik[m] = exp(mod->eik[m] - max);
			sum += mod->eik[m];
			++m;
//fprintf(stderr, " (%f; %f)\n", sum, max);
		}

		if (!isfinite(sum))
			exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "Read %u non-finite likelihood\n", i));

		/* compute previously log likelihood */
		ll += log(sum) + max;

		tmp = 0.;
		m = i * opt->n_states;
		for (k = 0; k < opt->n_states; ++k) {
			mod->eik[m] /= sum;
			tmp += mod->eik[m];
			++m;
		}
/* save time as this no longer triggers
*/
		if (fabs(tmp - 1.) > 1e-6) {
			mmessage(INFO_MSG, NO_ERROR, "eik[%u] sum != 1.0 (%e)\n", i, tmp = 1.);
			exit(0);
		}
	}
	return ll;
} /* e_step */

void m_step(data *dat, model *mod, options *opt) {
	size_t j1, j2;
	unsigned int i, j, k, m;
	int l;
	double sum, tmp1, tmp2;
	double epsilon = 1e-8;
	double *dptr;

	if (opt->estimate & ESTIMATE_PI)
		for (k = 0; k <= opt->n_indices; ++k)
			mod->npi[k] = epsilon;

	/* use pseudocounts to avoid boundaries */
	for (j = 0; j < opt->primer_len + opt->n_indices - 1; ++j)
		mod->ndelta[j] = 0.;

	dptr = mod->ngamma;
	for (i = 0; i < NUM_NUCLEOTIDES; ++i)
		for (j = 0; j < NUM_NUCLEOTIDES; ++j) {
			*dptr = 0.;
			dptr++;
		}
	
	for (i = 0; i < NUM_NUCLEOTIDES; ++i)
		mod->nbeta[i] = 0.;

	for (i = 0; i < NUM_NUCLEOTIDES; ++i)
		mod->neta[i] = 0.;

	for (i = 0; i < dat->fqd->n_reads; ++i) {

		/* update mixing proportions */
		sum = 0;
		j1 = i * opt->n_states;

		/* for all primer positions */
		for (k = 0; k < opt->n_indices; ++k) {
			/* and all possible primers */
			for (j2 = 0; j2 < opt->n_primers; ++j2, ++j1) {
//if (!k) fprintf(stderr, "Read %u assigned to group %u wp %f\n", i, k, mod->eik[j1]);
				if (opt->estimate & ESTIMATE_PI)
					mod->npi[k] += mod->eik[j1];
				sum += mod->eik[j1];

/* TESTING CODE
				if (!isfinite(mod->npi[k])) {
					fprintf(stderr, "pi[%u] += eik[%u][%u] = %f\n", k, i, k, mod->eik[j1]);
					exit(EXIT_FAILURE);
				}
*/
			}
		}

		/* and no primer */
		sum += mod->eik[j1];
		if (opt->estimate & ESTIMATE_PI)
			mod->npi[opt->n_indices] += mod->eik[j1];
/* save time as this no longer triggers
*/
		if (fabs(sum - 1.) > 1e-6)
			exit(mmessage(INFO_MSG, NO_ERROR, "eik[%u] sum != 1.0 (%f)\n", i, sum - 1.));
/* TESTING CODE
		if (!isfinite(mod->npi[opt->n_indices])) {
			fprintf(stderr, "pi[%u] += eik[%u][%u] = %f\n", opt->n_indices, i, opt->n_indices, mod->eik[j1]);
			exit(EXIT_FAILURE);
		}
*/

		/* primer site at one of these positions */
		m = i * opt->n_states;
		for (k = 0; k < opt->n_indices; ++k) {
			first_sequence(opt->cprimer, opt->primer, opt->primer_len);
			do {
				l = -k;
				for (j = 0; j < read_length(dat->fqd, i); ++j) {
					/* within primer */
					if (l >= 0 && (unsigned int) l < opt->primer_len) {
						/* informative site */
						if (popcnt[(int)opt->primer[l]] < NUM_NUCLEOTIDES) {
							if (!nuccmp(opt->cprimer[l], dat, i, j)) {
//if (!j) fprintf(stderr, "Read %u match at position %u %f %f\n", i, j, mod->eik[m], mod->ndelta[j]);
								mod->ndelta[j] += mod->eik[m] / popcnt[(int)opt->primer[l]] ;
								if (!isfinite(mod->ndelta[j]))
									exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "delta[%u] += eik[%u][%u] = %f\n",
										j, i, k, mod->eik[m]));
							} else {
//if (!j) fprintf(stderr, "Read %u nonmatch at position %u %f %f\n", i, j, mod->eik[m], mod->ndelta[j]);
								mod->ngamma[err_idx(opt->cprimer[l], dat, i, j)] += mod->eik[m] / popcnt[(int)opt->primer[l]];
								if (!isfinite(mod->ngamma[err_idx(opt->cprimer[l], dat, i, j)]))
									exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "gamma[%c][%c] += eik[%u][%u] = %f (%u)\n",
										iupac_to_char[(int)opt->cprimer[l]], nuc_char(dat, i, j), i, k, mod->eik[m], j));
							}
						/* barcode site */
						} else {
							mod->nbeta[nuc_idx(dat, i, j)] += mod->eik[m];
/* debugging
							if (!isfinite(mod->nbeta[nuc_idx(dat, i, j)]))
								exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "beta[%c] += eik[%u][%u] = %f (%u)\n",
									nuc_char(dat, i, j), i, k, mod->eik[m], j));
*/
						}
					/* outside primer */
					} else if (l < 0) {
						mod->nbeta[nuc_idx(dat, i, j)] += mod->eik[m];
					} else {
						mod->neta[nuc_idx(dat, i, j)] += mod->eik[m];
/* debugging
						if (!isfinite(mod->neta[nuc_idx(dat, i, j)]))
							exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "eta[%c] += eik[%u][%u] = %f (%u)\n",
								nuc_char(dat, i, j), i, k, mod->eik[m], j));
*/
					}
					l++;
				}
				m++;
			} while (next_sequence(opt->cprimer, opt->primer, opt->primer_len));
		}

		/* no primer site: read is iid~ eta */
		for (j = 0; j < read_length(dat->fqd, i); ++j) {
			mod->neta[nuc_idx(dat, i, j)] += mod->eik[m];
			if (!isfinite(mod->neta[nuc_idx(dat, i, j)]))
				exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "eta[%c] += eik[%u][%u] = %f (%u)\n",
					nuc_char(dat, i, j), i, k, mod->eik[m], j));
		}
	}

	if (opt->estimate & ESTIMATE_PI) {
		tmp1 = 0;
		for (k = 0; k <= opt->n_indices; ++k) {
			mod->npi[k] /= dat->fqd->n_reads + (opt->n_indices + 1)*epsilon;
			tmp1 += mod->npi[k];
/* save time as this no longer triggers
			if (mod->npi[k] < 0 || mod->npi[k] > 1) {
				mmessage(ERROR_MSG, INTERNAL_ERROR, "pi[%u]: %f not in [0,1]\n", k, mod->npi[k]);
				exit(EXIT_FAILURE);
			}
*/
		}

		if (fabs(tmp1 - 1.0) > 1e-6)
			exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "pi do not sum to 1 (%f)\n", tmp1 - 1.0));
	}

	if (opt->n_indices > opt->primer_len)
		exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "your programmer did "
			"not anticipate number of possible primer positions "
			"(%u) might be more than primer length (%u)\n",
			opt->n_indices, opt->primer_len));

	tmp1 = tmp2 = 0;
	for (j = 0; j < opt->primer_len + opt->n_indices - 1; ++j) {

		if (j < opt->n_indices) {
			tmp1 += mod->npi[j];// + (opt->n_indices + 1)*epsilon;
//fprintf(stderr, "Divide by %f (num=%f tmp1=%f nr=%f %e", tmp1 * (dat->fqd->n_reads + 2*epsilon), mod->ndelta[j], tmp1, (dat->fqd->n_reads + 2*epsilon), tmp1*(dat->fqd->n_reads + 2*epsilon) - tmp1*dat->fqd->n_reads);
			mod->ndelta[j] /= tmp1 * dat->fqd->n_reads;
//fprintf(stderr, " pi=%f -> %f %u)\n", mod->npi[j], mod->ndelta[j], j);
		} else if (j >= opt->primer_len) {
			tmp1 -= mod->npi[j - opt->primer_len];// + (opt->n_indices + 1)*epsilon;
			mod->ndelta[j] /= tmp1 * dat->fqd->n_reads;
//fprintf(stderr, "Divide by %f (%f %f %f %u)\n", tmp1 * (dat->fqd->n_reads + 2*epsilon), mod->ndelta[j], tmp1, mod->npi[j - opt->primer_len], j);
		} else {
			mod->ndelta[j] /= (1 - mod->npi[opt->n_indices])
				* dat->fqd->n_reads;
//fprintf(stderr, "Divide by %f (%f %u)\n", (1 - mod->npi[opt->n_indices]) * dat->fqd->n_reads, 1 - mod->npi[opt->n_indices], dat->fqd->n_reads);
		}
//fprintf(stderr, "delta[%u] - 1. = %e\n", j, mod->ndelta[j] - 1.);

		if (mod->ndelta[j] < epsilon) {
			if (mod->ndelta[j] < -epsilon)
				exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "delta[%u] < 0 (%e)", j, mod->ndelta[j] - 0.));
			mod->ndelta[j] = epsilon;
		}
		if (mod->ndelta[j] > 1 - epsilon) {
			if (mod->ndelta[j] > 1 + epsilon)
				exit(mmessage(ERROR_MSG, INTERNAL_ERROR, "delta[%u] > 1 (%e)", j, mod->ndelta[j] - 1.));
			mod->ndelta[j] = 1. - epsilon;
		}

	}
	tmp1 = tmp2 = 0;
	for (i = 0; i < NUM_NUCLEOTIDES; ++i) {
		tmp1 += mod->neta[i];
		tmp2 += mod->nbeta[i];
		sum = 0;
		for (j = 0; j < NUM_NUCLEOTIDES; ++j)
			if (i != j)
				sum += mod->ngamma[i*NUM_NUCLEOTIDES + j];
		for (j = 0; j < NUM_NUCLEOTIDES; ++j)
			if (i != j) {
				mod->ngamma[i*NUM_NUCLEOTIDES + j] /= sum;
/* no longer triggers anymore
*/
				if (mod->ngamma[i*NUM_NUCLEOTIDES + j] <= 0 || mod->ngamma[i*NUM_NUCLEOTIDES + j] >= 1) {
					message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG, INTERNAL_ERROR, "gamma[%u][%u]: %f not in (0,1)\n", i, j, mod->ngamma[i*NUM_NUCLEOTIDES + j]);
					exit(EXIT_FAILURE);
				}
			}
	}

	for (i = 0; i < NUM_NUCLEOTIDES; ++i) {
		mod->neta[i] /= tmp1;
		mod->nbeta[i] /= tmp2;
	}

} /* m_step */

void assign_index(data *dat, model *mod, options *opt, int which) {
	double max;
	size_t i, k, l = 0;

	for (i = 0; i < dat->fqd->n_reads; ++i) {
		max = -INFINITY;
		for (k = 0; k <= opt->n_indices; ++k)
			if (mod->eik[i * opt->n_states + k] > max) {
				max = mod->eik[i * opt->n_states + k];
				l = k;
			}

		/* three possible places to store result */
		if (!which)
			dat->index[i] = l;
		else if (which == BEST)
			dat->best_index[i] = l;
		else
			dat->optimal_index[i] = l;
	}
} /* assign_index */

char *em_error_message(int err_no) {
	if (err_no == EM_EXCEED_MAX_ITERATIONS)
		return "Exceed maximum allowed iterations";
	else if (err_no == EM_ASCENT_VIOLATION)
		return "Log likelihood decline";
	else if (err_no)
		return "unknown error";
	else
		return "";
} /* em_error_message */
