/**
 * @file cluster_lloyds.h
 * @author Yudi Zhang
 *
 * Header file for culster_lloyds
 *
 */

#ifndef CLUSTER_LLOYDS_H
#define CLUSTER_LLOYDS_H

#include <float.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <curses.h>
#include <math.h>

#include "constants.h"
#include "haplotype_data.h"
#include "haplotype_option.h"

#define LOG3 1.098612288668109782108

typedef int (*func1)(void *, data_t **, int, int, int, unsigned int *, unsigned int *, int);
typedef void (*func2)(void *, data_t **, int, int, int, unsigned int *, unsigned int *);
typedef double (*compute_rule)(void *d, data_t **, unsigned int *ic, double *, unsigned int K, unsigned int n, unsigned int p);

data_t *find_unique(data_t *mat, int n, int *n_unique);
void compute_pij(void *d, double **prob_t, double **prob_f);
int compare (const void * a, const void * b);

double compute_criterion_hap(void *d, data_t **seeds, unsigned int *ic, double *criterion, unsigned int K, unsigned int n, unsigned int p);
int fastq_lloyds_step1 (void *d, data_t **seeds, int K, int n, int p, unsigned int *nclass, unsigned int *ic1, int Is_Init);
void fastq_lloyds_step2 (void *d, data_t **seeds, int K, int n, int p, unsigned int *nclass, unsigned int *ic1);

double cluster_lloyds(void *d, data_t **seeds, unsigned int *nclass, unsigned int *ic1, unsigned int n, unsigned int p, unsigned int K, unsigned int max_iter, double *cost, int *ifault, unsigned int *iter, func1 step1, func2 step2, compute_rule compute_sum_cost);

#endif /* cluster_lloyds_h */




