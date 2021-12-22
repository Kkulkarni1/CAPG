/**
 * @file effici.h
 * @author Yudi Zhang
 *
 * Header file for effici
 *
 */

#ifndef EFFICI_H
#define EFFICI_H

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
#include "cluster_lloyds.h"
#include "myfun.h"

typedef void (*fun)(void *, data_t **, int, int, int, unsigned int *, unsigned int *, int);
typedef int (*fun_iter)(void *, data_t **, int, int, int, unsigned int *, unsigned int *);
typedef double (*compute_crit)(void *, unsigned int *, double *, unsigned int, unsigned int);
typedef void (*fun_init)(void *, data_t **, int, int, int, unsigned int *, unsigned int *, double *);

int fastq_lloyds_efficient_step1 (void *d, data_t **seeds, int K, int n, int p, unsigned int *nclass, unsigned int *ic1, int Init);
void fastq_lloyds_efficient_step2 (void *d, data_t **seeds, int K, int n, int p, unsigned int *nclass, unsigned int *ic1, int Init);
int fastq_macqueen_iter(void *d, data_t **seeds, int K, int n, int p, unsigned int *nclass, unsigned int *ic1);
void fastq_macqueen_ini(void *d, data_t **seeds, int K, int n, int p, unsigned int *nclass, unsigned int *ic1);
void hw_init(void *d, data_t **seeds, int K, int n, int p, unsigned int *nclass, unsigned int *ic, double *cost);
void optra_haplotype(void *d, data_t **seeds, int K, int n, int p, unsigned int *nclass, unsigned int *ic, double *cost, int i, int Is_Live, unsigned int *indx);

double compute_cost(void *d, unsigned int *ic, double *criterion, unsigned int K, unsigned int n);

double cluster_lloyds2 (void *d, data_t **seeds, unsigned int *nclass, unsigned int *ic1, unsigned int n, unsigned int p, unsigned int K, unsigned int max_iter, double *cost, int *ifault, unsigned int *iter, func1 step1, fun step2, compute_crit compute_sum_cost);

double cluster_macqueen(void *d, unsigned int n, unsigned int p, data_t **seeds,
			unsigned int K, unsigned int *ic1, unsigned int *nclass,
			unsigned int max_iter, double *cost, int *ifault,
			unsigned int *iter, func2 ini, fun_iter itr, compute_crit compute_sum_cost);

double cluster_hw (void *d, data_t **seeds, unsigned int *nclass, unsigned int *ic, unsigned int n, unsigned int p, unsigned int K, unsigned int max_iter, double *cost, int *ifault, unsigned int *iter, fun_init init);

#endif /* effici_h */
