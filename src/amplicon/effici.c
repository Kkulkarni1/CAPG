/**
 * @file effici.c
 * @author Yudi Zhang
 *
 * Functions for efficient algorithm applied to fastq data assuming literal
 * quality scores and uniform substitution probabilties.  It also contains
 * the general interface to Lloyd's algorithm that will presumably one
 * day replace the code in kmodes.c::kmodes_lloyd().
 *
 */

#include "effici.h"

/**
 * Compute cost criterion.
 *
 * @param d           data structure
 * @param ic          cluster assignments (n x 1)
 * @param criterion   cost per cluster
 * @param K           number of clusters
 * @param n           number of observations
 * @return            return sum of cluster costs
 */
double compute_cost(void *d, unsigned int *ic, double *criterion, unsigned int K, unsigned int n) {
	
	data *dp = (data *)d;
	
	SETZERO_VECTOR(criterion, K);
	
	for (int i = 0; i < n; ++i)
		criterion[ic[i]] += dp->v_ik[i][ic[i]];
	
	double sum = 0;
	for (int k = 0; k < K; ++k)
		sum += criterion[k];
	
	PRINT_VECTOR(criterion, K);
	
	return sum;
	
} /* compute_cost */

/**
 * Step 2 for lloyd applied to fastq data using uniform substitution
 * probabilities and literal quality scores.  This step updates the
 * haplotype nucleotides.
 *
 * @param d        void pointer
 * @param seeds    seeds
 * @param K        no. cluster
 * @param n        no. observation
 * @param p        no. coordinate
 * @param nclass   count in clsuters (1 x k)
 * @param ic1      cluster assignments (1 x n)
 * @param Init     if it is at the initial phase
 */
void fastq_lloyds_efficient_step2 (void *d, data_t **seeds, int K, int n, int p, unsigned int *nclass, unsigned int *ic1, int Init) {
	
	data *dp = (data *)d;
	unsigned int i, k, j, l, len_index;
	double e_kjN_diff;
	
	//Update e_kjn when doing the iteration
	if (!Init)
		for (i = 0; i < n; ++i) {
			if (dp->last_assign[i] == ic1[i])
				continue;
			for (j = 0; j < p; ++j) {
				for (l = 0; l < dp->tot_n_categories; ++l) {
					e_kjN_diff = dp->dmat[i][j] == dp->categories[l] ? dp->prob_t[i][j] : dp->prob_f[i][j];
					dp->e_kjN[dp->last_assign[i]][j][l] -= e_kjN_diff;
					dp->e_kjN[ic1[i]][j][l] += e_kjN_diff;
				}
			}
		}
	
	/* Store e_kjN; Second step: recompute the center */
	for (k = 0; k < K; ++k) {
		len_index = 0;
		
		if (Init)
			for (i = 0; i < n; i++)
				if (ic1[i] == k)
					dp->index[len_index++] = i;
		
		for (j = 0; j < p; ++j) {
			/* Store the current center for the following comparision */
			dp->last_cent[k][j] = seeds[k][j];
			
			double max_llj = -INFINITY;
			int max_id = -1;
			
			for (l = 0; l < dp->tot_n_categories; ++l) {
				if (Init)
					for (i = 0; i < len_index; ++i)
						dp->e_kjN[k][j][l] += dp->dmat[dp->index[i]][j] == dp->categories[l] ? dp->prob_t[dp->index[i]][j] : dp->prob_f[dp->index[i]][j];
				
				if (dp->e_kjN[k][j][l] > max_llj) {  /*find the largest e_kjN[k][j][l] for each site*/
					max_llj = dp->e_kjN[k][j][l];
					max_id = l;
				}
			}
			seeds[k][j] = dp->categories[max_id]; //Compute new Hk
		}
	}
}/* fastq_lloyds_efficient_step2 */

/**
 * Function to carry out step 1 for fastq data using uniform substitution
 * probabilities and literal quality scores.  This step reassigns observations
 * to clusters.
 *
 * @param d        void pointer to the data structure
 * @param seeds    seeds
 * @param K        no. clusters
 * @param n        no. observations
 * @param p        no. coordinates
 * @param nclass    cluster counts (1 x K)
 * @param ic1        cluster assignments (1 x n)
 * @param Init    is this an initialization step
 * @return keep_going
 */
int fastq_lloyds_efficient_step1 (void *d, data_t **seeds, int K, int n, int p, unsigned int *nclass, unsigned int *ic1, int Init) {
	
	unsigned int i, k, j;
	int keep_going = 0;
	double v_ik_diff;
	data *dp = (data *)d;
	
	// Ressiagn the reads to new Hk, compare if each site of the centers has changed or not
	for (i = 0; i < n; ++i) {        /* for each observation */
		if (!Init)
			dp->last_assign[i] = ic1[i];      // Store the last cluster assignment
		
		double max = -INFINITY;
		int max_id = -1;
		
		for (k = 0; k < K; ++k) {    /* compare to each haplotype */
			for (j = 0; j < p; ++j) {
				
				if (Init)
					dp->v_ik[i][k] += dp->dmat[i][j] == seeds[k][j] ? dp->prob_t[i][j] : dp->prob_f[i][j];
				
				else if (dp->last_cent[k][j] != seeds[k][j]) {
					
					v_ik_diff = dp->prob_t[i][j] - dp->prob_f[i][j];
					
					if (dp->dmat[i][j] == seeds[k][j])
						dp->v_ik[i][k] += v_ik_diff;
					else if (dp->dmat[i][j] == dp->last_cent[k][j])
						dp->v_ik[i][k] -= v_ik_diff;
				}
			}
			if (dp->v_ik[i][k] > max) { // Find the biggest v_ik for each read
				max = dp->v_ik[i][k];
				max_id = k;
			}
		}
		if (Init) {
			nclass[max_id]++;
			ic1[i] = max_id;
		} else if (dp->last_assign[i] != max_id) {
			keep_going = 1;
			nclass[max_id]++;
			nclass[dp->last_assign[i]]--;
			ic1[i] = max_id;
		}
	}
	
	return keep_going;
} /* fastq_lloyds_efficient_step1 */

/**
 * Effificient version of Lloyd's algorithm. Alternate between updating cluster assignments and cluster centers until convergence.
 *
 * @param d void pointer data structure: defined by user
 * @param seeds seeds
 * @param nclass count in clsuters (1 x k)
 * @param ic1 cluster assignments (1 x n)
 * @param n number of observations
 * @param p number of coordinates
 * @param K number of clusters
 * @param max_iter maximum allowed iterations
 * @param cost cost per cluster
 * @param ifault pointer to error status
 * @param iter pointer to number of iterations
 * @param step1 function pointer: user defined way to carry out the 1st step of lloyds
 * @param step2 function pointer: user defined way to carry out the 2nd step of lloyds
 * @return optimized criterion
 */
double cluster_lloyds2 (void *d, data_t **seeds, unsigned int *nclass, unsigned int *ic1, unsigned int n, unsigned int p, unsigned int K, unsigned int max_iter, double *cost, int *ifault, unsigned int *iter, func1 step1, fun step2, compute_crit compute_sum_cost)
{
	
	*ifault = KMODES_NO_ERROR;
	unsigned int m = 0;
	int keep_going = 1;
	
	/* user requests too few or too many clusters */
	if (K < 1 || K > n) {
		*ifault = KMODES_CALLER_INPUT_ERROR;
		/* [KSD, TODO] Let caller exit gracefully. */
		mmessage(ERROR_MSG, INVALID_USER_INPUT,
			 "Request too few or many cluasters");
		exit(EXIT_FAILURE); // return INFINITY;
	}
	
	/* initialize cluster counts and cost */
	SETZERO_VECTOR(nclass, K);
	
	step1(d, seeds, K, n, p, nclass, ic1, 1);
	
	step2(d, seeds, K, n, p, nclass, ic1, 1);
	
	while (m++ < max_iter && keep_going) {
		
		keep_going = step1(d, seeds, K, n, p, nclass, ic1, 0);
		
		if (!keep_going)
			break;
		
		step2(d, seeds, K, n, p, nclass, ic1, 0);
	}
	
	if (m > max_iter)
		*ifault = KMODES_EXCEED_ITER_WARNING;
	
	*iter = m;
	
	return compute_sum_cost(d, ic1, cost, K, n);
	
} /* fastq_lloyds_efficient */

/**
 * Modified MacQueen's algorithm. Initialization step is the same as Lloyd's efficient.
 
 * @param d void pointer data structure: defined by user
 * @param seeds seeds
 * @param nclass count in clsuters (1 x k)
 * @param ic1 cluster assignments (1 x n)
 * @param n number of observations
 * @param p number of coordinates
 * @param K number of clusters
 * @param max_iter maximum allowed iterations
 * @param cost cost per cluster
 * @param ifault pointer to error status
 * @param iter pointer to number of iterations
 * @return optimized criterion
 */
double cluster_macqueen(void *d, unsigned int n, unsigned int p, data_t **seeds,
			unsigned int K, unsigned int *ic1, unsigned int *nclass,
			unsigned int max_iter, double *cost, int *ifault,
			unsigned int *iter, func2 ini, fun_iter itr, compute_crit compute_sum_cost) {
	
	*ifault = KMODES_NO_ERROR;
	unsigned int k, m = 0;
	int keep_going = 1;
	
	/* user requests too few or too many clusters */
	if (K < 1 || K > n) {
		*ifault = KMODES_CALLER_INPUT_ERROR;
		mmessage(ERROR_MSG, INVALID_USER_INPUT,
			 "Request too few or many cluasters");
		exit(EXIT_FAILURE); // return INFINITY;
	}
	
	/* find the initial cluster assignment */
	SETZERO_VECTOR(nclass, K);
	
	ini(d, seeds, K, n, p, nclass, ic1);
	/* Iteration */
	while (m++ < max_iter && keep_going) {
		printf("%d\n", m);
		keep_going = itr(d, seeds, K, n, p, nclass, ic1);
		if (!keep_going)
			break;
	}
	m--;
	
	if (m > max_iter && max_iter)
		*ifault = KMODES_EXCEED_ITER_WARNING;
	
	/* error if initialization routine produces empty cluster */
	for (k = 0; k < K; ++k)
		if (nclass[k] == 0) {
			*ifault = KMODES_NULL_CLUSTER_ERROR;
			return INFINITY;
		}
	
	*iter = m;
	
	return compute_sum_cost(d, ic1, cost, K, n);
}

int fastq_macqueen_iter(void *d, data_t **seeds, int K, int n, int p, unsigned int *nclass, unsigned int *ic1) {
	
	unsigned int i, k, j, l;
	int keep_going = 0;
	double v_ik_diff, e_kjN_diff;
	data *dp = (data *)d;
	
	for (k = 0; k < K; ++k)
		for (j = 0; j < p; ++j)
			dp->last_cent[k][j] = seeds[k][j];
	
	for (i = 0; i < n; ++i) {
		
		dp->last_assign[i] = ic1[i];
		double max = -INFINITY;
		int max_id = -1;
		
		/* Update cluster assignment */
		for (k = 0; k < K; ++k)
			if (dp->v_ik[i][k] > max) { // Find the biggest v_ik for each read
				max = dp->v_ik[i][k];
				max_id = k;
			}
		if (dp->last_assign[i] != max_id) {
			
			keep_going = 1;
			
			/* update cluster counts */
			nclass[max_id]++;
			nclass[dp->last_assign[i]]--;
			ic1[i] = max_id;
			
			for (j = 0; j < p; ++j) {
				
				double max_lst_llj = -INFINITY, max_cur_llj = -INFINITY;
				int max_lst_id = -1, max_cur_id = -1;
				
				for (l = 0; l < dp->tot_n_categories; ++l) {
					
					e_kjN_diff = dp->dmat[i][j] == dp->categories[l] ? dp->prob_t[i][j] : dp->prob_f[i][j];
					dp->e_kjN[dp->last_assign[i]][j][l] -= e_kjN_diff;
					dp->e_kjN[ic1[i]][j][l] += e_kjN_diff;
					
					/*find the largest e_kjN[k][j][l] for each site*/
					if (dp->e_kjN[dp->last_assign[i]][j][l] > max_lst_llj) {
						max_lst_llj = dp->e_kjN[dp->last_assign[i]][j][l];
						max_lst_id = l;
					}
					
					if (dp->e_kjN[ic1[i]][j][l] > max_cur_llj) {
						max_cur_llj = dp->e_kjN[ic1[i]][j][l];
						max_cur_id = l;
					}
				}
				/* Update haplotypes */
				seeds[dp->last_assign[i]][j] = dp->categories[max_lst_id];
				seeds[ic1[i]][j] = dp->categories[max_cur_id];
				
				/* Update vik, the last center should be the one used to compute the last vik, so after interate i from 1 to n, we update the values of last_cent */
				if (dp->last_cent[dp->last_assign[i]][j] != seeds[dp->last_assign[i]][j]) {
					v_ik_diff = dp->prob_t[i][j] - dp->prob_f[i][j];
					if (dp->dmat[i][j] == seeds[dp->last_assign[i]][j])
						dp->v_ik[i][dp->last_assign[i]] += v_ik_diff;
					else if (dp->dmat[i][j] == dp->last_cent[dp->last_assign[i]][j])
						dp->v_ik[i][dp->last_assign[i]] -= v_ik_diff;
				}
				
				if (dp->last_cent[ic1[i]][j] != seeds[ic1[i]][j]) {
					v_ik_diff = dp->prob_t[i][j] - dp->prob_f[i][j];
					if (dp->dmat[i][j] == seeds[ic1[i]][j])
						dp->v_ik[i][ic1[i]] += v_ik_diff;
					else if (dp->dmat[i][j] != seeds[ic1[i]][j])
						dp->v_ik[i][ic1[i]] -= v_ik_diff;
				}
			}
		}
	}
	
	return keep_going;
}

void fastq_macqueen_ini(void *d, data_t **seeds, int K, int n, int p, unsigned int *nclass, unsigned int *ic1) {
	
	data *dp = (data *)d;
	unsigned int i, k, j, l;
	double v_ik_diff;
	
	for (i = 0; i < n; ++i) {
		/* identify closest haplotype */
		double maxll = -INFINITY;
		int max_id = -1;
		for (k = 0; k < K; ++k) {    /* compare to each haplotype */
			for (j = 0; j < p; ++j)
				dp->v_ik[i][k] += dp->dmat[i][j] == seeds[k][j] ? dp->prob_t[i][j] : dp->prob_f[i][j];
			if (dp->v_ik[i][k] > maxll) { // Find the biggest v_ik for each read
				maxll = dp->v_ik[i][k];
				max_id = k;
			}
		}
		nclass[max_id]++;
		ic1[i] = max_id;
		
		/* update log probability of each cluster, position, haplotype nucleotide */
		for (j = 0; j < p; ++j)
			for (l = 0; l < dp->tot_n_categories; ++l)
				dp->e_kjN[ic1[i]][j][l] += dp->dmat[i][j] == dp->categories[l] ? dp->prob_t[i][j] : dp->prob_f[i][j];
	}
	
	/* compute most likely haplotype nucleotides and update vik */
	for (k = 0; k < K; ++k) {
		for (j = 0; j < p; ++j) {
			dp->last_cent[k][j] = seeds[k][j];
			double max_llj = -INFINITY;
			int max_id = -1;
			for (l = 0; l < dp->tot_n_categories; ++l)
				if (dp->e_kjN[k][j][l] > max_llj) {  /*find the largest e_kjN[k][j][l] for each site*/
					max_llj = dp->e_kjN[k][j][l];
					max_id = l;
				}
			seeds[k][j] = dp->categories[max_id];
			if (seeds[k][j] != dp->last_cent[k][j]) {
				for (i = 0; i < n; ++i) {
					v_ik_diff = dp->prob_t[i][j] - dp->prob_f[i][j];
					if (dp->dmat[i][j] == seeds[k][j])
						dp->v_ik[i][k] += v_ik_diff;
					else if (dp->dmat[i][j] == dp->last_cent[k][j])
						dp->v_ik[i][k] -= v_ik_diff;
				}
			}
		}
	}
}

/**
 * Hartigan and Wong's algorithm for clustering haplotype

 * @param d void pointer data structure: defined by user
 * @param seeds seeds
 * @param nclass count in clsuters (1 x k)
 * @param ic cluster assignments (1 x n)
 * @param n number of observations
 * @param p number of coordinates
 * @param K number of clusters
 * @param max_iter maximum allowed iterations
 * @param cost cost per cluster
 * @param ifault pointer to error status
 * @param iter pointer to number of iterations
 * @param init function pointer: user defined way to do initialization
 * @return optimized criterion
 
 **/
double cluster_hw (void *d, data_t **seeds, unsigned int *nclass, unsigned int *ic, unsigned int n, unsigned int p, unsigned int K, unsigned int max_iter, double *cost, int *ifault, unsigned int *iter, fun_init init) {
	
	*ifault = KMODES_NO_ERROR;
	
	unsigned int i, k, len_index, indx = 0;
	
	data *dp = (data *)d;
	
	/* user requests too few or too many clusters */
	if (K < 1 || K > n) {
		*ifault = KMODES_CALLER_INPUT_ERROR;
		/* [KSD, TODO] Let caller exit gracefully. */
		mmessage(ERROR_MSG, INVALID_USER_INPUT,
			 "Request too few or many cluasters");
		exit(EXIT_FAILURE); // return INFINITY;
	}
	
	/* initialize cluster counts and cost */
	SETZERO_VECTOR(nclass, K);
	SETZERO_VECTOR(cost, K);
	
	/* Initialization */
	init(d, seeds, K, n, p, nclass, ic, cost);
	
//	PRINT_VECTOR(nclass, K);
//	PRINT_VECTOR(ic, n);
//	PRINT_VECTOR(cost, K);
//	PRINT_MATRIX(dp->v_ik, n, K);
//	PRINT_MATRIX(seeds, K, p);
//	PRINT_3DARRAY(dp->e_kjN, K, p, dp->tot_n_categories);
	
	/* Optimal transfer stage */
	for (int m = 0; m < max_iter; m++) {
		
		printf("m: %d\n", m);
		
		for(i = 0; i < n; ++i) {
			indx++;
			/* If in the live cluster, operate for all the clusters */
			printf("%d; live_situation:%d\n", i, dp->live[ic[i]]);
			
				if (i < dp->live[ic[i]]) {
					
					for (k = 0; k < K; ++k)
						dp->index[k] = k;
					
					printf("choosed k\n");
					PRINT_VECTOR(dp->index, K);
					
					printf("indx_bf_oper: %d\n", indx);
					optra_haplotype(d, seeds, K, n, p, nclass, ic, cost, i, 1, &indx);
					
				}
				/* If not in the live set, repeat the above steps only for the live sets */
				else {
					len_index = 0;
					
					/* Find the live set */
					for (k = 0; k < K; ++k)
						if (i < dp->live[k])
							dp->index[len_index++] = k;
					dp->index[len_index++] = ic[i];
					
					printf("choosed k\n");
					PRINT_VECTOR(dp->index, len_index);
					
					printf("indx_bf_oper: %d\n", indx);
					optra_haplotype(d, seeds, len_index, n, p, nclass, ic, cost, i, 0, &indx);
				}
				
			}
		
			for (k = 0; k < K; ++k)
				dp->live[k] -= n;
		
		printf("number in each cluster:\n");
		PRINT_VECTOR(nclass, K);
		
		PRINT_3DARRAY(dp->e_kjN, K, p, dp->tot_n_categories);
		
		PRINT_MATRIX(dp->v_ik, n, K);
		
		printf("indx: %d\n", indx);
		if (indx == n)
			break;
		
		
	}
	
	PRINT_VECTOR(nclass, K);
	PRINT_VECTOR(ic, n);
	PRINT_VECTOR(cost, K);
	PRINT_MATRIX(seeds, K, p);
	
	exit(4);
	
	double sum = 0;
	
	for (int k = 0; k < K; ++k)
		sum += cost[k];
	
	return sum;

}


/**
 * Initialization for HW

 * @param d void pointer data structure: defined by user
 * @param seeds seeds
 * @param nclass count in clsuters (1 x k)
 * @param ic cluster assignments (1 x n)
 * @param n number of observations
 * @param p number of coordinates
 * @param K number of clusters
 * @param cost cost per cluster
 
 */
void hw_init(void *d, data_t **seeds, int K, int n, int p, unsigned int *nclass, unsigned int *ic, double *cost) {
	
	data *dp = (data *)d;
	unsigned int i, k, j, l;
	double max, v_ik_diff;
	int max_id;
	
	/* Find the biggest and the second biggest v_ik for each read */
	for (i = 0; i < n; i++) {
		
		max = -INFINITY;
		max_id = -1;
		
		for (k = 0; k < K; ++k) {    /* compare to each haplotype */
			for (j = 0; j < p; ++j)
				dp->v_ik[i][k] += dp->dmat[i][j] == seeds[k][j] ? dp->prob_t[i][j] : dp->prob_f[i][j];
			if (dp->v_ik[i][k] > max) {
//				max2 = max;
				max = dp->v_ik[i][k];
//				max2_id = max_id;
				max_id = k;
			}
//			else if (dp->v_ik[i][k] > max2 && dp->v_ik[i][k] < max) {
//				max2 = dp->v_ik[i][k];
//				max2_id = k;
//			}
		}
		nclass[max_id]++;
		ic[i] = max_id;
		
		/* update log probability of each cluster, position, haplotype nucleotide */
		for (j = 0; j < p; ++j)
			for (l = 0; l < dp->tot_n_categories; ++l)
				dp->e_kjN[ic[i]][j][l] += dp->dmat[i][j] == dp->categories[l] ? dp->prob_t[i][j] : dp->prob_f[i][j];
	}
	
	/* Store e_kjN; Second step: recompute the center */
	for (k = 0; k < K; ++k) {
		for (j = 0; j < p; ++j) {
			dp->last_cent[k][j] = seeds[k][j];
			double max_llj = -INFINITY;
			int max_id = -1;
			for (l = 0; l < dp->tot_n_categories; ++l)
				if (dp->e_kjN[k][j][l] > max_llj) {  /*find the largest e_kjN[k][j][l] for each site*/
					max_llj = dp->e_kjN[k][j][l];
					max_id = l;
				}
			seeds[k][j] = dp->categories[max_id];
			/* Update vik */
			if (seeds[k][j] != dp->last_cent[k][j]) {
				for (i = 0; i < n; ++i) {
					v_ik_diff = dp->prob_t[i][j] - dp->prob_f[i][j];
					if (dp->dmat[i][j] == seeds[k][j])
						dp->v_ik[i][k] += v_ik_diff;
					else if (dp->dmat[i][j] == dp->last_cent[k][j])
						dp->v_ik[i][k] -= v_ik_diff;
				}
			}
		}
		dp->live[k] = n + 1;	/* all clusters start live */
	}
	
	for (i = 0; i < n; ++i)
	/* initialize the cost */
		cost[ic[i]] += dp->v_ik[i][ic[i]];
	
} /* hw_init */

/**
 * Optimal transfer stage

 @param d data pointer
 @param seeds centers
 @param K no. clusters to compare
 @param n no. observations
 @param p no. of coordinates
 @param nclass no. of observations in each cluster
 @param ic observation assignment
 @param cost cost
 @param i ith observation
 @param Is_Live to justify if it is in the live set or not
 @param indx first 0, number of consecutive observations not transferred
 */
void optra_haplotype(void *d, data_t **seeds, int K, int n, int p, unsigned int *nclass, unsigned int *ic, double *cost, int i, int Is_Live, unsigned int *indx) {
	
	data *dp = (data *)d;
	unsigned int i1, k, j, l;
	double max_llj, v_ik_diff, diff = 0;
	int max_id;
	
	PRINT_VECTOR(dp->live, K);
	
	/* Store the last assignment */
	dp->last_assign[i] = ic[i];
	
	for (k = 0; k < K; ++k) {
		
		/* Store the last cost */
		dp->last_cost[dp->index[k]] = cost[dp->index[k]];
		
		for (i1 = 0; i1 < n; ++i1)
		/* Store the last v_ik of all the ovservations for reset */
			dp->last_vik[i1][dp->index[k]] = dp->v_ik[i1][dp->index[k]];
		
		for (j = 0; j < p; ++j) {
			/* Store the current center for the following comparision */
			dp->last_cent[dp->index[k]][j] = seeds[dp->index[k]][j];
			/* Store the last ekjn for reset */
			for (l = 0; l < dp->tot_n_categories; ++l)
				dp->last_ekjn[dp->index[k]][j][l] = dp->e_kjN[dp->index[k]][j][l];
		}
	} //CHECKED!
	
	printf("last_cent:\n");
	PRINT_MATRIX(dp->last_cent, K, p);
	/* Update the centers if move i from cluster l to the others */
//	for (j = 0; j < p; ++j) {
//		for (l = 0; l < dp->tot_n_categories; ++l) {
//			printf("DATA: %d\t, cat: %d\n", dp->dmat[i][j], dp->categories[l]);
//			if (dp->dmat[i][j] == dp->categories[l]) {
//
//				for (k = 0; k < K; ++k) {
//					if (dp->index[k] == dp->last_assign[i]) {
//						dp->e_kjN[dp->index[k]][j][l] -= dp->prob_t[i][j];
//						printf("k: %d, lass: %d, prob_t: %0.2lf\n", dp->index[k], dp->last_assign[i], dp->prob_t[i][j]);
//					}
//
//					else dp->e_kjN[dp->index[k]][j][l] += dp->prob_t[i][j];
//				}
//			}
//			else {
//				for (k = 0; k < K; ++k) {
//					if (dp->index[k] == dp->last_assign[i]) {
//						dp->e_kjN[dp->index[k]][j][l] -= dp->prob_f[i][j];
//						printf("k: %d, lass: %d, prob_f: %0.2lf\n", dp->index[k], dp->last_assign[i], dp->prob_f[i][j]);
//
//					}
//					else dp->e_kjN[dp->index[k]][j][l] += dp->prob_f[i][j];
//				}
//			}
//		}
//	}
	
	for (k = 0; k < K; ++k) {
		if (dp->index[k] != dp->last_assign[i]) {
			for (j = 0; j < p; ++j) {
//				printf("k: %d, lass: %d\n", dp->index[k], dp->last_assign[i]);
				for (l = 0; l < dp->tot_n_categories; ++l) {
//					printf("cat: %d\n", dp->categories[l]);
					if (dp->dmat[i][j] == dp->categories[l]) {
						dp->e_kjN[dp->index[k]][j][l] += dp->prob_t[i][j];
//						printf("prob_t: %0.2lf\n", dp->prob_t[i][j]);
//						printf("schange_ekjn: %0.2lf\n", dp->e_kjN[dp->index[k]][j][l]);
					}
					else {
						dp->e_kjN[dp->index[k]][j][l] += dp->prob_f[i][j];
//					printf("prob_f: %0.2lf\n", dp->prob_f[i][j]);
//						printf("schange_ekjn: %0.2lf\n", dp->e_kjN[dp->index[k]][j][l]);
						
					}
				}
			}
		}
		else {
			for (j = 0; j < p; ++j) {
				for (l = 0; l < dp->tot_n_categories; ++l) {
					if (dp->dmat[i][j] == dp->categories[l]) {
						dp->e_kjN[dp->index[k]][j][l] -= dp->prob_t[i][j];
					}
					else dp->e_kjN[dp->index[k]][j][l] -= dp->prob_f[i][j];
				}
			}
			
		}
	}
	
	printf("\nekjn:\n");
	PRINT_3DARRAY(dp->e_kjN, K, p, dp->tot_n_categories);

	for (k = 0; k < K; ++k) {
		
		for (j = 0; j < p; ++j) {
			
			max_llj = -INFINITY;
			max_id = -1;
			
			for (l = 0; l < dp->tot_n_categories; ++l)
				
				if (dp->e_kjN[dp->index[k]][j][l] > max_llj) {  /*find the largest e_kjN[k][j][l] for each site*/
					max_llj = dp->e_kjN[dp->index[k]][j][l];
					max_id = l;
				}
			seeds[dp->index[k]][j] = dp->categories[max_id]; //Compute new Hk
			
			printf("seeds: %d, %d, %d\n", k, j, seeds[dp->index[k]][j]);
			
			if (seeds[dp->index[k]][j] != dp->last_cent[dp->index[k]][j])
				for (i1 = 0; i1 < n; ++i1) {
					
					v_ik_diff = dp->prob_t[i1][j] - dp->prob_f[i1][j];
					if (dp->dmat[i1][j] == seeds[dp->index[k]][j])
						dp->v_ik[i1][dp->index[k]] += v_ik_diff;
					else if (dp->dmat[i1][j] == dp->last_cent[dp->index[k]][j])
						dp->v_ik[i1][dp->index[k]] -= v_ik_diff;
				}
		}
	}
	
	printf("new centers:\n");
	PRINT_MATRIX(seeds, K, p);
	
	printf("new vik: \n");
	PRINT_MATRIX(dp->v_ik, n, K);
	
	/* recompute the cost, first initialize it */
	for (k = 0; k < K; ++k)
		cost[dp->index[k]] = dp->v_ik[i][dp->index[k]];
	cost[dp->last_assign[i]] = 0;
	
	printf("initialize cost\n");
	PRINT_VECTOR(cost, K);
	
	if (Is_Live == 1) {
		printf("!!");
		for (i1 = 0; i1 < n; ++i1)
			if (i1 != i)
				cost[ic[i1]] += dp->v_ik[i1][ic[i1]];
	}
	
	else {
		
		printf("??");
		/* Only consider the live sets */
		for (k = 0; k < K; ++k)
			for (i1 = 0; i1 < n; ++i1)
				if (i1 != i)
					cost[dp->index[k]] += dp->v_ik[i1][dp->index[k]];
	}
	
	printf("new cost:\n");
	PRINT_VECTOR(cost, K);
	
	printf("last cost:\n");
	PRINT_VECTOR(dp->last_cost, K);
	
	/* find the cluster that increases likelihood the most */
	double max_dif = 0;
	
	for (k = 0; k < K; ++k) {
		//					printf("cost[k]%.01lf\n", cost[k]);
		//					printf("last_cost[k]%.01lf\n", dp->last_cost[k]);
//							printf("cost[ic[i]]%.01lf\n", cost[ic[i]]);
		//					printf("last_cost[ic[i]]%.01lf\n", dp->last_cost[ic[i]]);
		if (dp->index[k] != dp->last_assign[i]) {
			diff = cost[dp->index[k]] - dp->last_cost[dp->index[k]] + cost[ic[i]] - dp->last_cost[ic[i]];
			printf("k: %d\n", dp->index[k]);
			printf("diff: %0.1lf\n", diff);
		}
		
		if (diff > max_dif) {
			max_dif = diff;
			ic[i] = dp->index[k];
		}
	}
	
	printf("indx_bf_oper_in: %d\n", *indx);
	
//	printf("updated vik:\n");
//	PRINT_MATRIX(dp->v_ik, n, K);
//
	if (dp->last_assign[i] != ic[i]) {
		
		*indx = 0;
		nclass[ic[i]]++;
		nclass[dp->last_assign[i]]--;
		
		printf("lastass: %d\n", dp->last_assign[i]);
		printf("ass:%d\n", ic[i]);
	
		for (k = 0; k < K; ++k) {
			
			if (dp->index[k] != ic[i] && dp->index[k] != dp->last_assign[i]) {
				
				/* Reset vik and ekjn to the last value for the clusters that are not involved in the transformation of observation i */
				for (i1 = 0; i1 < n; ++i1)
					dp->v_ik[i1][dp->index[k]] = dp->last_vik[i1][dp->index[k]];
				for (j = 0; j < p; ++j) {
					for (l = 0; l < dp->tot_n_categories; ++l)
						dp->e_kjN[dp->index[k]][j][l] = dp->last_ekjn[dp->index[k]][j][l];
					/* reset the centers of those unaffected clusters */
					seeds[dp->index[k]][j] = dp->last_cent[dp->index[k]][j];
				}
			}
			/* update the live set */
			else dp->live[dp->index[k]] = n + i + 1;
		}
		printf("indx_af_oper_in: %d\n", *indx);
	}
	else {
		/* Reset vik, ekjn and the centers to their last values */
		for (k = 0; k < K; ++k) {
			for (i1 = 0; i1 < n; ++i1)
				dp->v_ik[i1][dp->index[k]] = dp->last_vik[i1][dp->index[k]];
			for (j = 0; j < p; ++j) {
				for (l = 0; l < dp->tot_n_categories; ++l)
					dp->e_kjN[dp->index[k]][j][l] = dp->last_ekjn[dp->index[k]][j][l];
				seeds[dp->index[k]][j] = dp->last_cent[dp->index[k]][j];
				
			}
		}
	}
//	PRINT_MATRIX(dp->v_ik, n, K);
//	printf("old vik:\n");
//	PRINT_MATRIX(dp->last_vik, n, K);
//	
//	printf("old ekjn:\n");
//	PRINT_3DARRAY(dp->e_kjN, K, p, dp->tot_n_categories);
//	PRINT_3DARRAY(dp->last_ekjn, K, p, dp->tot_n_categories);
////
	PRINT_MATRIX(seeds, K, p);
	
}

