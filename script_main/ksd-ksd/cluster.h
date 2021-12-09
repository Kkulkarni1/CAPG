#ifndef _CLUSTER_H
#define _CLUSTER_H

#include <stdlib.h>
#include <stdio.h>
#include <time.h>




enum {	MUTUAL_INFORMATION,		/*!< see Romano2014 */
	NORMALIZED_MUTUAL_INFORMATION,	/*!< see Vinh2010: really scaled */
	ADJUSTED_MUTUAL_INFORMATION,	/*!< see Vinh2010 */
	STANDARDIZED_MUTUAL_INFORMATION,/*!< don't know */
	VARIATION_OF_INFORMATION,	/*!< see Meila2007 */
	NORMALIZED_VARIATION_OF_INFORMATION,	/*!< see Vinh2010 */
	NUMBER_OF_INFORMATION_MEASURES
};

enum {	NO_MI_SCALING,		/*!< no scaling of mutual information */
	MIN_SCALING,		/*!< min(h(a),h(b)) */
	GEOMETRIC_MEAN_SCALING,	/*!< sqrt(h(a)*h(b)) */
	MEAN_SCALING,		/*!< (h(a)*h(b))/2 */
	MAX_SCALING,		/*!< max(h(a),h(b)) */
	JOINT_ENTROPY_SCALING,	/*!< h(a,b) */
	NUMBER_OF_SCALING_METHODS
};

enum {	RAND_INDEX,
	ADJUSTED_RAND_INDEX,
	E_INDEX,
	NUMBER_OF_INDICES
};

double mutual_information(unsigned int *id_est, unsigned int *id_true, unsigned int n, unsigned int k_est, unsigned int k_true, int mi_variant, int normalization);
double biological_homogeneity_index(unsigned int *id_est, unsigned int *id_true, unsigned int n, unsigned int k_est);
double cluster_index(unsigned int *id_est, unsigned int *id_true, unsigned int n, unsigned int k_est, unsigned int k_true, int rand_type);
void free_cluster_statics();	/* WARNING: YOU MUST CALL THIS FUNCTION IF YOU CHANGE K, EITHER TRUE K OR ESTIMATED K */

#endif
