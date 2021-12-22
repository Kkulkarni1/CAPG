/**
 * @file stages.h
 * @ Author Xiyu Peng
 * 
 * Function for multi-stages Ampliclust
 */
#ifndef __H_STAGES__
#define __H_STAGES__

#include <stddef.h>

#include "fastq.h"
#include "error.h"
#include "ampliclust.h"
#include "run_ampliclust.h"

/**
 *  Set of the confidence levels
 *  Note they should be divided by 1000 when being used
 * */

enum {
    ALPHA_SIZE = 5,
    ALPHA_NU = 1000 
};

enum {
    ALPHA_1=50,
    ALPHA_2=25,
    ALPHA_3=10,
    ALPHA_4=5,
    ALPHA_5=1
};

enum{
    D_SET,
    SAMPLE
} ;

enum{
    STAGE_SAMPLE,
    FINAL_SAMPLE
};

typedef struct _stages stages;

/**
 * Data.  Store the fastq file.
 */

 struct _stages{
     /*
     I plan to use just one array D to store all index of D_s of stage s.
     At the s stage, if we know |D_s| and we just need to find the first |D_s|
     element in the array D. To achieve this, we just need to resort the whole array
     at the end of each stage s and let the index of D_{s+1} at
     the beginning of the array. Also, store the size of |D_{s+1}|.

     Also I think about keeping the sample_idx at every stage s.
     Maybe we could make the sample_idx of stage s at the end of
     dataset D_s at every stage.

     D_0 and D_1 is the same as dat->read_idx.
     */

     size_t *D_idx;       /*<! Index of the dataset for all stages */
     size_t *D_size;      /*<! sizes of D_s of every stages */
     size_t *sample_size; /*<! number of reads in the samples of every stages */

     unsigned int current_stage; /*<! current stage */

     /* for the sample S_s at every stage */
     size_t *sample_idx;         /*<! index of the reads in the current sample */

     /* for reads of different length and different trimmed length */
     /* I wander if we should set min_offset here */
     unsigned int *max_offset;	/*<! maximum offset in the sample at each stage */
     unsigned int *max_read_position;	/*<! maximum read length + offset in the sample at each stage*/
     unsigned int *max_read_length;	/*<! maximum observed read length in the sample at each stage */

     /* for PvalueGenerator() */
     double *pvalue; /*<! p-value for the set D_s */
     double *confi_level; /*< ! final confidence level chosen of every stages */
     unsigned char * simu_read; /*<! a simulated nucleotide sequence */
     unsigned char * simu_qual; /*<! a simulated quality sequence */
     double * simu_test_stat; /*<! an array to store all simulated test statistics at stage s */

     /* size of the below increase at every stage */
     unsigned int sum_K;             /*< Number of clusters identified so far */
     unsigned int *all_K;            /*<! an array to store number of clusters at every stage */
     unsigned char *all_haplotypes; /*<! all haplotypes identified so far */
     double *all_pi;                /*<! all pi identified so far */
     double *delta_s;              /*<! An array to store all pointers to delta at every stage*/
     double *gamma_s;              /*<! An array to store all pointers to gamma at every stage*/
     double *lambda0_s;              /*<! An array to store all pointers to lambda at every stage*/
     double *lambda1_s;              /*<! An array to store all pointers to lambda at every stage*/

    /* final estimates of K and pi by classification */
     unsigned int best_K ;              /*<! final K identified at the end */
     double *best_pi ;            /*<! final pi estimated at the end */
     unsigned char *best_haplotypes; /*<! final hyplotypes identified at the end */
     double *best_delta;              /*<! best delta estimated at the end*/
     double *best_gamma;              /*<! best gamma estimated at the end*/
     double *best_lambda0;              /*<! best lambda0 estimated at the end */
     double *best_lambda1;              /*<! best lambda1 estimated at the end*/

     double *eik;		/*<! expectations */
     double ll;         /*<! log likelihood */

 };/* stages */

 /* main function*/
int multi_stage(model *mod, data *dat, options *opt, initializer *ini, simulator *sim, run_info *ri, stages *stag);

/* stages object */
int make_stages(stages **stag, data *dat, options *opt); /* with fdata loaded, before data object */
void free_stages(stages *stag);
int update_stages(stages *stag,data *dat,model *mod,options *opt, unsigned int current);

/* cluster identification */
int reconstruct_sample(stages *stag,options *opt);
int reorder_D_idx(stages *stag, int type);
int pvalue_generator(data *dat, stages *stag,model *mod,options *opt);
double likelihood_ratio(model *mod, unsigned char *read,
	unsigned char *qual, unsigned int len, unsigned int off,unsigned int max_read_length);

/* reclassification */
int reestimate_parameters(data *dat,stages *stag,options *opt);
int final_assignment(data *dat,stages *stag, options *opt,run_info *ri);
int eliminate_nullcluster(stages *stag, data *dat,options *opt);
int generate_E(stages *stag, data *dat,unsigned int K,FILE *fp);
int reestimate_dispersion(stages *stag, data *dat);


#endif

