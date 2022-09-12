/*
 *  @file initialize_kmodes.c
*/

#include "initialize_kmodes.h"
#include "array.h"

static inline double hd(data_t *x, data_t *y, unsigned int p);

/**
 * Hamming distance function.
 *
 * @param x    observation (1 x p)
 * @param y    observation (1 x p)
 * @param p    number of coordinates
 * @return    (weighted) Hamming distance
 */
static inline double hd(data_t *x, data_t *y, unsigned int p)
{
    double d = 0;
    for (unsigned int j = 0; j < p; ++j)
        d += (x[j] != y[j]);
    return d;
} /* hd */

int kmodes_ini (data *dat, options* opt)
{
    
    /* [KSD] The remaining code should be shared with k-modes.  Since it is
     * currently in run_kmodes.c, the best solution is to move it out of that
     * file, say to a file kmodes_initialize.c.  Then you can link with that
     * file.  Later, or instead, you may also wish to link to initializers
     * from initialize.c, which are designed for fastq data.
     */
    
    if (opt->pfile) {
        FILE *fp = fopen(opt->pfile, "r");
        for (unsigned int i = 0; i < dat->n_observations; ++i) {
            if (fscanf(fp, "%u", &dat->cluster_id[i]) != 1) {
                fclose(fp);
                return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
                                opt->pfile);
            }
            if (dat->cluster_id[i] >= opt->K) {
                fclose(fp);
                return mmessage(ERROR_MSG, INVALID_USER_INPUT,
                                "Partition file assigns to more than %u"
                                " clusters.", opt->K);
            }
        }
        fclose(fp);
    }
    
    if (opt->sfile) {
        FILE *fp = fopen(opt->sfile, "r");
        if (!fp)
            return mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt->sfile);
        
        opt->n_seed_set = 0;
        do {
            char c = fgetc(fp);
            
            /* new line with content */
            if (c != '\n' && !feof(fp))
                ++opt->n_seed_set;
            
            /* fast-forward through line */
            while (c != '\n' && !feof(fp)) c = fgetc(fp);
        } while (!feof(fp));
        
        data_t **seeds = dat->seeds; /* Important! */
        /* allocate space for seed set */
        if (opt->n_seed_set > opt->K) {
            data_t *tmp = malloc(opt->n_seed_set
                                 * dat->n_coordinates * sizeof *tmp);
            opt->seed_set = malloc(opt->n_seed_set * sizeof
                                   *opt->seed_set);
            if (!tmp || !opt->seed_set)
                return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
                                "options:seed_set");
            for (unsigned int k = 0; k < opt->n_seed_set; ++k) {
                opt->seed_set[k] = tmp;
                tmp += dat->n_coordinates;
            }
            seeds = opt->seed_set;
            opt->init_method = KMODES_INIT_RANDOM_FROM_SET;
        } else if (opt->n_seed_set < opt->K)
            mmessage(WARNING_MSG, NO_ERROR,
                     "Requesting %u clusters, but only %u seeds in "
                     "seed file '%s'.  Will generate remaining seeds"
                     " with chosen initialization method.\n", opt->K,
                     opt->n_seed_set, opt->sfile);
        else {    /* exactly options::K seeds provided */
            opt->init_method = KMODES_INIT_USER_SEEDS;
            if (opt->n_init > 1)
                mmessage(WARNING_MSG, INVALID_USER_INPUT,
                         "Resetting to one initialization.\n");
            opt->n_init = 1;
            opt->n_inner_init = 1;
        }
        
        fprintf(stderr, "Found %u seeds\n", opt->n_seed_set);
        rewind(fp);
        for (unsigned int k = 0; k < opt->n_seed_set; ++k) {
            if (fscan_data_ts(fp, seeds[k], dat->n_coordinates)) {
                fclose(fp);
                return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
                                opt->sfile);
            }
            if (opt->subtract_one)
                for (unsigned int j = 0; j < dat->n_coordinates;
                     ++j) {
                    if (seeds[k][j] == 0)
                        return mmessage(ERROR_MSG,
                                        INVALID_USER_INPUT,
                                        "Command option -1 "
                                        "requested, but -i "
                                        "<istr> has 0-based "
                                        "data.\n");
                    else
                        seeds[k][j] -= 1;
                }
            if (dat->ini_seeds && opt->n_seed_set  < opt->K)
                memcpy(dat->ini_seeds[k], seeds[k],
                       dat->n_coordinates * sizeof **seeds);
        }
        fclose(fp);
    }
    
    return NO_ERROR;
}

/**
 Initialize seeds

 @param dat data structure
 @param opt option structure
 @return err
 */
int initialize(data *dat, options *opt)
{
    int err = NO_ERROR;
    double inner_total = INFINITY, dtmp;
    data_t **seeds = dat->seeds;
    unsigned int *seed_idx = dat->seed_idx;
    
    dat->use_ini = 0;
    for (unsigned int j = 0; j < opt->n_inner_init; ++j) {
        if (!opt->n_sd_idx && !opt->pfile) {
            
            if (opt->init_method
                == KMODES_INIT_USER_SEEDS) {
                unsigned int idx[opt->K];
                
                idx[0] = 1; idx[1] = 47; idx[2] = 95;
                
                for (int i = 0; i < opt->K; ++i)
                    COPY_1ARRAY(dat->seeds[i], dat->dmat[idx[i]], dat->n_coordinates);
            }  // For the test
            
            else if (opt->init_method
                == KMODES_INIT_RANDOM_FROM_PARTITION)
                kmodes_init_random_from_partition(dat->dmat,
                                                  dat->n_observations, dat->n_coordinates,
                                                  opt->K, seeds, seed_idx,
                                                  opt->true_cluster);
            else if (opt->init_method
                     == KMODES_INIT_RANDOM_FROM_SET)
                kmodes_init_random_from_set(opt->K,
                                            dat->n_coordinates, opt->n_seed_set,
                                            seeds, opt->seed_set);
            else if (opt->init_method
                     == KMODES_INIT_RANDOM_SEEDS) {
                unsigned int k = 0, l = 0, m;
                do {
                    /*
                     memcpy(seeds[k], dat->dmat[k],
                     dat->n_coordinates * **seeds);
                     */
                    seed_idx[k] = l;
                    for (m = 0; m < dat->n_coordinates; ++m)
                        seeds[k][m] = dat->dmat[l][m];
                    
                    for (m = 0; m < k; ++m)
                        if (!hd(seeds[k], seeds[m], dat->n_coordinates))
                            break;
                    if (!k || m == k) ++k;
                    ++l;
                } while (k < opt->K);
            } else
                kmodes_init(dat->dmat, dat->n_observations,
                            dat->n_coordinates, opt->K,
                            opt->n_seed_set,
                            seeds, seed_idx,
                            opt->init_method, opt->weight);
        } else if (opt->n_sd_idx) {    /* will only happen once */
            for (unsigned int k = 0; k < opt->K; ++k) {
                for (unsigned int j = 0; j < dat->n_coordinates;
                     ++j)
                    seeds[k][j] =
                    dat->dmat[seed_idx[k]][j];
                /* retaining: WHY DOES THIS NOT WORK???
                 memcpy(dat->seeds[k],
                 dat->dmat[seed_idx[k]],
                 dat->n_coordinates * *dat->seeds[k]);
                 */
            }
        } else if (opt->pfile) {    /* will only happen once */
            kmodes_init_from_partition(dat->dmat,
                                       dat->n_observations, dat->n_coordinates, opt->K,
                                       opt->weight, seeds, dat->cluster_id);
        }
//        if (opt->n_inner_init > 1) {
//            if (opt->kmodes_algorithm == KMODES_HUANG)
//                dtmp = kmodes_huang(dat->dmat,
//                                    dat->n_observations, dat->n_coordinates,
//                                    seeds, opt->K, dat->cluster_id,
//                                    dat->cluster_size, 0, dat->criterion,
//                                    &err, &dat->iter, opt->kopt);//opt->weight, opt->update_modes);
//            else if (opt->kmodes_algorithm == KMODES_HARTIGAN_WONG)
//                dtmp = kmodes_hw(dat->dmat, dat->n_observations,
//                                 dat->n_coordinates, seeds, opt->K,
//                                 dat->cluster_id, dat->cluster_size, 0,
//                                 dat->criterion, &err, &dat->iter,
//                                 opt->kopt); //opt->weight, opt->update_modes, opt->use_qtran);
//            else
//                return mmessage(ERROR_MSG, INVALID_USER_INPUT,
//                                "Invalid algorithm in ini-kmodes.\n");
//            
//            if (opt->quiet > MINIMAL)
//                fprintf(stdout, "Inner initialization %*u of "
//                        "%u: %f (%f)\n", (int)
//                        (log10(opt->n_inner_init) + 1), j,
//                        opt->n_inner_init, dtmp, inner_total);
//            
//            /* repeat if obtain null cluster */
//            if (err == KMODES_NULL_CLUSTER_ERROR) {
//                if (j) --j;
//                continue;
//            } else if (err)
//                return err;
//            
//            if (dtmp < inner_total) {
//                inner_total = dtmp;
//                seeds = seeds == dat->seeds
//                ? dat->ini_seeds : dat->seeds;
//                seed_idx = seed_idx == dat->seed_idx
//                ? dat->ini_seed_idx : dat->seed_idx;
//                dat->use_ini = !dat->use_ini;
//            }
//        }
    }
    
    if (opt->n_inner_init > 1)
        dat->use_ini = !dat->use_ini;
    
    return err;
} /* initialize */


