//
//  myfun.c
//  Test
//
//  Created by Yudi Zhang on 7/29/18.
//  Copyright © 2018 Yudi Zhang. All rights reserved.
//

#include "myfun.h"

void printdata(data *x)
{
    printf("READS:\n");
    PRINT_VECTOR(x->fdata->reads, x->n_observations * x->fdata->n_max_length);
    printf("\n");
    
    printf("QMAT:\n");
    PRINT_MATRIX(x->qmat, x->n_observations, x->fdata->n_max_length)
    ;
    
    printf("DMAT:\n");
    PRINT_MATRIX(x->dmat, x->n_observations, x->fdata->n_max_length)
    ;
}
