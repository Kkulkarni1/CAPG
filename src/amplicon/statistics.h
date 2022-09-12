#ifndef __STATISTICS_H__
#define __STATISTICS_H__

#include <math.h>
#include <stddef.h>

double aic(double ll, size_t k);
double bic(double ll, size_t k, size_t n);

#endif
