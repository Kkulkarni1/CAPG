#include "statistics.h"

double aic(double ll, size_t k) {
	return (2*k - 2*ll);
} /* aic */

double bic(double ll, size_t k, size_t n) {
	return (log(n)*k - 2*ll);
} /* bic */
