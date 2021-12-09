#include "io_kmodes.h"
#include "math.h"
#include "error.h"

void fprint_data_ts(FILE *fp, data_t *v, size_t n, int width, int newline)
{
	size_t i;
	for (i = 0; i < n; ++i)
		if (width)
			fprintf(fp, " %*" AMPLICLUST_PRIu_data_t, width, v[i]);
		else
			fprintf(fp, " %" AMPLICLUST_PRIu_data_t, v[i]);
	if (newline) fprintf(fp, "\n");
} /* fprint_data_ts */

#ifdef USE_CURSES
void wprint_data_ts(WINDOW *wp, data_t *v, size_t n, int width, int newline)
{
	size_t i;
	for (i = 0; i < n; ++i)
		if (width)
			wprintw(wp, " %*" AMPLICLUST_PRIu_data_t, width, v[i]);
		else
			wprintw(wp, " %" AMPLICLUST_PRIu_data_t, v[i]);
	if (newline) wprintw(wp, "\n");
} /* wprint_data_ts */
#endif

int fscan_data_ts(FILE *fp, data_t *v, size_t n)
{
	size_t i;
	for (i = 0; i < n; ++i)
		if (fscanf(fp, "%" AMPLICLUST_SCNu_data_t, &v[i]) != 1)
			return FILE_FORMAT_ERROR;
	return NO_ERROR;
} /* fscan_data_ts */
