#ifndef __IO_KMODES_H__
#define __IO_KMODES_H__

#include <stdio.h>	/* FILE, fprintf() */
//#include <stddef.h>
#ifdef USE_CURSES
#include <curses.h>	/* WINDOW, wprintw() */
#endif

#include "constants.h"

/* print vectors */
void fprint_data_ts(FILE *fp, data_t *v, size_t len, int width, int newline);
#ifdef USE_CURSES
void wprint_data_ts(WINDOW *fp, data_t *v, size_t len, int width, int newline);
#endif

/* read vectors */
int fscan_data_ts(FILE *fp, data_t *v, size_t len);

#endif
