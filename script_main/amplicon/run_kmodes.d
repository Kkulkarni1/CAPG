run_kmodes.o run_kmodes.d : run_kmodes.c /usr/local/include/Rmath.h \
  /usr/local/include/Rconfig.h run_kmodes.h kmodes.h constants.h \
  cmdline.h error.h io.h io_kmodes.h cluster.h timing_mach.h \
  matrix_exponential.h
