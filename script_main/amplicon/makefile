CC =	gcc

CFLAGS = -Wall -Wextra -pedantic -O3 -std=c11 -D_POSIX_C_SOURCE=200112L -D_XOPEN_SOURCE=500		# c11 with POSIX 20 112L and optimization

# -e DBG=1
ifdef DBG
	CFLAGS = -Wall -Wextra -pedantic -g -std=c11 -D_POSIX_C_SOURCE=200112L -D_XOPEN_SOURCE=500 -DDEBUGGING_ON		# debugging
endif
# -e PROF=1
ifdef PROF
	CFLAGS = -Wall -Wextra -pedantic -O3 -std=c11 -D_POSIX_C_SOURCE=200112L -D_XOPEN_SOURCE=500 -pg
endif
# -e CURSES=1
ifdef CURSES
	CFLAGS += -DUSE_CURSES
endif

IFLAGS = 
LDFLAGS = -lm

ifdef CURSES
	LDFLAGS += -lncurses
endif


# -e YUDI=1
ifdef YUDI
	CFLAGS += -DYUDI -I/usr/local/Cellar/r/3.6.0/include/ -I/usr/local/Cellar/openblas/0.3.5/include/
	LDFLAGS += -L/usr/local/Cellar/r/3.6.0/lib/ -L/usr/local/opt/openblas/lib
endif

# -e CENSOR=1
ifdef CENSOR
	CFLAGS += -DCENSOR_QUALITY_SCORES
endif



# Local variables
srcs = $(wildcard *.c)
hds = $(wildcard *.h)
objs = $(srcs:.c=.o)
deps = $(srcs:.c=.d)

# remove source code that is not to be compiled into certain executables 
kmodes_objs = run_kmodes.o error.o cmdline.o io.o io_kmodes.o kmodes.o math.o cluster.o statistics.o order.o util.o timing_mach.o matrix_exponential.o sample.o r8lib.o
khaplotype_objs = run_khaplotype.o haplotype_option.o haplotype_data.o initialize_kmodes.o cluster_lloyds.o fastq.o effici.o error.o cmdline.o io.o io_kmodes.o kmodes.o math.o cluster.o statistics.o order.o util.o timing_mach.o matrix_exponential.o sample.o r8lib.o align.o
ampliclust_objs = run_ampliclust.o ampliclust.o aecm.o fastq.o error.o cmdline.o io.o kmodes.o math.o cluster.o statistics.o align.o order.o util.o options.o model_options.o initialize_options.o simulate_options.o model.o data.o simulate.o initialize.o stages.o hash.o sample.o matrix_exponential.o r8lib.o nuc.o qual.o sequence.o brent.o
fqmorph_objs = fqmorph.o fastq.o align.o error.o cmdline.o order.o util.o nuc.o sequence.o qual.o
merger_objs = merger.c fastq.c error.c align.c nuc.c sequence.c qual.c
trim_objs = trim.o trim_em.o fastq.o align.o io.o error.o cmdline.o nuc.o qual.o sequence.o
sam_objs = run_sam.o sam.o sequence.o nuc.o qual.o error.o fastq.o align.o io.o cmdline.o
roshan_objs = roshan.o sam.o sequence.o nuc.o qual.o error.o fastq.o align.o io.o cmdline.o order.o math.o brent.o pgamma.o lgamma.o gamma.o lgammacor.o dpois.o pchisq.o dnorm.o pnorm.o fmax2.o cospi.o stirlerr.o bd0.o mlutils.o chebyshev.o r8lib.o vcf.o pick_reads.o

roshan:	$(roshan_objs)
	$(CC) -o roshan $(roshan_objs) $(CFLAGS) $(LDFLAGS) -lz

sam:	$(sam_objs)
	$(CC) -o sam $(sam_objs) $(CFLAGS) $(LDFLAGS) -lz

trim:	$(trim_objs)
	$(CC) -o trim $(trim_objs) $(CFLAGS) $(LDFLAGS)

kmodes:	$(kmodes_objs)
	$(CC) -o run_kmodes $(kmodes_objs) $(CFLAGS) $(LDFLAGS) -lRmath

khaplotype:	$(khaplotype_objs)
	$(CC) -o run_khaplotype $(khaplotype_objs) $(CFLAGS) $(LDFLAGS) -lRmath

ampliclust:	$(ampliclust_objs)
	$(CC) -o run_ampliclust $(ampliclust_objs) $(CFLAGS) $(LDFLAGS) -lRmath -llapack -lopenblas

fqmorph:	$(fqmorph_objs)
	$(CC) -o fqmorph $(fqmorph_objs) $(CFLAGS) $(LDFLAGS)

merger:		$(merger_objs)
	$(CC) -o merger $(merger_objs) $(CFLAGS) $(LDFLAGS)

include $(deps)

%.d : %.c
	-@$(SHELL) -ec '$(CC) -MM $(CFLAGS) $(IFLAGS) $< \
		| sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@'

.PHONY : archive clean clear dummy

archive: 
	tar czvf ksd.tar.gz $(srcs) $(hds) makefile

clean:	
	-rm $(objs) $(deps) roshan run_kmodes run_ampliclust run_khaplotype fqmorph merger sam trim 2>/dev/null

clear:
	rm $(objs)

dummy:
	@echo $(objs)
