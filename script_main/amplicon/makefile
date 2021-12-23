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
sam_objs = run_sam.o sam.o sequence.o nuc.o qual.o error.o fastq.o align.o io.o cmdline.o
roshan_objs = roshan.o sam.o sequence.o nuc.o qual.o error.o fastq.o align.o io.o cmdline.o order.o math.o brent.o pgamma.o lgamma.o gamma.o lgammacor.o dpois.o pchisq.o dnorm.o pnorm.o fmax2.o cospi.o stirlerr.o bd0.o mlutils.o chebyshev.o r8lib.o vcf.o pick_reads.o

capg_amp:	$(roshan_objs)
	$(CC) -o capg_amp $(roshan_objs) $(CFLAGS) $(LDFLAGS) -lz

sam:	$(sam_objs)
	$(CC) -o sam $(sam_objs) $(CFLAGS) $(LDFLAGS) -lz


include $(deps)

%.d : %.c
	-@$(SHELL) -ec '$(CC) -MM $(CFLAGS) $(IFLAGS) $< \
		| sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@'

.PHONY : archive clean clear dummy

archive: 
	tar czvf ksd.tar.gz $(srcs) $(hds) makefile

clean:	
	-rm $(objs) $(deps) roshan sam 2>/dev/null

clear:
	rm $(objs)

dummy:
	@echo $(objs)
