CFLAGS=-std=c11 -g -O3 -fopenmp -fstrict-aliasing -ftree-vectorize -march=native -mtune=native -fopt-info-vec-missed -fopt-info-vec-optimized -fopt-info-loop-missed -fopt-info-loop-optimized
CFLAGS_QUIET=-std=c11 -g -O3 -fopenmp -fstrict-aliasing -ftree-vectorize -march=native -mtune=native
LDFLAGS=-fopenmp
LDLIBS=-lquadmath -lm

CPPFLAGS=-Ivectorclass
CXXFLAGS=-g -O3 -fstrict-aliasing -ftree-vectorize -march=native -mtune=native -fopt-info-vec-optimized -fopt-info-vec-missed -fopt-info-loop-optimized -fopt-info-loop-missed

globalsums: globalsums.o do_sum.o do_sum_wdigittrunc.o do_sum_wbittrunc.o do_ldsum.o do_ldsum_wdigittrunc.o do_ldsum_wbittrunc.o \
            do_kahan_sum.o do_kahan_sum_v.o do_kahan_sum_gcc_v.o do_knuth_sum.o do_knuth_sum_v.o

globalsums.o: globalsums.c
	$(CC) $(CPPFLAGS) $(CFLAGS_QUIET) -c $<

clean:
	rm -f globalsums *.o
