CFLAGS=-fopenmp -std=c11 -g -O3 -ftree-vectorize
#CFLAGS=-fopenmp -std=c11 -g -O3 -ftree-vectorize -fopt-info-missed 
#CFLAGS=-fopenmp -std=c11 -g -Wall
LDFLAGS=-fopenmp
LDLIBS=-lquadmath

globalsums: globalsums.o

clean:
	rm -f globalsums *.o
