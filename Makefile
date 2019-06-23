CFLAGS=-fopenmp -std=c11 -g -O3 -ftree-vectorize -march=native -mtune=native
#CFLAGS=-fopenmp -std=c11 -g -O3 -ftree-vectorize -fopt-info-missed 
#CFLAGS=-fopenmp -std=c11 -g -Wall
LDFLAGS=-fopenmp
LDLIBS=-lquadmath -lm

globalsums: globalsums.o

clean:
	rm -f globalsums *.o
