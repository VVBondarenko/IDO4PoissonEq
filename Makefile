CC=gcc
CFLAGS=-I./inc/ -g -lm -fopenmp -fopenmp-simd -O3 -ffast-math
#-lgsl -lgslcblas

SRC = src/main.c \
	src/grid.c

HEADERS = inc/grid.h

OBJ = main.o \
	grid.o
	
%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

all:
	$(CC) -c $(SRC) $(CFLAGS); \
    $(CC) $(OBJ) $(CFLAGS) -o solver
build: $(OBJ)
	$(CC) $(OBJ) $(CFLAGS) -o solver


#hellomake: hellomake.o hellofunc.o 
#	gcc -o hellomake hellomake.o hellofunc.o -I.


