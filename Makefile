CC=gcc-5
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
cross:
	$(CC) -c src/grid.c $(CFLAGS); \
	$(CC) -c src/cross.c $(CFLAGS); \
	$(CC) -c src/main_cross.c $(CFLAGS); \
	$(CC) main_cross.o cross.o grid.o $(CFLAGS) -o cross_solver

gs-ido:
	$(CC) -c src/grid.c $(CFLAGS); \
    $(CC) -c src/cross.c $(CFLAGS); \
    $(CC) -c src/gs-ido.c $(CFLAGS); \
    $(CC) -c src/ido_original.c $(CFLAGS); \
    $(CC) -c src/main_test.c $(CFLAGS); \
    $(CC) main_test.o cross.o grid.o gs-ido.o ido_original.o $(CFLAGS) -o test_solver


CFD:
	$(CC) -c src/grid.c $(CFLAGS); \
	$(CC) -c src/cross.c $(CFLAGS); \
	$(CC) -c src/main_cross_CFD.c $(CFLAGS); \
	$(CC) main_cross_CFD.o cross.o grid.o $(CFLAGS) -o CFD_cross_solver

clean:
	rm *solver *.o

#hellomake: hellomake.o hellofunc.o 
#	gcc -o hellomake hellomake.o hellofunc.o -I.


