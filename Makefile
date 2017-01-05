CC=gcc
CFLAGS=-I./inc/ -g -lm 
#-lgsl -lgslcblas

SRC = src/main.c \
	src/grid.c

HEADERS = inc/grid.h

OBJ = main.o \
	grid.o
	
%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

all:
	$(CC) -c $(SRC) $(CFLAGS)

build: $(OBJ)
	$(CC) $(OBJ) $(CFLAGS) -o solver


#hellomake: hellomake.o hellofunc.o 
#	gcc -o hellomake hellomake.o hellofunc.o -I.


