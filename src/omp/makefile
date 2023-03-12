# Authors:	João Leitão		 93088
#		  	José Brito		 93106
#  		   	Xavier Fernandes 93202

CC = gcc
CFLAGS = -Wall -Wextra -pedantic -O3 -fopenmp
EXE = foxes-rabbits-omp 

all: $(EXE)

foxes-rabbits-omp: foxes-rabbits-omp.c
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm $(EXE)