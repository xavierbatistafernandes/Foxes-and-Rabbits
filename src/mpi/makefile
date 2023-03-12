# Authors:	João Leitão		 93088
#		  	José Brito		 93106
#  		   	Xavier Fernandes 93202

CC = mpicc
CFLAGS = -Wall -Wextra -pedantic -O3 -fopenmp
EXE = foxes-rabbits-mpi

all: $(EXE)

foxes-rabbits-mpi: foxes-rabbits-mpi.c
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm $(EXE)