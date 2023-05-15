# Foxes and Rabbits - Ecosystem Simulation
This project involves simulating an ecosystem with two species: the Iberian fox (Vulpes vulpes) and the Iberian rabbit (Oryctolagus cuniculus).
The goal is to explore the benefits of parallel programming on both shared-memory and distributed-memory systems using OpenMP and MPI, respectively.



## Program Versions
There are four versions of the program, each implementing the simulation program in a different way:

1. **Serial**: Runs the simulation sequentially without parallelization.
2. **OpenMP**: Utilizes OpenMP for parallelization on shared-memory systems.
3. **MPI**: Employs MPI for parallelization on distributed-memory systems.
4. **Hybrid**: Combines OpenMP and MPI, providing parallelization on both shared-memory and distributed-memory systems.



## Requirements

### 1. Open Multithread-Processing (OpenMP)
This project uses the OpenMP library. Through OpenMP it is possible to exploit parallelism within a shared-memory system 
(for instance, your local computer).

OpenMP provides directives for defining parallel regions (independent work loads) and work-sharing constructs. These enable 
the creation of multiple threads which can then execute simultaneously smaller tasks within the same memory system.

To specify the number of threads during execution, you need to change the `OMP_NUM_THREADS` environment variable.

```
$ export OMP_NUM_THREADS=<number of threads>
```

Replace `<number of threads>` with the desired numbers of threads to be launched during execution.

For OpenMP documentation, see:
- https://www.openmp.org/
- http://gcc.gnu.org/projects/gomp
- https://computing.llnl.gov/tutorials/openMP/


### 2. Open Message Passing Interface (OpenMPI)
This project also utilizes the Message Passing Interface (MPI) implementation called OpenMPI.

MPI enables parallelization not only in shared-memory systems but also in distributed-memory systems, such as clusters of computers.

To execute the simulation across multiple computers, a cluster setup is required, but setting it up is beyond the scope of this repository.

To run the simulation using OpenMPI on a single machine, you can use the following command:

```
$ mpirun -n <number of processes> ./foxes-rabbits-mpi <simulation arguments>
```

Where `<number of processes>` specifies the desired number of MPI processes, and `<simulation arguments>` represents the necessary arguments for the simulation.

If you have a host file specifying the machines in your cluster, you can use the following command:

```
$ mpirun -n <number of processes> -hostfile <host filepath> ./foxes-rabbits-mpi <simulation arguments>
```

Replace `<host filepath>` with the path to your host file.

For OpenMPI documentation, see:
- http://www.mpi-forum.org
- http://www.open-mpi.org
- http://www.mcs.anl.gov/research/projects/mpi/www/www3



## Simulation Rules

The simulation takes place on a square grid containing cells. At the start, some of the cells are
occupied by either a rabbit, a fox, or a rock, the rest are empty. The simulation consists of computing
how the population evolves over discrete time steps (generations) according to certain rules.

A simplified version, for context, of the rules are as follows:

- ==Rocks== don’t move and neither animal can occupy cells with rocks (they are too steep to climb).
- At each time step, a ==rabbit== tries to move to a neighboring empty cell. If no neighboring cell is empty, it stays.
- At each time step, a ==fox== moves to a cell containing a ==rabbit==, eating it. When there are no ==rabbits==, the ==foxes==
move to an empty cell.
- Both ==rabbits== and ==foxes== move up or down, left or right, but not diagonally.
- Both animals have a breeding age, which will leave a child behind when reached.
- ==Rabbits== do not suffer from starvation, while foxes need to eat a rabbit before dying from starvation.

The simulator progresses through a red-black scheme/organization:

![Figure 1](images/figure1.png)

This means that a single generation actually consists of two sub-generations. In the first sub-generation,
only the red cells are processed, and in the second sub-generation the black cells are processed.
The red-black scheme will help parallelize the computation of the cells. Cells of the same color
are diagonal to each other, while the animals can only move parallel to the axes.


For a more detailed version of the simulation rules, please consult 'project-description.pdf' inside the `docs` directory of this
repository.



## Program Execution

All four implementations are organized in this repository inside the `src` directory.
For simplicity, each implementation consists of a single source file and its respective makefile.
To compile a specific version, move to the corresponding directory inside the `src` directory and run the `make` command to generate the executable.

#### Input Data
The program takes ten command line parameters, all positive integers:

```
$ ./foxes-rabbits <# generations> <M> <N> <# rocks> <# rabbits> <rabbit breeding> <# foxes> <fox breeding> <fox starvation> <seed>
```

Remember that for parallel implementions using OpenMP, you need to change the `OMP_NUM_THREADS` environment variable.

#### Output Data
All implementations send to standard output, `stdout`, just three integers in one line separated
with one space in this order: the final number of ==rocks==, ==rabbits== and ==foxes==.


## Implementation Design

### Serial

### OpenMP

### MPI

### Hybrid (OpenMP + MPI)


## Perfomance Analysis


Table 2 contains the execution time of different tests ran on the lab machines to benchmark the performance of 
our MPI implementation of the project (checkerboard decomposition). The table shows the execution times and speedups 
obtained for 8, 16, 32 and 64 processes (inside the parenthesis are the total number of machines), for the following testing maps: (Note that no serial time is
provided for the final 3 tests as they are too big to fit on a single machine and run serially and so the shown
speedup is the speedup relative to the 8 process execution time)


![Table 1](images/table1.png)
![Table 2](images/table2.png)


