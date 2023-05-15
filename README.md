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

### Open Multithread-Processing (OpenMP)
This project uses the OpenMP library.

Through OpenMP it is possible to exploit parallelism within a shared-memory system (for instance, your local computer).

OpenMP provides directives for defining parallel regions (independent work loads) and work-sharing constructs. These enable 
the creation of multiple threads that execute smaller tasks simultaneously within the same memory system.

To specify the number of threads during execution, you need to change the `OMP_NUM_THREADS` environment variable.

```
$ export OMP_NUM_THREADS=<number of threads>
```

For OpenMP documentation, see:
- https://www.openmp.org/
- http://gcc.gnu.org/projects/gomp
- https://computing.llnl.gov/tutorials/openMP/


### Open Message Passing Interface (OpenMPI)
This project also utilizes the Message Passing Interface (MPI) implementation called OpenMPI.

MPI enables parallelization not only in shared-memory systems but also in distributed-memory systems, such as clusters of computers.

To execute the simulation across multiple computers, a cluster setup is required, which is beyond the scope of this repository.

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
how the population evolves over discrete time steps (generations) according to certain rules which are as follows:

(TODO)

A more detailed description of the simulation rules this project follows can be read in the 
'project-description.pdf' file, in the Simulation Rules chapter.



## Program Execution

All four implementations are organized in this repository inside the 'src' directory.
For simplicity, each implementation consists of a single source file and its respective makefile.
To compile a specific version, move to the corresponding directory inside the 'src' directory and run the 'make' command to generate the executable.

The program takes ten command line parameters, all positive integers:

```
$ ./foxes-rabbits <# generations> <M> <N> <# rocks> <# rabbits> <rabbit breeding> <# foxes> <fox breeding> <fox starvation> <seed>
```

Remember that for parallel implementions using OpenMP, you need to change the `OMP_NUM_THREADS` environment variable.



## Implementation Design

### Serial

### OpenMP

### MPI

### Hybrid (OpenMP + MPI)


## Perfomance Analysis

![Table 1](images/table1.png)
![Table 2](images/table2.png)


