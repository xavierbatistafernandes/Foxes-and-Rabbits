# Foxes and Rabbits
This repository contains a programming project regarding Parallel and Distributed Computing

```python
MPI_Init (&argc, &argv);

MPI_Comm_rank (MPI_COMM_WORLD, &id);
MPI_Comm_size (MPI_COMM_WORLD, &p);
```

## 1. Requirements

## 2. Simulation Rules

## 3. Input Data

The program takes ten command line parameters, all positive integers:

`$ foxes-rabbits <# generations> <M> <N> <# rocks> <# rabbits> <rabbit breeding> <# foxes> <fox breeding> <fox starvation> <seed>`

## 4. Implementations
Distint versions of the simulation have been implemented. For base reference, there is a serial version of the program. Then, we explore parallelism through OpenMP and later through MPI. Finally, OpenMP and MPI are combined to further improved the speed of the simulation.

### 4.1. Serial implementation
### 4.2. OpenMP implementation
### 4.3. MPI implementation
### 4.4. OpenMP + MPI implementation


