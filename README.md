# Foxes and Rabbits
This repository contains a programming project regarding Parallel and Distributed Computing

```python
MPI_Init (&argc, &argv);

MPI_Comm_rank (MPI_COMM_WORLD, &id);
MPI_Comm_size (MPI_COMM_WORLD, &p);
```


## Input Data

The program takes ten command line parameters, all positive integers:

`$ foxes-rabbits <# generations> <M> <N> <# rocks> <# rabbits> <rabbit breeding> <# foxes> <fox breeding> <fox starvation> <seed>`
