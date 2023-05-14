# Foxes and Rabbits
This repository contains a programming project regarding Parallel and Distributed Computing

```python
MPI_Init (&argc, &argv);

MPI_Comm_rank (MPI_COMM_WORLD, &id);
MPI_Comm_size (MPI_COMM_WORLD, &p);
```

## 1. Project Description
The purpose of this class project is to give students hands-on experience in parallel programming on
both shared-memory and distributed-memory systems, using OpenMP and MPI, respectively. For this
assignment you are to write a sequential and two parallel implementations of a program to simulate
an ecosystem with two species: the Iberian fox (Vulpes vulpes) and the Iberian rabbit (Oryctolagus
cuniculus).

The simulation takes place on a square grid containing cells. At the start, some of the cells are
occupied by either a rabbit, a fox, or a rock, the rest are empty. The simulation consists of computing
how the population evolves over discrete time steps (generations) according to certain rules described
next.

## 2. Requirements
(TODO)

## 3. Simulation Rules
Initially the grid is populated with rocks, rabbits and foxes. The world grid has a finite size and
animals cannot move outside its boundaries. The i-axis (vertical) and j-axis (horizontal) both start at
0 (in the upper left corner of the grid) and end at limit-1 (with dimensions supplied on the command
line for each simulation).
The animals move in this grid, and can breed and/or die. At the start of each new generation,
the age of each animal is incremented. If its age reaches its breeding period, then it reproduces as
described next. Foxes may starve if they don’t eat a rabbit within a specified number of generations.
The period for breeding and starving is defined at the start of the simulation. All animals are born
with age 0, when the simulation is started and a new generation is going to be computed they have all
been incremented so that during generation 1 all initial animals have age 1.

# 3.1 Rules for Rocks
Rocks don’t move and neither animal can occupy cells with rocks (they are too steep to climb).

# 3.2 Rules for Rabbits
At each time step, a rabbit tries to move to a neighboring empty cell, picked according to the method
described below if there are multiple empty neighbors. If no neighboring cell is empty, it stays.
Rabbits move up or down, left or right, but not diagonally (like Rooks not Bishops).
If the age of the rabbit reaches the breeding age, when it moves it breeds, leaving behind a rabbit
of age 0, and its age is reset to 0. A rabbit cannot breed if it doesn’t move (and its age doesn’t get
reset until it actually breeds).
Rabbits never starve.

# 3.3 Rules for Foxes
At each time step, if one of the neighboring cells has a rabbit, the fox moves to that cell eating the
rabbit. If multiple neighboring cells have rabbits, then one is chosen using the method described
below. If no neighbors have rabbits and if one of the neighboring cells is empty, the fox moves
there (picked using the method described below from empty neighbors if there is more than one).
Otherwise, it stays. Same as rabbits, foxes move up or down, left or right, but not diagonally.
If a fox reaches a breeding age, when it moves it breeds, leaving behind a fox of age 0 (both for
breeding and eating), and its breeding age is reset to 0. A fox cannot breed if it doesn’t move (and its
breeding age doesn’t get reset until it actually breeds).
Foxes only eat rabbits, not other foxes. If a fox reaches a starvation age (generations since last
having eaten) and doesn’t eat in the current generation it dies. Hence, death occurs at the end of the
current generation, namely after having reproduced.

# 2.4 Traversal Order Independence
Since we want the simulation to have the same result independently of the order that the individual
cells are processed, you should implement the simulator using a so-called red-black scheme (as in a
checkerboard, depicted in Figure 1) for updating the board for each generation.

Figure 1: Red-black cell organization.
(TODO)

This means that a actually generation consists of two sub-generations. In the first sub-generation,
only the red cells are processed, and in the second sub-generation the black cells are processed. In an
even numbered row red cells are the ones with an even column number, and in an odd numbered row
red cells are the ones with an odd column number.

The red-black scheme will help parallelize the computation of the cells. Cells of the same color
are diagonal to each other, while the animals can only move parallel to the axes. This allows the
evaluation of all the positions in a sub-generation, i.e., over the red or black cells, in parallel. The
red-black scheme allows you to think of each sub-generation as a separate parallel loop over the red
(or black) cells, with no dependencies between the iterations of the loop.

Note that in the red-black scheme a rabbit or fox could end up moving twice in a generation. You
need to prevent this, making sure that if the animal has moved in first sub-generation it will not move
again in the next.

The rules that follow about selecting cells and resolving conflicts apply for each sub-generation.

# 2.5 Rules for Selecting a Cell when Multiple Choices Are Possible
If multiple cells are open for either a rabbit or a fox to move to, or for a fox to move to eat a rabbit,
number the possible choices starting from 0, clockwise starting from the 12:00 position (i.e., up, right,
down, left). Note that only cells that are unoccupied (for moves) or occupied by rabbits (for foxes to
eat) should be numbered. Let p be the number of possible moves.
Determine C, the grid cell number where the animal being evaluated is positioned. If the position
of this cell is (i, j) in a world grid with (0, 0) as the grid origin, and M × N the grid size, the grid cell
number is given by C = i × N + j.

Then, the cell to select is then determined by C mod p. For example, if there are p = 3 possible
cells to choose from, say up, down and left, then if C mod p is 0 the selected cell is up from the
current cell; if it is 1 then select down; and if it is 2 then select left.

# 2.6 Resolving Conflicts
A cell may get updated multiple times during one generation due to rabbit or fox movements. If a cell
is updated multiple times in a single generation, the conflict is resolved as follows:
1. If two rabbits end up in a cell, then the cell ends up with the rabbit with the greatest current age
(closest to breeding – bigger rabbits win). The other rabbit disappears.
2. If two foxes end up in a cell, then, just like in the case of the rabbits, the cell ends up with the
fox with the greatest current age (closest to breeding – bigger foxes win). If the foxes have
the same age, the resulting fox gets the lowest starvation age of the tied foxes. The other fox
disappears. Note: If one of the foxes has starvation age 0, the resulting fox must have starvation
0 (at some point a rabbit moved to this cell and the winning fox ate it).
3. If a fox and a rabbit end up in a cell, of course the fox eats the rabbit and gets its starvation age
reset.

Multiple conflicts can be solved applying the above rules in any order.



## 4. Input Data

The program takes ten command line parameters, all positive integers:

`$ foxes-rabbits <# generations> <M> <N> <# rocks> <# rabbits> <rabbit breeding> <# foxes> <fox breeding> <fox starvation> <seed>`

## 5. Implementations
Distint versions of the simulation have been implemented. For base reference, there is a serial version of the program. Then, we explore parallelism through OpenMP and later through MPI. Finally, OpenMP and MPI are combined to further improved the speed of the simulation.

### 5.1. Serial implementation
### 5.2. OpenMP implementation
### 5.3. MPI implementation
### 5.4. OpenMP + MPI implementation


