/*****************************************
* Authors:	João Leitão		 93088       *
*		  	José Brito		 93106       *
*  		   	Xavier Fernandes 93202       *
*****************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>
#include <string.h>
#include <mpi.h>

#define COMM_UP 0
#define COMM_DOWN 1

#define ROCK    '*'
#define RABBIT  'R'
#define FOX     'F'
#define EMPTY   ' '

#define dir_UP          0
#define dir_RIGHT       1
#define dir_DOWN        2
#define dir_LEFT        3
#define dir_UP_LEFT     4
#define dir_UP_RIGHT    5
#define dir_DOWN_LEFT   6
#define dir_DOWN_RIGHT  7  

#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n) - 1)
#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH(id,p,n) - BLOCK_LOW(id,p,n) + 1)
#define BLOCK_OWNER(index,p,n) (((p)*((index)+1)-1)/(n))

#define MSG_TAG 12345
#define NUM_GHOST_ROWS 3


typedef struct Cell {   
    
    char animal;
    int breed_age;
    int starve_age;

} Cell;


/* FUNCTION: r4_uni
 * ARGUMENTS: uint32_t *seed - Seed for randomization
 * DESCRIPTION: Generate random uniform number to be used by "generate_element"
 * RETURN: Uniform number
 */
float r4_uni(uint32_t *seed) {
    int seed_input, sseed;
    seed_input = *seed;
    *seed ^= (*seed << 13);
    *seed ^= (*seed >> 17);
    *seed ^= (*seed << 5);
    sseed = *seed;
    return 0.5 + 0.2328306e-09 * (seed_input + sseed);
}

/* FUNCTION: generate_element
 * ARGUMENTS: Cell *curr_world     - Current world state
              int n                - Number of rocks, rabbits or foxes
              char atype           - Character representative of rocks, rabbits or foxes
              uint32_t *seed       - Seed for randomization
              int M                - Number of rows
              int N                - Number of columns
              int first_row_proc   - First row to be analized by the processor
              int last_row_proc    - Last row to be analized by the processor
              int first_col_proc   - First column to be analized by the processor
              int last_col_proc    - Last column to be analized by the processor
              int local_cols       - Number of columns (local)
 * DESCRIPTION: Randomly distributes elements (rocks, foxes, rabbits) throughout
                the world matrix
 * RETURN: N/A
 */
void generate_element(Cell *curr_world, int n, char atype, uint32_t *seed, int M, int N, int first_row_proc, int last_row_proc, int first_col_proc, int last_col_proc, int local_cols){
    
    int i, j, k;
    int c = 0;
    int local_i,local_j;

    for(k = 0; k < n; k++){      
        i = M * r4_uni(seed);
        j = N * r4_uni(seed);
        if(i >= first_row_proc && i <= last_row_proc && j <= last_col_proc && j >= first_col_proc)
        {   
            local_i = i + NUM_GHOST_ROWS - first_row_proc;
            local_j = j + NUM_GHOST_ROWS - first_col_proc;

            c = local_i*local_cols + local_j;

            if(curr_world[c].animal == EMPTY)
            {
                curr_world[c].animal = atype;
            
                if(atype != ROCK)
                {
                    curr_world[c].breed_age = 0;
                    curr_world[c].starve_age = 0;
                }
            }
        }
            
    }
}

/* FUNCTION: initialize_worlds
 * ARGUMENTS: Cell *curr_world  - Current world state
              Cell *next_world  - Next world state
              int M             - Number of rows
              int N             - Number of columns
 * DESCRIPTION: Initialize worlds
 * RETURN: N/A
 */
void initialize_worlds(Cell *curr_world, Cell *next_world, int M, int N){

    int c = 0;
    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){                /*Initializing all entries as empty cells*/
            curr_world[c].animal = EMPTY;
            next_world[c].animal = EMPTY;
            c++;
        }
    }
}

/* FUNCTION: calculate_direction
 * ARGUMENTS: int x      - Number of possible moves
              int *moves - Vector marking potential movement directions
              int c      - Cell number
 * DESCRIPTION: Calculates movement direction
 * RETURN: direction
 */
int calculate_direction(int x, int *moves, int c){

    int count = 0;                  /*Count variable*/
    int direction = -1;             /*Direction variable*/

    for(int l = 0; l < 4; l++){     /*Loop that iterates over a vector to find available moves*/
        if(moves[l] == 1){          /*If in a position, the vector has a value of 1*/
            if(count == c%x){       /*Check if count is equal to c mod x*/
                direction = l;      /*If so, then direction is equal to the vector index, which corresponds to a real direction*/
                break;              /*Break loop*/
            }
            count++;                /*Increment count variable*/
        }
    }
    return direction;               /*Return*/
}

/* FUNCTION: get_moves_rabbit
 * ARGUMENTS: Cell *curr_world  - Current world state
              int i             - Current cell's row number (local)
              int j             - Current cell's column number
              int local_rows    - Number of rows (local)
              int local_cols    - Number of columns (local)
              int *p            - Number of adjacent empty cells
              int *empty_cells  - Vector marking available adjacent empty cells
              int global_row    - Current cell's row number (global)
              int global_col    - Current cell's column number (global)
              int total_rows    - Number of rows of the entire world
              int total_cols    - Number of columns of the entire world
 * DESCRIPTION: Counts the number os adjacent empty cells and registers them in a vector
                to find the number of possible movements for a rabbit cell
 * RETURN: N/A
 */                   /*Calculate p*/
void get_moves_rabbit(Cell *curr_world, int i, int j, int local_rows, int local_cols, int *p, int *empty_cells, int global_row, int global_col, int total_rows, int total_cols){

    int c = i*local_cols + j;
    curr_world[c].breed_age++;                      /*Increment rabbit breeding age*/
    if((global_row - 1 >= 0) && (i - 1 >= 0)){      /*Check if position above is in bounds*/
        if(curr_world[c-local_cols].animal == EMPTY){        /*If position above is empty*/
            empty_cells[dir_UP] = 1;                /*Mark that in the vector*/
            *p+=1;                                  /*Increment number of empty cells available*/
        }   
    }
    if((global_col + 1 < total_cols) && (j + 1 < local_cols)){                                /*Check if position to the left is in bounds*/
        if(curr_world[c+1].animal == EMPTY){        /*If position to the right is empty*/
            empty_cells[dir_RIGHT] = 1;             /*Mark that in the vector*/  
            *p+=1;                                  /*Increment number of empty cells available*/
        }
    }
    if((global_row + 1 < total_rows) && (i + 1 < local_rows)){   /*Check if position below is in bounds*/
        if(curr_world[c+local_cols].animal == EMPTY){            /*If position below is empty*/
            empty_cells[dir_DOWN] = 1;                  /*Mark that in the vector*/ 
            *p+=1;                                      /*Increment number of empty cells available*/
        }   
    }
    if((global_col - 1 >= 0) && (j - 1 >= 0)){                               /*Check if position to the left is in bounds*/
        if(curr_world[c-1].animal == EMPTY){        /*If position to the left is empty*/
            empty_cells[dir_LEFT] = 1;              /*Mark that in the vector*/ 
            *p+=1;                                  /*Increment number of empty cells available*/
        }
    }
}


/* FUNCTION: get_moves_fox
 * ARGUMENTS: Cell *curr_world      - Current world state
              int i                 - Current cell's row number (local)
              int j                 - Current cell's column number
              int local_rows        - Number of rows for the local process
              int local_cols        - Number of columns for the local process
              int *p                - Number of adjacent empty cells
              int *adj_rabbits      - Number of adjacent rabbit cells
              int *empty_cells      - Vector marking available adjacent empty cells
              int *adj_rabbit_cells - Vector marking available adjacent rabbit cells
              int global_row        - Current cell's row number (global)
              int global_col        - Current cell's column number (global)
              int total_rows        - Number of rows of the entire world
              int total_cols        - Number of columns of the entire world
 * DESCRIPTION: Searches adjacent cells, prioritising rabbits and then empty cells, to find the number of possible movements for a fox cell
 * RETURN: N/A
 */
void get_moves_fox(Cell *curr_world, int i, int j, int local_rows, int local_cols, int *p, int *adj_rabbit, int *empty_cells, int *adj_rabbit_cells, int global_row, int global_col, int total_rows, int total_cols){

    int c = i*local_cols + j;
    curr_world[c].breed_age++;                      /*Increment fox breeding age*/
    curr_world[c].starve_age++;                     /*Increment fox starving age*/

    if((global_row - 1 >= 0) && (i - 1 >= 0)){                  /*Check if position above is in bounds*/
        if(curr_world[c-local_cols].animal == RABBIT){          /*If position above is a rabbit, mark it*/
            adj_rabbit_cells[dir_UP] = 1;           
            *adj_rabbit+=1;                         /*Increment number of adjacent rabbits*/
        }
        else if(curr_world[c-local_cols].animal == EMPTY)    /*Else check if it is empty*/
        {
            empty_cells[dir_UP] = 1;                /*Mark that in the vector*/ 
            *p+=1;                                  /*Increment number of empty cells available*/
        }
    }
    if((global_col + 1 < total_cols) && (j + 1 < local_cols)){                                  /*Check if position to the left is in bounds*/
        if(curr_world[c+1].animal == RABBIT){       /*If position to the right is a rabbit, mark it*/
            adj_rabbit_cells[dir_RIGHT] = 1;
            *adj_rabbit+=1;                         /*Increment number of adjacent rabbits*/
        }           
        else if(curr_world[c+1].animal == EMPTY)    /*Else check if it is empty*/
        {
            empty_cells[dir_RIGHT] = 1;             /*Mark that in the vector*/ 
            *p+=1;                                  /*Increment number of empty cells available*/
        }   
    }
    if((global_row + 1 < total_rows) && (i + 1 < local_rows)){      /*Check if position below is in bounds*/
        if(curr_world[c+local_cols].animal == RABBIT){              /*If position below is a rabbit, mark it*/
            adj_rabbit_cells[dir_DOWN] = 1;
            *adj_rabbit+=1;                         /*Increment number of adjacent rabbits*/
        }
        else if(curr_world[c+local_cols].animal == EMPTY)           /*Else check if it is empty*/
        {
            empty_cells[dir_DOWN] = 1;              /*Mark that in the vector*/ 
            *p+=1;                                  /*Increment number of empty cells available*/
        }
    }
    if((global_col - 1 >= 0) && (j - 1 >= 0)){                       /*Check if position to the left is in bounds*/
        if(curr_world[c-1].animal == RABBIT){       /*If position to the left is a rabbit, mark it*/
            adj_rabbit_cells[dir_LEFT] = 1;
            *adj_rabbit+=1;                         /*Increment number of adjacent rabbits*/
        }
        else if(curr_world[c-1].animal == EMPTY)    /*Else check if it is empty*/
        {
            empty_cells[dir_LEFT] = 1;              /*Mark that in the vector*/ 
            *p+=1;                                  /*Increment number of empty cells available*/
        }
    }
}

/* FUNCTION: solve_conflicts_and_breeding
 * ARGUMENTS: Cell *curr_world  - Current world state
              Cell *next_world  - Next world state
              int *updated      - Vector marking already updated cells
              int c1            - Cell's current position (local)
              int i2            - Cell's next row
              int j2            - Cell's next column
              int breed_age_f   - Breeding age for foxes
              int breed_age_r   - Breeding age for rabbits
              int local_cols    - Number of columns for the local process
 * DESCRIPTION: Solves conflicts between animal in a particular cell, while also determining whether breeding occurs
 * RETURN: N/A
 */
void solve_conflicts_and_breeding(Cell *curr_world, Cell *next_world, int *updated, int c1, int i2, int j2, int breed_age_f, int breed_age_r, int local_cols){
    
    int c2 = i2*local_cols + j2;
    int conflict_type = 0;                          /*Conflict type identifier variable*/
    int breeding_flag = 0;                          /*Flag informing if breeding occurs*/

    Cell tmp;                                       /*Auxiliar Cell variable*/

    next_world[c1].animal = EMPTY;                  /*Setting the previous entry to empty because animal moved*/

    tmp.animal = curr_world[c1].animal;             /*Saving cell that is moving in auxiliar variable*/
    tmp.breed_age = curr_world[c1].breed_age;
    tmp.starve_age = curr_world[c1].starve_age;

    if(curr_world[c1].animal == FOX && curr_world[c1].breed_age >= breed_age_f) {           /*Checks if a fox reached or surpassed breeding age*/
        tmp.breed_age = 0;                                                                  /*Resets breeding age*/
        breeding_flag = 1;                                                                  /*Marks flag to know breeding must happen*/
    }

    if(curr_world[c1].animal == RABBIT && curr_world[c1].breed_age >= breed_age_r) {        /*Checks if a rabbit reached or surpassed breeding age*/
        tmp.breed_age = 0;                                                                  /*Resets breeding age*/
        breeding_flag = 1;                                                                  /*Marks flag to know breeding must happen*/
    }
    
    if(breeding_flag) {                          /*Checks if breeding must happen*/
        next_world[c1].animal = tmp.animal;      /*Leaves the bred animal behind*/
        next_world[c1].breed_age = 0;            /*Bred animal begins with 0 breed age*/
        next_world[c1].starve_age = -1;          /*Rabbits do not have hunger*/

        if(curr_world[c1].animal == FOX)         
            next_world[c1].starve_age = 0;       /*If it's a fox than it begins with 0 starvation*/
    }
    
    /*Conflict type identification*/
    if(tmp.animal == FOX && next_world[c2].animal == FOX)                 /*Conflict between two foxes*/
        conflict_type = 1;
    else if (tmp.animal == RABBIT && next_world[c2].animal == RABBIT)     /*Conflict between two rabbits*/
        conflict_type = 2;
    else if (tmp.animal == RABBIT && next_world[c2].animal == FOX)        /*Conflict between a fox and a rabbit, rabbits bumps into fox*/
        conflict_type = 3;
    else if (tmp.animal == FOX && next_world[c2].animal == RABBIT)        /*Conflict between a fox and a rabbit, fox bumps into rabbit*/
        conflict_type = 4;
        
    
    /*Conflict solving*/
    switch (conflict_type) {
        case 0:                                                     /*No conflict occured, animal just moved*/              
            next_world[c2].animal = tmp.animal;                     /*Updates the new cell with the moving animal*/
            next_world[c2].breed_age = tmp.breed_age;
            next_world[c2].starve_age = tmp.starve_age;
            break;

        case 1:                                                     /*Conflict between two foxes*/
            if(tmp.breed_age > next_world[c2].breed_age) {          /*The fox with larger breeding age wins*/
                next_world[c2].breed_age = tmp.breed_age;       
                
                if(next_world[c2].starve_age != 0) {                /*Checks if one of the foxes has eatten a rabbit*/
                    next_world[c2].starve_age = tmp.starve_age;
                }
            }
            else if(tmp.breed_age == next_world[c2].breed_age)      /*If the foxes have the same breeding age*/
            {
                if(tmp.starve_age < next_world[c2].starve_age)      /*The one with least starvation wins*/
                {
                    next_world[c2].starve_age = tmp.starve_age;
                }
            }
            break;

        case 2:                                                     /*Conflict between two rabbits*/
            if(tmp.breed_age > next_world[c2].breed_age) {          /*The one with higher beeding age wins*/
                next_world[c2].breed_age = tmp.breed_age;
            }
            break;

        case 3:                                                     /*Conflict between a fox and a rabbit, rabbits bumps into fox*/
            next_world[c2].animal = FOX;                            /*The fox always wins*/
            next_world[c2].starve_age = 0;                          /*Resetting the fox's starvation*/
            break;
            
        case 4:                                                     /*Conflict between a fox and a rabbit, fox bumps into rabbit*/
            next_world[c2].animal = FOX;                            /*The fox always wins*/
            next_world[c2].breed_age = tmp.breed_age;               /*Saving the fox's breeding age*/
            next_world[c2].starve_age = 0;                          /*Resetting the fox's starvation*/
            break;
            
        default:
            break;

    }

    updated[c2] = 1;
}

/* FUNCTION: update_next_world
 * ARGUMENTS: Cell *curr_world      - Current world state
              Cell *next_world      - Next world state
              int * updated         - Array marking already updated cells
              int i                 - Current cell's row number
              int j                 - Current cell's column number
              int direction         - Direction the animal will move
              int breed_age_f       - Fox breeding age
              int breed_age_r       - Rabbit breeding age
              int local_cols        - Number of columns for the local process
 * DESCRIPTION: Moves the animal to a new cell, depending on the direction computed, solving any
                conflicts that may occur and breeding computation
 * RETURN: N/A
 */
void update_next_world(Cell *curr_world, Cell *next_world, int *updated, int i, int j, int direction, int breed_age_f, int breed_age_r, int local_cols){
    
    int c = i*local_cols + j;
    if(direction == dir_UP)                                                                                        /*If movement direction is upwards*/
    {
        solve_conflicts_and_breeding(curr_world, next_world, updated, c, i-1, j, breed_age_f, breed_age_r, local_cols);      /*Call solving breeding and conflicts funtion for destination [i-1][j]*/
    }
    else if(direction == dir_DOWN)                                                                                 /*If movement direction is downwards*/
    {
        solve_conflicts_and_breeding(curr_world, next_world, updated, c, i+1, j, breed_age_f, breed_age_r, local_cols);      /*Call solving breeding and conflicts funtion for destination [i+1][j]*/
    }
    else if(direction == dir_LEFT)                                                                                 /*If movement direction is to the left*/
    {
        solve_conflicts_and_breeding(curr_world, next_world, updated, c, i, j-1, breed_age_f, breed_age_r, local_cols);      /*Call solving breeding and conflicts funtion for destination [i][j-1]*/
    }   
    else if(direction == dir_RIGHT)                                                                                /*If movement direction is to the right*/
    { 
        solve_conflicts_and_breeding(curr_world, next_world, updated, c, i, j+1, breed_age_f, breed_age_r, local_cols);      /*Call solving breeding and conflicts funtion for destination [i][j+1]*/
    }  
    else                                                                                                           /*If no movement occurs*/
    {   
         
        next_world[c].animal = curr_world[c].animal;                                                               /*Maintain animal in same position, but write it to "new_world"*/
        next_world[c].breed_age = curr_world[c].breed_age;
        next_world[c].starve_age = curr_world[c].starve_age;
    }
}

/* FUNCTION: movement
 * ARGUMENTS: Cell *curr_world     - Current world state
              Cell *next_world     - Next world state
              int *updated         - Vector marking already updated cells
              int local_rows       - Number of rows for the local process
              int local_cols       - Number of columns for the local process
              int i                - Cell's current row
              int j                - Cell's current column
              int breed_age_f      - Breeding age for foxes
              int breed_age_r      - Breeding age for rabbits
              int global_row       - Global index of the row of the cell being processed
              int global_col       - Global index of the column of the cell being processed
              int total_rows       - Total number of rows of the map
              int total_cols       - Total number of columns of the map
 * DESCRIPTION: Calculates the overall movement of a cell, including breeding and conflict solving, writing the result to the "next_world" vector
 * RETURN: N/A
 */
void movement(Cell *curr_world, Cell *next_world, int *updated, int local_rows, int local_cols, \
              int i, int j, int breed_age_f, int breed_age_r,                                    \
              int global_row, int global_col, int total_rows, int total_cols) {
    
    int c = i*local_cols + j;                             /*Local cell number*/
    int c2 = global_row*total_cols + global_col;          /*Global cell number*/
    int p = 0;                          /*Number of adjacent empty cells*/
    int adj_rabbit = 0;                 /*Number of adjacent rabbit cells*/
    int direction = -1;                 /*Movement direction*/

    int *empty_cells;
    int *adj_rabbit_cells;

    empty_cells = (int*) malloc(4*sizeof(int));             /*Allocating a vector to mark adjacent and non-diagonal empty cells*/
    if(empty_cells == NULL) {
        exit(-1);
    }

    adj_rabbit_cells = (int*) malloc(4*sizeof(int));        /*Allocating a vector to mark adjacent and non-diagonal rabbit cells*/
    if(adj_rabbit_cells == NULL) {
        exit(-1);
    }

    empty_cells[dir_UP] = 0;                                /*Initialize empty cells direction marking vector*/
    empty_cells[dir_RIGHT] = 0;
    empty_cells[dir_DOWN] = 0;
    empty_cells[dir_LEFT] = 0;

    adj_rabbit_cells[dir_UP] = 0;                           /*Initialize adjacent rabbit cells direction marking vector*/
    adj_rabbit_cells[dir_RIGHT] = 0;
    adj_rabbit_cells[dir_DOWN] = 0;
    adj_rabbit_cells[dir_LEFT] = 0;

    if(curr_world[c].animal == RABBIT)  /*If the current cell is a rabbit*/
        get_moves_rabbit(curr_world, i, j, local_rows, local_cols, &p, &empty_cells[0], global_row, global_col, total_rows, total_cols);                                                 /*Calculate possible moves for rabbit (how many empty cells and in which directions)*/

    if(curr_world[c].animal == FOX)     /*If the current cell is a fox*/
        get_moves_fox(curr_world, i, j, local_rows, local_cols, &p, &adj_rabbit, &empty_cells[0], &adj_rabbit_cells[0], global_row, global_col, total_rows, total_cols);                   /*Calculate possible moves for fox (how many adjacent rabbits or empty cells and in which directions)*/
    
    
    if(adj_rabbit > 0)                  /*If the number of adjacent rabbits is greater than 0 (only applies for foxes, as "adj_rabbits" is unchanged for rabbits)*/                                  
    {
        direction = calculate_direction(adj_rabbit, adj_rabbit_cells, c2);                                       /*Calculate movement direction*/
        update_next_world(curr_world, next_world, updated, i, j, direction, breed_age_f, breed_age_r, local_cols);        /*Update the next world matrix for conflict resolution in next step*/ 
    }
    else if(p >= 0)                     /*If the number of empty cells is greater than 0*/
    {
        direction = calculate_direction(p, empty_cells, c2);                                                     /*Calculate movement direction*/
        update_next_world(curr_world, next_world, updated, i, j, direction, breed_age_f, breed_age_r, local_cols);        /*Update the next world matrix for conflict resolution in next step*/   
    }

    /*Free allocated memory*/
    free(empty_cells);
    free(adj_rabbit_cells);      
}

/* FUNCTION: communicate_ghost_points
 * ARGUMENTS: Cell *curr_world       - Current world state
              Cell *next_world       - Next world state
              Cell **ghost_snd       - Vectors of ghost points to be sent
              Cell **ghost_rcv       - Vectors of ghost points to be received
              int num_ghost_pts_ver  - Number of ghost points vertically
              int num_ghost_pts_hor  - Number of ghost points horrizontally
              int num_ghost_pts_dgn  - Number of ghost points diagonally
              int num_rows_proc      - Number of rows to be processed (local)
              int num_cols_proc      - Number of columns to be processed (local)
              MPI_Comm CART_COMM     - Represents a logical group of MPI processes
              MPI_Datatype MPI_Cell  - Custom MPI datatype
              int *dims              - Optimal division dimensions
              int *source_coords     - Coordinates of the processor (local)
 * DESCRIPTION: Handles communication between processes by sending and receiving necessary data for the computation of each processes's map
 * RETURN: N/A
 */
void communicate_ghost_points(Cell *curr_world, Cell *next_world, Cell **ghost_snd, Cell **ghost_rcv,  \
                                int num_ghost_pts_ver, int num_ghost_pts_hor, int num_ghost_pts_dgn,   \
                                int num_rows_proc, int num_cols_proc,                                  \
                                MPI_Comm CART_COMM, MPI_Datatype MPI_Cell,                             \
                                int *dims, int *source_coords) {

    int offset1, offset2;                       /*Offset variables*/
    int target_rank, target_coords[2];          /*Process rank and coordinates*/
    int i;                                      /*Iteration variable*/

    int comm_flags[8] = {0};                    /*Communication direction flags*/
    MPI_Request requests_snd[8];
    MPI_Request requests_rcv[8];

    /*Initialize MPI Request variables*/
    for(i = 0; i < 8; i++) {
        requests_snd[i] = MPI_REQUEST_NULL;
        requests_rcv[i] = MPI_REQUEST_NULL;
    }
    
    /* Cheking in which directions communication should be done */
    /* Analysing vertical communication */
    if (source_coords[0] > 0) {
        comm_flags[dir_UP] = 1;
    }

    if (source_coords[0] < dims[0]-1) {
        comm_flags[dir_DOWN] = 1;
    }

    /* Analysing horizontal communication */
    if (source_coords[1] > 0) {
        comm_flags[dir_LEFT] = 1;
    }

    if (source_coords[1] < dims[1]-1) {
        comm_flags[dir_RIGHT] = 1;
    }

    /* Analysing diagonal communication */
    if (comm_flags[dir_UP] && comm_flags[dir_LEFT]) {
        comm_flags[dir_UP_LEFT] = 1;
    }

    if (comm_flags[dir_UP] && comm_flags[dir_RIGHT]) {
        comm_flags[dir_UP_RIGHT] = 1;
    }

    if (comm_flags[dir_DOWN] && comm_flags[dir_LEFT]) {
        comm_flags[dir_DOWN_LEFT] = 1;
    }

    if (comm_flags[dir_DOWN] && comm_flags[dir_RIGHT]) {
        comm_flags[dir_DOWN_RIGHT] = 1;
    }


    /* Now we exchange data between marked processes */
    /* Sending data vertically */
    offset1 = (num_cols_proc + 2*NUM_GHOST_ROWS);               /* Offset is an entire row of curr_world */

    /* Up */
    if (comm_flags[dir_UP]) {
        offset2 = NUM_GHOST_ROWS*offset1 + NUM_GHOST_ROWS;
        for (i = 0; i < NUM_GHOST_ROWS; i++) {
            memcpy(ghost_snd[dir_UP] + i*num_cols_proc, curr_world + offset2 + i*offset1, num_cols_proc*sizeof(Cell));
        }

        target_coords[0] = source_coords[0] - 1;   
        target_coords[1] = source_coords[1];       
        MPI_Cart_rank(CART_COMM, target_coords, &target_rank);

        MPI_Isend(ghost_snd[dir_UP], num_ghost_pts_ver, MPI_Cell, target_rank, MSG_TAG, MPI_COMM_WORLD, &requests_snd[dir_UP]);

        MPI_Irecv(ghost_rcv[dir_UP], num_ghost_pts_ver, MPI_Cell, target_rank, MSG_TAG, MPI_COMM_WORLD, &requests_rcv[dir_UP]);
    }

    /* Down */
    if (comm_flags[dir_DOWN]) {
        offset2 = num_rows_proc*offset1 + NUM_GHOST_ROWS;
        for (i = 0; i < NUM_GHOST_ROWS; i++) {
            memcpy(ghost_snd[dir_DOWN] + i*num_cols_proc, curr_world + offset2 + i*offset1, num_cols_proc*sizeof(Cell));
        }

        target_coords[0] = source_coords[0] + 1;   
        target_coords[1] = source_coords[1];       
        MPI_Cart_rank(CART_COMM, target_coords, &target_rank);

        MPI_Isend(ghost_snd[dir_DOWN], num_ghost_pts_ver, MPI_Cell, target_rank, MSG_TAG, MPI_COMM_WORLD, &requests_snd[dir_DOWN]);

        MPI_Irecv(ghost_rcv[dir_DOWN], num_ghost_pts_ver, MPI_Cell, target_rank, MSG_TAG, MPI_COMM_WORLD, &requests_rcv[dir_DOWN]);
    }

    /* Sending data horizontally */
    /* Left */
    if (comm_flags[dir_LEFT]) {
        offset2 = NUM_GHOST_ROWS*offset1 + NUM_GHOST_ROWS;
        for (i = 0; i < num_rows_proc; i++) {
            memcpy(ghost_snd[dir_LEFT] + i*NUM_GHOST_ROWS, curr_world + offset2 + i*offset1, NUM_GHOST_ROWS*sizeof(Cell));
        }

        target_coords[0] = source_coords[0];   
        target_coords[1] = source_coords[1] - 1;       
        MPI_Cart_rank(CART_COMM, target_coords, &target_rank);

        MPI_Isend(ghost_snd[dir_LEFT], num_ghost_pts_hor, MPI_Cell, target_rank, MSG_TAG, MPI_COMM_WORLD, &requests_snd[dir_LEFT]);

        MPI_Irecv(ghost_rcv[dir_LEFT], num_ghost_pts_hor, MPI_Cell, target_rank, MSG_TAG, MPI_COMM_WORLD, &requests_rcv[dir_LEFT]);
    }
    
    /* Right */
    if (comm_flags[dir_RIGHT]) {
        offset2 = NUM_GHOST_ROWS*offset1 + num_cols_proc;
        for (i = 0; i < num_rows_proc; i++) {
            memcpy(ghost_snd[dir_RIGHT] + i*NUM_GHOST_ROWS, curr_world + offset2 + i*offset1, NUM_GHOST_ROWS*sizeof(Cell));
        }

        target_coords[0] = source_coords[0];   
        target_coords[1] = source_coords[1] + 1;       
        MPI_Cart_rank(CART_COMM, target_coords, &target_rank);

        MPI_Isend(ghost_snd[dir_RIGHT], num_ghost_pts_hor, MPI_Cell, target_rank, MSG_TAG, MPI_COMM_WORLD, &requests_snd[dir_RIGHT]);

        MPI_Irecv(ghost_rcv[dir_RIGHT], num_ghost_pts_hor, MPI_Cell, target_rank, MSG_TAG, MPI_COMM_WORLD, &requests_rcv[dir_RIGHT]);
    }

    /* Sending data diagonally */
    /* Up-Left */
    if (comm_flags[dir_UP_LEFT]) {
        offset2 = NUM_GHOST_ROWS*offset1 + NUM_GHOST_ROWS;
        for (i = 0; i < NUM_GHOST_ROWS; i++) {
            memcpy(ghost_snd[dir_UP_LEFT] + i*NUM_GHOST_ROWS, curr_world + offset2 + i*offset1, NUM_GHOST_ROWS*sizeof(Cell));
        }

        target_coords[0] = source_coords[0] - 1;   
        target_coords[1] = source_coords[1] - 1;       
        MPI_Cart_rank(CART_COMM, target_coords, &target_rank);

        MPI_Isend(ghost_snd[dir_UP_LEFT], num_ghost_pts_dgn, MPI_Cell, target_rank, MSG_TAG, MPI_COMM_WORLD, &requests_snd[dir_UP_LEFT]);

        MPI_Irecv(ghost_rcv[dir_UP_LEFT], num_ghost_pts_dgn, MPI_Cell, target_rank, MSG_TAG, MPI_COMM_WORLD, &requests_rcv[dir_UP_LEFT]);
    }

    /* Up-Right */
    if (comm_flags[dir_UP_RIGHT]) {
        offset2 = NUM_GHOST_ROWS*offset1 + num_cols_proc;
        for (i = 0; i < NUM_GHOST_ROWS; i++) {
            memcpy(ghost_snd[dir_UP_RIGHT] + i*NUM_GHOST_ROWS, curr_world + offset2 + i*offset1, NUM_GHOST_ROWS*sizeof(Cell));
        }

        target_coords[0] = source_coords[0] - 1;   
        target_coords[1] = source_coords[1] + 1;       
        MPI_Cart_rank(CART_COMM, target_coords, &target_rank);

        MPI_Isend(ghost_snd[dir_UP_RIGHT], num_ghost_pts_dgn, MPI_Cell, target_rank, MSG_TAG, MPI_COMM_WORLD, &requests_snd[dir_UP_RIGHT]);

        MPI_Irecv(ghost_rcv[dir_UP_RIGHT], num_ghost_pts_dgn, MPI_Cell, target_rank, MSG_TAG, MPI_COMM_WORLD, &requests_rcv[dir_UP_RIGHT]);
    }

    /* Down-Left */
    if (comm_flags[dir_DOWN_LEFT]) {
        offset2 = num_rows_proc*offset1 + NUM_GHOST_ROWS;
        for (i = 0; i < NUM_GHOST_ROWS; i++) {
            memcpy(ghost_snd[dir_DOWN_LEFT] + i*NUM_GHOST_ROWS, curr_world + offset2 + i*offset1, NUM_GHOST_ROWS*sizeof(Cell));
        }

        target_coords[0] = source_coords[0] + 1;   
        target_coords[1] = source_coords[1] - 1;       
        MPI_Cart_rank(CART_COMM, target_coords, &target_rank);

        MPI_Isend(ghost_snd[dir_DOWN_LEFT], num_ghost_pts_dgn, MPI_Cell, target_rank, MSG_TAG, MPI_COMM_WORLD, &requests_snd[dir_DOWN_LEFT]);

        MPI_Irecv(ghost_rcv[dir_DOWN_LEFT], num_ghost_pts_dgn, MPI_Cell, target_rank, MSG_TAG, MPI_COMM_WORLD, &requests_rcv[dir_DOWN_LEFT]);
    }

    /* Down-Right */
    if (comm_flags[dir_DOWN_RIGHT]) {
        offset2 = num_rows_proc*offset1 + num_cols_proc;
        for (i = 0; i < NUM_GHOST_ROWS; i++) {
            memcpy(ghost_snd[dir_DOWN_RIGHT] + i*NUM_GHOST_ROWS, curr_world + offset2 + i*offset1, NUM_GHOST_ROWS*sizeof(Cell));
        }

        target_coords[0] = source_coords[0] + 1;   
        target_coords[1] = source_coords[1] + 1;       
        MPI_Cart_rank(CART_COMM, target_coords, &target_rank);

        MPI_Isend(ghost_snd[dir_DOWN_RIGHT], num_ghost_pts_dgn, MPI_Cell, target_rank, MSG_TAG, MPI_COMM_WORLD, &requests_snd[dir_DOWN_RIGHT]);

        MPI_Irecv(ghost_rcv[dir_DOWN_RIGHT], num_ghost_pts_dgn, MPI_Cell, target_rank, MSG_TAG, MPI_COMM_WORLD, &requests_rcv[dir_DOWN_RIGHT]);
    }

    /* Wait for all messages to be sent and received*/
    MPI_Waitall(8, requests_rcv, MPI_STATUSES_IGNORE);
    MPI_Waitall(8, requests_snd, MPI_STATUSES_IGNORE);
    
    /* Copy all ghost points received to the current world */
    /* Up */
    if (comm_flags[dir_UP]) {
        offset2 = NUM_GHOST_ROWS;
        for (i = 0; i < NUM_GHOST_ROWS; i++) {
            memcpy(curr_world + offset2 + i*offset1, ghost_rcv[dir_UP] + i*num_cols_proc, num_cols_proc*sizeof(Cell));
            memcpy(next_world + offset2 + i*offset1, ghost_rcv[dir_UP] + i*num_cols_proc, num_cols_proc*sizeof(Cell));
        }
    }

    /* Down */
    if (comm_flags[dir_DOWN]) {
        offset2 = (num_rows_proc + NUM_GHOST_ROWS)*offset1 + NUM_GHOST_ROWS;
        for (i = 0; i < NUM_GHOST_ROWS; i++) {
            memcpy(curr_world + offset2 + i*offset1, ghost_rcv[dir_DOWN] + i*num_cols_proc, num_cols_proc*sizeof(Cell));
            memcpy(next_world + offset2 + i*offset1, ghost_rcv[dir_DOWN] + i*num_cols_proc, num_cols_proc*sizeof(Cell));
        }
    }

    /* Left */
    if (comm_flags[dir_LEFT]) {
        offset2 = NUM_GHOST_ROWS*offset1;
        for (i = 0; i < num_rows_proc; i++) {
            memcpy(curr_world + offset2 + i*offset1, ghost_rcv[dir_LEFT] + i*NUM_GHOST_ROWS, NUM_GHOST_ROWS*sizeof(Cell));
            memcpy(next_world + offset2 + i*offset1, ghost_rcv[dir_LEFT] + i*NUM_GHOST_ROWS, NUM_GHOST_ROWS*sizeof(Cell));
        }
    }
    
    /* Right */
    if (comm_flags[dir_RIGHT]) {
        offset2 = NUM_GHOST_ROWS*offset1 + num_cols_proc + NUM_GHOST_ROWS;
        for (i = 0; i < num_rows_proc; i++) {
            memcpy(curr_world + offset2 + i*offset1, ghost_rcv[dir_RIGHT] + i*NUM_GHOST_ROWS, NUM_GHOST_ROWS*sizeof(Cell));
            memcpy(next_world + offset2 + i*offset1, ghost_rcv[dir_RIGHT] + i*NUM_GHOST_ROWS, NUM_GHOST_ROWS*sizeof(Cell));
        }
    }

    /* Up-Left */
    if (comm_flags[dir_UP_LEFT]) {
        offset2 = 0;
        for (i = 0; i < NUM_GHOST_ROWS; i++) {
            memcpy(curr_world + offset2 + i*offset1, ghost_rcv[dir_UP_LEFT] + i*NUM_GHOST_ROWS, NUM_GHOST_ROWS*sizeof(Cell));
            memcpy(next_world + offset2 + i*offset1, ghost_rcv[dir_UP_LEFT] + i*NUM_GHOST_ROWS, NUM_GHOST_ROWS*sizeof(Cell));
        }
    }

    /* Up-Right */
    if (comm_flags[dir_UP_RIGHT]) {
        offset2 = num_cols_proc + NUM_GHOST_ROWS;
        for (i = 0; i < NUM_GHOST_ROWS; i++) {
            memcpy(curr_world + offset2 + i*offset1, ghost_rcv[dir_UP_RIGHT] + i*NUM_GHOST_ROWS, NUM_GHOST_ROWS*sizeof(Cell));
            memcpy(next_world + offset2 + i*offset1, ghost_rcv[dir_UP_RIGHT] + i*NUM_GHOST_ROWS, NUM_GHOST_ROWS*sizeof(Cell));
        }
    }

    /* Down-Left */
    if (comm_flags[dir_DOWN_LEFT]) {
        offset2 = (num_rows_proc + NUM_GHOST_ROWS)*offset1;
        for (i = 0; i < NUM_GHOST_ROWS; i++) {
            memcpy(curr_world + offset2 + i*offset1, ghost_rcv[dir_DOWN_LEFT] + i*NUM_GHOST_ROWS, NUM_GHOST_ROWS*sizeof(Cell));
            memcpy(next_world + offset2 + i*offset1, ghost_rcv[dir_DOWN_LEFT] + i*NUM_GHOST_ROWS, NUM_GHOST_ROWS*sizeof(Cell));
        }
    }

    /* Down-Right */
    if (comm_flags[dir_DOWN_RIGHT]) {
        offset2 = (num_rows_proc + NUM_GHOST_ROWS)*offset1 + num_cols_proc + NUM_GHOST_ROWS;
        for (i = 0; i < NUM_GHOST_ROWS; i++) {
            memcpy(curr_world + offset2 + i*offset1, ghost_rcv[dir_DOWN_RIGHT] + i*NUM_GHOST_ROWS, NUM_GHOST_ROWS*sizeof(Cell));
            memcpy(next_world + offset2 + i*offset1, ghost_rcv[dir_DOWN_RIGHT] + i*NUM_GHOST_ROWS, NUM_GHOST_ROWS*sizeof(Cell));
        }
    }
}

/* FUNCTION: print_results
 * ARGUMENTS: Cell *world       - World to count rabbits, foxes and rocks from
              int local_rows    - Number of rows on local process
              int local_cols    - Number of rocolumnsws on local process
              int id            - Process id
 * DESCRIPTION: Prints the final results to stdout
 * RETURN: N/A
 */
void print_results(Cell *world, int local_rows, int local_cols, int id){
    
    int i, j;                               /*Iteration variables*/
    int c = 0;
    int local_rabbits = 0;                  /*Final number of rabbits for local process*/
    int local_foxes = 0;                    /*Final number of foxes for local process*/
    int local_rocks = 0;                    /*Final number of rocks for local process*/
    int global_rabbits;                     /*Final number of rabbits*/
    int global_foxes;                       /*Final number of foxes*/
    int global_rocks;                       /*Final number of rocks*/

    /*Count foxes, rabbits and rocks for local process*/
    for(i = NUM_GHOST_ROWS; i < local_rows - NUM_GHOST_ROWS; i++){  

        for(j = NUM_GHOST_ROWS; j < local_cols - NUM_GHOST_ROWS; j++){

            c = i*local_cols + j;

            if(world[c].animal == FOX) {
                local_foxes++;
            }
            else if(world[c].animal == RABBIT) {
                local_rabbits++;
            }
            else if(world[c].animal == ROCK) {
                local_rocks++;
            }
        }
    }

    /*Get total number of rocks, rabbits and foxes*/
    MPI_Reduce(&local_rocks, &global_rocks, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_rabbits, &global_rabbits, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_foxes, &global_foxes, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    /*Process 0 prints results*/
    if(id == 0)
    {
        fprintf(stdout, "%d %d %d\n",global_rocks, global_rabbits, global_foxes);   /*Print result to stdout*/
    }
}

int main(int argc, char *argv[]){

    int generations = 0;                        /*Number of generations to simulate*/
    int M = 0;                                  /*Number of rows of the world matrix*/
    int N = 0;                                  /*Number of columns of the world matrix*/
    int num_rocks = 0;                          /*Number of rocks in the world matrix*/
    int num_rabbits = 0;                        /*Number of rabbits in the world matrix*/
    int num_foxes = 0;                          /*Number of foxes in the world matrix*/
    int breed_age_r = 0;                        /*Breeding age for rabbits*/
    int breed_age_f = 0;                        /*Breeding age for foxes*/
    int starve_age_f = 0;                       /*Starving age for foxes*/
    double exec_time = 0;                       /*Execution time of simulation*/
    double global_exec_time = 0;                /*Final execution time*/
    uint32_t seed = 0;                          /*Seed for randomization*/
    int i = 0;                                  /*Iteration variables*/
    int id, p;                                  /*Process id and number of processes*/

    int *updated;                               /*Matrix for marking already updated cells in a sub-generation*/
    Cell *curr_world;                           /*Matrix for storing current world state*/
    Cell *next_world;                           /*Matrix for storing next world state*/

    Cell **ghost_rcv;                           /*Buffer for receiving ghost rows and columns*/
    Cell **ghost_snd;                           /*Buffer for sending ghost rows and columns*/
    
    int first_row_proc;                         /*First row for a certain process*/
    int last_row_proc;                          /*Last row for a certain process*/
    int num_rows_proc;                          /*Number of rows per process*/

    int first_col_proc;                         /*First column for a certain process*/
    int last_col_proc;                          /*Last column for a certain process*/
    int num_cols_proc;                          /*Number of columns per process*/

    int local_rows, local_cols;                 /*Number of rows and columns for a local process*/
    int total_rows, total_cols;                 /*Total number of rows and columns*/
    
    int num_ghost_pts_ver;                      /*Number of ghost points to communicate in a vertical direction*/
    int num_ghost_pts_hor;                      /*Number of ghost points to communicate in a horizontal direction*/
    int num_ghost_pts_dgn;                      /*Number of ghost points to communicate in a diagonal direction*/
    

    MPI_Datatype MPI_Cell;                      /*Custom MPI datatype*/
    MPI_Datatype types[3] = {MPI_CHAR, MPI_INT, MPI_INT};

    MPI_Aint displacements[3];
    MPI_Aint base_address;

    Cell dummy_cell;
    int lengths[3] = {1, 1, 1};
    int dims[2] = {0, 0};

    if (argc != 11) {                        /*Verifies if number of args is correct*/
        exit(-1);
    }

    for(i = 1; i <= 10; i++) {               /*Verifies if all args are positive integers and assigns them to the correct variables*/
        if(atoi(argv[i]) < 0) {
            exit(-1);
        }
        
        switch (i) {
            case 1:
                generations = atoi(argv[i]);
                break;
            case 2:
                M = atoi(argv[i]);
                break;
            case 3:
                N = atoi(argv[i]);
                break;
            case 4:
                num_rocks = atoi(argv[i]);
                break;
            case 5:
                num_rabbits = atoi(argv[i]);
                break;
            case 6:
                breed_age_r = atoi(argv[i]);
                break;
            case 7:
                num_foxes = atoi(argv[i]);
                break;
            case 8:
                breed_age_f = atoi(argv[i]);
                break;
            case 9:
                starve_age_f = atoi(argv[i]);
                break;
            case 10:
                seed = atoi(argv[i]);
                break;
            default:
                exit(-1);
        }
    }

    MPI_Init(&argc, &argv);                              /*Initialize MPI*/
    MPI_Comm_rank(MPI_COMM_WORLD, &id);                  /*Get process ID*/
    MPI_Comm_size(MPI_COMM_WORLD, &p);                   /*Get number of processes*/

    /*Variables for creating process grid*/
    MPI_Comm CART_COMM;
    int reorder = 1;
    int periods[2] = {1, 1};
    int id_cart;
    int id_coords[2];

    MPI_Dims_create(p, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &CART_COMM);
    MPI_Comm_rank(CART_COMM, &id_cart);
    MPI_Cart_coords(CART_COMM, id_cart, 2, id_coords);                       

    first_row_proc = BLOCK_LOW(id_coords[0], dims[0], M);                 /*Get first row*/
    last_row_proc = BLOCK_HIGH(id_coords[0], dims[0], M);                 /*Get last row*/
    num_rows_proc = BLOCK_SIZE(id_coords[0], dims[0], M);                 /*Get number of rows for process*/

    first_col_proc = BLOCK_LOW(id_coords[1], dims[1], N);                 /*Get first column*/
    last_col_proc = BLOCK_HIGH(id_coords[1], dims[1], N);                 /*Get last column*/
    num_cols_proc = BLOCK_SIZE(id_coords[1], dims[1], N);                 /*Get number of columns for process*/

    local_rows = num_rows_proc + 2*NUM_GHOST_ROWS;                        /*Number of rows for local process including ghost rows*/
    local_cols = num_cols_proc + 2*NUM_GHOST_ROWS;                        /*Number of columns for local process including ghost columns*/

    total_rows = M;                                                       /*Copy M, N variables for easier understanding*/
    total_cols = N;

    /*Allocate worlds*/
    curr_world = (Cell *) malloc(local_rows * local_cols * sizeof(Cell));
    next_world = (Cell *) malloc(local_rows * local_cols * sizeof(Cell));
    updated = (int *) calloc(local_rows * local_cols, sizeof(int));        /*Create 'updated' matrix to mark already moved animal*/
    
    if(curr_world == NULL || next_world == NULL || updated == NULL) {
        exit(-1);
    }

    initialize_worlds(curr_world, next_world, local_rows, local_cols);    /*Initialize worlds*/
    
    /*Create MPI Datatype*/
    MPI_Get_address(&dummy_cell, &base_address);
    MPI_Get_address(&dummy_cell.animal, &displacements[0]);
    MPI_Get_address(&dummy_cell.breed_age, &displacements[1]);
    MPI_Get_address(&dummy_cell.starve_age, &displacements[2]);

    displacements[0] = MPI_Aint_diff(displacements[0], base_address);
    displacements[1] = MPI_Aint_diff(displacements[1], base_address);
    displacements[2] = MPI_Aint_diff(displacements[2], base_address);

    MPI_Type_create_struct(3, lengths, displacements, types, &MPI_Cell);
    MPI_Type_commit(&MPI_Cell);

    num_ghost_pts_ver = num_cols_proc * NUM_GHOST_ROWS;                 /*Number of ghost points to communicate in a vertical direction*/
    num_ghost_pts_hor = num_rows_proc * NUM_GHOST_ROWS;                 /*Number of ghost points to communicate in a horizontal direction*/
    num_ghost_pts_dgn = NUM_GHOST_ROWS * NUM_GHOST_ROWS;                /*Number of ghost points to communicate in a diagonal direction*/

    /* Allocating ghost point vectors */
    ghost_rcv = (Cell **) malloc(8 * sizeof(Cell *));
    ghost_snd = (Cell **) malloc(8 * sizeof(Cell *));

    /* Checking if memory allocations were successful */
    if (ghost_rcv == NULL || ghost_snd == NULL) {
        exit(-1);
    }

    ghost_rcv[dir_UP] = (Cell *) malloc(num_ghost_pts_ver * sizeof(Cell));
    ghost_snd[dir_UP] = (Cell *) malloc(num_ghost_pts_ver * sizeof(Cell));

    ghost_rcv[dir_DOWN] = (Cell *) malloc(num_ghost_pts_ver * sizeof(Cell));
    ghost_snd[dir_DOWN] = (Cell *) malloc(num_ghost_pts_ver * sizeof(Cell));

    ghost_rcv[dir_LEFT] = (Cell *) malloc(num_ghost_pts_hor * sizeof(Cell));
    ghost_snd[dir_LEFT] = (Cell *) malloc(num_ghost_pts_hor * sizeof(Cell));

    ghost_rcv[dir_RIGHT] = (Cell *) malloc(num_ghost_pts_hor * sizeof(Cell));
    ghost_snd[dir_RIGHT] = (Cell *) malloc(num_ghost_pts_hor * sizeof(Cell));

    ghost_rcv[dir_UP_LEFT] = (Cell *) malloc(num_ghost_pts_dgn * sizeof(Cell));
    ghost_snd[dir_UP_LEFT] = (Cell *) malloc(num_ghost_pts_dgn * sizeof(Cell));

    ghost_rcv[dir_UP_RIGHT] = (Cell *) malloc(num_ghost_pts_dgn * sizeof(Cell));
    ghost_snd[dir_UP_RIGHT] = (Cell *) malloc(num_ghost_pts_dgn * sizeof(Cell));

    ghost_rcv[dir_DOWN_LEFT] = (Cell *) malloc(num_ghost_pts_dgn * sizeof(Cell));
    ghost_snd[dir_DOWN_LEFT] = (Cell *) malloc(num_ghost_pts_dgn * sizeof(Cell));

    ghost_rcv[dir_DOWN_RIGHT] = (Cell *) malloc(num_ghost_pts_dgn * sizeof(Cell));
    ghost_snd[dir_DOWN_RIGHT] = (Cell *) malloc(num_ghost_pts_dgn * sizeof(Cell));

    /* Checking if memory allocations were successful */
    for (i = 0; i < 8; i++) {
        if(ghost_rcv[i] == NULL || ghost_snd[i] == NULL) {
            exit(-1);
        }
    }
    
    /*Generate elements to fill world*/
    generate_element(curr_world, num_rocks, ROCK, &seed, M, N, first_row_proc, last_row_proc, first_col_proc, last_col_proc, local_cols);
    generate_element(curr_world, num_rabbits, RABBIT, &seed, M, N, first_row_proc, last_row_proc, first_col_proc, last_col_proc, local_cols);
    generate_element(curr_world, num_foxes, FOX, &seed, M, N, first_row_proc, last_row_proc, first_col_proc, last_col_proc, local_cols);
    
    memcpy(next_world, curr_world, sizeof(Cell)*local_rows*local_cols);     /*Copy world*/
    
    exec_time = -omp_get_wtime();                                           /*Start timing simulation*/

#pragma omp parallel                                                        /*Create parallel region*/
{
    int i, j, k, g;
    int c = 0;                                          /*Cell number*/

    int num_threads = omp_get_num_threads();            /*Number of running threads*/     
    int num_tasks = 8*num_threads;                      /*Number of tasks to launch*/
    int L = (local_rows + num_tasks - 1)/num_tasks;     /*The amount of work (rows) for each task to process*/

    int start, end;                                     /*Variables that mark the first and last index of the rows each task should compute*/
    int global_row, global_col;

    int flags[num_tasks+1];                             /*Vector of flag variables to establish dependencies between tasks*/                           

    /*Generation loop*/
    for(g = 1; g <= generations; g++){

        /*Receive and send ghost rows*/
        #pragma omp single
        communicate_ghost_points(curr_world, next_world, ghost_snd, ghost_rcv, num_ghost_pts_ver, num_ghost_pts_hor, num_ghost_pts_dgn, \
                            num_rows_proc, num_cols_proc, CART_COMM, MPI_Cell, dims, id_coords);

        /*Red Sub-Generation*/
        #pragma omp single
        {
            /* Creating tasks with even id */
            for(k = 0; k < num_tasks; k+=2) {

                #pragma omp task depend(out:flags[k])                        /*In order to avoid race condition between consecutive task borders*/
                {                                             
                    (k == 0)? (start = 1): (start = k*L);
                    end = (k+1)*L;
                
                    for(i = start; i < end && i < local_rows-1; i++){        /*Iterate over all rows*/

                        global_row = i + first_row_proc - NUM_GHOST_ROWS;                      /*Calculate global row index*/                         
                        if (global_row < 0 || global_row > total_rows-1) {                     /*Check if it is inbounds*/
                            continue;
                        }  
                                                                            
                        for(j = (global_row + first_col_proc)%2 + 1; j < local_cols-1; j+=2){        /*Process row's red cells*/

                            global_col = j + first_col_proc - NUM_GHOST_ROWS;                        /*Calculate global column index*/ 
                            if (global_col < 0 || global_col > total_cols-1) {                       /*Check if it is inbounds*/
                                continue;
                            }

                            c = i*local_cols + j;                                                                   
                            if(curr_world[c].animal == FOX || curr_world[c].animal == RABBIT)        /*Only try to move if the cell has a fox or rabbit*/
                                movement(curr_world, next_world, updated, local_rows, local_cols, \
                                        i, j, breed_age_f, breed_age_r,                          \
                                        global_row, global_col, total_rows, total_cols);        
                        }
                    }
                }
            }

            /* Creating tasks with odd id */
            for(k = 1; k < num_tasks; k+=2) {

            #pragma omp task depend(in:flags[k-1], flags[k+1])                        /*In order to avoid race condition between*/
                {                                                                     /*consecutive tasks borders*/
                    start = k*L;
                    end = (k+1)*L;
                
                    for(i = start; i < end && i < local_rows-1; i++){                 /*Iterate over all rows*/

                        global_row = i + first_row_proc - NUM_GHOST_ROWS;             /*Calculate global row index*/                         
                        if (global_row < 0 || global_row > total_rows-1) {            /*Check if it is inbounds*/
                            continue;
                        }  
                                                                            
                        for(j = (global_row + first_col_proc)%2 + 1; j < local_cols-1; j+=2){        /*Process row's red cells*/

                            global_col = j + first_col_proc - NUM_GHOST_ROWS;                        /*Calculate global column index*/ 
                            if (global_col < 0 || global_col > total_cols-1) {                       /*Check if it is inbounds*/
                                continue;
                            }

                            c = i*local_cols + j;                                                                   
                            if(curr_world[c].animal == FOX || curr_world[c].animal == RABBIT)         /*Only try to move if the cell has a fox or rabbit*/
                                movement(curr_world, next_world, updated, local_rows, local_cols, \
                                        i, j, breed_age_f, breed_age_r,                          \
                                        global_row, global_col, total_rows, total_cols);        
                        }
                    }
                }
            }
        }

        /*Copy world*/
        #pragma omp single
        memcpy(curr_world, next_world, local_rows*local_cols*sizeof(Cell));


        /*Black Sub-Generation*/
        #pragma omp single
        {
            /* Creating tasks with even id */
            for(k = 0; k < num_tasks; k+=2) {

                #pragma omp task depend(out:flags[k])                                   /*In order to avoid race condition between*/
                {                                                                       /*consecutive tasks borders*/
                    (k == 0)? (start = 2): (start = k*L);
                    end = (k+1)*L;
                
                    for(i = start; i < end && i < local_rows-2; i++){                   /*Iterate over all rows*/

                        global_row = i + first_row_proc - NUM_GHOST_ROWS;                      /*Calculate global row index*/                         
                        if (global_row < 0 || global_row > total_rows-1) {                     /*Check if it is inbounds*/
                            continue;
                        }  
                                                                            
                        for(j = (global_row + first_col_proc)%2 + 2; j < local_cols-2; j+=2){        /*Process row's red cells*/

                            global_col = j + first_col_proc - NUM_GHOST_ROWS;                        /*Calculate global column index*/ 
                            if (global_col < 0 || global_col > total_cols-1) {                       /*Check if it is inbounds*/
                                continue;
                            }


                            c = i*local_cols + j; 
                            if(updated[c] == 1)                                         /*Check if cell was updated on previous sub-generation*/                                                    
                                    continue;                                           /*If so, skip*/

                            if(curr_world[c].animal == FOX || curr_world[c].animal == RABBIT)         /*Only try to move if the cell has a fox or rabbit*/
                                movement(curr_world, next_world, updated, local_rows, local_cols, \
                                        i, j, breed_age_f, breed_age_r,                          \
                                        global_row, global_col, total_rows, total_cols);        
                        }
                    }
                }
            }

            /* Creating tasks with odd id */
            for(k = 1; k < num_tasks; k+=2) {

                #pragma omp task depend(in:flags[k-1], flags[k+1])                      /*In order to avoid race condition between*/
                {                                                                       /*consecutive tasks borders*/
                    start = k*L;
                    end = (k+1)*L;
                
                    for(i = start; i < end && i < local_rows-2; i++){                   /*Iterate over all rows*/

                        global_row = i + first_row_proc - NUM_GHOST_ROWS;                      /*Calculate global row index*/                         
                        if (global_row < 0 || global_row > total_rows-1) {                     /*Check if it is inbounds*/
                            continue;
                        }  
                                                                            
                        for(j = (global_row + first_col_proc)%2 + 2; j < local_cols-2; j+=2){        /*Process row's red cells*/

                            global_col = j + first_col_proc - NUM_GHOST_ROWS;                        /*Calculate global column index*/ 
                            if (global_col < 0 || global_col > total_cols-1) {                       /*Check if it is inbounds*/
                                continue;
                            }

                            c = i*local_cols + j;
                            if(updated[c] == 1)                                         /*Check if cell was updated on previous sub-generation*/                             
                                    continue;                                           /*If so, skip*/

                            if(curr_world[c].animal == FOX || curr_world[c].animal == RABBIT)         /*Only try to move if the cell has a fox or rabbit*/
                                movement(curr_world, next_world, updated, local_rows, local_cols, \
                                        i, j, breed_age_f, breed_age_r,                          \
                                        global_row, global_col, total_rows, total_cols);        
                        }
                    }
                }
            }
        }
        
        /*Copy world and remove foxes with starving ages greater than the limit*/
        #pragma omp for
        for(i = 0; i < local_rows; i++){
            
            global_row = i + first_row_proc - NUM_GHOST_ROWS;              /*Calculate global row index*/                         
            if (global_row < 0 || global_row > total_rows-1) {             /*Check if it is inbounds*/
                continue;
            }  

            for(j = 0; j < local_cols; j++){

                global_col = j + first_col_proc - NUM_GHOST_ROWS;           /*Calculate global column index*/ 
                if (global_col < 0 || global_col > total_cols-1) {          /*Check if it is inbounds*/
                    continue;
                }

                c = i*local_cols + j;     
                if(next_world[c].animal == FOX) {                           /*If the current position is a fox, check its starve age*/
                    if(next_world[c].starve_age >= starve_age_f) {          /*If it is greater or equal to the limit*/
                        next_world[c].animal = EMPTY;                       /*Remove it*/
                    }
                }

                updated[c] = 0;                                             /*Reset updated vector*/
                curr_world[c].animal = next_world[c].animal;                /*Copy to curr_world*/
                curr_world[c].breed_age = next_world[c].breed_age;
                curr_world[c].starve_age = next_world[c].starve_age;
            }
        }
    }
}

    exec_time += omp_get_wtime();                                   /*Finish timing simulation*/

    /*Get the execution time from the process that took longest*/
    MPI_Reduce(&exec_time, &global_exec_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); 

    if(id == 0)
        fprintf(stderr, "%.5fs\n", global_exec_time);              /*Print execution time to stderr*/
    
    print_results(curr_world, local_rows, local_cols, id);         /*Print results to stdout*/

    /*Free allocated memory*/
    free(curr_world);
    free(next_world);
    free(updated);

    for (i = 0; i < 8; i++) {
        free(ghost_rcv[i]);
        free(ghost_snd[i]);
    }

    free(ghost_rcv);
    free(ghost_snd);
    
    /*Terminate MPI*/
    MPI_Finalize();

    return 0;
}