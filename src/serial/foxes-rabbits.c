/*****************************************
* Authors:	João Leitão		 93088       *
*		  	José Brito		 93106       *
*  		   	Xavier Fernandes 93202       *
*****************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>

#define ROCK    '*'
#define RABBIT  'R'
#define FOX     'F'
#define EMPTY   ' '

#define dir_UP      0
#define dir_RIGHT   1
#define dir_DOWN    2
#define dir_LEFT    3
#define dir_CENTER  4

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
 * ARGUMENTS: Cell **curr_world - Current world state
              int n             - Number of rocks, rabbits or foxes
              char atype        - Character representative of rocks, rabbits or foxes
              uint32_t *seed    - Seed for randomization
              int M             - Number of rows
              int N             - Number of columns
 * DESCRIPTION: Randomly distributes elements (rocks, foxes, rabbits) throughout
                the world matrix
 * RETURN: N/A
 */
void generate_element(Cell **curr_world, int n, char atype, uint32_t *seed, int M, int N){
    
    int i, j, k;
    for(k = 0; k < n; k++){
        i = M * r4_uni(seed);
        j = N * r4_uni(seed);
        if(curr_world[i][j].animal == ' ')
        {
            curr_world[i][j].animal = atype;
            
            if(atype != ROCK)
            {
                curr_world[i][j].breed_age = 0;
                if(atype == FOX){
                    curr_world[i][j].starve_age = 0;
                }
            }
        }
    }
}

/* FUNCTION: initialize_worlds
 * ARGUMENTS: Cell **curr_world - Current world state
              Cell **next_world - Next world state
              int M             - Number of rows
              int N             - Number of columns

 * DESCRIPTION: Initialize worlds
 * RETURN: N/A
 */
void initialize_worlds(Cell **curr_world, Cell **next_world, int M, int N) {

    for(int i = 0; i < M; i++){

        curr_world[i] = (Cell *) malloc(N * sizeof(Cell));  /*Allocating memory to represent the worlds*/
        if(curr_world[i] == NULL) {                         /*Current world is static, used to compute its next state*/
            exit(-1);                                       /*Next world will be updated with the new state*/
        }

        next_world[i] = (Cell *) malloc(N * sizeof(Cell));
        if(next_world[i] == NULL) {
            exit(-1);
        }

        for(int j = 0; j < N; j++) {                /*Initializing all entries as empty cells*/

            curr_world[i][j].animal = ' ';
            curr_world[i][j].breed_age = -1;
            curr_world[i][j].starve_age = -1;

            next_world[i][j].animal = ' ';
            next_world[i][j].breed_age = -1;
            next_world[i][j].starve_age = -1;
                
        }
    }
}


/* FUNCTION: calculate_direction
 * ARGUMENTS: int x      - Number of possible moves
              int *moves - Vector marking potential movement directions
              int c      - Cell number (i*N + j)

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
 * ARGUMENTS: Cell **curr_world - Current world state
              int i             - Current cell's row number
              int j             - Current cell's column number
              int M             - Number of rows
              int N             - Number of collumns
              int *p            - Number of adjacent empty cells
              int *empty_cells  - Vector marking available adjacent empty cells

 * DESCRIPTION: Counts the number os adjacent empty cells and registers them in a vector
                to find the number of possible movements for a rabbit cell
 * RETURN: N/A
 */
void get_moves_rabbit(Cell **curr_world, int i, int j, int M, int N,int *p, int *empty_cells){

    curr_world[i][j].breed_age++;                   /*Increment rabbit breeding age*/
    if(i - 1 >= 0){                                 /*Check if position above is in bounds*/
        if(curr_world[i-1][j].animal == EMPTY){     /*If position above is empty*/
            empty_cells[dir_UP] = 1;                /*Mark that in the vector*/
            *p+=1;                                  /*Increment number of empty cells available*/
        }   
    }
    if(j + 1 < N){                                  /*Check if position to the left is in bounds*/
        if(curr_world[i][j+1].animal == EMPTY){     /*If position to the right is empty*/
            empty_cells[dir_RIGHT] = 1;             /*Mark that in the vector*/  
            *p+=1;                                  /*Increment number of empty cells available*/
        }
    }
    if(i + 1 < M){                                  /*Check if position below is in bounds*/
        if(curr_world[i+1][j].animal == EMPTY){     /*If position below is empty*/
            empty_cells[dir_DOWN] = 1;              /*Mark that in the vector*/ 
            *p+=1;                                  /*Increment number of empty cells available*/
        }   
    }
    if(j - 1 >= 0){                                 /*Check if position to the left is in bounds*/
        if(curr_world[i][j-1].animal == EMPTY){     /*If position to the left is empty*/
            empty_cells[dir_LEFT] = 1;              /*Mark that in the vector*/ 
            *p+=1;                                  /*Increment number of empty cells available*/
        }
    }
}


/* FUNCTION: get_moves_fox
 * ARGUMENTS: Cell **curr_world     - Current world state
              int i                 - Current cell's row number
              int j                 - Current cell's column number
              int M                 - Number of rows
              int N                 - Number of columns
              int *p                - Number of adjacent empty cells
              int *adj_rabbits      - Number of adjacent rabbit cells
              int *empty_cells      - Vector marking available adjacent empty cells
              int *adj_rabbit_cells - Vector marking available adjacent rabbit cells

 * DESCRIPTION: Searches adjacent cells, prioritising rabbits and then empty cells, to find the number of possible movements for a fox cell
 * RETURN: N/A
 */
void get_moves_fox(Cell **curr_world, int i, int j, int M, int N,int *p, int *adj_rabbit, int *empty_cells,int *adj_rabbit_cells){
   
    curr_world[i][j].breed_age++;                   /*Increment fox breeding age*/
    curr_world[i][j].starve_age++;                  /*Increment fox starving age*/

    if(i - 1 >= 0){                                 /*Check if position above is in bounds*/
        if(curr_world[i-1][j].animal == RABBIT){    /*If position above is a rabbit, mark it*/
            adj_rabbit_cells[dir_UP] = 1;           
            *adj_rabbit+=1;                         /*Increment number of adjacent rabbits*/
        }
        else if(curr_world[i-1][j].animal == EMPTY) /*Else check if it is empty*/
        {
            empty_cells[dir_UP] = 1;                /*Mark that in the vector*/ 
            *p+=1;                                  /*Increment number of empty cells available*/
        }
    }
    if(j + 1 < N){                                  /*Check if position to the left is in bounds*/
        if(curr_world[i][j+1].animal == RABBIT){    /*If position to the right is a rabbit, mark it*/
            adj_rabbit_cells[dir_RIGHT] = 1;
            *adj_rabbit+=1;                         /*Increment number of adjacent rabbits*/
        }           
        else if(curr_world[i][j+1].animal == EMPTY) /*Else check if it is empty*/
        {
            empty_cells[dir_RIGHT] = 1;             /*Mark that in the vector*/ 
            *p+=1;                                  /*Increment number of empty cells available*/
        }   
    }
    if(i + 1 < M){                                  /*Check if position below is in bounds*/
        if(curr_world[i+1][j].animal == RABBIT){    /*If position below is a rabbit, mark it*/
            adj_rabbit_cells[dir_DOWN] = 1;
            *adj_rabbit+=1;                         /*Increment number of adjacent rabbits*/
        }
        else if(curr_world[i+1][j].animal == EMPTY) /*Else check if it is empty*/
        {
            empty_cells[dir_DOWN] = 1;              /*Mark that in the vector*/ 
            *p+=1;                                  /*Increment number of empty cells available*/
        }
    }
    if(j - 1 >= 0){                                 /*Check if position to the left is in bounds*/
        if(curr_world[i][j-1].animal == RABBIT){    /*If position to the left is a rabbit, mark it*/
            adj_rabbit_cells[dir_LEFT] = 1;
            *adj_rabbit+=1;                         /*Increment number of adjacent rabbits*/
        }
        else if(curr_world[i][j-1].animal == EMPTY) /*Else check if it is empty*/
        {
            empty_cells[dir_LEFT] = 1;              /*Mark that in the vector*/ 
            *p+=1;                                  /*Increment number of empty cells available*/
        }
    }
}

/* FUNCTION: solve_conflicts_and_breeding
 * ARGUMENTS: Cell **curr_world - Current world state
              Cell **next_world - Next world state
              int **updated     - Matrix marking already updated cells
              int i             - Cell's current row
              int j             - Cell's current column
              int i2            - Cell's next row
              int j2            - Cell's next column
              int breed_age_f   - Breeding age for foxes
              int breed_age_r   - Breeding age for rabbits

 * DESCRIPTION: Solves conflicts between animal in a particular cell, while also determining whether breeding occurs
 * RETURN: N/A
 */
void solve_conflicts_and_breeding(Cell **curr_world, Cell **next_world, int **updated, int i, int j, int i2, int j2, int breed_age_f, int breed_age_r){
    
    int conflict_type = 0;                  /*Conflict type identifier variable*/
    int breeding_flag = 0;                  /*Flag informing if breeding occurs*/

    Cell tmp;                               /*Auxiliar Cell variable*/

    next_world[i][j].animal = EMPTY;        /*Setting the previous entry to empty because animal moved*/
    next_world[i][j].breed_age = -1;        
    next_world[i][j].starve_age = -1;

    tmp.animal = curr_world[i][j].animal;           /*Saving cell that is moving in auxiliar variable*/
    tmp.breed_age = curr_world[i][j].breed_age;
    tmp.starve_age = curr_world[i][j].starve_age;

    if(curr_world[i][j].animal == FOX && curr_world[i][j].breed_age >= breed_age_f) {       /*Checks if a fox reached or surpassed breeding age*/
        tmp.breed_age = 0;                                                                  /*Resets breeding age*/
        breeding_flag = 1;                                                                  /*Marks flag to know breeding must happen*/
    }

    if(curr_world[i][j].animal == RABBIT && curr_world[i][j].breed_age >= breed_age_r) {    /*Checks if a rabbit reached or surpassed breeding age*/
        tmp.breed_age = 0;                                                                  /*Resets breeding age*/
        breeding_flag = 1;                                                                  /*Marks flag to know breeding must happen*/
    }
    
    if(breeding_flag) {                            /*Checks if breeding must happen*/
        next_world[i][j].animal = tmp.animal;      /*Leaves the bred animal behind*/
        next_world[i][j].breed_age = 0;            /*Bred animal begins with 0 breed age*/
        next_world[i][j].starve_age = -1;          /*Rabbits do not have hunger*/

        if(curr_world[i][j].animal == FOX)         
            next_world[i][j].starve_age = 0;       /*If it's a fox than it begins with 0 starvation*/
    }
    
    /*Conflict type identification*/
    if(tmp.animal == FOX && next_world[i2][j2].animal == FOX)                 /*Conflict between two foxes*/
        conflict_type = 1;
    else if (tmp.animal == RABBIT && next_world[i2][j2].animal == RABBIT)     /*Conflict between two rabbits*/
        conflict_type = 2;
    else if (tmp.animal == RABBIT && next_world[i2][j2].animal == FOX)        /*Conflict between a fox and a rabbit, rabbits bumps into fox*/
        conflict_type = 3;
    else if (tmp.animal == FOX && next_world[i2][j2].animal == RABBIT)        /*Conflict between a fox and a rabbit, fox bumps into rabbit*/
        conflict_type = 4;
        
    
    /*Conflict solving*/
    switch (conflict_type) {
        case 0:                                                         /*No conflict occured, animal just moved*/              
            next_world[i2][j2].animal = tmp.animal;                     /*Updates the new cell with the moving animal*/
            next_world[i2][j2].breed_age = tmp.breed_age;
            next_world[i2][j2].starve_age = tmp.starve_age;
            break;

        case 1:                                                         /*Conflict between two foxes*/
            if(tmp.breed_age > next_world[i2][j2].breed_age) {          /*The fox with larger breeding age wins*/
                next_world[i2][j2].breed_age = tmp.breed_age;       
                
                if(next_world[i2][j2].starve_age != 0) {                /*Checks if one of the foxes has eatten a rabbit*/
                    next_world[i2][j2].starve_age = tmp.starve_age;
                }
            }
            else if(tmp.breed_age == next_world[i2][j2].breed_age)      /*If the foxes have the same breeding age*/
            {
                if(tmp.starve_age < next_world[i2][j2].starve_age)      /*The one with least starvation wins*/
                {
                    next_world[i2][j2].starve_age = tmp.starve_age;
                }
            }
            break;

        case 2:                                                         /*Conflict between two rabbits*/
            if(tmp.breed_age > next_world[i2][j2].breed_age) {          /*The one with higher beeding age wins*/
                next_world[i2][j2].breed_age = tmp.breed_age;
            }
            break;

        case 3:                                                         /*Conflict between a fox and a rabbit, rabbits bumps into fox*/
            next_world[i2][j2].animal = FOX;                            /*The fox always wins*/
            next_world[i2][j2].starve_age = 0;                          /*Resetting the fox's starvation*/
            break;
            
        case 4:                                                         /*Conflict between a fox and a rabbit, fox bumps into rabbit*/
            next_world[i2][j2].animal = FOX;                            /*The fox always wins*/
            next_world[i2][j2].breed_age = tmp.breed_age;               /*Saving the fox's breeding age*/
            next_world[i2][j2].starve_age = 0;                          /*Resetting the fox's starvation*/
            break;
            
        default:
            break;

    }
    updated[i2][j2] = 1;
}


/* FUNCTION: update_next_world
 * ARGUMENTS: Cell **curr_world     - Current world state
              Cell **next_world     - Next world state
              int **updated         - Matrix marking already updated cells
              int i                 - Current cell's row number
              int j                 - Current cell's column number
              int direction         - Direction the animal will move
              int breed_age_f       - Fox breeding age
              int breed_age_r       - Rabbit breeding age

 * DESCRIPTION: Moves the animal to a new cell, depending on the direction computed, solving any
                conflicts that may occur and breeding computation
 * RETURN: N/A
 */
void update_next_world(Cell **curr_world, Cell **next_world, int **updated, int i, int j, int direction, int breed_age_f, int breed_age_r){
    
    if(direction == dir_UP)                                                                                         /*If movement direction is upwards*/
    {
        solve_conflicts_and_breeding(curr_world, next_world, updated, i, j, i-1, j, breed_age_f, breed_age_r);      /*Call solving breeding and conflicts funtion for destination [i-1][j]*/
    }
    else if(direction == dir_DOWN)                                                                                  /*If movement direction is downwards*/
    {
        solve_conflicts_and_breeding(curr_world, next_world, updated, i, j, i+1, j, breed_age_f, breed_age_r);      /*Call solving breeding and conflicts funtion for destination [i+1][j]*/
    }
    else if(direction == dir_LEFT)                                                                                  /*If movement direction is to the left*/
    {
        solve_conflicts_and_breeding(curr_world, next_world, updated, i, j, i, j-1, breed_age_f, breed_age_r);      /*Call solving breeding and conflicts funtion for destination [i][j-1]*/
    }   
    else if(direction == dir_RIGHT)                                                                                 /*If movement direction is to the right*/
    { 
        solve_conflicts_and_breeding(curr_world, next_world, updated, i, j, i, j+1, breed_age_f, breed_age_r);      /*Call solving breeding and conflicts funtion for destination [i][j+1]*/
    }  
    else                                                                                                            /*If no movement occurs*/
    {           
        next_world[i][j].animal = curr_world[i][j].animal;                                                          /*Maintain animal in same position, but write it to "new_world"*/
        next_world[i][j].breed_age = curr_world[i][j].breed_age;
        next_world[i][j].starve_age = curr_world[i][j].starve_age;
    }
}


/* FUNCTION: movement
 * ARGUMENTS: Cell **curr_world     - Current world state
              Cell **next_world     - Next world state
              int **updated         - Matrix marking already updated cells
              int M                 - Number of rows
              int N                 - Number of columns
              int i                 - Cell's current row
              int j                 - Cell's current column
              int *empty_cells      - Vector marking adjacent empty cells
              int *adj_rabbit_cells - Vector marking adjacent rabbit cells
              int breed_age_f       - Breeding age for foxes
              int breed_age_r       - Breeding age for rabbits

 * DESCRIPTION: Calculates the overall movement of a cell, including breeding and conflict solving, writing the result to the "next_world" matrix
 * RETURN: N/A
 */
void movement(Cell **curr_world, Cell **next_world, int **updated, int M, int N, int i, int j, int *empty_cells, int *adj_rabbit_cells, int breed_age_f, int breed_age_r){
    
    int c = i*N + j;                    /*Cell number*/
    int p = 0;                          /*Number of adjacent empty cells*/
    int adj_rabbit = 0;                 /*Number of adjacent rabbit cells*/
    int direction = -1;                 /*Movement direction*/
    
    empty_cells[dir_UP] = 0;            /*Initialize empty cells direction marking vector*/
    empty_cells[dir_RIGHT] = 0;
    empty_cells[dir_DOWN] = 0;
    empty_cells[dir_LEFT] = 0;

    adj_rabbit_cells[dir_UP] = 0;       /*Initialize adjacent rabbit cells direction marking vector*/
    adj_rabbit_cells[dir_RIGHT] = 0;
    adj_rabbit_cells[dir_DOWN] = 0;
    adj_rabbit_cells[dir_LEFT] = 0;

    if(curr_world[i][j].animal == RABBIT)                                               /*If the current cell is a rabbit*/
        get_moves_rabbit(curr_world,i,j,M,N,&p,&empty_cells[0]);                        /*Calculate possible moves for rabbit (how many empty cells and in which directions)*/

    if(curr_world[i][j].animal == FOX)                                                             /*If the current cell is a fox*/
        get_moves_fox(curr_world,i,j,M,N,&p,&adj_rabbit,&empty_cells[0],&adj_rabbit_cells[0]);     /*Calculate possible moves for fox (how many adjacent rabbits or empty cells and in which directions)*/


    if(adj_rabbit > 0)                  /*If the number of adjacent rabbits is greater than 0 (only applies for foxes, as "adj_rabbits" is unchanged for rabbits)*/                                  
    {
        direction = calculate_direction(adj_rabbit,adj_rabbit_cells,c);                                      /*Calculate movement direction*/
        update_next_world(curr_world, next_world, updated, i, j, direction, breed_age_f, breed_age_r);       /*Update the next world matrix for conflict resolution in next step*/
    } 
    else if(p >= 0)                     /*If the number of empty cells is greater than 0*/
    {
        direction = calculate_direction(p,empty_cells,c);                                                     /*Calculate movement direction*/
        update_next_world(curr_world, next_world, updated, i, j, direction, breed_age_f, breed_age_r);        /*Update the next world matrix for conflict resolution in next step*/   
    }               
}



/* FUNCTION: compute_generation
 * ARGUMENTS: Cell **curr_world     - Current world state
              Cell **next_world     - Next world state
              int *empty_cells      - Vector marking adjacent empty cells
              int *adj_rabbit_cells - Vector marking adjacent rabbit cells
              int **updated         - Matrix marking already updated cells
              int M                 - Number of rows
              int N                 - Number of columns
              int breed_age_f       - Breeding age for foxes
              int breed_age_r       - Breeding age for rabbits
              int starve_age_f      - Starving age for foxes

 * DESCRIPTION: Computes a full generation, divided in two sub-generations
 * RETURN: N/A
 */
void compute_generation(Cell **curr_world, Cell **next_world, int *empty_cells, int *adj_rabbit_cells, int **updated, int M, int N, \
                        int breed_age_f, int breed_age_r, int starve_age_f) {

    int i, j;                   /*Iteration variables*/
    int loop_start;             /*Variable to determine which column to start on when traversing a row*/

    /*Red Sub-Generation*/
    loop_start = 0;
    for(i = 0; i < M; i++){                                                                       /*Iterate over all rows*/                
        for(j = loop_start; j < N; j+=2){                                                         /*But only the red cells*/
            if(curr_world[i][j].animal == FOX || curr_world[i][j].animal == RABBIT)               /*Only try to move if the cell has a fox or rabbit*/
                movement(curr_world, next_world, updated, M, N, i, j, empty_cells, adj_rabbit_cells, breed_age_f, breed_age_r);
        }
        loop_start == 0? (loop_start = 1): (loop_start = 0);                                      /*Switch starting column for next row*/                         
    }

    /*Copy world*/
    for(i = 0; i < M; i++) {                                                                        
         for(j = 0; j < N; j++){
            curr_world[i][j].animal = next_world[i][j].animal;
            curr_world[i][j].breed_age = next_world[i][j].breed_age;
            curr_world[i][j].starve_age = next_world[i][j].starve_age;
        }
    }

    /*Black Sub-Generation*/
    loop_start = 1; 
    for(i = 0; i < M; i++) {                                                                      /*Iterate over all rows*/   
        for(j = loop_start; j < N; j+=2){                                                         /*But only the black cells*/
            if(updated[i][j] == 1) {       /*If an animal has moved to a black cell from a red cell, then update is set to 1 in that position, meaning in this loop it is skipped and resets it*/
                updated[i][j] = 0;
                continue;
            }
            if(curr_world[i][j].animal == FOX || curr_world[i][j].animal == RABBIT)               /*Only try to move if the cell has a fox or rabbit*/
                movement(curr_world, next_world, updated, M, N, i, j, empty_cells, adj_rabbit_cells, breed_age_f, breed_age_r);
        }
        loop_start == 0? (loop_start = 1) : (loop_start = 0);                                     /*Switch starting column for next row*/
    }
    
    /*Copy world and remove foxes with starving ages greater than the limit*/
    for(i = 0; i < M; i++) {                
        for(j = 0; j < N; j++) {
            if(next_world[i][j].animal == FOX) {
                if(next_world[i][j].starve_age >= starve_age_f) {
                    next_world[i][j].animal = EMPTY;
                    next_world[i][j].breed_age = -1;
                    next_world[i][j].starve_age = -1;
                }
            }
            curr_world[i][j].animal = next_world[i][j].animal;
            curr_world[i][j].breed_age = next_world[i][j].breed_age;
            curr_world[i][j].starve_age = next_world[i][j].starve_age; 
        }
    }
}

/* FUNCTION: print_results
 * ARGUMENTS: Cell **world - World to count rabbits, foxes and rocks from
              int M        - Number of rows
              int N        - Number of columns

 * DESCRIPTION: Prints the final results to stdout
 * RETURN: N/A
 */
void print_results(Cell **world, int M, int N) {
    
    int i, j;                               /*Iteration variables*/
    int final_rabbits = 0;                  /*Final number of rabbits*/
    int final_foxes = 0;                    /*Final number of foxes*/
    int final_rocks = 0;                    /*Final number of rocks*/
    
    for(i = 0; i < M; i++){                 /*Count foxes, rabbits and rocks*/
        for(j = 0; j < N; j++){
            if(world[i][j].animal == FOX) {
                final_foxes++;
            }
            else if(world[i][j].animal == RABBIT) {
                final_rabbits++;
            }
            else if(world[i][j].animal == ROCK) {
                final_rocks++;
            }
        }
    }
    
    fprintf(stdout, "%d %d %d\n",final_rocks, final_rabbits, final_foxes);   /*Print result to stdout*/
}


int main(int argc, char *argv[]) {

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
    uint32_t seed = 0;                          /*Seed for randomization*/
    int i, j;                                   /*Iteration variables*/
    
    int *empty_cells;                           /*Vector for marking adjacent empty cells*/
    int *adj_rabbit_cells;                      /*Vector for marking adjacent rabbit cells*/
    int **updated;                              /*Matrix for marking already updated cells in a sub-generation*/
    Cell **curr_world;                          /*Matrix for storing current world state*/
    Cell **next_world;                          /*Matrix for storing next world state*/
    
    if (argc != 11) {                           /*Verifies if number of args is correct*/
        exit(-1);
    }

    for(int i = 1; i <= 10; i++) {               /*Verifies if all args are positive integers and assigns them to the correct variables*/
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

    curr_world = (Cell **) malloc(M * sizeof(Cell));
    if(curr_world == NULL) {
        exit(-1);
    }

    next_world = (Cell **) malloc(M * sizeof(Cell));
    if(next_world == NULL) {
        exit(-1);
    }
    
    empty_cells = (int*) malloc(4*sizeof(int));
    if(empty_cells == NULL) {
        exit(-1);
    }

    adj_rabbit_cells = (int*) malloc(4*sizeof(int));
    if(adj_rabbit_cells == NULL) {
        exit(-1);
    }

    initialize_worlds(curr_world, next_world, M, N);    /*Initialize worlds*/

    
    updated = (int **) calloc(M, sizeof(int *));        /*Create 'updated' matrix to mark already moved animal*/
    if(updated == NULL) {
        exit(-1);
    }

    for(int i = 0; i < M; i++) {                    
        updated[i] = (int *) calloc(N, sizeof(int));
        if (updated[i] == NULL) {
            exit(-1);
        }
    }
    
    /*Generate elements to fill world*/
    generate_element(curr_world, num_rocks, ROCK, &seed, M, N);
    generate_element(curr_world, num_rabbits, RABBIT, &seed, M, N);
    generate_element(curr_world, num_foxes, FOX, &seed, M, N);


    for(i = 0; i < M; i++) {                            /*Copy current world to next world matrix*/
        for(j = 0; j < N; j++) {
            next_world[i][j].animal = curr_world[i][j].animal;
            next_world[i][j].breed_age = curr_world[i][j].breed_age;
            next_world[i][j].starve_age = curr_world[i][j].starve_age;
        }
    }

    exec_time = -omp_get_wtime();                       /*Start timing simulation*/
    
    /*Generation loop*/
    for(int i = 1; i <= generations; i++) {
        compute_generation(curr_world, next_world, empty_cells, adj_rabbit_cells, updated, M, N, breed_age_f, breed_age_r, starve_age_f); 
    }

    exec_time += omp_get_wtime();                       /*Finish timing simulation*/
    fprintf(stderr, "%.5fs\n", exec_time);              /*Print execution time to stderr*/

    print_results(curr_world, M, N);                    /*Print results to stdout*/

    /*Free allocated memory*/
    for(i = 0; i < M; i++) {        
        free(curr_world[i]);
        free(next_world[i]);
        free(updated[i]);
    }
    free(curr_world);
    free(next_world);
    free(updated);
    free(empty_cells);
    free(adj_rabbit_cells);

    return 0;
}