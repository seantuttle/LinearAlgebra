// main.c
#include "vector.h"
#include "matrix.h"
#include "equations.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define NUM_ROWS 4
#define NUM_COLS 4

int main(void)
{
    clock_t tic = clock();
    
    printf("Use this file to test out features of the Linear Algebra library!\n");


    clock_t toc = clock();
    double total_time = (double) (toc - tic) / CLOCKS_PER_SEC;
    printf("\nTime taken: %f seconds\n", total_time);

    return 0;
}
