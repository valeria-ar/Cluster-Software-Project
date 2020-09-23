#include "matrixB.h"
#include <stdio.h>
#include <stdlib.h>

/*
 * matrixB Summary:
 *
 * A matrixB is a matrix which is represented by 'matrixA' (sparse matrix),
 * vector 'f' which is the max sum of matrix B row (we need it to compute B hat),
 * and vector 'K' (ranks vector) - we save only one vector K for the original B and then use it for sub matrix B.
 *
 * allocateB  -      Creates a new matrixB, meaning, creates new sparse matrix with 'size' according to number of
 * 			    	 vertexes in this sub matrixB. The function allocates as well vector f.
 * freeB      -      Frees all resources used by 'matrixB'
 */

void freeB(struct _matrixB *B);
matrixB* allocateB(int size);

matrixB* allocateB(int size) {
    matrixB *B = (matrixB *) calloc(1, sizeof(matrixB));
    double *f = (double *) calloc(size, sizeof(double));
    spmat *matrixA;
    matrixA = spmat_allocate_list(size);
    if (B == NULL) {
        printf("error-couldn't allocate\n");
        exit(EXIT_FAILURE);
    }
    if (f == NULL) {
        printf("error-couldn't allocate\n");
        exit(EXIT_FAILURE);
    }

    B->size = size;
    B->matrixA = matrixA;
    B->f = f;
    B->freeB = freeB;
    return B;
}

/* Frees all resources used by 'B' */
void freeB(struct _matrixB *B) {
    B->matrixA->free(B->matrixA);
    free(B->f);

    free(B);
}





