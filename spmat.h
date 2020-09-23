#ifndef _SPMAT_H
#define _SPMAT_H
#include "listVertex.h"

typedef struct _spmat {
    /* Matrix size (n*n) */
    int n;

    /* Adds 'row' as the 'i'th row in the matrix */
    void (*add_row)(struct _spmat *matrixA, int *row, int i, int k);

    /* Frees all resources used by 'matrixA' */
    void (*free)(struct _spmat *matrixA);

    /* Multiplies matrix 'matrixA' by vector 'v' */
    void (*mult)(const struct _spmat *matrixA, const double *v, double *result);

    /* Copies to 'ABg' the rows and columns relevant to the vertexes in 'g', from 'AB' */
    void (*copy_spmat)(struct _spmat *ABg, struct _spmat *AB, struct _listVertex *g);

    /* Returns matrixA[i][j] */
    int (*Aij)(struct _spmat *matrixA, int i, int j);

    /* Returns the product of row 'numRow' from 'matrixA' with vector 's' */
    double (*multRowAVec)(struct _spmat *matrixA, int numRow, double *s);

    /* Updates vector 'result' to hold the sums of each row in matrix 'matrixA' */
    void (*sumRowsA)(struct _spmat *matrixA, double *result);

    /* Updates vector 'row' to hold only the columns from 'g' according to 'matrixA' */
    int (*ARowI)(struct _spmat *matrixA, int *row, int indexRow, listVertex *g);

    /* Updates 'sum' as it would have been if vertex 'i' was not in 'matrixA' */
    void (*deletefromCol)(struct _spmat *matrixA, int sign, double *sum, int i);

    /* Returns the norm of matrix "B" which consists of 'matrixA' and 'vecK'+'M' (see matrixB.c) */
    double (*normOfBabs)(struct _spmat *matrixA, int size, int M, int *vecK);

    /* Private field for inner implementation.
     * Should not be read or modified externally */
    void *private;
} spmat;

/* Creates a new spars-matrix of size 'n' */
struct _spmat* spmat_allocate_list(int n);

#endif
