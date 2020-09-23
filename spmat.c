#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "spmat.h"

/*
 * spmat Summary:
 *
 * A linked-list implementation of sparse matrices.

 * A linked-lists sparse matrix contains a set of linked lists, one for each row of the matrix.
 * Each cell contains its value and column index, as well as
 * a pointer to the next non-zero element in the same row (or NULL if it is the last non-zero cell in its row).
 *
 * spmat_allocate_list   -  Creates a new spars-matrix of size 'n'
 * add_row_list          -  Adds 'row' as the 'i'th row in the matrix
 * free_list             -  Frees all resources used by 'matrixA'
 * mult_list             -  Multiplies matrix 'matrixA' by vector 'v'
 * copy_spmat            -  Copies to 'ABg' the rows and columns relevant to the vertexes in 'g', from 'AB'
 * Aij                   -  Returns matrixA[i][j]
 * multRowAVec           -  Returns the product of row 'numRow' from 'matrixA' with vector 's'
 * sumRowsA              -  Updates vector 'result' to hold the sums of each row in matrix 'matrixA'
 * ARowI                 -  Updates vector 'row' to hold only the columns from 'g' according to 'matrixA'
 * deletefromCol         -  Updates 'sum' as it would have been if vertex 'i' was not in 'matrixA'
 * normOfBabs            -  Returns the norm of matrix "B" which consists of 'matrixA' and 'vecK'+'M' (see matrixB.c)
 */

spmat* spmat_allocate_list(int n);
void add_row_list(struct _spmat *matrixA, int *row, int i, int k);
void free_list(struct _spmat *matrixA);
void mult_list(const struct _spmat *matrixA, const double *v, double *result);
void copy_spmat(struct _spmat *ABg, struct _spmat *AB, struct _listVertex *g);
int Aij(struct _spmat *matrixA,int i, int j);
double multRowAVec(struct _spmat *matrixA, int numRow, double *s);
void sumRowsA(struct _spmat *matrixA, double *result);
int ARowI (struct _spmat *matrixA, int *row, int indexRow, listVertex *g);
void deletefromCol(struct _spmat* matrixA,int sign, double* sum,int i);
double normOfBabs(struct _spmat* matrixA, int size, int M, int* vecK);

typedef struct _list {
    int value;
    int col;
    struct _list *next;
} list;

/* Creates a new spars-matrix of size 'n' */
spmat* spmat_allocate_list(int n) {
    list *arr = (list *) calloc(n, sizeof(list));
    spmat *spmatList = (spmat *) calloc(1, sizeof(spmat));
    if (arr == NULL) {
        printf("Error: allocation failed\n");
        exit(EXIT_FAILURE);
    }
    if (spmatList == NULL) {
        printf("Error: allocation failed\n");
        exit(EXIT_FAILURE);
    }

    arr->col = -5;
    spmatList->n = n;
    spmatList->free = free_list;
    spmatList->mult = mult_list;
    spmatList->add_row = add_row_list;
    spmatList->Aij = Aij;
    spmatList->copy_spmat = copy_spmat;
    spmatList->multRowAVec = multRowAVec;
    spmatList->sumRowsA = sumRowsA;
    spmatList->private = arr;
    spmatList->ARowI = ARowI;
    spmatList->deletefromCol = deletefromCol;
    spmatList->normOfBabs = normOfBabs;

    return spmatList;
}

/*
* The function receives a pointer to a spars matrix 'matrixA', a pointer to a row 'row',
* and the row index 'i'.
*
* The function adds 'row' as the 'i'th row in the matrix. Called before any other call,
* exactly n times in order (i = 0 to n-1)
*/
void add_row_list(struct _spmat *matrixA, int *row, int i, int k) {
    int j = 0;
    list *parr, *tail = NULL, *arr = (list *) matrixA->private;

    parr = &arr[i];
    parr->col = -5;
    while (j < k) {
        if (tail != NULL) {
            tail->next = (list *) calloc(1, sizeof(list));
            if (tail->next == NULL) {
                printf("Error: allocation failed\n");
                exit(EXIT_FAILURE);
            }
            tail = tail->next;
            tail->value = 1;
            tail->col = *row;
        } else {/*if it's the first node in the row*/
            parr->value = 1;
            parr->col = *row;
            parr->next = NULL;
            tail = parr;
        }
        row++;
        j++;
    }
    if (k != 0) {
        tail->next = NULL;
    }
}

/* Frees all resources used by 'matrixA' */
void free_list(struct _spmat *matrixA) {
    int i = 0, n = matrixA->n;
    list *node, *temp, *arr = (list *) matrixA->private;

    while (i < n) {/*iterate each row in the matrix*/
        node = arr;
        if (node->value != 0) {/*the row has node to free*/
            node = node->next;
            while (node != NULL) {/*frees each cell in the row*/
                temp = node;
                node = node->next;
                free(temp);
            }
        }
        arr++;
        i++;
    }
    free(matrixA->private);
    free(matrixA);
}

/*
* The function receives a pointer to a spars matrix 'matrixA', a pointer to a start vector 'v',
* and a pointer to the final vector 'result'.
* The function multiplies matrix 'matrixA' by vector 'v', and saves the product to 'result'
* ('result' is pre-allocated).
*/
void mult_list(const struct _spmat *matrixA, const double *v, double *result) {
    int i = 0, n = matrixA->n;
    double sum, *presult;
    list *node, *arr = (list *) matrixA->private;

    presult = result;/*pointer for result vector*/
    while (i < n) {/*iterate each row in the matrix*/
        node = arr;
        sum = 0;
        if ((node->value != 0) && (node->col != -5)) {/*the row has at least one node*/
            while (node != NULL) {/*dot product row with given vector*/
                sum = sum + (double) (v[node->col]);
                node = node->next;
            }
        }
        *presult = sum;
        presult++;
        arr++;
        i++;
    }
}

/*
 * The function receives a large matrix 'AB' and copies to the smaller matrix 'ABg'
 * the rows and columns relevant to the vertexes in 'g'.
 */
void copy_spmat(struct _spmat *ABg, struct _spmat *AB, struct _listVertex *g) {
    int size, i = 0, numNeighbors, *row, *prow;
    struct _listVertex *pointer1;
    size = g->size(g);

    row = (int *) calloc(size, sizeof(int));/*array to save one row at a time*/
    if (row == NULL) {
        printf("Error: allocation failed\n");
        exit(EXIT_FAILURE);
    }
    pointer1 = g;
    while (pointer1 != NULL) {
        prow = row;
        numNeighbors = ARowI(AB, prow, pointer1->value, g);
        prow = row;
        ABg->add_row(ABg, prow, i, numNeighbors);/*add row to spars-matrix according to the implementation*/
        pointer1 = pointer1->next;
        i++;
    }
    free(row);
}

/* returns matrixA[i][j] */
int Aij(struct _spmat *matrixA,int i, int j) {
    list *parr, *arr = (list *) matrixA->private;
    if (i == j) {
        return 0;
    }
    parr = &arr[i];
    while ((parr != NULL) && (parr->col <= j)) {
        if (parr->col == j) {
            return 1;
        }
        parr = parr->next;
    }
    return 0;
}

/* Returns the product of row 'numRow' from 'matrixA' with vector 's' */
double multRowAVec(struct _spmat *matrixA, int numRow, double *s) {
    list *parr, *arr = (list *) matrixA->private;
    int sum = 0;
    double *ps;
    ps = s;
    parr = &arr[numRow];
    if ((parr->value != 0) && (parr->col != -5)) {/*the row has at least one node*/
        while ((parr != NULL)) {
            if (parr->col == numRow) {
                sum += ps[parr->col] * (-1);
            } else {
                sum += ps[parr->col];
            }
            parr = parr->next;
        }
    }
    return sum;
}

/* Updates vector 'result' to hold the sums of each row in matrix 'matrixA' */
void sumRowsA(struct _spmat *matrixA, double *result) {
    int i = 0, n = matrixA->n;
    double sum, *presult;
    list *node, *arr = (list *) matrixA->private;

    presult = result;/*pointer for result vector*/
    while (i < n) {/*iterate each row in the matrix*/
        node = arr;
        sum = 0;
        if ((node->value != 0) && (node->col != -5)) {/*the row has at least one node*/
            while (node != NULL) {/*dot product row with given vector*/
                sum++;
                node = node->next;
            }
        }
        *presult = sum;
        presult++;
        arr++;
        i++;
    }
}


/*
 * The function receives a matrix 'matrixA', an empty vector 'row' and a list of vertexes 'g'.
 * The function updates 'row' according to 'g' and the big matrix 'matrixA' so that we can use 'row' later on
 * to build a smaller spars matrix only with the vertexes in 'g'.
 */
int ARowI (struct _spmat *matrixA, int *row, int indexRow, listVertex *g) {
    list *arr = (list *) matrixA->private;
    list *node = &arr[indexRow];
    listVertex *Pg = g;
    int *Prow = row, numNeighbors = 0, j = 0;
    if ((node->value != 0) && (node->col != -5)) {/*the row has at least one node*/
        while ((node != NULL) && (Pg != NULL)) {
            if (Pg->value != indexRow) {
                if (node->col == Pg->value) {
                    /* we insert to 'row'*/
                    *Prow = j;
                    numNeighbors++;
                    Prow++;
                    node = node->next;
                    Pg = Pg->next;
                    j++;
                } else {
                    if (node->col > Pg->value) {
                        Pg = Pg->next;
                        j++;
                    } else {
                        node = node->next;
                    }
                }
            } else {
                Pg = Pg->next;
                j++;
            }
        }
    }
    return numNeighbors;
}

/*
 * The function receives a matrix 'matrixA', and a vector of sums 'sum'.
 * The function updates 'sum' as it would have been if vertex 'i' was not in 'matrixA'.
*/
void deletefromCol(struct _spmat* matrixA,int sign, double* sum,int i) {
    list *parr, *arr = (list *) matrixA->private;
    double *psum;
    parr = &arr[i];
    psum = sum;
    if ((parr->value != 0) && (parr->col != -5)) {/*the row has at least one node*/
        while ((parr != NULL)) {
            psum[parr->col] += 2 * sign;
            parr = parr->next;
        }
    }
}

/*
 * The function receives 'matrixA' and vector 'vecK' + 'M', which are the parts of a "B" matrix. (see matrixB.c)
 * The function returns the largest sum of all the row sums of B, if each cell in B was
 * taken as an absolute value.
 */
double normOfBabs(struct _spmat* matrixA, int size, int M, int* vecK) {
    list *node, *arr = (list *) matrixA->private;
    int i = 0, *ki, *kj, j;
    double max, temp;

    max = 0;
    ki = vecK;

    while (i < size) {/*iterate each row in the matrix*/
        node = arr;
        temp = 0;
        j = 0;
        kj = vecK;
        if ((node->value != 0) && (node->col != -5)) {/*the row has at least one node*/
            while (node != NULL) {/*dot product row with given vector*/
                while (j < node->col) {
                    temp += fabs((double) (((double) (*ki * *kj)) / M));
                    j++;
                    kj++;
                }
                temp += fabs((double) 1 - ((double) (((double) (*ki * *kj)) / M)));
                node = node->next;
                j++;
                kj++;
            }
        }
        while (j < size) {
            temp += fabs((double) (((double) (*ki * *kj)) / M));
            j++;
            kj++;
        }
        if (temp > max) {
            max = temp;
        }
        ki++;
        arr++;
        i++;
    }
    return max;
}

