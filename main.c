#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "listGroup.h"
#include "spmat.h"
#include "matrixB.h"

int adjacencyMatrix(FILE* input, spmat* matrixA, int size, int* K);
void S(double *eigenVec, double *vectorS, int size);
double modularityQ(double* vectorS, matrixB* B, int size, int M, int* vecK, listVertex* g);
void dotProductdeltafVec(double* deltaf, double* vec, double* result, int size);
void dotProductBvec(matrixB* B, double* vec, double* result, int size, int M, int* vecK, listVertex* g);
void dotProductMatAndVecM(double* PKs, double* vec, int* K, int M, int size, listVertex* g);
double dotProduct2Vec(double* v1, double* v2, int size);
double dotProductVecAndK(double* v1, int* v2, int size, listVertex* g);
void dominantAbsoluteEigenVector(matrixB* B, double *currV, double *nextV, int size, double norm, int M, int* vecK, listVertex* g, int sizeB);
double dominantAbsoluteEigenValue(matrixB* B, double *bk, int size, double norm, int M,int* vecK, listVertex* g);
double leadingEigenValue(matrixB* B, double *eigenVec, int size, double norm, int M, int* vecK, listVertex* g);
void algorithm4(double *s, matrixB* B, int size, listVertex* g, int M, int *vecK, int sizeB);
double BgRow(matrixB* B, int iA, int size, double *s, int M, int* vecK, int iValue, listVertex* g);
void Algorithm3(listGroup *P, int numV, matrixB* B, double norm, int M, FILE *output, int* vecK);
listGroup* insertingGroups(listVertex *g, listGroup *G);
void divideG(double* s, listVertex* g, matrixB *Bg, double norm, int M, int size, int* vecK, int sizeB);
matrixB* BG(listVertex* g, matrixB *B, int M, int size,int* vecK);
void oneGroup(double* s, int size);
void outputFile(FILE *output, listGroup *O, int sizeO);
void initialize(listVertex* vertexes, int* indexValue);
double denomDivision(double *nextV, double *currV, int size);
double multRowK(int *vecK, int numRow, double *vecS,int M, int size, int iValue, listVertex* g);
double sumK(int* k, int size, listVertex* g);
void sumRowsKM(double* f, int* K, int M, int size, listVertex* g);
void sumRowsB(matrixB* B, int size, int M, int* vecK, listVertex* g);
void Bsums(matrixB* B, listVertex *g, double *sums, double *s, int M, int* vecK, int size, int *indexValue);
void sumAfterOneChange(int M, int sign, double *sums, int j,matrixB *B, listVertex *g, int* vecK, int maxval, int* indexValue);
double sImprovment(int size,int *indices, int indexImproveMax,double *s, double improveMax);
void maximization(matrixB* B, int *indices,double *sums, listVertex *g, double *s,int *vecK,int *indexValue, int M, int size, listVertex *unmoved);
void divideS(int size, double *s,listVertex *g1,listVertex *g2,listVertex *g, int *sizeG1, int *sizeG2);
int conditionCheck(int sizeG1,int sizeG2,listGroup **O,listGroup **P,listVertex *g2,listVertex *g1,int sizeO);

int main(int argc, char* argv[]) {
    FILE *input, *output;
    int numOfVertex, i, *pointerK, M, *vecK;
    double norm;
    spmat *pointerA;
    matrixB *B;
    listGroup *P;
    listVertex *g;

    /*checking that we got 3 arguments or more*/
    if (argc < 3) {
        printf("Error: invalid input\n");
        exit(EXIT_FAILURE);
    }
    /*opening input and output*/
    input = fopen(argv[1], "rb");
    if (input == NULL) {
        printf("Error: invalid input\n");
        exit(EXIT_FAILURE);
    }
    output = fopen(argv[2], "wb");
    if (output == NULL) {
        printf("Error: invalid output\n");
        exit(EXIT_FAILURE);
    }

    i = fread(&numOfVertex, sizeof(int), 1, input);
    if (i != 1) {
        printf("Error: couldn't read input file\n");
        exit(EXIT_FAILURE);
    }
    B = allocateB(numOfVertex);
    vecK = (int *) calloc(numOfVertex, sizeof(int));
    if (vecK == NULL) {
        printf("Error: allocation failed\n");
        exit(EXIT_FAILURE);
    }
    pointerA = B->matrixA;
    pointerK = vecK;
    M = adjacencyMatrix(input, pointerA, numOfVertex, pointerK);
    if (M >= -0.00001 && M <= 0.00001) {
        printf("Error: M=0\n");
        exit(EXIT_FAILURE);
    }
    pointerK = vecK;
    norm = B->matrixA->normOfBabs(B->matrixA, numOfVertex, M, pointerK);

    i = 0;
    g = allocateV();
    while (i < numOfVertex) {/*inserting all vertexes to a single group 'P' (a trivial division)*/
        g->addVertex(g, i);
        i++;
    }
    P = allocateG(g);
    pointerK = vecK;
    Algorithm3(P, numOfVertex, B, norm, M, output, pointerK);

    B->freeB(B);
    free(vecK);
    fclose(input);
    fclose(output);
    return 0;
}

/*
 * The function recieves a sub-group of vertexes 'g', the original matrix 'B'
 * and the sum of ranks 'M'.
 *
 * The function calculates and return sub-matrix of 'B' according to 'g'.
 */
matrixB* BG(listVertex* g, matrixB *B, int M, int size, int* vecK) {
    matrixB *Bg;
    listVertex *pg;
    int *PvecK;
    size = g->size(g);
    Bg = allocateB(size);
    pg = g;
    Bg->matrixA->copy_spmat(Bg->matrixA, B->matrixA, pg);
    PvecK = vecK;
    sumRowsB(Bg, size, M, PvecK, pg);

    return Bg;
}

/*
 * The function receives a a matrix 'Bg' according to vertexes in 'g',
 * and updates vector 's' to be the division to 2 groups according to the eigenvector.
 */
void divideG(double* s, listVertex* g, matrixB *Bg, double norm, int M, int size, int* vecK, int sizeB) {
    double *b0, *Pb0, *eigenVec, eigenVal, Q, *Ps;
    int n, i, *pk;
    listVertex *Pg;
    size = Bg->size;
    b0 = (double *) calloc(size, sizeof(double));
    if (b0 == NULL) {
        printf("Error: allocation failed\n");
        exit(EXIT_FAILURE);
    }
    n = 0;
    while (n == 0) {
        i = 0;
        Pb0 = b0;
        while (i < size) {
            *Pb0 = rand();
            if (*Pb0 != 0) {
                n = 1;
            }
            Pb0++;
            i++;
        }
    }
    eigenVec = (double *) calloc(size, sizeof(double));
    if (eigenVec == NULL) {
        printf("Error: allocation failed\n");
        exit(EXIT_FAILURE);
    }
    Pg = g;
    pk = vecK;
    dominantAbsoluteEigenVector(Bg, b0, eigenVec, size, norm, M, pk, Pg, sizeB);
    Pg = g;
    pk = vecK;
    eigenVal = leadingEigenValue(Bg, eigenVec, size, norm, M, pk, Pg);
    Ps = s;
    S(eigenVec, Ps, size);
    Ps = s;
    Pg = g;
    pk = vecK;
    algorithm4(Ps, Bg, size, Pg, M, pk, sizeB);
    if (eigenVal <= -0.00001) {/*The group g is indivisible*/
        Ps = s;
        oneGroup(s, size);
    } else {
        Pg = g;
        pk = vecK;
        Q = modularityQ(s, Bg, size, M, pk, Pg);
        if (Q <= -0.00001) {/*The group g is indivisible*/
            Ps = s;
            oneGroup(s, size);
        }
    }
    free(b0);
    free(eigenVec);
}

/*
 * The function receives vector 's' and updates it to be a vector of 1
 * (means the vertexes are all in the same group)
 */
void oneGroup(double* s, int size) {
    int i;
    i = 0;
    while (i < size) {
        *s++ = 1;
        i++;
    }
}

/*
 * The function receives a pointer to the input file, pointer to an empty matrix,
 * the size of a single column and row in the matrix, and a pointer to an array of the vertex's ranks.
 *
 * The function updates the adjacency matrix 'matrixA' and returns the sum of ranks.
 */
int adjacencyMatrix(FILE* input, spmat* matrixA, int size, int* K) {
    int j = 0, i, sum = 0, *row, *prow, numNeighbors;

    row = (int *) calloc(size, sizeof(int));/*array to save one row at a time*/
    if (row == NULL) {
        printf("Error: allocation failed\n");
        exit(EXIT_FAILURE);
    }
    while (j < size) {
        i = fread(&numNeighbors, sizeof(int), 1, input);
        if (i != 1) {
            printf("Error: couldn't read input file\n");
            exit(EXIT_FAILURE);
        }
        *K = numNeighbors; /* i'th v-rank is numNeighbors */
        sum += numNeighbors; /* counting sum of ranks */
        prow = row;
        i = fread(prow, sizeof(int), numNeighbors, input);
        if (i != numNeighbors) {
            printf("Error: couldn't read input file\n");
            exit(EXIT_FAILURE);
        }
        prow = row;/*prow is a pointer to row*/
        matrixA->add_row(matrixA, prow, j, numNeighbors);/*add row to spars-matrix according to the implementation*/
        K++;
        j++;
    }
    free(row);
    return sum;
}

/*
 * The function receives the eigenvector and updates vector 's' to hold the division to 2 groups.
 * if: eigenvector[i] < 0 -> s[i] = -1
 * else: s[i] = 1
 */
void S(double *eigenVec, double *vectorS, int size) {
    int i = 0;
    while (i < size) {
        if (*eigenVec < -0.00001) {
            *vectorS = -1;
        } else {
            *vectorS = 1;
        }
        vectorS++;
        eigenVec++;
        i++;
    }
}

/*
 * The function receives matrix 'B' and vector 'vectorS'
 * calculates and returns 'deltaQ'
 */
double modularityQ(double* vectorS, matrixB* B, int size, int M, int* vecK, listVertex* g) {
    double Q, *result, *Presult, *Ps;
    int *Pk;
    listVertex *Pg;

    result = (double *) calloc(size, sizeof(double));
    if (result == NULL) {
        printf("Error: allocation failed\n");
        exit(EXIT_FAILURE);
    }
    Presult = result;
    Ps = vectorS;
    Pk = vecK;
    Pg = g;
    dotProductBvec(B, Ps, Presult, size, M, Pk, Pg);/*Bs*/
    Ps = vectorS;
    Presult = result;
    dotProductdeltafVec(B->f, Ps, Presult, size);/*Bhat*s = Bs-delta*f*s*/
    Q = (0.5) *
        dotProduct2Vec(vectorS, result, size); /* 0.5(s*Bhat*s) =  0.5(s(B*s-delta*f*s)) = 0.5(s*B*s-s*delta*f*s) */
    free(result);
    return Q;
}

/*
 * The function receives 'deltaf' which holds the (sum of B's rows)*I, and its size 'size'.
 * The function updates 'result' to hold the dot product of vector 'vec' and 'deltaf'.
 */
void dotProductdeltafVec(double* deltaf, double* vec, double* result, int size) {
    double *f, *Presult, *Pvec;
    int i = 0;
    f = deltaf;
    Presult = result;
    Pvec = vec;
    while (i < size) {
        *Presult++ -= *f++ * *Pvec++;
        i++;
    }
}

/*
 * The function updates vector 'result' to hold the dot product of matrix 'B' and vector 'vec'
 */
void dotProductBvec(matrixB* B, double* vec, double* result, int size, int M, int* vecK, listVertex* g) {
    double *Pvec, *Presult;
    int *PK;
    listVertex *Pg;
    Pvec = vec;
    Presult = result;
    B->matrixA->mult(B->matrixA, Pvec, Presult);/* dot product spars mat A and vector s */

    Presult = result;
    Pvec = vec;
    PK = vecK;
    Pg = g;
    dotProductMatAndVecM(Presult, Pvec, PK, M, size, Pg);/* dot product mat K/M and vector s */
    Presult = result;
}

/*
 * The function updates vector 'result' to hold the dot product of matrix 'K/M' with 'vec'
 */
void dotProductMatAndVecM(double* result, double* vec, int* K, int M, int size, listVertex* g) {
    int i = 0, *Pk;
    double KS, temp, *pvec, *presult;
    listVertex *Pg;
    Pk = K;
    pvec = vec;
    Pg = g;
    KS = dotProductVecAndK(pvec, Pk, size, Pg);
    temp = (double) KS / M;
    Pk = K;
    Pg = g;
    presult = result;
    while (i < size) {
        *presult -= Pk[Pg->value] * temp;
        Pg = Pg->next;
        presult++;
        i++;
    }
}

/*
 * Returns the dot product of double vectors 'v1' and 'v2'
 */
double dotProduct2Vec(double* v1, double* v2, int size) {
    double *pv1, *pv2, dot = 0;
    pv1 = v1;
    pv2 = v2;
    while (size > 0) {
        dot += *pv1++ * *pv2++;
        size--;
    }

    return dot;
}

/*
 * Returns the dot product of double vector 'v1' and vector 'K' (ranks)
 */
double dotProductVecAndK(double* v1, int* k, int size, listVertex* g) {
    double *pv1, dot = 0;
    int *pk;
    pv1 = v1;
    pk = k;
    while (size > 0) {
        dot += *pv1++ * pk[g->value];
        g = g->next;
        size--;
    }
    return dot;
}

/*
* The function receives the  matrix 'B', a pointer to the first vector 'currV',
* and size of matrix 'size'.
*
* In each iteration the function uses current vector to produce a new one.
* When done, the vector produced in the final iteration is the desired eigen vector.
* The function updates the eigen vector to vector 'currV'.
* If there is an infinite loop while finding the eigen vector, the function stops and exits.
*/
void dominantAbsoluteEigenVector(matrixB* B, double *currV, double *nextV, int size, double norm, int M, int* vecK, listVertex* g, int sizeB) {
    int i, k = 0, flag, *pk;
    double diff, *nextP, *currP, *infinite, *infP, *f;
    listVertex *Pg;
    infinite = (double *) calloc(size, sizeof(double));
    if (infinite == NULL) {
        printf("Error: allocation failed\n");
        exit(EXIT_FAILURE);
    }
    diff = 1;
    while (diff > 0.00001) {
        nextP = nextV;
        currP = currV;
        pk = vecK;
        Pg = g;
        dotProductBvec(B, currP, nextP, size, M, pk, Pg);/* nextP  = Bg*currP */
        nextP = nextV;
        currP = currV;
        f = B->f;
        dotProductdeltafVec(f, currP, nextP,
                            size);/* nextP = Bhat*currP = Bg*currP - delta*f*currP = currP(Bg - delta*f)   */

        i = 0;
        nextP = nextV;
        currP = currV;
        while (i < size) {/*shift Bghat*s = Bg*s + norm*I*s */
            *nextP = *nextP + norm * *currP;
            i++;
            nextP++;
            currP++;
        }
        nextP = nextV;
        currP = currV;
        diff = denomDivision(nextP, currP, size);
        /*currV = nextV*/
        i = 0;
        nextP = nextV;
        currP = currV;
        infP = infinite;
        flag = 0;
        while (i < size) {
            if ((k % 2 == 0) && (k >= size * 2)) {
                *infP = *currP;
            }
            if ((k % 2 == 1) && (k >= size * 2)) {
                if (*infP != *nextP) {
                    flag = 1;
                }
            }
            *currP++ = *nextP++;
            infP++;
            i++;
        }
        if ((k % 2 == 1) && (k >= size * 2)) {
            if (flag == 0) {
                printf("Error: infinite loop in power iteration\n");
                exit(EXIT_FAILURE);
            }
        }
        if (k > (12 * sizeB + 10000 * sizeB)) {
            printf("Error: infinite loop in power iteration\n");
            exit(EXIT_FAILURE);
        }
        k++;
    }
    free(infinite);
}

/*
 * The function receives a vector 'nextV' and 'currV',
 * calculates the denominator (as explained in the assignment)
 * and divides 'nextV' accordingly.
 *
 * The function returns the difference between the current and previous eigenvector.
 */
double denomDivision(double *nextV, double *currV, int size) {
    int i = 0;
    double *nextP = nextV, diff, *currP, denom = 0;
    while (i < size) {
        denom += *nextP * *nextP;
        i++;
        nextP++;
    }
    denom = (double) sqrt(denom);
    diff = 0;
    nextP = nextV;
    currP = currV;
    i = 0;
    while (i < size) {
        *nextP = (double) *nextP / denom;
        if (fabs(*currP - *nextP) > diff) {
            diff = fabs(*currP - *nextP);
        }
        i++;
        nextP++;
        currP++;
    }
    return diff;
}

/*
 * The function receives matrix 'B', the eigenvector 'bk',
 * and returns the matching eigenvalue of the shifted matrix 'B' hat.
 */
double dominantAbsoluteEigenValue(matrixB* B, double *bk, int size, double norm, int M,int* vecK, listVertex* g) {
    double *Pbk, *result, *Presult, up, down, *f;
    int i, *pk;
    listVertex *Pg;
    result = (double *) calloc(size, sizeof(double));
    if (result == NULL) {
        printf("Error: allocation failed\n");
        exit(EXIT_FAILURE);
    }
    Pbk = bk;
    Presult = result;
    Pg = g;
    pk = vecK;
    dotProductBvec(B, Pbk, Presult, size, M, pk, Pg);/* result  = Bg*bk */
    Pbk = bk;
    Presult = result;
    f = B->f;
    dotProductdeltafVec(f, Pbk, Presult,
                        size);/* result = BhatShift*bk = BgShift*bk - delta*f*bk = bk(BgShift - delta*f) */

    i = 0;
    Pbk = bk;
    Presult = result;
    while (i < size) {/*shift Bghat*s = Bg*s + norm*I*s */
        *Presult = *Presult + norm * *Pbk;
        i++;
        Presult++;
        Pbk++;
    }

    Presult = result;
    Pbk = bk;
    up = dotProduct2Vec(Pbk, Presult, size);/* bk*(BhatShift*bk) = bk*result */
    Pbk = bk;
    down = dotProduct2Vec(Pbk, Pbk, size); /* bk*bk */

    free(result);
    return (double) up / down;
}

/*
 * The function receives an eigenvector 'eigenVec' and returns the leading eigenvalue of matrix 'B'
 */
double leadingEigenValue(matrixB* B, double *eigenVec, int size, double norm, int M, int* vecK, listVertex* g) {
    double beta;
    int *pk;
    listVertex *pg;
    pk = vecK;
    pg = g;
    beta = dominantAbsoluteEigenValue(B, eigenVec, size, norm, M, pk, pg);
    beta = beta - norm;
    return beta;
}

/*
 * The function vector 's' which represents the division into 2 groups,
 * and uses algorithm 4 (as described in the assignment) to maximize the division.
 */
void algorithm4(double *s, matrixB* B, int size, listVertex* g, int M, int *vecK, int sizeB) {
    int *D, *indices, *Pindices, *indexValue, *PindexValue;
    double *Ps, *sums, *psums;
    listVertex *k, *unmoved, *pg;
    D = vecK;
    indices = (int *) calloc(size, sizeof(int));
    indexValue = (int *) calloc(sizeB, sizeof(int));
    sums = (double *) calloc(size, sizeof(double));
    if ((indices == NULL) || (indexValue == NULL)) {
        printf("Error: allocation failed\n");
        exit(EXIT_FAILURE);
    }

    unmoved = allocateV();
    unmoved->copy(g, unmoved, size);
    PindexValue = indexValue;
    k = unmoved;
    initialize(k, PindexValue);
    Pindices = indices;
    psums = sums;
    pg = g;
    Ps = s;
    PindexValue = indexValue;
    k = unmoved;
    maximization(B, Pindices, psums, pg, Ps, D, PindexValue, M, size, k);

    free(indices);
    free(indexValue);
    free(sums);
}

/*
 * The function maximize a division into 2 groups,
 * moving one node that brings the largest improvement to the other group.
 */
void maximization(matrixB* B, int *indices,double *sums, listVertex *g, double *s,int *vecK,int *indexValue, int M, int size, listVertex *unmoved) {
    int *Pindices, *D, *PindexValue, i, indexScoreMax, maxval, id, indexImproveMax;
    double Q = 1, improveMax, *psums, *Ps, scoreMax, checkScore, improvePrev;
    listVertex *pg, *k, *maxNode = NULL;
    while (Q > 0.00001) {
        Pindices = indices;
        improveMax = -1 * HUGE_VAL;
        psums = sums;
        pg = g;
        Ps = s;
        D = vecK;
        PindexValue = indexValue;
        Bsums(B, pg, psums, Ps, M, D, size, PindexValue);/* sums = B*s */
        PindexValue = indexValue;
        i = 0;
        while (i < size) {
            k = unmoved;
            scoreMax = -1 * HUGE_VAL;
            psums = sums;
            if (i != 0) {
                psums = sums;
                D = vecK;
                PindexValue = indexValue;
                sumAfterOneChange(M, s[indexScoreMax], psums, indexScoreMax, B, k, D, maxval, PindexValue);
            }
            k = unmoved;
            while (k != NULL) {
                Ps = s;
                psums = sums;
                id = PindexValue[k->value];
                D = vecK;
                checkScore = (-4) * Ps[id] * psums[id] + 4 * (double) ((D[k->value]) * (D[k->value])) / M;
                if (checkScore > scoreMax) {
                    scoreMax = checkScore;
                    indexScoreMax = id;
                    maxNode = k;
                    maxval = maxNode->value;
                }
                k = k->next;
            }
            s[indexScoreMax] = -1 * s[indexScoreMax];
            *Pindices = indexScoreMax;;
            if (i == 0) {
                improveMax = scoreMax;
                improvePrev = scoreMax;
                indexImproveMax = i;
            } else {
                improvePrev = improvePrev + scoreMax;
                if (improvePrev > improveMax) {
                    improveMax = improvePrev;
                    indexImproveMax = i;
                }
            }
            k = unmoved;
            unmoved = k->removeVertex(k, maxNode);
            Pindices++;
            i++;
        }
        Pindices = indices;
        Q = sImprovment(size, Pindices, indexImproveMax, s, improveMax);
        if (Q > 0.00001) {
            unmoved = allocateV();
            unmoved->copy(g, unmoved, size);
        }
    }
}

/*
 * The function finds the node that will bring the largest improvement in the modularity if we mve it, and moves it.
 * It returns the improvement 'deltaQ'
 */
double sImprovment(int size,int *indices, int indexImproveMax,double *s, double improveMax) {
    double Q, *ps;
    int i, j, *Pindices;
    i = size - 1;
    Pindices = indices;
    ps = s;
    while (i != indexImproveMax) {
        j = Pindices[i];
        ps[j] = -1 * ps[j];
        i--;
    }
    if (indexImproveMax == (size - 1)) {
        Q = 0;
    } else {
        Q = improveMax;
    }
    return Q;
}

/*
 * The function receives a groups of vertexes 'vertexes' and an array 'indexValue'
 * and updates 'indexValue' to hold the index of each vertex according to sub group 'vertexes'
 */
void initialize(listVertex* vertexes, int* indexValue) {
    int i = 0;
    listVertex *g;
    g = vertexes;
    while (g != NULL) {
        indexValue[g->value] = i;
        g = g->next;
        i++;
    }
}

/*
 * The function receives vector 'sums' that represent dot the product of matrix 'B' with vector s.
 * It updates vector 'sum' to be dot product of 'B' and vector s if
 * s[j] changes to be -s[j].
 */
void sumAfterOneChange(int M, int sign, double *sums, int j, matrixB *B, listVertex *g, int* vecK, int maxval, int* indexValue) {
    double *psum;
    int *k, kj, ki, *PindexValue, id;
    listVertex *pg;
    psum = sums;
    B->matrixA->deletefromCol(B->matrixA, sign, psum, j);
    psum = sums;
    pg = g;
    k = vecK;
    kj = k[maxval];
    PindexValue = indexValue;
    while (pg != NULL) {
        id = PindexValue[pg->value];
        ki = k[pg->value];
        psum[id] -= (double) 2 * (double) sign * (double) (((double) (ki * kj)) / M);
        pg = pg->next;
    }
}

/*
 * The function updates vector 'sums' to be the dot product of 'B' and vector 's'
 */
void Bsums(matrixB* B, listVertex *g, double *sums, double *s, int M, int* vecK, int size, int *indexValue) {
    int i = 0, id, *PindexValue;
    double *Ps, *sum;
    int *D;
    listVertex *pg, *pg2;
    sum = sums;
    pg = g;
    Ps = s;
    D = vecK;
    PindexValue = indexValue;
    while (i < size) {
        id = PindexValue[pg->value];
        pg2 = g;
        *sum = BgRow(B, id, size, Ps, M, D, pg->value, pg2);
        pg = pg->next;
        sum++;
        i++;
    }
}

/*
 * The function returns the dot product of the 'numRow'th row in 'B' with vector 's'
 */
double BgRow(matrixB* B, int numRow, int size, double *s, int M, int* vecK, int iValue, listVertex* g) {
    int Ais, *pk;
    double Kis, *Ps;
    listVertex *pg;
    Ps = s;
    Ais = B->matrixA->multRowAVec(B->matrixA, numRow, Ps); /* row i in matrix A * vec s when s[i] = -s[i] */
    pk = vecK;
    Ps = s;
    pg = g;
    Kis = multRowK(pk, numRow, Ps, M, size, iValue, pg);
    return (double) (Ais - Kis);
}

/*
 * The function returns the dot product of the 'numRow'th row in matrix K/M (ranks) with vector 's'
 */
double multRowK(int *vecK, int numRow, double *vecS,int M, int size, int iValue, listVertex* g) {
    int i = 0, *pk, ki, kj;
    double sum = 0, *ps;
    listVertex *pg;

    pk = vecK;
    ki = pk[iValue];
    ps = vecS;
    pg = g;
    while (i < size) {
        if (i == numRow) {
            sum += (double) -1 * (double) *ps * (double) (((double) (ki * ki)) / M);
        } else {
            kj = pk[pg->value];
            sum += (double) *ps * (double) (((double) (ki * kj)) / M);
        }
        pg = pg->next;
        ps++;
        i++;
    }
    return sum;
}

/*
 * The function updates vector f of the original matrix 'B' to hold the sum of its rows.
 */
void sumRowsB(matrixB* B, int size, int M, int* vecK, listVertex* g) {
    int *pk;
    double *f;
    listVertex *pg;
    f = B->f;
    B->matrixA->sumRowsA(B->matrixA, f);
    pk = vecK;
    pg = g;
    f = B->f;
    sumRowsKM(f, pk, M, size, pg);
}

/*
 * The function updates vector f of the original matrix 'B' to hold the sum of its rows
 * using the fact that B = A - K/M, it calculates K/M part
 */
void sumRowsKM(double* f, int* K, int M, int size, listVertex* g) {
    int i = 0, *Pk;
    double KS, temp, *vecf;
    listVertex *Pg;
    Pg = g;
    Pk = K;
    KS = sumK(Pk, size, Pg);
    temp = (double) KS / M;
    Pk = K;
    Pg = g;
    vecf = f;
    while (i < size) {
        *vecf -= Pk[Pg->value] * temp;
        Pg = Pg->next;
        vecf++;
        i++;
    }
}

/*
 * The function returns the sum of ranks according to sub group of vertexes 'g'
 */
double sumK(int* k, int size, listVertex* g) {
    double dot = 0;
    int *pk;
    pk = k;
    while (size > 0) {
        dot += pk[g->value];
        g = g->next;
        size--;
    }
    return dot;
}

/*
 * The function receives a pointer to a set of groups of vertexes 'P',
 * that holds one group with all the vertexes in the matrix.
 * Also, it receives a pointer to an empty set of groups of vertexes 'O'.
 * 'numV' is the number of vertexes.
 *
 * The function finishes when 'O' holds the best division of the vertexes to groups.
 */
void Algorithm3(listGroup *P, int numV, matrixB* B, double norm, int M, FILE *output, int* vecK) {
    listVertex *g1, *g2, *g, *Pg, *pg1, *pg2;
    listGroup *pointerO, *O;
    int sizeG1, sizeG2, i = 0, size, sizeO = 0, *pVecK; /* 'sizeO' is the number of groups in 'O' */
    double *s, *ps;
    matrixB *Bg;
    O = allocateG(NULL); /*output set*/

    while (P != NULL) {
        g1 = allocateV();
        g2 = allocateV();
        g = P->value;
        size = g->size(g);
        s = (double *) calloc(size, sizeof(double));
        if (s == NULL) {
            printf("Error: allocation failed\n");
            exit(EXIT_FAILURE);
        }
        if (size == numV) {
            pVecK = vecK;
            Pg = g;
            sumRowsB(B, size, M, pVecK, Pg);
            pVecK = vecK;
            divideG(s, g, B, norm, M, size, pVecK, numV);
        } else {
            pVecK = vecK;
            Bg = BG(g, B, M, size, pVecK);
            pVecK = vecK;
            divideG(s, g, Bg, norm, M, size, pVecK, numV);
        }

        sizeG1 = 0;
        sizeG2 = 0;
        i = 0;
        Pg = g;
        ps = s;
        pg1 = g1;
        pg2 = g2;
        divideS(size, ps, pg1, pg2, Pg, &sizeG1, &sizeG2);

        /* finished with vector s */
        free(s);
        P = P->removeG(P);
        pg1 = g1;
        pg2 = g2;
        sizeO = conditionCheck(sizeG1, sizeG2, &O, &P, pg2, pg1, sizeO);
        if (size != numV) {
            Bg->freeB(Bg);
        }
    }
    i = fwrite(&sizeO, sizeof(int), 1, output);/*1 int of size of group O for output*/
    if (i != 1) {
        printf("Error: Failed writing to output \n");
        exit(EXIT_FAILURE);
    }
    pointerO = O;
    outputFile(output, pointerO, sizeO);
    O->freeGroup(O);
}

/*
 * The function receives vector 's' that holds a division into 2 groups,
 * and adds each vertex to the correct sub group according to 's'.
 */
void divideS(int size, double *s,listVertex *g1,listVertex *g2,listVertex *g, int *sizeG1, int *sizeG2) {
    int i = 0;
    listVertex *Pg;
    Pg = g;

    /* initializing the 2 groups according to s */
    while (i < size) {
        if (s[i] == 1) {
            *sizeG1 += 1;
            g1->addVertex(g1, Pg->value);
        } else {
            if (s[i] == -1) {
                *sizeG2 += 1;
                g2->addVertex(g2, Pg->value);
            } else {
                printf("Error: invalid value in vector S\n");
                exit(EXIT_FAILURE);
            }
        }
        i++;
        Pg = Pg->next;
    }
}

/*
 * The function receives 2 groups and puts them back to group O or P
 * according to the conditions in the project.
 */
int conditionCheck(int sizeG1,int sizeG2,listGroup **O,listGroup **P,listVertex *g2,listVertex *g1,int sizeO) {
    /* CONDITION 1 */
    if (sizeG1 == 0) {
        /* inserting 'G2' to 'O' */
        *O = insertingGroups(g2, *O);
        sizeO++;
        g1->freeVertex(g1);

    } else if (sizeG2 == 0) {
        /* inserting 'G1' to 'O' */
        *O = insertingGroups(g1, *O);
        sizeO++;
        g2->freeVertex(g2);
    } else {
        /* CONDITION 2 */
        if (sizeG1 == 1) {
            /* inserting 'G1' to 'O' */
            *O = insertingGroups(g1, *O);
            sizeO++;
        }

        if (sizeG2 == 1) {
            /* inserting 'G2' to 'O' */
            *O = insertingGroups(g2, *O);
            sizeO++;
        }

        /* CONDITION 3 */
        if (sizeG1 > 1) {
            /* inserting 'G1' to 'P' */
            *P = insertingGroups(g1, *P);
        }

        if (sizeG2 > 1) {
            /* inserting 'G2' to 'P' */
            *P = insertingGroups(g2, *P);
        }
    }
    return sizeO;
}

/*
 * Inserting group of vertexes with the vertexes 'g' to 'G'
 */
listGroup* insertingGroups(listVertex *g, listGroup *G) {
    listGroup *TEMP, *pointerG;
    pointerG = G;
    if (pointerG == NULL) {
        G = allocateG(g);
    } else {
        if (pointerG->value == NULL) {
            pointerG->value = g;
        } else {
            TEMP = allocateG(g);
            TEMP->next = pointerG;
            G = TEMP;
        }
    }
    return G;
}

/*
 * The function writes to the output file (according to the format in the project),
 * the best division to groups of all the vertexes.
 */
void outputFile(FILE *output, listGroup *O, int sizeO) {
    int i = 0, j = 0, sizeGroup = 0, val, k;
    listVertex *g, *pointerg;
    while (i < sizeO) {
        g = O->value;
        pointerg = g;
        sizeGroup = pointerg->size(pointerg);
        k = fwrite(&sizeGroup, sizeof(int), 1, output);/*1 int of size of group i in O for output*/
        if (k != 1) {
            printf("Error: write to output failed\n");
            exit(EXIT_FAILURE);
        }
        j = 0;
        while (j < sizeGroup) {
            val = g->value;
            k = fwrite(&val, sizeof(int), 1, output);/*1 int of size of group O for output*/
            if (k != 1) {
                printf("Error: write to output failed\n");
                exit(EXIT_FAILURE);
            }
            g = g->next;
            j++;
        }
        O = O->next;
        i++;
    }
}