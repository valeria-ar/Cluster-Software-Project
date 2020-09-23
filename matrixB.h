#ifndef MATRIXB_H_
#define MATRIXB_H_
#include "spmat.h"

typedef struct _matrixB {
    struct _spmat *matrixA;
    double *f;
    int size;

    /* Frees all resources used by 'matrixB' */
    void (*freeB)(struct _matrixB *B);

} matrixB;

struct _matrixB* allocateB(int size);

#endif /* MATRIXB_H_ */
