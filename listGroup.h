#ifndef LISTGROUP_H_
#define LISTGROUP_H_
#include "listVertex.h"

typedef struct _listGroup {
    struct _listVertex *value;
    struct _listGroup *next;

    /* Frees all resources used by 'listG' */
    void (*freeGroup)(struct _listGroup *listG);

    /* Returns 1 if 'listG' holds no vertexes, and 0 otherwise */
    int (*isEmpty)(struct _listGroup *listG);

    /* Deletes 'listG' and returns a pointer to the next node */
    struct _listGroup *(*removeG)(struct _listGroup *listG);

} listGroup;

/* Creates a new listGroup node, which value is the list of vertexes 'value' */
listGroup* allocateG(struct _listVertex *value);

#endif /* LISTGROUP_H_ */
