#ifndef LISTVERTEX_H_
#define LISTVERTEX_H_

typedef struct _listVertex {
    int value;
    int head;
    struct _listVertex *next;
    struct _listVertex *prev;

    /* Frees all resources used by listV */
    void (*freeVertex)(struct _listVertex *listV);

    /* Creates a new listVertex node with the value 'val' and inserts it before 'g' in the list */
    void (*addVertex)(struct _listVertex *g, int val);

    /* Removes node 'g' from the list */
    struct _listVertex *(*removeVertex)(struct _listVertex *k, struct _listVertex *g);

    /* Returns the length of the list of which 'g' is the first node */
    int (*size)(struct _listVertex *g);

    /* Copies nodes from one list to another */
    void (*copy)(struct _listVertex *from, struct _listVertex *to, int size);

} listVertex;

/* Creates a new listVertex node */
struct _listVertex* allocateV();

#endif /* LISTVERTEX_H_ */
