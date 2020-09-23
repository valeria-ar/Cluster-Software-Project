#include <stdio.h>
#include <stdlib.h>
#include "listVertex.h"

/*
 * listVertex Summary:
 *
 * A listVertex is a linked-list node, which value is a vertex, as well as the pointer to the next node
 * in the list (might be NULL), and a pointer to the previous node.
 *
 * AllocateV    -    Creates a new listVertex node
 * freeVertex   -    Frees all resources used by 'listV'
 * addVertex    -    Creates a new listVertex node with the value 'val' and inserts it before 'g' in the list
 * removeVertex -    Removes node 'g' from the list
 * size         -    Returns the length of the list of which 'g' is the first node
 * copy         -    Copies nodes from one list to another
 */

listVertex* allocateV();
void addVertex(listVertex* g,int value);
void freeVertex(struct _listVertex *listV);
listVertex* removeVertex(listVertex* k, listVertex* g);
int size(listVertex* g);
void copy(listVertex* from, listVertex* to, int size);

listVertex* allocateV() {
    listVertex *listV = (listVertex *) calloc(1, sizeof(listVertex));

    if (listV == NULL) {
        printf("Error: couldn't allocate\n");
        exit(EXIT_FAILURE);
    }

    listV->next = NULL;
    listV->prev = NULL;
    listV->head = -5;

    listV->addVertex = addVertex;
    listV->copy = copy;
    listV->freeVertex = freeVertex;
    listV->removeVertex = removeVertex;
    listV->size = size;

    return listV;
}

/* Frees all resources used by listV */
void freeVertex(struct _listVertex *listV) {
    listVertex *temp, *node;
    temp = listV;
    while (temp != NULL) {
        node = temp;
        temp = temp->next;
        free(node);
    }
}

/*
 * Adds a new vertex with value 'val' before 'g'.
 * After the function, 'g' is the new vertex
 */
void addVertex(listVertex* g, int val) {
    listVertex *temp;
    if ((g->head == -5) || (g == NULL)) {
        g->head = val;
        g->value = val;
        g->prev = g;
    } else {
        temp = allocateV();
        temp->value = val;
        temp->prev = g->prev;
        g->prev->next = temp;
        g->prev = temp;
    }
}

/* Removes vertex 'g' from the list. 'k' is the first node in the list to which 'g' belongs */
listVertex* removeVertex(listVertex* k, listVertex* g) {
    listVertex *temp;
    if ((g->prev == g)) {/*g is 1 node, after removing it g will be null*/
        temp = g;
        g = NULL;
        k = g;
        free(temp);
    } else {/*there is more then 1 node*/
        if (g->value == g->head) {/*the node is the head, we will update the next node to be the head*/
            g->next->head = g->next->value;
            temp = g;
            g = g->next;
            g->prev = temp->prev;
            k = g;
        } else if (g->next == NULL) {/*the node is the last one*/
            temp = g;
            g = g->prev;
            k->prev = g;
            g->next = NULL;
        } else {/*the node is in the middle*/
            temp = g;
            g = g->prev;
            g->next = temp->next;
            (temp->next)->prev = g;
        }
        free(temp);
    }
    return k;
}

/*
 * The function receives a node 'g' and returns the number of nodes in the list.
 * 'g' is assumed to be the first node in the list.
 */
int size(listVertex* g) {
    int size = 0;
    if ((g == NULL) || (g->head == -5)) {
        return size;
    }
    while (g != NULL) {
        size++;
        g = g->next;
    }
    return size;
}

/* Copies 'size' nodes from list 'from' to list 'to' */
void copy(listVertex* from, listVertex* to, int size) {
    int i = 0;
    listVertex *Pfrom, *Pto;
    Pfrom = from;
    Pto = to;
    to->value = -5;
    while (i < size) {
        Pto->addVertex(Pto, Pfrom->value);
        Pfrom = Pfrom->next;
        i++;
    }
}
