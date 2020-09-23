#include "listGroup.h"
#include <stdio.h>
#include <stdlib.h>

/*
 * listGroup Summary:
 *
 * A listGroup is a linked-list node, which value is a group of vertexes (i.e listVertex), as well as
 * the pointer to the next node in the list (might be NULL).
 *
 * allocateG    -    Creates a new listGroup node, which value is the list of vertexes 'value'
 * freeGroup    -    Frees all resources used by 'listG'
 * isEmpty      -    Returns 1 if 'listG' holds no vertexes, and 0 otherwise
 * removeG      -    Deletes 'listG' and returns a pointer to the next node
 */

listGroup* allocateG(listVertex *value);
void freeGroup(struct _listGroup *listG);
listGroup* removeG(struct _listGroup *listG);
int isEmpty(struct _listGroup *listG);

listGroup* allocateG(listVertex *value) {
    listGroup *listG = (listGroup *) calloc(1, sizeof(listGroup));

    if (listG == NULL) {
        printf("error-couldn't allocate\n");
        exit(EXIT_FAILURE);
    }
    listG->value = value;
    listG->next = NULL;

    listG->freeGroup = freeGroup;
    listG->isEmpty = isEmpty;
    listG->removeG = removeG;

    return listG;

}

/* Frees all resources used by 'listG' */
void freeGroup(struct _listGroup *listG) {
    listGroup *temp;
    listVertex *node;
    if (listG != NULL) {
        if (listG->value == NULL) {
            free(listG);
        } else {
            while (listG != NULL && listG->value != NULL) {
                temp = listG;
                listG = listG->next;
                node = temp->value;
                (node)->freeVertex(node);
                free(temp);
            }
        }
    }
}

/* Returns 1 if 'listG' holds no vertexes, and 0 otherwise */
int isEmpty(struct _listGroup *listG) {
    if (listG->value == NULL) {
        return 1;
    } else
        return 0;
}

/*
 * The function receives a pointer to 'listG', and removes it from the list.
 * 'listG' is assumed to be the first node in the list.
 *
 * The function returns a pointer to the next node after 'listG'. (might be NULL)
 */
listGroup* removeG(struct _listGroup *listG) {
    listGroup *temp = listG,
            *nextG = listG->next;
    (temp->value)->freeVertex(temp->value);
    free(temp);
    listG = nextG;
    return listG;
}