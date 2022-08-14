/*
 * Copyright 2010 Volkan Yazıcı <volkan.yazici@gmail.com>
 * Copyright 2006-2010 The Apache Software Foundation
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not
 * use this file except in compliance with the License. You may obtain a copy of
 * the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
 * License for the specific language governing permissions and limitations under
 * the License.
*/
#include "../../config.h"
#include "../../globals.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include "ads/pqueue.h"


#define left(i)   ((i) << 1)
#define right(i)  (((i) << 1) + 1)
#define parent(i) ((i) >> 1)


pqueue_t *
pqueue_init(size_t n,
            pqueue_cmp_pri_f cmppri,
            pqueue_get_pri_f getpri,
            pqueue_set_pri_f setpri,
            pqueue_get_pos_f getpos,
            pqueue_set_pos_f setpos)
{
    pqueue_t *q;

    if (!(q = malloc(sizeof(pqueue_t))))
        return NULL;

    /* Need to allocate n+1 elements since element 0 isn't used. */
    if (!(q->d = malloc((n + 1) * sizeof(void *)))) {
        free(q);
        return NULL;
    }

    q->size = 1;
    q->avail = q->step = (n+1);  /* see comment above about n+1 */
    q->cmppri = cmppri;
    q->setpri = setpri;
    q->getpri = getpri;
    q->getpos = getpos;
    q->setpos = setpos;

    return q;
}


void
pqueue_free(pqueue_t *q)
{
    free(q->d);
    free(q);
}


size_t
pqueue_size(pqueue_t *q)
{
    /* queue element 0 exists but doesn't count since it isn't used. */
    return (q->size - 1);
}


static void
bubble_up(pqueue_t *q, size_t i)
{
    size_t parent_node;
    void *moving_node = q->d[i];
    pqueue_pri_t moving_pri = q->getpri(moving_node);

    for (parent_node = parent(i);
         ((i > 1) && q->cmppri(q->getpri(q->d[parent_node]), moving_pri));
         i = parent_node, parent_node = parent(i))
    {
        q->d[i] = q->d[parent_node];
        q->setpos(q->d[i], i);
    }

    q->d[i] = moving_node;
    q->setpos(moving_node, i);
}


static size_t
maxchild(pqueue_t *q, size_t i)
{
    size_t child_node = left(i);

    if (child_node >= q->size)
        return 0;

    if ((child_node+1) < q->size &&
        q->cmppri(q->getpri(q->d[child_node]), q->getpri(q->d[child_node+1])))
        child_node++; /* use right child instead of left */

    return child_node;
}


static void
percolate_down(pqueue_t *q, size_t i)
{
    size_t child_node;
    void *moving_node = q->d[i];
    pqueue_pri_t moving_pri = q->getpri(moving_node);

    while ((child_node = maxchild(q, i)) &&
           q->cmppri(moving_pri, q->getpri(q->d[child_node])))
    {
        q->d[i] = q->d[child_node];
        q->setpos(q->d[i], i);
        i = child_node;
    }

    q->d[i] = moving_node;
    q->setpos(moving_node, i);
}


int
pqueue_insert(pqueue_t *q, void *d)
{
    void *tmp;
    size_t i;
    size_t newsize;

    if (!q) return 1;

    /* allocate more memory if necessary */
    if (q->size >= q->avail) {
        newsize = q->size + q->step;
        if (!(tmp = realloc(q->d, sizeof(void *) * newsize)))
            return 1;
        q->d = tmp;
        q->avail = newsize;
    }

    /* insert item */
    i = q->size++;
    q->d[i] = d;
    bubble_up(q, i);

    return 0;
}


void
pqueue_change_priority(pqueue_t *q,
                       pqueue_pri_t new_pri,
                       void *d)
{
    size_t posn;
    pqueue_pri_t old_pri = q->getpri(d);

    q->setpri(d, new_pri);
    posn = q->getpos(d);
    if (q->cmppri(old_pri, new_pri))
        bubble_up(q, posn);
    else
        percolate_down(q, posn);
}


int
pqueue_remove(pqueue_t *q, void *d)
{
    size_t posn = q->getpos(d);
    q->d[posn] = q->d[--q->size];
    if (q->cmppri(q->getpri(d), q->getpri(q->d[posn])))
        bubble_up(q, posn);
    else
        percolate_down(q, posn);

    return 0;
}

int
pqueue_remove_n(pqueue_t *q, int n)
{
    //size_t posn = q->getpos(d);
    void *d=q->d[n];
    q->d[n] = q->d[--q->size];
    if (q->cmppri(q->getpri(d), q->getpri(q->d[n])))
        bubble_up(q, n);
    else
        percolate_down(q, n);

    return 0;
}

void *
pqueue_pop(pqueue_t *q)
{
    void *head;

    if (!q || q->size == 1)
        return NULL;

    head = q->d[1];
    q->d[1] = q->d[--q->size];
    percolate_down(q, 1);

    return head;
}

void *
pqueue_pop_n(pqueue_t *q, int n)
{
    void *head;

    if (!q || q->size == n)
        return NULL;

    head = q->d[n];
    q->d[n] = q->d[--q->size];
    percolate_down(q, n);


    return head;
}
void *
pqueue_peek(pqueue_t *q)
{
    void *d;
    if (!q || q->size == 1)
        return NULL;
    d = q->d[1];
    return d;
}

void *
pqueue_peek_n(pqueue_t *q, int d)
{
    //size_t posn = q->getpos(d);
    void *ds;
    if (!q || q->size == d)
        return NULL;
    ds = q->d[d];
    return ds;
}

void
pqueue_dump(pqueue_t *q,
            FILE *out,
            pqueue_print_entry_f print)
{
    int i;

    fprintf(stdout,"posn\tleft\tright\tparent\tmaxchild\t...\n");
    for (i = 1; i < q->size ;i++) {
        fprintf(stdout,
                "%d\t%d\t%d\t%d\t%ul\t",
                i,
                left(i), right(i), parent(i),
                (unsigned int)maxchild(q, i));
        print(out, q->d[i]);
    }
}


static void
set_pos(void *d, size_t val)
{
    /* do nothing */
}


static void
set_pri(void *d, pqueue_pri_t pri)
{
    /* do nothing */
}


void
pqueue_print(pqueue_t *q,
             FILE *out,
             pqueue_print_entry_f print)
{
    pqueue_t *dup;
	void *e;

    dup = pqueue_init(q->size,
                      q->cmppri, q->getpri, set_pri,
                      q->getpos, set_pos);
    dup->size = q->size;
    dup->avail = q->avail;
    dup->step = q->step;

    memcpy(dup->d, q->d, (q->size * sizeof(void *)));

    while ((e = pqueue_pop(dup)))
		print(out, e);

    pqueue_free(dup);
}


static int
subtree_is_valid(pqueue_t *q, int pos)
{
    if (left(pos) < q->size) {
        /* has a left child */
        if (q->cmppri(q->getpri(q->d[pos]), q->getpri(q->d[left(pos)])))
            return 0;
        if (!subtree_is_valid(q, left(pos)))
            return 0;
    }
    if (right(pos) < q->size) {
        /* has a right child */
        if (q->cmppri(q->getpri(q->d[pos]), q->getpri(q->d[right(pos)])))
            return 0;
        if (!subtree_is_valid(q, right(pos)))
            return 0;
    }
    return 1;
}


int
pqueue_is_valid(pqueue_t *q)
{
    return subtree_is_valid(q, 1);
}







 



 

/* Function to create an empty priority queue */
pqueue_bsf* pqueue_bsf_init(int k)
{

        pqueue_bsf *q;

    if (!(q = malloc(sizeof(pqueue_bsf))))
        return NULL;
    if (!(q->position = malloc(sizeof(long unsigned int)*k)))
        return NULL;
    if (!(q->knn = malloc(sizeof(float)*k)))
        return NULL;
    if (!(q->node = malloc(sizeof(isax_node *)*k)))
        return NULL;

    for (int i = 0; i < k; ++i)
    {
        q->knn[i]=FLT_MAX;
    }
    q->k=k;
    q->nowk=0;
    return q;
}

 

/* Function to insert value into priority queue */



 

/* Function to check priority and place element */

void pqueue_bsf_insert(pqueue_bsf *q,float data,long int position, isax_node *node)
{

    int i,j;

    for (i = 0; i < q->k; i++)
    {

        if (data <= q->knn[i])
        {
            //for (int w=0; w<q->k; w++) {
                //printf("[%d] %f\n", w, q->knn[w]);
            //}

            for (j = q->k-1; j > i; j--)
            {

                q->knn[j] =q->knn[j - 1];
                q->position[j]=q->position[j - 1];
                q->node[j] = q->node[j - 1];
            }

            q->knn[i] = data;
            q->position[i]=position;
            q->node[i] = node;

            //for (int w=0; w<q->k; w++) {
                //printf("[%d] %f\n", w, q->knn[w]);
            //}

            return;

        }

    }
}


void pqueue_bsfre_insert(pqueue_bsf *q,float data,long int position, isax_node *node)
{

    int i,j;
    if(q->nowk==q->k)
    {
        q->k=q->k*2;
        realloc(q->position,sizeof(long int)*(q->k));
        realloc(q->knn ,sizeof(float)*(q->k));
        

    }
    for (i = 0; i < q->nowk; i++)
    {
        if (data <= q->knn[i])
        {

            for (j = q->nowk; j > i; j--)
            {

                q->knn[j] =q->knn[j - 1];
                q->position[j]=q->position[j - 1];
            }
            q->knn[i] = data;
            q->position[i]=position;
            //for (int w=0; w<q->k; w++) {
                //printf("[%d] %f\n", w, q->knn[w]);
            //}
            q->nowk++;
            return;
        }
    }
    q->knn[q->nowk] = data;
    q->position[q->nowk]=position;
    q->nowk++;
}
 

/* Function to delete an element from queue */



 

/* Function to display queue elements */





