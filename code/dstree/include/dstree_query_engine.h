//
//  dstree_query_engine.h
//  ds-tree C version
//
//  Created by Karima Echihabi on 18/12/2016
//  Copyright 2012 Paris Descartes University. All rights reserved.
//  
//  Based on isax code written by Zoumpatianos on 3/12/12.
//  Copyright 2012 University of Trento. All rights reserved.
//


#ifndef al_dstree_dstree_query_engine_h
#define al_dstree_dstree_query_engine_h
#include "../config.h"
#include "../globals.h"

#include "dstree_index.h"
#include "dstree_node.h"

typedef struct query_result {
    ts_type distance;
    struct dstree_node *node;
    size_t pqueue_position;
};

/// Data structure for sorting the query.
typedef struct q_index
    {   double value;
        int    index;
    } q_index;

static int cmp_pri(double next, double curr)
{
	return (next > curr);
}


static double
get_pri(void *a)
{
	return (double) ((struct query_result *) a)->distance;
}


static void
set_pri(void *a, double pri)
{
	((struct query_result *) a)->distance = (float)pri;
}


static size_t
get_pos(void *a)
{
	return ((struct query_result *) a)->pqueue_position;
}


static void
set_pos(void *a, size_t pos)
{
	((struct query_result *) a)->pqueue_position = pos;
}



//struct query_result approximate_search (ts_type *ts, struct dstree_index *index);
//struct query_result exact_search (ts_type *ts, struct dstree_index *index,ts_type minimum_distance);

struct query_result exact_search (ts_type *query_ts, ts_type * query_reordered, int * query_order, unsigned int offset, struct dstree_index *index,ts_type minimum_distance);
struct query_result exact_search_leaves_file (ts_type *query_ts, ts_type * query_ts_reordered, int * query_order, unsigned int offset, struct dstree_index *index,ts_type minimum_distance);
struct query_result exact_serial_search_leaves_file (ts_type *query_ts, ts_type * query_ts_reordered, int * query_order, unsigned int offset, struct dstree_index *index,ts_type minimum_distance);


struct query_result approximate_search (ts_type *query_ts, ts_type * query_ts_reordered, int * query_order, unsigned int offset, ts_type bsf, struct dstree_index *index);
void dstree_calc_tlb (ts_type *query_ts, struct dstree_index *index, struct dstree_node * curr_node);


#endif
