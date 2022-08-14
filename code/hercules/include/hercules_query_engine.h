//
//  hercules_query_engine.h
//
//  Created by Karima Echihabi on 18/12/2018
//


#ifndef al_hercules_hercules_query_engine_h
#define al_hercules_hercules_query_engine_h
#include "../config.h"
#include "../globals.h"

#include "hercules_index.h"

#include "hercules_node.h"
#include "pqueue.h"

typedef struct query_result {
  ts_type distance;
  struct hercules_node *node;
  ts_type max_distance;
  size_t pqueue_position;
};

typedef struct bsf_snapshot {
  ts_type distance;
  double time;  
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

static double
get_max_pri(void *a)
{
	return (double) ((struct query_result *) a)->max_distance;
}

static void
set_pri(void *a, double pri)
{
	((struct query_result *) a)->distance = (float)pri;
}

static void
set_max_pri(void *a, double pri)
{
	((struct query_result *) a)->max_distance = (float)pri;
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


typedef struct localStack {
    struct hercules_node **val; 
    int top;
    int bottom;
}localStack;

typedef struct siss_query_worker_data
{
  ts_type *query_ts,*query_ts_reordered;
  int * query_order;
  unsigned int offset;
  ts_type *query_paa;  
  ts_type minimum_distance;
  ts_type epsilon;
  ts_type delta;
  unsigned int k;
  pqueue_t *pq;
  struct query_result * candidates;
  unsigned int *candidates_series_count;
  unsigned int *candidates_count;
  unsigned int *candidates_idx;
  ts_type * ts_list;
  struct hercules_index *index;  
  //pthread_mutex_t *lock_current_root_node;
  pthread_mutex_t *lock_queue;
  pthread_barrier_t *lock_barrier;
  pthread_rwlock_t *lock_bsf;
  struct query_result *knn_results;
  int * node_counter;
  int node_list_size;
  bool lockvalueq;
  bool lockvaluebsf;
  struct hercules_node ** node_list;
  int * curr_k_size;

  double thread_input_time;
  double thread_realdist_time;
  double thread_insert_node_time;
  double thread_random_io;
  double thread_sequential_io;
  
  unsigned long *load_point;
  unsigned long * label_number;
  float* minidisvector;
  int sum_of_lab;

} siss_query_worker_data;

typedef struct query_settings {
  int serial;
  int sims;
  int psq;
  int pmq;
  int num_queues;
  int num_threads;
  int num_read_threads;
  float approx_stop_condition;
  float exact_stop_condition;
  float eapca_threshold;
  float sax_threshold;
}query_settings;


//struct query_result approximate_search (ts_type *ts, struct hercules_index *index);
//struct query_result exact_search (ts_type *ts, struct hercules_index *index,ts_type minimum_distance);

struct query_result exact_search (ts_type *query_ts, ts_type * query_reordered, int * query_order, unsigned int offset, struct hercules_index *index,ts_type minimum_distance, ts_type epsilon, ts_type delta);

struct query_result approximate_search (ts_type *query_ts, ts_type * query_ts_reordered, int * query_order, unsigned int offset, ts_type bsf, struct hercules_index *index);

void approximate_knn_search (ts_type *query_ts, ts_type * query_ts_reordered,
			     int * query_order, unsigned int offset,
			     ts_type bsf, struct hercules_index *index,
			     struct query_result *knn_results,
			     unsigned int k, struct bsf_snapshot ** bsf_snapshots,
			     unsigned int * cur_bsf_snapshot,
			     unsigned int * curr_size, int serial);

void hercules_calc_tlb (ts_type *query_ts, struct hercules_index *index, struct hercules_node * curr_node);

void update_query_stats(struct hercules_index * index, unsigned int query_id, unsigned int found_knn, struct query_result bsf_result);
void print_query_stats(struct hercules_index * index, unsigned int query_num, unsigned int found_knn, char * queries);
void get_query_stats(struct hercules_index * index, unsigned int found_knn);
void print_bsf_snapshots(struct hercules_index * index, unsigned int query_num, unsigned int k,
			 char * queries, struct bsf_snapshot **bsf_snapshots, unsigned int cur_bsf_snapshot);

void  exact_knn_search_track_bsf (ts_type *query_ts, ts_type * query_ts_reordered, int * query_order,
				  unsigned int offset, struct hercules_index *index,ts_type minimum_distance,
				  ts_type epsilon, ts_type delta, unsigned int k,
				  unsigned int q_id, char * qfilename, 
				  struct bsf_snapshot ** bsf_snapshots, unsigned int * cur_bsf_snapshot,query_settings q_settings);
void  exact_knn_search_track_pruning (ts_type *query_ts, ts_type * query_ts_reordered,
				      int * query_order, unsigned int offset,
				      struct hercules_index *index,ts_type minimum_distance,
				      ts_type epsilon, ts_type delta,
				      unsigned int k, unsigned int q_id, char * qfilename, query_settings q_settings);

void  exact_knn_search_max_policy (ts_type *query_ts, ts_type * query_ts_reordered,
				   int * query_order, unsigned int offset,
				   struct hercules_index *index,ts_type minimum_distance,
				   ts_type epsilon, ts_type delta,
				   unsigned int k, unsigned int q_id, char * qfilename, query_settings q_settings);

void print_pruning_snapshots(struct hercules_node * node,
			     ts_type node_bsf,
			     ts_type node_mindist,
			     unsigned int k,
			     unsigned int query_num,
			     char * queries);

void dump_mindists (struct hercules_index *index,
		    struct hercules_node *node,
		    ts_type *query_ts);
ts_type get_node_QoS(struct hercules_index * index, struct hercules_node * node);


void  exact_incr_knn_search (ts_type *query_ts, ts_type * query_ts_reordered,
			int * query_order, unsigned int offset,
			struct hercules_index *index,ts_type minimum_distance,
			ts_type epsilon, ts_type r_delta,
			     unsigned int k, unsigned int q_id, char * qfilename, unsigned int nprobes, query_settings q_settings);

void  exact_de_knn_search (ts_type *query_ts, ts_type * query_ts_reordered,
			int * query_order, unsigned int offset,
			struct hercules_index *index,ts_type minimum_distance,
			ts_type epsilon, ts_type r_delta,
			   unsigned int k, unsigned int q_id, char * qfilename, query_settings q_settings);

void  exact_ng_knn_search (ts_type *query_ts, ts_type * query_ts_reordered,
			int * query_order, unsigned int offset,
			struct hercules_index *index,ts_type minimum_distance,
			   unsigned int k, unsigned int q_id, char * qfilename, unsigned int nprobes,query_settings q_settings);

void exact_de_knn_search_populate_queue(ts_type *query_ts,
					struct hercules_index *index,
					struct query_result * knn_results,
					int k,
					pqueue_t *pq,
					ts_type delta,
					ts_type epsilon,
					pqueue_t * candidate_list,
					int stop_condition);

void exact_de_knn_search_populate_queue_lf(ts_type *query_ts,
					   struct hercules_index *index,
					   struct query_result * knn_results,
					   int k,
					   pqueue_t *pq,
					   ts_type delta,
					   ts_type epsilon,
					   struct query_result * candidates,
					   unsigned int *candidates_count,
					   int stop_condition);

void* exact_de_knn_search_worker(void *qwdata);
void  exact_de_knn_search_psq(ts_type *query_ts, ts_type * query_ts_reordered,
			      int * query_order, unsigned int offset,
			      struct hercules_index *index,ts_type minimum_distance,
			      ts_type epsilon, double delta,
			      unsigned int k, unsigned int q_id, char * qfilename,
			      query_settings q_settings);
void search_worker_insert_node(ts_type *query_ts,
			       struct hercules_index *index,
			       struct hercules_node *node,
			       struct query_result * knn_results,
			       int k,
			       pqueue_t *pq,
			       ts_type delta,
			       ts_type epsilon,
			       pthread_mutex_t *lock_queue);

int queue_bounded_sorted_insert(struct  query_result *q, struct query_result d, unsigned int *cur_size, unsigned int k);

void* exact_de_knn_read_worker_lf(void *qwdata);
void* exact_de_knn_search_worker_lf(void *qwdata);
void  exact_de_knn_search_psq_lf (ts_type *query_ts, ts_type * query_ts_reordered,
				  int * query_order, unsigned int offset,
				  struct hercules_index *index,ts_type minimum_distance,
				  ts_type epsilon, double delta,
				  unsigned int k, unsigned int q_id, char * qfilename,
				  query_settings q_settings);

void search_worker_insert_node_lf(ts_type *query_ts,
				  struct hercules_index *index,
				  struct hercules_node *node,
				  struct query_result * knn_results,
				  int k,
				  struct query_result *candidates,
				  unsigned int * candidates_count,
				  unsigned int * candidates_series_count,
				  ts_type delta,
				  ts_type epsilon,
				  pthread_mutex_t *lock_queue);

void* exact_de_knn_read_worker_lf_sims(void *qwdata);
void  exact_de_knn_search_psq_lf_sims (ts_type *query_ts, ts_type * query_ts_reordered,
				       int * query_order, unsigned int offset, ts_type *query_paa,
				       struct hercules_index *index,ts_type minimum_distance,
				       ts_type epsilon, double delta,
				       unsigned int k, unsigned int q_id, char * qfilename,
				       query_settings q_settings);
void* exact_de_knn_mindist_worker_lf(void *qwdata);

void  exact_de_knn_search_psq_lf_sims_adaptive (ts_type *query_ts, ts_type * query_ts_reordered,
				       int * query_order, unsigned int offset, ts_type *query_paa,
				       struct hercules_index *index,ts_type minimum_distance,
				       ts_type epsilon, double delta,
				       unsigned int k, unsigned int q_id, char * qfilename,
				       query_settings q_settings);
void  exact_de_knn_search_noapprox (ts_type *query_ts, ts_type * query_ts_reordered,
			   int * query_order, unsigned int offset,
			   struct hercules_index *index,ts_type minimum_distance,
			   ts_type epsilon, ts_type r_delta,
			   struct query_result * knn_results, unsigned int k, unsigned int q_id, char * qfilename, query_settings q_settings,
 			   struct query_result * candidates, unsigned int candidates_count, unsigned int *candidates_idx,
               ts_type *ts_list, FILE *ifile);
void approximate_greedy_knn_search (ts_type *query_ts, ts_type * query_ts_reordered,
			     int * query_order, unsigned int offset,
			     ts_type bsf, struct hercules_index *index,
			     struct query_result *knn_results,
			     unsigned int k, struct bsf_snapshot ** bsf_snapshots,
			     unsigned int * cur_bsf_snapshot,
			     unsigned int * curr_size,
			     int serial,
				 pqueue_t *pq,
				 FILE *fp,
                 ts_type *ts_list,
				 ts_type epsilon,
                 ts_type delta,
                 float stop_condition);

#endif
