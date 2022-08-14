//
//  hercules_node.h
//
//  Created by Karima Echihabi on 18/12/2018
//

#ifndef herculeslib_hercules_node_h
#define herculeslib_hercules_node_h


#include "../config.h"
#include "../globals.h"    
#include "hercules_node_split.h"
#include "hercules_file_buffer.h"
#include "pqueue.h"
#include "hercules_query_engine.h"
#include <pthread.h>

struct hercules_node {

  struct node_segment_split_policy * node_segment_split_policies;
  short * node_points;
  short * hs_node_points;
  short split_segment;
  struct segment_sketch * node_segment_sketches;
  struct segment_sketch * hs_node_segment_sketches;
  struct node_split_policy * split_policy; 

  struct hercules_node *left_child;
  struct hercules_node *right_child;
  struct hercules_node *parent;

  struct hercules_file_buffer * file_buffer;

  char * filename;

  mean_stdev_range range;  
 
  int num_node_segment_split_policies;  
  short num_node_points;  //number of vertical split points
  short num_hs_node_points;  //number of horizontal split points

  int max_segment_length; 
  int max_value_length; 
  
  unsigned int node_size;
  unsigned int level;

  unsigned char is_leaf;
  boolean is_left;
  
  unsigned long long file_pos;
  unsigned char is_pruned;

  pthread_mutex_t lock_data;
  pthread_mutex_t lock_split;  
  unsigned char is_splitting;

  int proc_finished;
  int write_finished;

  int thread_id;
  
};

struct hercules_node * hercules_root_node_init();
struct hercules_node * hercules_leaf_node_init();
enum response node_init_segments(struct hercules_node * node, short * split_points, int segment_size);

enum response append_ts_to_node(struct hercules_index * index, struct hercules_node * node, ts_type * timeseries);  
enum response update_node_statistics(struct hercules_node * node, ts_type * timeseries);
enum response update_node_statistics_parallel(struct hercules_node * node, ts_type * timeseries,struct segment_sketch * series_segment_sketch);

enum response create_hercules_node_filename(struct hercules_index_settings *settings,
                                          struct hercules_node * node, struct hercules_node * parent_node);

//ts_type calculate_ts_in_node_distance (struct hercules_index *index, struct hercules_node *node, ts_type *query);
//ts_type calculate_node_distance (struct hercules_index *index, struct hercules_node *node, ts_type *query);
ts_type calculate_node_min_distance (struct hercules_index *index, struct hercules_node *node, ts_type *query);
ts_type calculate_ts_in_node_distance (struct hercules_index *index,
				       struct hercules_node *node,
				       ts_type *query_ts_reordered,
				       int * query_order,
				       unsigned int offset,
				       ts_type bound);
ts_type calculate_node_distance (struct hercules_index *index, struct hercules_node *node, ts_type *query_ts_reordered, int *query_order, unsigned int offset, ts_type bsf);

void calculate_node_knn_distance (struct hercules_index *index, struct hercules_node *node,
				  ts_type *query_ts,
				  ts_type *query_ts_reordered, int *query_order,
				  unsigned int offset, ts_type bsf,
				  unsigned int k,
				  struct query_result  *knn_results,
				  struct bsf_snapshot ** bsf_snapshots,
				  unsigned int *cur_bsf_snapshot,
				  unsigned int * cur_size,
				  int serial);

boolean calculate_node_knn_distance_psq (struct hercules_index *index, struct hercules_node *node,
					 ts_type *query_ts,
					 ts_type *query_ts_reordered,				  
					 int *query_order,
					 unsigned int offset, ts_type bsf,
					 unsigned int k,
					 struct query_result  *knn_results,
					 struct bsf_snapshot ** bsf_snapshots,
					 unsigned int *cur_bsf_snapshot,
					 unsigned int * cur_size,
					 void* qwdata);


enum response update_node_ancestors_statistics_for_non_split_segments(struct hercules_node * leaf_node);
enum response update_node_ancestors_statistics_for_split_segment(struct hercules_node * leaf_node, ts_type * timeseries);
  
enum response update_node_ancestors_statistics_for_non_split_segments_parallel(struct hercules_node * leaf_node);
enum response update_node_ancestors_statistics_for_split_segment_parallel(struct hercules_node * node, ts_type * timeseries,struct segment_sketch * series_segment_sketch);

enum response hercules_index_node_update_synopsis(struct hercules_index *index, struct hercules_node *node, int sims,void * tdata);
enum response hercules_index_node_update_synopsis_hsplit(struct hercules_node * leaf_node);
enum response hercules_index_node_update_synopsis_vsplit(struct hercules_node * node, ts_type * timeseries,struct segment_sketch * series_segment_sketch);



#endif
