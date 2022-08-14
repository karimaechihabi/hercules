//
//  hercules_node_split.h
//
//  Created by Karima Echihabi on 18/12/2018
//

#ifndef herculeslib_hercules_node_split_h
#define herculeslib_hercules_node_split_h
#include "../config.h"
#include "../globals.h"
#include "hercules_index.h"

struct segment_sketch{
  ts_type * indicators;
  int num_indicators;
  pthread_mutex_t * lock_synopsis;
}; 

struct node_segment_split_policy {
  ts_type indicator_split_value;
  int indicator_split_idx; //Set Idx = 0 for mean_based split and Idx = 1 for stdev_based split 
};

struct  node_split_policy {
  short split_from;
  short split_to;
  struct node_segment_split_policy curr_node_segment_split_policy;

  int indicator_split_idx;
  ts_type indicator_split_value;

  struct segment_sketch series_segment_sketch;
   
};

enum response split_node_create_children (struct hercules_index * index, struct hercules_node * node, short * child_node_points, int num_child_node_points);
enum response node_segment_sketch_update_sketch(struct segment_sketch * node_sketch, ts_type * series, int fromIdx, int toIdx);

enum response series_segment_sketch_do_sketch(struct segment_sketch * series_segment_sketch, ts_type * series, int fromIdx, int toIdx);

enum response mean_node_segment_split_policy_split(struct node_segment_split_policy * policy, struct segment_sketch sketch, struct segment_sketch * ret);
enum response stdev_node_segment_split_policy_split(struct node_segment_split_policy * policy, struct segment_sketch sketch, struct segment_sketch * ret);

boolean is_split_policy_mean(struct node_segment_split_policy policy);
boolean is_split_policy_stdev(struct node_segment_split_policy policy);
struct hercules_node *  create_child_node(struct hercules_node * parent);
ts_type range_calc(struct segment_sketch sketch, int len);
//boolean node_split_policy_route_to_left (struct hercules_node * node, ts_type * series);
boolean node_split_policy_route_to_left (struct hercules_node * node, ts_type * series, struct segment_sketch * series_segment_sketch);
short get_hs_split_point(short * points, short from, short to, int size_points);
struct segment_sketch * init_segment_sketches(int num_child_segments, int num_indicators);
void  split_node_distribute_data(struct hercules_index * index, struct hercules_node *node, struct segment_sketch * sketch);
void  split_node (struct hercules_index * index, struct hercules_node *node, struct segment_sketch * sketch, void * thread_data);
void split_node_evaluate_policies (struct hercules_index * index, struct hercules_node *node, int * hs_split_point);
enum response node_segment_sketch_update_sketch_parallel(struct segment_sketch * node_segment_sketch, ts_type * series, int fromIdx, int toIdx,struct segment_sketch * series_segment_sketch);

#endif
