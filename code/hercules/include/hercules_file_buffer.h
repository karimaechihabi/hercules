//
//  hercules_file_buffer.h
//
//  Created by Karima Echihabi on 18/12/2018
//

#ifndef hercules_hercules_file_buffer_h
#define hercules_hercules_file_buffer_h
#include "../config.h"
#include "../globals.h"
#include "ts.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hercules_node.h"
#include "hercules_index.h"


struct hercules_file_buffer {

  struct hercules_node * node; //the buffer points back to its node
  struct hercules_file_map * position_in_map; //the buffer points back to its position in file map

  ts_type ** buffered_list;
  unsigned int disk_count; //  by default
  int buffered_list_size;   //number of series currently stored in this buffer

  boolean in_disk; //false by default
  boolean do_not_flush;
};

enum response hercules_file_buffer_init(struct hercules_node *node);
enum response flush_buffer_to_disk(struct hercules_index *index, struct hercules_node *node);
enum response clear_file_buffer(struct hercules_index *index, struct hercules_node * node);
enum response delete_file_buffer(struct hercules_index * index,struct hercules_node * node);
//enum response get_all_time_series_in_node(struct hercules_index * index, struct hercules_node * node);
ts_type ** get_all_time_series_in_node(struct hercules_index * index, struct hercules_node * node,int serial);
void * get_all_time_series_in_node_parallel(struct hercules_index * index, struct hercules_node * node, int serial, void * tdata);

enum response flush_leaf_to_leaves_file(struct hercules_index *index, struct hercules_node *node, int sims);
enum response flush_leaf_to_leaves_file_update_stats_serial(struct hercules_index *index, struct hercules_node *node, int sims);

enum response flush_leaf_to_leaves_file_update_stats_parallel(struct hercules_index *index, struct hercules_node *node, int sims,void * tdata);
void * flush_leaf_worker(void *transferdata);
enum response hercules_index_flush_leaves(struct hercules_index *index, int num_threads);
void * hercules_index_flush_leaf_worker(void *transferdata);



#endif				   
