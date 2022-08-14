//
//  hercules_file_loaders.h
//
//  Created by Karima Echihabi on 18/12/2018
//


#ifndef hercules_hercules_file_loaders_h
#define hercules_hercules_file_loaders_h
#include "../config.h"
#include "../globals.h"
#include "ts.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hercules_index.h"
#include "hercules_query_engine.h"
#include "calc_utils.h"

enum response hercules_query_ascii_file(const char *ifilename, int q_num, const char delimiter, struct hercules_index *index, float minimum_distance, ts_type epsilon, ts_type delta);
enum response hercules_query_binary_file(const char *ifilename, int q_num,
				       struct hercules_index *index, float minimum_distance,
				       ts_type epsilon, ts_type delta);  
enum response hercules_knn_query_binary_file(const char *ifilename, int q_num, struct hercules_index *index,
					   float minimum_distance, ts_type epsilon, ts_type r_delta,
					   unsigned int k, boolean track_bsf, boolean track_pruning,
					   boolean dump_mindists, boolean max_policy,
					   unsigned int nprobes, unsigned char incremental,
					   query_settings q_settings,
					   int q_skip);
enum response hercules_index_binary_file(const char *ifilename, file_position_type ts_num, struct hercules_index *index);
enum response hercules_index_ascii_file(const char *ifilename, file_position_type ts_num, const char delimiter, struct hercules_index *index);
enum response reorder_query(ts_type * query_ts, ts_type * query_ts_reordered, int * query_order, int ts_length);
enum response hercules_tlb_binary_file(const char *ifilename, int q_num, struct hercules_index *index,float minimum_distance);
enum response hercules_index_binary_file_p(const char *ifilename, file_position_type ts_num,
					 struct hercules_index *index, int num_threads, int start_parallel);
void* indexbulkloadingworker(void *transferdata);

#endif
