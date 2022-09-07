//
//  hercules_query_engine.c
//
//  Created by Karima Echihabi on 18/12/2018
//


#include "../config.h"
#include "../globals.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/hercules_query_engine.h"
#include "../include/hercules_file_buffer.h"
#include "../include/hercules_file_buffer_manager.h"
#include "../include/hercules_index.h"
#include "../include/hercules_node.h"
#include "../include/sax.h"

#include "../include/pqueue.h"
#ifdef VALUES
#include <values.h>
#endif
#include <pthread.h>

struct query_result approximate_search (ts_type *query_ts, ts_type * query_ts_reordered, int * query_order, unsigned int offset, ts_type bsf, struct hercules_index *index) {

  struct query_result result;
  struct hercules_node * node = index->first_node;
  //ts_type bsf = FLT_MAX;  //no bsf known so far
    
  //Allocate memory for the series sketch
  struct segment_sketch timeseries_segment_sketch;
  timeseries_segment_sketch.indicators = NULL;
  timeseries_segment_sketch.indicators = malloc (sizeof(ts_type)*2);
  if(timeseries_segment_sketch.indicators == NULL) {
    fprintf(stderr,"Error in hercules_index.c: Could not allocate memory for" 
	    "series segment sketch indicators.\n");
  }

  timeseries_segment_sketch.num_indicators = 2;

  if (node != NULL) {
    // Traverse tree
    while (!node->is_leaf) {
      if(node_split_policy_route_to_left(node,query_ts,&timeseries_segment_sketch)) 
	{
	  node = node->left_child;
	}
      else
	{
	  node = node->right_child;
	}
    }

    result.distance = calculate_node_distance(index, node, query_ts_reordered, query_order,offset,bsf);
    result.node = node;
  }
  else {
    printf("Error: index is empty \n");        
    result.node = NULL;
    result.distance = MAXFLOAT;
  }
    
    
  return result;
}
void approximate_greedy_knn_search (ts_type *query_ts, ts_type * query_ts_reordered,
				    int * query_order, unsigned int offset,
				    ts_type bsf, struct hercules_index *index,
				    struct query_result *knn_results,
				    unsigned int k, struct bsf_snapshot ** bsf_snapshots,
				    unsigned int * cur_bsf_snapshot,
				    unsigned int * curr_size,
				    int serial,
				    pqueue_t *pq,
				    FILE *ifile,
				    ts_type *ts_list,
				    ts_type epsilon,
				    ts_type delta,
				    float approx_stop_condition) 
{

  curr_size = 0;
  bsf = FLT_MAX;
  ts_type kth_bsf = FLT_MAX;
  ts_type temp_bsf = FLT_MAX;
  ts_type old_kth_bsf = FLT_MAX;
  ts_type distance = FLT_MAX;

  fseek(ifile, 0, SEEK_SET);    
  unsigned int found_knn = 0;
 
  struct query_result * n;
  struct query_result temp;

  int ts_length = index->settings->timeseries_size;

  int max_visited_leaves =  (int)approx_stop_condition;
  while ((n = pqueue_pop(pq)))
    {
      kth_bsf =  knn_results[k-1].distance;
      if (n->distance > kth_bsf)
	{
          break;
	}

      if (n->node->is_leaf) // n is a leaf
        {
          COUNT_LOADED_NODE
	    COUNT_LOADED_TS(n->node->node_size)
	    COUNT_PARTIAL_INPUT_TIME_START       
	    fseek(ifile,
		  n->node->file_pos * ts_length * sizeof(ts_type),
		  SEEK_SET);
		  
	  fread(ts_list,sizeof(ts_type),
		ts_length*n->node->node_size,
		ifile);	    
	  COUNT_PARTIAL_INPUT_TIME_END       
		  
	    for (int idx = 0; idx < n->node->node_size; ++idx)
	      {
		kth_bsf =  knn_results[k-1].distance;
		distance = ts_euclidean_distance_SIMD(query_ts,
						      &ts_list[idx*ts_length],
						      ts_length,
						      kth_bsf);
		if (distance < kth_bsf)
		  {
		    struct query_result object_result;// =  malloc(sizeof(struct query_result));
		    object_result.node = n->node;
		    object_result.distance =  distance;				
		    queue_bounded_sorted_insert(knn_results, object_result, &curr_size, k);
		  }	    
	      }	
	  if (loaded_nodes_count >= max_visited_leaves )
	    break;      

        }
      // If it is an intermediate node calculate mindist for children
      // and push them in the queue
      else  //n is an internal node
        {
          temp = knn_results[k-1];
          kth_bsf =  temp.distance;
	  
          ts_type child_distance;
     	  child_distance = calculate_node_min_distance(index, n->node->left_child,query_ts);
	  
	  if ((child_distance < kth_bsf)) 
	    {
	      struct query_result * mindist_result_left = malloc(sizeof(struct query_result));
	      mindist_result_left->node = n->node->left_child;
	      mindist_result_left->distance =  child_distance;
	      pqueue_insert(pq, mindist_result_left);
	    }

	  child_distance = calculate_node_min_distance(index, n->node->right_child,query_ts);

	  if ((child_distance < kth_bsf) ) 	  
	    {
	      struct query_result * mindist_result_right = malloc(sizeof(struct query_result));
	      mindist_result_right->node = n->node->right_child;
	      mindist_result_right->distance =  child_distance;
	      pqueue_insert(pq, mindist_result_right);
	    }
        }
      free(n);
    }

}

void approximate_knn_search (ts_type *query_ts, ts_type * query_ts_reordered,
			     int * query_order, unsigned int offset,
			     ts_type bsf, struct hercules_index *index,
			     struct query_result *knn_results,
			     unsigned int k, struct bsf_snapshot ** bsf_snapshots,
			     unsigned int * cur_bsf_snapshot,
			     unsigned int * curr_size,
			     int serial) {

  struct query_result result;
  struct hercules_node * node = index->first_node;

  struct segment_sketch timeseries_segment_sketch;
  timeseries_segment_sketch.indicators = NULL;
  timeseries_segment_sketch.indicators = malloc (sizeof(ts_type)*2);
  if(timeseries_segment_sketch.indicators == NULL) {
    fprintf(stderr,"Error in hercules_index.c: Could not allocate memory for" 
	    "series segment sketch indicators.\n");
  }

  timeseries_segment_sketch.num_indicators = 2;
    
    
  if (node != NULL)
    {
      // Traverse tree
      while (!node->is_leaf) {
	if(node_split_policy_route_to_left(node,query_ts,&timeseries_segment_sketch)) 
	  {
	    node = node->left_child;
	  }
	else
	  {
	    node = node->right_child;
	  }
      }

      calculate_node_knn_distance(index, node, query_ts,query_ts_reordered, query_order,offset,
				  bsf,k,knn_results,bsf_snapshots, cur_bsf_snapshot,curr_size,serial);
    }
  else
    {
      printf("Error: index is empty \n");        
    }
  free(timeseries_segment_sketch.indicators);

}


struct query_result exact_search (ts_type *query_ts, ts_type * query_ts_reordered, int * query_order, unsigned int offset, struct hercules_index *index,ts_type minimum_distance, ts_type epsilon, ts_type delta) {
    
  ts_type bsf = FLT_MAX;


  struct query_result approximate_result = approximate_search(query_ts, query_ts_reordered, query_order, offset, bsf, index);

  struct query_result bsf_result = approximate_result;    

  COUNT_PARTIAL_TIME_END                 
    index->stats->query_filter_total_time  = partial_time;	
    
  index->stats->query_filter_input_time  = partial_input_time;
  index->stats->query_filter_output_time = partial_output_time;
  index->stats->query_filter_load_node_time = partial_load_node_time;
  index->stats->query_filter_cpu_time    = partial_time-partial_input_time-partial_output_time;
  index->stats->query_filter_seq_input_count   = partial_seq_input_count;
  index->stats->query_filter_seq_output_count  = partial_seq_output_count;
  index->stats->query_filter_rand_input_count  = partial_rand_input_count;
  index->stats->query_filter_rand_output_count = partial_rand_output_count;
    
  index->stats->query_filter_loaded_nodes_count = loaded_nodes_count;
  index->stats->query_filter_loaded_ts_count = loaded_ts_count;
  index->stats->query_filter_checked_nodes_count = checked_nodes_count;
  index->stats->query_filter_checked_ts_count = checked_ts_count;

  index->stats->query_approx_distance = sqrtf(approximate_result.distance);
  index->stats->query_approx_node_filename = approximate_result.node->filename;
  index->stats->query_approx_node_size = approximate_result.node->node_size;;
  index->stats->query_approx_node_level = approximate_result.node->level;
    
  index->stats->queries_filter_total_time  +=  index->stats->query_filter_total_time;
    
  index->stats->queries_filter_input_time  +=  index->stats->query_filter_input_time;
  index->stats->queries_filter_output_time +=  index->stats->query_filter_output_time;
  index->stats->queries_filter_load_node_time += index->stats->query_filter_load_node_time;
  index->stats->queries_filter_cpu_time    += index->stats->query_filter_cpu_time;

  index->stats->queries_filter_seq_input_count   += index->stats->query_filter_seq_input_count;
  index->stats->queries_filter_seq_output_count  += index->stats->query_filter_seq_output_count;
  index->stats->queries_filter_rand_input_count  += index->stats->query_filter_rand_input_count;
  index->stats->queries_filter_rand_output_count += index->stats->query_filter_rand_output_count;
    
    
  if(approximate_result.node != NULL) {  
    index->stats->query_approx_distance = sqrtf(approximate_result.distance);
    index->stats->query_approx_node_filename = approximate_result.node->filename;
    index->stats->query_approx_node_size = approximate_result.node->node_size;;
    index->stats->query_approx_node_level = approximate_result.node->level;
  }
 
  // Early termination...
  if (approximate_result.distance == 0) {
    index->stats->query_exact_distance = sqrtf(approximate_result.distance);
    index->stats->query_exact_node_filename = approximate_result.node->filename;
    index->stats->query_exact_node_size = approximate_result.node->node_size;;
    index->stats->query_exact_node_level = approximate_result.node->level; 
    return approximate_result;
  }

  RESET_QUERY_COUNTERS()
    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
 
    
   
    pqueue_t *pq = pqueue_init(index->first_node->node_size, 
                               cmp_pri, get_pri, set_pri, get_pos, set_pos);
	
  if(approximate_result.node != NULL) {
    // Insert approximate result in heap.
    pqueue_insert(pq, &approximate_result);
  }
    
  struct query_result *do_not_remove = &approximate_result;


  //Add the root to the priority queue

  struct query_result * root_pq_item = malloc(sizeof(struct query_result));
  root_pq_item->node = index->first_node;
  root_pq_item->distance = calculate_node_min_distance (index, index->first_node, query_ts);
    
  //initialize the lb distance to be the distance of the query to the root.

  pqueue_insert(pq, root_pq_item);    

  struct query_result * n;
  while ((n = pqueue_pop(pq)))
    {
      if (n->distance > bsf_result.distance/(1 + epsilon))
	{
          break;
        }
      if (n->node->is_leaf) // n is a leaf
        {
   
	  ts_type distance = calculate_node_distance(index, n->node, query_ts_reordered, query_order, offset, bsf_result.distance);

	  if (distance < bsf_result.distance)
	    {
	      bsf_result.distance = distance;
	      bsf_result.node = n->node;
	    }
        }
      // If it is an intermediate node calculate mindist for children
      // and push them in the queue
      else  //n is an internal node
        {
	        
	  ts_type child_distance;
	  child_distance = calculate_node_min_distance(index, n->node->left_child,query_ts);

	  if (child_distance < bsf_result.distance/(1 + epsilon) )
	    {
	      struct query_result * mindist_result_left = malloc(sizeof(struct query_result));
	      mindist_result_left->node = n->node->left_child;
	      mindist_result_left->distance =  child_distance;
	      pqueue_insert(pq, mindist_result_left);
	    }

	  child_distance = calculate_node_min_distance(index, n->node->right_child,query_ts);

	  if (child_distance < bsf_result.distance/(1 + epsilon) )
	    {
	      struct query_result * mindist_result_right = malloc(sizeof(struct query_result));
	      mindist_result_right->node = n->node->right_child;
	      mindist_result_right->distance =  child_distance;
	      pqueue_insert(pq, mindist_result_right);
	    }
        }
      // Free the node currently popped.
      if(n != do_not_remove)
	free(n);
    }
    
  // Free the nodes that were not popped.
  while ((n = pqueue_pop(pq)))
    {
      if(n != do_not_remove)
	free(n);
    }
  // Free the priority queue.
    
  pqueue_free(pq);


  COUNT_PARTIAL_TIME_END 


    index->stats->query_refine_total_time  = partial_time;	
        
  index->stats->query_refine_input_time  = partial_input_time;
  index->stats->query_refine_output_time = partial_output_time;
  index->stats->query_refine_load_node_time = partial_load_node_time;    
  index->stats->query_refine_cpu_time    = partial_time
    - partial_input_time
    - partial_output_time;
  index->stats->query_refine_seq_input_count   = partial_seq_input_count;
  index->stats->query_refine_seq_output_count  = partial_seq_output_count;
  index->stats->query_refine_rand_input_count  = partial_rand_input_count;
  index->stats->query_refine_rand_output_count = partial_rand_output_count;

  index->stats->query_total_time  = partial_time
    +  index->stats->query_filter_total_time;
  index->stats->query_total_input_time  = partial_input_time
    +  index->stats->query_filter_input_time;
  index->stats->query_total_output_time  = partial_output_time
    +  index->stats->query_filter_output_time;
  index->stats->query_total_load_node_time  = partial_load_node_time
    +  index->stats->query_filter_load_node_time;
  index->stats->query_total_cpu_time  =  index->stats->query_total_time
    -  index->stats->query_total_input_time
    -  index->stats->query_total_output_time;
    
  index->stats->query_total_seq_input_count   = partial_seq_input_count
    + index->stats->query_filter_seq_input_count;
  index->stats->query_total_seq_output_count   = partial_seq_output_count
    + index->stats->query_filter_seq_output_count;
  index->stats->query_total_rand_input_count   = partial_rand_input_count
    + index->stats->query_filter_rand_input_count;
  index->stats->query_total_rand_output_count   = partial_rand_output_count
    + index->stats->query_filter_rand_output_count;
    
  index->stats->query_exact_distance = sqrtf(bsf_result.distance);
  index->stats->query_exact_node_filename = bsf_result.node->filename;
  index->stats->query_exact_node_size = bsf_result.node->node_size;;
  index->stats->query_exact_node_level = bsf_result.node->level;

  index->stats->query_refine_loaded_nodes_count = loaded_nodes_count;
  index->stats->query_refine_loaded_ts_count = loaded_ts_count;
  index->stats->query_refine_checked_nodes_count = checked_nodes_count;
  index->stats->query_refine_checked_ts_count = checked_ts_count;
  
  index->stats->query_total_loaded_nodes_count = loaded_nodes_count 
    + index->stats->query_filter_loaded_nodes_count;
  index->stats->query_total_loaded_ts_count = loaded_ts_count
    + index->stats->query_filter_loaded_ts_count;
  index->stats->query_total_checked_nodes_count = checked_nodes_count
    + index->stats->query_filter_checked_nodes_count;
  index->stats->query_total_checked_ts_count = checked_ts_count
    + index->stats->query_filter_checked_ts_count;


  index->stats->queries_refine_total_time  += index->stats->query_refine_total_time;	
    
  index->stats->queries_refine_input_time  += partial_input_time;
  index->stats->queries_refine_output_time += partial_output_time;
  index->stats->queries_refine_load_node_time += partial_load_node_time;        
  index->stats->queries_refine_cpu_time    += partial_time
    - partial_input_time
    - partial_output_time;
  index->stats->queries_refine_seq_input_count   += partial_seq_input_count;
  index->stats->queries_refine_seq_output_count  += partial_seq_output_count;
  index->stats->queries_refine_rand_input_count  += partial_rand_input_count;
  index->stats->queries_refine_rand_output_count += partial_rand_output_count;


  index->stats->queries_total_input_time  = index->stats->queries_refine_input_time
    +index->stats->queries_filter_input_time;	
  index->stats->queries_total_output_time  = index->stats->queries_refine_output_time
    +index->stats->queries_filter_output_time;	
  index->stats->queries_total_load_node_time  = index->stats->queries_refine_load_node_time
    +index->stats->queries_filter_load_node_time;
  index->stats->queries_total_cpu_time  = index->stats->queries_refine_cpu_time
    +index->stats->queries_filter_cpu_time;
    
  index->stats->queries_total_time  = index->stats->queries_refine_total_time
    +index->stats->queries_filter_total_time;
    
  index->stats->queries_total_seq_input_count   = index->stats->queries_filter_seq_input_count
    + index->stats->queries_refine_seq_input_count;
  index->stats->queries_total_seq_output_count  = index->stats->queries_filter_seq_output_count
    + index->stats->queries_refine_seq_output_count;
  index->stats->queries_total_rand_input_count  = index->stats->queries_filter_rand_input_count
    + index->stats->queries_refine_rand_input_count;
  index->stats->queries_total_rand_output_count = index->stats->queries_filter_rand_output_count
    + index->stats->queries_refine_rand_output_count;

       
  //keep a running sum then divide by the total number of queries
  index->stats->queries_avg_checked_nodes_count += index->stats->query_total_checked_nodes_count;
  index->stats->queries_avg_checked_ts_count += index->stats->query_total_checked_ts_count;
  index->stats->queries_avg_loaded_nodes_count += index->stats->query_total_loaded_nodes_count;
  index->stats->queries_avg_loaded_ts_count += index->stats->query_total_loaded_ts_count; 

    
  //COUNT_TOTAL_TIME_START
  return bsf_result;
}

void  exact_de_knn_search_psq (ts_type *query_ts, ts_type * query_ts_reordered,
			       int * query_order, unsigned int offset,
			       struct hercules_index *index,ts_type minimum_distance,
			       ts_type epsilon, double delta,
			       unsigned int k, unsigned int q_id, char * qfilename,
			       query_settings q_settings)
{
      
  unsigned int curr_size = 0;
  ts_type bsf = FLT_MAX;
  ts_type kth_bsf = FLT_MAX;
  ts_type temp_bsf = FLT_MAX;
  int i = 0; 
  //the next NN found by incremental search
  unsigned int found_knn = 0;
  int node_counter = 0;
    
  //queue containing kNN results
  struct query_result * knn_results = calloc(k,sizeof(struct query_result));
  for (int idx = 0; idx < k; ++idx)
    {
      knn_results[idx].node = NULL;
      knn_results[idx].distance = FLT_MAX;      
    }
    
  //return k approximate results
  approximate_knn_search(query_ts, query_ts_reordered, query_order, offset, bsf,
			 index,knn_results,k, NULL, NULL,&curr_size, q_settings.serial);
    
    
  //set the approximate result to be the first item in the queue
  struct query_result approximate_result = knn_results[0];    
  //struct query_result bsf_result = approximate_result;    

  COUNT_PARTIAL_TIME_END                 
    index->stats->query_filter_total_time  = partial_time;	
    
  index->stats->query_filter_input_time  = partial_input_time;
  index->stats->query_filter_output_time = partial_output_time;
  index->stats->query_filter_load_node_time = partial_load_node_time;
  index->stats->query_filter_cpu_time    = partial_time-partial_input_time-partial_output_time;
  index->stats->query_filter_seq_input_count   = partial_seq_input_count;
  index->stats->query_filter_seq_output_count  = partial_seq_output_count;
  index->stats->query_filter_rand_input_count  = partial_rand_input_count;
  index->stats->query_filter_rand_output_count = partial_rand_output_count;
    
  index->stats->query_filter_loaded_nodes_count = loaded_nodes_count;
  index->stats->query_filter_loaded_ts_count = loaded_ts_count;
  index->stats->query_filter_checked_nodes_count = checked_nodes_count;
  index->stats->query_filter_checked_ts_count = checked_ts_count;

  index->stats->query_approx_distance = sqrtf(approximate_result.distance);
  index->stats->query_approx_node_filename = approximate_result.node->filename;
  index->stats->query_approx_node_size = approximate_result.node->node_size;;
  index->stats->query_approx_node_level = approximate_result.node->level;
    
  index->stats->queries_filter_total_time  +=  index->stats->query_filter_total_time;
    
  index->stats->queries_filter_input_time  +=  index->stats->query_filter_input_time;
  index->stats->queries_filter_output_time +=  index->stats->query_filter_output_time;
  index->stats->queries_filter_load_node_time += index->stats->query_filter_load_node_time;
  index->stats->queries_filter_cpu_time    += index->stats->query_filter_cpu_time;

  index->stats->queries_filter_seq_input_count   += index->stats->query_filter_seq_input_count;
  index->stats->queries_filter_seq_output_count  += index->stats->query_filter_seq_output_count;
  index->stats->queries_filter_rand_input_count  += index->stats->query_filter_rand_input_count;
  index->stats->queries_filter_rand_output_count += index->stats->query_filter_rand_output_count;
    
    
  if(approximate_result.node != NULL) {  
    index->stats->query_approx_distance = sqrtf(approximate_result.distance);
    index->stats->query_approx_node_filename = approximate_result.node->filename;
    index->stats->query_approx_node_size = approximate_result.node->node_size;;
    index->stats->query_approx_node_level = approximate_result.node->level;
  }
 
  // Early termination...
  if (approximate_result.distance == 0) {
    index->stats->query_exact_distance = sqrtf(approximate_result.distance);
    index->stats->query_exact_node_filename = approximate_result.node->filename;
    index->stats->query_exact_node_size = approximate_result.node->node_size;;
    index->stats->query_exact_node_level = approximate_result.node->level; 
    update_query_stats(index,q_id, found_knn, approximate_result);
    get_query_stats(index, found_knn);
    print_query_stats(index, q_id, found_knn,qfilename);	     
  }

    
  RESET_QUERY_COUNTERS()
    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
 
    struct query_result bsf_result = approximate_result;
    
  pqueue_t *initial_queue = pqueue_init(index->first_node->node_size, 
					cmp_pri, get_pri, set_pri, get_pos, set_pos);

  pqueue_t *candidate_list = pqueue_init(index->first_node->node_size, 
					 cmp_pri, get_pri, set_pri, get_pos, set_pos);

  //set up the threads
    
  //creating maxquerythread threads
  pthread_t threadid[q_settings.num_threads];
  struct siss_query_worker_data qwdata;

  //exclusive locks
  pthread_mutex_t lock_queue=PTHREAD_MUTEX_INITIALIZER;//,lock_current_root_node=PTHREAD_MUTEX_INITIALIZER;

  //reader/writer lock
  pthread_rwlock_t lock_bsf=PTHREAD_RWLOCK_INITIALIZER;

  //once all threads join barrier, the barrier unblock all threads
  pthread_barrier_t lock_barrier;
  pthread_barrier_init(&lock_barrier, NULL, q_settings.num_threads);
    
  qwdata.query_ts=query_ts;
  qwdata.query_ts_reordered=query_ts_reordered;
  qwdata.query_order=query_order;
  qwdata.offset=offset;
  qwdata.minimum_distance=minimum_distance;
  qwdata.delta=delta;
  qwdata.epsilon=epsilon;
  qwdata.pq=candidate_list;
  qwdata.node_counter = &node_counter;    
  qwdata.lockvalueq=false;
  qwdata.lock_queue=&lock_queue;
  qwdata.lock_bsf=&lock_bsf;
  qwdata.index=index;
  //qwdata.bsf_result = &bsf_result;
  qwdata.knn_results =  knn_results;
  qwdata.k=k;    
  qwdata.lock_barrier=&lock_barrier;
  qwdata.curr_k_size = &curr_size;
 
  //Add the root to the priority queue
  struct query_result * root_pq_item = malloc(sizeof(struct query_result));
  root_pq_item->node = index->first_node;
  root_pq_item->distance = calculate_node_min_distance (index, index->first_node, query_ts);
  pqueue_insert(initial_queue, root_pq_item);
    

  struct query_result * n;
  struct query_result temp;

  ts_type bsf_stop = -1;
  ts_type r_delta_stop = -1;
  struct hercules_node * leaf_stop = NULL;

  int number_nodes = 1;
    
  //At this point, the queue contains only the approximate answer, not the root
  exact_de_knn_search_populate_queue(query_ts,
				     index,
				     knn_results,
				     k,
				     initial_queue,
				     delta,
				     epsilon,
				     candidate_list,
				     q_settings.num_threads);

  
  //put pq nodes in an array to process in parallel
  int num_nodes = pqueue_size(initial_queue);
  struct hercules_node **node_list = calloc(num_nodes,sizeof(struct hercules_node *));

  for (i = 0 ; i < num_nodes; ++i)
    {
      n = pqueue_pop(initial_queue);
      node_list[i] = n->node;
      free(n);
    }
    
  qwdata.node_list_size = num_nodes;    
  qwdata.node_list = node_list;

  //create maxquerythread threads, each will call the function
  //exact_search_worker_inmemory and pass to it rfdata
  for (i = 0; i < q_settings.num_threads; i++)
    {
      pthread_create(&(threadid[i]),NULL,exact_de_knn_search_worker,(void*)&(qwdata));
    }


  //block until all threads have completed execution
  for (i = 0; i < q_settings.num_threads; i++)
    {
      pthread_join(threadid[i],NULL);
    }

  COUNT_PARTIAL_TIME_END
    //printf ("Finished_Search = %lf\n", partial_time);  
    COUNT_PARTIAL_TIME_START

    // Free the nodes that where not popped.
    while ((n = pqueue_pop(initial_queue)))
      {
	free(n);
      }
  // Free the priority queue.
  pthread_barrier_destroy(&lock_barrier);
  pqueue_free(initial_queue);

  while ((n = pqueue_pop(candidate_list)))
    {
      free(n);
    }
  pqueue_free(candidate_list);
    
  //report the elements that were not reported already
  for (unsigned int pos = found_knn; pos < k; ++pos)
    {
      bsf_result = knn_results[pos];
      found_knn = pos+1;
      COUNT_PARTIAL_TIME_END
        update_query_stats(index,q_id, found_knn, bsf_result);
      get_query_stats(index, found_knn);
      print_query_stats(index, q_id, found_knn,qfilename);	
      //report all results for found_knn - last_found_knn or print their results
      RESET_QUERY_COUNTERS()
        RESET_PARTIAL_COUNTERS()
        COUNT_PARTIAL_TIME_START	
	}
  free(knn_results);	
  free(node_list);


}

void  exact_de_knn_search_psq_lf_sims_adaptive (ts_type *query_ts, ts_type * query_ts_reordered,
						int * query_order, unsigned int offset, ts_type *query_paa,
						struct hercules_index *index,ts_type minimum_distance,
						ts_type epsilon, double delta,
						unsigned int k, unsigned int q_id, char * qfilename,
						query_settings q_settings)
{
  printf ("Adaptive with EAPCA_THRESHOLD = %g and SAX_THRESHOLD = %g\n", q_settings.eapca_threshold, q_settings.sax_threshold);
      
  unsigned int curr_size = 0;
  ts_type bsf = FLT_MAX;
  ts_type kth_bsf = FLT_MAX;
  ts_type temp_bsf = FLT_MAX;
  int i = 0; 
  //the next NN found by incremental search
  unsigned int found_knn = 0;
  int node_counter = 0;
  unsigned int candidates_series_count = 0;
  unsigned int candidates_count = 0;
  unsigned int candidates_idx = 0;
  int sum_of_lab = 0;
  //queue containing kNN results
  struct query_result * knn_results = NULL;
  knn_results = calloc(k,sizeof(struct query_result));
  for (int idx = 0; idx < k; ++idx)
    {
      knn_results[idx].node = NULL;
      knn_results[idx].distance = FLT_MAX;      
    }
    
  pqueue_t *initial_queue = pqueue_init(index->first_node->node_size, 
					cmp_pri, get_pri, set_pri, get_pos, set_pos);

  struct query_result * root_pq_item = malloc(sizeof(struct query_result));
  root_pq_item->node = index->first_node;
  root_pq_item->distance = calculate_node_min_distance (index, index->first_node, query_ts);
  pqueue_insert(initial_queue, root_pq_item);


  ts_type * ts_list  =  malloc(index->settings->max_leaf_size * index->settings->timeseries_size * sizeof(ts_type));

  const char *leaves_filename = malloc(sizeof(char) * (strlen(index->settings->root_directory) + 18));
  leaves_filename = strcpy(leaves_filename, index->settings->root_directory);
  leaves_filename = strcat(leaves_filename, "leaves_raw.idx\0");
    
  FILE *ifile= NULL;
  ifile = fopen(leaves_filename,"rb");

  approximate_greedy_knn_search(query_ts, query_ts_reordered, query_order, offset, bsf,
				index,knn_results,k, NULL, NULL,&curr_size, q_settings.serial,
				initial_queue, ifile, ts_list, epsilon, delta, q_settings.approx_stop_condition);



  //set the approximate result to be the first item in the queue
  struct query_result approximate_result = knn_results[0];    
  //struct query_result bsf_result = approximate_result;    

  COUNT_PARTIAL_TIME_END                 
    index->stats->query_filter_total_time  = partial_time;	
    
  index->stats->query_filter_input_time  = partial_input_time;
  index->stats->query_filter_output_time = partial_output_time;
  index->stats->query_filter_load_node_time = partial_load_node_time;
  index->stats->query_filter_cpu_time    = partial_time-partial_input_time-partial_output_time;
  index->stats->query_filter_seq_input_count   = partial_seq_input_count;
  index->stats->query_filter_seq_output_count  = partial_seq_output_count;
  index->stats->query_filter_rand_input_count  = partial_rand_input_count;
  index->stats->query_filter_rand_output_count = partial_rand_output_count;
    
  index->stats->query_filter_loaded_nodes_count = loaded_nodes_count;
  index->stats->query_filter_loaded_ts_count = loaded_ts_count;
  index->stats->query_filter_checked_nodes_count = checked_nodes_count;
  index->stats->query_filter_checked_ts_count = checked_ts_count;

  index->stats->query_approx_distance = sqrtf(approximate_result.distance);
  index->stats->query_approx_node_filename = approximate_result.node->filename;
  index->stats->query_approx_node_size = approximate_result.node->node_size;;
  index->stats->query_approx_node_level = approximate_result.node->level;
    
  index->stats->queries_filter_total_time  +=  index->stats->query_filter_total_time;
    
  index->stats->queries_filter_input_time  +=  index->stats->query_filter_input_time;
  index->stats->queries_filter_output_time +=  index->stats->query_filter_output_time;
  index->stats->queries_filter_load_node_time += index->stats->query_filter_load_node_time;
  index->stats->queries_filter_cpu_time    += index->stats->query_filter_cpu_time;

  index->stats->queries_filter_seq_input_count   += index->stats->query_filter_seq_input_count;
  index->stats->queries_filter_seq_output_count  += index->stats->query_filter_seq_output_count;
  index->stats->queries_filter_rand_input_count  += index->stats->query_filter_rand_input_count;
  index->stats->queries_filter_rand_output_count += index->stats->query_filter_rand_output_count;
    
    
  if(approximate_result.node != NULL) {  
    index->stats->query_approx_distance = sqrtf(approximate_result.distance);
    index->stats->query_approx_node_filename = approximate_result.node->filename;
    index->stats->query_approx_node_size = approximate_result.node->node_size;;
    index->stats->query_approx_node_level = approximate_result.node->level;
  }
 
  // Early termination...
  if (approximate_result.distance == 0) {
    index->stats->query_exact_distance = sqrtf(approximate_result.distance);
    index->stats->query_exact_node_filename = approximate_result.node->filename;
    index->stats->query_exact_node_size = approximate_result.node->node_size;;
    index->stats->query_exact_node_level = approximate_result.node->level; 
    update_query_stats(index,q_id, found_knn, approximate_result);
    get_query_stats(index, found_knn);
    print_query_stats(index, q_id, found_knn,qfilename);	     
  }

    
  RESET_QUERY_COUNTERS()
    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
 
    struct query_result bsf_result = approximate_result;
    
  struct query_result * candidates =  calloc(index->stats->leaves_counter,
					     sizeof(struct query_result));
  int non_pruned_count = 0;
    

  struct query_result * n;
  struct query_result temp;

  ts_type bsf_stop = -1;
  ts_type r_delta_stop = -1;
  struct hercules_node * leaf_stop = NULL;

  int number_nodes = 1;
  //start sequentially, then create threads once queue has enough nodes
  //n = pqueue_pop(pq);
    
  //At this point, we pick the search from where the approximate search ended

  while ((n = pqueue_pop(initial_queue)))
    {
      kth_bsf =  knn_results[k-1].distance;
      if (n->distance > kth_bsf/(1 + epsilon))
    	{
          break;
        }

      if (n->node->is_leaf) // n is a leaf
        {
	  candidates[candidates_count].node = n->node;
	  candidates[candidates_count].distance = n->distance;
          ++candidates_count;
        }
      else  //n is an internal node
        {
          temp = knn_results[k-1];
          kth_bsf =  temp.distance;
	  
          ts_type child_distance;
	  child_distance = calculate_node_min_distance(index, n->node->left_child,query_ts);
	  
      
	  if ((child_distance < kth_bsf/(1+epsilon)) &&
	      (n->node->left_child != approximate_result.node)) 
	    {

	      struct query_result * mindist_result_left = malloc(sizeof(struct query_result));
	      mindist_result_left->node = n->node->left_child;
	      mindist_result_left->distance =  child_distance;
	      pqueue_insert(initial_queue, mindist_result_left);
	    }

	  child_distance = calculate_node_min_distance(index, n->node->right_child,query_ts);

	  if ((child_distance < kth_bsf/(1+epsilon))  &&
	      (n->node->right_child != approximate_result.node))	  
	    {
	      struct query_result * mindist_result_right = malloc(sizeof(struct query_result));
	      mindist_result_right->node = n->node->right_child;
	      mindist_result_right->distance =  child_distance;
	      pqueue_insert(initial_queue, mindist_result_right);
	    }
        }
      free(n);
    }
  while ((n = pqueue_pop(initial_queue)))
    {
      free(n);
    }
  pqueue_free(initial_queue);

  qsort(candidates,
	candidates_count,
	sizeof(struct query_result),
	compare_leaf_pos);


  ts_type total_leaves = index->stats->leaves_counter * 1.0;
  ts_type total_series = index->first_node->node_size * 1.0;
  float first_level_pruning = 1 - (candidates_count*1.0/total_leaves);
  printf("1st LEVEL FILTERING NODES = %u, CANDIDATE_LEAVES = %g, PRUNING=%g\n", candidates_count,total_leaves,first_level_pruning);
   
  if(first_level_pruning != 1)
    {

      if (first_level_pruning < q_settings.eapca_threshold)
	{

	  printf("Continue sequentially at idx = %g\n", total_leaves-candidates_count);
       
	  exact_de_knn_search_noapprox(query_ts, query_ts_reordered, query_order, offset,
				       index, minimum_distance, epsilon, delta, knn_results,k,
				       q_id, qfilename, q_settings,
				       candidates,candidates_count, &candidates_idx, 
				       ts_list, ifile);
     
	}
      else
	{
	  pthread_t threadid[q_settings.num_threads];
	  struct siss_query_worker_data qwdata[q_settings.num_threads];

	  pthread_mutex_t lock_queue=PTHREAD_MUTEX_INITIALIZER;//,lock_current_root_node=PTHREAD_MUTEX_INITIALIZER;

	  pthread_rwlock_t lock_bsf=PTHREAD_RWLOCK_INITIALIZER;

	  pthread_barrier_t lock_barrier;
	  pthread_barrier_init(&lock_barrier, NULL, q_settings.num_threads);
	  int num_nodes = pqueue_size(initial_queue);

	  for (int i = 0; i < q_settings.num_threads; i++)
	    {
	      qwdata[i].query_ts=query_ts;
	      qwdata[i].query_ts_reordered=query_ts_reordered;
	      qwdata[i].query_order=query_order;
	      qwdata[i].offset=offset;
	      qwdata[i].query_paa=query_paa;      
	      qwdata[i].minimum_distance=minimum_distance;
	      qwdata[i].delta=delta;
	      qwdata[i].epsilon=epsilon;
	      qwdata[i].candidates=candidates;
	      qwdata[i].candidates_count=&candidates_count;
	      qwdata[i].candidates_idx=&candidates_idx;
	      qwdata[i].node_counter = &node_counter;    
	      qwdata[i].lockvalueq=false;
	      qwdata[i].lock_queue=&lock_queue;
	      qwdata[i].lock_bsf=&lock_bsf;
	      qwdata[i].index=index;
	      qwdata[i].knn_results =  knn_results;
	      qwdata[i].k=k;    
	      qwdata[i].lock_barrier=&lock_barrier;
	      qwdata[i].curr_k_size = &curr_size;

	      qwdata[i].sum_of_lab = 0;
	      qwdata[i].candidates_series_count=&candidates_series_count;
	    }    
	  for (i = 0; i < q_settings.num_threads; i++)
	    {
	      pthread_create(&(threadid[i]),NULL,exact_de_knn_mindist_worker_lf,(void*)&(qwdata[i]));
	    }

	  for (i = 0; i < q_settings.num_threads; i++)
	    {
	      pthread_join(threadid[i],NULL);
	      sum_of_lab+=qwdata[i].sum_of_lab;
	    }
	  float second_level_pruning =  1- (sum_of_lab*1.0/total_series);
	  printf("2nd LEVEL FILTERING CANDIDATES = %u, PRUNING=%.6g\n", sum_of_lab, second_level_pruning);

	  if (second_level_pruning < q_settings.sax_threshold)
	    {
	      printf("Continue sequentially at idx = %g\n", total_leaves-candidates_count);	     
	      exact_de_knn_search_noapprox(query_ts, query_ts_reordered, query_order, offset,
					   index, minimum_distance, epsilon, delta, knn_results,k,
					   q_id, qfilename, q_settings,
					   candidates,candidates_count, &candidates_idx, 
					   ts_list, ifile);
	    }
	  else
	    {
	      unsigned long * label_number=malloc(sizeof(unsigned long )*(sum_of_lab));
	      float* minidisvector=malloc(sizeof(float)*(sum_of_lab));
   
	      sum_of_lab=0;

	      for (i = 0; i < q_settings.num_threads; i++)
		{
		  memcpy(&(label_number[sum_of_lab]),qwdata[i].label_number,sizeof(unsigned long)*qwdata[i].sum_of_lab);
		  memcpy(&(minidisvector[sum_of_lab]),qwdata[i].minidisvector,sizeof(float)*qwdata[i].sum_of_lab);
		  free(qwdata[i].label_number);
		  free(qwdata[i].minidisvector);
		  sum_of_lab+=qwdata[i].sum_of_lab;
		}
    
	      unsigned int total_read_threads =  q_settings.num_threads * q_settings.num_read_threads;
	      pthread_t readthreadid[total_read_threads];
    
	      struct siss_query_worker_data read_qwdata[total_read_threads];

	      candidates_idx = 0;
	      for (int i = 0; i < total_read_threads; i++)
		{
		  read_qwdata[i].query_ts=query_ts;
		  read_qwdata[i].delta=delta;
		  read_qwdata[i].epsilon=epsilon;
		  read_qwdata[i].candidates=candidates;
		  read_qwdata[i].candidates_count=&candidates_count;
		  read_qwdata[i].candidates_idx=&candidates_idx;
		  read_qwdata[i].lockvalueq=false;
		  read_qwdata[i].lock_queue=&lock_queue;
		  read_qwdata[i].lock_bsf=&lock_bsf;
		  read_qwdata[i].index=index;
		  read_qwdata[i].knn_results = knn_results;
		  read_qwdata[i].k=k;    
		  read_qwdata[i].curr_k_size = &curr_size;
		  read_qwdata[i].thread_input_time = 0;
		  read_qwdata[i].thread_insert_node_time = 0;
		  read_qwdata[i].thread_realdist_time = 0;		
		  read_qwdata[i].thread_random_io = 0;
		  read_qwdata[i].thread_sequential_io = 0;
    
		  read_qwdata[i].minidisvector=minidisvector;
		  read_qwdata[i].sum_of_lab=sum_of_lab;
		  read_qwdata[i].load_point=label_number;
		}    

	      for (i = 0; i < total_read_threads; i++)
		{
		  pthread_create(&(readthreadid[i]),NULL,exact_de_knn_read_worker_lf_sims,(void*)&(read_qwdata[i]));
		}            
    
	      double temp_time = 0;
	      unsigned long temp_rio = 0;
	      unsigned long temp_sio = 0;
 
	      for (i = 0; i < total_read_threads; i++)
		{
		  pthread_join(readthreadid[i],NULL);
#if DETAILED_STATS == 1
		  temp_time = fmax(temp_time, read_qwdata[i].thread_input_time);
		  temp_rio = temp_rio + read_qwdata[i].thread_random_io ;
		  temp_sio = temp_sio + read_qwdata[i].thread_sequential_io ;
#endif
		}
	
	      COUNT_PARTIAL_TIME_END
		COUNT_PARTIAL_TIME_START

#if DETAILED_STATS == 1
		partial_input_time = temp_time;
	      loaded_ts_count = temp_rio;
	      loaded_nodes_count = temp_rio;
#endif
	      free(label_number);
	      free(minidisvector);
	    }
	}
    }

  for (unsigned int pos = found_knn; pos < k; ++pos)
    {
      bsf_result = knn_results[pos];
      found_knn = pos+1;
      COUNT_PARTIAL_TIME_END
	update_query_stats(index,q_id, found_knn, bsf_result);
      get_query_stats(index, found_knn);
      print_query_stats(index, q_id, found_knn,qfilename);	
      RESET_QUERY_COUNTERS()
	RESET_PARTIAL_COUNTERS()
	COUNT_PARTIAL_TIME_START	
	}

  free(candidates);
  free(knn_results);	
  fflush(stdout);

  free(ts_list);
  free(leaves_filename);
  fclose(ifile);

}

void  exact_de_knn_search_psq_lf_sims (ts_type *query_ts, ts_type * query_ts_reordered,
				       int * query_order, unsigned int offset, ts_type *query_paa,
				       struct hercules_index *index,ts_type minimum_distance,
				       ts_type epsilon, double delta,
				       unsigned int k, unsigned int q_id, char * qfilename,
				       query_settings q_settings)
{
      
  printf ("Non adaptive\n");
  unsigned int curr_size = 0;
  ts_type bsf = FLT_MAX;
  ts_type kth_bsf = FLT_MAX;
  ts_type temp_bsf = FLT_MAX;
  int i = 0; 
  //the next NN found by incremental search
  unsigned int found_knn = 0;
  int node_counter = 0;
  unsigned int candidates_series_count = 0;
  unsigned int candidates_count = 0;
  unsigned int candidates_idx = 0;
  int sum_of_lab = 0;
  //queue containing kNN results
  struct query_result * knn_results = NULL;
  knn_results = calloc(k,sizeof(struct query_result));
  for (int idx = 0; idx < k; ++idx)
    {
      knn_results[idx].node = NULL;
      knn_results[idx].distance = FLT_MAX;      
    }
    
  //return k approximate results
  approximate_knn_search(query_ts, query_ts_reordered, query_order, offset, bsf,
			 index,knn_results,k, NULL, NULL,&curr_size, q_settings.serial);
    
    
  //set the approximate result to be the first item in the queue
  struct query_result approximate_result = knn_results[0];    

  COUNT_PARTIAL_TIME_END                 
    index->stats->query_filter_total_time  = partial_time;	
    
  index->stats->query_filter_input_time  = partial_input_time;
  index->stats->query_filter_output_time = partial_output_time;
  index->stats->query_filter_load_node_time = partial_load_node_time;
  index->stats->query_filter_cpu_time    = partial_time-partial_input_time-partial_output_time;
  index->stats->query_filter_seq_input_count   = partial_seq_input_count;
  index->stats->query_filter_seq_output_count  = partial_seq_output_count;
  index->stats->query_filter_rand_input_count  = partial_rand_input_count;
  index->stats->query_filter_rand_output_count = partial_rand_output_count;
    
  index->stats->query_filter_loaded_nodes_count = loaded_nodes_count;
  index->stats->query_filter_loaded_ts_count = loaded_ts_count;
  index->stats->query_filter_checked_nodes_count = checked_nodes_count;
  index->stats->query_filter_checked_ts_count = checked_ts_count;

  index->stats->query_approx_distance = sqrtf(approximate_result.distance);
  index->stats->query_approx_node_filename = approximate_result.node->filename;
  index->stats->query_approx_node_size = approximate_result.node->node_size;;
  index->stats->query_approx_node_level = approximate_result.node->level;
    
  index->stats->queries_filter_total_time  +=  index->stats->query_filter_total_time;
    
  index->stats->queries_filter_input_time  +=  index->stats->query_filter_input_time;
  index->stats->queries_filter_output_time +=  index->stats->query_filter_output_time;
  index->stats->queries_filter_load_node_time += index->stats->query_filter_load_node_time;
  index->stats->queries_filter_cpu_time    += index->stats->query_filter_cpu_time;

  index->stats->queries_filter_seq_input_count   += index->stats->query_filter_seq_input_count;
  index->stats->queries_filter_seq_output_count  += index->stats->query_filter_seq_output_count;
  index->stats->queries_filter_rand_input_count  += index->stats->query_filter_rand_input_count;
  index->stats->queries_filter_rand_output_count += index->stats->query_filter_rand_output_count;
    
    
  if(approximate_result.node != NULL) {  
    index->stats->query_approx_distance = sqrtf(approximate_result.distance);
    index->stats->query_approx_node_filename = approximate_result.node->filename;
    index->stats->query_approx_node_size = approximate_result.node->node_size;;
    index->stats->query_approx_node_level = approximate_result.node->level;
  }
 
  // Early termination...
  if (approximate_result.distance == 0) {
    index->stats->query_exact_distance = sqrtf(approximate_result.distance);
    index->stats->query_exact_node_filename = approximate_result.node->filename;
    index->stats->query_exact_node_size = approximate_result.node->node_size;;
    index->stats->query_exact_node_level = approximate_result.node->level; 
    update_query_stats(index,q_id, found_knn, approximate_result);
    get_query_stats(index, found_knn);
    print_query_stats(index, q_id, found_knn,qfilename);	     
  }

    
  RESET_QUERY_COUNTERS()
    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
 
    struct query_result bsf_result = approximate_result;
    
  pqueue_t *initial_queue = pqueue_init(index->first_node->node_size, 
					cmp_pri, get_pri, set_pri, get_pos, set_pos);

  struct query_result * candidates =  calloc(index->stats->leaves_counter,
					     sizeof(struct query_result));
  int non_pruned_count = 0;
    

  pthread_t threadid[q_settings.num_threads];
  struct siss_query_worker_data qwdata[q_settings.num_threads];

  pthread_mutex_t lock_queue=PTHREAD_MUTEX_INITIALIZER;//,lock_current_root_node=PTHREAD_MUTEX_INITIALIZER;

  pthread_rwlock_t lock_bsf=PTHREAD_RWLOCK_INITIALIZER;

  pthread_barrier_t lock_barrier;
  pthread_barrier_init(&lock_barrier, NULL, q_settings.num_threads);

   
  //Add the root to the priority queue
  struct query_result * root_pq_item = malloc(sizeof(struct query_result));
  root_pq_item->node = index->first_node;
  root_pq_item->distance = calculate_node_min_distance (index, index->first_node, query_ts);
  pqueue_insert(initial_queue, root_pq_item);
    

  struct query_result * n;
  struct query_result temp;

  ts_type bsf_stop = -1;
  ts_type r_delta_stop = -1;
  struct hercules_node * leaf_stop = NULL;

  int number_nodes = 1;


  exact_de_knn_search_populate_queue_lf(query_ts,
					index,
					knn_results,
					k,
					initial_queue,
					delta,
					epsilon,
					candidates,
					&candidates_count,
					q_settings.num_threads);
    
  
  //put pq nodes in an array to process in parallel
  int num_nodes = pqueue_size(initial_queue);
  struct hercules_node **node_list= calloc(num_nodes,sizeof(struct hercules_node *));

  for (i = 0 ; i < num_nodes; ++i)
    {
      n = pqueue_pop(initial_queue);
      node_list[i] = n->node;
      free(n);
    }


  for (int i = 0; i < q_settings.num_threads; i++)
    {
      qwdata[i].query_ts=query_ts;
      qwdata[i].query_ts_reordered=query_ts_reordered;
      qwdata[i].query_order=query_order;
      qwdata[i].offset=offset;
      qwdata[i].query_paa=query_paa;      
      qwdata[i].minimum_distance=minimum_distance;
      qwdata[i].delta=delta;
      qwdata[i].epsilon=epsilon;
      qwdata[i].candidates=candidates;
      qwdata[i].candidates_count=&candidates_count;
      qwdata[i].candidates_idx=&candidates_idx;
      qwdata[i].node_counter = &node_counter;    
      qwdata[i].lockvalueq=false;
      qwdata[i].lock_queue=&lock_queue;
      qwdata[i].lock_bsf=&lock_bsf;
      qwdata[i].index=index;
      qwdata[i].knn_results =  knn_results;
      qwdata[i].k=k;    
      qwdata[i].lock_barrier=&lock_barrier;
      qwdata[i].curr_k_size = &curr_size;
      qwdata[i].node_list_size = num_nodes;    
      qwdata[i].node_list = node_list;

      qwdata[i].sum_of_lab = 0;
      qwdata[i].candidates_series_count=&candidates_series_count;
    }    

    
  for (i = 0; i < q_settings.num_threads; i++)
    {
      pthread_create(&(threadid[i]),NULL,exact_de_knn_search_worker_lf,(void*)&(qwdata[i]));
    }

    
  for (i = 0; i < q_settings.num_threads; i++)
    {
      pthread_join(threadid[i],NULL);
    }
 

  qsort(candidates,
	candidates_count,
	sizeof(struct query_result),
	compare_leaf_pos);

  for (i = 0; i < q_settings.num_threads; i++)
    {
      pthread_create(&(threadid[i]),NULL,exact_de_knn_mindist_worker_lf,(void*)&(qwdata[i]));
    }

  for (i = 0; i < q_settings.num_threads; i++)
    {
      pthread_join(threadid[i],NULL);
      sum_of_lab+=qwdata[i].sum_of_lab;
    }
  unsigned long * label_number=malloc(sizeof(unsigned long )*(sum_of_lab));
  float* minidisvector=malloc(sizeof(float)*(sum_of_lab));
   
  sum_of_lab=0;

  for (i = 0; i < q_settings.num_threads; i++)
    {
      memcpy(&(label_number[sum_of_lab]),qwdata[i].label_number,sizeof(unsigned long)*qwdata[i].sum_of_lab);
      memcpy(&(minidisvector[sum_of_lab]),qwdata[i].minidisvector,sizeof(float)*qwdata[i].sum_of_lab);
      free(qwdata[i].label_number);
      free(qwdata[i].minidisvector);
      sum_of_lab+=qwdata[i].sum_of_lab;
    }
    
  unsigned int total_read_threads =  q_settings.num_threads * q_settings.num_read_threads;
  pthread_t readthreadid[total_read_threads];
    
  struct siss_query_worker_data read_qwdata[total_read_threads];

  candidates_idx = 0;
  for (int i = 0; i < total_read_threads; i++)
    {
      read_qwdata[i].query_ts=query_ts;
      read_qwdata[i].delta=delta;
      read_qwdata[i].epsilon=epsilon;
      read_qwdata[i].candidates=candidates;
      read_qwdata[i].candidates_count=&candidates_count;
      read_qwdata[i].candidates_idx=&candidates_idx;
      read_qwdata[i].lockvalueq=false;
      read_qwdata[i].lock_queue=&lock_queue;
      read_qwdata[i].lock_bsf=&lock_bsf;
      read_qwdata[i].index=index;
      read_qwdata[i].knn_results = knn_results;
      read_qwdata[i].k=k;    
      read_qwdata[i].curr_k_size = &curr_size;
      read_qwdata[i].thread_input_time = 0;
      read_qwdata[i].thread_insert_node_time = 0;
      read_qwdata[i].thread_realdist_time = 0;		
      read_qwdata[i].thread_random_io = 0;
      read_qwdata[i].thread_sequential_io = 0;
    
      read_qwdata[i].minidisvector=minidisvector;
      read_qwdata[i].sum_of_lab=sum_of_lab;
      read_qwdata[i].load_point=label_number;
    }    

  for (i = 0; i < total_read_threads; i++)
    {
      pthread_create(&(readthreadid[i]),NULL,exact_de_knn_read_worker_lf_sims,(void*)&(read_qwdata[i]));
    }            
    
  double temp_time = 0;
  unsigned long temp_rio = 0;
  unsigned long temp_sio = 0;
 
  for (i = 0; i < total_read_threads; i++)
    {
      pthread_join(readthreadid[i],NULL);
#if DETAILED_STATS == 1
      temp_time = fmax(temp_time, read_qwdata[i].thread_input_time);
      temp_rio = temp_rio + read_qwdata[i].thread_random_io ;
      temp_sio = temp_sio + read_qwdata[i].thread_sequential_io ;
#endif
    }
	
  COUNT_PARTIAL_TIME_END
    COUNT_PARTIAL_TIME_START



#if DETAILED_STATS == 1
    partial_input_time = temp_time;
  loaded_ts_count = temp_sio;
  loaded_nodes_count = temp_rio;
#endif

  for (unsigned int pos = found_knn; pos < k; ++pos)
    {
      bsf_result = knn_results[pos];
      found_knn = pos+1;
      COUNT_PARTIAL_TIME_END
        update_query_stats(index,q_id, found_knn, bsf_result);
      get_query_stats(index, found_knn);
      print_query_stats(index, q_id, found_knn,qfilename);	
      RESET_QUERY_COUNTERS()
        RESET_PARTIAL_COUNTERS()
        COUNT_PARTIAL_TIME_START	
	}

  free(label_number);
  free(minidisvector);

  free(candidates);
  free(knn_results);	
  free(node_list);

}

void  exact_de_knn_search_psq_lf (ts_type *query_ts, ts_type * query_ts_reordered,
				  int * query_order, unsigned int offset,
				  struct hercules_index *index,ts_type minimum_distance,
				  ts_type epsilon, double delta,
				  unsigned int k, unsigned int q_id, char * qfilename,
				  query_settings q_settings)
{
      
  unsigned int curr_size = 0;
  ts_type bsf = FLT_MAX;
  ts_type kth_bsf = FLT_MAX;
  ts_type temp_bsf = FLT_MAX;
  int i = 0; 
  unsigned int found_knn = 0;
  int node_counter = 0;
  unsigned int candidates_count = 0;
  unsigned int candidates_idx = 0;
    
  struct query_result * knn_results = calloc(k,sizeof(struct query_result));
  for (int idx = 0; idx < k; ++idx)
    {
      knn_results[idx].node = NULL;
      knn_results[idx].distance = FLT_MAX;      
    }
    
  approximate_knn_search(query_ts, query_ts_reordered, query_order, offset, bsf,
			 index,knn_results,k, NULL, NULL,&curr_size, q_settings.serial);
    
    
  struct query_result approximate_result = knn_results[0];    

  COUNT_PARTIAL_TIME_END                 
    index->stats->query_filter_total_time  = partial_time;	
    
  index->stats->query_filter_input_time  = partial_input_time;
  index->stats->query_filter_output_time = partial_output_time;
  index->stats->query_filter_load_node_time = partial_load_node_time;
  index->stats->query_filter_cpu_time    = partial_time-partial_input_time-partial_output_time;
  index->stats->query_filter_seq_input_count   = partial_seq_input_count;
  index->stats->query_filter_seq_output_count  = partial_seq_output_count;
  index->stats->query_filter_rand_input_count  = partial_rand_input_count;
  index->stats->query_filter_rand_output_count = partial_rand_output_count;
    
  index->stats->query_filter_loaded_nodes_count = loaded_nodes_count;
  index->stats->query_filter_loaded_ts_count = loaded_ts_count;
  index->stats->query_filter_checked_nodes_count = checked_nodes_count;
  index->stats->query_filter_checked_ts_count = checked_ts_count;

  index->stats->query_approx_distance = sqrtf(approximate_result.distance);
  index->stats->query_approx_node_filename = approximate_result.node->filename;
  index->stats->query_approx_node_size = approximate_result.node->node_size;;
  index->stats->query_approx_node_level = approximate_result.node->level;
    
  index->stats->queries_filter_total_time  +=  index->stats->query_filter_total_time;
    
  index->stats->queries_filter_input_time  +=  index->stats->query_filter_input_time;
  index->stats->queries_filter_output_time +=  index->stats->query_filter_output_time;
  index->stats->queries_filter_load_node_time += index->stats->query_filter_load_node_time;
  index->stats->queries_filter_cpu_time    += index->stats->query_filter_cpu_time;

  index->stats->queries_filter_seq_input_count   += index->stats->query_filter_seq_input_count;
  index->stats->queries_filter_seq_output_count  += index->stats->query_filter_seq_output_count;
  index->stats->queries_filter_rand_input_count  += index->stats->query_filter_rand_input_count;
  index->stats->queries_filter_rand_output_count += index->stats->query_filter_rand_output_count;
    
    
  if(approximate_result.node != NULL) {  
    index->stats->query_approx_distance = sqrtf(approximate_result.distance);
    index->stats->query_approx_node_filename = approximate_result.node->filename;
    index->stats->query_approx_node_size = approximate_result.node->node_size;;
    index->stats->query_approx_node_level = approximate_result.node->level;
  }
 
  if (approximate_result.distance == 0) {
    index->stats->query_exact_distance = sqrtf(approximate_result.distance);
    index->stats->query_exact_node_filename = approximate_result.node->filename;
    index->stats->query_exact_node_size = approximate_result.node->node_size;;
    index->stats->query_exact_node_level = approximate_result.node->level; 
    update_query_stats(index,q_id, found_knn, approximate_result);
    get_query_stats(index, found_knn);
    print_query_stats(index, q_id, found_knn,qfilename);	     
  }

    
  RESET_QUERY_COUNTERS()
    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
 
    struct query_result bsf_result = approximate_result;
    
  pqueue_t *initial_queue = pqueue_init(index->first_node->node_size, 
					cmp_pri, get_pri, set_pri, get_pos, set_pos);

  struct query_result * candidates =  calloc(index->stats->leaves_counter,
					     sizeof(struct query_result));
  int non_pruned_count = 0;
    

  pthread_t threadid[q_settings.num_threads];
  struct siss_query_worker_data qwdata[q_settings.num_threads];

  pthread_mutex_t lock_queue=PTHREAD_MUTEX_INITIALIZER;//,lock_current_root_node=PTHREAD_MUTEX_INITIALIZER;

  pthread_rwlock_t lock_bsf=PTHREAD_RWLOCK_INITIALIZER;

  pthread_barrier_t lock_barrier;
  pthread_barrier_init(&lock_barrier, NULL, q_settings.num_threads);

   
  struct query_result * root_pq_item = malloc(sizeof(struct query_result));
  root_pq_item->node = index->first_node;
  root_pq_item->distance = calculate_node_min_distance (index, index->first_node, query_ts);
  pqueue_insert(initial_queue, root_pq_item);
    

  struct query_result * n;
  struct query_result temp;

  ts_type bsf_stop = -1;
  ts_type r_delta_stop = -1;
  struct hercules_node * leaf_stop = NULL;

  int number_nodes = 1;


  exact_de_knn_search_populate_queue_lf(query_ts,
					index,
					knn_results,
					k,
					initial_queue,
					delta,
					epsilon,
					candidates,
					&candidates_count,
					q_settings.num_threads);
    
  
  int num_nodes = pqueue_size(initial_queue);
  struct hercules_node **node_list = calloc(num_nodes,sizeof(struct hercules_node *));

  for (i = 0 ; i < num_nodes; ++i)
    {
      n = pqueue_pop(initial_queue);
      node_list[i] = n->node;
      free(n);
    }


  for (int i = 0; i < q_settings.num_threads; i++)
    {
      qwdata[i].query_ts=query_ts;
      qwdata[i].query_ts_reordered=query_ts_reordered;
      qwdata[i].query_order=query_order;
      qwdata[i].offset=offset;
      qwdata[i].minimum_distance=minimum_distance;
      qwdata[i].delta=delta;
      qwdata[i].epsilon=epsilon;
      qwdata[i].candidates=candidates;
      qwdata[i].candidates_count=&candidates_count;
      qwdata[i].candidates_idx=&candidates_idx;
      qwdata[i].node_counter = &node_counter;    
      qwdata[i].lockvalueq=false;
      qwdata[i].lock_queue=&lock_queue;
      qwdata[i].lock_bsf=&lock_bsf;
      qwdata[i].index=index;
      qwdata[i].knn_results =  knn_results;
      qwdata[i].k=k;    
      qwdata[i].lock_barrier=&lock_barrier;
      qwdata[i].curr_k_size = &curr_size;
      qwdata[i].node_list_size = num_nodes;    
      qwdata[i].node_list = node_list;      
    }    

    
  for (i = 0; i < q_settings.num_threads; i++)
    {
      pthread_create(&(threadid[i]),NULL,exact_de_knn_search_worker_lf,(void*)&(qwdata[i]));
    }

  for (i = 0; i < q_settings.num_threads; i++)
    {
      pthread_join(threadid[i],NULL);
    }
 

  qsort(candidates,
	candidates_count,
	sizeof(struct query_result),
	compare_leaf_pos);

  unsigned int total_read_threads =  q_settings.num_threads * q_settings.num_read_threads;
  pthread_t readthreadid[total_read_threads];
  struct siss_query_worker_data read_qwdata[total_read_threads];


  for (int i = 0; i < total_read_threads; i++)
    {
      read_qwdata[i].query_ts=query_ts;
      read_qwdata[i].delta=delta;
      read_qwdata[i].epsilon=epsilon;
      read_qwdata[i].candidates=candidates;
      read_qwdata[i].candidates_count=&candidates_count;
      read_qwdata[i].candidates_idx=&candidates_idx;
      read_qwdata[i].lockvalueq=false;
      read_qwdata[i].lock_queue=&lock_queue;
      read_qwdata[i].lock_bsf=&lock_bsf;
      read_qwdata[i].index=index;
      read_qwdata[i].knn_results =  knn_results;
      read_qwdata[i].k=k;    
      read_qwdata[i].curr_k_size = &curr_size;
      read_qwdata[i].ts_list = malloc (index->settings->max_leaf_size *
				       index->settings->timeseries_size * sizeof(ts_type));

      read_qwdata[i].thread_input_time = 0;
      read_qwdata[i].thread_insert_node_time = 0;
      read_qwdata[i].thread_realdist_time = 0;		

    }    

  for (i = 0; i < total_read_threads; i++)
    {
      pthread_create(&(readthreadid[i]),NULL,exact_de_knn_read_worker_lf,(void*)&(read_qwdata[i]));
    }            



  //block until all threads have completed execution
  for (i = 0; i < total_read_threads; i++)
    {
      pthread_join(readthreadid[i],NULL);
    }

  double temp_time = 0;
    
#if DETAILED_STATS == 1
  for (i = 0; i < total_read_threads; i++)
    {
      threads_total_input_time[i] +=  read_qwdata[i].thread_input_time;
      threads_total_realdist_time[i] +=  read_qwdata[i].thread_realdist_time;
      threads_total_insert_node_time[i] +=  read_qwdata[i].thread_insert_node_time;            
      temp_time = fmax(temp_time, read_qwdata[i].thread_input_time);
    }    
#endif
	
  COUNT_PARTIAL_TIME_END
    COUNT_PARTIAL_TIME_START

    // Free the nodes that where not popped.
    while ((n = pqueue_pop(initial_queue)))
      {
	free(n);
      }
  // Free the priority queue.
  pthread_barrier_destroy(&lock_barrier);
  pqueue_free(initial_queue);

  free(candidates);
    
#if DETAILED_STATS == 1
  partial_input_time = temp_time;
#endif

  for (unsigned int pos = found_knn; pos < k; ++pos)
    {
      bsf_result = knn_results[pos];
      found_knn = pos+1;
      COUNT_PARTIAL_TIME_END
        update_query_stats(index,q_id, found_knn, bsf_result);
      get_query_stats(index, found_knn);
      print_query_stats(index, q_id, found_knn,qfilename);	
      RESET_QUERY_COUNTERS()
        RESET_PARTIAL_COUNTERS()
        COUNT_PARTIAL_TIME_START	
	}

  free(knn_results);	
  free(node_list);


  for (int i = 0; i < total_read_threads; i++)
    {
      free(read_qwdata[i].ts_list);
    }


}


void exact_de_knn_search_populate_queue_lf(ts_type *query_ts,
					   struct hercules_index *index,
					   struct query_result * knn_results,
					   int k,
					   pqueue_t *pq,
					   ts_type delta,
					   ts_type epsilon,
					   struct query_result * candidates,
					   unsigned int * candidates_count,
					   int stop_condition)
{

  struct query_result *n, temp;
  int number_nodes = 0;
  ts_type kth_bsf;
  struct hercules_node *kth_node;

  while ((n = pqueue_pop(pq)))
    {

      kth_bsf =  knn_results[k-1].distance;
      
      if (n->distance > kth_bsf/(1 + epsilon))
	{
          break;
	}
      
      if (n->node->is_leaf) // n is a leaf
	{
	  struct query_result cand_result; 
	  cand_result.node = n->node;
	  cand_result.distance=n->distance;
	  candidates[*candidates_count] = cand_result;
	  ++(*candidates_count); 	  
	}
      else  //n is an internal node
	{
	  temp = knn_results[k-1];
          kth_bsf =  temp.distance;
          kth_node = temp.node;
	    
          ts_type child_distance;
	  child_distance = calculate_node_min_distance(index, n->node->left_child,query_ts);
	  
	  if ((child_distance < kth_bsf/(1+epsilon)) &
	      n->node->left_child != kth_node) //add epsilon
	    {
	      struct query_result * mindist_result_left = malloc(sizeof(struct query_result));
	      mindist_result_left->node = n->node->left_child;
	      mindist_result_left->distance =  child_distance;
	      pqueue_insert(pq, mindist_result_left);
	    }

	  child_distance = calculate_node_min_distance(index, n->node->right_child,query_ts);

	  if ((child_distance < kth_bsf/(1+epsilon)) &
	      n->node->right_child != kth_node) 
	    {
	      struct query_result * mindist_result_right = malloc(sizeof(struct query_result));
	      mindist_result_right->node = n->node->right_child;
	      mindist_result_right->distance =  child_distance;
	      pqueue_insert(pq, mindist_result_right);
	    }

	  if (pqueue_size(pq) >= stop_condition)
	    {
	      free(n);
	      break;
	    }
	}
      free(n);
    }
}
void exact_de_knn_search_populate_queue(ts_type *query_ts,
					struct hercules_index *index,
					struct query_result * knn_results,
					int k,
					pqueue_t *pq,
					ts_type delta,
					ts_type epsilon,
					pqueue_t * candidate_list,
					int stop_condition)
{

  struct query_result *n, temp;
  int number_nodes = 0;
  ts_type kth_bsf;
  struct hercules_node *kth_node;
  
  while ((n = pqueue_pop(pq)))
    {

      kth_bsf =  knn_results[k-1].distance;
      
      if (n->distance > kth_bsf/(1 + epsilon))
	{
          break;
	}
      
      if (n->node->is_leaf) // n is a leaf
	{
	  struct query_result * mindist_result = malloc(sizeof(struct query_result));
	  mindist_result->node = n->node;
	  mindist_result->distance=n->distance;
	  pqueue_insert(candidate_list, mindist_result);
	}
      else  //n is an internal node
	{
	  temp = knn_results[k-1];
          kth_bsf =  temp.distance;
          kth_node = temp.node;
	    
          ts_type child_distance;
	  child_distance = calculate_node_min_distance(index, n->node->left_child,query_ts);
	  
	  if ((child_distance < kth_bsf/(1+epsilon)) &
	      n->node->left_child != kth_node) 
	    {
	      struct query_result * mindist_result_left = malloc(sizeof(struct query_result));
	      mindist_result_left->node = n->node->left_child;
	      mindist_result_left->distance =  child_distance;
	      pqueue_insert(pq, mindist_result_left);
	    }

	  child_distance = calculate_node_min_distance(index, n->node->right_child,query_ts);

	  if ((child_distance < kth_bsf/(1+epsilon)) &
	      n->node->right_child != kth_node) 
	    {
	      struct query_result * mindist_result_right = malloc(sizeof(struct query_result));
	      mindist_result_right->node = n->node->right_child;
	      mindist_result_right->distance =  child_distance;
	      pqueue_insert(pq, mindist_result_right);
	    }

	  if (pqueue_size(pq) >= stop_condition)
	    {
	      free(n);
	      break;
	    }
	}
      free(n);
    }
}

void* exact_de_knn_search_worker(void *qwdata)
{

  
  struct hercules_node *current_node;
  struct query_result *n;
  struct hercules_index *index=((siss_query_worker_data*)qwdata)->index;
  int *query_order=((siss_query_worker_data*)qwdata)->query_order;
  ts_type *query_ts_reordered=((siss_query_worker_data*)qwdata)->query_ts_reordered;
  ts_type *query_ts=((siss_query_worker_data*)qwdata)->query_ts;
  pqueue_t *pq=((siss_query_worker_data*)qwdata)->pq;
  //struct query_result *do_not_remove = ((siss_query_worker_data*)qwdata)->bsf_result;
  unsigned int offset=((siss_query_worker_data*)qwdata)->offset;
  ts_type minimum_distance=((siss_query_worker_data*)qwdata)->minimum_distance;
  ts_type epsilon=((siss_query_worker_data*)qwdata)->epsilon;
  ts_type delta=((siss_query_worker_data*)qwdata)->delta;    
  unsigned int k = ((siss_query_worker_data*)qwdata)->k;  
  struct query_result *knn_results = (((siss_query_worker_data*)qwdata)->knn_results);
  ts_type kth_bsf=knn_results[k-1].distance;
  int current_node_number = -1;
  int * curr_size = ((siss_query_worker_data*)qwdata)->curr_k_size;
  while (1) 
    {

      current_node_number = __sync_fetch_and_add(((siss_query_worker_data*)qwdata)->node_counter,1);
      if(current_node_number >= ((siss_query_worker_data*)qwdata)->node_list_size)
	break;

      current_node = ((siss_query_worker_data*)qwdata)->node_list[current_node_number];
      
      search_worker_insert_node(query_ts,
				index,
				current_node,
				knn_results,
				k,
				pq,
				delta,
				epsilon,
				((siss_query_worker_data*)qwdata)->lock_queue);
    }



  
  pthread_barrier_wait(((siss_query_worker_data*)qwdata)->lock_barrier);

  COUNT_PARTIAL_TIME_END
    COUNT_PARTIAL_TIME_START
    
    while (1)
      {
	pthread_mutex_lock(((siss_query_worker_data*)qwdata)->lock_queue);
	n = pqueue_pop(pq);
	pthread_mutex_unlock(((siss_query_worker_data*)qwdata)->lock_queue);
	if(n==NULL)
	  break;
	kth_bsf = knn_results[k-1].distance;
	if (n->distance > kth_bsf/(1 + epsilon))
	  {
	    break;
	  }
	else 
	  {
	    if (n->node->is_leaf)
	      {
		calculate_node_knn_distance_psq(index, n->node, query_ts,query_ts_reordered,
						query_order, offset, kth_bsf,
						k,knn_results,NULL, NULL,curr_size,
						qwdata);
	      }
            
	  }
	free(n);
      }

}



//reading one series at a time
void* exact_de_knn_search_worker_lf(void *qwdata)
{

  
  struct hercules_node *current_node;
  struct query_result n;
  struct hercules_index *index=((siss_query_worker_data*)qwdata)->index;
  int *query_order=((siss_query_worker_data*)qwdata)->query_order;
  ts_type *query_ts_reordered=((siss_query_worker_data*)qwdata)->query_ts_reordered;
  ts_type *query_ts=((siss_query_worker_data*)qwdata)->query_ts;
  pqueue_t *pq=((siss_query_worker_data*)qwdata)->pq;
  struct query_result * candidates=((siss_query_worker_data*)qwdata)->candidates;  
  unsigned int * candidates_count=((siss_query_worker_data*)qwdata)->candidates_count;
  unsigned int offset=((siss_query_worker_data*)qwdata)->offset;
  ts_type minimum_distance=((siss_query_worker_data*)qwdata)->minimum_distance;
  ts_type epsilon=((siss_query_worker_data*)qwdata)->epsilon;
  ts_type delta=((siss_query_worker_data*)qwdata)->delta;    
  unsigned int k = ((siss_query_worker_data*)qwdata)->k;  
  struct query_result *knn_results = (((siss_query_worker_data*)qwdata)->knn_results);
  ts_type kth_bsf=knn_results[k-1].distance;
  int current_node_number = -1;
  int * curr_size = ((siss_query_worker_data*)qwdata)->curr_k_size;
  int ts_length = index->settings->timeseries_size;

  unsigned int idx = -1;

  while (1) 
    {

      current_node_number = __sync_fetch_and_add(((siss_query_worker_data*)qwdata)->node_counter,1);
      if(current_node_number >= ((siss_query_worker_data*)qwdata)->node_list_size)
	break;

      current_node = ((siss_query_worker_data*)qwdata)->node_list[current_node_number];
      
      search_worker_insert_node_lf(query_ts,
				   index,
				   current_node,
				   knn_results,
				   k,
				   candidates,				   
				   ((siss_query_worker_data*)qwdata)->candidates_count,
				   ((siss_query_worker_data*)qwdata)->candidates_series_count,
				   delta,
				   epsilon,
				   ((siss_query_worker_data*)qwdata)->lock_queue);
    }
  
}
void * exact_de_knn_read_worker_lf_sims(void *qwdata)
{
  struct hercules_index *index=((siss_query_worker_data*)qwdata)->index;
  FILE *raw_file = fopen(index->leaves_raw_filename, "rb");
  fseek(raw_file, 0, SEEK_SET);
  ts_type *ts_buffer = malloc(index->settings->ts_byte_size);
  ts_type *ts =((siss_query_worker_data*)qwdata)->query_ts;
  unsigned long  t=0;
  unsigned long p = 0;
  unsigned int k = ((siss_query_worker_data*)qwdata)->k;  
  unsigned long sum_of_lab=((siss_query_worker_data*)qwdata)->sum_of_lab;
  float *minidisvector=((siss_query_worker_data*)qwdata)->minidisvector;
  struct query_result *knn_results = (((siss_query_worker_data*)qwdata)->knn_results);
  ts_type kth_bsf=knn_results[k-1].distance;
  int * curr_size = ((siss_query_worker_data*)qwdata)->curr_k_size;
  ts_type epsilon=((siss_query_worker_data*)qwdata)->epsilon;
  ts_type delta=((siss_query_worker_data*)qwdata)->delta;    

  float bsf,dist;

#if DETAILED_STATS == 1
  double tS = 0;
  double tE = 0;
  struct timeval start_time;
  struct timeval end_time;
#endif
    
  while(1)
    {          
      kth_bsf= (((siss_query_worker_data*)qwdata)->knn_results[k-1].distance); 
      if (t>=sum_of_lab) 
        {    
	  break; 
        } 
        
      p=((siss_query_worker_data*)qwdata)->load_point[t];
      if (minidisvector[t] <= kth_bsf)        
	{
#if DETAILED_STATS == 1
	  gettimeofday(&start_time, NULL);
	  fseek(raw_file, p * index->settings->ts_byte_size, SEEK_SET); 
	  fread(ts_buffer, index->settings->ts_byte_size, 1, raw_file);
	  gettimeofday(&end_time, NULL);
	  tS = start_time.tv_sec*1000000 + (start_time.tv_usec);
	  tE = end_time.tv_sec*1000000 + (end_time.tv_usec);
	  ((siss_query_worker_data*)qwdata)->thread_input_time += (tE - tS);
	  ++(((siss_query_worker_data*)qwdata)->thread_random_io);
#else
	  fseek(raw_file, p * index->settings->ts_byte_size, SEEK_SET); 
	  fread(ts_buffer, index->settings->ts_byte_size, 1, raw_file);
#endif
	  dist = ts_euclidean_distance_SIMD(ts, ts_buffer, index->settings->timeseries_size, kth_bsf); 
	  if(dist < kth_bsf)  
	    {  
	      struct query_result object_result;
	      object_result.node = NULL;
	      object_result.distance =  dist;		
	      pthread_rwlock_wrlock(((siss_query_worker_data*)qwdata)->lock_bsf);
	      kth_bsf= (((siss_query_worker_data*)qwdata)->knn_results[k-1].distance); 
	      if (dist < kth_bsf)
		{
		  queue_bounded_sorted_insert(knn_results, object_result, curr_size, k);
		}
	      pthread_rwlock_unlock(((siss_query_worker_data*)qwdata)->lock_bsf);
	    } 
        } 
    }

  free(ts_buffer);
  fclose(raw_file);
} 

 
void* exact_de_knn_read_worker_lf(void *qwdata)
{

  struct hercules_index *index=((siss_query_worker_data*)qwdata)->index;
  FILE * raw_file = fopen(index->leaves_raw_filename, "rb");
  	  
  fseek(raw_file, 0, SEEK_SET);

  struct hercules_node *current_node;
  struct query_result n;
  ts_type *query_ts=((siss_query_worker_data*)qwdata)->query_ts;
  struct query_result * candidates=((siss_query_worker_data*)qwdata)->candidates;  
  unsigned int * candidates_count=((siss_query_worker_data*)qwdata)->candidates_count;
  ts_type epsilon=((siss_query_worker_data*)qwdata)->epsilon;
  ts_type delta=((siss_query_worker_data*)qwdata)->delta;    
  unsigned int k = ((siss_query_worker_data*)qwdata)->k;  
  struct query_result *knn_results = (((siss_query_worker_data*)qwdata)->knn_results);
  ts_type kth_bsf=knn_results[k-1].distance;
  int * curr_size = ((siss_query_worker_data*)qwdata)->curr_k_size;
  ts_type * ts_list = (((siss_query_worker_data*)qwdata)->ts_list);
  int ts_length = index->settings->timeseries_size;
 
  unsigned int idx = -1;

  ts_type distance = FLT_MAX;

  double tS = 0;
  double tE = 0;
  struct timeval start_time;
  struct timeval end_time;

  while(1)
    { 
      if(*((siss_query_worker_data*)qwdata)->candidates_idx >=*candidates_count)
	break;

      kth_bsf= (((siss_query_worker_data*)qwdata)->knn_results[k-1].distance); 

      n = ((siss_query_worker_data*)qwdata)->candidates[__sync_fetch_and_add(((siss_query_worker_data*)qwdata)->candidates_idx,1)];

      if  (n.distance <= kth_bsf/(1 + epsilon))
	{
#if DETAILED_STATS == 1 
	  pthread_mutex_lock(((siss_query_worker_data*)qwdata)->lock_queue);     
	  gettimeofday(&start_time, NULL);
	  fseek(index->leaves_raw_file,
		n.node->file_pos * ts_length * sizeof(ts_type),
		SEEK_SET);
	  fread(ts_list,sizeof(ts_type),
		ts_length*n.node->node_size,
		index->leaves_raw_file);
	  gettimeofday(&end_time, NULL);
	  pthread_mutex_unlock(((siss_query_worker_data*)qwdata)->lock_queue);     
	  tS = start_time.tv_sec*1000000 + (start_time.tv_usec);
	  tE = end_time.tv_sec*1000000 + (end_time.tv_usec);
	  ((siss_query_worker_data*)qwdata)->thread_input_time += (tE - tS);
#else
	  fseek(raw_file,
		n.node->file_pos * ts_length * sizeof(ts_type),
		SEEK_SET);
	  fread(ts_list,sizeof(ts_type),
		ts_length*n.node->node_size,
		raw_file);
#endif
	      
#if DETAILED_STATS == 1
	  gettimeofday(&start_time, NULL);
#endif

	      
	  for (int idx = 0; idx < n.node->node_size; ++idx)
	    {
	      kth_bsf =  knn_results[k-1].distance;
	      distance = ts_euclidean_distance_SIMD(query_ts,
						    &ts_list[idx*ts_length],
						    ts_length,
						    kth_bsf);
	      if (distance < kth_bsf)
		{
		  struct query_result object_result;
		  object_result.node = n.node;
		  object_result.distance =  distance;		
		  pthread_rwlock_wrlock(((siss_query_worker_data*)qwdata)->lock_bsf);
		  kth_bsf= (((siss_query_worker_data*)qwdata)->knn_results[k-1].distance); 
		  if (distance < kth_bsf)
		    {
		      queue_bounded_sorted_insert(knn_results, object_result, curr_size, k);
		    }

		  pthread_rwlock_unlock(((siss_query_worker_data*)qwdata)->lock_bsf);
		}	    
	    }	  
        
#if DETAILED_STATS == 1
	  gettimeofday(&end_time, NULL);
	  tS = start_time.tv_sec*1000000 + (start_time.tv_usec);
	  tE = end_time.tv_sec*1000000 + (end_time.tv_usec);
	  ((siss_query_worker_data*)qwdata)->thread_realdist_time += (tE - tS);
#endif

	}
	  
	    
    }           
  fclose(raw_file);
}
 
void* exact_de_knn_mindist_worker_lf(void *qwdata)
{

  struct hercules_index *index=((siss_query_worker_data*)qwdata)->index;

  struct hercules_node *current_node;
  struct query_result n;
  ts_type *query_ts=((siss_query_worker_data*)qwdata)->query_ts;
  struct query_result * candidates=((siss_query_worker_data*)qwdata)->candidates;  
  unsigned int * candidates_count=((siss_query_worker_data*)qwdata)->candidates_count;
  ts_type epsilon=((siss_query_worker_data*)qwdata)->epsilon;
  ts_type delta=((siss_query_worker_data*)qwdata)->delta;    
  unsigned int k = ((siss_query_worker_data*)qwdata)->k;  
  struct query_result *knn_results = (((siss_query_worker_data*)qwdata)->knn_results);
  ts_type kth_bsf=knn_results[k-1].distance;
  int * curr_size = ((siss_query_worker_data*)qwdata)->curr_k_size;
  int ts_length = index->settings->timeseries_size;
  ts_type *paa=((siss_query_worker_data*)qwdata)->query_paa;
  
  ((siss_query_worker_data*)qwdata)->label_number=malloc(sizeof(unsigned long)*10000);
  ((siss_query_worker_data*)qwdata)->minidisvector=malloc(sizeof(float)*10000);

  unsigned int idx = -1;
  unsigned long max_number=10000;

  ts_type distance = FLT_MAX;
  float mindist;
  unsigned int  candidates_idx;

  // fseek(index->leaves_raw_file, 0, SEEK_SET);


  while(1)
    { 
      candidates_idx=__sync_fetch_and_add(((siss_query_worker_data*)qwdata)->candidates_idx,1);
         
      if(candidates_idx >=*candidates_count)
	break;
	
      n = ((siss_query_worker_data*)qwdata)->candidates[candidates_idx];

      if  (n.distance <= kth_bsf)	
	{
	  for (int i = 0; i < n.node->node_size; ++i )
	    {	    
	      unsigned long idx = n.node->file_pos + i;
	      sax_type *sax = &index->sax_cache[(unsigned long long) (idx *index->settings->paa_segments)];
	      mindist = minidist_paa_to_isax_raw_SIMD((float*)paa, sax, index->settings->max_sax_cardinalities,
						      index->settings->sax_bit_cardinality,
						      index->settings->sax_alphabet_cardinality,
						      index->settings->paa_segments, MINVAL, MAXVAL,
						      (float)index->settings->mindist_sqrt);
	      if(mindist <= kth_bsf) 
		{
		  if ( ((siss_query_worker_data*)qwdata)->sum_of_lab>=max_number)
		    {
		      max_number=(max_number+10000);
		      unsigned long* change_lab=((siss_query_worker_data*)qwdata)->label_number;
		      unsigned long* change_minivec=((siss_query_worker_data*)qwdata)->minidisvector;
		      ((siss_query_worker_data*)qwdata)->label_number=malloc(sizeof(unsigned long)*(max_number+10000));
		      ((siss_query_worker_data*)qwdata)->minidisvector=malloc(sizeof(float)*(max_number+10000));
		      memcpy(((siss_query_worker_data*)qwdata)->label_number,change_lab,sizeof(unsigned long)*max_number);
		      memcpy(((siss_query_worker_data*)qwdata)->minidisvector,change_minivec,sizeof(float)*max_number);
		      free(change_lab);
		      free(change_minivec);
		  
		    }
		  ((siss_query_worker_data*)qwdata)->label_number[((siss_query_worker_data*)qwdata)->sum_of_lab]= idx;
		  ((siss_query_worker_data*)qwdata)->minidisvector[((siss_query_worker_data*)qwdata)->sum_of_lab]=mindist;
		  ((siss_query_worker_data*)qwdata)->sum_of_lab++;
		}
	    }
	}	 	    
    }           
}
void  exact_de_knn_search (ts_type *query_ts, ts_type * query_ts_reordered,
			   int * query_order, unsigned int offset,
			   struct hercules_index *index,ts_type minimum_distance,
			   ts_type epsilon, ts_type r_delta,
			   unsigned int k, unsigned int q_id, char * qfilename, query_settings q_settings)
{
    
  unsigned int curr_size = 0;
  ts_type bsf = FLT_MAX;
  ts_type kth_bsf = FLT_MAX;
  ts_type temp_bsf = FLT_MAX;

  ts_type * ts_list  =  malloc(index->settings->max_leaf_size * index->settings->timeseries_size * sizeof(ts_type));
  const char *leaves_filename = malloc(sizeof(char) * (strlen(index->settings->root_directory) + 18));
  leaves_filename = strcpy(leaves_filename, index->settings->root_directory);
  leaves_filename = strcat(leaves_filename, "leaves_raw.idx\0");
    

  FILE *ifile= NULL;
  ifile = fopen(leaves_filename,"rb");
  fseek(ifile, 0, SEEK_SET);    
  unsigned int found_knn = 0;
    
  struct query_result * knn_results = calloc(k,sizeof(struct query_result));
  for (int idx = 0; idx < k; ++idx)
    {
      knn_results[idx].node = NULL;
      knn_results[idx].distance = FLT_MAX;      
    }

  //return k approximate results
  approximate_knn_search(query_ts, query_ts_reordered, query_order, offset, bsf,
			 index,knn_results,k, NULL, NULL,&curr_size, q_settings.serial);
    
    
  //set the approximate result to be the first item in the queue
  struct query_result approximate_result = knn_results[0];    

  COUNT_PARTIAL_TIME_END                 
    index->stats->query_filter_total_time  = partial_time;	
    
  index->stats->query_filter_input_time  = partial_input_time;
  index->stats->query_filter_output_time = partial_output_time;
  index->stats->query_filter_load_node_time = partial_load_node_time;
  index->stats->query_filter_cpu_time    = partial_time-partial_input_time-partial_output_time;
  index->stats->query_filter_seq_input_count   = partial_seq_input_count;
  index->stats->query_filter_seq_output_count  = partial_seq_output_count;
  index->stats->query_filter_rand_input_count  = partial_rand_input_count;
  index->stats->query_filter_rand_output_count = partial_rand_output_count;
    
  index->stats->query_filter_loaded_nodes_count = loaded_nodes_count;
  index->stats->query_filter_loaded_ts_count = loaded_ts_count;
  index->stats->query_filter_checked_nodes_count = checked_nodes_count;
  index->stats->query_filter_checked_ts_count = checked_ts_count;

  index->stats->query_approx_distance = sqrtf(approximate_result.distance);
  index->stats->query_approx_node_filename = approximate_result.node->filename;
  index->stats->query_approx_node_size = approximate_result.node->node_size;;
  index->stats->query_approx_node_level = approximate_result.node->level;
    
  index->stats->queries_filter_total_time  +=  index->stats->query_filter_total_time;
    
  index->stats->queries_filter_input_time  +=  index->stats->query_filter_input_time;
  index->stats->queries_filter_output_time +=  index->stats->query_filter_output_time;
  index->stats->queries_filter_load_node_time += index->stats->query_filter_load_node_time;
  index->stats->queries_filter_cpu_time    += index->stats->query_filter_cpu_time;

  index->stats->queries_filter_seq_input_count   += index->stats->query_filter_seq_input_count;
  index->stats->queries_filter_seq_output_count  += index->stats->query_filter_seq_output_count;
  index->stats->queries_filter_rand_input_count  += index->stats->query_filter_rand_input_count;
  index->stats->queries_filter_rand_output_count += index->stats->query_filter_rand_output_count;
    
    
  if(approximate_result.node != NULL) {  
    index->stats->query_approx_distance = sqrtf(approximate_result.distance);
    index->stats->query_approx_node_filename = approximate_result.node->filename;
    index->stats->query_approx_node_size = approximate_result.node->node_size;;
    index->stats->query_approx_node_level = approximate_result.node->level;
  }
 
  // Early termination...
  if (approximate_result.distance == 0) {
    index->stats->query_exact_distance = sqrtf(approximate_result.distance);
    index->stats->query_exact_node_filename = approximate_result.node->filename;
    index->stats->query_exact_node_size = approximate_result.node->node_size;;
    index->stats->query_exact_node_level = approximate_result.node->level; 
    update_query_stats(index,q_id, found_knn, approximate_result);
  }

    
  RESET_QUERY_COUNTERS()
    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
 
    struct query_result bsf_result = approximate_result;

  struct query_result * non_pruned_leaves =  calloc(index->stats->leaves_counter,
						    sizeof(struct query_result));
  int non_pruned_count = 0;
    
  pqueue_t *pq = pqueue_init(index->first_node->node_size, 
			     cmp_pri, get_pri, set_pri, get_pos, set_pos);
	

  struct query_result * root_pq_item = malloc(sizeof(struct query_result));
  root_pq_item->node = index->first_node;
  root_pq_item->distance = calculate_node_min_distance (index, index->first_node, query_ts);
    
  //initialize the lb distance to be the distance of the query to the root.

  pqueue_insert(pq, root_pq_item);    

  struct query_result * n;
  struct query_result temp;

    

  while ((n = pqueue_pop(pq)))
    {
      kth_bsf =  knn_results[k-1].distance;
      if (n->distance > kth_bsf/(1 + epsilon))
	{
          break;
        }


      if (n->node->is_leaf) // n is a leaf
        {
	  non_pruned_leaves[non_pruned_count].node = n->node;
	  non_pruned_leaves[non_pruned_count].distance = n->distance;
          ++non_pruned_count;
        }
      else  //n is an internal node
        {
          temp = knn_results[k-1];
          kth_bsf =  temp.distance;
	  
          ts_type child_distance;
	  child_distance = calculate_node_min_distance(index, n->node->left_child,query_ts);
	  
	  if ((child_distance < kth_bsf/(1+epsilon)) &&
	      (n->node->left_child != approximate_result.node)) 
	    {
	      struct query_result * mindist_result_left = malloc(sizeof(struct query_result));
	      mindist_result_left->node = n->node->left_child;
	      mindist_result_left->distance =  child_distance;
	      pqueue_insert(pq, mindist_result_left);
	    }

	  child_distance = calculate_node_min_distance(index, n->node->right_child,query_ts);
	  if ((child_distance < kth_bsf/(1+epsilon))  &&
	      (n->node->right_child != approximate_result.node)) 
	    {
	      struct query_result * mindist_result_right = malloc(sizeof(struct query_result));
	      mindist_result_right->node = n->node->right_child;
	      mindist_result_right->distance =  child_distance;
	      pqueue_insert(pq, mindist_result_right);
	    }
        }
      free(n);
    }
  while ((n = pqueue_pop(pq)))
    {
      free(n);
    }
  pqueue_free(pq);

  qsort(non_pruned_leaves,
	non_pruned_count,
	sizeof(struct query_result),
	compare_leaf_pos);

  ts_type distance = FLT_MAX;
  int ts_length = index->settings->timeseries_size;
    
    
  for (int i = 0; i < non_pruned_count; ++i)
    {

      if (non_pruned_leaves[i].distance <= kth_bsf/(1 + epsilon))
	{
	  COUNT_LOADED_NODE      
	    COUNT_LOADED_TS(non_pruned_leaves[i].node->node_size)

	    COUNT_PARTIAL_INPUT_TIME_START
	    fseek(ifile,
		  non_pruned_leaves[i].node->file_pos * ts_length * sizeof(ts_type),
		  SEEK_SET);
	  
	  fread(ts_list,sizeof(ts_type),
		ts_length*non_pruned_leaves[i].node->node_size,
		ifile);	    
	  COUNT_PARTIAL_INPUT_TIME_END       
	  
	    for (int idx = 0; idx < non_pruned_leaves[i].node->node_size; ++idx)
	      {
		kth_bsf =  knn_results[k-1].distance;
		distance = ts_euclidean_distance_SIMD(query_ts,
						      &ts_list[idx*ts_length],
						      ts_length,
						      kth_bsf);
		if (distance < kth_bsf)
		  {
		    struct query_result object_result;// =  malloc(sizeof(struct query_result));
		    object_result.node = non_pruned_leaves[i].node;
		    object_result.distance =  distance;				
		    queue_bounded_sorted_insert(knn_results, object_result, &curr_size, k);
		  }	    
	      }	  
	}
    }
    
  for (unsigned int pos = found_knn; pos < k; ++pos)
    {
      bsf_result = knn_results[pos];
      found_knn = pos+1;
      COUNT_PARTIAL_TIME_END
        update_query_stats(index,q_id, found_knn, bsf_result);
      get_query_stats(index, found_knn);
      print_query_stats(index, q_id, found_knn,qfilename);	
      RESET_QUERY_COUNTERS()
        RESET_PARTIAL_COUNTERS()
        COUNT_PARTIAL_TIME_START	
	}
    
      
  //free the results, eventually do something with them!!

  free(non_pruned_leaves);
  free(knn_results);
  free(ts_list);
  free(leaves_filename);
  fclose(ifile);
}

void  exact_de_knn_search_noapprox (ts_type *query_ts, ts_type * query_ts_reordered,
				    int * query_order, unsigned int offset,
				    struct hercules_index *index,ts_type minimum_distance,
				    ts_type epsilon, ts_type r_delta,
				    struct query_result * knn_results, unsigned int k, unsigned int q_id, char * qfilename, query_settings q_settings,
				    struct query_result * candidates, unsigned int candidates_count, unsigned int *candidates_idx,
				    ts_type *ts_list, FILE *ifile)
{
    
  unsigned int curr_size = k;
  ts_type kth_bsf = knn_results[k-1].distance;
  ts_type temp_bsf = FLT_MAX;

  fseek(ifile, 0, SEEK_SET);    
  unsigned int found_knn = 0;

  struct query_result bsf_result = knn_results[k-1];


  ts_type distance = FLT_MAX;
  int ts_length = index->settings->timeseries_size;
    
  for (int i = 0; i < candidates_count; ++i)
    {
      if (candidates[i].distance <= kth_bsf/(1 + epsilon))
	{
	  COUNT_LOADED_NODE      
	    COUNT_LOADED_TS(candidates[i].node->node_size)

	    COUNT_PARTIAL_INPUT_TIME_START
	    fseek(ifile,
		  candidates[i].node->file_pos * ts_length * sizeof(ts_type),
		  SEEK_SET);
	  
	  fread(ts_list,sizeof(ts_type),
		ts_length*candidates[i].node->node_size,
		ifile);	    
	  COUNT_PARTIAL_INPUT_TIME_END       
	  
	    for (int idx = 0; idx < candidates[i].node->node_size; ++idx)
	      {
		kth_bsf =  knn_results[k-1].distance;
		distance = ts_euclidean_distance_SIMD(query_ts,
						      &ts_list[idx*ts_length],
						      ts_length,
						      kth_bsf);
		if (distance < kth_bsf)
		  {

		    struct query_result object_result;// =  malloc(sizeof(struct query_result));
		    object_result.node = candidates[i].node;
		    object_result.distance =  distance;				
		    queue_bounded_sorted_insert(knn_results, object_result, &curr_size, k);
		  }	    
	      }	  
	}
    }
}

void search_worker_insert_node(ts_type *query_ts,
			       struct hercules_index *index,
			       struct hercules_node *node,
			       struct query_result * knn_results,
			       int k,
			       pqueue_t *pq,
			       ts_type delta,
			       ts_type epsilon,
			       pthread_mutex_t *lock_queue)
{
  
  ts_type distance = calculate_node_min_distance(index,
						 node,
						 query_ts);
  ts_type kth_bsf;
  kth_bsf =  knn_results[k-1].distance;
      
  if(distance < kth_bsf/(1+epsilon))
    {
      if (node->is_leaf) 
	{   
	  struct query_result * mindist_result = malloc(sizeof(struct query_result));
	  mindist_result->node = node;
	  mindist_result->distance=distance;
	  pthread_mutex_lock(lock_queue);
	  pqueue_insert(pq, mindist_result);
	  pthread_mutex_unlock(lock_queue);
	}
      else
	{   
	  search_worker_insert_node(query_ts,
				    index,
				    node->left_child,
				    knn_results,
				    k,
				    pq,
				    delta,
				    epsilon,
				    lock_queue);

	  search_worker_insert_node(query_ts,
				    index,
				    node->right_child,
				    knn_results,
				    k,
				    pq,
				    delta,
				    epsilon,
				    lock_queue);
	}
    }
}
void search_worker_insert_node_lf(ts_type *query_ts,
				  struct hercules_index *index,
				  struct hercules_node *node,
				  struct query_result * knn_results,
				  int k,
				  struct query_result *candidates,
				  unsigned int *candidates_count,
				  unsigned int *candidates_series_count,
				  ts_type delta,
				  ts_type epsilon,
				  pthread_mutex_t *lock_queue)
{
  
  ts_type distance = calculate_node_min_distance(index,
						 node,
						 query_ts);
  ts_type kth_bsf;
  kth_bsf =  knn_results[k-1].distance;
  
  if(distance < kth_bsf/(1+epsilon))
    {
      if (node->is_leaf) 
	{   
	  struct query_result cand_result;
	  cand_result.node = node;
	  cand_result.distance=distance;

	  pthread_mutex_lock(lock_queue);
	  candidates[*candidates_count] = cand_result;
	  ++(*candidates_count);
	  (*candidates_series_count) = 	  (*candidates_series_count) + node->node_size;
	  pthread_mutex_unlock(lock_queue);
	}
      else
	{   
	  search_worker_insert_node_lf(query_ts,
				       index,
				       node->left_child,
				       knn_results,
				       k,
				       candidates,
				       candidates_count,
				       candidates_series_count,
				       delta,
				       epsilon,
				       lock_queue);
	  
	  search_worker_insert_node_lf(query_ts,
				       index,
				       node->right_child,
				       knn_results,
				       k,
				       candidates,
				       candidates_count,
				       candidates_series_count,
				       delta,
				       epsilon,
				       lock_queue);
	}
    }
}

void  exact_ng_knn_search (ts_type *query_ts, ts_type * query_ts_reordered,
			   int * query_order, unsigned int offset,
			   struct hercules_index *index,ts_type minimum_distance,
			   unsigned int k, unsigned int q_id, char * qfilename, unsigned int nprobes, query_settings q_settings)
{
    
  unsigned int curr_size = 0;
  ts_type bsf = FLT_MAX;
  ts_type kth_bsf = FLT_MAX;
  ts_type temp_bsf = FLT_MAX;
  unsigned int cur_probes = 0;
    
  unsigned int found_knn = 0;
    
  struct query_result * knn_results = calloc(k,sizeof(struct query_result));
  for (int idx = 0; idx < k; ++idx)
    {
      knn_results[idx].node = NULL;
      knn_results[idx].distance = FLT_MAX;      
    }

  approximate_knn_search(query_ts, query_ts_reordered, query_order, offset, bsf,
			 index,knn_results,k, NULL, NULL,&curr_size,q_settings.serial);
    
  ++cur_probes;
    
  struct query_result approximate_result = knn_results[0];    

  COUNT_PARTIAL_TIME_END                 
    index->stats->query_filter_total_time  = partial_time;	
    
  index->stats->query_filter_input_time  = partial_input_time;
  index->stats->query_filter_output_time = partial_output_time;
  index->stats->query_filter_load_node_time = partial_load_node_time;
  index->stats->query_filter_cpu_time    = partial_time-partial_input_time-partial_output_time;
  index->stats->query_filter_seq_input_count   = partial_seq_input_count;
  index->stats->query_filter_seq_output_count  = partial_seq_output_count;
  index->stats->query_filter_rand_input_count  = partial_rand_input_count;
  index->stats->query_filter_rand_output_count = partial_rand_output_count;
    
  index->stats->query_filter_loaded_nodes_count = loaded_nodes_count;
  index->stats->query_filter_loaded_ts_count = loaded_ts_count;
  index->stats->query_filter_checked_nodes_count = checked_nodes_count;
  index->stats->query_filter_checked_ts_count = checked_ts_count;

  index->stats->query_approx_distance = sqrtf(approximate_result.distance);
  index->stats->query_approx_node_filename = approximate_result.node->filename;
  index->stats->query_approx_node_size = approximate_result.node->node_size;;
  index->stats->query_approx_node_level = approximate_result.node->level;
    
  index->stats->queries_filter_total_time  +=  index->stats->query_filter_total_time;
    
  index->stats->queries_filter_input_time  +=  index->stats->query_filter_input_time;
  index->stats->queries_filter_output_time +=  index->stats->query_filter_output_time;
  index->stats->queries_filter_load_node_time += index->stats->query_filter_load_node_time;
  index->stats->queries_filter_cpu_time    += index->stats->query_filter_cpu_time;

  index->stats->queries_filter_seq_input_count   += index->stats->query_filter_seq_input_count;
  index->stats->queries_filter_seq_output_count  += index->stats->query_filter_seq_output_count;
  index->stats->queries_filter_rand_input_count  += index->stats->query_filter_rand_input_count;
  index->stats->queries_filter_rand_output_count += index->stats->query_filter_rand_output_count;
    
    
  if(approximate_result.node != NULL) {  
    index->stats->query_approx_distance = sqrtf(approximate_result.distance);
    index->stats->query_approx_node_filename = approximate_result.node->filename;
    index->stats->query_approx_node_size = approximate_result.node->node_size;;
    index->stats->query_approx_node_level = approximate_result.node->level;
  }
 
  // Early termination...
  if (approximate_result.distance == 0) {
    index->stats->query_exact_distance = sqrtf(approximate_result.distance);
    index->stats->query_exact_node_filename = approximate_result.node->filename;
    index->stats->query_exact_node_size = approximate_result.node->node_size;;
    index->stats->query_exact_node_level = approximate_result.node->level; 
    update_query_stats(index,q_id, found_knn, approximate_result);
    get_query_stats(index, found_knn);
    print_query_stats(index, q_id, found_knn,qfilename);	     
  }

    
  RESET_QUERY_COUNTERS()
    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
 
    struct query_result bsf_result = approximate_result;
    
   
  pqueue_t *pq = pqueue_init(index->first_node->node_size, 
			     cmp_pri, get_pri, set_pri, get_pos, set_pos);
	

  struct query_result * root_pq_item = malloc(sizeof(struct query_result));
  root_pq_item->node = index->first_node;
  root_pq_item->distance = calculate_node_min_distance (index, index->first_node, query_ts);
    

  pqueue_insert(pq, root_pq_item);    

  struct query_result * n;
  struct query_result temp;


    
  while ((n = pqueue_pop(pq)) &&  (cur_probes < nprobes))
    {
      
      if (n->distance > bsf_result.distance)
	{
          break;
        }


      if (n->node->is_leaf) // n is a leaf
        {
 	  calculate_node_knn_distance(index, n->node, query_ts,query_ts_reordered,
				      query_order, offset, bsf_result.distance,
				      k,knn_results,NULL, NULL,&curr_size,q_settings.serial);
          
	  ++cur_probes;
	  
        }
      else  //n is an internal node
        {
          temp = knn_results[k-1];
          kth_bsf =  temp.distance;
	  
          ts_type child_distance;
	  child_distance = calculate_node_min_distance(index, n->node->left_child,query_ts);
	  
	  if ((child_distance < kth_bsf) &&
	      (n->node->left_child != approximate_result.node)) 
	    {
	      struct query_result * mindist_result_left = malloc(sizeof(struct query_result));
	      mindist_result_left->node = n->node->left_child;
	      mindist_result_left->distance =  child_distance;
	      pqueue_insert(pq, mindist_result_left);
	    }

	  child_distance = calculate_node_min_distance(index, n->node->right_child,query_ts);

	  if ((child_distance < kth_bsf)  &&
	      (n->node->right_child != approximate_result.node))
	    {
	      struct query_result * mindist_result_right = malloc(sizeof(struct query_result));
	      mindist_result_right->node = n->node->right_child;
	      mindist_result_right->distance =  child_distance;
	      pqueue_insert(pq, mindist_result_right);
	    }
        }
      free(n);
    }
  while ((n = pqueue_pop(pq)))
    {
      free(n);
    }
  pqueue_free(pq);
    
  for (unsigned int pos = found_knn; pos < k; ++pos)
    {
      bsf_result = knn_results[pos];
      found_knn = pos+1;
      COUNT_PARTIAL_TIME_END
        update_query_stats(index,q_id, found_knn, bsf_result);
      get_query_stats(index, found_knn);
      print_query_stats(index, q_id, found_knn,qfilename);	
      RESET_QUERY_COUNTERS()
        RESET_PARTIAL_COUNTERS()
        COUNT_PARTIAL_TIME_START	
	}


  free(knn_results);	
}



void  exact_knn_search_max_policy (ts_type *query_ts, ts_type * query_ts_reordered,
				   int * query_order, unsigned int offset,
				   struct hercules_index *index,ts_type minimum_distance,
				   ts_type epsilon, ts_type delta,
				   unsigned int k, unsigned int q_id, char * qfilename,query_settings q_settings)
{
    
  unsigned int curr_size = 0;
  ts_type bsf = FLT_MAX;
  ts_type kth_bsf = FLT_MAX;
  ts_type temp_bsf = FLT_MAX;
    
  unsigned int found_knn = 0;
    
  struct query_result * knn_results = calloc(k,sizeof(struct query_result));
  for (int idx = 0; idx < k; ++idx)
    {
      knn_results[idx].node = NULL;
      knn_results[idx].distance = FLT_MAX;      
    }
      

  approximate_knn_search(query_ts, query_ts_reordered, query_order, offset, bsf,
			 index,knn_results,k, NULL, NULL,&curr_size,q_settings.serial);
    

  struct query_result approximate_result = knn_results[0];    

  COUNT_PARTIAL_TIME_END                 
    index->stats->query_filter_total_time  = partial_time;	
    
  index->stats->query_filter_input_time  = partial_input_time;
  index->stats->query_filter_output_time = partial_output_time;
  index->stats->query_filter_load_node_time = partial_load_node_time;
  index->stats->query_filter_cpu_time    = partial_time-partial_input_time-partial_output_time;
  index->stats->query_filter_seq_input_count   = partial_seq_input_count;
  index->stats->query_filter_seq_output_count  = partial_seq_output_count;
  index->stats->query_filter_rand_input_count  = partial_rand_input_count;
  index->stats->query_filter_rand_output_count = partial_rand_output_count;
    
  index->stats->query_filter_loaded_nodes_count = loaded_nodes_count;
  index->stats->query_filter_loaded_ts_count = loaded_ts_count;
  index->stats->query_filter_checked_nodes_count = checked_nodes_count;
  index->stats->query_filter_checked_ts_count = checked_ts_count;

  index->stats->query_approx_distance = sqrtf(approximate_result.distance);
  index->stats->query_approx_node_filename = approximate_result.node->filename;
  index->stats->query_approx_node_size = approximate_result.node->node_size;;
  index->stats->query_approx_node_level = approximate_result.node->level;
    
  index->stats->queries_filter_total_time  +=  index->stats->query_filter_total_time;
    
  index->stats->queries_filter_input_time  +=  index->stats->query_filter_input_time;
  index->stats->queries_filter_output_time +=  index->stats->query_filter_output_time;
  index->stats->queries_filter_load_node_time += index->stats->query_filter_load_node_time;
  index->stats->queries_filter_cpu_time    += index->stats->query_filter_cpu_time;

  index->stats->queries_filter_seq_input_count   += index->stats->query_filter_seq_input_count;
  index->stats->queries_filter_seq_output_count  += index->stats->query_filter_seq_output_count;
  index->stats->queries_filter_rand_input_count  += index->stats->query_filter_rand_input_count;
  index->stats->queries_filter_rand_output_count += index->stats->query_filter_rand_output_count;
    
    
  if(approximate_result.node != NULL) {  
    index->stats->query_approx_distance = sqrtf(approximate_result.distance);
    index->stats->query_approx_node_filename = approximate_result.node->filename;
    index->stats->query_approx_node_size = approximate_result.node->node_size;;
    index->stats->query_approx_node_level = approximate_result.node->level;
  }
 
  if (approximate_result.distance == 0) {
    index->stats->query_exact_distance = sqrtf(approximate_result.distance);
    index->stats->query_exact_node_filename = approximate_result.node->filename;
    index->stats->query_exact_node_size = approximate_result.node->node_size;;
    index->stats->query_exact_node_level = approximate_result.node->level; 
    update_query_stats(index,q_id, found_knn, approximate_result);
    get_query_stats(index, found_knn);
    print_query_stats(index, q_id, found_knn,qfilename);	     
  }

  RESET_QUERY_COUNTERS()
    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
 
    struct query_result bsf_result = approximate_result;
    
   
  pqueue_t *pq = pqueue_init(index->first_node->node_size, 
			     cmp_pri, get_max_pri, set_max_pri, get_pos, set_pos);
	

  struct query_result * root_pq_item = malloc(sizeof(struct query_result));
  root_pq_item->node = index->first_node;
  root_pq_item->distance = calculate_node_min_distance (index, index->first_node, query_ts);
  root_pq_item->max_distance = calculate_node_max_distance (index, index->first_node, query_ts);;
    

  pqueue_insert(pq, root_pq_item);    

  struct query_result * n;
  struct query_result temp;
    
  while ((n = pqueue_pop(pq)))
    {
      

      kth_bsf = knn_results[k-1].distance;
	
      if (n->distance > kth_bsf/(1+epsilon)) 
	{
	  continue;
	}

      if (n->node->is_leaf) // n is a leaf
	{
	  calculate_node_knn_distance(index, n->node, query_ts,query_ts_reordered,
				      query_order, offset, bsf_result.distance,
				      k,knn_results,NULL, NULL,&curr_size,q_settings.serial);
	  
	}
      else  //n is an internal node
	{
	  temp = knn_results[k-1];
	  kth_bsf =  temp.distance;
	  
	  ts_type child_distance;
	  ts_type child_max_distance;	  
	  child_distance = calculate_node_min_distance(index, n->node->left_child,query_ts);
	  child_max_distance = calculate_node_max_distance(index, n->node->left_child,query_ts);
	  if ((child_distance < kth_bsf/(1+epsilon)) &&
	      (n->node->left_child != approximate_result.node)) //add epsilon
	    {
	      struct query_result * mindist_result_left = malloc(sizeof(struct query_result));
	      mindist_result_left->node = n->node->left_child;
	      mindist_result_left->distance =  child_distance;
	      mindist_result_left->max_distance = child_max_distance;	      
	      pqueue_insert(pq, mindist_result_left);
	    }

	  child_distance = calculate_node_min_distance(index, n->node->right_child,query_ts);
	  child_max_distance = calculate_node_max_distance(index, n->node->right_child,query_ts);
	  if ((child_distance < kth_bsf/(1+epsilon))  &&
	      (n->node->right_child != approximate_result.node)) //add epsilon	  
	    {
	      struct query_result * mindist_result_right = malloc(sizeof(struct query_result));
	      mindist_result_right->node = n->node->right_child;
	      mindist_result_right->distance =  child_distance;
	      mindist_result_right->max_distance = child_max_distance;	      	      
	      pqueue_insert(pq, mindist_result_right);
	    }
	}
      free(n);
    }

  for (unsigned int pos = 1; pos <= k; ++pos)
    {
      bsf_result = knn_results[pos-1];
      COUNT_PARTIAL_TIME_END
	update_query_stats(index,q_id, pos, bsf_result);
      get_query_stats(index, pos);
      print_query_stats(index, q_id, pos,qfilename);	
      RESET_QUERY_COUNTERS()
	RESET_PARTIAL_COUNTERS()
        COUNT_PARTIAL_TIME_START	
	}
    
      
  while ((n = pqueue_pop(pq)))
    {
      free(n);
    }
  pqueue_free(pq);
  free(knn_results);	
}

void  exact_knn_search_track_pruning (ts_type *query_ts, ts_type * query_ts_reordered,
				      int * query_order, unsigned int offset,
				      struct hercules_index *index,ts_type minimum_distance,
				      ts_type epsilon, ts_type delta,
				      unsigned int k, unsigned int q_id, char * qfilename, query_settings q_settings)
{
  unsigned int curr_size = 0;    
  ts_type bsf = FLT_MAX;
  ts_type kth_bsf = FLT_MAX;
  ts_type temp_bsf = FLT_MAX;
    
  unsigned int found_knn = 0;
    
  struct query_result * knn_results = calloc(k,sizeof(struct query_result));
  for (int idx = 0; idx < k; ++idx)
    {
      knn_results[idx].node = NULL;
      knn_results[idx].distance = FLT_MAX;      
    }
      
  approximate_knn_search(query_ts, query_ts_reordered, query_order, offset, bsf,
			 index,knn_results,k, NULL, NULL,&curr_size, q_settings.serial);
    

  struct query_result approximate_result = knn_results[0];    

  COUNT_PARTIAL_TIME_END                 
    index->stats->query_filter_total_time  = partial_time;	
    
  index->stats->query_filter_input_time  = partial_input_time;
  index->stats->query_filter_output_time = partial_output_time;
  index->stats->query_filter_load_node_time = partial_load_node_time;
  index->stats->query_filter_cpu_time    = partial_time-partial_input_time-partial_output_time;
  index->stats->query_filter_seq_input_count   = partial_seq_input_count;
  index->stats->query_filter_seq_output_count  = partial_seq_output_count;
  index->stats->query_filter_rand_input_count  = partial_rand_input_count;
  index->stats->query_filter_rand_output_count = partial_rand_output_count;
    
  index->stats->query_filter_loaded_nodes_count = loaded_nodes_count;
  index->stats->query_filter_loaded_ts_count = loaded_ts_count;
  index->stats->query_filter_checked_nodes_count = checked_nodes_count;
  index->stats->query_filter_checked_ts_count = checked_ts_count;

  index->stats->query_approx_distance = sqrtf(approximate_result.distance);
  index->stats->query_approx_node_filename = approximate_result.node->filename;
  index->stats->query_approx_node_size = approximate_result.node->node_size;;
  index->stats->query_approx_node_level = approximate_result.node->level;
    
  index->stats->queries_filter_total_time  +=  index->stats->query_filter_total_time;
    
  index->stats->queries_filter_input_time  +=  index->stats->query_filter_input_time;
  index->stats->queries_filter_output_time +=  index->stats->query_filter_output_time;
  index->stats->queries_filter_load_node_time += index->stats->query_filter_load_node_time;
  index->stats->queries_filter_cpu_time    += index->stats->query_filter_cpu_time;

  index->stats->queries_filter_seq_input_count   += index->stats->query_filter_seq_input_count;
  index->stats->queries_filter_seq_output_count  += index->stats->query_filter_seq_output_count;
  index->stats->queries_filter_rand_input_count  += index->stats->query_filter_rand_input_count;
  index->stats->queries_filter_rand_output_count += index->stats->query_filter_rand_output_count;
    
    
  if(approximate_result.node != NULL) {  
    index->stats->query_approx_distance = sqrtf(approximate_result.distance);
    index->stats->query_approx_node_filename = approximate_result.node->filename;
    index->stats->query_approx_node_size = approximate_result.node->node_size;;
    index->stats->query_approx_node_level = approximate_result.node->level;
  }
 
  if (approximate_result.distance == 0) {
    index->stats->query_exact_distance = sqrtf(approximate_result.distance);
    index->stats->query_exact_node_filename = approximate_result.node->filename;
    index->stats->query_exact_node_size = approximate_result.node->node_size;;
    index->stats->query_exact_node_level = approximate_result.node->level; 
    update_query_stats(index,q_id, found_knn, approximate_result);
    get_query_stats(index, found_knn);
  }

  RESET_QUERY_COUNTERS()
    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
 
    struct query_result bsf_result = approximate_result;
    
   
  pqueue_t *pq = pqueue_init(index->first_node->node_size, 
			     cmp_pri, get_pri, set_pri, get_pos, set_pos);
	

  struct query_result * root_pq_item = malloc(sizeof(struct query_result));
  root_pq_item->node = index->first_node;
  root_pq_item->distance = calculate_node_min_distance (index, index->first_node, query_ts);
    

  pqueue_insert(pq, root_pq_item);    

  struct query_result * n;
  struct query_result temp;

  while ((n = pqueue_pop(pq)))
    {
      for (unsigned int pos = found_knn; pos < k; ++pos)
	{
	  bsf_result = knn_results[pos];
	
	  if (n->distance > bsf_result.distance/(1+epsilon)) //add epsilon+1
	    {
	      found_knn = pos+1;
	      COUNT_PARTIAL_TIME_END

		update_query_stats(index,q_id, found_knn, bsf_result);
	      get_query_stats(index, found_knn);
	      if (found_knn < k){	    
		bsf_result = knn_results[found_knn];
	      }
	  
	      RESET_QUERY_COUNTERS()
		RESET_PARTIAL_COUNTERS()
		COUNT_PARTIAL_TIME_START
		}

	}

      if (found_knn == k)
	{
	  break;
	}

      if (n->node->is_leaf) // n is a leaf
        {
 	  calculate_node_knn_distance(index, n->node, query_ts,query_ts_reordered,
				      query_order, offset, bsf_result.distance,
				      k,knn_results,NULL, NULL,&curr_size,q_settings.serial);
	  
        }
      else  //n is an internal node
        {
          temp = knn_results[k-1];
          kth_bsf =  temp.distance;
	  
          ts_type child_distance;
	  child_distance = calculate_node_min_distance(index, n->node->left_child,query_ts);

	  if ((child_distance < kth_bsf/(1+epsilon)) &&
	      (n->node->left_child != approximate_result.node)) 
	    {
	      struct query_result * mindist_result_left = malloc(sizeof(struct query_result));
	      mindist_result_left->node = n->node->left_child;
	      mindist_result_left->distance =  child_distance;
	      pqueue_insert(pq, mindist_result_left);
	    }
	  else if ((n->node->left_child != approximate_result.node))
	    {
	      print_pruning_snapshots(n->node->left_child,
				      kth_bsf,
				      child_distance,
				      found_knn+1,
				      q_id,
				      qfilename);
	    }
	  child_distance = calculate_node_min_distance(index, n->node->right_child,query_ts);

	  if ((child_distance < kth_bsf/(1+epsilon))  &&
	      (n->node->right_child != approximate_result.node)) 
	    {
	      struct query_result * mindist_result_right = malloc(sizeof(struct query_result));
	      mindist_result_right->node = n->node->right_child;
	      mindist_result_right->distance =  child_distance;
	      pqueue_insert(pq, mindist_result_right);
	    }
	  else if ((n->node->right_child != approximate_result.node))
	    {
	      print_pruning_snapshots(n->node->right_child,
				      kth_bsf,
				      child_distance,
				      found_knn+1,
				      q_id,
				      qfilename);
	    }
			     
        }
      free(n);
    }

  for (unsigned int pos = found_knn; pos < k; ++pos)
    {
      bsf_result = knn_results[pos];
      found_knn = pos+1;
      COUNT_PARTIAL_TIME_END
        update_query_stats(index,q_id, found_knn, bsf_result);
      get_query_stats(index, found_knn);
      RESET_QUERY_COUNTERS()
        RESET_PARTIAL_COUNTERS()
        COUNT_PARTIAL_TIME_START	
	}
    
  while ((n = pqueue_pop(pq)))
    {
      free(n);
    }

  pqueue_free(pq);
  free(knn_results);	
}

void dump_mindists (struct hercules_index *index,
		    struct hercules_node *node,
		    ts_type *query_ts)
{

  ts_type distance;
  ts_type QoS;
  
  distance = calculate_node_min_distance(index, node, query_ts);
  QoS = get_node_QoS(index,node);
  
  printf("%*s%lf\t%d\t%d\t%lf\n",
	 node->level,
	 "",
	 sqrtf(distance),
	 node->num_node_points,
	 node->level,
	 sqrtf(QoS)
	 );

  if(!node->is_leaf)
    {
      dump_mindists(index, node->left_child, query_ts);
      dump_mindists(index, node->right_child, query_ts);		
    }

}

	
void  exact_knn_search_track_bsf (ts_type *query_ts, ts_type * query_ts_reordered, int * query_order,
				  unsigned int offset, struct hercules_index *index,ts_type minimum_distance,
				  ts_type epsilon, ts_type delta, unsigned int k,
				  unsigned int q_id, char * qfilename, 
				  struct bsf_snapshot ** bsf_snapshots, unsigned int * cur_bsf_snapshot, query_settings q_settings)
{
    
  unsigned int curr_size = 0;
  ts_type bsf = FLT_MAX;
  ts_type kth_bsf = FLT_MAX;
  ts_type temp_bsf = FLT_MAX;
    
  unsigned int found_knn = 0;
  unsigned int last_found_knn = 0;
    
  struct query_result * knn_results = calloc(k,sizeof(struct query_result));
  for (int idx = 0; idx < k; ++idx)
    {
      knn_results[idx].node = NULL;
      knn_results[idx].distance = FLT_MAX;      
    }
    
  approximate_knn_search(query_ts, query_ts_reordered, query_order, offset, bsf, index,knn_results,k, bsf_snapshots, cur_bsf_snapshot,&curr_size,q_settings.serial);
    

  struct query_result approximate_result = knn_results[0];    

  COUNT_PARTIAL_TIME_END

    for (int idx = 0; idx < k; ++idx)
      {
	bsf_snapshots[idx][*cur_bsf_snapshot].distance = knn_results[idx].distance;
	bsf_snapshots[idx][*cur_bsf_snapshot].time = partial_time;
      }
  ++(*cur_bsf_snapshot);
      
  index->stats->query_filter_total_time  = partial_time;	
    
  index->stats->query_filter_input_time  = partial_input_time;
  index->stats->query_filter_output_time = partial_output_time;
  index->stats->query_filter_load_node_time = partial_load_node_time;
  index->stats->query_filter_cpu_time    = partial_time-partial_input_time-partial_output_time;
  index->stats->query_filter_seq_input_count   = partial_seq_input_count;
  index->stats->query_filter_seq_output_count  = partial_seq_output_count;
  index->stats->query_filter_rand_input_count  = partial_rand_input_count;
  index->stats->query_filter_rand_output_count = partial_rand_output_count;
    
  index->stats->query_filter_loaded_nodes_count = loaded_nodes_count;
  index->stats->query_filter_loaded_ts_count = loaded_ts_count;
  index->stats->query_filter_checked_nodes_count = checked_nodes_count;
  index->stats->query_filter_checked_ts_count = checked_ts_count;

  index->stats->query_approx_distance = sqrtf(approximate_result.distance);
  index->stats->query_approx_node_filename = approximate_result.node->filename;
  index->stats->query_approx_node_size = approximate_result.node->node_size;;
  index->stats->query_approx_node_level = approximate_result.node->level;
    
  index->stats->queries_filter_total_time  +=  index->stats->query_filter_total_time;
    
  index->stats->queries_filter_input_time  +=  index->stats->query_filter_input_time;
  index->stats->queries_filter_output_time +=  index->stats->query_filter_output_time;
  index->stats->queries_filter_load_node_time += index->stats->query_filter_load_node_time;
  index->stats->queries_filter_cpu_time    += index->stats->query_filter_cpu_time;

  index->stats->queries_filter_seq_input_count   += index->stats->query_filter_seq_input_count;
  index->stats->queries_filter_seq_output_count  += index->stats->query_filter_seq_output_count;
  index->stats->queries_filter_rand_input_count  += index->stats->query_filter_rand_input_count;
  index->stats->queries_filter_rand_output_count += index->stats->query_filter_rand_output_count;
    
    
  if(approximate_result.node != NULL) {  
    index->stats->query_approx_distance = sqrtf(approximate_result.distance);
    index->stats->query_approx_node_filename = approximate_result.node->filename;
    index->stats->query_approx_node_size = approximate_result.node->node_size;;
    index->stats->query_approx_node_level = approximate_result.node->level;
  }
 
  if (approximate_result.distance == 0) {
    index->stats->query_exact_distance = sqrtf(approximate_result.distance);
    index->stats->query_exact_node_filename = approximate_result.node->filename;
    index->stats->query_exact_node_size = approximate_result.node->node_size;;
    index->stats->query_exact_node_level = approximate_result.node->level; 
    update_query_stats(index,q_id, found_knn, approximate_result);
    get_query_stats(index, found_knn);
    print_query_stats(index, q_id, found_knn,qfilename);	     
  }


  RESET_QUERY_COUNTERS()
    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
 
    struct query_result bsf_result = approximate_result;
    
   
  pqueue_t *pq = pqueue_init(index->first_node->node_size, 
			     cmp_pri, get_pri, set_pri, get_pos, set_pos);
	
 
  struct query_result * root_pq_item = malloc(sizeof(struct query_result));
  root_pq_item->node = index->first_node;
  root_pq_item->distance = calculate_node_min_distance (index, index->first_node, query_ts);
    
  //initialize the lb distance to be the distance of the query to the root.

  pqueue_insert(pq, root_pq_item);    

  struct query_result * n;
  struct query_result temp;
    
  while ((n = pqueue_pop(pq)))
    {
      for (unsigned int pos = found_knn; pos < k; ++pos)
	{
	  bsf_result = knn_results[pos];
	
	  if (n->distance > bsf_result.distance/(1+epsilon)) //add epsilon+1
	    {
	      last_found_knn = found_knn;
	      found_knn = pos+1;
	      COUNT_PARTIAL_TIME_END

		update_query_stats(index,q_id, found_knn, bsf_result);
	      get_query_stats(index, found_knn);
	      print_query_stats(index, q_id, found_knn,qfilename);
	  
	      for (int idx = 0; idx < k; ++idx)
		{
		  bsf_snapshots[idx][*cur_bsf_snapshot].distance = knn_results[idx].distance;
		  bsf_snapshots[idx][*cur_bsf_snapshot].time = partial_time;
		}
	      ++(*cur_bsf_snapshot);
	  
	      if (found_knn < k){	    
		bsf_result = knn_results[found_knn];
	      }
	  
	      RESET_QUERY_COUNTERS()
		RESET_PARTIAL_COUNTERS()
		COUNT_PARTIAL_TIME_START
		}

	}

      if (found_knn == k)
	{
	  break;
	}
      if (n->node->is_leaf) // n is a leaf
        {
 	  calculate_node_knn_distance(index, n->node, query_ts,query_ts_reordered,
				      query_order, offset, bsf_result.distance,
				      k,knn_results,bsf_snapshots,cur_bsf_snapshot,&curr_size,q_settings.serial);
        }
      else  //n is an internal node
        {
          temp = knn_results[k-1];
          kth_bsf =  temp.distance;
	  
          ts_type child_distance;
	  child_distance = calculate_node_min_distance(index, n->node->left_child,query_ts);

	  if ((child_distance < kth_bsf/(1+epsilon)) &&
	      (n->node->left_child != approximate_result.node)) 
	    {
	      struct query_result * mindist_result_left = malloc(sizeof(struct query_result));
	      mindist_result_left->node = n->node->left_child;
	      mindist_result_left->distance =  child_distance;
	      pqueue_insert(pq, mindist_result_left);
	    }

	  child_distance = calculate_node_min_distance(index, n->node->right_child,query_ts);

	  if ((child_distance < kth_bsf/(1+epsilon))  &&
	      (n->node->right_child != approximate_result.node)) 	  
	    {
	      struct query_result * mindist_result_right = malloc(sizeof(struct query_result));
	      mindist_result_right->node = n->node->right_child;
	      mindist_result_right->distance =  child_distance;
	      pqueue_insert(pq, mindist_result_right);
	    }
        }
      free(n);
    }

  for (unsigned int pos = found_knn; pos < k; ++pos)
    {
      bsf_result = knn_results[pos];
      found_knn = pos+1;
      COUNT_PARTIAL_TIME_END
        update_query_stats(index,q_id, found_knn, bsf_result);
      get_query_stats(index, found_knn);
	
      print_query_stats(index, q_id, found_knn,qfilename);	
      RESET_QUERY_COUNTERS()
        RESET_PARTIAL_COUNTERS()
        COUNT_PARTIAL_TIME_START	
	}
  for (int idx = 0; idx < k; ++idx)
    {
      bsf_snapshots[idx][*cur_bsf_snapshot].distance = knn_results[idx].distance;
      bsf_snapshots[idx][*cur_bsf_snapshot].time = partial_time;
    }
  ++(*cur_bsf_snapshot);

  print_bsf_snapshots(index, q_id,k,qfilename,bsf_snapshots, *cur_bsf_snapshot);	    
      
  while ((n = pqueue_pop(pq)))
    {
      free(n);
    }
  pqueue_free(pq);
  free(knn_results);


}

void hercules_calc_tlb (ts_type *query_ts, struct hercules_index *index, struct hercules_node * curr_node) {
    
  ts_type curr_lb_dist = 0;
  ts_type curr_exact_dist = 0;
    
  if (curr_node == NULL)
    {
      return;
    }
  
  if (!curr_node->is_leaf)
    {		
      hercules_calc_tlb(query_ts,index,curr_node->left_child);
      hercules_calc_tlb(query_ts,index,curr_node->right_child);
    }
  else {
    curr_lb_dist = calculate_node_min_distance (index, curr_node, query_ts);
    if (curr_node->file_buffer->buffered_list_size == 0) 
      {	
	curr_node->file_buffer->buffered_list = get_all_time_series_in_node(index, curr_node,0);
	curr_node->file_buffer->buffered_list_size = curr_node->file_buffer->disk_count;

	if (curr_node->file_buffer->buffered_list == NULL)
	  {
	    fprintf(stderr,
		    "Error in hercules_index.c:  Could not retrieve all time series for node %s.\n",
		    curr_node->filename);
	  }	
      }
    total_ts_count = total_ts_count + curr_node->file_buffer->buffered_list_size;
    ++leaf_nodes_count;
    
    for (int idx = 0; idx < curr_node->file_buffer->buffered_list_size; ++idx)
      {  
	curr_exact_dist = ts_euclidean_distance(query_ts,
						curr_node->file_buffer->buffered_list[idx],
						index->settings->timeseries_size);
	total_tlb += sqrtf(curr_lb_dist / curr_exact_dist);		
      }	   
  }

}



void get_query_stats(struct hercules_index * index, unsigned int found_knn)
{

  if (total_ts_count != 0)
    {
      index->stats->query_pruning_ratio = 1.0 - ((double)index->stats->query_total_checked_ts_count/
						 index->stats->total_ts_count);
    }
    
  if (found_knn == 1)
    {
      if (index->stats->query_exact_distance != 0)
	{
	  index->stats->query_eff_epsilon =   (index->stats->query_approx_distance
					       -index->stats->query_exact_distance
					       )/ index->stats->query_exact_distance;	  
	}

    }

}


void print_bsf_snapshots(struct hercules_index * index, unsigned int query_num, unsigned int k,
			 char * queries, struct bsf_snapshot **bsf_snapshots, unsigned int cur_bsf_snapshot)
{

  for (unsigned int i = 0; i < k; ++i)
    {
      for (unsigned int j = 0; j < cur_bsf_snapshot; ++j)
	{
	  printf("Query_bsf_snapshot_time_secs\t%lf\t%s\t%u\t%u\n",
		 bsf_snapshots[i][j].time/1000000,	        
		 queries,
		 query_num,
		 i+1
		 );
	  printf("Query_bsf_snapshot_distance\t%lf\t%s\t%u\t%u\n",
		 sqrtf(bsf_snapshots[i][j].distance), 
		 queries,
		 query_num,
		 i+1
		 );      
	}
    }
}

void print_pruning_snapshots(struct hercules_node * node,
			     ts_type node_bsf,
			     ts_type node_mindist,
			     unsigned int k,
			     unsigned int query_num,
			     char * queries)			     
{
  
  printf("Query_pruning_snapshot_node_filename\t%s\t%s\t%u\t%u\n",
	 node->filename,	        
	 queries,
	 query_num,
	 k
	 );

  printf("Query_pruning_snapshot_node_level\t%u\t%s\t%u\t%u\n",
	 node->level,	        
	 queries,
	 query_num,
	 k
	 );  
  
  printf("Query_pruning_snapshot_node_bsf\t%lf\t%s\t%u\t%u\n",
	 node_bsf,	        
	 queries,
	 query_num,
	 k
	 );

  printf("Query_pruning_snapshot_node_mindist\t%lf\t%s\t%u\t%u\n",
	 node_mindist,	        
	 queries,
	 query_num,
	 k
	 );  
  

}


void print_query_stats(struct hercules_index * index, unsigned int query_num, unsigned int found_knn, char * queries)
{

  if (found_knn == 1)
    {
      printf("Query_filter_input_time_secs\t%lf\t%s\t%u\t%u\n",
	     index->stats->query_filter_input_time/1000000,	       
	     queries,
	     query_num,
	     found_knn
	     );

      printf("Query_filter_output_time_secs\t%lf\t%s\t%u\t%u\n",
	     index->stats->query_filter_output_time/1000000,	       
	     queries,
	     query_num,
	     found_knn
	     );

      printf("Query_filter_load_node_time_secs\t%lf\t%s\t%u\t%u\n",
	     index->stats->query_filter_load_node_time/1000000,	       
	     queries,
	     query_num,
	     found_knn
	     );

      printf("Query_filter_cpu_time_secs\t%lf\t%s\t%u\t%u\n",
	     index->stats->query_filter_cpu_time/1000000,	       
	     queries,
	     query_num,
	     found_knn
	     );
	
      printf("Query_filter_total_time_secs\t%lf\t%s\t%u\t%u\n",
	     index->stats->query_filter_total_time/1000000,	       
	     queries,	       
	     query_num,
	     found_knn
	     );

      printf("Query_filter_seq_input_count\t%llu\t%s\t%u\t%u\n",
	     index->stats->query_filter_seq_input_count,
	     queries,
	     query_num,
	     found_knn
	     );
	
      printf("Query_filter_seq_output_count\t%llu\t%s\t%u\t%u\n",
	     index->stats->query_filter_seq_output_count,
	     queries,
	     query_num,
	     found_knn
	     );
	
      printf("Query_filter_rand_input_count\t%llu\t%s\t%u\t%u\n",
	     index->stats->query_filter_rand_input_count,
	     queries,
	     query_num,
	     found_knn	       
	     );
	
      printf("Query_filter_rand_output_count\t%llu\t%s\t%u\t%u\n",
	     index->stats->query_filter_rand_output_count,
	     queries,
	     query_num,
	     found_knn	       
	     );
      printf("Query_filter_checked_nodes_count\t%u\t%s\t%u\t%u\n",
	     index->stats->query_filter_checked_nodes_count,
	     queries,
	     query_num,
	     found_knn	       
	     );	

      printf("Query_filter_checked_ts_count\t%llu\t%s\t%u\t%u\n",
	     index->stats->query_filter_checked_ts_count,
	     queries,
	     query_num,
	     found_knn	       
	     );	

      printf("Query_filter_loaded_nodes_count\t%u\t%s\t%u\t%u\n",
	     index->stats->query_filter_loaded_nodes_count,
	     queries,
	     query_num,
	     found_knn	       
	     );	

      printf("Query_filter_loaded_ts_count\t%llu\t%s\t%u\t%u\n",
	     index->stats->query_filter_loaded_ts_count,
	     queries,
	     query_num,
	     found_knn	       
	     );
      printf("Query_approx_distance\t%f\t%s\t%u\t%u\n",
	     index->stats->query_approx_distance,
	     queries,
	     query_num,
	     found_knn	       
	     );	

      printf("Query_approx_node_filename\t%s\t%s\t%u\t%u\n",
	     index->stats->query_approx_node_filename,
	     queries,
	     query_num,
	     found_knn	       
	     );
	
      printf("Query_approx_node_size\t%u\t%s\t%u\t%u\n",
	     index->stats->query_approx_node_size,
	     queries,
	     query_num,
	     found_knn	       
	     );	

      printf("Query_approx_node_level\t%u\t%s\t%u\t%u\n",
	     index->stats->query_approx_node_level,
	     queries,
	     query_num,
	     found_knn
	     );

      printf("Query_eff_epsilon\t%f\t%s\t%u\t%u\n",
	     index->stats->query_eff_epsilon,	       
	     queries,
	     query_num,
	     found_knn	       
	     );		

	
    }

  printf("Query_refine_input_time_secs\t%lf\t%s\t%u\t%u\n",
	 index->stats->query_refine_input_time/1000000,	       
	 queries,
	 query_num,
	 found_knn	       
	 );

  printf("Query_refine_output_time_secs\t%lf\t%s\t%u\t%u\n",
	 index->stats->query_refine_output_time/1000000,	       
	 queries,
	 query_num,
	 found_knn	       
	 );

  printf("Query_refine_load_node_time_secs\t%lf\t%s\t%u\t%u\n",
	 index->stats->query_refine_load_node_time/1000000,	       
	 queries,
	 query_num,
	 found_knn	       
	 );

  printf("Query_refine_cpu_time_secs\t%lf\t%s\t%u\t%u\n",
	 index->stats->query_refine_cpu_time/1000000,	       
	 queries,
	 query_num,
	 found_knn	       
	 );
	
  printf("Query_refine_total_time_secs\t%lf\t%s\t%u\t%u\n",
	 index->stats->query_refine_total_time/1000000,	       
	 queries,
	 query_num,
	 found_knn	       
	 );

  printf("Query_refine_seq_input_count\t%llu\t%s\t%u\t%u\n",
	 index->stats->query_refine_seq_input_count,
	 queries,
	 query_num,
	 found_knn	       
	 );
	
  printf("Query_refine_seq_output_count\t%llu\t%s\t%u\t%u\n",
	 index->stats->query_refine_seq_output_count,
	 queries,
	 query_num,
	 found_knn	       
	 );
	
  printf("Query_refine_rand_input_count\t%llu\t%s\t%u\t%u\n",
	 index->stats->query_refine_rand_input_count,
	 queries,
	 query_num,
	 found_knn	       
	 );
	
  printf("Query_refine_rand_output_count\t%llu\t%s\t%u\t%u\n",
	 index->stats->query_refine_rand_output_count,
	 queries,
	 query_num,
	 found_knn	       
	 );

	
  printf("Query_refine_checked_nodes_count\t%u\t%s\t%u\t%u\n",
	 index->stats->query_refine_checked_nodes_count,
	 queries,
	 query_num,
	 found_knn	       
	 );	

  printf("Query_refine_checked_ts_count\t%llu\t%s\t%u\t%u\n",
	 index->stats->query_refine_checked_ts_count,
	 queries,
	 query_num,
	 found_knn	       
	 );	

  printf("Query_refine_loaded_nodes_count\t%u\t%s\t%u\t%u\n",
	 index->stats->query_refine_loaded_nodes_count,
	 queries,
	 query_num,
	 found_knn	       
	 );	

  printf("Query_refine_loaded_ts_count\t%llu\t%s\t%u\t%u\n",
	 index->stats->query_refine_loaded_ts_count,
	 queries,
	 query_num,
	 found_knn	       
	 );

  printf("Query_total_input_time_secs\t%lf\t%s\t%u\t%u\n",
	 index->stats->query_total_input_time/1000000,	       
	 queries,
	 query_num,
	 found_knn	       
	 );

  printf("Query_total_output_time_secs\t%lf\t%s\t%u\t%u\n",
	 index->stats->query_total_output_time/1000000,	       
	 queries,
	 query_num,
	 found_knn	       
	 );

  printf("Query_total_load_node_time_secs\t%lf\t%s\t%u\t%u\n",
	 index->stats->query_total_load_node_time/1000000,	       
	 queries,
	 query_num,
	 found_knn	       
	 );

  printf("Query_total_cpu_time_secs\t%lf\t%s\t%u\t%u\n",
	 index->stats->query_total_cpu_time/1000000,	       
	 queries,
	 query_num,
	 found_knn	       
	 );
	
  printf("Query_total_time_secs\t%lf\t%s\t%u\t%u\n",
	 index->stats->query_total_time/1000000,	       
	 queries,
	 query_num,
	 found_knn	       
	 );

  printf("Query_total_seq_input_count\t%llu\t%s\t%u\t%u\n",
	 index->stats->query_total_seq_input_count,
	 queries,
	 query_num,
	 found_knn	       
	 );
	
  printf("Query_total_seq_output_count\t%llu\t%s\t%u\t%u\n",
	 index->stats->query_total_seq_output_count,
	 queries,
	 query_num,
	 found_knn	       
	 );
	
  printf("Query_total_rand_input_count\t%llu\t%s\t%u\t%u\n",
	 index->stats->query_total_rand_input_count,
	 queries,
	 query_num,
	 found_knn	       
	 );
	
  printf("Query_total_rand_output_count\t%llu\t%s\t%u\t%u\n",
	 index->stats->query_total_rand_output_count,
	 queries,
	 query_num,
	 found_knn	       
	 );


  printf("Query_total_checked_nodes_count\t%u\t%s\t%u\t%u\n",
	 index->stats->query_total_checked_nodes_count,
	 queries,
	 query_num,
	 found_knn	       
	 );	

  printf("Query_total_checked_ts_count\t%llu\t%s\t%u\t%u\n",
	 index->stats->query_total_checked_ts_count,
	 queries,
	 query_num,
	 found_knn	       
	 );	

  printf("Query_total_loaded_nodes_count\t%u\t%s\t%u\t%u\n",
	 index->stats->query_total_loaded_nodes_count,
	 queries,
	 query_num,
	 found_knn	       
	 );	

  printf("Query_total_loaded_ts_count\t%llu\t%s\t%u\t%u\n",
	 index->stats->query_total_loaded_ts_count,
	 queries,
	 query_num,
	 found_knn	       
	 );

  printf("Query_exact_distance\t%f\t%s\t%u\t%u\n",
	 index->stats->query_exact_distance,
	 queries,
	 query_num,
	 found_knn	       
	 );	

  printf("Query_exact_node_filename\t%s\t%s\t%u\t%u\n",
	 index->stats->query_exact_node_filename,
	 queries,
	 query_num,
	 found_knn	       
	 );
	
  printf("Query_exact_node_size\t%u\t%s\t%u\t%u\n",
	 index->stats->query_exact_node_size,
	 queries,
	 query_num,
	 found_knn	       
	 );	

  printf("Query_exact_node_level\t%u\t%s\t%u\t%u\n",
	 index->stats->query_exact_node_level,
	 queries,
	 query_num,
	 found_knn	       
	 );


  printf("Query_exact_node_file_pos\t%u\t%s\t%u\t%u\n",
	 index->stats->query_exact_node_file_pos,
	 queries,
	 query_num,
	 found_knn	       
	 );
	
	
  printf("Query_pruning_ratio_level\t%f\t%s\t%u\t%u\n",
	 index->stats->query_pruning_ratio,	       
	 queries,
	 query_num,
	 found_knn	       
	 );

}

void update_query_stats(struct hercules_index * index,
			unsigned int query_id,
			unsigned int found_knn,
			struct query_result bsf_result)
{
  index->stats->query_refine_total_time  = partial_time;	
        
  index->stats->query_refine_input_time  = partial_input_time;
  index->stats->query_refine_output_time = partial_output_time;
  index->stats->query_refine_load_node_time = partial_load_node_time;    
  index->stats->query_refine_cpu_time    = partial_time
    - partial_input_time
    - partial_output_time;
  index->stats->query_refine_seq_input_count   = partial_seq_input_count;
  index->stats->query_refine_seq_output_count  = partial_seq_output_count;
  index->stats->query_refine_rand_input_count  = partial_rand_input_count;
  index->stats->query_refine_rand_output_count = partial_rand_output_count;

  index->stats->query_refine_loaded_nodes_count = loaded_nodes_count;
  index->stats->query_refine_loaded_ts_count = loaded_ts_count;
  index->stats->query_refine_checked_nodes_count = checked_nodes_count;
  index->stats->query_refine_checked_ts_count = checked_ts_count;
  
  if (found_knn == 1)
    {
      index->stats->query_total_time  = partial_time
	+  index->stats->query_filter_total_time;
      index->stats->query_total_input_time  = partial_input_time
	+  index->stats->query_filter_input_time;
      index->stats->query_total_output_time  = partial_output_time
	+  index->stats->query_filter_output_time;
      index->stats->query_total_load_node_time  = partial_load_node_time
	+  index->stats->query_filter_load_node_time;
      index->stats->query_total_cpu_time  =  index->stats->query_total_time
	-  index->stats->query_total_input_time
	-  index->stats->query_total_output_time;
    
      index->stats->query_total_seq_input_count   = partial_seq_input_count
	+ index->stats->query_filter_seq_input_count;
      index->stats->query_total_seq_output_count   = partial_seq_output_count
	+ index->stats->query_filter_seq_output_count;
      index->stats->query_total_rand_input_count   = partial_rand_input_count
	+ index->stats->query_filter_rand_input_count;
      index->stats->query_total_rand_output_count   = partial_rand_output_count
	+ index->stats->query_filter_rand_output_count;

      index->stats->query_total_loaded_nodes_count = loaded_nodes_count 
	+ index->stats->query_filter_loaded_nodes_count;
      index->stats->query_total_loaded_ts_count = loaded_ts_count
	+ index->stats->query_filter_loaded_ts_count;
      index->stats->query_total_checked_nodes_count = checked_nodes_count
	+ index->stats->query_filter_checked_nodes_count;
      index->stats->query_total_checked_ts_count = checked_ts_count
	+ index->stats->query_filter_checked_ts_count;      
    }
  else
    {
      index->stats->query_total_time  = partial_time;
      index->stats->query_total_input_time  = partial_input_time;
      index->stats->query_total_output_time  = partial_output_time;
      index->stats->query_total_load_node_time  = partial_load_node_time;

      index->stats->query_total_cpu_time  =  index->stats->query_total_time
	-  index->stats->query_total_input_time
	-  index->stats->query_total_output_time;
    
      index->stats->query_total_seq_input_count   = partial_seq_input_count;
      index->stats->query_total_seq_output_count   = partial_seq_output_count;
      index->stats->query_total_rand_input_count   = partial_rand_input_count;
      index->stats->query_total_rand_output_count   = partial_rand_output_count;

      index->stats->query_total_loaded_nodes_count = loaded_nodes_count; 
      index->stats->query_total_loaded_ts_count = loaded_ts_count;
      index->stats->query_total_checked_nodes_count = checked_nodes_count;
      index->stats->query_total_checked_ts_count = checked_ts_count;
      
    }
  
  index->stats->query_exact_distance = sqrtf(bsf_result.distance);

  if (bsf_result.node != NULL)
    {
      index->stats->query_exact_node_filename = bsf_result.node->filename;
      index->stats->query_exact_node_size = bsf_result.node->node_size;;
      index->stats->query_exact_node_level = bsf_result.node->level;
      index->stats->query_exact_node_file_pos = bsf_result.node->file_pos;
    }
}



ts_type get_node_QoS(struct hercules_index * index, struct hercules_node * node)
{

  ts_type node_range_value = 0;
  for (int i=0; i <  node->num_node_points;++i)
    { 
      struct segment_sketch curr_node_segment_sketch = node->node_segment_sketches[i];

      //This is the QoS of this segment. QoS is the estimation quality evaluated as =
      //QoS = segment_length * (max_mean-min_mean) * ((max_mean-min_mean) +
      //     (max_stdev * max_stdev))
      //The smaller the QoS, the more effective the bounds are for similarity
      //estimation
	
      node_range_value += range_calc(curr_node_segment_sketch,
				     get_segment_length(node->node_points, i));

    }

  return node_range_value;
  
}

int queue_bounded_sorted_insert(struct  query_result *q, struct query_result d, unsigned int *cur_size, unsigned int k)
{
  struct query_result  temp;
  size_t i;
  size_t newsize;

  bool is_duplicate = false;
  for (unsigned int itr = 0 ; itr < *cur_size ; ++itr)
    {
      if (q[itr].distance == d.distance)
	is_duplicate = true;
    }
   
  if (!is_duplicate)
    {
    
      /* the queue is full, overwrite last element*/
      if (*cur_size == k) {      
	q[k-1].distance = d.distance;
	q[k-1].node = d.node;
      }
      else
	{
	  q[*cur_size].distance = d.distance;
	  q[*cur_size].node = d.node;      
	  ++(*cur_size);
	}

      unsigned int idx,j;

      idx = 1;
    
      while (idx < *cur_size)
	{
	  j = idx;
	  while ( j > 0 &&  ( (q[j-1]).distance > q[j].distance)) 
	    {
	      temp = q[j];
	      q[j].distance = q[j-1].distance;
	      q[j].node = q[j-1].node;	
	      q[j-1].distance = temp.distance;
	      q[j-1].node = temp.node;		
	      --j;
	    }
	  ++idx;
	}
    }
  return 0;
}

