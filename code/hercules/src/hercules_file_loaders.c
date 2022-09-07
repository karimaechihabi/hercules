//
//  hercules_file_loaders.c
//
//  Created by Karima Echihabi on 18/12/2018
//

#include "../config.h"
#include "../globals.h"
#include <stdio.h>
#include "../include/hercules_file_loaders.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include "../include/hercules_query_engine.h"
#include <pthread.h>


enum response hercules_query_binary_file(const char *ifilename, int q_num, struct hercules_index *index,
					 float minimum_distance, ts_type epsilon,
					 ts_type r_delta )
{
  RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START

    FILE * ifile;
  COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_INPUT_TIME_START    
    ifile = fopen (ifilename,"rb");
  COUNT_PARTIAL_INPUT_TIME_END    
    if (ifile == NULL) {
      fprintf(stderr, "File %s not found!\n", ifilename);
      exit(-1);
    }
    
  COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_INPUT_TIME_START
    fseek(ifile, 0L, SEEK_END);
  file_position_type sz = (file_position_type) ftell(ifile);
  fseek(ifile, 0L, SEEK_SET);
  COUNT_PARTIAL_INPUT_TIME_END        
    unsigned int ts_length = index->settings->timeseries_size;
  file_position_type total_records = sz/ts_length * sizeof(ts_type);
  unsigned int offset = 0;

  if (total_records < q_num) {
    fprintf(stderr, "File %s has only %llu records!\n", ifilename, total_records);
    exit(-1);
  }

  unsigned int q_loaded = 0;
  ts_type * query_ts = calloc(1,sizeof(ts_type) * ts_length);
  ts_type * query_ts_reordered = calloc(1,sizeof(ts_type) * ts_length);
  int * query_order = calloc(1,sizeof(int) * ts_length);
  if( query_order == NULL)
    return FAILURE;
  while (q_loaded < q_num)
    {
  
      RESET_QUERY_COUNTERS ()
	  
        COUNT_PARTIAL_SEQ_INPUT      
        COUNT_PARTIAL_INPUT_TIME_START
	fread(query_ts, sizeof(ts_type), ts_length, ifile);
      COUNT_PARTIAL_INPUT_TIME_END


	reorder_query(query_ts,query_ts_reordered,query_order,ts_length);        

      
      struct query_result result = exact_search(query_ts, query_ts_reordered, query_order, offset, index, minimum_distance, epsilon, r_delta);	

      q_loaded++;
 
      get_query_stats(index,1);	
      print_query_stats(index,q_loaded, 1, ifilename);
      RESET_PARTIAL_COUNTERS()
	COUNT_PARTIAL_TIME_START
	  
	}
  COUNT_PARTIAL_TIME_END
    RESET_PARTIAL_COUNTERS()
      
    free(query_ts);
  free(query_ts_reordered);
  free(query_order);

  if(fclose(ifile))
    {   
      fprintf(stderr, "Error in hercules_file_loaders.c: Could not close the query filename %s", ifilename);
      return FAILURE;
    }

  return SUCCESS;
}

enum response hercules_knn_query_binary_file(const char *ifilename,
					     int q_num,
					     struct hercules_index *index,
					     float minimum_distance,
					     ts_type epsilon,
					     ts_type r_delta,
					     unsigned int k,
					     boolean track_bsf,
					     boolean track_pruning,
					     boolean all_mindists,
					     boolean max_policy,
					     unsigned int nprobes,
					     unsigned char incremental,					   
					     query_settings q_settings,
					     int q_skip)					   
{

  struct bsf_snapshot ** bsf_snapshots = NULL;
  unsigned int max_bsf_snapshots;
  unsigned int cur_bsf_snapshot;
  if(track_bsf)
    {
      max_bsf_snapshots = 10000;
      cur_bsf_snapshot = 0;

      bsf_snapshots = calloc(k, sizeof(struct bsf_snapshot*));      
      for (unsigned int i = 0; i < k; ++i)
	{
	  bsf_snapshots[i] = calloc(max_bsf_snapshots, sizeof(struct bsf_snapshot));
	  for (unsigned int j = 0; j < max_bsf_snapshots; ++j)
	    {
	      bsf_snapshots[i][j].distance = FLT_MAX;
	      bsf_snapshots[i][j].time = FLT_MAX;
	    }      
	}
    }

  
  RESET_PARTIAL_COUNTERS()

    COUNT_PARTIAL_TIME_START

    FILE * ifile;
  COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_INPUT_TIME_START    
    ifile = fopen (ifilename,"rb");
  COUNT_PARTIAL_INPUT_TIME_END    
    if (ifile == NULL) {
      fprintf(stderr, "File %s not found!\n", ifilename);
      exit(-1);
    }
    
  COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_INPUT_TIME_START
    fseek(ifile, 0L, SEEK_END);
  file_position_type sz = (file_position_type) ftell(ifile);
  fseek(ifile, 0L, SEEK_SET);
  COUNT_PARTIAL_INPUT_TIME_END        
    unsigned int ts_length = index->settings->timeseries_size;
  file_position_type total_records = sz/ts_length * sizeof(ts_type);
  unsigned int offset = 0;

  if (total_records < q_num) {
    fprintf(stderr, "File %s has only %llu records!\n", ifilename, total_records);
    exit(-1);
  }


    
  ts_type * query_ts = calloc(1,sizeof(ts_type) * ts_length);
  ts_type * query_ts_reordered = calloc(1,sizeof(ts_type) * ts_length);
  int * query_order = calloc(1,sizeof(int) * ts_length);
  if( query_order == NULL)
    return FAILURE;
  ts_type * query_paa = NULL;
  query_paa = calloc(1,sizeof(ts_type) * index->settings->paa_segments);

    
  fseek(ifile,
	q_skip * ts_length * sizeof(ts_type),
	SEEK_SET);
	
  unsigned int q_loaded = q_skip;
  q_num = q_num + q_skip;
    
  while (q_loaded < q_num)
    {
  
      RESET_QUERY_COUNTERS ()
	  
        COUNT_PARTIAL_SEQ_INPUT      
        COUNT_PARTIAL_INPUT_TIME_START
	fread(query_ts, sizeof(ts_type), ts_length, ifile);
      COUNT_PARTIAL_INPUT_TIME_END


	reorder_query(query_ts,query_ts_reordered,query_order,ts_length);        

      q_loaded++;
      if (track_bsf)
	{
	  cur_bsf_snapshot = 0;
	  exact_knn_search_track_bsf(query_ts, query_ts_reordered, query_order, offset,
				     index, minimum_distance, epsilon, r_delta, k, q_loaded, ifilename,
				     bsf_snapshots, &cur_bsf_snapshot,q_settings);
	  
	  for (unsigned int i = 0; i < k; ++i)
	    {
	      for (unsigned int j = 0; j < max_bsf_snapshots; ++j)
		{
		  bsf_snapshots[i][j].distance = FLT_MAX;
		  bsf_snapshots[i][j].time = FLT_MAX;	      
		}      
	    }
	}
      else if (track_pruning)
	{
	  exact_knn_search_track_pruning(query_ts, query_ts_reordered, query_order, offset,
					 index, minimum_distance, epsilon, r_delta, k, q_loaded, ifilename,q_settings);
	  if (all_mindists)
	    {
	      printf("\nmindist\tnum_segments\tlevel\tQoS\n");
	      dump_mindists(index,index->first_node, query_ts);
	    }	    
	}
      else if (max_policy)
	{
	  exact_knn_search_max_policy(query_ts, query_ts_reordered, query_order, offset,
				      index, minimum_distance, epsilon, r_delta, k, q_loaded, ifilename,q_settings);		  
	}	
      else if (nprobes)
	{
	  exact_ng_knn_search(query_ts, query_ts_reordered, query_order, offset,
			      index, minimum_distance, k,
			      q_loaded, ifilename, nprobes,q_settings);		  
	}
      else 
	{
       
	  if (q_settings.num_threads == 1)
	    {
	      exact_de_knn_search(query_ts, query_ts_reordered, query_order, offset,
				  index, minimum_distance, epsilon, r_delta, k,
				  q_loaded, ifilename, q_settings);
	    }
	  else
	    {
	      if (q_settings.sims)
		{
		  paa_from_ts(query_ts, query_paa, index->settings->paa_segments, 
			      index->settings->ts_values_per_paa_segment);	 
 
		  exact_de_knn_search_psq_lf_sims_adaptive(query_ts, query_ts_reordered, query_order, offset, query_paa,
							   index, minimum_distance, epsilon, r_delta, k,
							   q_loaded, ifilename, q_settings);	      	      
		}
	      else
		{
		  exact_de_knn_search_psq_lf(query_ts, query_ts_reordered, query_order, offset,
					     index, minimum_distance, epsilon, r_delta, k,
					     q_loaded, ifilename, q_settings);	      	      		
		}

	    }
	  
	}
  
      RESET_PARTIAL_COUNTERS()
	COUNT_PARTIAL_TIME_START
	  
	}
  COUNT_PARTIAL_TIME_END
    RESET_PARTIAL_COUNTERS()
      
    free(query_paa);
  free(query_ts);
  free(query_ts_reordered);
  free(query_order);

  if(fclose(ifile))
    {   
      fprintf(stderr, "Error in hercules_file_loaders.c: Could not close the query filename %s", ifilename);
      return FAILURE;
    }

  if (track_bsf)
    {
      for (unsigned int i = 0; i < k; ++i)
	{
	  free(bsf_snapshots[i]);
	}
      free(bsf_snapshots);
    }
  return SUCCESS;
}







enum response hercules_index_insert_parallel(struct hercules_index *index, struct hercules_node *node,  ts_type * timeseries, void *thread_data)
{

  struct hercules_node *subtree = NULL;
  boolean lock_subtree =false;
  boolean lock_node =false;
  
  struct segment_sketch timeseries_segment_sketch;
  timeseries_segment_sketch.indicators = NULL;
  timeseries_segment_sketch.indicators = malloc (sizeof(ts_type)*2);
  if(timeseries_segment_sketch.indicators == NULL) {
    fprintf(stderr,"Error in hercules_index.c: Could not allocate memory for" 
	    "series segment sketch indicators.\n");
  }

  timeseries_segment_sketch.num_indicators = 2;
  double tS = 0;
  double tE = 0;
  struct timeval start_time;
  struct timeval end_time;
#if DETAILED_STATS == 1 
  gettimeofday(&start_time, NULL);
#endif

 
  while (!node->is_leaf )//& !node->is_splitting)
    {
      if (node_split_policy_route_to_left(node,timeseries,&timeseries_segment_sketch))
	node = node->left_child;
      else
	node = node->right_child;
    }


  pthread_mutex_lock(&(node->lock_data));
  
  while (!node->is_leaf )//& !node->is_splitting)
    {
      pthread_mutex_unlock(&(node->lock_data));
      while (!node->is_leaf )//& !node->is_splitting)
	{
	  if (node_split_policy_route_to_left(node,timeseries,&timeseries_segment_sketch))
	    node = node->left_child;
	  else
	    node = node->right_child;
	}
      pthread_mutex_lock(&(node->lock_data));
    }

#if DETAILED_STATS == 1
  gettimeofday(&end_time, NULL);
  tS = start_time.tv_sec*1000000 + (start_time.tv_usec);
  tE = end_time.tv_sec*1000000 + (end_time.tv_usec);
  ((index_buffer_data*)thread_data)->thread_traverse_total_time += (tE - tS);
  gettimeofday(&start_time, NULL); 
#endif	      

  if(!update_node_statistics_parallel(node, timeseries,&timeseries_segment_sketch))
    {
      fprintf(stderr,"Error in hercules_index.c: could not update \
                        statistics at node %s\n", node->filename);
      return FAILURE;
    }

  if(!append_ts_to_node_parallel(index,node, timeseries, thread_data))
    {
      fprintf(stderr,"Error in hercules_index.c: could not append \
                        time series to node %s\n", node->filename);
      return FAILURE;
    }
    
#if DETAILED_STATS == 1
  gettimeofday(&end_time, NULL);
  tS = start_time.tv_sec*1000000 + (start_time.tv_usec);
  tE = end_time.tv_sec*1000000 + (end_time.tv_usec);
  ((index_buffer_data*)thread_data)->thread_append_total_time += (tE - tS);
  gettimeofday(&start_time, NULL); 
#endif

  if (node->node_size >= index->settings->max_leaf_size)
    {
      split_node_detailed(index,node, &timeseries_segment_sketch,thread_data);
      node->is_leaf = 0;
    } 

  pthread_mutex_unlock(&(node->lock_data));

#if DETAILED_STATS == 1
  gettimeofday(&end_time, NULL);
  tS = start_time.tv_sec*1000000 + (start_time.tv_usec);
  tE = end_time.tv_sec*1000000 + (end_time.tv_usec);
  ((index_buffer_data*)thread_data)->thread_split_total_time += (tE - tS);    
  gettimeofday(&start_time, NULL); 
#endif

  free(timeseries_segment_sketch.indicators);    
  return SUCCESS;
}

enum response hercules_index_binary_file(const char *ifilename, file_position_type ts_num, struct hercules_index *index)
{
  double parse_time = 0;
  int ts_length = index->settings->timeseries_size;
  ts_type * ts = NULL;
    
  ts = NULL;
  ts = malloc (sizeof(ts_type) * ts_length); 

  if(ts == NULL)
    {
      fprintf(stderr,"Error in hercules_file_loaders.c: Could not allocate memory for ts.\n");
      return FAILURE;	
    }
    
  FILE * ifile; 
  COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_INPUT_TIME_START
    ifile = fopen (ifilename,"rb");
  COUNT_PARTIAL_INPUT_TIME_END
    if (ifile == NULL) {
      fprintf(stderr, "Error in hercules_file_loaders.c: File %s not found!\n", ifilename);
      return FAILURE;
    }
  COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_INPUT_TIME_START    
    fseek(ifile, 0L, SEEK_END);
  file_position_type sz = (file_position_type) ftell(ifile);
  COUNT_PARTIAL_INPUT_TIME_END    
    file_position_type total_records = sz/index->settings->timeseries_size * sizeof(ts_type);

  COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_INPUT_TIME_START
    fseek(ifile, 0L, SEEK_SET);
  COUNT_PARTIAL_INPUT_TIME_END        
    if (total_records < ts_num) {
      fprintf(stderr, "File %s has only %llu records!\n", ifilename, total_records);
      return FAILURE;
    }
    	
  file_position_type ts_loaded = 0;    
    
  int percentage = 100;
  if (percentage > 100)
    {
      percentage = (int) (ts_num / (file_position_type) 100);
    }

  while (ts_loaded<ts_num)
    {
#ifndef DEBUG
#if VERBOSE_LEVEL == 2
      printf("\r\x1b[32mLoading: \x1b[36m%d\x1b[0m",(ts_loaded + 1));
#endif
#endif
      COUNT_PARTIAL_SEQ_INPUT
	COUNT_PARTIAL_INPUT_TIME_START
        fread(ts, sizeof(ts_type), index->settings->timeseries_size, ifile);
      COUNT_PARTIAL_INPUT_TIME_END

	if (!hercules_index_insert(index,ts))
	  {
	    fprintf(stderr, "Error in hercules_file_loaders.c:  Could not \
                           add the time series to the index.\n");
	    return FAILURE;              
	  }
      COUNT_PARTIAL_TIME_END

	index->stats->idx_building_total_time  += partial_time;	
      index->stats->idx_building_input_time  += partial_input_time;
      index->stats->idx_building_output_time += partial_output_time;

      index->stats->idx_building_seq_input_count  += partial_seq_input_count;
      index->stats->idx_building_seq_output_count += partial_seq_output_count;
      index->stats->idx_building_rand_input_count  += partial_rand_input_count;
      index->stats->idx_building_rand_output_count += partial_rand_output_count;	    
 
      RESET_PARTIAL_COUNTERS()
	COUNT_PARTIAL_TIME_START

        ts_loaded++;

     
      if(ts_loaded % percentage == 0)
        {
	  float distance = 0;
        }

    }

  free(ts);
  COUNT_PARTIAL_INPUT_TIME_START    
    if(fclose(ifile))
      {   
        fprintf(stderr, "Error in hercules_file_loaders.c: Could not close the filename %s", ifilename);
        return FAILURE;
      }
  COUNT_PARTIAL_INPUT_TIME_END    

    COUNT_PARTIAL_TIME_END
    index->stats->idx_building_total_time  += partial_time;	
  index->stats->idx_building_input_time  += partial_input_time;
  index->stats->idx_building_output_time += partial_output_time;
  index->stats->idx_building_seq_input_count  += partial_seq_input_count;
  index->stats->idx_building_seq_output_count += partial_seq_output_count;
  index->stats->idx_building_rand_input_count  += partial_rand_input_count;
  index->stats->idx_building_rand_output_count += partial_rand_output_count;	    
 
  RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START

    return SUCCESS;      
}

enum response reorder_query(ts_type * query_ts, ts_type * query_ts_reordered, int * query_order, int ts_length)
{
  
  q_index *q_tmp = malloc(sizeof(q_index) * ts_length);
  int i;
	
  if( q_tmp == NULL )
    return FAILURE;

  for( i = 0 ; i < ts_length ; i++ )
    {
      q_tmp[i].value = query_ts[i];
      q_tmp[i].index = i;
    }
	
  qsort(q_tmp, ts_length, sizeof(q_index),znorm_comp);

  for( i=0; i<ts_length; i++)
    {
          
      query_ts_reordered[i] = q_tmp[i].value;
      query_order[i] = q_tmp[i].index;
    }
  free(q_tmp);

  return SUCCESS;
}


int znorm_comp(const void *a, const void* b)
{
  q_index* x = (q_index*)a;
  q_index* y = (q_index*)b;


  if (fabsf(y->value) > fabsf(x->value) )
    return 1;
  else if (fabsf(y->value) == fabsf(x->value))
    return 0;
  else
    return -1;

}


