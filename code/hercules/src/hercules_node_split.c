//
//  hercules_node_split.c
//  
//
//  Created by Karima Echihabi on 18/12/2018
//


#include "../config.h"
#include "../globals.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/hercules_node.h"
#include "../include/hercules_index.h"
#include "../include/hercules_node_split.h"
#include "../include/calc_utils.h"
#include "float.h"


enum response mean_node_segment_split_policy_split(struct node_segment_split_policy * policy, struct segment_sketch sketch, struct segment_sketch * ret)
{

  const int num_splits=2; //default split into 2 node
  
  ts_type max_mean = sketch.indicators[0]; 
  ts_type min_mean = sketch.indicators[1];
  
  policy->indicator_split_idx = 0 ;  //stdev based split
  policy->indicator_split_value = (max_mean + min_mean) / 2;  //the mean value is split value

  ret[0].num_indicators = sketch.num_indicators;
  ret[1].num_indicators = sketch.num_indicators;
    
  for (int i = 0; i < num_splits; i++) {
    for(int j= 0; j<ret[i].num_indicators; ++j){
      ret[i].indicators[j] = sketch.indicators[j];
    }
  }

  ret[0].indicators[1] = policy->indicator_split_value;
  ret[1].indicators[0] = policy->indicator_split_value;

  return SUCCESS;
}

enum response stdev_node_segment_split_policy_split(struct node_segment_split_policy * policy, struct segment_sketch sketch, struct segment_sketch *ret){

  const int num_splits=2; //default split into 2 node
  
  ts_type max_stdev = sketch.indicators[2]; 
  ts_type min_stdev = sketch.indicators[3];

  policy->indicator_split_idx = 1 ;  //stdev based split
  policy->indicator_split_value = (ts_type)(max_stdev + min_stdev) / 2;  //the mean of stdev value is split value

  ret[0].num_indicators = sketch.num_indicators;
  ret[1].num_indicators = sketch.num_indicators;
    
  for (int i = 0; i < num_splits; i++) {
    for(int j= 0; j<ret[i].num_indicators; ++j){
      ret[i].indicators[j] = sketch.indicators[j];
    }
  }

  ret[0].indicators[2] = policy->indicator_split_value;
  ret[1].indicators[3] = policy->indicator_split_value;

  return SUCCESS;
}

enum response series_segment_sketch_do_sketch(struct segment_sketch * series_segment_sketch, ts_type * series, int fromIdx, int toIdx)
{
  
#ifdef __SSE__
  calc_mean_stdev_SIMD(series, fromIdx, toIdx, &series_segment_sketch->indicators[0], &series_segment_sketch->indicators[1]);
#elif
  calc_mean_stdev(series, fromIdx, toIdx, &series_segment_sketch->indicators[0], &series_segment_sketch->indicators[1]);
#endif
  

  return SUCCESS;
}

enum response node_segment_sketch_update_sketch_parallel(struct segment_sketch * node_segment_sketch, ts_type * series, int fromIdx, int toIdx,struct segment_sketch * series_segment_sketch)
{
  if(!series_segment_sketch_do_sketch (series_segment_sketch, series, fromIdx, toIdx))
    {
      fprintf(stderr,"Error in node_init_segments(): Could not calculate the series segment  \
                    sketch.\n");
      return FAILURE;	
    }
    

  if (node_segment_sketch->indicators == NULL)
    {
      node_segment_sketch->indicators =  malloc(sizeof(ts_type) * 4);
      if(node_segment_sketch->indicators == NULL) {
        fprintf(stderr,"Error in hercules_node_split.c: Could not allocate memory for node \
                        segment sketch indicators.\n");
        return FAILURE;
      }     
      node_segment_sketch->indicators[0] = -FLT_MAX ; //for max mean
      node_segment_sketch->indicators[1] = FLT_MAX; //for min mean
      node_segment_sketch->indicators[2] = -FLT_MAX; //for max stdev
      node_segment_sketch->indicators[3] = FLT_MAX; //for min stdev
      node_segment_sketch->num_indicators = 4;                     
    }
   
   
  node_segment_sketch->indicators[0] = fmaxf(node_segment_sketch->indicators[0], series_segment_sketch->indicators[0]);
  node_segment_sketch->indicators[1] = fminf(node_segment_sketch->indicators[1], series_segment_sketch->indicators[0]);
  node_segment_sketch->indicators[2] = fmaxf(node_segment_sketch->indicators[2], series_segment_sketch->indicators[1]);
  node_segment_sketch->indicators[3] = fminf(node_segment_sketch->indicators[3], series_segment_sketch->indicators[1]);
  node_segment_sketch->num_indicators = 4;                     


  return SUCCESS;
}

enum response node_segment_sketch_update_sketch(struct segment_sketch * node_segment_sketch, ts_type * series, int fromIdx, int toIdx)
{

  struct segment_sketch series_segment_sketch;

  series_segment_sketch.indicators = NULL;
  series_segment_sketch.indicators = malloc (sizeof(ts_type)*2);
  if(series_segment_sketch.indicators == NULL) {
    fprintf(stderr,"Error in hercules_node_split.c: Could not allocate memory for series \
                        segment sketch indicators.\n");
    return FAILURE;
  }
  series_segment_sketch.num_indicators = 2;  
   
  if(!series_segment_sketch_do_sketch (&series_segment_sketch, series, fromIdx, toIdx))
    {
      fprintf(stderr,"Error in node_init_segments(): Could not calculate the series segment  \
                    sketch.\n");
      return FAILURE;	
    }
    

  if (node_segment_sketch->indicators == NULL)
    {
      node_segment_sketch->indicators =  malloc(sizeof(ts_type) * 4);
      if(node_segment_sketch->indicators == NULL) {
        fprintf(stderr,"Error in hercules_node_split.c: Could not allocate memory for node \
                        segment sketch indicators.\n");
        return FAILURE;
      }     
      node_segment_sketch->indicators[0] = -FLT_MAX ; //for max mean
      node_segment_sketch->indicators[1] = FLT_MAX; //for min mean
      node_segment_sketch->indicators[2] = -FLT_MAX; //for max stdev
      node_segment_sketch->indicators[3] = FLT_MAX; //for min stdev
      node_segment_sketch->num_indicators = 4;                     
    }
   
   
  node_segment_sketch->indicators[0] = fmaxf(node_segment_sketch->indicators[0], series_segment_sketch.indicators[0]);
  node_segment_sketch->indicators[1] = fminf(node_segment_sketch->indicators[1], series_segment_sketch.indicators[0]);
  node_segment_sketch->indicators[2] = fmaxf(node_segment_sketch->indicators[2], series_segment_sketch.indicators[1]);
  node_segment_sketch->indicators[3] = fminf(node_segment_sketch->indicators[3], series_segment_sketch.indicators[1]);
  node_segment_sketch->num_indicators = 4;                     

  free(series_segment_sketch.indicators);
  series_segment_sketch.indicators = NULL;

  return SUCCESS;
}

boolean node_split_policy_route_to_left (struct hercules_node * node, ts_type * series, struct segment_sketch * series_segment_sketch)
{
  
  struct node_split_policy *policy = node->split_policy;
  boolean route_to_left = false;
  
  if(!series_segment_sketch_do_sketch (series_segment_sketch, series, policy->split_from, policy->split_to))
    {
      fprintf(stderr,"Error in hercules_node_split.c: Could not calculate the series  \
                        segment sketch .\n");
    }

     
  
  if (series_segment_sketch->indicators[policy->indicator_split_idx] < policy->indicator_split_value)
    {
      route_to_left = true;
    }

  return route_to_left;

}



enum response split_node_create_children (struct hercules_index * index, struct hercules_node * node, short * child_node_points, int num_child_node_points)
{
  if (!node->is_leaf)
    {
	
      fprintf(stderr,"Error in hercules_node_split.c: Trying to split a node  \
                        that is not a leaf %s.\n", node->filename);
      return FAILURE;
    }
  else
    {
      
      node->left_child = create_child_node(node);
      if (node->left_child == NULL)
	{
	  fprintf(stderr,"Error in hercules_node_split.c: Left child not \
                        initialized properly.\n");
	  return FAILURE;
	}
    
      if (!node_init_segments(node->left_child,child_node_points,num_child_node_points))
	{
	  fprintf(stderr,"Error in hercules_node_split.c: Could not initialize       \
                        segments for left child of node %s.\n", node->filename);
	  return FAILURE;
	}
    
      node->left_child->is_leaf = 1;
      node->left_child->is_left = true;

      if (!create_hercules_node_filename(index->settings, node->left_child, node))
	{
	  fprintf(stderr,"Error in hercules_node_split.c: Could not create        \
                        filename for left child of node %s.\n", node->filename);
	  return FAILURE;
	}
  
      node->right_child = create_child_node(node);
      if (node->right_child == NULL)
	{
	  fprintf(stderr,"Error in hercules_node_split.c: Right child not \
                        initialized properly.\n");
	  return FAILURE;
	}
    
	
      if (!node_init_segments(node->right_child,child_node_points,num_child_node_points))
	{
	  fprintf(stderr,"Error in hercules_node_split.c: Could not initialize       \
                        segments for right child of node %s.\n", node->filename);
	  return FAILURE;
	}
    
      node->right_child->is_leaf = 1;

      if (!create_hercules_node_filename(index->settings, node->right_child, node))
	{
	  fprintf(stderr,"Error in hercules_node_split.c: Could not create        \
                        filename for right child of node %s.\n", node->filename);
	  return FAILURE;
	}

    }

  return SUCCESS;
}
	
	
struct hercules_node *  create_child_node(struct hercules_node * parent) 
{

    
  struct hercules_node * child_node = hercules_leaf_node_init();

  if (child_node == NULL)
    {
      fprintf(stderr,"Error in hercules_node_split.c: Could not initialize        \
                        child node of parent %s.\n", parent->filename);
      return NULL;
    }    


  child_node->node_segment_split_policies = parent->node_segment_split_policies;
  child_node->num_node_segment_split_policies = parent->num_node_segment_split_policies;
  
  child_node->range = parent->range;

  child_node->parent = parent;

  child_node->level = parent->level + 1;

  child_node->is_leaf = 1;



  return child_node;
}


short get_hs_split_point(short * points, short from, short to, int size_points) 
{

  short * res;
  res= (short *) bsearch(&to, points, size_points,
			 sizeof(short), compare_short);

  

  if(res != NULL) //key found,
    {
      return from;
    }
  else
    {
      return to;
    }
}



boolean is_split_policy_mean(struct node_segment_split_policy policy)
{
  
  if (policy.indicator_split_idx == 0)
    return true;
  else
    return false;
}

boolean is_split_policy_stdev(struct node_segment_split_policy policy)
{
  
  if (policy.indicator_split_idx == 1)
    return true;
  else
    return false;
  

}

ts_type range_calc(struct segment_sketch sketch, int len)
{
  ts_type mean_width = sketch.indicators[0] -sketch.indicators[1];
  ts_type stdev_upper = sketch.indicators[2];
  ts_type stdev_lower = sketch.indicators[3];

  return len * (mean_width * mean_width + stdev_upper * stdev_upper);
  
}


enum response calc_split_points (short * points, int ts_length, int segment_size )
{

  int avg_length = ts_length / segment_size;

  for (int i = 0; i < segment_size; ++i)
    {
      points[i] = (short) ((i+1) * avg_length);    
    }

  //set the last one
  points [segment_size-1] = (short) ts_length;

  return SUCCESS;
}

/*
  This function uses the vertical split points to
  initialize the the horizontal split points in hs_split_points
  The number of hs_split_points is double that of the vertical ones
*/


enum response calc_hs_split_points(short * hs_split_points,int * num_hs_split_points, short * split_points, int segment_size, int min_length)
{
  int c = 0;

  for (int i = 0; i < segment_size; i++) {
    short length = split_points[i]; //i==0
    if (i > 0)
      {
	length = (short) (split_points[i] - split_points[i - 1]);
      }
    if (length >= min_length * 2)
      {
	int start = 0;
	if (i > 0)
	  {
	    start = split_points[i - 1];
	  }
	hs_split_points[c++] = (short) (start + (length / 2));
      }
    hs_split_points[c++] = split_points[i];
  }

  *num_hs_split_points = c;

  return SUCCESS;
}


void  split_node_distribute_data(struct hercules_index * index, struct hercules_node *node, struct segment_sketch * sketch)
{
  ts_type ** ts_list = NULL;
  ts_list =     get_all_time_series_in_node(index, node,0);
      
  for (int idx=0; idx < node->node_size ;++idx)      
    {
      if (node_split_policy_route_to_left(node, ts_list[idx],sketch))
        {
	  if(!update_node_statistics(node->left_child, ts_list[idx]))
	    {
	      fprintf(stderr,"Error in hercules_index.c: could not update \
                             statistics at left child of\
                             node %s\n", node->filename);
	      return FAILURE;
	    }

	  if(!append_ts_to_child_node_parallel(index, node->left_child,ts_list[idx]))
	    {
	      fprintf(stderr,"Error in hercules_index.c: could not append \
                           time series to left child of \
                           node %s\n", node->filename);
	      return FAILURE;
	    }

        }
      else
        {
	  if(!update_node_statistics(node->right_child, ts_list[idx]))
	    {
	      fprintf(stderr,"Error in hercules_index.c: could not update \
                             statistics at right child of\
                             node %s\n", node->filename);
	      return FAILURE;
	    }

	  if(!append_ts_to_child_node_parallel(index, node->right_child,ts_list[idx]))
	    {
	      fprintf(stderr,"Error in hercules_index.c: could not append \
                           time series to right child of \
                           node %s\n", node->filename);
	      return FAILURE;
	    }	  
        }
    }

     
  for (int i = 0 ; i < index->settings->max_leaf_size; ++i)
    {
      free(ts_list[i]);
    }

  free(ts_list);
      
}

void  split_node_distribute_data_parallel(struct hercules_index * index, struct hercules_node *node, struct segment_sketch * sketch, void * tdata)
{
  get_all_time_series_in_node_parallel(index, node,0, tdata);
  ts_type ** ts_list = ((index_buffer_data*)tdata)->split_node_data;	
      
  for (int idx=0; idx < node->node_size ;++idx)      
    {
      if (node_split_policy_route_to_left(node, ts_list[idx],sketch))
        {
	  if(!update_node_statistics(node->left_child, ts_list[idx]))
	    {
	      fprintf(stderr,"Error in hercules_index.c: could not update \
                             statistics at left child of\
                             node %s\n", node->filename);
	      return FAILURE;
	    }

	  if(!append_ts_to_child_node_parallel(index, node->left_child,ts_list[idx]))
	    {
	      fprintf(stderr,"Error in hercules_index.c: could not append \
                           time series to left child of \
                           node %s\n", node->filename);
	      return FAILURE;
	    }

        }
      else
        {
	  if(!update_node_statistics(node->right_child, ts_list[idx]))
	    {
	      fprintf(stderr,"Error in hercules_index.c: could not update \
                             statistics at right child of\
                             node %s\n", node->filename);
	      return FAILURE;
	    }

	  if(!append_ts_to_child_node_parallel(index, node->right_child,ts_list[idx]))
	    {
	      fprintf(stderr,"Error in hercules_index.c: could not append \
                           time series to right child of \
                           node %s\n", node->filename);
	      return FAILURE;
	    }	  
        }
    }
}

void  split_node_detailed (struct hercules_index * index, struct hercules_node *node, struct segment_sketch * sketch, void * thread_data)
{
  node->is_splitting = true;
  struct node_split_policy curr_node_split_policy;
  ts_type max_diff_value = (FLT_MAX * (-1));
  ts_type avg_children_range_value = 0;
  short hs_split_point = -1;
  short * child_node_points;
  int num_child_node_points = 0;
  const int num_child_segments = 2; //by default split to two subsegments
  short split_segment = -1;
      
  double tS = 0;
  double tE = 0;
  struct timeval start_time;
  struct timeval end_time;
  for (int i=0; i <  node->num_node_points;++i)
    { 
      struct segment_sketch curr_node_segment_sketch = node->node_segment_sketches[i];

      //This is the QoS of this segment. QoS is the estimation quality evaluated as =
      //QoS = segment_length * (max_mean_min_mean) * ((max_mean_min_mean) +
      //     (max_stdev * max_stdev))
      //The smaller the QoS, the more effective the bounds are for similarity
      //estimation
	
      ts_type node_range_value = range_calc(curr_node_segment_sketch,
					    get_segment_length(node->node_points, i));

      //for every split policy
      for (int j = 0; j < node->num_node_segment_split_policies;++j)
	{
	  struct node_segment_split_policy curr_node_segment_split_policy =
	    node->node_segment_split_policies[j]; 	  


	  //to hold the two child segments
	  struct segment_sketch * child_node_segment_sketches = NULL;
	  child_node_segment_sketches = init_segment_sketches(num_child_segments,curr_node_segment_sketch.num_indicators);	  

	  if (is_split_policy_mean(curr_node_segment_split_policy))
            mean_node_segment_split_policy_split(&curr_node_segment_split_policy,
						 curr_node_segment_sketch,
						 child_node_segment_sketches);	    
	  else if (is_split_policy_stdev(curr_node_segment_split_policy))
            stdev_node_segment_split_policy_split(&curr_node_segment_split_policy,
						  curr_node_segment_sketch,
						  child_node_segment_sketches);	    
	  else{
            fprintf(stderr,"Error in hercules_index.c: Split policy was not \
                            set properly for node %s\n", node->filename);
	    return FAILURE; 
          }
	  
	  ts_type range_values[num_child_segments];
	  for (int k = 0; k < num_child_segments; ++k)
	    {
	      struct segment_sketch child_node_segment_sketch  = child_node_segment_sketches[k];
	      range_values[k] = range_calc(child_node_segment_sketch,
					   get_segment_length(node->node_points, i));
	    }

	  //diff_value represents the splitting benefit
	  //B = QoS(N) - (QoS_leftNode + QoS_rightNode)/2
	  //the higher the diff_value, the better is the splitting
	  
	  avg_children_range_value = calc_mean(range_values, 0, num_child_segments);
	  ts_type diff_value = node_range_value - avg_children_range_value;
	  
          if (diff_value > max_diff_value) {
            max_diff_value = diff_value;
            curr_node_split_policy.split_from = get_segment_start(node->node_points, i);
            curr_node_split_policy.split_to = get_segment_end(node->node_points, i);
            curr_node_split_policy.indicator_split_idx = curr_node_segment_split_policy.indicator_split_idx;	    	    
            curr_node_split_policy.indicator_split_value = curr_node_segment_split_policy.indicator_split_value;
            curr_node_split_policy.curr_node_segment_split_policy = curr_node_segment_split_policy;
          }
	  for (int k = 0; k< num_child_segments;++k)
	    {
	      free(child_node_segment_sketches[k].indicators);
	      free(child_node_segment_sketches[k].lock_synopsis);	    
	    }
	  free(child_node_segment_sketches);
	} 
    }

  //add trade-off for horizontal split
  max_diff_value = max_diff_value *2;

  //we want to test every possible split policy for each horizontal segment
  for (int i=0; i <  node->num_hs_node_points;++i)
    { 
      struct segment_sketch curr_hs_node_segment_sketch = node->hs_node_segment_sketches[i];
      ts_type node_range_value = range_calc(curr_hs_node_segment_sketch,
					    get_segment_length(node->hs_node_points, i));

      //for every split policy
      for (int j = 0; j < node->num_node_segment_split_policies;++j)
	{
	  struct node_segment_split_policy curr_hs_node_segment_split_policy = node->node_segment_split_policies[j]; 	  

	  struct segment_sketch * child_node_segment_sketches = NULL; //to hold the two child segments
	  child_node_segment_sketches = init_segment_sketches(num_child_segments, curr_hs_node_segment_sketch.num_indicators);
	  
	  if (is_split_policy_mean(curr_hs_node_segment_split_policy))
            mean_node_segment_split_policy_split(&curr_hs_node_segment_split_policy,
						 curr_hs_node_segment_sketch,
						 child_node_segment_sketches);	    
	  else if (is_split_policy_stdev(curr_hs_node_segment_split_policy))
            stdev_node_segment_split_policy_split(&curr_hs_node_segment_split_policy,
						  curr_hs_node_segment_sketch,
						  child_node_segment_sketches);	    
	  else
	    printf("split policy not initialized properly\n");

	  ts_type range_values[num_child_segments];
	  for (int k = 0; k < num_child_segments; ++k)
	    {
	      struct segment_sketch child_node_segment_sketch  = child_node_segment_sketches[k];
	      range_values[k] = range_calc(child_node_segment_sketch,
					   get_segment_length(node->hs_node_points, i));
	    }
	  
	  avg_children_range_value = calc_mean(range_values, 0, num_child_segments);

	  ts_type diff_value = node_range_value - avg_children_range_value;
	    
          if (diff_value > max_diff_value) {
            max_diff_value = diff_value;
            curr_node_split_policy.split_from = get_segment_start(node->hs_node_points, i);
            curr_node_split_policy.split_to = get_segment_end(node->hs_node_points, i);
            curr_node_split_policy.indicator_split_idx = curr_hs_node_segment_split_policy.indicator_split_idx;	    	    
            curr_node_split_policy.indicator_split_value = curr_hs_node_segment_split_policy.indicator_split_value;
            curr_node_split_policy.curr_node_segment_split_policy = curr_hs_node_segment_split_policy;
	    hs_split_point = get_hs_split_point(node->node_points,
						curr_node_split_policy.split_from,
						curr_node_split_policy.split_to,
						node->num_node_points);
          }
	  
	  for (int k = 0; k< num_child_segments;++k)
	    {
	      free(child_node_segment_sketches[k].indicators);
	      free(child_node_segment_sketches[k].lock_synopsis);	    
	    }
	  
          free(child_node_segment_sketches);
	} 
    }
      
  node->split_policy = NULL;
  node->split_policy = malloc (sizeof(struct node_split_policy));
  if(node->split_policy == NULL)
    {
      fprintf(stderr,"Error in hercules_index.c: could not allocate memory \
                        for the split policy of node  %s\n", node->filename);
      return FAILURE;
    }
  node->split_policy->split_from = curr_node_split_policy.split_from;      	
  node->split_policy->split_to = curr_node_split_policy.split_to;
  node->split_policy->indicator_split_idx = curr_node_split_policy.indicator_split_idx;
  node->split_policy->indicator_split_value = curr_node_split_policy.indicator_split_value;
  node->split_policy->curr_node_segment_split_policy = curr_node_split_policy.curr_node_segment_split_policy;


  //when hs_split_point stays less than 0, it means that
  //considering splitting a vertical segment is not worth it
  //according to the QoS heuristic

  if (hs_split_point < 0) 
    {
      num_child_node_points = node->num_node_points;
      child_node_points = NULL;
      child_node_points = malloc(sizeof(short) * num_child_node_points);
      if(child_node_points == NULL)
        {
	  fprintf(stderr,"Error in hercules_index.c: could not allocate memory \
                          for the child node segment points node  %s\n", node->filename);
	  return FAILURE;
        }	
      //children will have the same number of segments as parent
      for(int i = 0; i < num_child_node_points;++i)
	{
	  child_node_points[i] = node->node_points[i];  
	}
    }
  else
    {
      num_child_node_points = node->num_node_points + 1;	
      child_node_points = NULL;
      child_node_points = malloc(sizeof(short) * num_child_node_points);
      if(child_node_points == NULL)
        {
	  fprintf(stderr,"Error in hercules_index.c: could not allocate memory \
                          for the child node segment points node  %s\n", node->filename);
	  return FAILURE;
        }
      //children will have one additional segment than the parent
      for(int i = 0; i < (num_child_node_points-1);++i)
	{
	  child_node_points[i] = node->node_points[i];  
	}
      child_node_points[num_child_node_points-1] = hs_split_point; //initialize newly added point

      qsort(child_node_points,num_child_node_points,sizeof(short),compare_short);

    }
      
  //this will put the time series of this node in the file_buffer->buffered_list aray
  //it will include the time series in disk and those in memory
      
  if(!split_node_create_children(index,node, child_node_points, num_child_node_points))
    {
      fprintf(stderr,"Error in hercules_index.c: could not split node %s.\n", node->filename);
      return FAILURE;
    }

  free(child_node_points);

  node->file_buffer->do_not_flush = true;

      
  ts_type ** ts_list;
  get_all_time_series_in_node_parallel(index, node,0,thread_data);
  ts_list = ((index_buffer_data*)thread_data)->split_node_data;
	
      
  for (int idx=0; idx < node->node_size ;++idx)      
    {
      if (node_split_policy_route_to_left(node, ts_list[idx],sketch))
        {
	  if(!update_node_statistics_parallel(node->left_child, ts_list[idx],sketch ))
	    {
	      fprintf(stderr,"Error in hercules_index.c: could not update \
                             statistics at left child of\
                             node %s\n", node->filename);
	      return FAILURE;
	    }

	  if(!append_ts_to_child_node_parallel(index, node->left_child,ts_list[idx], thread_data))
	    {
	      fprintf(stderr,"Error in hercules_index.c: could not append \
                           time series to left child of \
                           node %s\n", node->filename);
	      return FAILURE;
	    }

        }
      else
        {
	  if(!update_node_statistics_parallel(node->right_child, ts_list[idx], sketch))
	    {
	      fprintf(stderr,"Error in hercules_index.c: could not update \
                             statistics at right child of\
                             node %s\n", node->filename);
	      return FAILURE;
	    }

	  if(!append_ts_to_child_node_parallel(index, node->right_child,ts_list[idx], thread_data))
	    {
	      fprintf(stderr,"Error in hercules_index.c: could not append \
                           time series to right child of \
                           node %s\n", node->filename);
	      return FAILURE;
	    }	  
        }
    }

      
  if(!delete_file_buffer_parallel(index,node))
    {
      fprintf(stderr,"Error in hercules_index.c: could not delete file buffer for \
                           node %s\n", node->filename);
      return FAILURE;
    }
  node->is_splitting = false;
  __sync_fetch_and_add(&(index->stats->leaf_nodes_count),1);

}


void  split_node (struct hercules_index * index, struct hercules_node *node, struct segment_sketch * sketch, void * thread_data)
{

  //printf("thread = %lu splitted_node = %p level = %d BEGIN\n", pthread_self(),(void *)node, node->level);      
  node->is_splitting = true;

  short * child_node_points;
  int num_child_node_points = 0;
  short hs_split_point = -1;
  
  short split_segment = -1;

  double tS = 0;
  double tE = 0;
  struct timeval start_time;
  struct timeval end_time;

  //we want to test every possible split policy for each segment
      
  //for each segment
  split_node_evaluate_policies(index, node,&hs_split_point);
  if (hs_split_point < 0) 
    {
      num_child_node_points = node->num_node_points;
      child_node_points = NULL;
      child_node_points = calloc(num_child_node_points,sizeof(short));
      if(child_node_points == NULL)
        {
	  fprintf(stderr,"Error in hercules_index.c: could not allocate memory \
                          for the child node segment points node  %s\n", node->filename);
	  return FAILURE;
        }	
      //children will have the same number of segments as parent
      for(int i = 0; i < num_child_node_points;++i)
	{
	  child_node_points[i] = node->node_points[i];  
	}
    }
  else
    {
      num_child_node_points = node->num_node_points + 1;	
      child_node_points = NULL;
      child_node_points = calloc(num_child_node_points,sizeof(short));
      if(child_node_points == NULL)
        {
	  fprintf(stderr,"Error in hercules_index.c: could not allocate memory \
                          for the child node segment points node  %s\n", node->filename);
	  return FAILURE;
        }
      //children will have one additional segment than the parent
      for(int i = 0; i < num_child_node_points-1;++i)
	{
	  child_node_points[i] = node->node_points[i];  
	}
      child_node_points[num_child_node_points-1] = hs_split_point; //initialize newly added point

      qsort(child_node_points,num_child_node_points,sizeof(short),compare_short);

    }

  //this will put the time series of this node in the file_buffer->buffered_list aray
  //it will include the time series in disk and those in memory

#if DETAILED_STATS == 1
  gettimeofday(&end_time, NULL);
  tS = start_time.tv_sec*1000000 + (start_time.tv_usec);
  tE = end_time.tv_sec*1000000 + (end_time.tv_usec);
  ((index_buffer_data*)thread_data)->thread_evaluate_total_time += (tE - tS);
  gettimeofday(&start_time, NULL); 
#endif
  if(!split_node_create_children(index,node, child_node_points, num_child_node_points))
    {
      fprintf(stderr,"Error in hercules_index.c: could not split node %s.\n", node->filename);
      return FAILURE;
    }

  free(child_node_points);

  node->file_buffer->do_not_flush = true;

     
  if (!get_file_buffer_parallel(index, node))
    { 
      fprintf(stderr,"Error in hercules_index.c: could not get the file \
                           buffer for node %s.\n", node->filename);
      return FAILURE;              
    }

  //split_node_distribute_data_parallel(index,node, sketch, thread_data);
  split_node_distribute_data(index,node, sketch);
  if(!delete_file_buffer_parallel(index,node))
    {
      fprintf(stderr,"Error in hercules_index.c: could not delete file buffer for \
                           node %s\n", node->filename);
      return FAILURE;
    }
  node->is_splitting = false;

}

void split_node_evaluate_policies (struct hercules_index * index, struct hercules_node *node, int * hs_split_point)
{
  ts_type max_diff_value = (FLT_MAX * (-1));
  ts_type avg_children_range_value = 0;

  const int num_child_segments = 2; //by default split to two subsegments
  struct node_split_policy curr_node_split_policy;
  
  for (int i=0; i <  node->num_node_points;++i)
    { 
      struct segment_sketch curr_node_segment_sketch = node->node_segment_sketches[i];

      //This is the QoS of this segment. QoS is the estimation quality evaluated as =
      //QoS = segment_length * (max_mean_min_mean) * ((max_mean_min_mean) +
      //     (max_stdev * max_stdev))
      //The smaller the QoS, the more effective the bounds are for similarity
      //estimation
	
      ts_type node_range_value = range_calc(curr_node_segment_sketch,
					    get_segment_length(node->node_points, i));

      //for every split policy
      for (int j = 0; j < node->num_node_segment_split_policies;++j)
	{
	  struct node_segment_split_policy curr_node_segment_split_policy =
	    node->node_segment_split_policies[j]; 	  


	  //to hold the two child segments
	  struct segment_sketch * child_node_segment_sketches = NULL;
	  child_node_segment_sketches = init_segment_sketches(num_child_segments,curr_node_segment_sketch.num_indicators);	  

	  if (is_split_policy_mean(curr_node_segment_split_policy))
            mean_node_segment_split_policy_split(&curr_node_segment_split_policy,
						 curr_node_segment_sketch,
						 child_node_segment_sketches);	    
	  else if (is_split_policy_stdev(curr_node_segment_split_policy))
            stdev_node_segment_split_policy_split(&curr_node_segment_split_policy,
						  curr_node_segment_sketch,
						  child_node_segment_sketches);	    
	  else{
            fprintf(stderr,"Error in hercules_index.c: Split policy was not \
                            set properly for node %s\n", node->filename);
	    return FAILURE; 
          }
	  
	  ts_type range_values[num_child_segments];
	  for (int k = 0; k < num_child_segments; ++k)
	    {
	      struct segment_sketch child_node_segment_sketch  = child_node_segment_sketches[k];
	      range_values[k] = range_calc(child_node_segment_sketch,
					   get_segment_length(node->node_points, i));
	    }

	  //diff_value represents the splitting benefit
	  //B = QoS(N) - (QoS_leftNode + QoS_rightNode)/2
	  //the higher the diff_value, the better is the splitting
	  
	  avg_children_range_value = calc_mean(range_values, 0, num_child_segments);
	  ts_type diff_value = node_range_value - avg_children_range_value;
	  
          if (diff_value > max_diff_value) {
            max_diff_value = diff_value;
            curr_node_split_policy.split_from = get_segment_start(node->node_points, i);
            curr_node_split_policy.split_to = get_segment_end(node->node_points, i);
            curr_node_split_policy.indicator_split_idx = curr_node_segment_split_policy.indicator_split_idx;	    	    
            curr_node_split_policy.indicator_split_value = curr_node_segment_split_policy.indicator_split_value;
            curr_node_split_policy.curr_node_segment_split_policy = curr_node_segment_split_policy;
          }
	  for (int k = 0; k< num_child_segments;++k)
	    {
	      free(child_node_segment_sketches[k].indicators);
	      free(child_node_segment_sketches[k].lock_synopsis);	    
	    }
	  free(child_node_segment_sketches);
	} 
    }

  //add trade-off for horizontal split
  max_diff_value = max_diff_value *2;

  //we want to test every possible split policy for each horizontal segment
  for (int i=0; i <  node->num_hs_node_points;++i)
    { 
      struct segment_sketch curr_hs_node_segment_sketch = node->hs_node_segment_sketches[i];
      ts_type node_range_value = range_calc(curr_hs_node_segment_sketch,
					    get_segment_length(node->hs_node_points, i));

      //for every split policy
      for (int j = 0; j < node->num_node_segment_split_policies;++j)
	{
	  struct node_segment_split_policy curr_hs_node_segment_split_policy = node->node_segment_split_policies[j]; 	  

	  struct segment_sketch * child_node_segment_sketches = NULL; //to hold the two child segments
	  child_node_segment_sketches = init_segment_sketches(num_child_segments, curr_hs_node_segment_sketch.num_indicators);
	  
	  if (is_split_policy_mean(curr_hs_node_segment_split_policy))
            mean_node_segment_split_policy_split(&curr_hs_node_segment_split_policy,
						 curr_hs_node_segment_sketch,
						 child_node_segment_sketches);	    
	  else if (is_split_policy_stdev(curr_hs_node_segment_split_policy))
            stdev_node_segment_split_policy_split(&curr_hs_node_segment_split_policy,
						  curr_hs_node_segment_sketch,
						  child_node_segment_sketches);	    
	  else
	    printf("split policy not initialized properly\n");

	  ts_type range_values[num_child_segments];
	  for (int k = 0; k < num_child_segments; ++k)
	    {
	      struct segment_sketch child_node_segment_sketch  = child_node_segment_sketches[k];
	      range_values[k] = range_calc(child_node_segment_sketch,
					   get_segment_length(node->hs_node_points, i));
	    }
	  
	  avg_children_range_value = calc_mean(range_values, 0, num_child_segments);

	  ts_type diff_value = node_range_value - avg_children_range_value;
	    
          if (diff_value > max_diff_value) {
            max_diff_value = diff_value;
            curr_node_split_policy.split_from = get_segment_start(node->hs_node_points, i);
            curr_node_split_policy.split_to = get_segment_end(node->hs_node_points, i);
            curr_node_split_policy.indicator_split_idx = curr_hs_node_segment_split_policy.indicator_split_idx;	    	    
            curr_node_split_policy.indicator_split_value = curr_hs_node_segment_split_policy.indicator_split_value;
            curr_node_split_policy.curr_node_segment_split_policy = curr_hs_node_segment_split_policy;
	    hs_split_point = get_hs_split_point(node->node_points,
						curr_node_split_policy.split_from,
						curr_node_split_policy.split_to,
						node->num_node_points);
          }
	  
	  for (int k = 0; k< num_child_segments;++k)
	    {
	      free(child_node_segment_sketches[k].indicators);
	      free(child_node_segment_sketches[k].lock_synopsis);	    
	    }
	  
          free(child_node_segment_sketches);
	} 
    }
      
  node->split_policy = NULL;
  node->split_policy = malloc (sizeof(struct node_split_policy));
  if(node->split_policy == NULL)
    {
      fprintf(stderr,"Error in hercules_index.c: could not allocate memory \
                        for the split policy of node  %s\n", node->filename);
      return FAILURE;
    }
  node->split_policy->split_from = curr_node_split_policy.split_from;      	
  node->split_policy->split_to = curr_node_split_policy.split_to;
  node->split_policy->indicator_split_idx = curr_node_split_policy.indicator_split_idx;
  node->split_policy->indicator_split_value = curr_node_split_policy.indicator_split_value;
  node->split_policy->curr_node_segment_split_policy = curr_node_split_policy.curr_node_segment_split_policy;



}

void  hercules_index_split_node (struct hercules_index * index, struct hercules_node *node, struct segment_sketch * sketch, void * thread_data)
{
  //printf("thread = %lu splitted_node = %p level = %d BEGIN\n", pthread_self(),(void *)node, node->level);      
  node->is_splitting = true;
  struct node_split_policy curr_node_split_policy;
  ts_type max_diff_value = (FLT_MAX * (-1));
  ts_type avg_children_range_value = 0;
  short hs_split_point = -1;
  short * child_node_points;
  int num_child_node_points = 0;
  const int num_child_segments = 2; //by default split to two subsegments
  short split_segment = -1;
      
  //we want to test every possible split policy for each segment
  double tS = 0;
  double tE = 0;
  struct timeval start_time;
  struct timeval end_time;
  //for each segment 
  for (int i=0; i <  node->num_node_points;++i)
    { 
      struct segment_sketch curr_node_segment_sketch = node->node_segment_sketches[i];

      //This is the QoS of this segment. QoS is the estimation quality evaluated as =
      //QoS = segment_length * (max_mean_min_mean) * ((max_mean_min_mean) +
      //     (max_stdev * max_stdev))
      //The smaller the QoS, the more effective the bounds are for similarity
      //estimation
	
      ts_type node_range_value = range_calc(curr_node_segment_sketch,
					    get_segment_length(node->node_points, i));

      //for every split policy
      for (int j = 0; j < node->num_node_segment_split_policies;++j)
	{
	  struct node_segment_split_policy curr_node_segment_split_policy =
	    node->node_segment_split_policies[j]; 	  


	  //to hold the two child segments
	  struct segment_sketch * child_node_segment_sketches = NULL;
	  child_node_segment_sketches = init_segment_sketches(num_child_segments,curr_node_segment_sketch.num_indicators);	  

	  if (is_split_policy_mean(curr_node_segment_split_policy))
            mean_node_segment_split_policy_split(&curr_node_segment_split_policy,
						 curr_node_segment_sketch,
						 child_node_segment_sketches);	    
	  else if (is_split_policy_stdev(curr_node_segment_split_policy))
            stdev_node_segment_split_policy_split(&curr_node_segment_split_policy,
						  curr_node_segment_sketch,
						  child_node_segment_sketches);	    
	  else{
            fprintf(stderr,"Error in hercules_index.c: Split policy was not \
                            set properly for node %s\n", node->filename);
	    return FAILURE; 
          }
	  
	  ts_type range_values[num_child_segments];
	  for (int k = 0; k < num_child_segments; ++k)
	    {
	      struct segment_sketch child_node_segment_sketch  = child_node_segment_sketches[k];
	      range_values[k] = range_calc(child_node_segment_sketch,
					   get_segment_length(node->node_points, i));
	    }

	  //diff_value represents the splitting benefit
	  //B = QoS(N) - (QoS_leftNode + QoS_rightNode)/2
	  //the higher the diff_value, the better is the splitting
	  
	  avg_children_range_value = calc_mean(range_values, 0, num_child_segments);
	  ts_type diff_value = node_range_value - avg_children_range_value;
	  
          if (diff_value > max_diff_value) {
            max_diff_value = diff_value;
            curr_node_split_policy.split_from = get_segment_start(node->node_points, i);
            curr_node_split_policy.split_to = get_segment_end(node->node_points, i);
            curr_node_split_policy.indicator_split_idx = curr_node_segment_split_policy.indicator_split_idx;	    	    
            curr_node_split_policy.indicator_split_value = curr_node_segment_split_policy.indicator_split_value;
            curr_node_split_policy.curr_node_segment_split_policy = curr_node_segment_split_policy;
          }
	  for (int k = 0; k< num_child_segments;++k)
	    {
	      free(child_node_segment_sketches[k].indicators);
	      free(child_node_segment_sketches[k].lock_synopsis);	    
	    }
	  free(child_node_segment_sketches);
	} 
    }

  //add trade-off for horizontal split
  max_diff_value = max_diff_value *2;

  //we want to test every possible split policy for each horizontal segment
  for (int i=0; i <  node->num_hs_node_points;++i)
    { 
      struct segment_sketch curr_hs_node_segment_sketch = node->hs_node_segment_sketches[i];
      ts_type node_range_value = range_calc(curr_hs_node_segment_sketch,
					    get_segment_length(node->hs_node_points, i));

      //for every split policy
      for (int j = 0; j < node->num_node_segment_split_policies;++j)
	{
	  struct node_segment_split_policy curr_hs_node_segment_split_policy = node->node_segment_split_policies[j]; 	  

	  struct segment_sketch * child_node_segment_sketches = NULL; //to hold the two child segments
	  child_node_segment_sketches = init_segment_sketches(num_child_segments, curr_hs_node_segment_sketch.num_indicators);
	  
	  if (is_split_policy_mean(curr_hs_node_segment_split_policy))
            mean_node_segment_split_policy_split(&curr_hs_node_segment_split_policy,
						 curr_hs_node_segment_sketch,
						 child_node_segment_sketches);	    
	  else if (is_split_policy_stdev(curr_hs_node_segment_split_policy))
            stdev_node_segment_split_policy_split(&curr_hs_node_segment_split_policy,
						  curr_hs_node_segment_sketch,
						  child_node_segment_sketches);	    
	  else
	    printf("split policy not initialized properly\n");

	  ts_type range_values[num_child_segments];
	  for (int k = 0; k < num_child_segments; ++k)
	    {
	      struct segment_sketch child_node_segment_sketch  = child_node_segment_sketches[k];
	      range_values[k] = range_calc(child_node_segment_sketch,
					   get_segment_length(node->hs_node_points, i));
	    }
	  
	  avg_children_range_value = calc_mean(range_values, 0, num_child_segments);

	  ts_type diff_value = node_range_value - avg_children_range_value;
	    
          if (diff_value > max_diff_value) {
            max_diff_value = diff_value;
            curr_node_split_policy.split_from = get_segment_start(node->hs_node_points, i);
            curr_node_split_policy.split_to = get_segment_end(node->hs_node_points, i);
            curr_node_split_policy.indicator_split_idx = curr_hs_node_segment_split_policy.indicator_split_idx;	    	    
            curr_node_split_policy.indicator_split_value = curr_hs_node_segment_split_policy.indicator_split_value;
            curr_node_split_policy.curr_node_segment_split_policy = curr_hs_node_segment_split_policy;
	    hs_split_point = get_hs_split_point(node->node_points,
						curr_node_split_policy.split_from,
						curr_node_split_policy.split_to,
						node->num_node_points);
          }
	  
	  for (int k = 0; k< num_child_segments;++k)
	    {
	      free(child_node_segment_sketches[k].indicators);
	      free(child_node_segment_sketches[k].lock_synopsis);	    
	    }
	  
          free(child_node_segment_sketches);
	} 
    }
      
  node->split_policy = NULL;
  node->split_policy = malloc (sizeof(struct node_split_policy));
  if(node->split_policy == NULL)
    {
      fprintf(stderr,"Error in hercules_index.c: could not allocate memory \
                        for the split policy of node  %s\n", node->filename);
      return FAILURE;
    }
  node->split_policy->split_from = curr_node_split_policy.split_from;      	
  node->split_policy->split_to = curr_node_split_policy.split_to;
  node->split_policy->indicator_split_idx = curr_node_split_policy.indicator_split_idx;
  node->split_policy->indicator_split_value = curr_node_split_policy.indicator_split_value;
  node->split_policy->curr_node_segment_split_policy = curr_node_split_policy.curr_node_segment_split_policy;


  //when hs_split_point stays less than 0, it means that
  //considering splitting a vertical segment is not worth it
  //according to the QoS heuristic

  if (hs_split_point < 0) 
    {
      num_child_node_points = node->num_node_points;
      child_node_points = NULL;
      child_node_points = malloc(sizeof(short) * num_child_node_points);
      if(child_node_points == NULL)
        {
	  fprintf(stderr,"Error in hercules_index.c: could not allocate memory \
                          for the child node segment points node  %s\n", node->filename);
	  return FAILURE;
        }	
      //children will have the same number of segments as parent
      for(int i = 0; i < num_child_node_points;++i)
	{
	  child_node_points[i] = node->node_points[i];  
	}
    }
  else
    {
      num_child_node_points = node->num_node_points + 1;	
      child_node_points = NULL;
      child_node_points = malloc(sizeof(short) * num_child_node_points);
      if(child_node_points == NULL)
        {
	  fprintf(stderr,"Error in hercules_index.c: could not allocate memory \
                          for the child node segment points node  %s\n", node->filename);
	  return FAILURE;
        }
      //children will have one additional segment than the parent
      for(int i = 0; i < (num_child_node_points-1);++i)
	{
	  child_node_points[i] = node->node_points[i];  
	}
      child_node_points[num_child_node_points-1] = hs_split_point; //initialize newly added point

      qsort(child_node_points,num_child_node_points,sizeof(short),compare_short);

    }
      
  //this will put the time series of this node in the file_buffer->buffered_list aray
  //it will include the time series in disk and those in memory
  /*
    #if DETAILED_STATS == 1
    gettimeofday(&end_time, NULL);
    tS = start_time.tv_sec*1000000 + (start_time.tv_usec);
    tE = end_time.tv_sec*1000000 + (end_time.tv_usec);
    ((index_buffer_data*)thread_data)->thread_evaluate_total_time += (tE - tS);
    gettimeofday(&start_time, NULL); 
    #endif
  */
      
  if(!split_node_create_children(index,node, child_node_points, num_child_node_points))
    {
      fprintf(stderr,"Error in hercules_index.c: could not split node %s.\n", node->filename);
      return FAILURE;
    }

  free(child_node_points);

  node->file_buffer->do_not_flush = true;

  /*
    if (!get_file_buffer_parallel(index, node))
    { 
    fprintf(stderr,"Error in hercules_index.c: could not get the file \
    buffer for node %s.\n", node->filename);
    return FAILURE;              
    }
  */
     
      
  ts_type ** ts_list;
  //ts_list = get_all_time_series_in_node(index, node,0);
  get_all_time_series_in_node_parallel(index, node,0,thread_data);
  ts_list = ((index_thread_data*)thread_data)->split_node_data;
	
  //copying the contents of the the node being split
  //in case it gets flushed from memory to disk
      
  for (int idx=0; idx < node->node_size ;++idx)      
    {
      if (node_split_policy_route_to_left(node, ts_list[idx],sketch))
        {
	  if(!update_node_statistics_parallel(node->left_child, ts_list[idx],sketch ))
	    {
	      fprintf(stderr,"Error in hercules_index.c: could not update \
                             statistics at left child of\
                             node %s\n", node->filename);
	      return FAILURE;
	    }

	  if(!append_ts_to_child_node_parallel(index, node->left_child,ts_list[idx], thread_data))
	    {
	      fprintf(stderr,"Error in hercules_index.c: could not append \
                           time series to left child of \
                           node %s\n", node->filename);
	      return FAILURE;
	    }

        }
      else
        {
	  if(!update_node_statistics_parallel(node->right_child, ts_list[idx], sketch))
	    {
	      fprintf(stderr,"Error in hercules_index.c: could not update \
                             statistics at right child of\
                             node %s\n", node->filename);
	      return FAILURE;
	    }

	  if(!append_ts_to_child_node_parallel(index, node->right_child,ts_list[idx], thread_data))
	    {
	      fprintf(stderr,"Error in hercules_index.c: could not append \
                           time series to right child of \
                           node %s\n", node->filename);
	      return FAILURE;
	    }	  
        }
    }

      
  if(!delete_file_buffer_parallel(index,node))
    {
      fprintf(stderr,"Error in hercules_index.c: could not delete file buffer for \
                           node %s\n", node->filename);
      return FAILURE;
    }
  node->is_splitting = false;
  __sync_fetch_and_add(&(index->stats->leaf_nodes_count),1);

}

