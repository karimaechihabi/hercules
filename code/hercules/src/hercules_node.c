//
//  hercules_node.c
//
//  Created by Karima Echihabi on 18/12/2018
//



#include <stdio.h>
#include <stdlib.h>
#include "../config.h"
#include "../globals.h"
#include "../include/hercules_node.h"
#include "../include/hercules_index.h"
#include <math.h>
#include <pqueue.h>
#include <float.h>
#include <limits.h>
#include "../include/pqueue.h"
#include "../include/hercules_query_engine.h"
#include "immintrin.h"

/**
 This function initializes a hercules root node.
 */

struct hercules_node * hercules_root_node_init(struct hercules_index_settings * settings) 
{

    struct hercules_node * node = hercules_leaf_node_init();
    int ts_length = settings->timeseries_size;
    int segment_size = settings->init_segments;
    
    node->node_segment_split_policies = malloc (sizeof (struct node_segment_split_policy) * 2 );

    if(node->node_segment_split_policies == NULL) {
     fprintf(stderr,"Error in hercules_node.c: Could not allocate memory for root node segment \
                     split policies.\n");
     return FAILURE;	
    }
    
    node->node_segment_split_policies[0].indicator_split_idx = 0;  //Mean  
    node->node_segment_split_policies[1].indicator_split_idx = 1;  //Stdev
    node->num_node_segment_split_policies = 2;  //Stdev and Mean
          
    //calc the split points by segmentSize
    short * split_points = NULL;
    split_points = malloc (sizeof(short) * segment_size);

    if(split_points == NULL) {
      fprintf(stderr,"Error in hercules_node.c: Could not \
                     allocate memory for root split points.\n");
      return FAILURE;	
    }

    if(!calc_split_points(split_points, ts_length, segment_size)) 
    {
      fprintf(stderr,"Error in hercules_node.c: Could not \
                     calculate the split points for the root.\n");
      return FAILURE;	
    }
    
    if(!node_init_segments(node, split_points, segment_size))
    {
      fprintf(stderr,"Error in hercules_node.c: Could not \
                     initialize the segments for the root.\n");
      return FAILURE;	
    }
    
    if(!create_hercules_node_filename(settings, node,NULL))
    {
      fprintf(stderr,"Error in hercules_node.c: Could not \
                     create a filename for the root node.\n");
      return FAILURE;	
    }

    free(split_points);
    return node;
}

/**
 This function initalizes a hercules leaf node.
 */

struct hercules_node * hercules_leaf_node_init(void) 
{
    COUNT_NEW_NODE
 
    struct hercules_node *node = malloc(sizeof(struct hercules_node));
    if(node == NULL) {
        fprintf(stderr,"error: could not allocate memory for new node.\n");
        return NULL;
    }
    node->right_child = NULL;
    node->left_child = NULL;
    node->parent = NULL;

    node->filename = NULL;
    
    node->node_segment_split_policies=NULL;
    node->num_node_segment_split_policies=2;
    
    node->range=0;  

    node->level = 0;   
    node->is_left = false;
    node->is_leaf = true;
    node->is_pruned = true;    //a leaf is pruned by default
    
    node->split_policy=NULL; 
    node->node_points=NULL;   
    node->hs_node_points=NULL;
    node->num_node_points=0;   
    node->num_hs_node_points=0;
    node->split_segment = -1;
    
    node->node_segment_sketches=NULL;
    node->hs_node_segment_sketches=NULL;
    
    node->max_segment_length = 2; 
    node->max_value_length = 10; 

    node->file_buffer = NULL;
      
    node->node_size = 0;
    
    node->proc_finished = 0;
    node->write_finished = 0;

    node->thread_id = -1;

    pthread_mutex_init(&(node->lock_data), NULL);
    pthread_mutex_init(&(node->lock_split), NULL);    

    hercules_file_buffer_init(node);
    
    return node;
}



enum response node_init_segments(struct hercules_node * node, short * split_points, int segment_size)
{     

  node->num_node_points = segment_size;

  node->node_points = NULL;
  node->node_points = malloc (sizeof(short) * segment_size);
  if(node->node_points == NULL) {
     fprintf(stderr,"Error in node_init_segments(): Could not allocate memory for node points.\n");
     return FAILURE;	
  }
    
  for (int i = 0; i < segment_size; ++i)
  {
    node->node_points[i] = split_points[i];
  }
  
  node->hs_node_points = NULL;
  node->hs_node_points = malloc (sizeof(short) * segment_size * 2);
  if(node->hs_node_points == NULL) {
     fprintf(stderr,"Error in node_init_segments(): Could not allocate memory for hs node points.\n");
     return FAILURE;	
  }
  
  int min_length = 1; //mininum length of new segment = 1
  int num_hs_split_points = 0;
    
  calc_hs_split_points(node->hs_node_points, &node->num_hs_node_points,node->node_points, segment_size, min_length);

  //allocate mem for array of sketches


  node->node_segment_sketches = NULL;
  node->node_segment_sketches = malloc (sizeof(struct segment_sketch) * segment_size);
  if(node->node_segment_sketches == NULL) {
     fprintf(stderr,"Error in node_init_segments(): Could not allocate memory for node segment sketches.\n");
     return FAILURE;	
  }
  
  node->hs_node_segment_sketches = NULL;  
  node->hs_node_segment_sketches = malloc (sizeof(struct segment_sketch) * node->num_hs_node_points);
  if(node->hs_node_segment_sketches == NULL) {
     fprintf(stderr,"Error in node_init_segments(): Could not allocate memory for hs node segment sketches.\n");
     return FAILURE;	
  }
  
  //allocate memory for vertical indicators
  for (int i = 0; i < segment_size; ++i)
  {
    node->node_segment_sketches[i].indicators = NULL;    
    node->node_segment_sketches[i].indicators = malloc (sizeof(ts_type) * 4); //node segment has 4 indicators

    if(node->node_segment_sketches[i].indicators == NULL) {
     fprintf(stderr,"Error in node_init_segments(): Could not allocate memory for node segment \
                     sketch indicator.\n");
     return FAILURE;	
    }
    
    node->node_segment_sketches[i].indicators[0]= -FLT_MAX;
    node->node_segment_sketches[i].indicators[1]= FLT_MAX;
    node->node_segment_sketches[i].indicators[2]= -FLT_MAX;
    node->node_segment_sketches[i].indicators[3]= FLT_MAX;
    node->node_segment_sketches[i].num_indicators= 4;
  }

  //allocate memory for horizontal indicators
  for (int i = 0; i < node->num_hs_node_points; ++i)
  {
    node->hs_node_segment_sketches[i].indicators = NULL;
    node->hs_node_segment_sketches[i].indicators = malloc (sizeof(ts_type) * 4); //node segment has 4 indicators
    if(node->hs_node_segment_sketches[i].indicators == NULL) {
     fprintf(stderr,"Error in node_init_segments(): Could not allocate memory for hs node segment \
                     sketch indicator.\n");
     return FAILURE;	
    }
    node->hs_node_segment_sketches[i].indicators[0]= -FLT_MAX;
    node->hs_node_segment_sketches[i].indicators[1]= FLT_MAX;
    node->hs_node_segment_sketches[i].indicators[2]= -FLT_MAX;
    node->hs_node_segment_sketches[i].indicators[3]= FLT_MAX;    
    node->hs_node_segment_sketches[i].num_indicators= 4;
  }  

  return SUCCESS;

}

enum response append_ts_to_node(struct hercules_index * index,
				struct hercules_node * node,
				ts_type * timeseries)
{

  if (!get_file_buffer(index, node))
  {
    fprintf(stderr, "Error in hercules_index.c:  Could not get the \
                     file buffer for this node.\n");
    return FAILURE;              
  }

  if (node->file_buffer == NULL)
  {
    fprintf(stderr, "Error in hercules_index.c:  Null file buffer for \
                     this node after creating it.\n");
    return FAILURE;              
  }

  int idx = node->file_buffer->buffered_list_size;  

  int ts_length = index->settings->timeseries_size;
  int max_leaf_size = index->settings->max_leaf_size;

  if (idx == 0)
  {
      node->file_buffer->buffered_list = NULL;
      node->file_buffer->buffered_list = malloc(sizeof(struct ts_type *) * max_leaf_size);
      
      if (node->file_buffer->buffered_list == NULL)
      {
        fprintf(stderr, "Error in hercules_index.c:  Could not \
                         allocate memory for the buffered list. \n");
        return FAILURE;                  
      }
  }
  /*
  node->file_buffer->buffered_list[idx] = NULL;
  node->file_buffer->buffered_list[idx] = malloc(sizeof(ts_type) * ts_length);
  */

  node->file_buffer->buffered_list[idx] = (ts_type *) index->buffer_manager->current_record;
  index->buffer_manager->current_record += sizeof(ts_type) * ts_length;
  index->buffer_manager->current_record_index++;
    
  if (node->file_buffer->buffered_list[idx] == NULL)
  {
        fprintf(stderr, "Error in hercules_index.c:  Could not \
                         allocate memory for the time series in\
                         the buffer.\n");
        return FAILURE;                  
  }
  
  for(int i=0; i<ts_length; ++i)
  {
        node->file_buffer->buffered_list[idx][i] = timeseries[i];
  }
    
  ++node->file_buffer->buffered_list_size;
  index->buffer_manager->current_count += ts_length; 
    
  return SUCCESS;
}
/*
enum response append_ts_to_node_parallel(struct hercules_index * index,
					 struct hercules_node * node,
					 ts_type * timeseries)

{
  char *cr;

  
  if (!get_file_buffer_parallel(index, node))
  {
    fprintf(stderr, "Error in hercules_index.c:  Could not get the \
                     file buffer for this node.\n");
    return FAILURE;              
  }

  if (node->file_buffer == NULL)
  {
    fprintf(stderr, "Error in hercules_index.c:  Null file buffer for \
                     this node after creating it.\n");
    return FAILURE;              
  }

  int idx = node->file_buffer->buffered_list_size;  

  int ts_length = index->settings->timeseries_size;
  int max_leaf_size = index->settings->max_leaf_size;

  if (idx == 0)
  {
      node->file_buffer->buffered_list = NULL;
      node->file_buffer->buffered_list = malloc(sizeof(struct ts_type *) * max_leaf_size);
      
      if (node->file_buffer->buffered_list == NULL)
      {
        fprintf(stderr, "Error in hercules_index.c:  Could not \
                         allocate memory for the buffered list. \n");
        return FAILURE;                  
      }
  }
  

  //index->buffer_manager->current_record_index++;
  __sync_fetch_and_add(&(index->buffer_manager->current_record_index),1);
  cr = __sync_fetch_and_add(&(index->buffer_manager->current_record),sizeof(ts_type) * ts_length);
  node->file_buffer->buffered_list[idx] = (ts_type *) cr;
  memcpy((void *) cr, (void *) timeseries, sizeof(ts_type) * ts_length);

  
  if (node->file_buffer->buffered_list[idx] == NULL)
  {
        fprintf(stderr, "Error in hercules_index.c:  Could not \
                         allocate memory for the time series in\
                         the buffer.\n");
        return FAILURE;                  
  }


  //pthread_mutex_lock(&(index->buffer_manager->lock_file_buffer)); 

  //pthread_mutex_unlock(&(index->buffer_manager->lock_file_buffer));

  //NORMALLY NOT NECESSARY TO SYNC, DOUBLE CHECK
  //__sync_fetch_and_add(&(node->file_buffer->buffered_list_size),1);
  ++node->file_buffer->buffered_list_size; 

  
  return SUCCESS;
}
*/

enum response append_ts_to_node_parallel(struct hercules_index * index,
					 struct hercules_node * node,
					 ts_type * timeseries,
					 void * thread_data)

{

  char *cr;

  int idx = node->file_buffer->buffered_list_size;  

  int ts_length = index->settings->timeseries_size;
  int max_leaf_size = index->settings->max_leaf_size;

  if (idx == 0)
  {
      node->file_buffer->buffered_list = NULL;
      node->file_buffer->buffered_list = malloc(sizeof(struct ts_type *) * max_leaf_size);
      
      if (node->file_buffer->buffered_list == NULL)
      {
        fprintf(stderr, "Error in hercules_index.c:  Could not \
                         allocate memory for the buffered list. \n");
        return FAILURE;                  
      }
  }
  

  cr = ((index_thread_data*)thread_data)->buffer_offset + ((index_thread_data*)thread_data)->buffer_counter * (sizeof(ts_type) * ts_length);  
  node->file_buffer->buffered_list[idx] = (ts_type *) cr;
  memcpy((void *) cr, (void *) timeseries, sizeof(ts_type) * ts_length);
  ++(((index_thread_data*)thread_data)->buffer_counter);
  //printf ("hercules_node.c: inserted %g\n",   node->file_buffer->buffered_list[idx][0] );
  
  if (node->file_buffer->buffered_list[idx] == NULL)
  {
        fprintf(stderr, "Error in hercules_index.c:  Could not \
                         allocate memory for the time series in\
                         the buffer.\n");
        return FAILURE;                  
  }

  ++node->file_buffer->buffered_list_size; 
  
  return SUCCESS;
}

enum response append_ts_to_child_node(struct hercules_index * index,
				struct hercules_node * node,
				ts_type * timeseries)
{

  // fprintf(stderr, "IN APPEND TS TO CHILD NODE.\n");
  if (!get_file_buffer(index, node))
  {
    fprintf(stderr, "Error in hercules_index.c:  Could not get the \
                     file buffer for this node.\n");
    return FAILURE;              
  }

  if (node->file_buffer == NULL)
  {
    fprintf(stderr, "Error in hercules_index.c:  Null file buffer for \
                     this node after creating it.\n");
    return FAILURE;              
  }

  int idx = node->file_buffer->buffered_list_size;  

  int ts_length = index->settings->timeseries_size;
  int max_leaf_size = index->settings->max_leaf_size;

  if (idx == 0)
  {
      node->file_buffer->buffered_list = NULL;
      node->file_buffer->buffered_list = malloc(sizeof(struct ts_type *) * max_leaf_size);
      
      if (node->file_buffer->buffered_list == NULL)
      {
        fprintf(stderr, "Error in hercules_index.c:  Could not \
                         allocate memory for the buffered list. \n");
        return FAILURE;                  
      }
  }


  node->file_buffer->buffered_list[idx] = (ts_type *) index->buffer_manager->current_record;
  index->buffer_manager->current_record += sizeof(ts_type) * ts_length;
  index->buffer_manager->current_record_index++;

  
  if (node->file_buffer->buffered_list[idx] == NULL)
  {
        fprintf(stderr, "Error in hercules_index.c:  Could not \
                         allocate memory for the time series in\
                         the buffer.\n");
        return FAILURE;                  
  }
  
  for(int i=0; i<ts_length; ++i)
  {
        node->file_buffer->buffered_list[idx][i] = timeseries[i];
  }
    
  ++node->file_buffer->buffered_list_size; 
    
  return SUCCESS;
}

enum response append_ts_to_child_node_parallel(struct hercules_index * index,
					       struct hercules_node * node,
					       ts_type * timeseries,
					       void * thread_data)
{
  char *cr;
  

  int idx = node->file_buffer->buffered_list_size;  

  int ts_length = index->settings->timeseries_size;
  int max_leaf_size = index->settings->max_leaf_size;

  if (idx == 0)
  {
      node->file_buffer->buffered_list = NULL;
      node->file_buffer->buffered_list = malloc(sizeof(struct ts_type *) * max_leaf_size);
      
      if (node->file_buffer->buffered_list == NULL)
      {
        fprintf(stderr, "Error in hercules_index.c:  Could not \
                         allocate memory for the buffered list. \n");
        return FAILURE;                  
      }
  }

  
  cr = ((index_thread_data*)thread_data)->buffer_offset + ((index_thread_data*)thread_data)->buffer_counter * (sizeof(ts_type) * ts_length);  
  node->file_buffer->buffered_list[idx] = (ts_type *) cr;
  memcpy((void *) cr, (void *) timeseries, sizeof(ts_type) * ts_length);
  ++(((index_thread_data*)thread_data)->buffer_counter);
  

  if (node->file_buffer->buffered_list[idx] == NULL)
  {
        fprintf(stderr, "Error in hercules_index.c:  Could not \
                         allocate memory for the time series in\
                         the buffer.\n");
        return FAILURE;                  
  }


   ++node->file_buffer->buffered_list_size;  
 
  return SUCCESS;
}

enum response create_hercules_node_filename(struct hercules_index_settings *settings,
                                          struct hercules_node * node,
                                          struct hercules_node * parent_node)
{
    int i;
    
    node->filename = malloc(sizeof(char) * (settings->max_filename_size));

    int l = 0;
    l += sprintf(node->filename+l ,"%02d", node->num_node_points);
    
    // If this has level other than 0 then it is not a root node and as such it does have some 
    // split data on its parent.

    if (node->level) {
      if(node->is_left)
      {
  	 l += sprintf(node->filename+l ,"%s", "_L");
      }
      else
      {
  	 l += sprintf(node->filename+l ,"%s", "_R");	
      }

      //Add parent split policy, just the code 0 for mean and 1 for stdev
      struct node_split_policy *policy;
      policy = parent_node->split_policy;
      
      if(policy->indicator_split_idx)
      {
  	 l += sprintf(node->filename+l ,"%s", "_1");
      }
      else
      {
  	 l += sprintf(node->filename+l ,"%s", "_0");	
      }
      
      l += sprintf(node->filename+l ,"_(%d,%d,%g)",policy->split_from, policy->split_to,
		   policy->indicator_split_value);
      

    }
    
    // If its level is 0 then it is a root
    l += sprintf(node->filename+l ,"_%d", node->level);	      

    return SUCCESS;
}


enum response update_node_statistics(struct hercules_node * node, ts_type * timeseries)
{

  //update vertical node_segment_sketch
  for (int i = 0; i < node->num_node_points; i++) {
    int from = 0;
    int to = 0;

    from = get_segment_start(node->node_points, i);
    to = get_segment_end(node->node_points, i);

    if (!node_segment_sketch_update_sketch(&node->node_segment_sketches[i], timeseries, from,to))
    {
        fprintf(stderr, "Error in hercules_index.c:  Could not update vertical sketch for node segment.\n");
        return FAILURE;                  
    }
          
  }
  
  //update horizontal node_segment_sketch
  for (int i = 0; i < node->num_hs_node_points; i++) {
    if (!node_segment_sketch_update_sketch(&node->hs_node_segment_sketches[i], timeseries, get_segment_start(node->hs_node_points, i), get_segment_end(node->hs_node_points, i)))
    {
        fprintf(stderr, "Error in hercules_index.c:  Could not update horizontal sketch for node segment.\n");
        return FAILURE;                  
    }
    
  }

    ++node->node_size;

    return SUCCESS;
}

enum response update_node_statistics_parallel(struct hercules_node * node, ts_type * timeseries,struct segment_sketch * series_segment_sketch)
{

  //update vertical node_segment_sketch
  for (int i = 0; i < node->num_node_points; i++) {
    int from = 0;
    int to = 0;

    from = get_segment_start(node->node_points, i);
    to = get_segment_end(node->node_points, i);

    if (!node_segment_sketch_update_sketch_parallel(&node->node_segment_sketches[i], timeseries, from,to, series_segment_sketch))
    {
        fprintf(stderr, "Error in hercules_index.c:  Could not update vertical sketch for node segment.\n");
        return FAILURE;                  
    }
          
  }
  
  //update horizontal node_segment_sketch
  for (int i = 0; i < node->num_hs_node_points; i++) {
    if (!node_segment_sketch_update_sketch_parallel(&node->hs_node_segment_sketches[i], timeseries, get_segment_start(node->hs_node_points, i), get_segment_end(node->hs_node_points, i),series_segment_sketch))
    {
        fprintf(stderr, "Error in hercules_index.c:  Could not update horizontal sketch for node segment.\n");
        return FAILURE;                  
    }
    
  }

    ++node->node_size;

    return SUCCESS;
}

enum response hercules_index_node_update_synopsis_vsplit(struct hercules_node * node, ts_type * timeseries,struct segment_sketch * series_segment_sketch)
{

  short split_segment;
  while (node != NULL)
    {
        split_segment = node->split_segment;
	if (split_segment != -1)
	  {
	    int from = 0;
	    int to = 0;

	    from = get_segment_start(node->node_points, split_segment);
	    to = get_segment_end(node->node_points, split_segment);
      
        pthread_mutex_lock(&(node->lock_data));
	    if (!node_segment_sketch_update_sketch_parallel(&node->node_segment_sketches[split_segment], timeseries, from,to, series_segment_sketch))
	      {
		fprintf(stderr, "Error in hercules_index.c:  Could not update vertical sketch for node segment.\n");
		return FAILURE;                  
	      }
        pthread_mutex_unlock(&(node->lock_data));

       }
      node = node->parent;
    }
    return SUCCESS;
}


enum response update_node_ancestors_statistics_for_split_segment(struct hercules_node * node, ts_type * timeseries)
{

  short split_segment;
  while (node != NULL)
    {
        split_segment = node->split_segment;
	if (split_segment != -1)
	  {
	    int from = 0;
	    int to = 0;

	    from = get_segment_start(node->node_points, split_segment);
	    to = get_segment_end(node->node_points, split_segment);

	    if (!node_segment_sketch_update_sketch(&node->node_segment_sketches[split_segment], timeseries, from,to))
	      {
		fprintf(stderr, "Error in hercules_index.c:  Could not update vertical sketch for node segment.\n");
		return FAILURE;                  
	      }
          }
      node = node->parent;
    }
    return SUCCESS;
}

enum response update_node_ancestors_statistics_for_non_split_segments(struct hercules_node * leaf_node)
{

  struct hercules_node * currP = NULL;
  struct hercules_node * currN = NULL;  
  currP = leaf_node -> parent;
  currN = leaf_node;
  
  int from = -1;
  int to = -1;
  short split_segment;
  int idx = -1;
  
  //we only keep track of one vertically split segment for each ancestor
  //if the parent has the same vertical segmentation as the child (node->split_segment = -1)
  //the parent's stats can be deduced from those of its children

  
  while (currP != NULL)
  {
    split_segment = currP->split_segment;
    //get the vertically split segment
    for (int i = 0; i < currP->num_node_points; i++) {
      if (i != split_segment)
      {	  

	if((i < split_segment) | (split_segment == -1))	
	  idx = i;
	else if (i > split_segment)
	  idx = i+1;
	
	currP->node_segment_sketches[i].indicators[0] = fmaxf(currP->node_segment_sketches[i].indicators[0], currN->node_segment_sketches[idx].indicators[0]);
	currP->node_segment_sketches[i].indicators[1] = fminf(currP->node_segment_sketches[i].indicators[1], currN->node_segment_sketches[idx].indicators[1]);      
	currP->node_segment_sketches[i].indicators[2] = fmaxf(currP->node_segment_sketches[i].indicators[2], currN->node_segment_sketches[idx].indicators[2]);
	currP->node_segment_sketches[i].indicators[3] = fminf(currP->node_segment_sketches[i].indicators[3], currN->node_segment_sketches[idx].indicators[3]);
      }
    }

    //reset to zero the first time the parent is processed

    //currP->node_size += currN->node_size;
    
    currN = currP;
    currP = currP->parent;
  }
  
    return SUCCESS;
}

enum response hercules_index_node_update_synopsis_hsplit(struct hercules_node * leaf_node)
{

  struct hercules_node * currP = NULL;
  struct hercules_node * currN = NULL;  
  currP = leaf_node -> parent;
  currN = leaf_node;
  
  int from = -1;
  int to = -1;
  short split_segment;
  int idx = -1;
  
  //we only keep track of one vertically split segment for each ancestor
  //if the parent has the same vertical segmentation as the child (node->split_segment = -1)
  //the parent's stats can be deduced from those of its children

  
  while (currP != NULL)
  {
    split_segment = currP->split_segment;
    //get the vertically split segment
    for (int i = 0; i < currP->num_node_points; i++) {
      if (i != split_segment)
      {	  

	if((i < split_segment) | (split_segment == -1))	
	  idx = i;
	else if (i > split_segment)
	  idx = i+1;
      pthread_mutex_lock(&(currP->lock_data));
	
	currP->node_segment_sketches[i].indicators[0] = fmaxf(currP->node_segment_sketches[i].indicators[0], currN->node_segment_sketches[idx].indicators[0]);
	currP->node_segment_sketches[i].indicators[1] = fminf(currP->node_segment_sketches[i].indicators[1], currN->node_segment_sketches[idx].indicators[1]);      
	currP->node_segment_sketches[i].indicators[2] = fmaxf(currP->node_segment_sketches[i].indicators[2], currN->node_segment_sketches[idx].indicators[2]);
	currP->node_segment_sketches[i].indicators[3] = fminf(currP->node_segment_sketches[i].indicators[3], currN->node_segment_sketches[idx].indicators[3]);
      pthread_mutex_unlock(&(currP->lock_data));

      }
    }

    //reset to zero the first time the parent is processed

    //currP->node_size += currN->node_size;
    
    currN = currP;
    currP = currP->parent;
  }
  
    return SUCCESS;
}

/*

  This function calculates the euclidean distance of query to a 
  given node.
  
  The hercules and isax both load all the time series in the node 
  and compare each to the query. 

  The function returns the smallest distance between the query
  and the time series in the node.

 */
void calculate_node_knn_distance (struct hercules_index *index, struct hercules_node *node,
				  ts_type *query_ts,
				  ts_type *query_ts_reordered, int *query_order,
				  unsigned int offset, ts_type bsf,
				  unsigned int k,
				  struct query_result  *knn_results,
				  struct bsf_snapshot ** bsf_snapshots,
				  unsigned int *cur_bsf_snapshot,
				  unsigned int * cur_size, int serial)
{
     ts_type distance = FLT_MAX;

     //get the k-th distance from the results queue
     ts_type kth_bsf = FLT_MAX;

     
    //count the number of leaves and the number of time series that were checked
    //this is different from the count_loaded_node and counted_loaded_ts
    //which count the number of leaves and time series that were not found in
    //memory and had to be retrieved from disk
    //checked_nodes = loaded_nodes + nodes_in_memory    
    COUNT_CHECKED_NODE
    COUNT_CHECKED_TS(node->node_size)		  

    
    //TEST THAT DATA IS FULLY IN MEM
    //If the leaf's data is in disk, load it
    if (node->file_buffer->buffered_list_size == 0) 
    {
      COUNT_LOADED_NODE      
      COUNT_LOADED_TS(node->node_size)
      COUNT_PARTIAL_LOAD_NODE_TIME_START
 
      //printf("loaded_node = %u, file_pos = %llu\n", loaded_nodes_count,node->file_pos);
      node->file_buffer->buffered_list = get_all_time_series_in_node(index, node,serial);
      node->file_buffer->buffered_list_size = node->file_buffer->disk_count;


      if (node->file_buffer->buffered_list == NULL)
       {
            fprintf(stderr, "Error in hercules_index.c:  Could not retrieve all time series for node %s.\n", node->filename);
       }
       COUNT_PARTIAL_LOAD_NODE_TIME_END      
    }
    //If the leaf's data is in memory, proceed. A leaf's data is either fully in disk or in memory 
    double tS_bsf;
    double tE_bsf;
    struct timeval current_time_bsf;
	
    for (int idx = 0; idx < node->file_buffer->buffered_list_size; ++idx)
    {  

      struct query_result result;
      result =  knn_results[k-1];
      
      kth_bsf = result.distance;

      #ifdef __SSE__
      distance = ts_euclidean_distance_SIMD(query_ts,
					   node->file_buffer->buffered_list[idx],
					   index->settings->timeseries_size,
					   kth_bsf);
      #elif
      
      distance = ts_euclidean_distance_reordered(query_ts_reordered,
						 node->file_buffer->buffered_list[idx],
						 offset,  //offset is 0 for whole matching
						 index->settings->timeseries_size,
						 kth_bsf,
						 query_order);

      #endif
      if (distance < kth_bsf)
      {
	struct query_result object_result;// =  malloc(sizeof(struct query_result));
	object_result.node = node;
	object_result.distance =  distance;
	

	queue_bounded_sorted_insert(knn_results, object_result, cur_size, k);
	if (cur_bsf_snapshot != NULL)
	{
	  gettimeofday(&current_time_bsf, NULL);
	  tS_bsf = partial_time_start.tv_sec*1000000 + (partial_time_start.tv_usec); 
	  tE_bsf = current_time_bsf.tv_sec*1000000  + (current_time_bsf.tv_usec);	
	  for (int j = 0; j < k; ++j)
	    {
	      bsf_snapshots[j][*cur_bsf_snapshot].distance = knn_results[j].distance;
	      bsf_snapshots[j][*cur_bsf_snapshot].time =tE_bsf - tS_bsf;
	    }
	  ++(*cur_bsf_snapshot);
	}
	
      }

    }

    //clearing the data for this node
     for (int i = 0 ; i < index->settings->max_leaf_size; ++i)
     {
       free(node->file_buffer->buffered_list[i]);
     }

     free(node->file_buffer->buffered_list);

     node->file_buffer->buffered_list = NULL;
     node->file_buffer->buffered_list_size = 0;
    
}


ts_type calculate_node_distance (struct hercules_index *index, struct hercules_node *node, ts_type *query_ts_reordered, int *query_order, unsigned int offset, ts_type bsf)
{
//    ts_type bsf = FLT_MAX;
     ts_type dist = FLT_MAX;

    //count the number of leaves and the number of time series that were checked
    //this is different from the count_loaded_node and counted_loaded_ts
    //which count the number of leaves and time series that were not found in
    //memory and had to be retrieved from disk
    //checked_nodes = loaded_nodes + nodes_in_memory    
    COUNT_CHECKED_NODE
    COUNT_CHECKED_TS(node->node_size)		  


    //TEST THAT DATA IS FULLY IN MEM
    //If the leaf's data is in disk, load it
    if (node->file_buffer->buffered_list_size == 0) 
    {
      COUNT_LOADED_NODE      
      COUNT_LOADED_TS(node->node_size)
      COUNT_PARTIAL_LOAD_NODE_TIME_START
 
	node->file_buffer->buffered_list = get_all_time_series_in_node(index, node,0);
      node->file_buffer->buffered_list_size = node->file_buffer->disk_count;


      if (node->file_buffer->buffered_list == NULL)
       {
            fprintf(stderr, "Error in hercules_index.c:  Could not retrieve all time series for node %s.\n", node->filename);
       }
       COUNT_PARTIAL_LOAD_NODE_TIME_END      
    }
    //If the leaf's data is in memory, proceed. A leaf's data is either fully in disk or in memory 

    dist = calculate_ts_in_node_distance(index, node, query_ts_reordered, query_order, offset,bsf); 

    //clearing the data for this node
     for (int i = 0 ; i < index->settings->max_leaf_size; ++i)
     {
       free(node->file_buffer->buffered_list[i]);
     }

     free(node->file_buffer->buffered_list);

    node->file_buffer->buffered_list = NULL;
    node->file_buffer->buffered_list_size = 0;
    

    return dist;
}

ts_type calculate_ts_in_node_distance (struct hercules_index *index,
				       struct hercules_node *node,
				       ts_type *query_ts_reordered,
				       int * query_order,
				       unsigned int offset,
				       ts_type bound)
{
    ts_type bsf = bound;
    ts_type temp_dist;
    
    for (int idx = 0; idx < node->file_buffer->buffered_list_size; ++idx)
    {  

      temp_dist = ts_euclidean_distance_reordered(query_ts_reordered,
						node->file_buffer->buffered_list[idx],
						offset,  //offset is 0 for whole matching
						index->settings->timeseries_size,
						bsf,
						query_order);     
      if (temp_dist < bsf)
      {
        bsf = temp_dist;	
      }
    }

    return bsf;
}



ts_type calculate_node_min_distance (struct hercules_index *index, struct hercules_node *node, ts_type *query)
{
    ts_type sum = 0;
    short *points = node->node_points;
    int num_points = (int) node->num_node_points;

    
    ts_type * mean_per_segment = malloc(sizeof(ts_type) * num_points);
    ts_type * stdev_per_segment= malloc(sizeof(ts_type) * num_points);
    /*
    calc_mean_per_segment(query, points,mean_per_segment, num_points);
    calc_stdev_per_segment(query, points,stdev_per_segment, num_points);
    */
    
    #ifdef __SSE__
    calc_mean_stdev_per_segment_SIMD(query,points,mean_per_segment,stdev_per_segment,num_points);
    #elif
    calc_mean_stdev_per_segment(query,points,mean_per_segment,stdev_per_segment,num_points);
    #endif
    
    ts_type temp_dist = 0;
    ts_type temp = 0;

   for (int i=0; i < num_points; ++i)
    {
      //use mean and standard deviation to estimate the distance
        temp_dist = 0;
        temp = 0;
        if ((stdev_per_segment[i] - node->node_segment_sketches[i].indicators[2]) *
	    (stdev_per_segment[i] - node->node_segment_sketches[i].indicators[3]) > 0)
	{
          temp = fmin(fabs(stdev_per_segment[i] - node->node_segment_sketches[i].indicators[2]),
				fabs(stdev_per_segment[i] - node->node_segment_sketches[i].indicators[3]));
          temp_dist += temp * temp;
        }

        if ((mean_per_segment[i] - node->node_segment_sketches[i].indicators[0]) *
	    (mean_per_segment[i] - node->node_segment_sketches[i].indicators[1]) > 0)
	{
          temp = fmin(fabs(mean_per_segment[i] - node->node_segment_sketches[i].indicators[0]),
				fabs(mean_per_segment[i] - node->node_segment_sketches[i].indicators[1]));
          temp_dist += temp * temp;
        }
        sum += temp_dist * get_segment_length(points, i);
    }

    //sum = sqrt(sum);
    
    free(mean_per_segment);
    free(stdev_per_segment);
    return sum;
}

ts_type calculate_node_max_distance (struct hercules_index *index, struct hercules_node *node, ts_type *query)
{
    ts_type sum = 0;
    short *points = node->node_points;
    int num_points = (int) node->num_node_points;

    
    ts_type * mean_per_segment = malloc(sizeof(ts_type) * num_points);
    ts_type * stdev_per_segment= malloc(sizeof(ts_type) * num_points);

    calc_mean_per_segment(query, points,mean_per_segment, num_points);
    calc_stdev_per_segment(query, points,stdev_per_segment, num_points);


    for (int i=0; i < num_points; ++i)
      {
	//use mean and standard deviation to estimate the distance
	ts_type temp_dist = 0;
	
	temp_dist += pow(stdev_per_segment[i] + node->node_segment_sketches[i].indicators[2], 2);
	    
	ts_type ub_threshold = (node->node_segment_sketches[i].indicators[0] + node->node_segment_sketches[i].indicators[1]) / 2.0;
	
        if (mean_per_segment[i] <= ub_threshold)
	  {
	    temp_dist += pow(fabs(mean_per_segment[i] - node->node_segment_sketches[i].indicators[0]), 2);
	  }
	  else
	    {
	      temp_dist += pow(fabs(mean_per_segment[i] - node->node_segment_sketches[i].indicators[1]), 2);	  
	    }
	  sum += temp_dist * get_segment_length(points, i);
      }
    //sum = sqrt(sum);
    
    free(mean_per_segment);
    free(stdev_per_segment);
    
    if (sum == 0)
      return FLT_MAX;
    else
      return sum;
}



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
					 void* qwdata)

{
     ts_type distance = FLT_MAX;

     //get the k-th distance from the results queue
     ts_type kth_bsf = FLT_MAX;

     
    //count the number of leaves and the number of time series that were checked
    //this is different from the count_loaded_node and counted_loaded_ts
    //which count the number of leaves and time series that were not found in
    //memory and had to be retrieved from disk
    //checked_nodes = loaded_nodes + nodes_in_memory    
    COUNT_CHECKED_NODE
    COUNT_CHECKED_TS(node->node_size)		  

    double tS = 0;
    double tE = 0;
    struct timeval start_time;
    struct timeval end_time;
    
    //TEST THAT DATA IS FULLY IN MEM
    //If the leaf's data is in disk, load it
    if (node->file_buffer->buffered_list_size == 0) 
    {
      COUNT_LOADED_NODE      
      COUNT_LOADED_TS(node->node_size)
      COUNT_PARTIAL_LOAD_NODE_TIME_START
	
      #if DETAILED_STATS == 1 
      gettimeofday(&start_time, NULL);
      node->file_buffer->buffered_list = get_all_time_series_in_node(index, node,0);
      gettimeofday(&end_time, NULL);
      tS = start_time.tv_sec*1000000 + (start_time.tv_usec);
      tE = end_time.tv_sec*1000000 + (end_time.tv_usec);
      ((siss_query_worker_data*)qwdata)->thread_input_time += (tE - tS);
      #else
      node->file_buffer->buffered_list = get_all_time_series_in_node(index, node,1);
      #endif

      node->file_buffer->buffered_list_size = node->file_buffer->disk_count;

      if (node->file_buffer->buffered_list == NULL)
       {
            fprintf(stderr, "Error in hercules_index.c:  Could not retrieve all time series for node %s.\n", node->filename);
       }
       COUNT_PARTIAL_LOAD_NODE_TIME_END      
    }
    //If the leaf's data is in memory, proceed. A leaf's data is either fully in disk or in memory 
    double tS_bsf;
    double tE_bsf;
    struct timeval current_time_bsf;

    #if DETAILED_STATS == 1
    gettimeofday(&start_time, NULL);
    #endif
    for (int idx = 0; idx < node->file_buffer->buffered_list_size; ++idx)
    {  

      struct query_result result;
      result =  knn_results[k-1];
      
      kth_bsf = result.distance;
      /*
      distance = ts_euclidean_distance_reordered(query_ts_reordered,_psq
						 node->file_buffer->buffered_list[idx],
						 offset,  //offset is 0 for whole matching
						 index->settings->timeseries_size,
						 kth_bsf,
						 query_order);

      */
     #ifdef __SSE__
     distance = ts_euclidean_distance_SIMD(query_ts,
					   node->file_buffer->buffered_list[idx],
					   index->settings->timeseries_size,
					   kth_bsf);
     #elif
     distance = ts_euclidean_distance_reordered(query_ts_reordered,_psq
						node->file_buffer->buffered_list[idx],
						offset,  //offset is 0 for whole matching
						index->settings->timeseries_size,
						kth_bsf,
						query_order);
     #endif
     if (distance < kth_bsf)
      {
	struct query_result object_result;// =  malloc(sizeof(struct query_result));
	object_result.node = node;
	object_result.distance =  distance;
	
	//no need to lock during approximate search
     if (qwdata != NULL)
	{
	  pthread_rwlock_wrlock(((siss_query_worker_data*)qwdata)->lock_bsf);
	  queue_bounded_sorted_insert(knn_results, object_result, cur_size, k);
	  pthread_rwlock_unlock(((siss_query_worker_data*)qwdata)->lock_bsf);
	}
	else
	{
	  queue_bounded_sorted_insert(knn_results, object_result, cur_size, k);	  
	}



	if (cur_bsf_snapshot != NULL)
	{
	  gettimeofday(&current_time_bsf, NULL);
	  tS_bsf = partial_time_start.tv_sec*1000000 + (partial_time_start.tv_usec); 
	  tE_bsf = current_time_bsf.tv_sec*1000000  + (current_time_bsf.tv_usec);	
	  for (int j = 0; j < k; ++j)
	    {
	      bsf_snapshots[j][*cur_bsf_snapshot].distance = knn_results[j].distance;
	      bsf_snapshots[j][*cur_bsf_snapshot].time =tE_bsf - tS_bsf;
	    }
	  ++(*cur_bsf_snapshot);
	}
	
      }

    }
    #if DETAILED_STATS == 1
      gettimeofday(&end_time, NULL);
      tS = start_time.tv_sec*1000000 + (start_time.tv_usec);
      tE = end_time.tv_sec*1000000 + (end_time.tv_usec);
      ((siss_query_worker_data*)qwdata)->thread_realdist_time += (tE - tS);
    #endif
      
      //clearing the data for this node
     for (int i = 0 ; i < index->settings->max_leaf_size; ++i)
     {
       free(node->file_buffer->buffered_list[i]);
     }

     free(node->file_buffer->buffered_list);

     node->file_buffer->buffered_list = NULL;
     node->file_buffer->buffered_list_size = 0;
    
}
struct segment_sketch * init_segment_sketches(int num_child_segments, int num_indicators)
{
  struct segment_sketch * child_node_segment_sketches = malloc (sizeof(struct segment_sketch) *
					num_child_segments);	  


  if(child_node_segment_sketches == NULL)
    {
      fprintf(stderr,"Error in hercules_node_split.c: could not allocate "
                            "memory for the child node segment sketches");
      return FAILURE;
    }


  for (int k = 0; k< num_child_segments;++k)
    {
      child_node_segment_sketches[k].indicators = NULL;
      child_node_segment_sketches[k].indicators =malloc (sizeof(ts_type) * num_indicators);
      if(child_node_segment_sketches[k].indicators == NULL)
	{
	  fprintf(stderr,"Error in hercules_node_split.c: could not allocate "
                               "memory for the child node segment sketches");
	  return FAILURE;
	}
      child_node_segment_sketches[k].indicators[0] = -FLT_MAX ; //for max mean
      child_node_segment_sketches[k].indicators[1] = FLT_MAX; //for min mean
      child_node_segment_sketches[k].indicators[2] = -FLT_MAX; //for max stdev
      child_node_segment_sketches[k].indicators[3] = FLT_MAX; //for min stdev
      child_node_segment_sketches[k].num_indicators = num_indicators;
      child_node_segment_sketches[k].lock_synopsis = malloc(sizeof(pthread_mutex_t)*4);
      for(int l =0; l < 4 ;l++)
	{
	  pthread_mutex_init(&(child_node_segment_sketches[k].lock_synopsis[l]), NULL);
	}	    
	    
    }
  return child_node_segment_sketches;
}

enum response hercules_index_node_update_synopsis(struct hercules_index *index, struct hercules_node *node, int sims,void * tdata)
{

  struct segment_sketch timeseries_segment_sketch;
  timeseries_segment_sketch.indicators = NULL;
  timeseries_segment_sketch.indicators = malloc (sizeof(ts_type)*2);
  if(timeseries_segment_sketch.indicators == NULL) {
        fprintf(stderr,"Error in hercules_index.c: Could not allocate memory for" 
                        "series segment sketch indicators.\n");
  }

 timeseries_segment_sketch.num_indicators = 2;

 int full_size = strlen(index->settings->root_directory) + strlen(node->filename)+1;


    int ts_length = index->settings->timeseries_size;
    unsigned int node_size;    
    node_size = node->node_size;

    ts_type ** ts_list;
    sax_type ** sax_list;

    int num_segments = index->settings->paa_segments;    

    sax_list = ((index_thread_data*)tdata)->node_sax_data;	
    get_all_time_series_in_node_parallel(index, node,0,tdata);
    ts_list = ((index_thread_data*)tdata)->split_node_data;


    // ++index->leaves_approx_pos;       
    //COUNT_PARTIAL_OUTPUT_TIME_START
    for (int idx=0; idx < node_size;++idx)
    {
       if (sims)
     	 {
	        sax_from_ts(ts_list[idx], sax_list[idx], index->settings->ts_values_per_paa_segment,
  		       index->settings->paa_segments, index->settings->sax_alphabet_cardinality,
		       index->settings->sax_bit_cardinality);
	   }

       hercules_index_node_update_synopsis_vsplit(node, ts_list[idx],&timeseries_segment_sketch);       
    }

     hercules_index_node_update_synopsis_hsplit(node);       


   //COUNT_PARTIAL_OUTPUT_TIME_END
   

    
	    const char *full_filename = malloc(sizeof(char) * full_size);
	    full_filename = strcpy(full_filename, index->settings->root_directory);
	    full_filename = strcat(full_filename, node->filename);
	    full_filename = strcat(full_filename, "\0");
            remove(full_filename);
        free(full_filename);

    free(timeseries_segment_sketch.indicators);   
     //printf("thread = %d leaf %d first_value %g\n", ((index_thread_data*)tdata)->thread_id, ((index_thread_data*)tdata)->current_leaf, ts_list[0][0]);

    return SUCCESS; 
    
}


