//
//
//  hercules_file_buffer.c
//  Created by Karima Echihabi on 18/12/2018
//

#include "../config.h"
#include "../globals.h"
#include <stdio.h>
#include <stdlib.h>
#include "../include/hercules_file_buffer.h"
#include "../include/hercules_file_buffer_manager.h"
#include "../include/hercules_index.h"
#include "../include/calc_utils.h"
#include <pthread.h>


enum response hercules_file_buffer_init(struct hercules_node *node)
{

  node->file_buffer = NULL;
  node->file_buffer = malloc(sizeof(struct hercules_file_buffer));

  if(node->file_buffer == NULL) {
    fprintf(stderr,"Error in hercules_file_buffer.c: Could not allocate memory for file buffer.\n");    
    return FAILURE;	
  }
  
  node->file_buffer->in_disk = false;
  node->file_buffer->disk_count = 0;
  
  node->file_buffer->buffered_list = NULL;
  node->file_buffer->buffered_list_size = 0;

  node->file_buffer->node = node;
  node->file_buffer->position_in_map = NULL;
  
  node->file_buffer->do_not_flush = false;
  
  return SUCCESS;
  
}
/* This function has to be improved by using mmap and reading the full block of data series once instead of a loop*/
ts_type ** get_all_time_series_in_node(struct hercules_index * index, struct hercules_node * node, int serial) 
{

  int ts_length = index->settings->timeseries_size;
  int max_leaf_size = index->settings->max_leaf_size;
  ts_type ** ret = NULL;
  ret = calloc(max_leaf_size, sizeof(ts_type *));


  if (node->file_buffer != NULL)
    {
      
      if (node->file_buffer->disk_count > 0)
	{
      
	  if(node->filename == NULL)
	    {   
	      fprintf(stderr, "Error in hercules_file_buffer.c: This node has data on disk but"
		      "could not get its filename.\n");
	      //return FAILURE;
	    }
      
	  if (!serial){    
	    int full_size = strlen(index->settings->root_directory) + strlen(node->filename)+1;
	
	    const char *full_filename = malloc(sizeof(char) * full_size);
	    full_filename = strcpy(full_filename, index->settings->root_directory);
	    full_filename = strcat(full_filename, node->filename);
	    full_filename = strcat(full_filename, "\0");
	
	    COUNT_PARTIAL_RAND_INPUT
	      COUNT_PARTIAL_INPUT_TIME_START
	      FILE *ts_file = fopen(full_filename, "r");
	    COUNT_PARTIAL_INPUT_TIME_END
       
	      if(ts_file == NULL)
		{   
		  fprintf(stderr, "Error in hercules_file_buffer.c: Could not open"
			  "the filename %s. Reason = %s\n", full_filename, strerror(errno));
		}
     
	    for (int i=0; i<node->file_buffer->disk_count;++i ) 
	      {       
		ret[i] = calloc(ts_length, sizeof(ts_type));
	    	    
		COUNT_PARTIAL_SEQ_INPUT
		  COUNT_PARTIAL_INPUT_TIME_START      
		  fread(ret[i],
			sizeof(ts_type),ts_length, ts_file);
		COUNT_PARTIAL_INPUT_TIME_END       
		  }
	    COUNT_PARTIAL_INPUT_TIME_START
	      if(fclose(ts_file))
		{   
		  fprintf(stderr, "Error in hercules_file_buffer.c: Could not close"
			  "the filename %s. Reason= %s.\n", full_filename, strerror(errno));
		} 
	    COUNT_PARTIAL_INPUT_TIME_END    
	      free(full_filename);
	  }
	  else
	    {
	      index->leaves_raw_file = fopen(index->leaves_raw_filename, "rb");
	  	  
	      fseek(index->leaves_raw_file, node->file_pos * index->settings->timeseries_size * sizeof(ts_type), SEEK_SET);

	      for (int i=0; i< node->file_buffer->disk_count;++i ) 
		{       
		  ret[i] = calloc(ts_length, sizeof(ts_type));
	    
		  COUNT_PARTIAL_SEQ_INPUT
		    COUNT_PARTIAL_INPUT_TIME_START      
		    fread(ret[i],sizeof(ts_type),ts_length, index->leaves_raw_file);
		  COUNT_PARTIAL_INPUT_TIME_END       
		    }       
	      fclose(index->leaves_raw_file);
	    }
       
  
	  int idx = node->file_buffer->buffered_list_size;
 
	  for (int i = 0 ; i < idx; ++i)
	    {
	      ret[i+(node->file_buffer->disk_count)] = calloc(ts_length, sizeof(ts_type));
	      for(int j=0; j<ts_length; ++j)
		{
		  ret[i+(node->file_buffer->disk_count)][j] = node->file_buffer->buffered_list[i][j];
		}
	    }
      

	}
      else
	{
	  int idx = node->file_buffer->buffered_list_size;

	  for (int i = 0 ; i < idx; ++i)
	    {
	      ret[i] = calloc(ts_length, sizeof(ts_type));
	      for(int j=0; j<ts_length; ++j)
		{
		  ret[i][j] = node->file_buffer->buffered_list[i][j];
		}
	    }    
	}
    }

  return ret;
}


void * get_all_time_series_in_node_parallel(struct hercules_index * index, struct hercules_node * node, int serial, void * tdata) 
{

  int ts_length = index->settings->timeseries_size;
  int max_leaf_size = index->settings->max_leaf_size;
  
  ts_type ** ret = ((index_thread_data*)tdata)->split_node_data;
  
  double tS = 0;
  double tE = 0;
  struct timeval start_time;
  struct timeval end_time;
 
  if (node->file_buffer != NULL)
    {
      if (node->file_buffer->disk_count > 0)
	{
#if DETAILED_STATS == 1 
	  gettimeofday(&start_time, NULL);
#endif


	  if(node->filename == NULL)
	    {   
	      fprintf(stderr, "Error in hercules_file_buffer.c: This node has data on disk but"
		      "could not get its filename.\n");
	    }
      
	  if (!serial){    
	    int full_size = strlen(index->settings->root_directory) + strlen(node->filename)+1;
	
	    const char *full_filename = malloc(sizeof(char) * full_size);
	    full_filename = strcpy(full_filename, index->settings->root_directory);
	    full_filename = strcat(full_filename, node->filename);
	    full_filename = strcat(full_filename, "\0");
	
            FILE *ts_file = fopen(full_filename, "r");
       
	    if(ts_file == NULL)
	      {   
		fprintf(stderr, "Error in hercules_file_buffer.c: Could not open"
			"the filename %s. Reason = %s\n", full_filename, strerror(errno));
	      }
	      
	    for (int i=0; i<node->file_buffer->disk_count;++i ) 
	      {       
		fread(ret[i],
		      sizeof(ts_type),ts_length, ts_file);
	      }
	    if (fclose(ts_file))
	      {   
		fprintf(stderr, "Error in hercules_file_buffer.c: Could not close"
			"the filename %s. Reason= %s.\n", full_filename, strerror(errno));
	      } 
	    free(full_filename);
	  }
	  else
	    {
	      index->leaves_raw_file = fopen(index->leaves_raw_filename, "rb");
	  	  
	      fseek(index->leaves_raw_file, node->file_pos * index->settings->timeseries_size * sizeof(ts_type), SEEK_SET);
	
	      for (int i=0; i< node->file_buffer->disk_count;++i ) 
		{       
	    
		  COUNT_PARTIAL_SEQ_INPUT
		    COUNT_PARTIAL_INPUT_TIME_START      
		    fread(ret[i],sizeof(ts_type),ts_length, index->leaves_raw_file);
		  COUNT_PARTIAL_INPUT_TIME_END       
		    }       
	      fclose(index->leaves_raw_file);
	    }
	  
#if DETAILED_STATS == 1
	  gettimeofday(&end_time, NULL);
	  tS = start_time.tv_sec*1000000 + (start_time.tv_usec);
	  tE = end_time.tv_sec*1000000 + (end_time.tv_usec);
	  ((index_thread_data*)tdata)->thread_split_input_time += (tE - tS);
	  gettimeofday(&start_time, NULL); 
#endif       
  
	  int idx = node->file_buffer->buffered_list_size;
 
	  for (int i = 0 ; i < idx; ++i)
	    {
	      for(int j=0; j<ts_length; ++j)
		{
		  ret[i+(node->file_buffer->disk_count)][j] = node->file_buffer->buffered_list[i][j];
		}
	    }
	}
      else
	{
	  int idx = node->file_buffer->buffered_list_size;

	  for (int i = 0 ; i < idx; ++i)
	    {
	      for(int j=0; j<ts_length; ++j)
		{
		  ret[i][j] = node->file_buffer->buffered_list[i][j];
		}
	    }    
	}
    }
}



enum response flush_buffer_to_disk(struct hercules_index *index, struct hercules_node *node)
{

  if (node->file_buffer != NULL)
    {
      if (node->file_buffer->buffered_list_size > 0 )
	{

	  if(node->filename == NULL)
	    {   
	      fprintf(stderr, "Error in hercules_file_buffer.c: Cannot flush the node to disk. "
		      "It does not have a filename. \n");
	      return FAILURE;
	    }
  
	  int full_size = strlen(index->settings->root_directory) + strlen(node->filename)+1;
     
	  const char *full_filename = malloc(sizeof(char) * full_size);
	  full_filename = strcpy(full_filename, index->settings->root_directory);
	  full_filename = strcat(full_filename, node->filename);
	  full_filename = strcat(full_filename, "\0");

	  COUNT_PARTIAL_RAND_OUTPUT
	    COUNT_PARTIAL_OUTPUT_TIME_START
	    FILE *ts_file = fopen(full_filename, "a");
	  COUNT_PARTIAL_OUTPUT_TIME_END
       
	    if(ts_file == NULL)
	      {   
		fprintf(stderr, "Error in hercules_file_buffer.c: Flushing node to disk.."
			"Could not open the filename %s. Reason= %s\n", node->filename, strerror(errno));
		return SUCCESS;
	      }

	  int num_ts =  node->file_buffer->buffered_list_size;

	  COUNT_PARTIAL_OUTPUT_TIME_START	
	    for (int idx = 0; idx < num_ts;++idx )
	      {
		COUNT_PARTIAL_SEQ_OUTPUT
		  if(!fwrite(node->file_buffer->buffered_list[idx], sizeof(ts_type), index->settings->timeseries_size, ts_file))
		    {   
		      fprintf(stderr, "Error in hercules_file_buffer.c: Could not "
			      "write the timeseries to file %s.\n", full_filename);
		      return FAILURE;
		    }

	      }

	  if(fclose(ts_file))
	    {   
	      fprintf(stderr, "Error in hercules_file_buffer.c: Flushing node to disk.. "
		      "Could not close the filename %s. Reason = %s.\n", full_filename, strerror(errno));
	      return FAILURE;
	    }
	  COUNT_PARTIAL_OUTPUT_TIME_END
      
	    node->file_buffer->disk_count += num_ts;

	  if (!clear_file_buffer(index, node))
	    {
	      fprintf(stderr, "Error in hercules_file_buffer.c: Flushing node to disk.. "
		      "Could not clear the buffer for %s.\n", full_filename);
	      return FAILURE;      
	    }
    
	  node->file_buffer->in_disk = true;
    
	  free(full_filename);
	}
    }
  return SUCCESS; 
    
}

enum response flush_buffer_to_disk_parallel(struct hercules_index *index, struct hercules_node *node)
{

  if (node->file_buffer != NULL)
    {
      if (node->file_buffer->buffered_list_size > 0 )
	{

	  if(node->filename == NULL)
	    {   
	      fprintf(stderr, "Error in hercules_file_buffer.c: Cannot flush the node to disk. "
		      "It does not have a filename. \n");
	      return FAILURE;
	    }
  
	  int full_size = strlen(index->settings->root_directory) + strlen(node->filename)+1;
     
	  const char *full_filename = malloc(sizeof(char) * full_size);
	  full_filename = strcpy(full_filename, index->settings->root_directory);
	  full_filename = strcat(full_filename, node->filename);
	  full_filename = strcat(full_filename, "\0");

	  FILE *ts_file = fopen(full_filename, "a");
       
	  if(ts_file == NULL)
	    {   
	      fprintf(stderr, "Error in hercules_file_buffer.c: Flushing node to disk.."
		      "Could not open the filename %s. Reason= %s\n", node->filename, strerror(errno));
	      return SUCCESS;
	    }

	  int num_ts =  node->file_buffer->buffered_list_size;

	  for (int idx = 0; idx < num_ts;++idx )
	    {
	      if(!fwrite(node->file_buffer->buffered_list[idx], sizeof(ts_type), index->settings->timeseries_size, ts_file))
		{   
		  fprintf(stderr, "Error in hercules_file_buffer.c: Could not "
			  "write the timeseries to file %s.\n", full_filename);
		  return FAILURE;
		}

	    }

	  if(fclose(ts_file))
	    {   
	      fprintf(stderr, "Error in hercules_file_buffer.c: Flushing node to disk.. "
		      "Could not close the filename %s. Reason = %s.\n", full_filename, strerror(errno));
	      return FAILURE;
	    }
      
	  node->file_buffer->disk_count += num_ts;

	  if (!clear_file_buffer(index, node))
	    {
	      fprintf(stderr, "Error in hercules_file_buffer.c: Flushing node to disk.. "
		      "Could not clear the buffer for %s.\n", full_filename);
	      return FAILURE;      
	    }
    
	  node->file_buffer->in_disk = true;
    
	  free(full_filename);
	}
    }
  return SUCCESS; 
    
}

enum response clear_file_buffer(struct hercules_index *index, struct hercules_node * node)
{

  if ((node->file_buffer) == NULL )
    {
      fprintf(stderr, "Error in hercules_file_buffer.c: Cannot clear a NULL buffer.\n");
      return FAILURE;          
    }
  else
    {
      if (node->file_buffer->buffered_list != NULL)
	{
	  free(node->file_buffer->buffered_list);
	}
      
      node->file_buffer->buffered_list = NULL;
      node->file_buffer->buffered_list_size = 0;
    }

  return SUCCESS;  
}


enum response delete_file_buffer(struct hercules_index * index,struct hercules_node * node)
{
  pthread_mutex_lock(&(index->buffer_manager->lock_file_map));    
  if (node->file_buffer->in_disk) //delete file if in disk
    {    
      int full_size = strlen(index->settings->root_directory) + strlen(node->filename)+1;
     
      const char *full_filename = malloc(sizeof(char) * full_size);
      full_filename = strcpy(full_filename, index->settings->root_directory);
      full_filename = strcat(full_filename, node->filename);
      full_filename = strcat(full_filename, "\0");

      if(!remove(full_filename)) //file deleted successfully
	{
	  node->file_buffer->disk_count = 0;
	  node->file_buffer->in_disk = false;
	}
      else 
	{
	  fprintf(stderr, "Error in hercules_file_buffer.c: Error deleting filename %s.\n", full_filename);
	  return FAILURE;          
	}
      free(full_filename);  
    }
  
  struct hercules_file_map * res = node->file_buffer->position_in_map;
    
  if (res != NULL)
    {
      if (res->prev == NULL) //first element in file map 
	{ 
	  index->buffer_manager->file_map = res->next;
	  if(res->next != NULL) //deleting the first and there are others elements in map
	    {
	      res->next->prev = NULL;
	    }
	  else  //deleting first and only element
	    {
	      index->buffer_manager->file_map_tail = NULL;
	    }
	}
      else if (res->next == NULL) //deleting the last element in the map
	{
	  res->prev->next = NULL;
	  index->buffer_manager->file_map_tail = res->prev;   
	}
      else
	{
	  res->prev->next = res->next;
	  res->next->prev = res->prev;      
	}

      free(res);
      res = NULL;
      --index->buffer_manager->file_map_size;
    
      if (!clear_file_buffer(index, node))
	{
	  fprintf(stderr, "Error in hercules_file_buffer.c: Deleting node.. "
		  "Could not clear the buffer for %s.\n", node->filename);
	  return FAILURE;      
	}

    }

  free(node->filename);
  node->filename = NULL;  

  free(node->file_buffer);
  node->file_buffer = NULL;
  pthread_mutex_unlock(&(index->buffer_manager->lock_file_map));
  return SUCCESS;  
}


enum response delete_file_buffer_parallel(struct hercules_index * index,struct hercules_node * node)
{

  if (node->file_buffer->in_disk) //delete file if in disk
    {    
      int full_size = strlen(index->settings->root_directory) + strlen(node->filename)+1;
     
      const char *full_filename = malloc(sizeof(char) * full_size);
      full_filename = strcpy(full_filename, index->settings->root_directory);
      full_filename = strcat(full_filename, node->filename);
      full_filename = strcat(full_filename, "\0");

      if(!remove(full_filename)) //file deleted successfully
	{
	  node->file_buffer->disk_count = 0;
	  node->file_buffer->in_disk = false;
	}
      else 
	{
	  fprintf(stderr, "Error in hercules_file_buffer.c: Error deleting filename %s.\n", full_filename);
	  return FAILURE;          
	}
      free(full_filename);  
    }

  free(node->filename);
  node->filename = NULL;  

  free(node->file_buffer);
  node->file_buffer = NULL;

  return SUCCESS;  
}

enum response flush_leaf_to_leaves_file_update_stats_serial(struct hercules_index *index, struct hercules_node *node, int sims)
{

  ts_type ** ts_list;

  unsigned int num_segments = index->settings->paa_segments;    
  sax_type *sax = NULL;
  sax = malloc(sizeof(sax_type) * num_segments);

  ts_list = get_all_time_series_in_node(index, node,0);

  for (int idx=0; idx < node->node_size;++idx)
    {
      COUNT_PARTIAL_SEQ_OUTPUT
	COUNT_PARTIAL_OUTPUT_TIME_START
	fwrite(ts_list[idx], sizeof(ts_type), index->settings->timeseries_size, index->leaves_raw_file);
      COUNT_PARTIAL_OUTPUT_TIME_END
	++index->leaves_raw_pos;
      if (sims)
	{
	  sax_from_ts(ts_list[idx], sax, index->settings->ts_values_per_paa_segment,
		      index->settings->paa_segments, index->settings->sax_alphabet_cardinality,
		      index->settings->sax_bit_cardinality);
	  COUNT_PARTIAL_OUTPUT_TIME_START
	    fwrite(sax, sizeof(sax_type), num_segments, index->leaves_sims_file);	   
	  COUNT_PARTIAL_OUTPUT_TIME_END
	    }

      update_node_ancestors_statistics_for_split_segment(node, ts_list[idx]);       
    }

  update_node_ancestors_statistics_for_non_split_segments(node);       

  
  for (int i = 0 ; i < index->settings->max_leaf_size; ++i)
    {
      free(ts_list[i]);
    }

  free(ts_list);
  free(sax);
    
  return SUCCESS; 
    
}


enum response hercules_index_flush_leaves(struct hercules_index *index, int num_threads)
{

  unsigned int max_leaf_size = index->settings->max_leaf_size;
  int sanity_counter= 1;
  int fin_number= index->stats->leaf_nodes_count;
  int i;
 
         
  index_thread_data *input_data=malloc(sizeof(index_thread_data)*(num_threads-1));

  pthread_t threadid[num_threads-1];
  pthread_barrier_t flush_barrier;
  pthread_barrier_init(&flush_barrier, NULL, num_threads); 
       
  fseek( index->leaves_raw_file, 0, SEEK_SET);
  fseek( index->leaves_sims_file, 0, SEEK_SET);


  for (i = 0; i < (num_threads-1); i++)
    {
      input_data[i].index=index;
      input_data[i].flush_barrier=&flush_barrier;
      input_data[i].finished=false;
      input_data[i].fin_number=fin_number;
      input_data[i].sanity_counter=&sanity_counter;
      input_data[i].threads_data = input_data;
      input_data[i].thread_id = i;
      input_data[i].node_sax_data = calloc(max_leaf_size, sizeof(sax_type *));
      input_data[i].split_node_data = calloc(max_leaf_size, sizeof(ts_type *));
      input_data[i].current_leaf = -1;

      for (int j =0; j < max_leaf_size ; ++j)
	{ 
	  input_data[i].node_sax_data[j] =  calloc(index->settings->paa_segments, sizeof(sax_type));
	  input_data[i].split_node_data[j] =  calloc(index->settings->timeseries_size, sizeof(ts_type));
	}
    }

  for (i = 0; i < (num_threads-1); i++)
    {
      pthread_create(&(threadid[i]),NULL,hercules_index_flush_leaf_worker,(void*)&(input_data[i]));
    }
       
       
  for (int i =0; i< fin_number;++i)
    {
      while (1)
	{  
	  if (__sync_fetch_and_add(&(index->leaves[i]->proc_finished),0)){
	    COUNT_PARTIAL_OUTPUT_TIME_START
	      for (int idx =0; idx < index->leaves[i]->node_size;++idx)
		{
		  fwrite(input_data[index->leaves[i]->thread_id].split_node_data[idx], 
			 sizeof(ts_type), index->settings->timeseries_size, index->leaves_raw_file);
		  fwrite(input_data[index->leaves[i]->thread_id].node_sax_data[idx], 
			 sizeof(sax_type), index->settings->paa_segments, index->leaves_sims_file);
		} 

	    COUNT_PARTIAL_OUTPUT_TIME_END
	      __sync_fetch_and_add(&(index->leaves[i]->write_finished),1);
	    break;
	  }

	}
    }
       
  pthread_barrier_wait(&flush_barrier);

  for (i = 0; i < (num_threads-1); i++)
    {
      for (int j =0; j < max_leaf_size ; ++j)
	{ 
	  free(input_data[i].node_sax_data[j]);
	  free(input_data[i].split_node_data[j]);
	}
      
      free(input_data[i].node_sax_data);
      free(input_data[i].split_node_data);

      pthread_join(threadid[i],NULL);
    }

  free(input_data);

  pthread_barrier_destroy(&flush_barrier);

}

void * hercules_index_flush_leaf_worker(void *transferdata)
{


  int num_leaves = ((index_thread_data*)transferdata)->fin_number;
  int thread_id = ((index_thread_data*)transferdata)->thread_id;
  struct hercules_index * index = ((index_thread_data*)transferdata)->index;
  struct hercules_node * node = NULL;

  int current_leaf;
  int cnt = 0;

  while (1)
    {
      ((index_thread_data*)transferdata)->current_leaf = __sync_fetch_and_add(((index_thread_data*)transferdata)->sanity_counter,1) - 1 ; 
      if( ((index_thread_data*)transferdata)->current_leaf   >= num_leaves){
	break;
      }	    
      node =   index->leaves[((index_thread_data*)transferdata)->current_leaf];
      index->leaves[((index_thread_data*)transferdata)->current_leaf]->thread_id = thread_id;
      hercules_index_node_update_synopsis(index,node,1, transferdata);
      
      __sync_fetch_and_add(&(node->proc_finished),1);

      while (1)
	{  

	  if (__sync_fetch_and_add(&(node->write_finished),0)){
	    break;
	  }
	}
      

    }
  pthread_barrier_wait((((index_thread_data*)transferdata)->flush_barrier));

}



