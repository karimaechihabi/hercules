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

//CORRECT STATS FOR ASCII
enum response hercules_query_ascii_file(const char *ifilename, int q_num, 
                  const char delimiter, struct hercules_index *index,
		  float minimum_distance, ts_type epsilon, ts_type r_delta)
{
    FILE * ifile;
    COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_INPUT_TIME_START
    ifile = fopen (ifilename,"r");
    COUNT_PARTIAL_INPUT_TIME_END    
    if (ifile == NULL) {
        fprintf(stderr, "Error in hercules_file_loaders.c: Could not open file %s!\n", ifilename);
        return FAILURE;
    }
    
    char *ts_str = NULL; 
    size_t linecap = 0;
    ssize_t linelen;
    unsigned int q_loaded = 0;
    ts_type * ts = malloc(sizeof(ts_type) * index->settings->timeseries_size);
    if (ts == NULL) {
        fprintf(stderr, "Error in hercules_file_loaders.c: Querying..\
                         Could not allocate memory for time series!\n");
        return FAILURE;
    }

    
    COUNT_PARTIAL_INPUT_TIME_START
    while ((linelen = getline(&ts_str, &linecap, ifile)) > 0 && q_loaded < q_num)
    {
        COUNT_PARTIAL_INPUT_TIME_END
	COUNT_PARTIAL_SEQ_INPUT
        if(ts_str == NULL)
        {
          fprintf(stderr,"Error in hercules_file_loaders.c: Querying..\
                         Could not get the time series from file %s.\n", ifilename);
          return FAILURE;	
        }
        if (!ts_parse_str(ts_str, ts, index->settings->timeseries_size, &delimiter))
        { 
           fprintf(stderr, "Error in hercules_file_loaders.c:  Querying..Could not parse the time series.\n");
           return FAILURE;              
        }
        
	  printf("\n\n");

        q_loaded++;
        COUNT_PARTIAL_INPUT_TIME_START    
    }



    if(fclose(ifile))
    {   
        fprintf(stderr, "Error in hercules_file_loaders.c: Could not close the query filename %s", ifilename);
        return FAILURE;
    }
    COUNT_PARTIAL_INPUT_TIME_END
      
    free(ts);
    free(ts_str);    
    return SUCCESS;    
}

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

        //answer each kNN query incrementally
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
	else if (incremental)
	{
	  exact_incr_knn_search(query_ts, query_ts_reordered, query_order, offset,
			   index, minimum_distance, epsilon, r_delta, k,
				q_loaded, ifilename, nprobes, q_settings);		  
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
	      /*
	      exact_de_knn_search_parallel_mqueue(query_ts, query_ts_reordered, query_order, offset,
						  index, minimum_distance, epsilon, r_delta, k,
						  q_loaded, ifilename, MAXQUERYTHREADS,N_PQUEUE);	      
	      */
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

enum response hercules_tlb_binary_file(const char *ifilename, int q_num, struct hercules_index *index,
						   float minimum_distance)
{

    FILE * ifile;

    ifile = fopen (ifilename,"rb");

    if (ifile == NULL) {
        fprintf(stderr, "File %s not found!\n", ifilename);
        exit(-1);
    }
    
    fseek(ifile, 0L, SEEK_END);
    file_position_type sz = (file_position_type) ftell(ifile);
    fseek(ifile, 0L, SEEK_SET);
    unsigned int ts_length = index->settings->timeseries_size;
    file_position_type total_records = sz/ts_length * sizeof(ts_type);
    unsigned int offset = 0;

    if (total_records < q_num) {
        fprintf(stderr, "File %s has only %llu records!\n", ifilename, total_records);
        exit(-1);
    }

    unsigned int q_loaded = 0;
    ts_type * query_ts = calloc(1,sizeof(ts_type) * ts_length);
      
    while (q_loaded < q_num)
    {
        total_tlb = 0;
	total_ts_count = 0;
	leaf_nodes_count = 0;

	fread(query_ts, sizeof(ts_type), ts_length, ifile);
      
	hercules_calc_tlb(query_ts, index, index->first_node);
	
        q_loaded++;
        print_tlb_stats(index,q_loaded, ifilename);

    }
      
     free(query_ts);

     if(fclose(ifile))
       {   
	 fprintf(stderr, "Error in hercules_file_loaders.c: Could not close the query filename %s", ifilename);
	 return FAILURE;
       }
     
     return SUCCESS;
}

enum response hercules_index_ascii_file(const char *ifilename, file_position_type ts_num, 
                           const char delimiter,struct  hercules_index *index)
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
    COUNT_INPUT_TIME_START
    ifile = fopen (ifilename,"r");
    COUNT_INPUT_TIME_END
    if (ifile == NULL) {
        fprintf(stderr, "File %s not found!\n", ifilename);
        exit(-1);
    }
    
	
    char *ts_str = NULL; //= malloc(sizeof(char) * 2000);
    size_t linecap = 0;
    ssize_t linelen;
    
    file_position_type ts_loaded = 0;    
    

    int percentage = 100;
    if (percentage > 100)
    {
       percentage = (int) (ts_num / (file_position_type) 100);
    }


    COUNT_INPUT_TIME_START
    while ((linelen = getline(&ts_str, &linecap, ifile)) > 0 && ts_loaded<ts_num)
    {
        COUNT_INPUT_TIME_END
#ifndef DEBUG
#if VERBOSE_LEVEL == 2
        printf("\r\x1b[32mLoading: \x1b[36m%d\x1b[0m",(ts_loaded + 1));
#endif
#endif
        COUNT_PARSE_TIME_START
        ts_parse_str(ts_str, ts, index->settings->timeseries_size, &delimiter);
        COUNT_PARSE_TIME_END


	if (!hercules_index_insert(index,ts))
	{
           fprintf(stderr, "Error in hercules_file_loaders.c:  Could not \
                           add the time series to the index.\n");
           return FAILURE;              
        }
    
        ts_loaded++;

	if(ts_loaded % percentage == 0)
            {
                float distance = 0;
            }

        COUNT_INPUT_TIME_START  
     }
    
    free(ts_str);
    free(ts);
    COUNT_INPUT_TIME_START
	fclose(ifile);
    COUNT_INPUT_TIME_END

      return SUCCESS;
}



 /*
enum response hercules_index_binary_file_p(const char *ifilename, file_position_type ts_num,
					 struct hercules_index *index, int num_threads, int start_parallel)

{
    double parse_time = 0;
    int ts_length = index->settings->timeseries_size;
    int max_leaf_size = index->settings->max_leaf_size;
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

    index_buffer_data *input_data=malloc(sizeof(index_buffer_data)*(num_threads-1));
    file_position_type *pos = malloc(sizeof(file_position_type));
    pthread_t threadid[num_threads-1];

    int j,conter=0;
    long long int i;
    int prev_flush_time=0,now_flush_time=0;
 

    
    int nodecounter=0;
    int sanity_counter=0;
    
    pthread_mutex_t lock_record=PTHREAD_MUTEX_INITIALIZER,lockfbl=PTHREAD_MUTEX_INITIALIZER,lock_index=PTHREAD_MUTEX_INITIALIZER,
                    lock_firstnode=PTHREAD_MUTEX_INITIALIZER,lock_disk=PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t *lockcbl;
    pthread_barrier_t lock_barrier1, lock_barrier2, lock_barrier3;
    pthread_barrier_init(&lock_barrier1, NULL, num_threads);
    pthread_barrier_init(&lock_barrier2, NULL, num_threads);
    pthread_barrier_init(&lock_barrier3, NULL, num_threads-1);

    lockcbl= malloc(sizeof(pthread_mutex_t)*LOCK_SIZE);

    double temp_traverse_total_time = 0;
    double temp_evaluate_total_time = 0;
    double temp_append_input_time = 0;
    double temp_append_total_time = 0;
    double temp_split_input_time = 0;
    double temp_split_total_time = 0;    


    
    for(i=0; i<LOCK_SIZE ;i++)
    {
        pthread_mutex_init(&lockcbl[i], NULL);
    }
    
    for (i = 0; i < (num_threads-1); i++)
    {
        input_data[i].index=index;
        input_data[i].lock_fbl=&lockfbl;
        input_data[i].lock_record=&lock_record;
        input_data[i].lock_cbl=lockcbl;
        input_data[i].lock_firstnode =&lock_firstnode;
        input_data[i].lock_index=&lock_index;
        input_data[i].blocid=i*READ_BLOCK_LENGTH;
        input_data[i].lock_disk=&lock_disk;
        //input_data[i].bufferpresize=bufferpresize;
        input_data[i].lock_barrier1=&lock_barrier1;
        input_data[i].lock_barrier2=&lock_barrier2;
        input_data[i].lock_barrier3=&lock_barrier3;
        input_data[i].finished=false;
        input_data[i].nodecounter=&nodecounter;
        input_data[i].sanity_counter=&sanity_counter;
        input_data[i].thread_traverse_total_time = 0;
        input_data[i].thread_evaluate_total_time = 0;
        input_data[i].thread_append_input_time = 0;
        input_data[i].thread_append_total_time = 0;
        input_data[i].thread_split_input_time = 0;
        input_data[i].thread_split_total_time = 0;
	input_data[i].split_node_data = calloc(max_leaf_size, sizeof(ts_type *));
	for (int j =0; j < max_leaf_size ; ++j)
	{
	  input_data[i].split_node_data[j] =  calloc(ts_length, sizeof(ts_type));
	}
    }

    ts     = malloc(sizeof(ts_type) * index->settings->timeseries_size*READ_BLOCK_LENGTH*(num_threads-1));
    ts_type * ts1    = malloc(sizeof(ts_type) * index->settings->timeseries_size*READ_BLOCK_LENGTH*(num_threads-1));
    ts_type * ts2;
    
        
    *pos = ftell(ifile);

    if(ts_num>READ_BLOCK_LENGTH*(num_threads-1))
      {
	//printf("first pos = %lu\n",(*pos)/1024);
	//printf("thread = %lu reading\n", pthread_self());	
        fread(ts1, sizeof(ts_type), index->settings->timeseries_size*READ_BLOCK_LENGTH*(num_threads-1), ifile);
	ts2   =ts;
        ts    = ts1;  
        ts1   =ts2;
	//memset(ts1, 0, sizeof(ts_type) * index->settings->timeseries_size*READ_BLOCK_LENGTH*(num_threads-1));
	  
        for ( j = 0; j < (num_threads-1); j++)
	  {   
            input_data[j].pos        = 0;
            input_data[j].ts         = ts;
            //input_data[j].saxv       =saxv;
            input_data[j].fin_number = READ_BLOCK_LENGTH*(num_threads-1);
	  }
	for (j = 0; j < (num_threads-1); j++)
	  {
            pthread_create(&(threadid[j]),NULL,indexbulkloadingworker,(void*)&(input_data[j]));
	  }
        for (i = READ_BLOCK_LENGTH*(num_threads-1)*2; i < ts_num; i+=READ_BLOCK_LENGTH*(num_threads-1))
	  {
            *pos = ftell(ifile);
	    //printf("next pos  = %lu\n",(*pos)/1024);
	    //printf("thread = %lu reading\n", pthread_self());
	    fread(ts1, sizeof(ts_type), index->settings->timeseries_size*READ_BLOCK_LENGTH*(num_threads-1), ifile);
            pthread_barrier_wait(&lock_barrier1);
	    sanity_counter = 0;
            //wait for the finish of other threads
            ts2   =ts;
            ts    =ts1;  
            ts1   =ts2;
            //__sync_fetch_and_add(&(index->buffer_manager->current_record_index),READ_BLOCK_LENGTH*(num_threads-1));	    
	
            for ( j = 0; j < (num_threads-1); j++)
            {   
	      input_data[j].pos        = 0;
	      input_data[j].ts         = ts;
	      //input_data[j].saxv       =saxv;
	      input_data[j].fin_number = READ_BLOCK_LENGTH*(num_threads-1);
            }
            //nodecounter=0;
            pthread_barrier_wait(&lock_barrier2);	    
	  }
	
        pthread_barrier_wait(&lock_barrier1);
        for (j = 0; j < (num_threads-1); j++)
	  {
            input_data[j].finished=true;
	  }
        //    __sync_fetch_and_add(&(index->buffer_manager->current_record_index),READ_BLOCK_LENGTH*(num_threads-1));	    	
        pthread_barrier_wait(&lock_barrier2);
        for (j = 0; j < (num_threads-1); j++)
        {
            pthread_join(threadid[j],NULL);
            #if DETAILED_STATS == 1
	    temp_traverse_total_time = fmax(temp_traverse_total_time, input_data[j].thread_traverse_total_time);
	    temp_evaluate_total_time = fmax(temp_evaluate_total_time, input_data[j].thread_evaluate_total_time);	
	    temp_append_total_time = fmax(temp_append_total_time, input_data[j].thread_append_total_time);
	    temp_append_input_time = fmax(temp_append_input_time, input_data[j].thread_append_input_time);
	    temp_split_total_time = fmax(temp_split_total_time, input_data[j].thread_split_total_time);
	    temp_split_input_time = fmax(temp_split_input_time, input_data[j].thread_split_input_time);		
            #endif	    
        }	
      }

    *pos = ftell(ifile);    
    //printf("next pos in rest of data= %lu\n",(*pos)/1024);
    COUNT_INPUT_TIME_START
    fread(ts1, sizeof(ts_type), index->settings->timeseries_size*(ts_num%(READ_BLOCK_LENGTH*(num_threads-1))), ifile);
    COUNT_INPUT_TIME_END
    sanity_counter = 0;
    ts2   =ts;
    ts    =ts1;
    ts1   =ts2;
    

    //handle the rest data
    int conter_ts_number=ts_num%(READ_BLOCK_LENGTH*(num_threads-1));
    //sax_save_number=conter_ts_number;
    for ( j = 0; j < (num_threads-1); j++)
    {   
      input_data[j].pos = 0;
      input_data[j].ts=ts;
      input_data[j].fin_number=fmin(conter_ts_number,READ_BLOCK_LENGTH * (num_threads-1));
      //input_data[j].saxv=saxv;
    }
    
    for ( j = 0; j < (num_threads-1); j++)
    { 
        input_data[j].finished=false;
        pthread_create(&(threadid[j]),NULL,indexbulkloadingworker,(void*)&(input_data[j]));
    }
    
    //printf("thread = %lu returned from bulkloading rest of the data\n", pthread_self());
    
    pthread_barrier_wait(&lock_barrier1);   
    for (j = 0; j < (num_threads-1); j++)
      {
	input_data[j].finished=true;
      }
    pthread_barrier_wait(&lock_barrier2);
    for (j = 0; j < (num_threads-1); j++)
      {
	pthread_join(threadid[j],NULL);
        #if DETAILED_STATS == 1
	temp_traverse_total_time = fmax(temp_traverse_total_time, input_data[j].thread_traverse_total_time);
	temp_evaluate_total_time = fmax(temp_evaluate_total_time, input_data[j].thread_evaluate_total_time);	
	temp_append_total_time = fmax(temp_append_total_time, input_data[j].thread_append_total_time);
	temp_append_input_time = fmax(temp_append_input_time, input_data[j].thread_append_input_time);
	temp_split_total_time = fmax(temp_split_total_time, input_data[j].thread_split_total_time);
	temp_split_input_time = fmax(temp_split_input_time, input_data[j].thread_split_input_time);		
        #endif
	for (int i = 0 ; i < max_leaf_size; ++i)
	  {
	    free(input_data[j].split_node_data[i]);
	  }
	
	free(input_data[j].split_node_data);

      }
     
    free(ts);
    free(ts1);
    free(input_data);
    free(lockcbl);
    free(pos);
    pthread_barrier_destroy(&lock_barrier1);
    pthread_barrier_destroy(&lock_barrier2);
    pthread_barrier_destroy(&lock_barrier3);

    COUNT_PARTIAL_INPUT_TIME_START
    if(fclose(ifile))
    {   
        fprintf(stderr, "Error in hercules_file_loaders.c: Could not close the filename %s", ifilename);
        return FAILURE;
    }
    COUNT_PARTIAL_INPUT_TIME_END    

    COUNT_PARTIAL_TIME_END
      
	      
    index->stats->idx_traverse_tree_total_time  = temp_traverse_total_time;
    index->stats->idx_evaluate_split_policies_total_time  = temp_evaluate_total_time;
    index->stats->idx_append_ts_to_leaf_total_time  = temp_append_total_time;
    index->stats->idx_append_ts_to_leaf_input_time  = temp_append_input_time;
    index->stats->idx_split_node_total_time  = temp_split_total_time;
    index->stats->idx_split_node_input_time  = temp_split_input_time;        
    


    index->stats->idx_building_total_time  += partial_time;	
    index->stats->idx_building_input_time  += partial_input_time;
    index->stats->idx_building_output_time += partial_output_time;
    index->stats->idx_building_seq_input_count  += partial_seq_input_count;
    index->stats->idx_building_seq_output_count += partial_seq_output_count;
    index->stats->idx_building_rand_input_count  += partial_rand_input_count;
    index->stats->idx_building_rand_output_count += partial_rand_output_count;	    
 
    //RESET_PARTIAL_COUNTERS()
    //COUNT_PARTIAL_TIME_START

      return SUCCESS;      

}

void* indexbulkloadingworker(void *transferdata)
{
    int fin_number=((index_buffer_data*)transferdata)->fin_number;
    file_position_type *pos = malloc(sizeof(file_position_type));
    //sax_type * saxv;
    //int offset_saxv=((index_buffer_data*)transferdata)->blocid*((index_buffer_data*)transferdata)->index->settings->paa_segments;
    int paa_segments=((index_buffer_data*)transferdata)->index->settings->paa_segments;
    struct hercules_index *index= ((index_buffer_data*)transferdata)->index;
    //int i=0,j;
    int i,j;
    pthread_barrier_t *lock_barrier1=((index_buffer_data*)transferdata)->lock_barrier1;
    pthread_barrier_t *lock_barrier2=((index_buffer_data*)transferdata)->lock_barrier2;
    pthread_barrier_t *lock_barrier3=((index_buffer_data*)transferdata)->lock_barrier3;
    int current_ts = -1;
    
    while(!((index_buffer_data*)transferdata)->finished)
    {   
      //saxv=(((index_buffer_data*)transferdata)->saxv);
      //        for (i=0;i<fin_number;i++)
      while (1)
	  {
	    
	    current_ts = __sync_fetch_and_add(((index_buffer_data*)transferdata)->sanity_counter,1);
	    if(current_ts >= fin_number){
	      //printf ("Break\n");
	      break;
            }	    
	    *pos= current_ts*index->settings->timeseries_size;
	    
	    //printf("thread = %lu inserting series %lu out of %d, ts[0] = %g\n", pthread_self(),*pos/256,fin_number, (((ts_type*)(((index_buffer_data*)transferdata)->ts)))[*pos] );
	    if (!hercules_index_insert_parallel(index, index->first_node,
					      (((ts_type*)(((index_buffer_data*)transferdata)->ts)+*pos)) , (void*)transferdata))	    
	      {
		fprintf(stderr, "Error in hercules_file_loaders.c:  Could not \
                           add the time series to the index.\n");
	      return FAILURE;              
	      }
        }
	    
      //printf("finished bulk loading\n");
	  // pthread_barrier_wait(lock_barrier3);

	pthread_barrier_wait(lock_barrier1);
        pthread_barrier_wait(lock_barrier2);
    }
    free(pos);
    //free(sax);
}

 enum response hercules_index_insert_parallel(struct hercules_index *index, struct hercules_node *node,  ts_type * timeseries, void *thread_data)
{

  struct hercules_node *subtree = NULL;
  boolean lock_subtree =false;
  boolean lock_node =false;
  
  //traverse the index tree to find the appropriate node
  //struct hercules_node * node = index->first_node;

  //Allocate memory for the series sketch
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


  //printf("thread = %lu locked_node = %p level = %d is_leaf = %d\n", pthread_self(),(void *)node, node->level,node->is_leaf);
  //pthread_mutex_lock(&(init_node->lock_data));

 
  while (!node->is_leaf )//& !node->is_splitting)
  {
    if (node_split_policy_route_to_left(node,timeseries,&timeseries_segment_sketch))
      node = node->left_child;
    else
      node = node->right_child;

   }

 #if DETAILED_STATS == 1
 gettimeofday(&end_time, NULL);
 tS = start_time.tv_sec*1000000 + (start_time.tv_usec);
 tE = end_time.tv_sec*1000000 + (end_time.tv_usec);
 ((index_buffer_data*)thread_data)->thread_traverse_total_time += (tE - tS);
 gettimeofday(&start_time, NULL); 
 #endif	      

 //if (node->parent == NULL)
 //   pthread_mutex_lock(&(node->lock_data));

  //printf("thread = %lu locking node = %p level = %d BEGIN\n", pthread_self(),(void *)node, node->level);    
  pthread_mutex_lock(&(node->lock_data));
   
  if (node->is_leaf)
  {
  
    //printf("thread = %lu processed_node = %p level = %d BEGIN\n", pthread_self(),(void *)node, node->level);    
    if(!update_node_statistics_parallel(node, timeseries,&timeseries_segment_sketch))
    {
        fprintf(stderr,"Error in hercules_index.c: could not update \
                        statistics at node %s\n", node->filename);
        return FAILURE;
    }

	//pthread_mutex_lock(&(index->buffer_manager->lock_file_buffer));
    if(!append_ts_to_node_parallel(index,node, timeseries))
    {
        fprintf(stderr,"Error in hercules_index.c: could not append \
                        time series to node %s\n", node->filename);
        return FAILURE;
    }
	//pthread_mutex_unlock(&(index->buffer_manager->lock_file_buffer));
    
    #if DETAILED_STATS == 1
    gettimeofday(&end_time, NULL);
    tS = start_time.tv_sec*1000000 + (start_time.tv_usec);
    tE = end_time.tv_sec*1000000 + (end_time.tv_usec);
    ((index_buffer_data*)thread_data)->thread_append_total_time += (tE - tS);
    gettimeofday(&start_time, NULL); 
    #endif

    //if split needed, split the node and refresh curr_node
    //if (node->parent == NULL)
    //  pthread_mutex_lock(&(index->buffer_manager->lock_file_map));

    if (node->node_size >= index->settings->max_leaf_size)
    {
      split_node_detailed(index,node, &timeseries_segment_sketch,thread_data);
      node->is_leaf = 0;
    } //end if_split_node
    //if (node->parent == NULL)
    //  pthread_mutex_unlock(&(index->buffer_manager->lock_file_map));
    
    #if DETAILED_STATS == 1
    gettimeofday(&end_time, NULL);
    tS = start_time.tv_sec*1000000 + (start_time.tv_usec);
    tE = end_time.tv_sec*1000000 + (end_time.tv_usec);
    ((index_buffer_data*)thread_data)->thread_split_total_time += (tE - tS);
    #endif
    //pthread_mutex_unlock(&(node->lock_data));
  } //end if_node_is_leaf
  else{
    //printf("No longer a leaf level = %d!\n", node->level);
 	//if (node->parent == NULL)
    //printf("thread = %lu unlocking node = %p level = %d BEGIN\n", pthread_self(),(void *)node, node->level);    
    //pthread_mutex_unlock(&(node->lock_data));    
    //pthread_mutex_destroy(&(node->lock_data));
    hercules_index_insert_parallel(index,node,timeseries, thread_data);
  }

 //if (node->parent == NULL)
    //printf("thread = %lu unlocking node = %p level = %d BEGIN\n", pthread_self(),(void *)node, node->level);    
   pthread_mutex_unlock(&(node->lock_data));
  free(timeseries_segment_sketch.indicators);    
  return SUCCESS;
}
*/

enum response hercules_index_binary_file_p(const char *ifilename, file_position_type ts_num,
					 struct hercules_index *index, int num_threads, int start_parallel)

{
    double parse_time = 0;
    int ts_length = index->settings->timeseries_size;
    int max_leaf_size = index->settings->max_leaf_size;
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

    index_buffer_data *input_data=malloc(sizeof(index_buffer_data)*(num_threads-1));
    file_position_type *pos = malloc(sizeof(file_position_type));
    pthread_t threadid[num_threads-1];

    int j,conter=0;
    long long int i;
    int prev_flush_time=0,now_flush_time=0;
 

    
    int nodecounter=0;
    int sanity_counter=0;
    
    pthread_mutex_t lock_record=PTHREAD_MUTEX_INITIALIZER,lockfbl=PTHREAD_MUTEX_INITIALIZER,lock_index=PTHREAD_MUTEX_INITIALIZER,
                    lock_firstnode=PTHREAD_MUTEX_INITIALIZER,lock_disk=PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t *lockcbl;
    pthread_barrier_t lock_barrier1, lock_barrier2, lock_barrier3, lock_barrier4; //barrier1 is the double buffer barrier and barrier2 is the flush barrier
    pthread_barrier_init(&lock_barrier1, NULL, num_threads); //dbBarrier
    pthread_barrier_init(&lock_barrier2, NULL, num_threads-1); //flushBarrier
    pthread_barrier_init(&lock_barrier3, NULL, num_threads);
    pthread_barrier_init(&lock_barrier4, NULL, num_threads-1); //bwBarrier (busy workers)   
    
    lockcbl= malloc(sizeof(pthread_mutex_t)*LOCK_SIZE);

    double temp_traverse_total_time = 0;
    double temp_evaluate_total_time = 0;
    double temp_append_input_time = 0;
    double temp_append_total_time = 0;
    double temp_split_input_time = 0;
    double temp_split_total_time = 0;
    double temp_flush_output_time = 0;
    double temp_flush_worker_time = 0;
    double temp_flush_coordinator_time = 0;
    double temp_reinsert_time = 0;
    double temp_sanity_counter_time = 0;
    double temp_insert_to_node_time = 0;
    double temp_bulk_loading_time = 0;
    

    int t_buffer_size = index->buffer_manager->max_record_index / (num_threads - 1);
    bool flush = false;
    bool handshake = false;
    volatile int cnt = 0;
    int flush_size = 0;
    volatile int flush_counter = 0;
    //bool flush_order = false;
    volatile int flush_order = 0;
    int finished_batch = 0;
    volatile bool barrier4_reached = false;    
    //int busy_insert_workers = num_threads - 1;
    	 
    for(i=0; i<LOCK_SIZE ;i++)
    {
        pthread_mutex_init(&lockcbl[i], NULL);
    }
    
    for (i = 0; i < (num_threads-1); i++)
    {
        input_data[i].index=index;
        input_data[i].lock_fbl=&lockfbl;
        input_data[i].lock_record=&lock_record;
        input_data[i].lock_cbl=lockcbl;
        input_data[i].lock_firstnode =&lock_firstnode;
        input_data[i].lock_index=&lock_index;
        input_data[i].blocid=i*READ_BLOCK_LENGTH;
        input_data[i].lock_disk=&lock_disk;
        //input_data[i].bufferpresize=bufferpresize;
        input_data[i].lock_barrier1=&lock_barrier1;
        input_data[i].lock_barrier2=&lock_barrier2;
        input_data[i].lock_barrier3=&lock_barrier3;
        input_data[i].lock_barrier4=&lock_barrier4;		
        input_data[i].finished=false;
        input_data[i].nodecounter=&nodecounter;
        input_data[i].sanity_counter=&sanity_counter;
        input_data[i].thread_traverse_total_time = 0;
        input_data[i].thread_evaluate_total_time = 0;
        input_data[i].thread_append_input_time = 0;
        input_data[i].thread_append_total_time = 0;
        input_data[i].thread_split_input_time = 0;
        input_data[i].thread_split_total_time = 0;
        input_data[i].thread_flush_output_time = 0;	
        input_data[i].thread_flush_worker_time = 0;	
        input_data[i].thread_flush_coordinator_time = 0;	
        input_data[i].thread_reinsert_time = 0;	
        input_data[i].thread_insert_to_node_time = 0;	
        input_data[i].thread_sanity_counter_time = 0;	
        input_data[i].thread_bulk_loading_time = 0;	
     	input_data[i].split_node_data = calloc(max_leaf_size, sizeof(ts_type *));
        input_data[i].buffer_max=t_buffer_size-max_leaf_size;
        input_data[i].buffer_offset= index->buffer_manager->current_record + (sizeof(ts_type) * ts_length) * t_buffer_size * i ;
        input_data[i].buffer_counter=0;
        input_data[i].flush_handshake=0;
        input_data[i].flush_order=&flush_order;
	input_data[i].is_flusher = false;
        input_data[i].flush_counter=&flush_counter;
	input_data[i].threads_data = input_data;
	input_data[i].num_insert_workers = num_threads-1;
	//input_data[i].busy_insert_workers = &busy_insert_workers;
	input_data[i].barrier4_reached = false;		
	input_data[i].is_busy = false;	
    	input_data[i].thread_id = i;
	input_data[i].is_full = false;
	
	for (int j =0; j < max_leaf_size ; ++j)
	{
	  input_data[i].split_node_data[j] =  calloc(ts_length, sizeof(ts_type));
	}
    }

    input_data[0].is_flusher = true;
    ts     = malloc(sizeof(ts_type) * index->settings->timeseries_size*READ_BLOCK_LENGTH*(num_threads-1));
    ts_type * ts1    = malloc(sizeof(ts_type) * index->settings->timeseries_size*READ_BLOCK_LENGTH*(num_threads-1));
    ts_type * ts2;

    *pos = ftell(ifile);

    if(ts_num>READ_BLOCK_LENGTH*(num_threads-1))
      {
	//printf("first pos = %lu\n",(*pos)/1024);
	//printf("thread = %lu reading\n", pthread_self());	
        fread(ts1, sizeof(ts_type), index->settings->timeseries_size*READ_BLOCK_LENGTH*(num_threads-1), ifile);
	ts2   =ts;
        ts    = ts1;  
        ts1   =ts2;
	//memset(ts1, 0, sizeof(ts_type) * index->settings->timeseries_size*READ_BLOCK_LENGTH*(num_threads-1));
	  
        for ( j = 0; j < (num_threads-1); j++)
	  {   
            input_data[j].pos        = 0;
            input_data[j].ts         = ts;
            //input_data[j].saxv       =saxv;
            input_data[j].fin_number = READ_BLOCK_LENGTH*(num_threads-1);
	  }
	for (j = 0; j < (num_threads-1); j++)
	  {
            pthread_create(&(threadid[j]),NULL,indexbulkloadingworker,(void*)&(input_data[j]));
	  }
        for (i = READ_BLOCK_LENGTH*(num_threads-1)*2; i <= ts_num; i+=READ_BLOCK_LENGTH*(num_threads-1))
	  {
            *pos = ftell(ifile);
	    //printf("next pos  = %lu\n",(*pos)/1024);
	    //printf("thread = %lu reading\n", pthread_self());
	    fread(ts1, sizeof(ts_type), index->settings->timeseries_size*READ_BLOCK_LENGTH*(num_threads-1), ifile);
	    int checked_batch = 0;
	    /*
            while (1)
	    {
	      pthread_mutex_lock(&lockfbl);
	      checked_batch = finished_batch;
	      pthread_mutex_unlock(&lockfbl);	      
	      if (checked_batch >= (num_threads - 1))
	      {
		break;	
	      }
	      if (__sync_fetch_and_add(&flush_counter,0) >= FLUSH_THRESHOLD)
	      {
   		    flush_coordinator(input_data);
	      }	   
	    }
	    pthread_mutex_lock(&lockfbl);
	    finished_batch = 0;
	    pthread_mutex_unlock(&lockfbl);	      	    

            while (__sync_fetch_and_add(&finished_batch,0) < (num_threads-1) )
	    {
	      if (__sync_fetch_and_add(&flush_counter,0) >= FLUSH_THRESHOLD)
	      {
   		    flush_coordinator(input_data);
	      }	   
	    }
	    finished_batch = 0;	    
	    */
            pthread_barrier_wait(&lock_barrier1);
            //busy_insert_workers = num_threads-1 ;
	    flush_counter = 0;
	    sanity_counter = 0;
	    pthread_barrier_wait(&lock_barrier1);
  	    //printf("DB coordinator unlocked both db barrier\n", pthread_self());	
	    //printf("DB coordinator ready for next batch \n", pthread_self());

            //wait for the finish of other threads
            ts2   =ts;
            ts    =ts1;  
            ts1   =ts2;
            //__sync_fetch_and_add(&(index->buffer_manager->current_record_index),READ_BLOCK_LENGTH*(num_threads-1));	    
            for ( j = 0; j < (num_threads-1); j++)
            {   
	      input_data[j].pos        = 0;
	      input_data[j].ts         = ts;
	      //input_data[j].saxv       =saxv;
	      input_data[j].fin_number = READ_BLOCK_LENGTH*(num_threads-1);
            }
            pthread_barrier_wait(&lock_barrier3);	    
	    //printf("DB coordinator unlocked barrier3 \n", pthread_self());	

	  }
        pthread_barrier_wait(&lock_barrier1);
	//busy_insert_workers = num_threads-1 ;
	flush_counter = 0;
	input_data[0].is_flusher = true;
	pthread_barrier_wait(&lock_barrier1);
	//printf("DB coordinator unlocked both db barrier to finish\n", pthread_self());	
	
        for (j = 0; j < (num_threads-1); j++)
	  {
            input_data[j].finished=true;
	  }
        pthread_barrier_wait(&lock_barrier3);
	//printf("DB coordinator unlocked barrier3 to finish\n", pthread_self());	

	    
        for (j = 0; j < (num_threads-1); j++)
        {
            pthread_join(threadid[j],NULL);
	       //printf("thread %lu joined\n", pthread_self());
            #if DETAILED_STATS == 1
	    temp_traverse_total_time = fmax(temp_traverse_total_time, input_data[j].thread_traverse_total_time);
	    temp_evaluate_total_time = fmax(temp_evaluate_total_time, input_data[j].thread_evaluate_total_time);	
	    temp_append_total_time = fmax(temp_append_total_time, input_data[j].thread_append_total_time);
	    temp_append_input_time = fmax(temp_append_input_time, input_data[j].thread_append_input_time);
	    //temp_append_input_time += input_data[j].thread_append_input_time;	
	    temp_split_total_time = fmax(temp_split_total_time, input_data[j].thread_split_total_time);
	    temp_split_input_time = fmax(temp_split_input_time, input_data[j].thread_split_input_time);
	    //temp_split_input_time += input_data[j].thread_split_input_time;
	    temp_flush_output_time = fmax(temp_flush_output_time, input_data[j].thread_flush_output_time);
	    temp_flush_worker_time = fmax(temp_flush_worker_time, input_data[j].thread_flush_worker_time);
	    temp_flush_coordinator_time = fmax(temp_flush_coordinator_time, input_data[j].thread_flush_coordinator_time);
	    temp_reinsert_time = fmax(temp_reinsert_time, input_data[j].thread_reinsert_time);
	    temp_insert_to_node_time = fmax(temp_insert_to_node_time, input_data[j].thread_insert_to_node_time);
	    temp_sanity_counter_time = fmax(temp_sanity_counter_time, input_data[j].thread_sanity_counter_time);
	    temp_bulk_loading_time = fmax(temp_bulk_loading_time, input_data[j].thread_bulk_loading_time);
		//temp_sanity_counter_time += input_data[j].thread_sanity_counter_time;
            #endif	    
        }	
      }

    
    *pos = ftell(ifile);    
    printf("next pos in rest of data= %lu\n",(*pos)/1024);
    COUNT_INPUT_TIME_START
    fread(ts1, sizeof(ts_type), index->settings->timeseries_size*(ts_num%(READ_BLOCK_LENGTH*(num_threads-1))), ifile);
    COUNT_INPUT_TIME_END
	    
    sanity_counter = 0;
    //busy_insert_workers = num_threads-1 ;
    flush_counter = 0;
    
    ts2   =ts;
    ts    =ts1;
    ts1   =ts2;
    

     //handle the rest data
    int conter_ts_number=ts_num%(READ_BLOCK_LENGTH*(num_threads-1));
    //sax_save_number=conter_ts_number;
    for ( j = 0; j < (num_threads-1); j++)
    {   
      input_data[j].pos = 0;
      input_data[j].ts=ts;
      input_data[j].fin_number=fmin(conter_ts_number,READ_BLOCK_LENGTH * (num_threads-1));
      /*
      if (input_data[j].buffer_counter >= input_data[j].buffer_max)
	++flush_counter;
      */
    }
    /*
    if (flush_counter >= FLUSH_THRESHOLD)
    {
      char * curr_time;
      curr_time= NULL;
      curr_time = malloc (sizeof(char) *26);
      get_current_time(curr_time);
      printf ("%s, batch remove ! %lu  \n", curr_time, flush_size);
      free(curr_time);
      
      flush_index_to_disk(index,index->first_node);
  
      memset(index->buffer_manager->mem_array,0,index->buffer_manager->max_record_index);
      index->buffer_manager->current_record_index = 0;
      index->buffer_manager->current_record = index->buffer_manager->mem_array;
     }	      
    */
    
    pthread_barrier_init(&lock_barrier1, NULL, num_threads);
    pthread_barrier_init(&lock_barrier2, NULL, num_threads-1);
    pthread_barrier_init(&lock_barrier3, NULL, num_threads);
    pthread_barrier_init(&lock_barrier4, NULL, num_threads-1);        
    for ( j = 0; j < (num_threads-1); j++)
    { 
        input_data[j].finished=false;
	pthread_create(&(threadid[j]),NULL,indexbulkloadingworker,(void*)&(input_data[j]));
    }
        
    //printf("thread = %lu returned from bulkloading rest of the data\n", pthread_self());

    //printf("DB coordinator unlocked db barrier reading rest of data\n", pthread_self());	    
    pthread_barrier_wait(&lock_barrier1);
    //busy_insert_workers = num_threads-1 ;
    pthread_barrier_wait(&lock_barrier1);
    
    for (j = 0; j < (num_threads-1); j++)
      {
	input_data[j].finished=true;
      }
    pthread_barrier_wait(&lock_barrier3);       
    for (j = 0; j < (num_threads-1); j++)
      {
	pthread_join(threadid[j],NULL);
        #if DETAILED_STATS == 1
	    temp_traverse_total_time = fmax(temp_traverse_total_time, input_data[j].thread_traverse_total_time);
	    temp_evaluate_total_time = fmax(temp_evaluate_total_time, input_data[j].thread_evaluate_total_time);	
	    temp_append_total_time = fmax(temp_append_total_time, input_data[j].thread_append_total_time);
	    temp_append_input_time = fmax(temp_append_input_time, input_data[j].thread_append_input_time);
	    //temp_append_input_time += input_data[j].thread_append_input_time;	
	    temp_split_total_time = fmax(temp_split_total_time, input_data[j].thread_split_total_time);
	    temp_split_input_time = fmax(temp_split_input_time, input_data[j].thread_split_input_time);
	    //temp_split_input_time += input_data[j].thread_split_input_time;
	    temp_flush_output_time = fmax(temp_flush_output_time, input_data[j].thread_flush_output_time);	
	    temp_flush_coordinator_time = fmax(temp_flush_coordinator_time, input_data[j].thread_flush_coordinator_time);
	    temp_flush_worker_time = fmax(temp_flush_worker_time, input_data[j].thread_flush_worker_time);
	    temp_insert_to_node_time = fmax(temp_insert_to_node_time, input_data[j].thread_insert_to_node_time);
	    temp_reinsert_time = fmax(temp_reinsert_time, input_data[j].thread_reinsert_time);
	    temp_sanity_counter_time = fmax(temp_sanity_counter_time, input_data[j].thread_sanity_counter_time);
	    temp_bulk_loading_time = fmax(temp_bulk_loading_time, input_data[j].thread_bulk_loading_time);
	    temp_split_input_time = fmax(temp_split_input_time, input_data[j].thread_split_input_time);
        #endif
	for (int i = 0 ; i < max_leaf_size; ++i)
	  {
	    free(input_data[j].split_node_data[i]);
	  }
	
	free(input_data[j].split_node_data);
 
    }

    printf ("finished\n");
    free(ts);
    free(ts1);
    free(input_data);
    free(lockcbl);
    free(pos);
    pthread_barrier_destroy(&lock_barrier1);
    pthread_barrier_destroy(&lock_barrier2);
    pthread_barrier_destroy(&lock_barrier3);
    pthread_barrier_destroy(&lock_barrier4);

    COUNT_PARTIAL_INPUT_TIME_START
    if(fclose(ifile))
    {   
        fprintf(stderr, "Error in hercules_file_loaders.c: Could not close the filename %s", ifilename);
        return FAILURE;
    }
    COUNT_PARTIAL_INPUT_TIME_END    

    COUNT_PARTIAL_TIME_END
      
	      
    index->stats->idx_traverse_tree_total_time  = temp_traverse_total_time;
    index->stats->idx_evaluate_split_policies_total_time  = temp_evaluate_total_time;
    index->stats->idx_append_ts_to_leaf_total_time  = temp_append_total_time;
    index->stats->idx_append_ts_to_leaf_input_time  = temp_append_input_time;
    index->stats->idx_split_node_total_time  = temp_split_total_time;
    index->stats->idx_split_node_input_time  = temp_split_input_time;    
    index->stats->idx_split_node_cpu_time  = temp_split_total_time - temp_split_input_time;
    index->stats->idx_flush_worker_time  = temp_flush_worker_time;    
    index->stats->idx_flush_coordinator_time  = temp_flush_coordinator_time;    
    index->stats->idx_insert_to_node_time  = temp_insert_to_node_time;    
    index->stats->idx_reinsert_time  = temp_reinsert_time;    
    index->stats->idx_sanity_counter_time  = temp_sanity_counter_time;    
    index->stats->idx_bulk_loading_time  = temp_bulk_loading_time;    

    /*
    index->stats->idx_building_total_time  += partial_time    
      +  index->stats->idx_traverse_tree_total_time 
      +  index->stats->idx_append_ts_to_leaf_total_time 
      +  index->stats->idx_evaluate_split_policies_total_time 
      +  index->stats->idx_split_node_total_time;
    */
    index->stats->idx_building_total_time  += partial_time;
    index->stats->idx_building_output_time = temp_flush_output_time + partial_output_time;
    
    index->stats->idx_building_input_time  += partial_input_time
      +  index->stats->idx_traverse_tree_input_time 
      +  index->stats->idx_append_ts_to_leaf_input_time 
      +  index->stats->idx_evaluate_split_policies_input_time 
      +  index->stats->idx_split_node_input_time 
      ;
    /*
    index->stats->idx_building_output_time += partial_output_time
      +  index->stats->idx_traverse_tree_output_time 
      +  index->stats->idx_append_ts_to_leaf_output_time 
      +  index->stats->idx_evaluate_split_policies_output_time 
      +  index->stats->idx_split_node_output_time 
      ;
    */
    
    index->stats->idx_building_cpu_time    += index->stats->idx_building_total_time 
      -  index->stats->idx_building_input_time 					 
      -  index->stats->idx_building_output_time 					 
      ;
    //index->stats->idx_building_total_time  += partial_time;	
    //index->stats->idx_building_input_time  += partial_input_time;
    //index->stats->idx_building_output_time += partial_output_time;
    //index->stats->idx_building_seq_input_count  += partial_seq_input_count;
    //index->stats->idx_building_seq_output_count += partial_seq_output_count;
    //index->stats->idx_building_rand_input_count  += partial_rand_input_count;
    //index->stats->idx_building_rand_output_count += partial_rand_output_count;	    
 
    //RESET_PARTIAL_COUNTERS()
    //COUNT_PARTIAL_TIME_START

      return SUCCESS;      

}
/*
void* indexbulkloadingworker(void *transferdata)
{
    int fin_number=((index_buffer_data*)transferdata)->fin_number;
    file_position_type *pos = malloc(sizeof(file_position_type));
    //sax_type * saxv;
    //int offset_saxv=((index_buffer_data*)transferdata)->blocid*((index_buffer_data*)transferdata)->index->settings->paa_segments;
    int paa_segments=((index_buffer_data*)transferdata)->index->settings->paa_segments;
    struct hercules_index *index= ((index_buffer_data*)transferdata)->index;
    //int i=0,j;
    int i,j;


    int current_ts = -1;

    
    while(!((index_buffer_data*)transferdata)->finished)
    {   
      //saxv=(((index_buffer_data*)transferdata)->saxv);
      //        for (i=0;i<fin_number;i++)
      pthread_barrier_t *lock_barrier1=((index_buffer_data*)transferdata)->lock_barrier1;
      pthread_barrier_t *lock_barrier3=((index_buffer_data*)transferdata)->lock_barrier3;      
      while (1)
	  {	    
	    current_ts = __sync_fetch_and_add(((index_buffer_data*)transferdata)->sanity_counter,1);
	    if(current_ts >= fin_number){
	      //printf("thread = %lu inserted %d series out of %d\n", pthread_self(), (((index_buffer_data*)transferdata)->buffer_counter) , fin_number);
	      //printf("thread = %d inserted %d series out of %d\n", ((index_buffer_data*)transferdata)->thread_id, (((index_buffer_data*)transferdata)->buffer_counter) , fin_number);
	      break;
            }	    
	    *pos= current_ts*index->settings->timeseries_size;
	    
	    //printf("thread = %lu inserting series #%d\n", pthread_self(), (((index_buffer_data*)transferdata)->buffer_counter) );
	    //printf("thread = %lu inserting series %lu out of %d, ts[0] = %g\n", pthread_self(),*pos/256,fin_number, (((ts_type*)(((index_buffer_data*)transferdata)->ts)))[*pos] );
	    if (!hercules_index_insert_parallel(index, index->first_node,
					      (((ts_type*)(((index_buffer_data*)transferdata)->ts)+*pos)) , (void*)transferdata))	    
	      {
		fprintf(stderr, "Error in hercules_file_loaders.c:  Could not \
                           add the time series to the index.\n");
	         return FAILURE;              
	      }
	    if ((((index_buffer_data*)transferdata)->buffer_counter) >= (((index_buffer_data*)transferdata)->buffer_max)  ||
		(__sync_fetch_and_add(((index_buffer_data*)transferdata)->flush_order,0)  != 0))
	      {
		//printf("thread = %lu ready to be flushed\n", pthread_self());
		printf("thread = %d ready to be flushed\n",((index_buffer_data*)transferdata)->thread_id );	     
		flush_worker(transferdata);
		pthread_barrier_wait(((index_buffer_data*)transferdata)->lock_barrier2);
		printf("thread = %d unlocked flush barrier \n", ((index_buffer_data*)transferdata)->thread_id );	     			
		((index_buffer_data*)transferdata)->flush_handshake = false;      
		((index_buffer_data*)transferdata)->buffer_counter = 0;		
	      }
	  }
      //pthread_mutex_lock(((index_buffer_data*)transferdata)->lock_fbl);
      //++(*(((index_buffer_data*)transferdata)->finished_batch));
      //pthread_mutex_unlock(((index_buffer_data*)transferdata)->lock_fbl);
      
           __sync_fetch_and_add(((index_buffer_data*)transferdata)->finished_batch,1);
      printf("thread = %d finished its batch\n",((index_buffer_data*)transferdata)->thread_id );	     	    
      pthread_barrier_wait(lock_barrier1);
      //printf("thread = %d finished its batch\n",((index_buffer_data*)transferdata)->thread_id );	     	          
      pthread_barrier_wait(lock_barrier3);
    }
    
    free(pos);
    //free(sax);
}
*/


void* indexbulkloadingworker(void *transferdata)
{
    int fin_number=((index_buffer_data*)transferdata)->fin_number;
    file_position_type *pos = malloc(sizeof(file_position_type));
    //sax_type * saxv;
    //int offset_saxv=((index_buffer_data*)transferdata)->blocid*((index_buffer_data*)transferdata)->index->settings->paa_segments;
    int paa_segments=((index_buffer_data*)transferdata)->index->settings->paa_segments;
    struct hercules_index *index= ((index_buffer_data*)transferdata)->index;
    //int i=0,j;
    int i,j,k;
    int num_threads = ((index_buffer_data*)transferdata)->num_insert_workers;
    
    int volatile cnt = 0;
    
    int current_ts = -1;

    int buffer_counter =  ((index_buffer_data*)transferdata)->buffer_counter;
    
     double tS = 0;
     double tE = 0;
     struct timeval start_time;
     struct timeval end_time;

	    #if DETAILED_STATS == 1
         gettimeofday(&start_time, NULL);
	    #endif


    while(!((index_buffer_data*)transferdata)->finished)
    {   

      function_while_1(transferdata, pos,fin_number,current_ts,index);   

      function_after_while_1(transferdata, index);

  	  ((index_buffer_data*)transferdata)->barrier4_reached = false;
  	   // printf("thread = %d resetting barrier4_reached to FALSE\n",((index_buffer_data*)transferdata)->thread_id );	 
    	    	 
	  pthread_barrier_wait(((index_buffer_data*)transferdata)->lock_barrier1);
	  pthread_barrier_wait(((index_buffer_data*)transferdata)->lock_barrier1);	 
	 //printf("thread = %d unlocking both db barrier\n",((index_buffer_data*)transferdata)->thread_id );	     	    

	  pthread_barrier_wait(((index_buffer_data*)transferdata)->lock_barrier3);
	 //printf("thread = %d unlocking barrier3\n",((index_buffer_data*)transferdata)->thread_id );	     	    

    }

	    #if DETAILED_STATS == 1
         gettimeofday(&end_time, NULL);
    	tS = start_time.tv_sec*1000000 + (start_time.tv_usec);
	    tE = end_time.tv_sec*1000000 + (end_time.tv_usec);
	    ((index_buffer_data*)transferdata)->thread_bulk_loading_time += (tE - tS);    
	    #endif
    
    free(pos);
    //free(sax);
}

void function_while_1(void *transferdata, file_position_type * pos, int fin_number, int current_ts, struct hercules_index * index)  
{
       double tS = 0;
     double tE = 0;
     struct timeval start_time;
     struct timeval end_time;


  while (1)
	  {	    
	    if (!(((index_buffer_data*)transferdata)->is_flusher))
	    {
	      if ((((index_buffer_data*)transferdata)->buffer_counter) >= (((index_buffer_data*)transferdata)->buffer_max))		
	      {
			((index_buffer_data*)transferdata)->is_full = true;	      
		    //__sync_fetch_and_add(((index_buffer_data*)transferdata)->flush_counter,1);
		    *(((index_buffer_data*)transferdata)->flush_counter) = *(((index_buffer_data*)transferdata)->flush_counter) + 1;
		    //printf("thread = %d is full and waiting to be flushed, busy_workers = %d\n",((index_buffer_data*)transferdata)->thread_id, *   	(((index_buffer_data*)transferdata)->busy_insert_workers));
   			break;
	      }
	      //if (((__sync_fetch_and_add(((index_buffer_data*)transferdata)->flush_order,0))))
		  if ( *(((index_buffer_data*)transferdata)->flush_order) )	      
          {
		     //printf("thread = %lu ready to be flushed\n", pthread_self());
		     //printf("thread = %d received order to be flushed in while 1\n",((index_buffer_data*)transferdata)->thread_id );
		  
   		     ((index_buffer_data*)transferdata)->barrier4_reached = true;
		     pthread_barrier_wait(((index_buffer_data*)transferdata)->lock_barrier4);
  		     //printf("thread = %d unlocked busy workers barrier in flush_worker\n", ((index_buffer_data*)transferdata)->thread_id );	     			      
 		     flush_worker(transferdata);
		      //printf("thread = %d returned from flushing in while 1\n",((index_buffer_data*)transferdata)->thread_id );		  
	       }
	    }
	    else
	    {
	      //printf("flusher is thread = %d\n",((index_buffer_data*)transferdata)->thread_id );	      
	      //if ( (((__sync_fetch_and_add(((index_buffer_data*)transferdata)->flush_counter,0))) >= index->settings->flush_threshold) ||
  	      if ( *(((index_buffer_data*)transferdata)->flush_counter) >=  index->settings->flush_threshold	       ||
		   ((((index_buffer_data*)transferdata)->buffer_counter) >= (((index_buffer_data*)transferdata)->buffer_max))
		   )
		  {
		    // printf("thread = %d will issue flush order\n",((index_buffer_data*)transferdata)->thread_id );
		     *(((index_buffer_data*)transferdata)->flush_order) = 1;
  		     ((index_buffer_data*)transferdata)->barrier4_reached = true;      
		     pthread_barrier_wait(((index_buffer_data*)transferdata)->lock_barrier4);
		     //printf("thread = %d unlocked busy workers barrier in flush_coordinator\n", ((index_buffer_data*)transferdata)->thread_id );	    			        
		     flush_coordinator(transferdata);
		     //printf("thread = %d (flusher) returned from flushing in while 1\n",((index_buffer_data*)transferdata)->thread_id );		  		  
		  }
	    }

	    #if DETAILED_STATS == 1
         gettimeofday(&start_time, NULL);
	    #endif

   	     current_ts = __sync_fetch_and_add(((index_buffer_data*)transferdata)->sanity_counter,1);

	    #if DETAILED_STATS == 1
         gettimeofday(&end_time, NULL);
    	tS = start_time.tv_sec*1000000 + (start_time.tv_usec);
	    tE = end_time.tv_sec*1000000 + (end_time.tv_usec);
	    ((index_buffer_data*)transferdata)->thread_sanity_counter_time += (tE - tS);    
	    #endif

	    if(current_ts >= fin_number){
			  //printf("thread = %lu inserted %d series out of %d\n", pthread_self(), (((index_buffer_data*)transferdata)->buffer_counter) , fin_number);
			  //printf("thread = %d inserted %d series out of %d\n", ((index_buffer_data*)transferdata)->thread_id, (((index_buffer_data*)transferdata)->buffer_counter) , fin_number);
			  /*
			  if (!((index_buffer_data*)transferdata)->flush_halted_insert)
			{
			  ((index_buffer_data*)transferdata)->flush_handshake = true;
			  pthread_barrier_wait(((index_buffer_data*)transferdata)->lock_barrier2 );
			  printf("thread = %d unlocked barrier 2 in break\n", ((index_buffer_data*)transferdata)->thread_id );	     
			  }
			  */
			//printf ("Break\n");
   	        break;
         }	    
	    *pos= current_ts*index->settings->timeseries_size;
	    #if DETAILED_STATS == 1
         gettimeofday(&start_time, NULL);
	    #endif

	    if (!hercules_index_insert_parallel(index, index->first_node,
					      (((ts_type*)(((index_buffer_data*)transferdata)->ts)+*pos)) , (void*)transferdata))	    
	    {
		fprintf(stderr, "Error in hercules_file_loaders.c:  Could not \
                           add the time series to the index.\n");
	         return FAILURE;              
	    }
	    #if DETAILED_STATS == 1
         gettimeofday(&end_time, NULL);
    	tS = start_time.tv_sec*1000000 + (start_time.tv_usec);
	    tE = end_time.tv_sec*1000000 + (end_time.tv_usec);
	    ((index_buffer_data*)transferdata)->thread_insert_to_node_time += (tE - tS);    
	    #endif

	    //printf("thread = %lu inserting series #%d\n", pthread_self(), (((index_buffer_data*)transferdata)->buffer_counter) );
	    //printf("thread = %lu inserting series %lu out of %d, ts[0] = %g\n", pthread_self(),*pos/256,fin_number, (((ts_type*)(((index_buffer_data*)transferdata)->ts)))[*pos] );
	  }

}

/*
void function_after_while_1(void *transferdata, struct hercules_index * index)
{
  //while((__sync_fetch_and_add(((index_buffer_data*)transferdata)->busy_insert_workers,0)) > 0)
    //while(*(((index_buffer_data*)transferdata)->busy_insert_workers) > 0)
  int num_threads = ((index_buffer_data*)transferdata)->num_insert_workers;  
  while(are_workers_busy(transferdata, num_threads))
       {
	 //printf("thread = %d inside while(busy), busy_insert_workers = %d \n",((index_buffer_data*)transferdata)->thread_id,  *(((index_buffer_data*)transferdata)->busy_insert_workers));
	 if ( ((index_buffer_data*)transferdata)->is_flusher)
	   {
	     function_after_while_1_flusher(transferdata, index);
	     //if ( (((__sync_fetch_and_add(((index_buffer_data*)transferdata)->flush_counter,0))) >= index->settings->flush_threshold))
	   }
	   else
	   {
	     function_after_while_1_non_flusher(transferdata, index);	     
	      //if (((__sync_fetch_and_add(((index_buffer_data*)transferdata)->flush_order,0))))
	   }
       }
}
*/

void function_after_while_1(void *transferdata, struct hercules_index * index)
{
  if (!(((index_buffer_data*)transferdata)->barrier4_reached))
  {
    pthread_barrier_wait(((index_buffer_data*)transferdata)->lock_barrier4);
    //printf("thread = %d unlocked busy workers barrier after while 1\n", ((index_buffer_data*)transferdata)->thread_id );	     			  
  }
  else
  {
    //printf("thread = %d has already unlocked busy workers barrier after while 1\n", ((index_buffer_data*)transferdata)->thread_id );	   			  
  }
  
  if ( ((index_buffer_data*)transferdata)->is_flusher)
    {
      function_after_while_1_flusher(transferdata, index);
    }
  else
    {
      function_after_while_1_non_flusher(transferdata, index);
    }
}

 /*
void function_after_while_1(void *transferdata, struct hercules_index * index)
{
  //while((__sync_fetch_and_add(((index_buffer_data*)transferdata)->busy_insert_workers,0)) > 0)
    //while(*(((index_buffer_data*)transferdata)->busy_insert_workers) > 0)
  int num_threads = ((index_buffer_data*)transferdata)->num_insert_workers;
  //pthread_barrier_wait(((index_buffer_data*)transferdata)->lock_barrier4);
  while(are_workers_busy(transferdata, num_threads))
       {
	 //printf("thread = %d inside while(busy), busy_insert_workers = %d \n",((index_buffer_data*)transferdata)->thread_id,  *(((index_buffer_data*)transferdata)->busy_insert_workers));
	 if ( ((index_buffer_data*)transferdata)->is_flusher)
	   {
	     function_after_while_1_flusher(transferdata, index);
	     //if ( (((__sync_fetch_and_add(((index_buffer_data*)transferdata)->flush_counter,0))) >= index->settings->flush_threshold))
	   }
	   else
	   {
	     function_after_while_1_non_flusher(transferdata, index);	     
	      //if (((__sync_fetch_and_add(((index_buffer_data*)transferdata)->flush_order,0))))
	   }
       }
}
 */



void function_after_while_1_flusher(void *transferdata, struct hercules_index * index)
{
  if ( *(((index_buffer_data*)transferdata)->flush_counter) >= index->settings->flush_threshold)
	       flush_coordinator(transferdata);
}

void function_after_while_1_non_flusher(void *transferdata, struct hercules_index * index)
{
  if (*(((index_buffer_data*)transferdata)->flush_order) )
  {
    //printf("thread = %lu ready to be flushed\n", pthread_self());
    printf("thread = %d received order to be flushed after while 1\n",((index_buffer_data*)transferdata)->thread_id );	     
    flush_worker(transferdata);
    printf("thread = %d returned from flushing after while 1\n",((index_buffer_data*)transferdata)->thread_id );		        
	  
  }
}
void flush_worker(void *transferdata)
{

  int s;
       double tS = 0;
     double tE = 0;
     struct timeval start_time;
     struct timeval end_time;

  #if DETAILED_STATS == 1 
  gettimeofday(&start_time, NULL);
  #endif

  //printf("thread = %lu unlocked barrier 2\n", pthread_self());
  //printf("thread = %d unlocking barrier 2\n", ((index_buffer_data*)transferdata)->thread_id );	    

 __sync_fetch_and_add(&(((index_buffer_data*)transferdata)->flush_handshake),1);
 /*
  if (!(((index_buffer_data*)transferdata)->barrier4_reached))
  {
    ((index_buffer_data*)transferdata)->barrier4_reached = true;
    pthread_barrier_wait(((index_buffer_data*)transferdata)->lock_barrier4);
    printf("thread = %d unlocked busy workers barrier in flush_worker\n", ((index_buffer_data*)transferdata)->thread_id );	     			      
  }
 */
  //printf("thread = %d unlocking flush barrier \n", ((index_buffer_data*)transferdata)->thread_id );	     			
  s = pthread_barrier_wait(((index_buffer_data*)transferdata)->lock_barrier2);

  //printf("thread = %d unlocked flush barrier \n", ((index_buffer_data*)transferdata)->thread_id );	     			
   #if DETAILED_STATS == 1
         gettimeofday(&end_time, NULL);
    	tS = start_time.tv_sec*1000000 + (start_time.tv_usec);
	    tE = end_time.tv_sec*1000000 + (end_time.tv_usec);
	    ((index_buffer_data*)transferdata)->thread_flush_worker_time += (tE - tS);    
   #endif

}


void flush_coordinator(void *transferdata)
{
  int flush_size = 0;
  struct hercules_index * index = ((index_buffer_data*)transferdata)->index;
  //index_buffer_data* input_data = ((index_buffer_data*)transferdata)->threads_data;
  int num_threads = ((index_buffer_data*)transferdata)->num_insert_workers;
  int j = 0;
  volatile int cnt = 0;
  int s;

  double tS = 0;
  double tE = 0;
  struct timeval start_time;
  struct timeval start_time2;

  struct timeval end_time;
  struct timeval end_time2;
  int worker_id = ((index_buffer_data*)transferdata)->thread_id;
  #if DETAILED_STATS == 1 
  gettimeofday(&start_time2, NULL);
  #endif

    
  // __sync_fetch_and_add(((index_buffer_data*)transferdata)->flush_order,1);
  //  *(((index_buffer_data*)transferdata)->flush_order) = 1;
  //printf("thread = %d set flush order bit in flush_coordinator\n", ((index_buffer_data*)transferdata)->thread_id );	     			        
   /*
  if (!(((index_buffer_data*)transferdata)->barrier4_reached))
    {
      ((index_buffer_data*)transferdata)->barrier4_reached = true;      
      pthread_barrier_wait(((index_buffer_data*)transferdata)->lock_barrier4);
      printf("thread = %d unlocked busy workers barrier in flush_coordinator\n", ((index_buffer_data*)transferdata)->thread_id );	     			        
    }
   */
   
  __sync_fetch_and_add(&(((index_buffer_data*)transferdata)->threads_data[worker_id].flush_handshake),1);
  
  for ( j = 0; j < num_threads; j++) 
    {   
      flush_size += ((index_buffer_data*)transferdata)->threads_data[j].buffer_counter;
      printf("Flusher sent order to thread %d\n", j);	
      //check that all workers have received the order
      while (  __sync_fetch_and_add(&(((index_buffer_data*)transferdata)->threads_data[j].flush_handshake),0) != 1) 
	{
	  //for (int tmp = 0; tmp < BUSYWAIT_THRESHOLD; tmp++) cnt++;
	  //for (int tmp = 0; tmp < 1000; tmp++) cnt++;
	  //printf("thread = %lu  (coordinator) busy waiting\n", pthread_self());
	}
      printf("Flusher received handshake from thread %d\n", j);		  
    }

  printf("Flusher checked all handshakes\n");
  
  char * curr_time;
  curr_time= NULL;
  curr_time = malloc (sizeof(char) *26);
  get_current_time(curr_time);
  printf ("%s, batch remove ! %lu  \n", curr_time, flush_size);
  free(curr_time);

  #if DETAILED_STATS == 1 
  gettimeofday(&start_time, NULL);
  #endif

  flush_index_to_disk(index,index->first_node);
 
  #if DETAILED_STATS == 1
  gettimeofday(&end_time, NULL);
  tS = start_time.tv_sec*1000000 + (start_time.tv_usec);
  tE = end_time.tv_sec*1000000 + (end_time.tv_usec);
  ((index_buffer_data*)transferdata)->thread_flush_output_time += (tE - tS);
  #endif
  
  memset(index->buffer_manager->mem_array,0,index->buffer_manager->max_record_index);
  index->buffer_manager->current_record_index = 0;
  index->buffer_manager->current_record = index->buffer_manager->mem_array;
  
  //__sync_fetch_and_sub(((index_buffer_data*)transferdata)->flush_order,1);
   *(((index_buffer_data*)transferdata)->flush_order) = 0;
   printf("Flusher reset order to FALSE \n");
      
  for ( j = 0; j < num_threads; j++)
    {
      __sync_fetch_and_sub(&(((index_buffer_data*)transferdata)->threads_data[j].flush_handshake),1);            	
      ((index_buffer_data*)transferdata)->threads_data[j].buffer_counter = 0;
      ((index_buffer_data*)transferdata)->is_full = false;	            
    }
  
  *(((index_buffer_data*)transferdata)->flush_counter) = 0;
  
  //printf("thread = %d unlocking barrier 2\n", ((index_buffer_data*)transferdata)->thread_id );
  //printf("thread = %lu (coordinator) unlocking barrier 2\n", pthread_self());
  //pthread_barrier_wait(&lock_barrier2);
  //printf("thread = %d (flusher) unlocking flush barrier \n", ((index_buffer_data*)transferdata)->thread_id );	       
  s =   pthread_barrier_wait(((index_buffer_data*)transferdata)->lock_barrier2);
      // s =   pthread_barrier_wait(lock_barrier4);
  //printf("thread = %d (flusher) unlocked flush barrier \n", ((index_buffer_data*)transferdata)->thread_id );	     
  //	printf("thread = %lu (coordinator) unlocked flush  barrier \n", pthread_self());	
  #if DETAILED_STATS == 1
  gettimeofday(&end_time2, NULL);
  tS = start_time2.tv_sec*1000000 + (start_time2.tv_usec);
  tE = end_time2.tv_sec*1000000 + (end_time2.tv_usec);
  ((index_buffer_data*)transferdata)->thread_flush_coordinator_time += (tE - tS);
  #endif

}



enum response hercules_index_insert_parallel(struct hercules_index *index, struct hercules_node *node,  ts_type * timeseries, void *thread_data)
{

  struct hercules_node *subtree = NULL;
  boolean lock_subtree =false;
  boolean lock_node =false;
  
  //traverse the index tree to find the appropriate node
  //struct hercules_node * node = index->first_node;

  //Allocate memory for the series sketch
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


  //printf("thread = %lu locking node = %p level = %d BEGIN\n", pthread_self(),(void *)node, node->level);    
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

      //printf("thread = %lu processed_node = %p level = %d BEGIN\n", pthread_self(),(void *)node, node->level);    
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

   //printf("thread = %lu unlocking node = %p level = %d BEGIN\n", pthread_self(),(void *)node, node->level);    

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
	   //PRINT_STATS(distance)
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


