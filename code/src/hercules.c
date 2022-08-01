//
//  main.c
//  Created by Karima Echihabi on 18/12/2018
//


#include "../config.h"
#include "../globals.h"
#include "../include/systemutils.h"


#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include "../include/ts.h"
#include "../include/hercules_file_loaders.h"
#include "../include/hercules_index.h"
#include "../include/hercules_node.h"
#include "../include/hercules_file_buffer.h"
#include "../include/hercules_file_buffer_manager.h"
#include "../include/hercules_query_engine.h"

#ifdef VALUES
#include <values.h>
#endif
		
int main (int argc, char **argv)
{
    INIT_STATS()
    COUNT_TOTAL_TIME_START

    static char * dataset = "/home/karima/myDisk/data/Cgenerator/data_current.txt";
    static char * queries = "/home/karima/myDisk/data/Cgenerator/query_current.txt";
    static char * dataset_hists = "/home/karima/myDisk/data/Cgenerator/data_current_hists.txt";    
    
    static char * index_path = "out/";
    static unsigned int dataset_size = 1000;
    static unsigned int queries_size = 5;
    static unsigned int queries_skip = 0;    
    static unsigned int time_series_size = 256;
    static unsigned int init_segments = 1;  
    static unsigned int leaf_size = 100;
    static double buffered_memory_size = 6439.2; 
    static int  use_ascii_input = 0;
    static int mode = 0;
    static float minimum_distance = FLT_MAX;
    static float epsilon = 0;  //by default perform exact search 
    static float delta = 1;  //by default perform exact search, delta = 0 means ng
    static unsigned int k = 1;
    static unsigned int nprobes = 0;       
    struct hercules_index * index = NULL;
    boolean is_index_new = 1;
    boolean track_bsf = 0;
    boolean track_pruning = 0;
    boolean all_mindists = 0;
    boolean max_policy = 0;
    unsigned char incremental = 0;
    unsigned char in_memory = 0;
    int flush_threshold = 1;    
    unsigned int initial_db_size = 20000 * 24;    

    int num_write_threads=1;
    
    query_settings q_settings;    
    q_settings.serial = 0;
    q_settings.sims = 0;
    q_settings.psq = 0;
    q_settings.pmq = 0;
    q_settings.num_queues = 1;
    q_settings.num_query_threads = 1;    
    q_settings.num_read_threads = 1;
    q_settings.approx_stop_condition = 80; 
    q_settings.exact_stop_condition = 1; 
    q_settings.eapca_threshold = 0.25; 
    q_settings.sax_threshold = 0.9; 

    /*Adding support for SIMS on leaves*/

    static int paa_segments = 16;
    static int sax_cardinality = 8;

    while (1)
    {    
        static struct option long_options[] =  {
            {"ascii-input", required_argument, 0, 'a'},
            {"buffer-size", required_argument, 0, 'b'},	    
	    {"epsilon", required_argument, 0, 'c'},
            {"dataset", required_argument, 0, 'd'},
	    {"delta", required_argument, 0, 'e'},
            {"queries-size", required_argument, 0, 'f'},
            {"track-bsf", required_argument, 0, 'g'},
            {"incremental", no_argument, 0, 'h'},
            {"track-pruning", required_argument, 0, 'i'},
            {"all-mindists", required_argument, 0, 'j'},
            {"k", required_argument, 0, 'k'},
            {"leaf-size", required_argument, 0, 'l'},	    
            {"max-policy", required_argument, 0, 'm'},
            {"dataset-hists", required_argument, 0, 'n'},
            {"nprobes", required_argument, 0, 'o'},
            {"index-path", required_argument, 0, 'p'},
            {"queries", required_argument, 0, 'q'},
            {"serial",required_argument, 0, 'r'},
	    {"minimum-distance", required_argument, 0, 's'},
            {"timeseries-size", required_argument, 0, 't'},
            {"psq", no_argument, 0, 'u'},
            {"pmq", no_argument, 0, 'v'},
            {"sax-cardinality", required_argument, 0, 'w'},
            {"paa-segments", required_argument, 0, 'x'},	    
            {"num-queues", required_argument, 0, 'A'},
            {"num-query-threads", required_argument, 0, 'B'},
            {"num-read-threads", required_argument, 0, 'C'},	    	    
            {"approx-stop-condition", required_argument, 0, 'D'},	    	    
            {"exact-stop-condition", required_argument, 0, 'E'},	    	    
            {"in-memory", required_argument, 0, 'M'},	    	    
            {"mode", required_argument, 0, 'y'},
            {"dataset-size", required_argument, 0, 'z'},
            {"sims", no_argument, 0, 'S'},
            {"queries-skip", required_argument, 0, 'Q'},
            {"flush-threshold", required_argument, 0, 'F'},	    
            {"num-write-threads", required_argument, 0, 'G'},	    	    	    
            {"eapca-threshold", required_argument, 0, 'H'},	    
            {"sax-threshold", required_argument, 0, 'I'},	
            {"initial-db-size", required_argument, 0, 'J'},	    	    	    
            {"help", no_argument, 0, '?'}
        };
        
        /* getopt_long stores the option index here. */
        int option_index = 0;
        
        int c = getopt_long (argc, argv, "",
                             long_options, &option_index);
        if (c == -1)
            break;
        switch (c)
        {
            case 'q':
                queries = optarg;
                break;
						
            case 's':
                 minimum_distance = atof(optarg);
                 break;

            case 'c':
                 epsilon = atof(optarg);
                 break;

           case 'e':
                 delta = atof(optarg);
                 break;

	    case 'b':
                 buffered_memory_size = atof(optarg);
                 break;
		 
            case 'f':
                queries_size = atoi(optarg);
                if (queries_size < 1)
		{
		  fprintf(stderr,"Please change the queries size to be greater than 0.\n");
		  exit(-1);
		}
                break;

	    case 'Q':
                queries_skip = atoi(optarg);
                break;
		
            case 'k':
                k = atoi(optarg);
                if (k < 1)
		{
		  fprintf(stderr,"Please change the k to be greater than 0.\n");
		  exit(-1);
		}
                break;		

            case 'g':
  	      track_bsf = atoi(optarg);
	      break;

            case 'i':
  	      track_pruning = atoi(optarg);
	      break;

            case 'j':
  	      all_mindists = atoi(optarg);
	      break;	      	      

            case 'm':
  	      max_policy = atoi(optarg);
	      break;	      	      
	      
            case 'd':
                dataset = optarg;
                break;	       

  	    case 'n':
	        dataset_hists = optarg;
	        break;	       

            case 'p':
                index_path = optarg;
                break;
		
            case 'y':
                mode = atoi(optarg);
                break;
		
            case 'z':
                dataset_size = atoi(optarg);
                if (dataset_size < 1)
		{
		  fprintf(stderr,"Please change the dataset size to be greater than 0.\n");
		  exit(-1);
		}
                break;
		
            case 't':
                time_series_size = atoi(optarg);
                break;

            case 'o':
                nprobes = atoi(optarg);
                break;		

            case 'r':
                q_settings.serial = atoi(optarg);
                break;

            case 'A':
                q_settings.num_queues = atoi(optarg);
                break;

            case 'B':
                q_settings.num_query_threads = atoi(optarg);
                break;

            case 'C':
                q_settings.num_read_threads = atoi(optarg);
                break;

            case 'D':
                q_settings.approx_stop_condition = atof(optarg);
                break;

            case 'E':
                q_settings.exact_stop_condition = atof(optarg);
                break;
		
            case 'u':
                q_settings.psq = atoi(optarg);
                break;
            case 'v':
                q_settings.pmq = atoi(optarg);
                break;				
		
            case 'H':
                q_settings.eapca_threshold = atof(optarg);
                break;
            case 'I':
                q_settings.sax_threshold = atof(optarg);
                break;				
            case 'J':
                initial_db_size = atoi(optarg);
                break;				
            case 'l':
                leaf_size = atoi(optarg);
                if (leaf_size <= 1)
		{
		  fprintf(stderr,"Please change the leaf size to be greater than 1.\n");
		  exit(-1);
		}
                break;
		
            case 'w':
                sax_cardinality = atoi(optarg);
                break;

            case 'x':
                paa_segments = atoi(optarg);
                break;
            case 'F':
                flush_threshold = atoi(optarg);
                break;		
            case 'G':
                num_write_threads = atoi(optarg);
                break;		
		
            case '?':
                printf("Usage:\n\
                       \t--Queries and Datasets should be single precision\n\
                       \t--floating points. They can be binary of ascii.\n\
                       \t--However performance is faster with binary files\n\
                       \t--dataset XX \t\t\tThe path to the dataset file\n\
                       \t--queries XX \t\t\tThe path to the queries file\n\
                       \t--dataset-size XX \t\tThe number of time series to load\n\
                       \t--queries-size XX \t\tThe number of queries to run\n\
                       \t--mode: 0=index, 1=query, 2=index & query  3=calc_tlb\t\t\n\
                       \t--index-path XX \t\tThe path of the output folder\n\
                       \t--buffer-size XX \t\tThe size of the buffer memory in MB\n\
                       \t--timeseries-size XX\t\tThe size of each time series\n\
                       \t--ascii-input X \t\t\0 for ascii files and 1 for binary files\n\
                       \t--leaf-size XX\t\t\tThe maximum size of each leaf\n\
                       \t--help\n\n\
                       \t--**********************EXAMPLES**********************\n\n\
                       \t--*********************INDEX MODE*********************\n\n\
                       \t--bin/hercules --dataset XX --dataset-size XX             \n\n\
                       \t--          --index-path XX --timeseries-size XX --mode 0\n\n\
                       \t--*********************QUERY MODE*********************\n\n\
                       \t--bin/hercules --queries XX --queries-size XX             \n\n\
                       \t--           --index-path XX --mode 1                 \n\n\
                       \t--*****************INDEX AND QUERY MODE***************\n\n\
                       \t--bin/hercules --dataset XX --dataset-size XX             \n\n\
                       \t--          --timeseries-size XX --index-path XX      \n\n\
                       \t--           --queries XX --queries-size XX --mode 2  \n\n\
                       \t--****************************************************\n\n");
  
                return 0;
                break;
            case 'a':
  	        use_ascii_input = atoi(optarg);
                break;
            case 'h':
	      incremental = atoi(optarg);
                break;
            case 'S':
  	        q_settings.sims = 1;
                break;
   	    default:
                exit(-1);
                break;
        }
    }

    minimum_distance = FLT_MAX;
    if (mode == 0)  //only build and store the index
    {

        RESET_PARTIAL_COUNTERS()
        COUNT_PARTIAL_TIME_START
	  struct hercules_index_settings * index_settings = hercules_index_settings_init(index_path,
										     time_series_size,   
										     init_segments,
										     paa_segments,
										     sax_cardinality, 		  
										     leaf_size,          
										     buffered_memory_size,
										     is_index_new,
										     q_settings.serial,
										     q_settings.sims,
										     flush_threshold);  

       if (index_settings == NULL)
       { 
         fprintf(stderr, "Error main.c:  Could not initialize the index settings.\n");
         return -1;              
       }
    
       index = hercules_index_init(index_settings);
       index->first_node = hercules_root_node_init(index->settings);

       if (index == NULL)
       { 
          fprintf(stderr, "Error main.c:  Could not initialize the index.\n");
          return -1;              
        }
        if (!use_ascii_input) {

	  printf ("Using num_query_threads = %d\n", q_settings.num_query_threads);
	  if (q_settings.num_query_threads == 1)
	    {
	    if (!hercules_index_binary_file(dataset, dataset_size, index))
            { 
               fprintf(stderr, "Error main.c:  Could not build the index.\n");
               return -1;              
            }
            COUNT_PARTIAL_TIME_END
	      index->stats->idx_building_total_time  += partial_time
	      +  index->stats->idx_traverse_tree_total_time 
	      +  index->stats->idx_append_ts_to_leaf_total_time 
	      +  index->stats->idx_evaluate_split_policies_total_time 
	      +  index->stats->idx_split_node_total_time 
	      ;	
	    index->stats->idx_building_input_time  += partial_input_time
	      +  index->stats->idx_traverse_tree_input_time 
	      +  index->stats->idx_append_ts_to_leaf_input_time 
	      +  index->stats->idx_evaluate_split_policies_input_time 
	      +  index->stats->idx_split_node_input_time 
	      ;	
	    index->stats->idx_building_output_time += partial_output_time
	      +  index->stats->idx_traverse_tree_output_time 
	      +  index->stats->idx_append_ts_to_leaf_output_time 
	      +  index->stats->idx_evaluate_split_policies_output_time 
	      +  index->stats->idx_split_node_output_time 
	      ;	
	    index->stats->idx_building_cpu_time    += index->stats->idx_building_total_time 
	      -  index->stats->idx_building_input_time 					 
	      -  index->stats->idx_building_output_time 					 
	      ;	    
	    index->stats->idx_building_seq_input_count   += partial_seq_input_count
	      +  index->stats->idx_traverse_tree_seq_input_count 
	      +  index->stats->idx_append_ts_to_leaf_seq_input_count 
	      +  index->stats->idx_evaluate_split_policies_seq_input_count 
	      +  index->stats->idx_split_node_seq_input_count 
	      ;	
	    index->stats->idx_building_seq_output_count  += partial_seq_output_count
	      +  index->stats->idx_traverse_tree_seq_output_count 
	      +  index->stats->idx_append_ts_to_leaf_seq_output_count 
	      +  index->stats->idx_evaluate_split_policies_seq_output_count 
	      +  index->stats->idx_split_node_seq_output_count 
	      ;	
	    index->stats->idx_building_rand_input_count  += partial_rand_input_count
	      +  index->stats->idx_traverse_tree_rand_input_count 
	      +  index->stats->idx_append_ts_to_leaf_rand_input_count 
	      +  index->stats->idx_evaluate_split_policies_rand_input_count 
	      +  index->stats->idx_split_node_rand_input_count 
	      ;	

	    index->stats->idx_building_rand_output_count += partial_rand_output_count	    
	      +  index->stats->idx_traverse_tree_rand_output_count 
	      +  index->stats->idx_append_ts_to_leaf_rand_output_count 
	      +  index->stats->idx_evaluate_split_policies_rand_output_count 
	      +  index->stats->idx_split_node_rand_output_count 
	      ;		    
	  }
	  else
	  {

	    if (!hercules_index_build(dataset, dataset_size, index,q_settings.num_query_threads,initial_db_size))
            { 
               fprintf(stderr, "Error main.c:  Could not build the index.\n");
               return -1;              
            }	    
        /*	    		
	    if (!hercules_index_binary_file_p(dataset, dataset_size, index,q_settings.num_query_threads,initial_db_size))
            { 
               fprintf(stderr, "Error main.c:  Could not build the index.\n");
               return -1;              
            }	    	    		
        */
	  }

	    RESET_PARTIAL_COUNTERS()
	    COUNT_PARTIAL_TIME_START
        
  	    if (num_write_threads == 1)
	    { 
			hercules_index_write_sequential(index);
        }
        else
        {
			hercules_index_write(index, num_write_threads);
        }
            COUNT_PARTIAL_TIME_END
	      //COUNT_TOTAL_TIME_END	      
	    index->stats->idx_writing_total_time  = partial_time;	
	    index->stats->idx_writing_input_time  = partial_input_time;
	    index->stats->idx_writing_output_time = partial_output_time;
	    index->stats->idx_writing_cpu_time    = partial_time
	                                          - partial_input_time
	                                          - partial_output_time;
	    index->stats->idx_writing_seq_input_count   = partial_seq_input_count;
	    index->stats->idx_writing_seq_output_count  = partial_seq_output_count;
	    index->stats->idx_writing_rand_input_count  = partial_rand_input_count;
	    index->stats->idx_writing_rand_output_count = partial_rand_output_count;	    
        
          hercules_get_index_stats (index);
          hercules_print_index_stats(index, dataset);
		
	    //COUNT_TOTAL_TIME_START	    
        }
        else 
        {
  	   if (!hercules_index_ascii_file(dataset, dataset_size, DELIMITER, index))
           { 
              fprintf(stderr, "Error main.c:  Could not build the index.\n");
              return -1;              
           }
	 
   	   if (!hercules_index_write_sequential(index)) //call it finalize index
           {  
              fprintf(stderr, "Error main.c:  Could not save the index to disk.\n");
              return -1;              
           }		
       }
    
	if (q_settings.serial)
	{
	  fclose (index->leaves_raw_file);
	}      	
	if (q_settings.sims)
	{
	  fclose (index->leaves_sims_file);
	} 
        	      
    }
    else if (mode == 1)  //read an existing index and execute queries
    {
            RESET_PARTIAL_COUNTERS()
            COUNT_PARTIAL_TIME_START
	    is_index_new = 0;
	    index = hercules_index_read(index_path);
  	    if (index == NULL) 
            { 
               fprintf(stderr, "Error main.c:  Could not read the index from disk.\n");
               return -1;              
            }
	    index->settings->serial = q_settings.serial;
            COUNT_PARTIAL_TIME_END
	      //COUNT_TOTAL_TIME_END	      
	    index->stats->idx_reading_total_time  = partial_time;	
	    index->stats->idx_reading_input_time  = partial_input_time;
	    index->stats->idx_reading_output_time = partial_output_time;
	    index->stats->idx_reading_cpu_time    = partial_time
	                                          - partial_input_time
	                                          - partial_output_time;
	    index->stats->idx_reading_seq_input_count   = partial_seq_input_count;
	    index->stats->idx_reading_seq_output_count  = partial_seq_output_count;
	    index->stats->idx_reading_rand_input_count  = partial_rand_input_count;
	    index->stats->idx_reading_rand_output_count = partial_rand_output_count;


            hercules_get_index_stats (index);
            hercules_print_index_stats(index,dataset);
	    //COUNT_TOTAL_TIME_START
	      
	    fprintf(stderr, ">>> Index loaded successfully from: %s\n", index_path);
	    /*
	    if (use_ascii_input)
	    {
	       if (!hercules_query_ascii_file(queries, queries_size, DELIMITER, index, minimum_distance, epsilon, delta))
               { 
                  fprintf(stderr, "Error main.c:  Could not execute the query.\n");
                  return -1;              		  
               }
	    }
	    else
	    {
            if (!hercules_knn_query_binary_file(queries, queries_size, index,
					      minimum_distance, epsilon, delta,k, track_bsf))
               { 
                  fprintf(stderr, "Error main.c:  Could not execute the query.\n");
                  return -1;              
               }
	       }*/

	    //calculate r_delta
	    ts_type r_delta = FLT_MAX;
	    
	    if ( delta > 0 && delta < 1)
	    {
	      FILE *fp;
	      int bins;
	      double *x, *y;
	      int i,j;
	      if((fp=fopen(dataset_hists, "r"))!=NULL)
		{
		  fscanf(fp, "%d", &bins);
		
		  x = (double *) calloc(bins, sizeof(double));
		  y = (double *) calloc(bins, sizeof(double));

		  for(int j=0; j<bins; j++)
		    fscanf(fp, "%lf%lf", &x[j], &y[j]);
		  fclose(fp);		
		}
	      else printf("Error opening data distribution file.\n");

	      for (j=0; (j<bins)&&(y[j]<(1-delta)); j++);
	      j--;
	      fprintf(stderr,"Using histogram bin %lf.\n", y[j]);
			
	      r_delta=x[j]+((((1-delta)-y[j])*(x[j+1]-x[j]))/(y[j+1]-y[j]));
	      fprintf(stderr,"Using r_delta = %lf.\n", r_delta);
           
	    }
	    if (q_settings.sims)
	    {		cache_sax_file(index);
	    }
           //answer delta-epsilon kNN queries, including k=1
	    hercules_knn_query_binary_file(queries, queries_size, index,
					 minimum_distance, epsilon, r_delta,
					 k,track_bsf,track_pruning, all_mindists,
					 max_policy,nprobes, incremental,q_settings, queries_skip);
	    if (q_settings.sims)
	    {
	      fclose(index->leaves_sims_file);
	    }
            /*
	    if (q_settings.serial
	      {
		fclose (index->leaves_raw_file);
	      }      	    
	    */
	    
    }
    else if (mode == 2) //build the index, execute queries and store the index
    { 
       RESET_PARTIAL_COUNTERS()
       COUNT_PARTIAL_TIME_START
	 struct hercules_index_settings * index_settings = hercules_index_settings_init(index_path,
										    time_series_size,   
										    init_segments,
										    paa_segments,
										    sax_cardinality,
										    leaf_size,						 
										    buffered_memory_size,
										    is_index_new,
										    q_settings.serial,
										    q_settings.sims,
										    flush_threshold);  

       if (index_settings == NULL)
       { 
         fprintf(stderr, "Error main.c:  Could not initialize the index settings.\n");
         return -1;              
       }
    
       index = hercules_index_init(index_settings);
       index->first_node = hercules_root_node_init(index->settings);
    
       if (index == NULL)
       { 
          fprintf(stderr, "Error main.c:  Could not initialize the index.\n");
          return -1;              
        }
        if (!use_ascii_input) {
	  
	  if (!hercules_index_binary_file(dataset, dataset_size, index))
            { 
               fprintf(stderr, "Error main.c:  Could not build the index.\n");
               return -1;              
            }
            COUNT_PARTIAL_TIME_END
	      //COUNT_TOTAL_TIME_END	      
	    index->stats->idx_building_total_time  = partial_time;	
	    index->stats->idx_building_input_time  = partial_input_time;
	    index->stats->idx_building_output_time = partial_output_time;
	    index->stats->idx_building_cpu_time    = partial_time
	                                           - partial_input_time
	                                           - partial_output_time;
	    index->stats->idx_building_seq_input_count   = partial_seq_input_count;
	    index->stats->idx_building_seq_output_count  = partial_seq_output_count;
	    index->stats->idx_building_rand_input_count  = partial_rand_input_count;
	    index->stats->idx_building_rand_output_count = partial_rand_output_count;


            hercules_get_index_stats (index);
            hercules_print_index_stats(index,dataset);
	    COUNT_TOTAL_TIME_START
	      
            if (!hercules_query_binary_file(queries, queries_size, index, minimum_distance, epsilon, delta))
            { 
               fprintf(stderr, "Error main.c:  Could not execute the query.\n");
               return -1;              
            }
	       
            RESET_PARTIAL_COUNTERS()
	    COUNT_PARTIAL_TIME_START
  	    if (!hercules_index_write_sequential(index)) 
            { 
               fprintf(stderr, "Error main.c:  Could not save the index to disk.\n");
               return -1;              
            }
            //COUNT_PARTIAL_TIME_END
	    index->stats->idx_writing_total_time  = partial_time;	
	    index->stats->idx_writing_input_time  = partial_input_time;
	    index->stats->idx_writing_output_time = partial_output_time;
	    index->stats->idx_writing_cpu_time    = partial_time
	                                          - partial_input_time
	                                          - partial_output_time;
	    index->stats->idx_writing_seq_input_count   = partial_seq_input_count;
	    index->stats->idx_writing_seq_output_count  = partial_seq_output_count;
	    index->stats->idx_writing_rand_input_count  = partial_rand_input_count;
	    index->stats->idx_writing_rand_output_count = partial_rand_output_count;	    
	    //ADD TIME COUNTERS	    
        }
        else 
        {
	  //ADD TIME COUNTERS HERE
  	   if (!hercules_index_ascii_file(dataset, dataset_size, DELIMITER, index))
           { 
              fprintf(stderr, "Error main.c:  Could not build the index.\n");
              return -1;              
           }
	 
	   if (!hercules_query_ascii_file(queries, queries_size, DELIMITER, index, minimum_distance, epsilon, delta))
           { 
              fprintf(stderr, "Error main.c:  Could not execute the query.\n");
              return -1;              
           }

   	   if (!hercules_index_write_sequential(index)) 
           {  
              fprintf(stderr, "Error main.c:  Could not save the index to disk.\n");
              return -1;              
           }		
       }
    }
    else if (mode == 3)  //read an existing index and execute queries
    {
	    is_index_new = 0;
	    index = hercules_index_read(index_path);
  	    if (index == NULL) 
            { 
               fprintf(stderr, "Error main.c:  Could not read the index from disk.\n");
               return -1;              
            }
	      
	    fprintf(stderr, ">>> Index loaded successfully from: %s\n", index_path);
	    if (!use_ascii_input)
	    {
	       if (!hercules_tlb_binary_file(queries, queries_size, index, minimum_distance))
               { 
                  fprintf(stderr, "Error main.c:  Could not execute the query.\n");
                  return -1;              
               }
	    }	    
    }
    
    else
    {
      fprintf(stderr, "Please use a valid mode. run hercules --help for more information. \n");
      return -1;              
    }
      
    COUNT_TOTAL_TIME_END
    fprintf(stderr,"Sanity check: combined indexing and querying times should be less than: %f secs \n",
      	   total_time/1000000);
    
    hercules_index_destroy(index, index->first_node, is_index_new);

    free(index->stats);
    
    free(index->settings->bit_masks);
    free(index->settings->max_sax_cardinalities);
    
    
    free(index->settings);
    free(index);

    //printf("\n");

    //malloc_stats_print(NULL, NULL, NULL);    
    return 0;
}

