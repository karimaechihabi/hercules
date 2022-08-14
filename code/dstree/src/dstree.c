//
//  main.c
//  ds-tree C version
//  Created by Karima Echihabi on 18/12/2016
//  Copyright 2016 Paris Descartes University. All rights reserved.
//  
//  Based on isax code written by Zoumpatianos on 3/12/12.
//  Copyright 2012 University of Trento. All rights reserved.
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
#include "../include/dstree_file_loaders.h"
#include "../include/dstree_index.h"
#include "../include/dstree_node.h"
#include "../include/dstree_file_buffer.h"
#include "../include/dstree_file_buffer_manager.h"

#ifdef VALUES
#include <values.h>
#endif
		
int main (int argc, char **argv)
{
    INIT_STATS()
    COUNT_TOTAL_TIME_START

    static char * dataset = "/home/karima/myDisk/data/Cgenerator/data_current.txt";
    static char * queries = "/home/karima/myDisk/data/Cgenerator/query_current.txt";
    
    static char * index_path = "out/";
    static unsigned int dataset_size = 1000;
    static unsigned int queries_size = 5;
    static unsigned int time_series_size = 256;
    static unsigned int init_segments = 1;  
    static unsigned int leaf_size = 100;
    static double buffered_memory_size = 6439.2; 
    static int  use_ascii_input = 0;
    static int mode = 0;
    static float minimum_distance = FLT_MAX;
    struct dstree_index * index = NULL;
    boolean is_index_new = 1;
    static int restructure = 0;
    static int serial = 0;
    static int skip_index = 0;        
    static int update_stats = 0;
    
    //printf("new code\n");
    
    while (1) 
    {    
        static struct option long_options[] =  {
            {"ascii-input", required_argument, 0, 'a'},
            {"buffer-size", required_argument, 0, 'b'},
            {"restructure", required_argument, 0, 'c'},
            {"dataset", required_argument, 0, 'd'},
            {"serial", required_argument, 0, 'e'},
            {"skip_index", required_argument, 0, 'f'},
            {"update_stats", required_argument, 0, 'g'},	    	    
            {"queries-size", required_argument, 0, 'k'},
            {"leaf-size", required_argument, 0, 'l'},
            {"index-path", required_argument, 0, 'p'},
            {"queries", required_argument, 0, 'q'},
	    {"minimum-distance", required_argument, 0, 's'},
            {"timeseries-size", required_argument, 0, 't'},
            {"mode", required_argument, 0, 'x'},
            {"dataset-size", required_argument, 0, 'z'},
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

	    case 'b':
                 buffered_memory_size = atof(optarg);
                 break;
		 
            case 'k':
                queries_size = atoi(optarg);
                if (queries_size < 1)
		{
		  fprintf(stderr,"Please change the queries size to be greater than 0.\n");
		  exit(-1);
		}
                break;
		
            case 'd':
                dataset = optarg;
                break;
		
            case 'p':
                index_path = optarg;
                break;
		
            case 'x':
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
		
            case 'l':
                leaf_size = atoi(optarg);
                if (leaf_size <= 1)
		{
		  fprintf(stderr,"Please change the leaf size to be greater than 1.\n");
		  exit(-1);
		}
                break;

            case 'c':
                restructure = atoi(optarg);
                break;
		
            case 'e':
                serial = atoi(optarg);
                break;
            case 'f':
                skip_index = atoi(optarg);
                break;
            case 'g':
                update_stats = atoi(optarg);
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
                       \t--bin/dstree --dataset XX --dataset-size XX             \n\n\
                       \t--          --index-path XX --timeseries-size XX --mode 0\n\n\
                       \t--*********************QUERY MODE*********************\n\n\
                       \t--bin/dstree --queries XX --queries-size XX             \n\n\
                       \t--           --index-path XX --mode 1                 \n\n\
                       \t--*****************INDEX AND QUERY MODE***************\n\n\
                       \t--bin/dstree --dataset XX --dataset-size XX             \n\n\
                       \t--          --timeseries-size XX --index-path XX      \n\n\
                       \t--           --queries XX --queries-size XX --mode 2  \n\n\
                       \t--****************************************************\n\n");
  
                return 0;
                break;
            case 'a':
  	        use_ascii_input = atoi(optarg);
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
       struct dstree_index_settings * index_settings = dstree_index_settings_init(index_path,
                                                                    time_series_size,   
                                                                    init_segments,       
                                                                    leaf_size,          
								    buffered_memory_size,
								    is_index_new);  

       if (index_settings == NULL)
       { 
         fprintf(stderr, "Error main.c:  Could not initialize the index settings.\n");
         return -1;              
       }
    
       index = dstree_index_init(index_settings);
       index->first_node = dstree_root_node_init(index->settings);

       if (index == NULL)
       { 
          fprintf(stderr, "Error main.c:  Could not initialize the index.\n");
          return -1;              
        }
        if (!use_ascii_input) {

	    if (!dstree_index_binary_file(dataset, dataset_size, index))
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


	    RESET_PARTIAL_COUNTERS()
	    COUNT_PARTIAL_TIME_START
  	    if (!dstree_index_write(index)) 
            { 
               fprintf(stderr, "Error main.c:  Could not save the index to disk.\n");
               return -1;              
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

            dstree_get_index_stats (index);
            dstree_print_index_stats(index, dataset);
	    //COUNT_TOTAL_TIME_START	    
        }
        else 
        {
  	   if (!dstree_index_ascii_file(dataset, dataset_size, DELIMITER, index))
           { 
              fprintf(stderr, "Error main.c:  Could not build the index.\n");
              return -1;              
           }
	 
   	   if (!dstree_index_write(index)) //call it finalize index
           {  
              fprintf(stderr, "Error main.c:  Could not save the index to disk.\n");
              return -1;              
           }		
       fclose (index->leaves_raw_file);

       }
      
    }
    else if (mode == 1)  //read an existing index and execute queries
    {
            RESET_PARTIAL_COUNTERS()
            COUNT_PARTIAL_TIME_START
	    is_index_new = 0;
	    index = dstree_index_read(index_path);
  	    if (index == NULL) 
            { 
               fprintf(stderr, "Error main.c:  Could not read the index from disk.\n");
               return -1;              
            }
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

            printf("BEFORE RESTRUCTURE\n");
            dstree_get_index_stats (index);
            dstree_print_index_stats(index,dataset);
	    //COUNT_TOTAL_TIME_START
	    fprintf(stderr, ">>> Index loaded successfully from: %s\n", index_path);

  	    if (restructure) //call it finalize index
            {
   	       if (!dstree_index_restructure_verify(index->first_node))
                { 
                  fprintf(stderr, "Error main.c:  Could not execute the query.\n");
                  return -1;              
                }
 	       if (!dstree_index_restructure_finalize(index,index->first_node))
  		    { 
		      fprintf(stderr, "Error main.c:  Could not execute the query.\n");
		      return -1;              
		   }	       
	    }
        else {   
           /*
	       if (!dstree_index_restructure_finalize(index,index->first_node))
           { 
                  fprintf(stderr, "Error main.c:  Could not execute the query.\n");
                  return -1;              
           }
     	   if (!dstree_index_write(index)) //call it finalize index
           {  
                 fprintf(stderr, "Error main.c:  Could not save the index to disk.\n");
                 return -1;              
           }		
  	       index = dstree_index_read(index_path);
           printf("AFTER RESTRUCTURE\n");
	       dstree_get_index_stats (index);
           dstree_print_index_stats(index,dataset);
          */
	   if (serial) //call it finalize index
           {
	     if (!skip_index)
                dstree_query_leaves_binary_file(queries, queries_size, index, minimum_distance);
	     else
	       dstree_serial_query_leaves_binary_file(queries, queries_size, index, minimum_distance);

	   }   

           if (!serial)
           {
 	     if (use_ascii_input)
	     {
	       if (!dstree_query_ascii_file(queries, queries_size, DELIMITER, index, minimum_distance))
               { 
                  fprintf(stderr, "Error main.c:  Could not execute the query.\n");
                  return -1;              		  
               }
	     }
	     else
	     {
	       if (!dstree_query_binary_file(queries, queries_size, index, minimum_distance))
               { 
                  fprintf(stderr, "Error main.c:  Could not execute the query.\n");
                  return -1;              
               }
	     }	    
          }
       }
	   //write the index back with stats about coaccesses
	   if (restructure || update_stats)
	   {
	     if (!dstree_index_write(index)) //call it finalize index
	       {  
		 fprintf(stderr, "Error main.c:  Could not save the index to disk.\n");
		 return -1;              
	       }
       }
	   
	   //later open and close only if serial
           fclose (index->leaves_raw_file);	     
	   
    }
    else if (mode == 2) //build the index, execute queries and store the index
    { 
       RESET_PARTIAL_COUNTERS()
       COUNT_PARTIAL_TIME_START
       struct dstree_index_settings * index_settings = dstree_index_settings_init(index_path,
                                                                    time_series_size,   
                                                                    init_segments,       
                                                                    leaf_size,          
								    buffered_memory_size,
								    is_index_new);  

       if (index_settings == NULL)
       { 
         fprintf(stderr, "Error main.c:  Could not initialize the index settings.\n");
         return -1;              
       }
    
       index = dstree_index_init(index_settings);
       index->first_node = dstree_root_node_init(index->settings);
    
       if (index == NULL)
       { 
          fprintf(stderr, "Error main.c:  Could not initialize the index.\n");
          return -1;              
        }
        if (!use_ascii_input) {
	  
	  if (!dstree_index_binary_file(dataset, dataset_size, index))
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


            dstree_get_index_stats (index);
            dstree_print_index_stats(index,dataset);
	    COUNT_TOTAL_TIME_START
	      
            if (!dstree_query_binary_file(queries, queries_size, index, minimum_distance))
            { 
               fprintf(stderr, "Error main.c:  Could not execute the query.\n");
               return -1;              
            }
	       
            RESET_PARTIAL_COUNTERS()
	    COUNT_PARTIAL_TIME_START
  	    if (!dstree_index_write(index)) 
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
  	   if (!dstree_index_ascii_file(dataset, dataset_size, DELIMITER, index))
           { 
              fprintf(stderr, "Error main.c:  Could not build the index.\n");
              return -1;              
           }
	 
	   if (!dstree_query_ascii_file(queries, queries_size, DELIMITER, index, minimum_distance))
           { 
              fprintf(stderr, "Error main.c:  Could not execute the query.\n");
              return -1;              
           }

   	   if (!dstree_index_write(index)) 
           {  
              fprintf(stderr, "Error main.c:  Could not save the index to disk.\n");
              return -1;              
           }		
       }
    }
    else if (mode == 3)  //read an existing index and execute queries
    {
	    is_index_new = 0;
	    index = dstree_index_read(index_path);
  	    if (index == NULL) 
            { 
               fprintf(stderr, "Error main.c:  Could not read the index from disk.\n");
               return -1;              
            }
	      
	    fprintf(stderr, ">>> Index loaded successfully from: %s\n", index_path);
	    if (!use_ascii_input)
	    {
	       if (!dstree_tlb_binary_file(queries, queries_size, index, minimum_distance))
               { 
                  fprintf(stderr, "Error main.c:  Could not execute the query.\n");
                  return -1;              
               }
	    }	    
    }
    
    else
    {
      fprintf(stderr, "Please use a valid mode. run dstree --help for more information. \n");
      return -1;              
    }
      
    COUNT_TOTAL_TIME_END
    fprintf(stderr,"Sanity check: combined indexing and querying times should be less than: %f secs \n",
      	   total_time/1000000);
    
    dstree_index_destroy(index, index->first_node, is_index_new);

    free(index->stats->leaves_heights);
    free(index->stats->leaves_sizes);
    free(index->stats->leaves_accesses);
    free(index->stats->leaves_coaccesses);
    free(index->stats);
    free(index->settings);
    free(index);
     
    //printf("\n");

    malloc_stats_print(NULL, NULL, NULL);    
    return 0;
}

