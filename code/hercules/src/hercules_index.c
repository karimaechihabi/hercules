//
//  hercules_index.c
//  Created by Karima Echihabi on 18/12/2018
//

#include "../config.h"
#include "../globals.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef VALUES
#include <values.h>
#endif

#include <sys/stat.h>
#include <float.h>
#include <math.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>
#include "../include/hercules_index.h"
#include "../include/hercules_node.h"
#include "../include/calc_utils.h"
#include "../include/hercules_node_split.h"
#include "../include/hercules_file_buffer_manager.h"
#include "../include/hercules_file_buffer.h"


/**
 This function initializes the settings of a hercules index
 */
struct hercules_index_settings * hercules_index_settings_init(const char * root_directory,
							  unsigned int timeseries_size, 
							  unsigned int init_segments,
							  int paa_segments,
							  int sax_bit_cardinality,				   	  
							  unsigned int max_leaf_size,
							  double buffered_memory_size,
							  boolean is_index_new,
							  int serial,
							  int sims,
							  int flush_threshold)
{

  if (is_index_new)
  {
    if(chdir(root_directory) == 0)
    {
        fprintf(stderr, "WARNING! Target index directory already exists. Please delete or choose a new one.\n");
        exit(-1);
    }
  }
    mkdir(root_directory, 0777);
  
    struct hercules_index_settings *settings = malloc(sizeof(struct hercules_index_settings));
    if(settings == NULL) {
        fprintf(stderr,"Error in hercules_index.c: could not allocate memory for index settings.\n");
        return NULL;
    }


    settings->root_directory = root_directory;
    settings->timeseries_size = timeseries_size;
    settings->init_segments = init_segments;
    settings->max_leaf_size = max_leaf_size;
    settings->buffered_memory_size = buffered_memory_size;
    settings->serial = serial;
    settings->sims = sims;
   
    /* Each leaf node has a file called: 
       numVerticalSegments_LeftorRight_SplitIndicator_(idxFrom,idxTo,splitValue)_level
             
       numVerticalSegments cannot exceed 2 characters (max_segment_length)
       LeftorRight is one character : L or R
       SplitIndicator is 0 (mean) or 1(stdev)
       idxFrom: Start of segment that was split 
       idxTo :  End of segment that was split 
       splitValue: Value used to split segment, (cannot exceed max_value_length characters) 
       level :  int
       number of punctuation marks (underscores:4, parentheses:2, commas:2): total 8
       null character: 1 
    */

     float segment_ends_size = ceil(log10(SHRT_MAX));
     float split_value_size = ceil(log10(INT_MAX)+1);

     settings->max_filename_size = 2+1+1+
                                   2*(segment_ends_size)+
                                   10 + split_value_size+8+1;          


    /* Adding Support for SIMS on leaves*/    
    
    if(paa_segments > (int)(8 * (int)sizeof(root_mask_type))){
        fprintf(stderr,"error: Too many paa segments. The maximum value is %zu.\n", 
                8 * sizeof(root_mask_type));
        return NULL;
    }
        
    settings->paa_segments = paa_segments;
    settings->ts_values_per_paa_segment = (int) timeseries_size/paa_segments;

    /*check if we need to cast to unsigned int*/
    settings->sax_byte_size = (sizeof(sax_type) * paa_segments);
    settings->ts_byte_size = (sizeof(ts_type) * timeseries_size);
    settings->position_byte_size = sizeof(file_position_type);
    
    settings->full_record_size = settings->sax_byte_size 
                                 + settings->position_byte_size 
                                 + settings->ts_byte_size;
    settings->partial_record_size = settings->sax_byte_size 
                                    + settings->position_byte_size;
    
    settings->sax_bit_cardinality = sax_bit_cardinality;
    settings->sax_alphabet_cardinality = pow(2, sax_bit_cardinality);
    
    settings->max_sax_cardinalities = malloc(sizeof(sax_type) * settings->paa_segments);
   
   for(int i=0; i<settings->paa_segments;i++)
     settings->max_sax_cardinalities[i] = settings->sax_bit_cardinality;
	
    //settings->mindist_sqrt = sqrtf((float) settings->timeseries_size /
    //                               (float) settings->paa_segments);
    settings->mindist_sqrt = (settings->timeseries_size /
                                   settings->paa_segments);
    
    // SEGMENTS * (CARDINALITY)
    float c_size = ceil(log10(settings->sax_alphabet_cardinality + 1));
    settings->max_filename_size = settings->paa_segments * 
                                  ((c_size * 2) + 2)
                                  + 5 + strlen(root_directory);
    
    
    if(paa_segments > sax_bit_cardinality)
    {
        settings->bit_masks = malloc(sizeof(root_mask_type) * (paa_segments+1));
        if(settings->bit_masks == NULL) {
            fprintf(stderr,"error: could not allocate memory for bit masks.\n");
            return NULL;
        }
        
        for (; paa_segments>=0; paa_segments--)
        {
            settings->bit_masks[paa_segments] = pow(2, paa_segments);
        }
    }
    else
    {
        settings->bit_masks = malloc(sizeof(root_mask_type) * (sax_bit_cardinality+1));
        if(settings->bit_masks == NULL) {
            fprintf(stderr,"error: could not allocate memory for bit masks.\n");
            return NULL;
        }
        
        for (; sax_bit_cardinality>=0; sax_bit_cardinality--)
        {
            settings->bit_masks[sax_bit_cardinality] = pow(2, sax_bit_cardinality);
        }
    }
    settings->flush_threshold = flush_threshold;
    
    return settings;
}



/**
 This function initializes an isax index
 @param isax_index_settings *settings
 @return isax_index
 */
struct hercules_index * hercules_index_init(struct hercules_index_settings *settings)
{
    struct hercules_index *index = malloc(sizeof(struct hercules_index));
    if(index == NULL) {
        fprintf(stderr,"Error in hercules_index.c: Could not allocate memory for index structure.\n");
        return NULL;
    }

    
    index->settings = settings;
    index->total_records = 0;
    index->sax_cache = NULL;
  
    hercules_init_stats(index);
    
    if (!init_file_buffer_manager(index))
    { 
      fprintf(stderr, "Error in hercules_index.c:  Could not initialize the \
                       file buffer manager for this index.\n");
      return NULL;              
    }
	
    return index;
}

enum response hercules_init_stats(struct hercules_index * index)
{
    index->stats = malloc(sizeof(struct stats_info));
    if(index->stats == NULL) {
        fprintf(stderr,"Error in hercules_index.c: Could not allocate memory for stats structure.\n");
        return FAILURE;
    }

    /*INDEX STATISTICS*/
    index->stats->idx_traverse_tree_input_time = 0;  
    index->stats->idx_traverse_tree_output_time = 0;
    index->stats->idx_traverse_tree_cpu_time = 0;
    index->stats->idx_traverse_tree_total_time = 0;

    index->stats->idx_traverse_tree_seq_input_count = 0;
    index->stats->idx_traverse_tree_seq_output_count = 0;
    index->stats->idx_traverse_tree_rand_input_count = 0;
    index->stats->idx_traverse_tree_rand_output_count = 0;    

    index->stats->idx_append_ts_to_leaf_input_time = 0;  
    index->stats->idx_append_ts_to_leaf_output_time = 0;
    index->stats->idx_append_ts_to_leaf_cpu_time = 0;
    index->stats->idx_append_ts_to_leaf_total_time = 0;

    index->stats->idx_append_ts_to_leaf_seq_input_count = 0;
    index->stats->idx_append_ts_to_leaf_seq_output_count = 0;
    index->stats->idx_append_ts_to_leaf_rand_input_count = 0;
    index->stats->idx_append_ts_to_leaf_rand_output_count = 0;    

    index->stats->idx_evaluate_split_policies_input_time = 0;  
    index->stats->idx_evaluate_split_policies_output_time = 0;
    index->stats->idx_evaluate_split_policies_cpu_time = 0;
    index->stats->idx_evaluate_split_policies_total_time = 0;

    index->stats->idx_evaluate_split_policies_seq_input_count = 0;
    index->stats->idx_evaluate_split_policies_seq_output_count = 0;
    index->stats->idx_evaluate_split_policies_rand_input_count = 0;
    index->stats->idx_evaluate_split_policies_rand_output_count = 0;    

    index->stats->idx_split_node_input_time = 0;  
    index->stats->idx_split_node_output_time = 0;
    index->stats->idx_split_node_cpu_time = 0;
    index->stats->idx_split_node_total_time = 0;

    index->stats->idx_split_node_seq_input_count = 0;
    index->stats->idx_split_node_seq_output_count = 0;
    index->stats->idx_split_node_rand_input_count = 0;
    index->stats->idx_split_node_rand_output_count = 0;    

    index->stats->idx_building_input_time = 0;  
    index->stats->idx_building_output_time = 0;
    index->stats->idx_building_cpu_time = 0;
    index->stats->idx_building_total_time = 0;

    index->stats->idx_building_seq_input_count = 0;
    index->stats->idx_building_seq_output_count = 0;
    index->stats->idx_building_rand_input_count = 0;
    index->stats->idx_building_rand_output_count = 0;    

    index->stats->idx_writing_input_time = 0;  
    index->stats->idx_writing_output_time = 0;
    index->stats->idx_writing_cpu_time = 0;
    index->stats->idx_writing_total_time = 0;

    index->stats->idx_writing_seq_input_count = 0;
    index->stats->idx_writing_seq_output_count = 0;
    index->stats->idx_writing_rand_input_count = 0;
    index->stats->idx_writing_rand_output_count = 0;

    index->stats->idx_reading_input_time = 0;  
    index->stats->idx_reading_output_time = 0;
    index->stats->idx_reading_cpu_time = 0;
    index->stats->idx_reading_total_time = 0;
  
    index->stats->idx_reading_seq_input_count = 0;
    index->stats->idx_reading_seq_output_count = 0;
    index->stats->idx_reading_rand_input_count = 0;
    index->stats->idx_reading_rand_output_count = 0;
    
    index->stats->idx_total_input_time = 0;
    index->stats->idx_total_output_time = 0;
    index->stats->idx_total_cpu_time = 0;
    index->stats->idx_total_time = 0;

    index->stats->idx_total_seq_input_count = 0;
    index->stats->idx_total_seq_output_count = 0;
    index->stats->idx_total_rand_input_count = 0;
    index->stats->idx_total_rand_output_count = 0;    

    index->stats->total_nodes_count = 0;
    index->stats->leaf_nodes_count = 0;
    index->stats->empty_leaf_nodes_count = 0;
  
    index->stats->idx_size_bytes = 0;
    index->stats->idx_size_blocks = 0;
    
    index->stats->min_fill_factor = FLT_MAX;
    index->stats->max_fill_factor = 0;
    index->stats->sum_fill_factor = 0;
    index->stats->sum_squares_fill_factor = 0;
    index->stats->avg_fill_factor = 0;
    index->stats->sd_fill_factor = 0;

    index->stats->min_height = FLT_MAX;
    index->stats->max_height = 0;
    index->stats->sum_height = 0;
    index->stats->sum_squares_height = 0;    
    index->stats->avg_height = 0;
    index->stats->sd_height = 0;
    
    /*PER QUERY STATISTICS*/    
    index->stats->query_filter_input_time = 0;
    index->stats->query_filter_output_time = 0;
    index->stats->query_filter_load_node_time = 0;
    index->stats->query_filter_cpu_time = 0;    
    index->stats->query_filter_total_time = 0;    

    index->stats->query_filter_seq_input_count = 0;
    index->stats->query_filter_seq_output_count = 0;
    index->stats->query_filter_rand_input_count = 0;
    index->stats->query_filter_rand_output_count = 0;    

    index->stats->query_filter_loaded_nodes_count = 0;
    index->stats->query_filter_checked_nodes_count = 0;    
    index->stats->query_filter_loaded_ts_count = 0;
    index->stats->query_filter_checked_ts_count = 0;

    index->stats->query_refine_input_time = 0;
    index->stats->query_refine_output_time = 0;
    index->stats->query_refine_load_node_time = 0;
    index->stats->query_refine_cpu_time = 0;
    index->stats->query_refine_total_time = 0;    

    index->stats->query_refine_seq_input_count = 0;
    index->stats->query_refine_seq_output_count = 0;
    index->stats->query_refine_rand_input_count = 0;
    index->stats->query_refine_rand_output_count = 0;    

    index->stats->query_refine_loaded_nodes_count = 0;
    index->stats->query_refine_checked_nodes_count = 0;    
    index->stats->query_refine_loaded_ts_count = 0;
    index->stats->query_refine_checked_ts_count = 0;

    index->stats->query_total_input_time = 0;
    index->stats->query_total_output_time = 0;
    index->stats->query_total_load_node_time = 0;
    index->stats->query_total_cpu_time = 0;
    index->stats->query_total_time = 0;    

    index->stats->query_total_seq_input_count = 0;
    index->stats->query_total_seq_output_count = 0;
    index->stats->query_total_rand_input_count = 0;
    index->stats->query_total_rand_output_count = 0;
    
    index->stats->query_total_loaded_nodes_count = 0;
    index->stats->query_total_checked_nodes_count = 0;    
    index->stats->query_total_loaded_ts_count = 0;
    index->stats->query_total_checked_ts_count = 0;
    
    index->stats->query_exact_distance = 0;
    index->stats->query_exact_node_filename = NULL;
    index->stats->query_exact_node_size = 0;    
    index->stats->query_exact_node_level = 0;
    index->stats->query_exact_node_file_pos = 0;
    
    index->stats->query_approx_distance = 0;
    index->stats->query_approx_node_filename = NULL;
    index->stats->query_approx_node_size = 0;
    index->stats->query_approx_node_level = 0;

    index->stats->query_lb_distance = 0;

    index->stats->query_tlb = 0;        
    index->stats->query_eff_epsilon = 0;    
    index->stats->query_pruning_ratio = 0;

    
    /*SUMMARY STATISTICS FOR ALL QUERIES*/        
    index->stats->queries_refine_input_time = 0;
    index->stats->queries_refine_output_time = 0;
    index->stats->queries_refine_load_node_time = 0;
    index->stats->queries_refine_cpu_time = 0;
    index->stats->queries_refine_total_time = 0;    

    index->stats->queries_refine_seq_input_count = 0;
    index->stats->queries_refine_seq_output_count = 0;
    index->stats->queries_refine_rand_input_count = 0;
    index->stats->queries_refine_rand_output_count = 0;        
    
    index->stats->queries_filter_input_time = 0;
    index->stats->queries_filter_output_time = 0;
    index->stats->queries_filter_load_node_time = 0;
    index->stats->queries_filter_cpu_time = 0;
    index->stats->queries_filter_total_time = 0;    

    index->stats->queries_filter_seq_input_count = 0;
    index->stats->queries_filter_seq_output_count = 0;
    index->stats->queries_filter_rand_input_count = 0;
    index->stats->queries_filter_rand_output_count = 0;    

    index->stats->queries_total_input_time = 0;
    index->stats->queries_total_output_time = 0;
    index->stats->queries_total_load_node_time = 0;    
    index->stats->queries_total_cpu_time = 0;
    index->stats->queries_total_time = 0;    

    index->stats->queries_total_seq_input_count = 0;
    index->stats->queries_total_seq_output_count = 0;
    index->stats->queries_total_rand_input_count = 0;
    index->stats->queries_total_rand_output_count = 0;        

    index->stats->queries_min_eff_epsilon = FLT_MAX;
    index->stats->queries_max_eff_epsilon = 0;
    index->stats->queries_sum_eff_epsilon = 0;
    index->stats->queries_sum_squares_eff_epsilon = 0;
    index->stats->queries_avg_eff_epsilon = 0;
    index->stats->queries_sd_eff_epsilon = 0;
    
    index->stats->queries_min_pruning_ratio =  FLT_MAX;
    index->stats->queries_max_pruning_ratio = 0;
    index->stats->queries_sum_pruning_ratio = 0;
    index->stats->queries_sum_squares_pruning_ratio = 0;
    index->stats->queries_avg_pruning_ratio = 0;
    index->stats->queries_sd_pruning_ratio = 0;

    index->stats->queries_min_tlb =  FLT_MAX;
    index->stats->queries_max_tlb = 0;
    index->stats->queries_sum_tlb = 0;
    index->stats->queries_sum_squares_tlb = 0;
    index->stats->queries_avg_tlb = 0;
    index->stats->queries_sd_tlb = 0;    


    index->stats->tlb_ts_count = 0;
    index->stats->eff_epsilon_queries_count = 0;

    
//    index->stats->total_queries_count = 0;
    
    /*COMBINED STATISTICS FOR INDEXING AND QUERY WORKLOAD*/            
    index->stats->total_input_time = 0;
    index->stats->total_output_time = 0;
    index->stats->total_load_node_time = 0;
    index->stats->total_cpu_time = 0;
    index->stats->total_time = 0;
    index->stats->total_time_sanity_check = 0;
    
    index->stats->total_seq_input_count = 0;
    index->stats->total_seq_output_count = 0;
    index->stats->total_rand_input_count = 0;
    index->stats->total_rand_output_count = 0;
    
    index->stats->total_parse_time = 0;	
    index->stats->total_ts_count = 0;	
    
    return SUCCESS;
}

void hercules_get_index_stats(struct hercules_index * index)
{
  /*
  index->stats->total_seq_input_count = index->stats->idx_building_seq_input_count
                                      + index->stats->idx_writing_seq_input_count
                                      + index->stats->idx_reading_seq_input_count;
                                      //    + index->stats->queries_total_seq_input_count;
  index->stats->total_seq_output_count = index->stats->idx_building_seq_output_count
                                      + index->stats->idx_writing_seq_output_count
                                      + index->stats->idx_reading_seq_output_count;
                                      //+ index->stats->queries_total_seq_output_count;
  index->stats->total_rand_input_count = index->stats->idx_building_rand_input_count
                                      + index->stats->idx_writing_rand_input_count
                                      + index->stats->idx_reading_rand_input_count;
                                      // + index->stats->queries_total_rand_input_count;
  index->stats->total_rand_output_count = index->stats->idx_building_rand_output_count
                                      + index->stats->idx_writing_rand_output_count
                                      + index->stats->idx_reading_rand_output_count;
                                         //+ index->stats->queries_total_rand_output_count;
					 */  
  index->stats->total_input_time = index->stats->idx_building_input_time
                                 + index->stats->idx_writing_input_time
                                  + index->stats->idx_reading_input_time;
                                 //   + index->stats->queries_total_input_time;
  index->stats->total_output_time = index->stats->idx_building_output_time
                                  + index->stats->idx_writing_output_time
                                  + index->stats->idx_reading_output_time;
                                  //    + index->stats->queries_total_output_time;  
  index->stats->total_cpu_time    = index->stats->idx_building_cpu_time
                                  + index->stats->idx_writing_cpu_time
                                  + index->stats->idx_reading_cpu_time;
                                 //    + index->stats->queries_total_cpu_time;

  index->stats->total_time    = index->stats->idx_building_total_time
                              + index->stats->idx_writing_total_time
                              + index->stats->idx_reading_total_time;

  //index->stats->total_time_sanity_check = total_time;
  
  //index->stats->load_node_time = load_node_time;
  index->stats->total_parse_time = total_parse_time;

  //index->stats->loaded_nodes_count = loaded_nodes_count;
  index->stats->leaf_nodes_count = leaf_nodes_count;
  index->stats->empty_leaf_nodes_count = empty_leaf_nodes_count;
  
  //index->stats->checked_nodes_count = checked_nodes_count;
  index->stats->total_nodes_count = total_nodes_count;
  index->stats->total_ts_count = total_ts_count;

  hercules_get_index_footprint(index);

}


void hercules_get_index_footprint(struct hercules_index * index)
{

    const char *filename = malloc(sizeof(char) * (strlen(index->settings->root_directory) + 9));
    filename = strcpy(filename, index->settings->root_directory);
    filename = strcat(filename, "root.idx\0");

    struct stat st;
    unsigned int  count_leaves;

    if (stat(filename, &st) == 0)
    {
      index->stats->idx_size_bytes = (long long) st.st_size;
      index->stats->idx_size_blocks = (long long) st.st_blksize;
    }

    count_leaves = index->stats->leaf_nodes_count;
    
    index->stats->avg_fill_factor =  ((double) index->stats->sum_fill_factor) / count_leaves  ;
    index->stats->sum_squares_fill_factor -= (pow( index->stats->sum_fill_factor,2) / count_leaves);
    index->stats->sd_fill_factor =  sqrt(((double) index->stats->sum_squares_fill_factor) / count_leaves);

    index->stats->avg_height     =  ((double) index->stats->sum_height) / count_leaves;
    index->stats->sum_squares_height -= (pow( index->stats->sum_height,2) / count_leaves);
    index->stats->sd_height =  sqrt(((double) index->stats->sum_squares_height) / count_leaves);

    free(filename);
}

  
void hercules_print_index_stats(struct hercules_index * index, char * dataset)
{
        /*
        printf("------------------------ \n");    
        printf("INDEX SUMMARY STATISTICS \n");
        printf("------------------------ \n");
        */	
      //  id = -1 for index and id = query_id for queries
        int id = -1;
        printf("Index_traverse_tree_input_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_traverse_tree_input_time/1000000,
	       dataset,
	       id,
	       id); 
        printf("Index_traverse_tree_output_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_traverse_tree_output_time/1000000,
	       dataset,
	       id,
	       id);
        printf("Index_traverse_tree_cpu_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_traverse_tree_cpu_time/1000000,
	       dataset,
	       id,
	       id);	
        printf("Index_traverse_tree_total_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_traverse_tree_total_time/1000000,
	       dataset,
	       id,
	       id);
	
        printf("Index_traverse_tree_seq_input_count\t%llu\t%s\t%d\n",
	       index->stats->idx_traverse_tree_seq_input_count,
	       dataset,
	       id,
	       id);
	
        printf("Index_traverse_tree_seq_output_count\t%llu\t%s\t%d\n",
	       index->stats->idx_traverse_tree_seq_output_count,
	       dataset,
	       id,
	       id);
	
        printf("Index_traverse_tree_rand_input_count\t%llu\t%s\t%d\n",
	       index->stats->idx_traverse_tree_rand_input_count,
	       dataset,
	       id,
	       id);
	
        printf("Index_traverse_tree_rand_output_count\t%llu\t%s\t%d\n",
	       index->stats->idx_traverse_tree_rand_output_count,
	       dataset,
	       id,
	       id); 

        printf("Index_append_ts_to_leaf_input_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_append_ts_to_leaf_input_time/1000000,
	       dataset,
	       id,
	       id); 
        printf("Index_append_ts_to_leaf_output_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_append_ts_to_leaf_output_time/1000000,
	       dataset,
	       id,
	       id);
        printf("Index_append_ts_to_leaf_cpu_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_append_ts_to_leaf_cpu_time/1000000,
	       dataset,
	       id,
	       id);	
        printf("Index_append_ts_to_leaf_total_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_append_ts_to_leaf_total_time/1000000,
	       dataset,
	       id,
	       id);
	
        printf("Index_append_ts_to_leaf_seq_input_count\t%llu\t%s\t%d\n",
	       index->stats->idx_append_ts_to_leaf_seq_input_count,
	       dataset,
	       id,
	       id);
	
        printf("Index_append_ts_to_leaf_seq_output_count\t%llu\t%s\t%d\n",
	       index->stats->idx_append_ts_to_leaf_seq_output_count,
	       dataset,
	       id,
	       id);
	
        printf("Index_append_ts_to_leaf_rand_input_count\t%llu\t%s\t%d\n",
	       index->stats->idx_append_ts_to_leaf_rand_input_count,
	       dataset,
	       id,
	       id);
	
        printf("Index_append_ts_to_leaf_rand_output_count\t%llu\t%s\t%d\n",
	       index->stats->idx_append_ts_to_leaf_rand_output_count,
	       dataset,
	       id,
	       id); 


        printf("Index_evaluate_split_policies_input_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_evaluate_split_policies_input_time/1000000,
	       dataset,
	       id,
	       id); 
        printf("Index_evaluate_split_policies_output_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_evaluate_split_policies_output_time/1000000,
	       dataset,
	       id,
	       id);
        printf("Index_evaluate_split_policies_cpu_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_evaluate_split_policies_cpu_time/1000000,
	       dataset,
	       id,
	       id);	
        printf("Index_evaluate_split_policies_total_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_evaluate_split_policies_total_time/1000000,
	       dataset,
	       id,
	       id);
	
        printf("Index_evaluate_split_policies_seq_input_count\t%llu\t%s\t%d\n",
	       index->stats->idx_evaluate_split_policies_seq_input_count,
	       dataset,
	       id,
	       id);
	
        printf("Index_evaluate_split_policies_seq_output_count\t%llu\t%s\t%d\n",
	       index->stats->idx_evaluate_split_policies_seq_output_count,
	       dataset,
	       id,
	       id);
	
        printf("Index_evaluate_split_policies_rand_input_count\t%llu\t%s\t%d\n",
	       index->stats->idx_evaluate_split_policies_rand_input_count,
	       dataset,
	       id,
	       id);
	
        printf("Index_evaluate_split_policies_rand_output_count\t%llu\t%s\t%d\n",
	       index->stats->idx_evaluate_split_policies_rand_output_count,
	       dataset,
	       id,
	       id); 



        printf("Index_split_node_input_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_split_node_input_time/1000000,
	       dataset,
	       id,
	       id); 
        printf("Index_split_node_output_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_split_node_output_time/1000000,
	       dataset,
	       id,
	       id);
        printf("Index_split_node_cpu_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_split_node_cpu_time/1000000,
	       dataset,
	       id,
	       id);	
        printf("Index_split_node_total_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_split_node_total_time/1000000,
	       dataset,
	       id,
	       id);
	
        printf("Index_split_node_seq_input_count\t%llu\t%s\t%d\n",
	       index->stats->idx_split_node_seq_input_count,
	       dataset,
	       id,
	       id);
	
        printf("Index_split_node_seq_output_count\t%llu\t%s\t%d\n",
	       index->stats->idx_split_node_seq_output_count,
	       dataset,
	       id,
	       id);
	
        printf("Index_split_node_rand_input_count\t%llu\t%s\t%d\n",
	       index->stats->idx_split_node_rand_input_count,
	       dataset,
	       id,
	       id);
	
        printf("Index_split_node_rand_output_count\t%llu\t%s\t%d\n",
	       index->stats->idx_split_node_rand_output_count,
	       dataset,
	       id,
	       id); 


        printf("Index_flush_worker_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_flush_worker_time/1000000,
	       dataset,
	       id,
	       id); 

        printf("Index_flush_coordinator_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_flush_coordinator_time/1000000,
	       dataset,
	       id,
	       id); 

        printf("Index_insert_to_node_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_insert_to_node_time/1000000,
	       dataset,
	       id,
	       id); 

        printf("Index_reinsert_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_reinsert_time/1000000,
	       dataset,
	       id,
	       id); 

        printf("Index_sanity_counter_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_sanity_counter_time/1000000,
	       dataset,
	       id,
	       id); 

        printf("Index_bulk_loading_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_bulk_loading_time/1000000,
	       dataset,
	       id,
	       id); 

        printf("Index_building_input_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_building_input_time/1000000,
	       dataset,
	       id,
	       id); 
        printf("Index_building_output_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_building_output_time/1000000,
	       dataset,
	       id,
	       id);
        printf("Index_building_cpu_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_building_cpu_time/1000000,
	       dataset,
	       id,
	       id);	
        printf("Index_building_total_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_building_total_time/1000000,
	       dataset,
	       id,
	       id);
	
        printf("Index_building_seq_input_count\t%llu\t%s\t%d\n",
	       index->stats->idx_building_seq_input_count,
	       dataset,
	       id,
	       id);
	
        printf("Index_building_seq_output_count\t%llu\t%s\t%d\n",
	       index->stats->idx_building_seq_output_count,
	       dataset,
	       id,
	       id);
	
        printf("Index_building_rand_input_count\t%llu\t%s\t%d\n",
	       index->stats->idx_building_rand_input_count,
	       dataset,
	       id,
	       id);
	
        printf("Index_building_rand_output_count\t%llu\t%s\t%d\n",
	       index->stats->idx_building_rand_output_count,
	       dataset,
	       id,
	       id); 

        printf("Index_writing_input_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_writing_input_time/1000000,
	       dataset,
	       id,
	       id);

        printf("Index_writing_output_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_writing_output_time/1000000,
	       dataset,
	       id,
	       id);

        printf("Index_writing_cpu_time_secs\t%lf\t%s\t%d\n",
               index->stats->idx_writing_cpu_time/1000000,
	       dataset,
	       id,
	       id);
	
        printf("Index_writing_total_time_secs\t%lf\t%s\t%d\n",
               index->stats->idx_writing_total_time/1000000,
	       dataset,
	       id,
	       id);	

        printf("Index_writing_seq_input_count\t%llu\t%s\t%d\n",
	       index->stats->idx_writing_seq_input_count,
	       dataset,
	       id,
	       id);

        printf("Index_writing_seq_output_count\t%llu\t%s\t%d\n",
	       index->stats->idx_writing_seq_output_count,
	       dataset,
	       id,
	       id);

        printf("Index_writing_rand_input_count\t%llu\t%s\t%d\n",
	       index->stats->idx_writing_rand_input_count,
	       dataset,
	       id,
	       id);

        printf("Index_writing_rand_output_count\t%llu\t%s\t%d\n",
	       index->stats->idx_writing_rand_output_count,
	       dataset,
	       id,
	       id);	

        printf("Index_reading_input_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_reading_input_time/1000000,
	       dataset,
	       id,
	       id
	     ); 
        printf("Index_reading_output_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_reading_output_time/1000000,
	       dataset,
	       id,
	       id);
        printf("Index_reading_cpu_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_reading_cpu_time/1000000,
	       dataset,
	       id,
	       id);	
        printf("Index_reading_total_time_secs\t%lf\t%s\t%d\n",
	       index->stats->idx_reading_total_time/1000000,
	       dataset,
	       id,
	       id);
	
        printf("Index_reading_seq_input_count\t%llu\t%s\t%d\n",
	       index->stats->idx_reading_seq_input_count,
	       dataset,
	       id,
	       id);
	
        printf("Index_reading_seq_output_count\t%llu\t%s\t%d\n",
	       index->stats->idx_reading_seq_output_count,
	       dataset,
	       id,
	       id);
	
        printf("Index_reading_rand_input_count\t%llu\t%s\t%d\n",
	       index->stats->idx_reading_rand_input_count,
	       dataset,
	       id,
	       id);
	
        printf("Index_reading_rand_output_count\t%llu\t%s\t%d\n",
	       index->stats->idx_reading_rand_output_count,
	       dataset,
	       id,
	       id);
	
        printf("Index_total_input_time_secs\t%lf\t%s\t%d\n",
	       index->stats->total_input_time/1000000,
	       dataset,
	       id,
	       id
	     ); 
        printf("Index_total_output_time_secs\t%lf\t%s\t%d\n",
	       index->stats->total_output_time/1000000,
	       dataset,
	       id,
	       id);
        printf("Index_total_cpu_time_secs\t%lf\t%s\t%d\n",
	       index->stats->total_cpu_time/1000000,
	       dataset,
	       id,
	       id);	
        printf("Index_total_time_secs\t%lf\t%s\t%d\n",
	       index->stats->total_time/1000000,
	       dataset,
	       id,
	       id);
	
        printf("Index_total_seq_input_count\t%llu\t%s\t%d\n",
	       index->stats->total_seq_input_count,
	       dataset,
	       id,
	       id);
	
        printf("Index_total_seq_output_count\t%llu\t%s\t%d\n",
	       index->stats->total_seq_output_count,
	       dataset,
	       id,
	       id);
	
        printf("Index_total_rand_input_count\t%llu\t%s\t%d\n",
	       index->stats->total_rand_input_count,
	       dataset,
	       id,
	       id);
	
        printf("Index_total_rand_output_count\t%llu\t%s\t%d\n",
	       index->stats->total_rand_output_count,
	       dataset,
	       id,
	       id); 

        printf("Internal_nodes_count\t%lu\t%s\t%d\n",
	       (index->stats->total_nodes_count - index->stats->leaf_nodes_count),
	       dataset,
	       id,
	       id);

        printf("Leaf_nodes_count\t%lu\t%s\t%d\n",
	       index->stats->leaf_nodes_count,
	       dataset,
	       id,
	       id);

	printf("Empty_leaf_nodes_count\t%lu\t%s\t%d\n",
               index->stats->empty_leaf_nodes_count,	       
	       dataset,
	       id,
	       id);

	printf("Total_nodes_count\t%lu\t%s\t%d\n",
	       index->stats->total_nodes_count,
	       dataset,
	       id,
	       id);

	double size_MB =  (index->stats->idx_size_bytes)*1.0/(1024*1024);

	printf("Index_size_MB\t%lf\t%s\t%d\n",
	       size_MB,
	       dataset,
	       id,
	       id);

	printf("Minimum_fill_factor\t%f\t%s\t%d\n",
	       index->stats->min_fill_factor,	       
	       dataset,
	       id,
	       id);	

	printf("Maximum_fill_factor\t%f\t%s\t%d\n",
	       index->stats->max_fill_factor,	       
	       dataset,
	       id,
	       id);

	printf("Average_fill_factor\t%f\t%s\t%d\n",
	       index->stats->avg_fill_factor,	       
	       dataset,
	       id,
	       id);	

	printf("SD_height\t%f\t%s\t%d\n",
	       index->stats->sd_height,	       
	       dataset,
	       id,
	       id);		

	printf("Minimum_height\t%u\t%s\t%d\n",
	       index->stats->min_height,	       
	       dataset,
	       id,
	       id);	

	printf("Maximum_height\t%u\t%s\t%d\n",
	       index->stats->max_height,	       
	       dataset,
	       id,
	       id);

	printf("Average_height\t%f\t%s\t%d\n",
	       index->stats->avg_height,	       
	       dataset,
	       id,
	       id);	

	printf("SD_height\t%f\t%s\t%d\n",
	       index->stats->sd_height,	       
	       dataset,
	       id,
	       id);			

	printf("Total_ts_count\t%u\t%s\t%d\n",
	       index->stats->total_ts_count,	       
	       dataset,
	       id,
	       id);

	for (int i = 0; i < index->stats->leaves_counter; ++i)	  
	{
 	    double fill_factor = ((double) index->stats->leaves_sizes[i])/index->settings->max_leaf_size;
  	    printf("Leaf_report_node_%d \t Height  %d  \t%s\t%d\n",
		   (i+1),
		   index->stats->leaves_heights[i],
	           dataset,
	           id,
		   id);
  	    printf("Leaf_report_node_%d \t Fill_Factor  %f \t%s\t%d\n",
		   (i+1),
		   fill_factor,	       
	           dataset,
	           id,
		   id);	    
	}

	

}




void print_tlb_stats(struct hercules_index * index, unsigned int query_num, char * queries)
{

        printf("Query_avg_tlb\t%lf\t%s\t%u\n",
	       total_tlb/total_ts_count,	       
	       queries,
	       query_num
	     );
	printf("Leaf_nodes_count\t%u\t%s\t%d\n",
	       leaf_nodes_count,	       
	       queries,
	       query_num);	
	printf("Total_ts_count\t%u\t%s\t%d\n",
	       total_ts_count,	       
	       queries,
	       query_num);
	

}

void flush_index_to_disk(struct hercules_index *index, struct hercules_node * node)
{
   if(!node->is_leaf) 
   {
     flush_index_to_disk(index, node->left_child);
     flush_index_to_disk(index, node->right_child);
   }
   else
   {
     flush_buffer_to_disk_parallel(index, node);
   }
}


void hercules_index_destroy(struct hercules_index *index, struct hercules_node *node, boolean is_index_new)
{

   if (node->level == 0) //root
   {
     if (node->node_segment_split_policies != NULL)
         free(node->node_segment_split_policies);
     if (index->buffer_manager != NULL)
	 destroy_buffer_manager(index);
     if (index->sax_cache != NULL)
       free(index->sax_cache);
     if (index->leaves != NULL)
       free (index->leaves);
     if (index->stats->leaves_heights != NULL)
       free (index->stats->leaves_heights);
     if (index->stats->leaves_sizes != NULL)
       free (index->stats->leaves_sizes);     
     if (index->leaves_sims_filename != NULL)
       free (index->leaves_sims_filename);
     if (index->leaves_raw_filename != NULL)
       free (index->leaves_raw_filename);

   }
   
   if(!node->is_leaf) {
     hercules_index_destroy(index, node->right_child, is_index_new);
     hercules_index_destroy(index, node->left_child,is_index_new);
   }
	
   if(node->split_policy != NULL)
   {
       free(node->split_policy);
   }
	
   if(node->filename != NULL)
   {
      free(node->filename);
   }
	
   if(node->file_buffer != NULL)
   {
       free(node->file_buffer->buffered_list);
       node->file_buffer->buffered_list = NULL;
       node->file_buffer->buffered_list_size = 0;
       free(node->file_buffer);
    }
	
    if(node->node_segment_sketches != NULL)
    {
  	   //free indicators first
       for (int i = 0; i < node->num_node_points; ++i)
       { 
	       free(node->node_segment_sketches[i].indicators);
       }

       free(node->node_segment_sketches);

     }

     if (node->node_points != NULL)
     {
	   free(node->node_points);	 
     }

     if (is_index_new)
       {
     if(node->hs_node_segment_sketches != NULL)
     {
        //free indicators first
        for (int i = 0; i < node->num_hs_node_points; ++i)
	{ 
	   free(node->hs_node_segment_sketches[i].indicators);
	}

	free(node->hs_node_segment_sketches);
     }


     if (node->hs_node_points != NULL)
     {
	free(node->hs_node_points);
     }
   }
    pthread_mutex_destroy(&(node->lock_data));
    free(node);


}

void destroy_buffer_manager(struct hercules_index *index)
{

  if(index->buffer_manager != NULL)
  {
    struct hercules_file_map *currP;
    struct hercules_file_map *temp;

    temp = NULL;
    currP = index->buffer_manager->file_map; 

    while(currP != NULL)
    {
      temp = currP;
      currP = currP->next;      
      free(temp);
    }

    free(index->buffer_manager->mem_array);     

    free(index->buffer_manager);
  }
 
}

enum response hercules_index_insert(struct hercules_index *index,  ts_type * timeseries)
{

  //traverse the index tree to find the appropriate node
  struct hercules_node * node = index->first_node;

  //Allocate memory for the series sketch
  struct segment_sketch timeseries_segment_sketch;
  timeseries_segment_sketch.indicators = NULL;
  timeseries_segment_sketch.indicators = malloc (sizeof(ts_type)*2);
  if(timeseries_segment_sketch.indicators == NULL) {
        fprintf(stderr,"Error in hercules_index.c: Could not allocate memory for" 
                        "series segment sketch indicators.\n");
  }

  timeseries_segment_sketch.num_indicators = 2;
  

  while (!node->is_leaf)
  {
    /*
    if(!update_node_statistics(node, timeseries))
    {
	fprintf(stderr,"Error in hercules_index.c: could not update \
                        statistics at node %s\n", node->filename);
        return  FAILURE;
    }
    */
    
    if (node_split_policy_route_to_left(node,timeseries,&timeseries_segment_sketch))
      node = node->left_child;
    else
      node = node->right_child;
  }
 
            COUNT_PARTIAL_TIME_END
	    index->stats->idx_traverse_tree_total_time  += partial_time;	
	    index->stats->idx_traverse_tree_input_time  += partial_input_time;
	    index->stats->idx_traverse_tree_output_time += partial_output_time;
	    index->stats->idx_traverse_tree_cpu_time    += partial_time
	                                           - partial_input_time
	                                           - partial_output_time;
	    index->stats->idx_traverse_tree_seq_input_count   += partial_seq_input_count;
	    index->stats->idx_traverse_tree_seq_output_count  += partial_seq_output_count;
	    index->stats->idx_traverse_tree_rand_input_count  += partial_rand_input_count;
	    index->stats->idx_traverse_tree_rand_output_count += partial_rand_output_count;	    
 
	    RESET_PARTIAL_COUNTERS()
	    COUNT_PARTIAL_TIME_START
   

  if (node->is_leaf)
  {

    if(!update_node_statistics(node, timeseries))
    {
        fprintf(stderr,"Error in hercules_index.c: could not update \
                        statistics at node %s\n", node->filename);
        return FAILURE;
    }

    if(!append_ts_to_node(index,node, timeseries))
    {
        fprintf(stderr,"Error in hercules_index.c: could not append \
                        time series to node %s\n", node->filename);
        return FAILURE;
    }
    

            COUNT_PARTIAL_TIME_END
	    index->stats->idx_append_ts_to_leaf_total_time  += partial_time;	
	    index->stats->idx_append_ts_to_leaf_input_time  += partial_input_time;
	    index->stats->idx_append_ts_to_leaf_output_time += partial_output_time;
	    index->stats->idx_append_ts_to_leaf_cpu_time    += partial_time
	                                           - partial_input_time
	                                           - partial_output_time;
	    index->stats->idx_append_ts_to_leaf_seq_input_count  += partial_seq_input_count;
	    index->stats->idx_append_ts_to_leaf_seq_output_count += partial_seq_output_count;
	    index->stats->idx_append_ts_to_leaf_rand_input_count  += partial_rand_input_count;
	    index->stats->idx_append_ts_to_leaf_rand_output_count += partial_rand_output_count;	    
 
	    RESET_PARTIAL_COUNTERS()
	    COUNT_PARTIAL_TIME_START

    //if split needed, split the node and refresh curr_node
    if (node->node_size >= index->settings->max_leaf_size)
    {

      struct node_split_policy curr_node_split_policy;
      ts_type max_diff_value = (FLT_MAX * (-1));
      ts_type avg_children_range_value = 0;
      short hs_split_point = -1;
      short * child_node_points;
      int num_child_node_points = 0;
      const int num_child_segments = 2; //by default split to two subsegments
      //we want to test every possible split policy for each segment
      
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

          child_node_segment_sketches = malloc (sizeof(struct segment_sketch) *
					        num_child_segments);

          if(child_node_segment_sketches == NULL)
          {
             fprintf(stderr,"Error in hercules_index.c: could not allocate \
                            memory for the child node segment sketches for \
                            node  %s\n", node->filename);
             return FAILURE;
          }
    
	  for (int k = 0; k< num_child_segments;++k)
	  {
	    child_node_segment_sketches[k].indicators = NULL;
	    child_node_segment_sketches[k].indicators =malloc (sizeof(ts_type) *
						       curr_node_segment_sketch.num_indicators);
            if(child_node_segment_sketches[k].indicators == NULL)
            {
               fprintf(stderr,"Error in hercules_index.c: could not allocate\
                               memory for the child node segment sketches \
                               indicators for node  %s\n", node->filename);
               return FAILURE;
            }
	    
	  }


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
          child_node_segment_sketches = malloc (sizeof(struct segment_sketch) *
						num_child_segments);
          if(child_node_segment_sketches == NULL)
          {
             fprintf(stderr,"Error in hercules_index.c: could not allocate memory \
                            for the horizontal child node segment sketches for \
                            node  %s\n", node->filename);
             return FAILURE;
          }
	  
	  for (int k = 0; k< num_child_segments;++k)
	  {
	    child_node_segment_sketches[k].indicators = NULL;
	    child_node_segment_sketches[k].indicators = malloc (sizeof(ts_type) *curr_hs_node_segment_sketch.num_indicators);
            if(child_node_segment_sketches[k].indicators == NULL)
            {
               fprintf(stderr,"Error in hercules_index.c: could not allocate memory \
                             for the horizontal child node segment sketches indicators \
                             for node  %s\n", node->filename);
               return FAILURE;
            }	    
	  }

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

	    //printf("from = %hd, to = %hd, num_node_points = %d,num_hs_node_points = %d, horizontal segment  = %hd, hs_split_point = %hd, hs_index = %d\n", curr_node_split_policy.split_from, curr_node_split_policy.split_to , node->num_node_points, node->num_hs_node_points,split_segment, hs_split_point,i);
          }
	  
	  for (int k = 0; k< num_child_segments;++k)
	  {
	    free(child_node_segment_sketches[k].indicators);
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

            COUNT_PARTIAL_TIME_END
	    index->stats->idx_evaluate_split_policies_total_time  += partial_time;	
	    index->stats->idx_evaluate_split_policies_input_time  += partial_input_time;
	    index->stats->idx_evaluate_split_policies_output_time += partial_output_time;
	    index->stats->idx_evaluate_split_policies_cpu_time    += partial_time
	                                           - partial_input_time
	                                           - partial_output_time;
	    index->stats->idx_evaluate_split_policies_seq_input_count   += partial_seq_input_count;
	    index->stats->idx_evaluate_split_policies_seq_output_count  += partial_seq_output_count;
	    index->stats->idx_evaluate_split_policies_rand_input_count  += partial_rand_input_count;
	    index->stats->idx_evaluate_split_policies_rand_output_count += partial_rand_output_count;	    
 
	    RESET_PARTIAL_COUNTERS()
	    COUNT_PARTIAL_TIME_START



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
	if (hs_split_point < node->node_points[0])
	{
	  node->split_segment = 0;
	}
        for(int i = 0; i < (num_child_node_points-1);++i)
	{
	  child_node_points[i] = node->node_points[i];
	  //printf ("node_point[%d] = %hd, ", i, node->node_points[i]);
	  if (hs_split_point > node->node_points[i])
	  {
	    node->split_segment = i+1;
	  }
	}
	//printf ("\n");
	child_node_points[num_child_node_points-1] = hs_split_point; //initialize newly added point

	qsort(child_node_points,num_child_node_points,sizeof(short),compare_short);
      }
      //printf("split_segment  = %hd,  hs_split_point = %hd\n", split_segment,hs_split_point);      
      //this will put the time series of this node in the file_buffer->buffered_list aray
      //it will include the time series in disk and those in memory

      if(!split_node_create_children(index,node, child_node_points, num_child_node_points))
      {
           fprintf(stderr,"Error in hercules_index.c: could not split node %s.\n", node->filename);
           return FAILURE;
      }

      free(child_node_points);

     node->file_buffer->do_not_flush = true;

     if (!get_file_buffer(index, node))
      { 
           fprintf(stderr,"Error in hercules_index.c: could not get the file \
                           buffer for node %s.\n", node->filename);
           return FAILURE;              
      }

     //resetting the parent node size to 0, the actual size will be updated during index materialization to disk

      
      ts_type ** ts_list;
      ts_list = get_all_time_series_in_node(index, node,0);

      //copying the contents of the the node being split
      //in case it gets flushed from memory to disk

      //printf ("splitting node %s with size %d", node->filename, node->file_buffer->buffered_list_size);
      for (int idx=0; idx <index->settings->max_leaf_size;++idx)
      {
        if (node_split_policy_route_to_left(node, ts_list[idx],&timeseries_segment_sketch))
        {
  	   if(!update_node_statistics(node->left_child, ts_list[idx]))
           {
             fprintf(stderr,"Error in hercules_index.c: could not update \
                             statistics at left child of\
                             node %s\n", node->filename);
             return FAILURE;
           }

           if(!append_ts_to_child_node(index, node->left_child,ts_list[idx]))
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

           if(!append_ts_to_child_node(index, node->right_child,ts_list[idx]))
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

     //printf("splitted_node = %p right_child = %d left_child = %d\n", (void *)node, node->right_child->node_size,node->left_child->node_size);
      
      if(!delete_file_buffer(index,node))
      {
        fprintf(stderr,"Error in hercules_index.c: could not delete file buffer for \
                           node %s\n", node->filename);
        return FAILURE;
      }
            COUNT_PARTIAL_TIME_END
	    index->stats->idx_split_node_total_time  += partial_time;	
	    index->stats->idx_split_node_input_time  += partial_input_time;
	    index->stats->idx_split_node_output_time += partial_output_time;
	    index->stats->idx_split_node_cpu_time    += partial_time
	                                           - partial_input_time
	                                           - partial_output_time;
	    index->stats->idx_split_node_seq_input_count   += partial_seq_input_count;
	    index->stats->idx_split_node_seq_output_count  += partial_seq_output_count;
	    index->stats->idx_split_node_rand_input_count  += partial_rand_input_count;
	    index->stats->idx_split_node_rand_output_count += partial_rand_output_count;	    
 
	    RESET_PARTIAL_COUNTERS()
	    COUNT_PARTIAL_TIME_START

      
    } //end if_split_node 
  } //end if_node_is_leaf

  return SUCCESS;
}

enum response hercules_index_write_sequential(struct hercules_index *index)
{
	fprintf(stderr, ">>> Storing index: %s\n", index->settings->root_directory);
	const char *filename = malloc(sizeof(char) * (strlen(index->settings->root_directory) + 9));
	filename = strcpy(filename, index->settings->root_directory);
	filename = strcat(filename, "root.idx\0");

        COUNT_PARTIAL_RAND_OUTPUT
        COUNT_PARTIAL_OUTPUT_TIME_START	  
	FILE *file = fopen(filename, "wb");
        COUNT_PARTIAL_OUTPUT_TIME_END	  

	free(filename);

	if (index->settings->sims)
	{
	  const char *sims_filename = malloc(sizeof(char) * (strlen(index->settings->root_directory) + 19));
	  strcpy(sims_filename, index->settings->root_directory);
	  strcat(sims_filename, "leaves_sims.idx\0");
	  index->leaves_sims_file = fopen(sims_filename, "wb");
	  
	  index->leaves_sims_pos = 0;	  
	  free(sims_filename);
	}
	
	if (index->settings->serial)
	{
	  const char *leaves_filename = malloc(sizeof(char) * (strlen(index->settings->root_directory) + 18));
	  filename = strcpy(leaves_filename, index->settings->root_directory);
	  filename = strcat(leaves_filename, "leaves_raw.idx\0");
	  index->leaves_raw_file = fopen(leaves_filename, "wb");
	  
	  index->leaves_raw_pos = 0;
	  
	  free(leaves_filename);
	}

	
        if(file == NULL)
        {   
           fprintf(stderr, "Error in hercules_index.c: Could not open"
		" the index file. Reason = %s\n", strerror(errno));
           return FAILURE;
        }	
	
	unsigned int timeseries_size = index->settings->timeseries_size;
	unsigned int max_leaf_size = index->settings->max_leaf_size;
	unsigned int init_segments = index->settings->init_segments;
	double buffered_memory_size = index->settings->buffered_memory_size;
	int sims = index->settings->sims;
        int serial =  index->settings->serial;
	int paa_segments = index->settings->paa_segments;
	int sax_bit_cardinality = index->settings->sax_bit_cardinality;	
	
	// SETTINGS DATA
        COUNT_PARTIAL_SEQ_OUTPUT
        COUNT_PARTIAL_SEQ_OUTPUT
        COUNT_PARTIAL_SEQ_OUTPUT
        COUNT_PARTIAL_SEQ_OUTPUT	  
        COUNT_PARTIAL_OUTPUT_TIME_START
	fwrite(&leaf_nodes_count, sizeof(unsigned long), 1, file);	  
	fwrite(&buffered_memory_size, sizeof(double), 1, file);
	fwrite(&timeseries_size, sizeof(unsigned int), 1, file);
	fwrite(&init_segments, sizeof(unsigned int), 1, file);
	fwrite(&paa_segments, sizeof(int), 1, file);
	fwrite(&sax_bit_cardinality, sizeof(int), 1, file);			
	fwrite(&max_leaf_size, sizeof(unsigned int), 1, file);
	fwrite(&serial, sizeof(int), 1, file);
	fwrite(&sims, sizeof(int), 1, file);	
        COUNT_PARTIAL_OUTPUT_TIME_END


	// NODES AND FILE BUFFERS
	hercules_node_write(index, index->first_node, file);
        COUNT_PARTIAL_OUTPUT_TIME_START
        fseek(file, 0L, SEEK_SET);
        fwrite(&leaf_nodes_count, sizeof(unsigned long), 1, file);	  	  	  
	fclose(file);
        COUNT_PARTIAL_OUTPUT_TIME_END
	  
	return SUCCESS;
}

enum response hercules_index_write(struct hercules_index *index, int num_write_threads)
{
	fprintf(stderr, ">>> Storing index in parallel: %s\n", index->settings->root_directory);
	const char *filename = malloc(sizeof(char) * (strlen(index->settings->root_directory) + 9));
	filename = strcpy(filename, index->settings->root_directory);
	filename = strcat(filename, "root.idx\0");

        COUNT_PARTIAL_RAND_OUTPUT
        COUNT_PARTIAL_OUTPUT_TIME_START	  
	FILE *file = fopen(filename, "wb");
        COUNT_PARTIAL_OUTPUT_TIME_END	  

	free(filename);

	if (index->settings->sims)
	{
	  const char *sims_filename = malloc(sizeof(char) * (strlen(index->settings->root_directory) + 19));
	  strcpy(sims_filename, index->settings->root_directory);
	  strcat(sims_filename, "leaves_sims.idx\0");
	  index->leaves_sims_file = fopen(sims_filename, "wb");
	  
	  index->leaves_sims_pos = 0;
	  
	  free(sims_filename);
	}
	
	if (index->settings->serial)
	{
	  const char *leaves_filename = malloc(sizeof(char) * (strlen(index->settings->root_directory) + 18));
	  filename = strcpy(leaves_filename, index->settings->root_directory);
	  filename = strcat(leaves_filename, "leaves_raw.idx\0");
	  index->leaves_raw_file = fopen(leaves_filename, "wb");
	  
	  index->leaves_raw_pos = 0;
	  
	  free(leaves_filename);
	}

	
        if(file == NULL)
        {   
           fprintf(stderr, "Error in hercules_index.c: Could not open"
		" the index file. Reason = %s\n", strerror(errno));
           return FAILURE;
        }	
	
	unsigned int timeseries_size = index->settings->timeseries_size;
	unsigned int max_leaf_size = index->settings->max_leaf_size;
	unsigned int init_segments = index->settings->init_segments;
	double buffered_memory_size = index->settings->buffered_memory_size;
	int sims = index->settings->sims;
        int serial =  index->settings->serial;
	int paa_segments = index->settings->paa_segments;
	int sax_bit_cardinality = index->settings->sax_bit_cardinality;	
	
	// SETTINGS DATA
        COUNT_PARTIAL_SEQ_OUTPUT
        COUNT_PARTIAL_SEQ_OUTPUT
        COUNT_PARTIAL_SEQ_OUTPUT
        COUNT_PARTIAL_SEQ_OUTPUT	  
        COUNT_PARTIAL_OUTPUT_TIME_START
	fwrite(&leaf_nodes_count, sizeof(unsigned long), 1, file);	  
	fwrite(&buffered_memory_size, sizeof(double), 1, file);
	fwrite(&timeseries_size, sizeof(unsigned int), 1, file);
	fwrite(&init_segments, sizeof(unsigned int), 1, file);
	fwrite(&paa_segments, sizeof(int), 1, file);
	fwrite(&sax_bit_cardinality, sizeof(int), 1, file);			
	fwrite(&max_leaf_size, sizeof(unsigned int), 1, file);
	fwrite(&serial, sizeof(int), 1, file);
	fwrite(&sims, sizeof(int), 1, file);	
        COUNT_PARTIAL_OUTPUT_TIME_END


    //add the root as it was not counted
    ++index->stats->leaf_nodes_count;

    index->leaves = calloc (index->stats->leaf_nodes_count, sizeof(struct hercules_node *));

        //node->file_pos = index->leaves_raw_pos;

    get_leaves(index, index->first_node);

    if (index->settings->sims){
      hercules_index_flush_leaves(index, num_write_threads);
    }

	// NODES AND FILE BUFFERS
	hercules_index_node_write(index, index->first_node, file);

    //flush_leaf_to_leaves_file_update_stats_parallel(index,node, index->settings->sims);
    //reset leaves_counter to 0 for backward compatibility
    index->stats->leaves_counter = 0;

        COUNT_PARTIAL_OUTPUT_TIME_START
        fseek(file, 0L, SEEK_SET);
        fwrite(&leaf_nodes_count, sizeof(unsigned long), 1, file);	  	  	  
	fclose(file);
        COUNT_PARTIAL_OUTPUT_TIME_END
	  
	return SUCCESS;
}

struct hercules_index * hercules_index_read(const char* root_directory)
{
        if(chdir(root_directory) != 0)
        {
            fprintf(stderr, "The index directory does not exist. "
                             "Please provide a valid directory.\n");
            exit (-1);
        }

	fprintf(stderr, ">>> Loading index: %s\n", root_directory);

	const char *filename = malloc(sizeof(char) * (strlen(root_directory) + 9));
	filename = strcpy(filename, root_directory);
	filename = strcat(filename, "root.idx\0");
	
	//printf ("Index file: %s\n",filename);
        COUNT_PARTIAL_RAND_INPUT
        COUNT_PARTIAL_INPUT_TIME_START
	FILE *file = fopen(filename, "rb");
        COUNT_PARTIAL_INPUT_TIME_END

	free(filename);
	
        if(file == NULL)
        {   
           fprintf(stderr, "Error in hercules_index.c: Could not open"
		" the index file. Reason = %s\n", strerror(errno));
           return FAILURE;
        }

	unsigned long count_leaves = 0;	
	unsigned int timeseries_size = 0;
	unsigned int max_leaf_size = 0;
	unsigned int init_segments = 0;
        double buffered_memory_size = 0;
	boolean is_index_new = 0;
        int serial = 0;
	int sims = 0;
	int paa_segments = 0;
	int sax_bit_cardinality = 0;	
	
	
        COUNT_PARTIAL_SEQ_INPUT
        COUNT_PARTIAL_SEQ_INPUT
        COUNT_PARTIAL_SEQ_INPUT
        COUNT_PARTIAL_SEQ_INPUT	  
        COUNT_PARTIAL_INPUT_TIME_START
        fread(&count_leaves, sizeof(unsigned long), 1, file);
	fread(&buffered_memory_size, sizeof(double), 1, file);
	fread(&timeseries_size, sizeof(unsigned int), 1, file);
	fread(&init_segments, sizeof(unsigned int), 1, file);
	fread(&paa_segments, sizeof(int), 1, file);
	fread(&sax_bit_cardinality, sizeof(int), 1, file);
	fread(&max_leaf_size, sizeof(unsigned int), 1, file);
	fread(&serial, sizeof(int), 1, file);	
	fread(&sims, sizeof(int), 1, file);
        COUNT_PARTIAL_INPUT_TIME_END

	  printf ("serial is %d",serial);
        struct hercules_index_settings * index_settings =
	                       hercules_index_settings_init(root_directory,
							  timeseries_size,
							  init_segments,
							  paa_segments,
							  sax_bit_cardinality,							  
                                                          max_leaf_size, 
							  buffered_memory_size,
							  is_index_new,
							  serial,
							  sims,
							  1);
    
        struct hercules_index * index = hercules_index_init(index_settings);

	index->stats->leaves_heights = calloc(count_leaves, sizeof(int));
	index->stats->leaves_sizes = calloc(count_leaves, sizeof(int));
	index->stats->leaves_counter = 0;

        if(index_settings->serial)
	{
	  const char *leaves_filename = malloc(sizeof(char) * (strlen(index->settings->root_directory) + 18));
	  leaves_filename = strcpy(leaves_filename, index->settings->root_directory);
	  leaves_filename = strcat(leaves_filename, "leaves_raw.idx\0");
	  index->leaves_raw_filename = leaves_filename;

	  //index->leaves_raw_file = fopen(leaves_filename, "rb");
	  
	  index->leaves_raw_pos = 0;	 
	  //an array of pointers to all the leaves, needed by the serial scan
	  index->leaves = calloc (count_leaves, sizeof(struct hercules_node *));
	}

        if(index_settings->sims)
	{
	  const char *sims_filename = malloc(sizeof(char) * (strlen(index->settings->root_directory) + 19));
	  strcpy(sims_filename, index->settings->root_directory);
	  strcat(sims_filename, "leaves_sims.idx\0");
	  index->leaves_sims_filename = sims_filename;
	  index->leaves_sims_file = fopen(sims_filename, "rb");
	  
	  //index->leaves_raw_file = fopen(leaves_filename, "rb");
	  
	  index->leaves_sims_pos = 0;
	}
       
	
	index->first_node = hercules_node_read(index, file);
        COUNT_PARTIAL_INPUT_TIME_START	
	fclose(file);
        COUNT_PARTIAL_INPUT_TIME_END
	  
	return index;
}


enum response hercules_index_set_stats(struct hercules_index *index, struct hercules_node *node)
{

  unsigned int height = node->level +1;
  unsigned int threshold = index->settings->max_leaf_size;
  //at this point all time series are on disk
  unsigned int node_size = node->node_size;

  double node_fill_factor = (node_size * 100.0)/threshold;  
    
  if (node_fill_factor < index->stats->min_fill_factor)
  {
    index->stats->min_fill_factor = node_fill_factor;
  }
  if (node_fill_factor > index->stats->max_fill_factor)
  {
    index->stats->max_fill_factor = node_fill_factor;
  }
  if (height < index->stats->min_height)
  {
    index->stats->min_height = height;    
  }
  if (height > index->stats->max_height)
  {
    index->stats->max_height = height;
  }
  
  if(node_size == 0)
  {
    COUNT_EMPTY_LEAF_NODE
   }
  index->stats->sum_fill_factor += node_fill_factor ;
  index->stats->sum_squares_fill_factor += pow(node_fill_factor,2) ;

  index->stats->sum_height += height;
  index->stats->sum_squares_height += pow(height,2) ;
     
  COUNT_TOTAL_TS(node_size)
  
}
  
enum response get_leaves(struct hercules_index *index, struct hercules_node *node)
{
	if(!node->is_leaf)
	{
		get_leaves(index, node->left_child);
		get_leaves(index, node->right_child);
    }
    else
    {
       node->file_pos = index->leaves_raw_pos;
	   index->leaves_raw_pos = index->leaves_raw_pos + node->node_size;  
 	   index->leaves[index->stats->leaves_counter] = node;
	   ++(index->stats->leaves_counter);
    }			


}

enum response hercules_node_write(struct hercules_index *index, struct hercules_node *node, FILE *file)
{

        COUNT_PARTIAL_SEQ_OUTPUT
	  //COUNT_PARTIAL_SEQ_OUTPUT	
        COUNT_PARTIAL_SEQ_OUTPUT  
	COUNT_PARTIAL_OUTPUT_TIME_START	  
        fwrite(&(node->is_leaf), sizeof(unsigned char), 1, file);
	//fwrite(&(node->node_size), sizeof(unsigned int), 1, file);
        fwrite(&(node->level), sizeof(unsigned int), 1, file);
        COUNT_PARTIAL_OUTPUT_TIME_END

	int idx = -1;
	short  split_segment;

	if(!node->is_leaf)
	{
		hercules_node_write(index, node->left_child, file);
		hercules_node_write(index, node->right_child, file);
		
		
		//update the parent's node_size
		
                COUNT_PARTIAL_SEQ_OUTPUT
                COUNT_PARTIAL_SEQ_OUTPUT
                COUNT_PARTIAL_SEQ_OUTPUT
		COUNT_PARTIAL_OUTPUT_TIME_START
		fwrite(&(node->node_size), sizeof(unsigned int), 1, file);		  
		fwrite(&(node->split_segment), sizeof(short), 1, file);		  
		fwrite(node->split_policy, sizeof(struct node_split_policy), 1, file);
		fwrite(&(node->num_node_points), sizeof(short), 1, file);
		fwrite(node->node_points, sizeof(short), node->num_node_points, file);
		//printf("node_size = %d, node_split_segment = %hd, left=%d, level=%d\n", node->node_size, node->split_segment, node->is_left, node->level);
		split_segment = node->split_segment;
                for (int i =0; i< node->num_node_points; ++i)
		{
		  //update the non-split segment sketches before writing to disk
		  /*
		  if (i != split_segment)
		  {
		      if((i < split_segment) | (split_segment == -1))	
			idx = i;
		      else if (i > split_segment)
			idx = i+1;
		      
		      node->node_segment_sketches[i].indicators[0] = fmaxf(node->left_child->node_segment_sketches[idx].indicators[0], node->right_child->node_segment_sketches[idx].indicators[0]);
		      node->node_segment_sketches[i].indicators[1] = fminf(node->left_child->node_segment_sketches[idx].indicators[1], node->right_child->node_segment_sketches[idx].indicators[1]);
		      node->node_segment_sketches[i].indicators[2] = fmaxf(node->left_child->node_segment_sketches[idx].indicators[2], node->right_child->node_segment_sketches[idx].indicators[2]);
		      node->node_segment_sketches[i].indicators[3] = fminf(node->left_child->node_segment_sketches[idx].indicators[3], node->right_child->node_segment_sketches[idx].indicators[3]);		
		  }
		  */
		    
                   COUNT_PARTIAL_SEQ_OUTPUT
                   COUNT_PARTIAL_SEQ_OUTPUT		     
                   fwrite(&(node->node_segment_sketches[i].num_indicators),
			  sizeof(int),
			  1,
			  file);
  	           fwrite(node->node_segment_sketches[i].indicators,
			  sizeof(ts_type),
			  node->node_segment_sketches[i].num_indicators,
			  file);
		   /*
		   printf("max_mean = %g, min_mean = %g, max_sd=%g, min_sd=%g\t",
			  node->node_segment_sketches[i].indicators[0],
			  node->node_segment_sketches[i].indicators[1],
			  node->node_segment_sketches[i].indicators[2],
			  node->node_segment_sketches[i].indicators[3]);			  
		   */			  
		}

		//printf("\n\n");
		COUNT_PARTIAL_OUTPUT_TIME_END
		
	}
	else
	{
	  fwrite(&(node->node_size), sizeof(unsigned int), 1, file);		  
	  if(node->filename != NULL) {
			int filename_size = strlen(node->filename);
                        COUNT_PARTIAL_SEQ_OUTPUT	
                        COUNT_PARTIAL_SEQ_OUTPUT  
                  	COUNT_PARTIAL_OUTPUT_TIME_START
               	        fwrite(&filename_size, sizeof(int), 1, file);
			fwrite(node->filename, sizeof(char), filename_size, file);

                        COUNT_PARTIAL_SEQ_OUTPUT
                        COUNT_PARTIAL_SEQ_OUTPUT
 		        fwrite(&(node->num_node_points), sizeof(short), 1, file);
		        fwrite(node->node_points, sizeof(short), node->num_node_points, file);
			//printf("node_size = %d, node_split_segment = %hd, left=%d, level=%d (leaf)\n", node->node_size, node->split_segment, node->is_left, node->level);
                        for (int i =0; i< node->num_node_points; ++i)
	        	{
                            COUNT_PARTIAL_SEQ_OUTPUT
                            COUNT_PARTIAL_SEQ_OUTPUT		     
                           fwrite(&(node->node_segment_sketches[i].num_indicators),
	  		           sizeof(int),
			           1,
			            file);
       	                   fwrite(node->node_segment_sketches[i].indicators,
			         sizeof(ts_type),
			         node->node_segment_sketches[i].num_indicators,
			         file);
			   /*
			   printf("max_mean = %g, min_mean = %g, max_sd=%g, min_sd=%g\t",
				  node->node_segment_sketches[i].indicators[0],
				  node->node_segment_sketches[i].indicators[1],
				  node->node_segment_sketches[i].indicators[2],
				  node->node_segment_sketches[i].indicators[3]);			  			  
			   */

    		        }

			//printf("\n\n");
		
                         COUNT_PARTIAL_OUTPUT_TIME_END
		

			//flush_buffer_to_disk(index,node);

			 if (index->settings->serial)
			 {
			   node->file_pos = index->leaves_raw_pos;
			   fwrite(&(node->file_pos), sizeof(unsigned long long), 1, file);
			   
			   //if (index->settings->is_new)
			   flush_leaf_to_leaves_file_update_stats_serial(index,node, index->settings->sims);
               //The leaves_raw_pos now means the leaf id so we can reuse it for the sax file
			   //++index->leaves_raw_pos;  
         
			   //index->leaves_raw_pos = index->leaves_raw_pos + node->node_size;  
			 }
			 
                        COUNT_LEAF_NODE
                         //collect stats while traversing the index
	    int full_size = strlen(index->settings->root_directory) + strlen(node->filename)+1;
	
	    const char *full_filename = malloc(sizeof(char) * full_size);
	    full_filename = strcpy(full_filename, index->settings->root_directory);
	    full_filename = strcat(full_filename, node->filename);
	    full_filename = strcat(full_filename, "\0");
        remove(full_filename);
        free(full_filename);
			hercules_index_set_stats(index, node);

		}
		else {
			int filename_size = 0;
                        COUNT_PARTIAL_SEQ_OUTPUT  
                  	COUNT_PARTIAL_OUTPUT_TIME_START
			fwrite(&filename_size, sizeof(int), 1, file);
                        COUNT_PARTIAL_OUTPUT_TIME_END						
		}
	  
	}

	
	return SUCCESS;
}

enum response hercules_index_node_write(struct hercules_index *index, struct hercules_node *node, FILE *file)
{

        COUNT_PARTIAL_SEQ_OUTPUT
	  //COUNT_PARTIAL_SEQ_OUTPUT	
        COUNT_PARTIAL_SEQ_OUTPUT  
	COUNT_PARTIAL_OUTPUT_TIME_START	  
        fwrite(&(node->is_leaf), sizeof(unsigned char), 1, file);
	//fwrite(&(node->node_size), sizeof(unsigned int), 1, file);
        fwrite(&(node->level), sizeof(unsigned int), 1, file);
        COUNT_PARTIAL_OUTPUT_TIME_END

	int idx = -1;
	short  split_segment;

	if(!node->is_leaf)
	{
		hercules_index_node_write(index, node->left_child, file);
		hercules_index_node_write(index, node->right_child, file);
		
		node->node_size = node->left_child->node_size + node->right_child->node_size;
		
		//update the parent's node_size
		
                COUNT_PARTIAL_SEQ_OUTPUT
                COUNT_PARTIAL_SEQ_OUTPUT
                COUNT_PARTIAL_SEQ_OUTPUT
		COUNT_PARTIAL_OUTPUT_TIME_START
		fwrite(&(node->node_size), sizeof(unsigned int), 1, file);		  
		fwrite(&(node->split_segment), sizeof(short), 1, file);		  
		fwrite(node->split_policy, sizeof(struct node_split_policy), 1, file);
		fwrite(&(node->num_node_points), sizeof(short), 1, file);
		fwrite(node->node_points, sizeof(short), node->num_node_points, file);
		//printf("node_size = %d, node_split_segment = %hd, left=%d, level=%d\n", node->node_size, node->split_segment, node->is_left, node->level);
		split_segment = node->split_segment;
                for (int i =0; i< node->num_node_points; ++i)
		{
		    
                   COUNT_PARTIAL_SEQ_OUTPUT
                   COUNT_PARTIAL_SEQ_OUTPUT		     
                   fwrite(&(node->node_segment_sketches[i].num_indicators),
			  sizeof(int),
			  1,
			  file);
  	           fwrite(node->node_segment_sketches[i].indicators,
			  sizeof(ts_type),
			  node->node_segment_sketches[i].num_indicators,
			  file);
		}

		//printf("\n\n");
		COUNT_PARTIAL_OUTPUT_TIME_END
		
	}
	else
	{
	  fwrite(&(node->node_size), sizeof(unsigned int), 1, file);		  
	  if(node->filename != NULL) {
			int filename_size = strlen(node->filename);
                        COUNT_PARTIAL_SEQ_OUTPUT	
                        COUNT_PARTIAL_SEQ_OUTPUT  
                  	COUNT_PARTIAL_OUTPUT_TIME_START
               	        fwrite(&filename_size, sizeof(int), 1, file);
			fwrite(node->filename, sizeof(char), filename_size, file);

                        COUNT_PARTIAL_SEQ_OUTPUT
                        COUNT_PARTIAL_SEQ_OUTPUT
 		        fwrite(&(node->num_node_points), sizeof(short), 1, file);
		        fwrite(node->node_points, sizeof(short), node->num_node_points, file);
			//printf("node_size = %d, node_split_segment = %hd, left=%d, level=%d (leaf)\n", node->node_size, node->split_segment, node->is_left, node->level);
                        for (int i =0; i< node->num_node_points; ++i)
	        	{
                            COUNT_PARTIAL_SEQ_OUTPUT
                            COUNT_PARTIAL_SEQ_OUTPUT		     
                           fwrite(&(node->node_segment_sketches[i].num_indicators),
	  		           sizeof(int),
			           1,
			            file);
       	                   fwrite(node->node_segment_sketches[i].indicators,
			         sizeof(ts_type),
			         node->node_segment_sketches[i].num_indicators,
			         file);

    		        }

			//printf("\n\n");
		
                         COUNT_PARTIAL_OUTPUT_TIME_END
		
			//flush_buffer_to_disk(index,node);

			 if (index->settings->serial)
			 {
			   //node->file_pos = index->leaves_raw_pos;

			   fwrite(&(node->file_pos), sizeof(unsigned long long), 1, file);
	           //index->leaves_raw_pos = index->leaves_raw_pos + node->node_size;  
 			   //index->leaves[index->stats->leaves_counter] = node;
			   //++(index->stats->leaves_counter);
			 }
			 
                        COUNT_LEAF_NODE
                         //collect stats while traversing the index
	
			hercules_index_set_stats(index, node);

		}
		else {
			int filename_size = 0;
                        COUNT_PARTIAL_SEQ_OUTPUT  
                  	COUNT_PARTIAL_OUTPUT_TIME_START
			fwrite(&filename_size, sizeof(int), 1, file);
                        COUNT_PARTIAL_OUTPUT_TIME_END						
		}
	  
	}

	
	return SUCCESS;
}

//IMPORTANT!!! remember to change THIS to adapt to the new format

struct hercules_node * hercules_node_read(struct hercules_index *index, FILE *file) {

        struct hercules_node *node = NULL;
        //to initialize node values for leaf and internal nodes
	node = hercules_leaf_node_init(sizeof(struct hercules_node));
	
	unsigned char is_leaf = 0;
        COUNT_PARTIAL_SEQ_INPUT
	  //COUNT_PARTIAL_SEQ_INPUT		     
        COUNT_PARTIAL_SEQ_INPUT
	COUNT_PARTIAL_INPUT_TIME_START	
	fread(&is_leaf, sizeof(unsigned char), 1, file);
	//fread(&(node->node_size), sizeof(unsigned int), 1, file);
        fread(&(node->level), sizeof(unsigned int), 1, file);
        COUNT_PARTIAL_INPUT_TIME_END		
	if(!is_leaf) {
		node->left_child = hercules_node_read(index, file);
		node->right_child = hercules_node_read(index, file);
                node->split_policy = NULL;
		node->split_policy = malloc(sizeof(struct node_split_policy));
                COUNT_PARTIAL_SEQ_INPUT	
                COUNT_PARTIAL_SEQ_INPUT
	        COUNT_PARTIAL_INPUT_TIME_START
		fread(&(node->node_size), sizeof(unsigned int), 1, file);		  
		fread(&(node->split_segment), sizeof(short), 1, file);		  
	        fread(node->split_policy, sizeof(struct node_split_policy), 1, file);
		fread(&(node->num_node_points), sizeof(short), 1, file);
	        COUNT_PARTIAL_INPUT_TIME_END
		  
                node->node_points = NULL;
		node->node_points = malloc(sizeof(short) * node->num_node_points);

                COUNT_PARTIAL_SEQ_INPUT		  		
	        COUNT_PARTIAL_INPUT_TIME_START
		fread(node->node_points, sizeof(short), node->num_node_points, file);
	        COUNT_PARTIAL_INPUT_TIME_END
		  
                node->node_segment_sketches = NULL;                           
  		node->node_segment_sketches = malloc(sizeof(struct segment_sketch)*
						     node->num_node_points);
                for (int i =0; i< node->num_node_points; ++i)
		{
                    COUNT_PARTIAL_SEQ_INPUT		  		
   	            COUNT_PARTIAL_INPUT_TIME_START		  
		    fread(&(node->node_segment_sketches[i].num_indicators),
			  sizeof(int),
			  1,
			  file);
   	            COUNT_PARTIAL_INPUT_TIME_END		  		    

		    node->node_segment_sketches[i].indicators = NULL;
		    node->node_segment_sketches[i].indicators = 
		      malloc(sizeof(ts_type) * node->node_segment_sketches[i].num_indicators );
                    COUNT_PARTIAL_SEQ_INPUT		  		
   	            COUNT_PARTIAL_INPUT_TIME_START		  		    
		    fread(node->node_segment_sketches[i].indicators,
			  sizeof(ts_type),
			  node->node_segment_sketches[i].num_indicators,
			  file);
   	            COUNT_PARTIAL_INPUT_TIME_END		  		    		    
		}
		node->file_buffer = NULL;
		node->is_leaf = 0;
		node->filename = NULL;
		
	}
	else {
		node->is_leaf = 1;
		hercules_file_buffer_init(node);

		int filename_size = 0;

                COUNT_PARTIAL_SEQ_INPUT
	        COUNT_PARTIAL_INPUT_TIME_START
		fread(&(node->node_size), sizeof(unsigned int), 1, file);		  		  
		fread(&filename_size, sizeof(int), 1, file);
                COUNT_PARTIAL_INPUT_TIME_END				

		node->file_buffer->disk_count = node->node_size;

		if(filename_size > 0)
		{
			node->filename = malloc(sizeof(char) * (filename_size + 1));
			
                        COUNT_PARTIAL_SEQ_INPUT
	                COUNT_PARTIAL_INPUT_TIME_START
			fread(node->filename, sizeof(char), filename_size, file);
                        COUNT_PARTIAL_INPUT_TIME_END							

			node->filename[filename_size] = '\0';

                        COUNT_PARTIAL_SEQ_INPUT
                        COUNT_PARTIAL_SEQ_INPUT
	                COUNT_PARTIAL_INPUT_TIME_START			  
 		        fread(&(node->num_node_points), sizeof(short), 1, file);
	                COUNT_PARTIAL_INPUT_TIME_END
			node->node_points = NULL;
         		node->node_points = malloc(sizeof(short) * node->num_node_points);
	                COUNT_PARTIAL_INPUT_TIME_START			  			
		        fread(node->node_points, sizeof(short), node->num_node_points, file);
	                COUNT_PARTIAL_INPUT_TIME_END

                        node->node_segment_sketches = NULL;                           
          		node->node_segment_sketches = malloc(sizeof(struct segment_sketch)*
						     node->num_node_points);
		
                        for (int i =0; i< node->num_node_points; ++i)
	        	{
                            COUNT_PARTIAL_SEQ_INPUT
                            COUNT_PARTIAL_SEQ_INPUT
        	            COUNT_PARTIAL_INPUT_TIME_START			      
                            fread(&(node->node_segment_sketches[i].num_indicators),
	  		           sizeof(int),
			           1,
			            file);
   	                    COUNT_PARTIAL_INPUT_TIME_END
		           node->node_segment_sketches[i].indicators = NULL;
        	           node->node_segment_sketches[i].indicators = 
		           malloc(sizeof(ts_type) * node->node_segment_sketches[i].num_indicators );

			   COUNT_PARTIAL_INPUT_TIME_START			      

       	                   fread(node->node_segment_sketches[i].indicators,
			         sizeof(ts_type),
			         node->node_segment_sketches[i].num_indicators,
			         file);
        	            COUNT_PARTIAL_INPUT_TIME_END			      
			   

    		        }

			COUNT_LEAF_NODE
			  
   	                index->stats->leaves_heights[index->stats->leaves_counter] = node->level + 1;
			index->stats->leaves_sizes[index->stats->leaves_counter] = node->node_size;

			if (index->settings->serial)
			  {
			    fread(&(node->file_pos), sizeof(unsigned long long), 1, file);
			    
			    index->leaves[index->stats->leaves_counter] = node;
			  }
			++(index->stats->leaves_counter);
			
			hercules_index_set_stats(index, node);			  
                 }
		else
		{
			node->filename = NULL;
			node->node_size = 0;
		}
		//get all timeseries for this node into its this file buffer

		//get_all_time_series_in_node(index,node);
		//node->file_buffer->disk_count = 0;
	}


	return node;
}

void cache_sax_file(struct hercules_index *index) {
    fseek(index->leaves_sims_file, 0L, SEEK_END);
    unsigned long size = ftell(index->leaves_sims_file);
    fseek(index->leaves_sims_file, 0L, SEEK_SET);
    index->sax_cache = malloc(size);
    index->sax_cache_size = size / index->settings->sax_byte_size;
    COUNT_INPUT_TIME_START
    fread(index->sax_cache, index->settings->sax_byte_size, index->sax_cache_size, index->leaves_sims_file);
    COUNT_INPUT_TIME_END
}


enum response hercules_index_build(const char *ifilename, file_position_type data_size,
					 struct hercules_index *index, int num_threads, unsigned int initial_db_size)

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
    if (total_records < data_size) {
        fprintf(stderr, "File %s has only %llu records!\n", ifilename, total_records);
        return FAILURE;
    }
    	
    file_position_type ts_loaded = 0;    
    
    int percentage = 100;
    if (percentage > 100)
    {
       percentage = (int) (data_size / (file_position_type) 100);
    }

    index_thread_data *input_data=malloc(sizeof(index_thread_data)*(num_threads-1));
    pthread_t threadid[num_threads-1];

    int j = 0;
    unsigned long long int i = 0;
    
    pthread_barrier_t db_barrier, continue_barrier, flush_barrier; //barrier1 is the double buffer barrier and barrier2 is the flush barrier
    pthread_barrier_init(&db_barrier, NULL, num_threads); //dbBarrier
    pthread_barrier_init(&continue_barrier, NULL, num_threads-1); //flushBarrier
    pthread_barrier_init(&flush_barrier, NULL, num_threads-1);
    

    ts_type ** double_buffer;

	unsigned int db_size[2]; 
	//unsigned int * db_counter;
	unsigned int db_counter[2];
	bool  finished[2];
    unsigned int flush_counter;
	unsigned int flush_order;

    static int toggle = 0;

    db_size[0] = 0;
    db_size[1] = 0;
    db_counter[0] = 0;

    db_counter[0] = 0;
    db_counter[1] = 0;
	//db_counter = NULL;

	double_buffer = NULL;

	finished[0] = false;
	finished[1] = false;
	flush_counter = 0;
	flush_order = 0;


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
    bool handshake = false;
    volatile int cnt = 0;
    int flush_size = 0;
    //volatile int flush_order = 0;
    	 
	db_size[toggle] = (unsigned int) ull_min((unsigned long long)initial_db_size, data_size);

    //This will hold both regions of the double buffer
    //double_buffer     = malloc(sizeof(ts_type) * index->settings->timeseries_size*initial_db_size*2);
	double_buffer     = calloc(2,sizeof(ts_type*));
    double_buffer[0]  = malloc(sizeof(ts_type) * index->settings->timeseries_size*initial_db_size);
    double_buffer[1]  = malloc(sizeof(ts_type) * index->settings->timeseries_size*initial_db_size);

    //double_buffer[1-toggle]   = malloc(sizeof(ts_type) * index->settings->timeseries_size*initial_db_size);

    //db_counter     = calloc(2,sizeof(unsigned int));

    
    for (i = 0; i < (num_threads-1); i++)
    {
		input_data[i].threads_data = input_data;

        input_data[i].index=index;

        input_data[i].buffer_counter=0;
        input_data[i].buffer_max=t_buffer_size-max_leaf_size -db_size[toggle];
        input_data[i].buffer_offset= index->buffer_manager->current_record + (sizeof(ts_type) * ts_length) * t_buffer_size * i ;

        input_data[i].continue_handshake=0;
		input_data[i].num_insert_workers = num_threads-1;
		input_data[i].thread_id = i;

 		input_data[i].is_flusher = false;
        input_data[i].flush_order=&flush_order;
        input_data[i].flush_counter=&flush_counter;
 		input_data[i].db_counter = db_counter;
 		input_data[i].db_size = db_size;
 		input_data[i].double_buffer = double_buffer;
        input_data[i].finished=finished;
        input_data[i].initial_db_size=initial_db_size;

        input_data[i].db_barrier=&db_barrier;
        input_data[i].continue_barrier=&continue_barrier;
        input_data[i].flush_barrier=&flush_barrier;

	

     	input_data[i].split_node_data = calloc(max_leaf_size, sizeof(ts_type *));
		for (int j =0; j < max_leaf_size ; ++j)
		{
		  input_data[i].split_node_data[j] =  calloc(ts_length, sizeof(ts_type));
		}

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
 
    }

    input_data[0].is_flusher = true;

    ts_type *timeseries;
    printf("read next %llu series\n", db_size[toggle]);	

    //timeseries = double_buffer+(toggle*initial_db_size);
	timeseries = double_buffer[toggle];
    fread(double_buffer[toggle], sizeof(ts_type), index->settings->timeseries_size*db_size[toggle], ifile);
   	toggle = 1-toggle;

	for (j = 0; j < (num_threads-1); j++)
	  {
            pthread_create(&(threadid[j]),NULL,hercules_index_insert_worker,(void*)&(input_data[j]));
		    //printf("created thread %d \n", j);	
	  }

//	for (i = db_size[1-toggle]; i <= data_size && (db_size[1-toggle] !=0); i+=db_size[1-toggle])
	for (i = db_size[1-toggle]; i < data_size; i+=db_size[toggle])
	  {
		    printf("read next %llu series\n", i);	
			db_size[toggle] = (unsigned int) ull_min((unsigned long long)initial_db_size, data_size-i);
		    //printf("db_size[toggle] = %llu\n", db_size[toggle]);	

		    //timeseries = double_buffer+(toggle*initial_db_size);
			timeseries = double_buffer[toggle];
    		fread(timeseries, sizeof(ts_type), index->settings->timeseries_size*db_size[toggle], ifile);

		    //fread(double_buffer[toggle], sizeof(ts_type), index->settings->timeseries_size*db_size[toggle], ifile);
			db_counter [toggle] = 0;

			toggle = 1 - toggle;
		    //printf("coordinator releasing db_barrier\n");	
			pthread_barrier_wait(&db_barrier);
		    //printf("coordinator released db_barrier\n");	
		    //printf("after barrier i =  %llu, toggle = %d \n",  i, toggle);	
		    //printf("after barrier db_size[toggle] = %llu\n", db_size[toggle]);	
	  }

    for (j = 0; j < (num_threads-1); j++)
	  {
            input_data[j].finished[toggle]=true;
	  }

	pthread_barrier_wait(&db_barrier);

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

			for (int i = 0 ; i < max_leaf_size; ++i)
			  {
			    free(input_data[j].split_node_data[i]);
			  }
			free(input_data[j].split_node_data);
      }	

    printf ("Finished\n");
    free(double_buffer[0]);
    free(double_buffer[1]);
    free(double_buffer);
    free(input_data);

    pthread_barrier_destroy(&db_barrier);
    pthread_barrier_destroy(&continue_barrier);
    pthread_barrier_destroy(&flush_barrier);

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

void* hercules_index_insert_worker(void *transferdata)
{
	int toggle = 0;
	int pos = 0;
    // int fin_number=((index_thread_data*)transferdata)->fin_number;
    struct hercules_index *index= ((index_thread_data*)transferdata)->index;
    //int i=0,j;
    int i,j,k;
    int num_threads = ((index_thread_data*)transferdata)->num_insert_workers;
    
    int volatile cnt = 0;
 
    unsigned int ts_id = 0;
    ts_type * ts = NULL;
    //file_position_type *offset = malloc(sizeof(file_position_type));
    file_position_type offset;

    int buffer_counter =  ((index_thread_data*)transferdata)->buffer_counter;
    unsigned int initial_db_size =  ((index_thread_data*)transferdata)->initial_db_size;
	

     double tS = 0;
     double tE = 0;
	 ts_type * timeseries;
	 struct timeval start_time;
     struct timeval end_time;

	    #if DETAILED_STATS == 1
         gettimeofday(&start_time, NULL);
	    #endif


    while(!((index_thread_data*)transferdata)->finished[toggle])
    {   

       if ((((index_thread_data*)transferdata)->buffer_counter) < (((index_thread_data*)transferdata)->buffer_max))
       {
	        //ts_id = __sync_fetch_and_add((((index_thread_data*)transferdata)->db_counter+toggle),1);
			//db_counter =  ((index_thread_data*)transferdata)->db_counter[toggle];
	        //ts_id = __sync_fetch_and_add(&db_counter,1);
			ts_id = __sync_fetch_and_add( &(((index_thread_data*)transferdata)->db_counter[toggle]),1);

   	        //printf("thread = %d has db_size[%d] =  %u \n",((index_thread_data*)transferdata)->thread_id,toggle, ((index_thread_data*)transferdata)->db_size[toggle]);	

		   while(ts_id < ((index_thread_data*)transferdata)->db_size[toggle])   
			{
  			   //offset =  (toggle*initial_db_size)+ts_id*index->settings->timeseries_size ;
				offset =  ts_id*index->settings->timeseries_size ;

 			    timeseries = (((ts_type*) ((index_thread_data*)transferdata)->double_buffer[toggle]) + offset);
 			   //printf("thread = %d inserting series # %u, series[0] = %g \n",((index_thread_data*)transferdata)->thread_id, ts_id,timeseries[0]);	

	           hercules_index_insert_series_to_node_parallel(index, index->first_node,
									    	//(((ts_type*)(((index_thread_data*)transferdata)->double_buffer[toggle])+offset)), 
											timeseries,
                                            (void*)transferdata);
    
			  //printf("thread = %d inserted series # %u \n",((index_thread_data*)transferdata)->thread_id, ts_id);	
  	    	  //db_counter =   ((index_thread_data*)transferdata)->db_counter[toggle];
			  //ts_id = __sync_fetch_and_add(&db_counter,1);

			   ts_id = __sync_fetch_and_add( &(((index_thread_data*)transferdata)->db_counter[toggle]),1);
			}   
      }
      else{
          //printf("thread = %lu, buffer_counter  = %d, buffer_max  = %d \n", 
		  //	pthread_self(), 
		  //	((index_thread_data*)transferdata)->buffer_counter,
		  //	((index_thread_data*)transferdata)->buffer_max
          //);
      }
      pthread_barrier_wait(((index_thread_data*)transferdata)->db_barrier);
	  //printf("thread = %d released db barrier\n",((index_thread_data*)transferdata)->thread_id);	     	    
 
      if (((index_thread_data*)transferdata)->is_flusher)
      {
		hercules_index_flush_coordinator(transferdata);
      }
      else
      {
		hercules_index_flush_worker(transferdata);	
      }
	  
      toggle = 1 - toggle;

    }
    //free(offset) ;
}



void hercules_index_flush_coordinator(void *transferdata)
{
  int flush_size = 0;
  struct hercules_index * index = ((index_thread_data*)transferdata)->index;
  //index_thread_data* input_data = ((index_thread_data*)transferdata)->threads_data;
  int num_threads = ((index_thread_data*)transferdata)->num_insert_workers;
  int j = 0;
  volatile int cnt = 0;
  int s;

  double tS = 0;
  double tE = 0;
  struct timeval start_time;
  struct timeval start_time2;

  struct timeval end_time;
  struct timeval end_time2;
  int worker_id = ((index_thread_data*)transferdata)->thread_id;
  #if DETAILED_STATS == 1 
  gettimeofday(&start_time2, NULL);
  #endif

  int buffer_counter = 0;

  int num_full_buffers = 0 ;

  __sync_fetch_and_add(&(((index_thread_data*)transferdata)->threads_data[worker_id].continue_handshake),1);
  
  for ( j = 0; j < num_threads; j++) 
    {   
      flush_size += (((index_thread_data*)transferdata)->threads_data[j].buffer_counter);
      //check that all workers have received the order
      while (  __sync_fetch_and_add(&(((index_thread_data*)transferdata)->threads_data[j].continue_handshake),0) != 1) 
  	  {
	  //for (int tmp = 0; tmp < BUSYWAIT_THRESHOLD; tmp++) cnt++;
	  //for (int tmp = 0; tmp < 1000; tmp++) cnt++;
	  //printf("thread = %lu  (coordinator) busy waiting\n", pthread_self());
	  }
      //printf("Flusher received handshake from thread %d\n", j);		  
    }

  //printf("Flusher checked all handshakes\n");
   
  if ( (((__sync_fetch_and_add(((index_thread_data*)transferdata)->flush_counter,0))) >= index->settings->flush_threshold) ||
         ((((index_thread_data*)transferdata)->buffer_counter) >= (((index_thread_data*)transferdata)->buffer_max))
       )
  {
    __sync_fetch_and_add(((index_thread_data*)transferdata)->flush_order,1);
    //printf("Flusher sent order\n");	
  }
  
  num_full_buffers = *(((index_thread_data*)transferdata)->flush_counter);
  __sync_fetch_and_sub(((index_thread_data*)transferdata)->flush_counter,num_full_buffers);
  //printf("Flusher reset flush_counter to %d (should be zero) \n", *(((index_thread_data*)transferdata)->flush_counter));		  

  s = pthread_barrier_wait(((index_thread_data*)transferdata)->continue_barrier);
  //printf("thread = %lu  (coordinator) released continue_barrier \n", pthread_self());
  __sync_fetch_and_sub(&(((index_thread_data*)transferdata)->threads_data[worker_id].continue_handshake),1);
    
  if (__sync_fetch_and_add(((index_thread_data*)transferdata)->flush_order,0)  != 0)
  {
    char * curr_time;
    curr_time= NULL;
    curr_time = malloc (sizeof(char) *26);
    get_current_time(curr_time);
    //printf ("%s, batch remove ! %lu  \n", curr_time, flush_size);
	printf ("%s, flushing started ! %lu  \n", curr_time);
    free(curr_time);

#if DETAILED_STATS == 1 
    gettimeofday(&start_time, NULL);
#endif

    flush_index_to_disk(index,index->first_node);
 
#if DETAILED_STATS == 1
    gettimeofday(&end_time, NULL);
    tS = start_time.tv_sec*1000000 + (start_time.tv_usec);
    tE = end_time.tv_sec*1000000 + (end_time.tv_usec);
    ((index_thread_data*)transferdata)->thread_flush_output_time += (tE - tS);
#endif
  
    memset(index->buffer_manager->mem_array,0,index->buffer_manager->max_record_index);
    index->buffer_manager->current_record_index = 0;
    index->buffer_manager->current_record = index->buffer_manager->mem_array;
    
    __sync_fetch_and_sub((((index_thread_data*)transferdata)->flush_order),1);            	
    //printf("Flusher reset flush_order to FALSE \n");

    s = pthread_barrier_wait(((index_thread_data*)transferdata)->flush_barrier);
    //printf("thread = %lu  (coordinator) released flush_barrier, number of full buffers = %d, flush order = %d \n", pthread_self(),
    //*(((index_thread_data*)transferdata)->flush_counter), *(((index_thread_data*)transferdata)->flush_order) );    

    buffer_counter = ((index_thread_data*)transferdata)->buffer_counter;
    __sync_fetch_and_sub(& (((index_thread_data*)transferdata)->buffer_counter) , buffer_counter);

    printf ("%s, flushing completed ! %lu  \n", curr_time);

    //for ( j = 0; j < num_threads; j++) 
    //{   
     //printf("buffer_counter[%d] = %d \n", j, ((index_thread_data*)transferdata)->threads_data[j].buffer_counter );
    //}
  }
}


void hercules_index_flush_worker(void *transferdata)
{

  int s;
  double tS = 0;
     double tE = 0;
     struct timeval start_time;
     struct timeval end_time;

  #if DETAILED_STATS == 1 
  gettimeofday(&start_time, NULL);
  #endif
  unsigned int buffer_counter = 0;
  //printf("thread = %lu unlocked barrier 2\n", pthread_self());
  //printf("thread = %d unlocking barrier 2\n", ((index_thread_data*)transferdata)->thread_id );	    

  if ((((index_thread_data*)transferdata)->buffer_counter) >= (((index_thread_data*)transferdata)->buffer_max))		
  {
      __sync_fetch_and_add(((index_thread_data*)transferdata)->flush_counter,1);

  }
  __sync_fetch_and_add(&(((index_thread_data*)transferdata)->continue_handshake),1);
  s = pthread_barrier_wait(((index_thread_data*)transferdata)->continue_barrier);
  //printf("thread = %lu   released continue_barrier \n", pthread_self());
  __sync_fetch_and_sub(&(((index_thread_data*)transferdata)->continue_handshake),1);
  
  if (__sync_fetch_and_add(((index_thread_data*)transferdata)->flush_order,0)  != 0)
  {
    buffer_counter = ((index_thread_data*)transferdata)->buffer_counter;
    s = pthread_barrier_wait(((index_thread_data*)transferdata)->flush_barrier);
    //printf("thread = %lu   released flush_barrier \n", pthread_self());    
    __sync_fetch_and_sub(& (((index_thread_data*)transferdata)->buffer_counter) , buffer_counter);
    //printf("thread = %lu  after flush buffer_counter  = %d \n", pthread_self(), ((index_thread_data*)transferdata)->buffer_counter );
  }

}


enum response hercules_index_insert_series_to_node_parallel(struct hercules_index *index, struct hercules_node *node,  ts_type * timeseries, void *thread_data)
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
 ((index_thread_data*)thread_data)->thread_traverse_total_time += (tE - tS);
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
    ((index_thread_data*)thread_data)->thread_append_total_time += (tE - tS);
    gettimeofday(&start_time, NULL); 
    #endif

    if (node->node_size >= index->settings->max_leaf_size)
    {
      hercules_index_split_node(index,node, &timeseries_segment_sketch,thread_data);
      node->is_leaf = 0;
    } 

    pthread_mutex_unlock(&(node->lock_data));

    #if DETAILED_STATS == 1
    gettimeofday(&end_time, NULL);
    tS = start_time.tv_sec*1000000 + (start_time.tv_usec);
    tE = end_time.tv_sec*1000000 + (end_time.tv_usec);
    ((index_thread_data*)thread_data)->thread_split_total_time += (tE - tS);    
    gettimeofday(&start_time, NULL); 
    #endif

   //printf("thread = %lu unlocking node = %p level = %d BEGIN\n", pthread_self(),(void *)node, node->level);    

  free(timeseries_segment_sketch.indicators);    
  return SUCCESS;
}


