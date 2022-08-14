//
//  hercules_index.h
//
//  Created by Karima Echihabi on 18/12/2018
//


#ifndef herculeslib_hercules_index_h
#define herculeslib_hercules_index_h
#include "../config.h"
#include "../globals.h"
#include "calc_utils.h"
#include "hercules_node.h"
#include "ts.h"

typedef struct index_buffer_data index_buffer_data;
typedef struct index_thread_data index_thread_data;

struct hercules_index_settings {
  const char* root_directory;  
  double buffered_memory_size;
  unsigned int timeseries_size;
  unsigned int init_segments;    
  unsigned int max_leaf_size;
  
  int ts_values_per_paa_segment;
  int paa_segments;

  sax_type * max_sax_cardinalities;
  
  unsigned int sax_byte_size;
  unsigned int position_byte_size;
  unsigned int ts_byte_size;
  unsigned int dataset_size;
  
  unsigned int full_record_size;
  unsigned int partial_record_size;
    
  // ALWAYS: SAX_ALPHABET_CARDINALITY = 2^SAX_BIT_CARDINALITY
  int sax_bit_cardinality;
  root_mask_type * bit_masks;
  int sax_alphabet_cardinality;
    
  ts_type mindist_sqrt;
  
  unsigned int max_filename_size;
  int serial;
  int sims;
  int flush_threshold;
};


struct data_record {
  ts_type * timeseries;
  sax_type  * sax;
};

struct stats_info {


	double idx_traverse_tree_total_time ;	
	double idx_traverse_tree_input_time;
	double idx_traverse_tree_output_time;
	double idx_traverse_tree_cpu_time;

	unsigned long idx_traverse_tree_seq_input_count;
	unsigned long idx_traverse_tree_seq_output_count;
	unsigned long idx_traverse_tree_rand_input_count ;
	unsigned long idx_traverse_tree_rand_output_count ;

	double idx_append_ts_to_leaf_total_time ;	
	double idx_append_ts_to_leaf_input_time;
	double idx_append_ts_to_leaf_output_time;
	double idx_append_ts_to_leaf_cpu_time;

	unsigned long idx_append_ts_to_leaf_seq_input_count;
	unsigned long idx_append_ts_to_leaf_seq_output_count;
	unsigned long idx_append_ts_to_leaf_rand_input_count ;
	unsigned long idx_append_ts_to_leaf_rand_output_count ;

	double idx_evaluate_split_policies_total_time ;	
	double idx_evaluate_split_policies_input_time;
	double idx_evaluate_split_policies_output_time;
	double idx_evaluate_split_policies_cpu_time;

	unsigned long idx_evaluate_split_policies_seq_input_count;
	unsigned long idx_evaluate_split_policies_seq_output_count;
	unsigned long idx_evaluate_split_policies_rand_input_count ;
	unsigned long idx_evaluate_split_policies_rand_output_count ;

	double idx_split_node_total_time ;	
	double idx_split_node_input_time;
	double idx_split_node_output_time;
	double idx_split_node_cpu_time;

	unsigned long idx_split_node_seq_input_count;
	unsigned long idx_split_node_seq_output_count;
	unsigned long idx_split_node_rand_input_count ;
	unsigned long idx_split_node_rand_output_count ;

	double idx_flush_worker_time;
	double idx_flush_coordinator_time;
	double idx_insert_to_node_time;
	double idx_reinsert_time;
	double idx_sanity_counter_time;
	double idx_bulk_loading_time;

        double idx_building_input_time;  
        double idx_building_output_time;
        double idx_building_cpu_time;
        double idx_building_total_time;

        unsigned long long idx_building_seq_input_count;
        unsigned long long idx_building_seq_output_count;
        unsigned long long idx_building_rand_input_count;
        unsigned long long idx_building_rand_output_count;

        double idx_writing_input_time;  
        double idx_writing_output_time;
        double idx_writing_cpu_time;
        double idx_writing_total_time;

        double idx_reading_input_time;  
        double idx_reading_output_time;
        double idx_reading_cpu_time;
        double idx_reading_total_time;
  
        double idx_total_input_time;
        double idx_total_output_time;
        double idx_total_cpu_time;
        double idx_total_time;
  

        unsigned long long idx_writing_seq_input_count;
        unsigned long long idx_writing_seq_output_count;
        unsigned long long idx_writing_rand_input_count;
        unsigned long long idx_writing_rand_output_count;

        unsigned long long idx_reading_seq_input_count;
        unsigned long long idx_reading_seq_output_count;
        unsigned long long idx_reading_rand_input_count;
        unsigned long long idx_reading_rand_output_count;
  
        unsigned long long idx_total_seq_input_count;
        unsigned long long idx_total_seq_output_count;
        unsigned long long idx_total_rand_input_count;
        unsigned long long idx_total_rand_output_count;

        double query_filter_input_time;
        double query_filter_output_time;
        double query_filter_load_node_time;
        double query_filter_cpu_time;
        double query_filter_total_time;    

        double queries_filter_input_time;
        double queries_filter_output_time;
        double queries_filter_load_node_time;
        double queries_filter_cpu_time;
        double queries_filter_total_time;

        double query_refine_input_time;
        double query_refine_output_time;
        double query_refine_load_node_time;
        double query_refine_cpu_time;
        double query_refine_total_time;    

        double queries_refine_input_time;
        double queries_refine_output_time;
        double queries_refine_load_node_time;
        double queries_refine_cpu_time;
        double queries_refine_total_time;

        double query_total_input_time;
        double query_total_output_time;
        double query_total_load_node_time;
        double query_total_cpu_time;
        double query_total_time;    

        double queries_total_input_time;
        double queries_total_output_time;
        double queries_total_load_node_time;
        double queries_total_cpu_time;
        double queries_total_time;

        double total_input_time;
        double total_output_time;
        double total_load_node_time;
        double total_cpu_time;
        double total_time;    
        double total_time_sanity_check;
  
        unsigned long long query_refine_seq_input_count; //refine for refinement step
        unsigned long long query_refine_seq_output_count;
        unsigned long long query_refine_rand_input_count;
        unsigned long long query_refine_rand_output_count;

        unsigned long long queries_refine_seq_input_count;
        unsigned long long queries_refine_seq_output_count;
        unsigned long long queries_refine_rand_input_count;
        unsigned long long queries_refine_rand_output_count;

        unsigned long long query_filter_seq_input_count; //filtering step
        unsigned long long query_filter_seq_output_count;
        unsigned long long query_filter_rand_input_count;
        unsigned long long query_filter_rand_output_count;

        unsigned long long queries_filter_seq_input_count;
        unsigned long long queries_filter_seq_output_count;
        unsigned long long queries_filter_rand_input_count;
        unsigned long long queries_filter_rand_output_count;

        unsigned long long query_total_seq_input_count;
        unsigned long long query_total_seq_output_count;
        unsigned long long query_total_rand_input_count;
        unsigned long long query_total_rand_output_count;

        unsigned long long queries_total_seq_input_count;
        unsigned long long queries_total_seq_output_count;
        unsigned long long queries_total_rand_input_count;
        unsigned long long queries_total_rand_output_count;

        unsigned long long total_seq_input_count;
        unsigned long long total_seq_output_count;
        unsigned long long total_rand_input_count;
        unsigned long long total_rand_output_count;
  
        unsigned int total_nodes_count;
        unsigned int leaf_nodes_count;
        unsigned int empty_leaf_nodes_count;  
  

        long long idx_size_bytes;
        long long idx_size_blocks; 


        int *  leaves_sizes;  
        double min_fill_factor;
        double max_fill_factor;
        double sum_fill_factor;
        double sum_squares_fill_factor;
        double avg_fill_factor;
        double sd_fill_factor;

        int *  leaves_heights;  
        unsigned int min_height;
        unsigned int max_height;
        unsigned int sum_height;
        unsigned int sum_squares_height;  
        double avg_height;
        double sd_height;
  
        int leaves_counter;
  
        //file_position_type loaded_records_count;
  /*
        unsigned long long query_refine_ts_disk_count;
        unsigned long long query_refine_ts_mem_count;
        unsigned long long query_filter_ts_disk_count;
        unsigned long long query_filter_ts_mem_count;
  */
  
        unsigned int query_loaded_nodes_count;
        unsigned int query_checked_nodes_count;
	unsigned long long query_loaded_ts_count;
        unsigned long long query_checked_ts_count;
  
        unsigned long long total_ts_count;


        unsigned int query_filter_loaded_nodes_count;
        unsigned int query_filter_checked_nodes_count;
	unsigned long long query_filter_loaded_ts_count;
        unsigned long long query_filter_checked_ts_count;

        unsigned int query_refine_loaded_nodes_count;
        unsigned int query_refine_checked_nodes_count;
	unsigned long long query_refine_loaded_ts_count;
        unsigned long long query_refine_checked_ts_count;

        unsigned int query_total_loaded_nodes_count;
        unsigned int query_total_checked_nodes_count;
	unsigned long long query_total_loaded_ts_count;
        unsigned long long query_total_checked_ts_count;
  




        double query_approx_distance;
        char * query_approx_node_filename;  
        unsigned int query_approx_node_size;
        unsigned int query_approx_node_level;  

        double query_exact_distance;
        char * query_exact_node_filename;
        unsigned int query_exact_node_size;
        unsigned int query_exact_node_level;
        unsigned long long query_exact_node_file_pos;
  
        double query_lb_distance;
  
        double query_eff_epsilon;
        double query_pruning_ratio;
        double query_tlb;  
        double query_sum_tlb;
  

        double queries_min_pruning_ratio;
        double queries_max_pruning_ratio;
        double queries_sum_pruning_ratio;
        double queries_sum_squares_pruning_ratio;
        double queries_avg_pruning_ratio;
        double queries_sd_pruning_ratio;
  
        double queries_min_eff_epsilon;
        double queries_max_eff_epsilon;
        double queries_sum_eff_epsilon;
        double queries_sum_squares_eff_epsilon;
        double queries_avg_eff_epsilon;
        double queries_sd_eff_epsilon;

        double queries_min_tlb;
        double queries_max_tlb;
        double queries_sum_tlb;
        double queries_sum_squares_tlb;
        double queries_avg_tlb;
        double queries_sd_tlb;
  
        double queries_total_checked_ts_count;
        double queries_total_checked_nodes_count;  
        double queries_avg_checked_ts_count;
        double queries_avg_checked_nodes_count;
        double queries_total_loaded_ts_count;
        double queries_total_loaded_nodes_count;  
        double queries_avg_loaded_ts_count;
        double queries_avg_loaded_nodes_count;

        unsigned long long tlb_ts_count; //this could be different from the total_ts_count. It includes only the ts for whiche d(Q,S) != 0
        unsigned long long eff_epsilon_queries_count; //the number of queries for which  d(Q,S) != 0

  
        unsigned int total_queries_count; //this could be different from the input queries to eliminate cases d(Q,S) = 0

        double total_parse_time;
 };


struct  hercules_index{
  unsigned long long total_records;
  struct hercules_node *first_node;
  struct hercules_index_settings *settings;
  struct hercules_file_buffer_manager * buffer_manager;
  struct stats_info * stats;
  
  const char * leaves_raw_filename;
  FILE * leaves_raw_file;
  unsigned long long leaves_raw_pos;
  struct hercules_node ** leaves;

  const char * leaves_sims_filename;
  FILE * leaves_sims_file;
  unsigned long long leaves_sims_pos;
  sax_type *sax_cache;
  unsigned long sax_cache_size;
  
};

typedef struct index_thread_data
{
  struct hercules_index *index;
  pthread_barrier_t *db_barrier;
  pthread_barrier_t *continue_barrier;
  pthread_barrier_t *flush_barrier;



  double thread_traverse_total_time;
  double thread_evaluate_total_time;  
  double thread_append_input_time;
  double thread_append_total_time;  
  double thread_split_input_time;
  double thread_split_total_time;
  double thread_flush_output_time;
  double thread_flush_worker_time;
  double thread_flush_coordinator_time;
  double thread_insert_to_node_time;
  double thread_reinsert_time;
  double thread_sanity_counter_time;
  double thread_bulk_loading_time;

  unsigned int buffer_max;
  unsigned int buffer_counter;
  char * buffer_offset;

  int continue_handshake;
  //bool * flush_order;
  bool is_flusher;
  index_thread_data *threads_data;
  int num_insert_workers;
  int thread_id;
  
  //ts_type * double_buffer;
  ts_type ** double_buffer;
  ts_type * ts;
  unsigned int * db_size; 
  unsigned int * db_counter;
  bool *finished;
  int * flush_counter;
  int * flush_order;  
  unsigned int initial_db_size;

  ts_type **split_node_data;
  sax_type **node_sax_data;

  int current_leaf;


  int *sanity_counter;
  int fin_number;


}index_thread_data;


typedef struct index_buffer_data
{
  //int ts_num,ts_loaded;
  file_position_type pos;
  struct hercules_index *index;
  ts_type * ts;
  sax_type *saxv;
  int fin_number,blocid;
  pthread_mutex_t *lock_record;
  pthread_mutex_t *lock_fbl;
  pthread_mutex_t *lock_index;
  pthread_mutex_t *lock_cbl;
  pthread_mutex_t *lock_firstnode;
  pthread_mutex_t *lock_nodeconter;
  pthread_mutex_t *lock_disk;
  file_position_type *fbl;
  int *bufferpresize;
  int workernumber;
  int total_workernumber;
  int *nodecounter;
  int *sanity_counter;
 
  pthread_barrier_t *lock_barrier1, *lock_barrier2, *lock_barrier3, *lock_barrier4;
  bool finished;
  double thread_traverse_total_time;
  double thread_evaluate_total_time;  
  double thread_append_input_time;
  double thread_append_total_time;  
  double thread_split_input_time;
  double thread_split_total_time;
  double thread_flush_output_time;
  double thread_flush_worker_time;
  double thread_flush_coordinator_time;
  double thread_insert_to_node_time;
  double thread_reinsert_time;
  double thread_sanity_counter_time;
  double thread_bulk_loading_time;

  int buffer_max;
  char * buffer_offset;
  int buffer_counter;
  int flush_handshake;
  //bool * flush_order;
  int * flush_order;  
  bool is_flusher;
  int *flush_counter;
  index_buffer_data *threads_data;
  int num_insert_workers;
  //int * busy_insert_workers;  
  bool is_busy;
  int thread_id;
  bool is_full;
  bool barrier4_reached;  
  
  ts_type **split_node_data;
  sax_type **node_sax_data;

  int current_leaf;
}index_buffer_data;


struct hercules_index * hercules_index_init(struct hercules_index_settings *settings);
struct hercules_index * hercules_index_read(const char* root_directory); 
enum response hercules_index_write_sequential(struct hercules_index *index);
enum response hercules_index_write_parallel(struct hercules_index *index, int num_write_threads);

struct hercules_node * hercules_node_read(struct hercules_index *index, FILE *file);
enum response hercules_node_write(struct hercules_index *index, struct hercules_node *node, FILE *file);
enum response hercules_node_write_parallel(struct hercules_index *index, struct hercules_node *node, FILE *file);


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
							  int flush_threshold);

void hercules_index_destroy(struct hercules_index *index, struct hercules_node *node, boolean is_index_new);
void destroy_buffer_manager(struct hercules_index *index);
enum response hercules_index_insert(struct hercules_index *index,  ts_type * timeseries);
enum response append_ts_to_buffered_list(struct hercules_index * index, struct hercules_node * node, ts_type * timeseries);

void hercules_print_stats (struct hercules_index *index);
enum response hercules_update_index_stats(struct hercules_index *index, struct hercules_node *node);
enum response hercules_init_stats(struct hercules_index * index);
void print_tlb_stats(struct hercules_index * index, unsigned int query_num, char * queries);
void cache_sax_file(struct hercules_index *index);
enum response hercules_index_insert_parallel(struct hercules_index *index, struct hercules_node *node,  ts_type * timeseries, void * data);  
//struct hercules_ts_buffer * copy_ts_in_mem(struct hercules_index * index, struct hercules_node * node);
void flush_index_to_distk(struct hercules_index *index, struct hercules_node * node);
void flush_worker(void *transferdata);
void flush_coordinator(void *transferdata);
enum response flush_leaves_parallel(struct hercules_index *index, int num_threads);
enum response get_leaves(struct hercules_index *index, struct hercules_node *node);


enum response hercules_index_write(struct hercules_index *index, int num_write_threads);
enum response hercules_index_build(const char *ifilename, file_position_type data_size,struct hercules_index *index, int num_threads, unsigned int initial_db_size);
void* hercules_index_insert_worker(void *transferdata);
void hercules_index_flush_coordinator(void *transferdata);
void hercules_index_flush_worker(void *transferdata);
enum response hercules_index_insert_series_to_node_parallel(struct hercules_index *index, struct hercules_node *node,  ts_type * timeseries, void *thread_data);
enum response hercules_index_set_stats(struct hercules_index *index, struct hercules_node *node);
enum response hercules_index_node_write(struct hercules_index *index, struct hercules_node *node, FILE *file);



#endif
