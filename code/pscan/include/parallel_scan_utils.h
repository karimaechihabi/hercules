//
//  parallel_scan_utils.h
//
//  Created by Karima Echihabi on 18/01/2022

#ifndef pscanlib_parallel_scan_utils_h
#define pscanlib_parallel_scan_utils_h

#include "../config.h"
#include "../globals.h"
#include "math.h"
#include "ts.h"

typedef struct thread_data thread_data;


/// Data structure for sorting the query.
typedef struct q_index
    {   double value;
        int    index;
    } q_index;

typedef struct query_result {
    ts_type distance;
    size_t file_position;
};

struct stats_info {


        double query_total_input_time;
        double query_total_output_time;
        double query_total_load_node_time;
        double query_total_cpu_time;
        double query_total_dist1_time;
        double query_total_dist2_time;
        double query_total_time;    

        double queries_total_input_time;
        double queries_total_output_time;
        double queries_total_load_node_time;
        double queries_total_cpu_time;

        double queries_total_time;

        unsigned long long query_total_seq_input_count;
        unsigned long long query_total_seq_output_count;
        unsigned long long query_total_rand_input_count;
        unsigned long long query_total_rand_output_count;

        unsigned long long queries_total_seq_input_count;
        unsigned long long queries_total_seq_output_count;
        unsigned long long queries_total_rand_input_count;
        unsigned long long queries_total_rand_output_count;
        unsigned long long total_ts_count;

        double query_exact_distance;
        ts_type query_exact_file_position;
  
        double query_lb_distance;  
  
        double query_tlb;
        double query_sum_tlb;
  
        double queries_min_tlb;
        double queries_max_tlb;
        double queries_sum_tlb;
        double queries_sum_squares_tlb;
        double queries_avg_tlb;
        double queries_sd_tlb;
  
  
        unsigned long long tlb_ts_count; //this could be different from the total_ts_count. It includes only the ts for which d(Q,S) != 0
    
}; 

typedef struct thread_data
{
  pthread_barrier_t *DBarrier;
  ts_type * *query_ts;
  
  pthread_rwlock_t *rwl_bsf;
  pthread_mutex_t *lock_bsf;
  
  int * finished;

  thread_data *threads_data;
  int num_query_workers;
  int thread_id;

  struct query_result *knn_results;
  unsigned int * curr_k_size;
  unsigned int k;
  unsigned int ts_length;

  ts_type ** d_buffer;
  int *db_size;
  int *db_counter;
  int toggle;
  int * global_counter;
  
}thread_data;


void update_query_stats(struct stats_info *stats,unsigned int query_id,
			unsigned int found_knn, struct query_result bsf_result);
void print_query_stats(struct stats_info *stats, unsigned int query_num,unsigned int found_knn, char * queries);

int znorm_comp(const void *a, const void* b);
enum response reorder_query(ts_type * query_ts, int * query_order, int ts_length);

enum response parallel_scan(const char * dataset, int dataset_size,const char * queries, int queries_size, unsigned int ts_length, float minimum_distance, struct stats_info *stats, unsigned int k, int num_query_threads,  int dbsize);
enum response parallel_scan_process_query(FILE *dataset_file, file_position_type ts_num, unsigned int ts_length,
					  ts_type *query_ts, unsigned int k,
					  struct query_result * knn_results,
					  unsigned int num_threads, unsigned int initial_db_size);
void* queryworker(void *transferdata);
int max(int a, int b);
int min(int a, int b);
  
#endif

