//
//  cals_utils.h
//
//  Created by Karima Echihabi on 18/12/2018
//

#ifndef herculeslib_hercules_calc_utils_h
#define herculeslib_hercules_calc_utils_h

#include "../config.h"
#include "../globals.h"
#include "math.h"
#include "hercules_file_buffer_manager.h"
#include "immintrin.h"

void calc_mean_stdev (ts_type * series, int start, int end, ts_type * mean, ts_type * stdev);
void calc_mean_stdev_per_segment (ts_type * series, short * segments, ts_type *means, ts_type *stdevs, int size);

ts_type calc_mean (ts_type *, int start, int end);
ts_type calc_stdev (ts_type *, int start, int end);
int compare(const void *a, const void *b);
short compare_short (const void * a, const void * b);
//short compare_file_map_entries (const void * a, const void * b);
void get_current_time(char * time_buf);
short compare_file_buffer (const void * a, const void * b);
ts_type calculate_mean_std_dev_range(struct segment_sketch sketch, int len);
int get_segment_start(short * points, int idx);
int get_segment_end(short * points, int idx);
int get_segment_length(short * points, int i);
int znorm_comp(const void *a, const void* b);
short compare_leaf_pos (const void * a, const void * b);
unsigned long long ull_min (unsigned long long a, unsigned long long b);

#ifdef __SSE__
void calc_mean_stdev_SIMD (ts_type * series, int start, int end, ts_type * mean, ts_type * stdev);
__m128 masked_read (int d, const float *x);
__m256 masked_read_8 (int d, const float *x);
void calc_mean_stdev_per_segment_SIMD (ts_type * series, short * segments, ts_type *means, ts_type *stdevs, int size);
#endif
#endif
