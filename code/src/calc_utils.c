//
//  calc_utils.c
//
//  Created by Karima Echihabi on 18/12/2018
//



#include "../config.h"
#include "../globals.h"
#include "../include/calc_utils.h"
#include "math.h"
#include "immintrin.h"
#include "assert.h"

void calc_mean_stdev (ts_type * series, int start, int end, ts_type * mean, ts_type * stdev)
{
  ts_type sum_x_squares=0, sum_x=0; //sum of x's and sum of x squares
  int i, count_x;

  *stdev = 0;
  *mean = 0;
  
  if (start >= end)
  {
    printf ("error in stdev start >= end\n");
  }
  else
  {
    count_x = end-start; //size of the series
  
    for (int i=start; i<end; i++) 
    {
     sum_x += series[i];
     sum_x_squares += series[i] * series[i];
    }

    *mean = sum_x/count_x; 
    sum_x_squares -= ((sum_x * sum_x) / count_x);

    //DO WE REALLY NEED SQRT???
    *stdev = sqrt(sum_x_squares/count_x);
  }
  
}
/*
void calc_mean_stdev_SIMD (ts_type * series, int start, int end, ts_type * mean, ts_type * stdev)
{
  float sum_x_squares=0, sum_x=0; //sum of x's and sum of x squares
  int i, count_x;
  __m256 v_t,v_s,v_d,v_sum, v_sum_squares;
  __m128 v_t_128,v_s_128,v_d_128,v_sum_128, v_sum_squares_128;
  
  float sumf[8], sum_squaresf[8];

  *stdev = 0;
  *mean = 0;
  
  if (start >= end)
  {
    printf ("error in stdev start >= end\n");
  }
  else
  {
    count_x = end-start; //size of the series

    int d = end-start;
    ts_type * x = series + start;

    while (d >= 8)
      {
	v_t=_mm256_loadu_ps (x); x = x + 8;
      
	v_sum = _mm256_hadd_ps (v_t, v_t);
	v_sum = _mm256_hadd_ps (v_sum, v_sum);
	_mm256_storeu_ps (sumf ,v_sum);
	sum_x +=sumf[0]+sumf[4];
      
	v_s=_mm256_mul_ps (v_t,v_t);
	v_sum_squares = _mm256_hadd_ps (v_s, v_s);
	v_sum_squares = _mm256_hadd_ps (v_sum_squares, v_sum_squares);
	_mm256_storeu_ps (sum_squaresf ,v_sum_squares);
	sum_x_squares +=sum_squaresf[0]+sum_squaresf[4];
  
	d -= 8;  
      }
    
    if (d >= 4)
      {

	v_t_128 = _mm128_loadu_ps (x);
	x = x + 4;
      
	v_sum_128 = _mm128_hadd_ps (v_t_128, v_t_128);
	v_sum_128 = _mm128_hadd_ps (v_sum_128, v_sum_128);
	_mm128_storeu_ps (sumf ,v_sum_128);
	sum_x +=sumf[0];
      
	v_s_128=_mm128_mul_ps (v_t_128,v_t_128);
	v_sum_squares_128 = _mm128_hadd_ps (v_s_128, v_s_128);
	v_sum_squares_128 = _mm128_hadd_ps (v_sum_squares_128, v_sum_squares_128);
	_mm128_storeu_ps (sum_squaresf ,v_sum_squares_128);
	sum_x_squares +=sum_squaresf[0];
  
	d -= 4;  
      }
    
    if (d > 0)
      {
	v_t_128 = masked_read(d,x);
	v_sum_128 =  _mm128_hadd_ps (v_t_128, v_t_128);
	v_sum_128 = _mm128_hadd_ps (v_sum_128, v_sum_128);
	_mm128_storeu_ps (sumf ,v_sum_128);
	sum_x +=sumf[0];

	v_s_128 = _mm128_mul_ps (v_t_128,v_t_128);
	v_sum_squares_128 = _mm128_hadd_ps (v_s_128, v_s_128);
	v_sum_squares_128 = _mm128_hadd_ps (v_sum_squares_128, v_sum_squares_128);
	_mm128_storeu_ps (sum_squaresf ,v_sum_squares_128);
	sum_x_squares +=sum_squaresf[0];
      }
    
    *mean = sum_x/count_x; 
    sum_x_squares -= ((sum_x * sum_x) / count_x);

    *stdev = sqrt(sum_x_squares/count_x);
  }
  
}
*/

#ifdef __SSE__
void calc_mean_stdev_SIMD (ts_type * series, int start, int end, ts_type * mean, ts_type * stdev)
{
  float sum_x_squares=0, sum_x=0; //sum of x's and sum of x squares
  int i, count_x;
  __m256 v_t,v_s,v_d,v_sum, v_sum_squares;
  __m128 v_t_128,v_s_128,v_d_128,v_sum_128, v_sum_squares_128;
  
  float sumf[8], sum_squaresf[8];

  *stdev = 0;
  *mean = 0;
  
  if (start >= end)
  {
    printf ("error in stdev start >= end\n");
  }
  else
  {
    count_x = end-start; //size of the series

    int d = end-start;
    ts_type * x = series + start;

    while (d >= 8)
      {
	v_t=_mm256_loadu_ps (x); x = x + 8;
      
	v_sum = _mm256_hadd_ps (v_t, v_t);
	v_sum = _mm256_hadd_ps (v_sum, v_sum);
	_mm256_storeu_ps (sumf ,v_sum);
	sum_x +=sumf[0]+sumf[4];
      
	v_s=_mm256_mul_ps (v_t,v_t);
	v_sum_squares = _mm256_hadd_ps (v_s, v_s);
	v_sum_squares = _mm256_hadd_ps (v_sum_squares, v_sum_squares);
	_mm256_storeu_ps (sum_squaresf ,v_sum_squares);
	sum_x_squares +=sum_squaresf[0]+sum_squaresf[4];
  
	d -= 8;  
      }
    
    if (d > 0)
      {
	v_t = masked_read_8 (d,x);
      
	v_sum = _mm256_hadd_ps (v_t, v_t);
	v_sum = _mm256_hadd_ps (v_sum, v_sum);
	_mm256_storeu_ps (sumf ,v_sum);
	sum_x +=sumf[0]+sumf[4];
      
	v_s=_mm256_mul_ps (v_t,v_t);
	v_sum_squares = _mm256_hadd_ps (v_s, v_s);
	v_sum_squares = _mm256_hadd_ps (v_sum_squares, v_sum_squares);
	_mm256_storeu_ps (sum_squaresf ,v_sum_squares);
	sum_x_squares +=sum_squaresf[0]+sum_squaresf[4];
	
      }
    
    *mean = sum_x/count_x; 
    sum_x_squares -= ((sum_x * sum_x) / count_x);

    *stdev = sqrt(sum_x_squares/count_x);
  }
  
}
 
__m128 masked_read (int d, const float *x)
{
    assert (0 <= d && d < 4);
    __attribute__((__aligned__(16))) float buf[4] = {0, 0, 0, 0};
    switch (d) {
      case 3:
        buf[2] = x[2];
      case 2:
        buf[1] = x[1];
      case 1:
        buf[0] = x[0];
    }
    return _mm_load_ps (buf);
    // cannot use AVX2 _mm_mask_set1_epi32
}

__m256 masked_read_8 (int d, const float *x)
{
    assert (0 <= d && d < 8);
    __attribute__((__aligned__(32))) float buf[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    switch (d) {      
      case 7:
        buf[6] = x[6];      
      case 6:
        buf[5] = x[5];
      case 5:
        buf[4] = x[4];
      case 4:
        buf[3] = x[3];
      case 3:
        buf[2] = x[2];
      case 2:
        buf[1] = x[1];
      case 1:
        buf[0] = x[0];
    }
    return _mm256_load_ps (buf);
    // cannot use AVX2 _mm_mask_set1_epi32
}

void calc_mean_stdev_per_segment_SIMD (ts_type * series, short * segments, ts_type *means, ts_type *stdevs, int size)
{
  int start=0, end;
  for (int i=0; i< size; i++) 
  {
    end = segments[i];
    calc_mean_stdev_SIMD (series, start, end,&means[i],&stdevs[i]);
    start = end;
  }
}

#endif

ts_type calc_mean (ts_type * series, int start, int end)
{

  ts_type mean = 0;

  if (start >= end)
  {
    int j = 0;
    j++;
    printf("error start > end \n");
  }
  else
  {
   for (int i=start; i < end; i++) 
   {
     mean += series[i];
   }
  
   mean /= (end-start); 
  }

  return mean;
  
}

void calc_mean_stdev_per_segment (ts_type * series, short * segments, ts_type *means, ts_type *stdevs, int size)
{
  int start=0, end;
  for (int i=0; i< size; i++) 
  {
    end = segments[i];
    calc_mean_stdev (series, start, end,&means[i],&stdevs[i]);
    start = end;
  }
}





/*
  Using the stdev computational formula.

*/

ts_type calc_stdev (ts_type * series, int start, int end)
{

  ts_type sum_x_squares=0, sum_x=0, stdev = 0; //sum of x's and sum of x squares
  int i, count_x;

  if (start >= end)
  {
    printf ("error in stdev start >= end\n");
  }
  else
  {
    count_x = end-start; //size of the series
  
    for (int i=start; i<end; i++) 
    {
     sum_x += series[i];
     sum_x_squares += pow(series[i],2);
    }
  
    sum_x_squares -= (pow(sum_x,2) / count_x);

    stdev = sqrt(sum_x_squares/count_x);
  }
  
  return stdev;
  
}

ts_type calc_mean_per_segment (ts_type * series, short * segments, ts_type *means, int size)
{

  int start=0, end;
  for (int i=0; i< size; i++) 
  {
    end = segments[i];
    means[i] = calc_mean (series, start, end);
    start = end;
  }
  
}

ts_type calc_stdev_per_segment (ts_type * series, short * segments, ts_type *stdevs, int size)
{

  int start=0, end;
  for (int i=0; i< size; i++) 
  {
    end = segments[i];
    stdevs[i] = calc_stdev (series, start, end);
    start = end;
  }
  
}
/*
short compare_leaf_pos (const void * a, const void * b)
{
 const struct hercules_node *const *leaf_a = a;
  const struct hercules_node *const *leaf_b = b;

  if ((*leaf_a)->file_pos < (*leaf_b)->file_pos )
	{
      printf ("a = %llu, b = %llu \n",leaf_a->file_pos,leaf_b->file_pos);
    return -1;
	}
  else if  ((*leaf_a)->file_pos == (*leaf_b)->file_pos)
    return 0;
  else
    return 1;
}

*/
short compare_leaf_pos (const void * a, const void * b)
{
  struct query_result leaf_a = *((struct query_result*) a);
  struct query_result leaf_b = *((struct query_result*) b);

  if (leaf_a.node->file_pos < leaf_b.node->file_pos )
    return -1;
  else if  (leaf_a.node->file_pos == leaf_b.node->file_pos)
    return 0;
  else
    return 1;
}


/* 
 This is the compare function used in the binary search code
 */

  

ts_type compare_ts_type (const void * a, const void * b)
{
  return ( *(ts_type*)a - *(ts_type*)b );
}

  

short compare_short (const void * a, const void * b)
{
  if (*(short*)a < *(short*)b )
    return -1;
  else if (*(short*)a == *(short*)b )
    return 0;
  else
    return 1;

}

unsigned long long ull_min (unsigned long long a, unsigned long long b)
{
  if (a < b )
    return a;
  else
    return b;
}

/*
short compare_file_map_entries (const void * a, const void * b)
{
  char * entry_a = (char *) a;
  struct hercules_file_map *entry_b = (struct hercules_file_map*) b;

  return ( strcmp(entry_a, entry_b->filename));

}
*/
short compare_file_buffer (const void * a, const void * b)
{
  struct hercules_file_buffer * entry_a = *((struct hercules_file_buffer**) a);
  struct hercules_file_buffer * entry_b = *((struct hercules_file_buffer**) b);

  if (entry_a->buffered_list_size < entry_b->buffered_list_size )
    return 1;
  else if  (entry_a->buffered_list_size == entry_b->buffered_list_size )
    return 0;
  else
    return -1;
}

/*
  returns the current time in string format.
*/

void get_current_time(char * time_buf)
{
    time_t timer;
    
    struct tm* tm_info;

    time(&timer);
    tm_info = localtime(&timer);

    strftime(time_buf, 26, "%Y-%m-%d %H:%M:%S", tm_info);
}

ts_type calculate_mean_std_dev_range(struct segment_sketch sketch, int len){

  ts_type mean_width = sketch.indicators[0]-sketch.indicators[1];

  ts_type stdev_upper = sketch.indicators[2];
  ts_type stdev_lower = sketch.indicators[3];

  return (len * (mean_width * mean_width + stdev_upper * stdev_upper));
  
}


int get_segment_start(short * points, int idx)
{
  if (idx == 0 )
    return 0;
  else
    return points[idx-1];
}

int get_segment_end(short * points, int idx)
{
  return points[idx];
}

int get_segment_length(short * points, int i)
{

  if (i == 0)
    return points[i];
  else
    return points[i] - points[i-1];

}


