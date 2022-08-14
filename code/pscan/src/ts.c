//
//  ts.c
//  isaxlib
//
//  Created by Kostas Zoumpatianos on 3/7/12.
//  Copyright 2012 University of Trento. All rights reserved.
//
#include "../config.h"
#include "../globals.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "../include/ts.h"
#include "immintrin.h"

/**
 This function converts a string of floats seperated by a delimeter into a ts 
 record of a size ts_size.
 @param char ts_str[]
 @param int ts_size
 @param const char * delims
 @return *ts
 */
enum response ts_parse_str(char ts_str[], ts_type * ts_out, int ts_size, const char * delims)
{
    int index=0;
    char *result = strtok( ts_str, delims );
	while( result != NULL ) {
   	        ts_out[index] =  atof(result);
		result = strtok( NULL, delims );
#ifdef SANITY_CHECK
        if (index >= ts_size)
        {
            fprintf(stderr, "Error in ts.c: Time series bigger than limit of %d.\n", ts_size);
            return FAILURE; 
        }
#endif
        index++;
	}
    free(result);
    return SUCCESS;
}

ts_type ts_euclidean_distance(ts_type * q, ts_type * t, int size, ts_type bsf) {

   ts_type distance = 0;
    while (size > 0 && distance < bsf) {
        size--;
        distance += (t[size] - q[size]) * (t[size] - q[size]);
    }
//    distance = sqrtf(distance);
    
    return distance;
}


/// Main function for calculating ED distance between the query, Q, and current data, T.
/// Note that Q is already sorted by absolute z-normalization value, |z_norm(Q[i])|
//ADD later std and mean for early abandoning z-normalization
//ts_type ts_euclidian_distance_reordered(ts_type * q, ts_type * t , int j , int  m , ts_type  mean , ts_type std , int* order, ts_type bsf)

/*ts_type ts_euclidean_distance_reordered(ts_type * q, ts_type * t , int j , int  size ,ts_type bsf, int * order)
{
    ts_type x;
    ts_type distance = 0;
    while (size > 0 && distance < bsf) {
        size--;
        x = t[order[size]];	
        distance += (x - q[size]) * (x - q[size]);
    }

    return distance;
}
*/

ts_type ts_euclidean_distance_reordered(ts_type * q, ts_type * t , int j , int  size ,ts_type bsf, int * order)
{
    int i;
    ts_type sum = 0;
    for ( i = 0 ; i < size && sum < bsf ; i++ )
    {
       //ts_type x = (T[(order[i]+j)]-mean)/std;
      ts_type x = t[order[i]];
      sum += (x-q[i])*(x-q[i]);      
    }
    return sum;
}




/** 
 This function prints a ts record of a size.
 @param ts *ts
 @param int size
*/
void ts_print(ts_type *ts, int size) 
{
    int i;
    for (i=0; i < size; i++) {
        printf("%lf", ts[i]);
    }
    printf("\n");
}


ts_type ts_euclidean_distance_SIMD(ts_type * t, ts_type * s, int size, ts_type bound) {
float distance = 0;
    int i =0;
float distancef[8];

__m256 v_t,v_s,v_d,distancev;
    while (size > 0 && distance <bound) {
        v_t=_mm256_loadu_ps (&t[i]);
        v_s=_mm256_loadu_ps (&s[i]);
        
        v_d= _mm256_sub_ps (v_t,v_s);

        v_d=_mm256_mul_ps (v_d,v_d);
        size-=8;

        i=i+8;
        distancev = _mm256_hadd_ps (v_d, v_d);
        distancev = _mm256_hadd_ps (distancev, distancev);
        _mm256_storeu_ps (distancef ,distancev);
        distance +=distancef[0]+distancef[4];

    }

//    distance = sqrtf(distance);
    
    return (ts_type) distance;
}
