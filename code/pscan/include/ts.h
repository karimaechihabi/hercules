//
//  ts.h
//  isaxlib
//
//  Created by Kostas Zoumpatianos on 3/7/12.
//  Copyright 2012 University of Trento. All rights reserved.
//

#ifndef dstreelib_ts_h
#define dstreelib_ts_h
#include "../config.h"
#include "../globals.h"
#include <math.h>

enum response ts_parse_str(char ts_str[], ts_type *ts_out, int ts_size, const char * delims);
void ts_print(ts_type *ts, int size);
ts_type ts_euclidean_distance_reordered(ts_type * q, ts_type * t , int j , int  m ,ts_type bsf, int* order);

float ts_euclidean_distance(ts_type * t, ts_type * s, int size, ts_type bsf);
ts_type ts_euclidean_distance_SIMD(ts_type * t, ts_type * s, int size, ts_type bound) ;

#endif