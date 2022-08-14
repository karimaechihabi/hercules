//
//  ts.h
//  isaxlib
//
//  Created by Kostas Zoumpatianos on 3/7/12.
//  Copyright 2012 University of Trento. All rights reserved.
//

#ifndef isaxlib_ts_h
#define isaxlib_ts_h
#include "../../../config.h"
#include "../../../globals.h"
void ts_parse_str(char ts_str[], ts_type *ts_out, int ts_size, const char * delims);
void ts_print(ts_type *ts, int size);
float ts_euclidean_distance(ts_type * t, ts_type * s, int size, float bound);
float ts_euclidean_distance_SIMD(ts_type * t, ts_type * s, int size, float bound);
float ts_euclidean_distance_neSIMD(ts_type * t, ts_type * s, int size, float bound);
#endif
