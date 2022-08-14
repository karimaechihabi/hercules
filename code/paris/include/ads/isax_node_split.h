//
//  isax_node_split.h
//  isaxlib
//
//  Created by Kostas Zoumpatianos on 3/12/12.
//  Copyright 2012 University of Trento. All rights reserved.
//

#ifndef isaxlib_isax_node_split_h
#define isaxlib_isax_node_split_h
#include "../../config.h"
#include "../../globals.h"
#include "isax_index.h"
#include "isax_node.h"

int simple_split_decision (isax_node_split_data * split_data, 
                                  isax_index_settings * settings);

int informed_split_decision (isax_node_split_data * split_data, 
                             isax_index_settings * settings,
                             isax_node_record * records_buffer,
                             int records_buffer_size);

void split_node(isax_index *index, isax_node *node);
void split_node_inmemory(isax_index *index, isax_node *node);

#endif
