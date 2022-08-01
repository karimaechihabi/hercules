//
//  hercules_file_buffer_manager.h
//
//  Created by Karima Echihabi on 18/12/2018
//

#ifndef hercules_hercules_file_buffer_manager_h
#define hercules_hercules_file_buffer_manager_h
#include "../config.h"
#include "../globals.h"
#include "ts.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hercules_index.h"
#include "hercules_file_buffer.h"

struct hercules_file_buffer_manager {
  struct hercules_file_map * file_map;
  struct hercules_file_map * file_map_tail;  
  pthread_mutex_t lock_file_map;
  pthread_mutex_t lock_file_buffer;
  
  unsigned long max_buffered_size;
  long current_count;
  long batch_remove_size;

  char *mem_array;
  char *current_record;
  int current_record_index;
  int max_record_index;
  
  int file_map_size;

};


struct hercules_file_map{
  struct hercules_file_buffer * file_buffer;
  struct hercules_file_map * next;
  struct hercules_file_map * prev;  
};

enum response init_file_buffer_manager(struct hercules_index *index);
enum response set_buffered_memory_size(struct hercules_index * index);
enum response get_file_buffer(struct hercules_index *index, struct hercules_node *node);
enum response save_all_buffers_to_disk(struct hercules_index *index);
enum response add_file_buffer_to_map(struct hercules_index * index, struct hercules_node *node);  

#endif
