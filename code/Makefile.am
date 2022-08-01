AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS =  ${ACLOCAL_FLAGS}
AM_CFLAGS="-g -O0"

lib_LIBRARIES= lib/libhercules.a 

lib_libhercules_a_SOURCES = src/pqueue.c src/hercules_file_loaders.c  src/hercules_index.c src/hercules_node.c src/hercules_node_split.c  src/ts.c  src/sax.c src/calc_utils.c src/hercules_query_engine.c  src/hercules_file_buffer.c src/hercules_file_buffer_manager.c  src/sax.c 

lib_libhercules_a_CFLAGS=-I/opt/local/include -Iinclude  -march=native -mavx -mavx2 -msse3 -fopenmp  

bin_PROGRAMS = bin/hercules

bin_hercules_SOURCES=src/hercules.c
bin_hercules_CFLAGS = -I/opt/local/include -Iinclude -I/usr/local/include/   -march=native
bin_hercules_LDADD = -lreadline -L$/lib -lhercules -lm  -lblas -lpthread
bin_hercules_LDFLAGS = -L/opt/local/lib -Llib/ -Wl,-rpath,-r/lib  -mavx -mavx2 -msse3 -fopenmp 
