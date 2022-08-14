AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS =  ${ACLOCAL_FLAGS}
AM_CFLAGS=""

lib_LIBRARIES= lib/libpscan.a 

lib_libpscan_a_SOURCES =  src/ts.c src/parallel_scan_utils.c 

lib_libpscan_a_CFLAGS=-I/opt/local/include -Iinclude  -march=native -mavx -mavx2 -msse3 -fopenmp

bin_PROGRAMS = bin/pscan

bin_pscan_SOURCES=src/parallel_scan.c
bin_pscan_CFLAGS = -I/opt/local/include -Iinclude -I/usr/local/include/  -march=native
bin_pscan_LDADD = -lreadline -L$/lib  -lpscan -lm -lblas -lpthread
bin_pscan_LDFLAGS = -L/opt/local/lib -Llib/ -Wl,-rpath,-r/lib -mavx -mavx2 -msse3 -fopenmp

