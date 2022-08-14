//
//  main.c
//  isaxlib
//
//  Created by Kostas Zoumpatianos on 3/12/12.
//  Copyright 2012 University of Trento. All rights reserved.
//
#define _GNU_SOURCE

#define PRODUCT "----------------------------------------------\
\nThis is the Adaptive Leaf iSAX index.\n\
Copyright (C) 2011-2014 University of Trento.\n\
----------------------------------------------\n\n"
#ifdef VALUES
#include <values.h>
#endif

#include "../../config.h"
#include "../../globals.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <unistd.h>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include <float.h>
#include <sched.h>

#include "ads/sax/sax.h"
#include "ads/sax/ts.h"
#include "ads/isax_visualize_index.h"
#include "ads/isax_file_loaders.h"
#include "ads/isax_visualize_index.h"
#include "ads/isax_first_buffer_layer.h"
#include "ads/isax_query_engine.h"
#include "ads/inmemory_query_engine.h"
#include "ads/parallel_inmemory_query_engine.h"
#include "ads/inmemory_index_engine.h"
#include "ads/parallel_query_engine.h"
#include "ads/parallel_index_engine.h"
#include "ads/inmemory_topk_engine.h"
#include "ads/DTWfunction.h"
//#define PROGRESS_CALCULATE_THREAD_NUMBER 4
//#define PROGRESS_FLUSH_THREAD_NUMBER 4
//#define QUERIES_THREAD_NUMBER 4
//#define DISK_BUFFER_SIZE 32

isax_index *idx;
void     INThandler(int);


int main (int argc, char **argv)
{
	signal(SIGINT, INThandler);

#ifndef BENCHMARK
    printf(PRODUCT);

#if VERBOSE_LEVEL == 0
    printf("Executing in silent mode. Please wait.\n");
#endif
#endif
    static char * dataset = "/home/botao/document/";
    static char * queries = "/home/botao/document/";
    static char * index_path = "/home/botao/document/myexperiment/";
    static char * labelset="/home/botao/document/myexperiment/";
    //static int dataset_size = 20100;//1000000;
    //static int dataset_size =   1000000;//simple2
    static long int dataset_size = 6000000;//testbench
    static int queries_size = 10;
    static int time_series_size = 256;
    static int paa_segments = 16;
    static int sax_cardinality = 8;
    static int leaf_size = 2000;
    static int min_leaf_size = 10;
    static unsigned long long initial_lbl_size = 2000;
    static unsigned long long flush_limit = 200000;
    static unsigned long long initial_fbl_size = 100;
    static char use_index = 0;
    static int complete_type = 0;
    static int total_loaded_leaves = 1;
    static int tight_bound = 0;
    static int aggressive_check = 0;
    static float minimum_distance = FLT_MAX;
    static int serial_scan = 0;
    static char knnlabel = 0;
    static int min_checked_leaves = -1;
    static int cpu_control_type = 81;
    static char inmemory_flag=0;
    int calculate_thread=8;
    int  function_type =0;
    N_PQUEUE =1;
    maxreadthread=5;
    read_block_length=20000;
    int k_size=0;
    long int labelsize=1;
    int topk=0;
    int dtwwindowsize=0;
    while (1)
    {
        static struct option long_options[] =  {
            {"use-index", no_argument, 0, 'a'},
            {"initial-lbl-size", required_argument, 0, 'b'},
            {"complete-type", required_argument, 0, 'c'},
            {"dataset", required_argument, 0, 'd'},
            {"total-loaded-leaves", required_argument, 0, 'e'},
            {"flush-limit", required_argument, 0, 'f'},
            {"aggressive-check", no_argument, 0, 'g'},
            {"help", no_argument, 0, 'h'},
            {"initial-fbl-size", required_argument, 0, 'i'},
            {"serial", no_argument, 0, 'j'},
            {"queries-size", required_argument, 0, 'k'},
            {"leaf-size", required_argument, 0, 'l'},
            {"min-leaf-size", required_argument, 0, 'm'},
            {"tight-bound", no_argument, 0, 'n'},
            {"read-thread", required_argument, 0, 'o'},
            {"index-path", required_argument, 0, 'p'},
            {"queries", required_argument, 0, 'q'},
            {"read-block", required_argument, 0, 'r'},
            {"minimum-distance", required_argument, 0, 's'},
            {"timeseries-size", required_argument, 0, 't'},
            {"min-checked-leaves", required_argument, 0, 'u'},
            {"in-memory", no_argument, 0, 'v'},
            {"cpu-type", required_argument, 0, 'w'},
            {"sax-cardinality", required_argument, 0, 'x'},
            {"function-type", required_argument, 0, 'y'},
            {"dataset-size", required_argument, 0, 'z'},
            {"k-size", required_argument, 0, '0'},
            {"knn-label-set", required_argument, 0, '1'},
            {"knn-label-size", required_argument, 0, '2'},
            {"knn", no_argument, 0, '3'},
            {"topk", no_argument, 0, '4'},
            {"dtwwindowsize", required_argument, 0, '5'},
            {"queue-number", required_argument, 0, '6'},
            {NULL, 0, NULL, 0}
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;
        int c = getopt_long (argc, argv, "",
                             long_options, &option_index);
        if (c == -1)
            break;
        switch (c)
        {
        	case 'j':
        		serial_scan = 1;
        		break;
            case 'g':
                aggressive_check = 1;
                break;

            case 's':
                minimum_distance = atof(optarg);
                break;

            case 'n':
                tight_bound = 1;
                break;

            case 'e':
                total_loaded_leaves = atoi(optarg);
                break;

            case 'c':
                complete_type = atoi(optarg);
                break;

            case 'q':
                queries = optarg;
                break;

            case 'k':
                queries_size = atoi(optarg);
                break;

            case 'd':
                dataset = optarg;
                break;

            case 'p':
                index_path = optarg;
                break;

            case 'z':
                dataset_size = atoi(optarg);
                break;

            case 't':
                time_series_size = atoi(optarg);
                break;

            case 'x':
                sax_cardinality = atoi(optarg);
                break;

            case 'l':
                leaf_size = atoi(optarg);
                break;

            case 'm':
                min_leaf_size = atoi(optarg);
                break;

            case 'b':
                initial_lbl_size = atoll(optarg);
                break;

            case 'f':
                flush_limit = atoll(optarg);
                break;

            case 'u':
            	min_checked_leaves = atoi(optarg);
            	break;
            case 'w':
                cpu_control_type = atoi(optarg);
                break;

            case 'y':
                function_type = atoi(optarg);
                break;
            case 'i':
                initial_fbl_size = atoll(optarg);
                break;
            case 'o':
                maxreadthread = atoi(optarg);
                break;
            case 'r':
                read_block_length = atoi(optarg);
                break;
            case '0':
                k_size = atoi(optarg); 
                
                break;
            case '1':
                labelset = optarg;
                break;
            case '2':
                labelsize =  atoi(optarg);
               
            case '3':
               knnlabel=1;
               break;
            case '4':
               topk=1;
               break;
            case '5':
                dtwwindowsize =  atoi(optarg);
                break;
            case '6':
                N_PQUEUE =  atoi(optarg);
                break;
            case 'h':
#ifdef BENCHMARK
                printf(PRODUCT);
#endif
                printf("Usage:\n\
                \t--dataset XX \t\t\tThe path to the dataset file\n\
                \t--queries XX \t\t\tThe path to the queries file\n\
                \t--dataset-size XX \t\tThe number of time series to load\n\
                \t--queries-size XX \t\tThe number of queries to do\n\
                \t--minimum-distance XX\t\tThe minimum distance we search (MAX if not set)\n\
                \t--use-index  \t\t\tSpecifies that an input index will be used\n\
                \t--index-path XX \t\tThe path of the output folder\n\
                \t--timeseries-size XX\t\tThe size of each time series\n\
                \t--sax-cardinality XX\t\tThe maximum sax cardinality in number of bits (power of two).\n\
                \t--leaf-size XX\t\t\tThe maximum size of each leaf\n\
                \t--min-leaf-size XX\t\tThe minimum size of each leaf\n\
                \t--initial-lbl-size XX\t\tThe initial lbl buffer size for each buffer.\n\
                \t--flush-limit XX\t\tThe limit of time series in memory at the same time\n\
                \t--initial-fbl-size XX\t\tThe initial fbl buffer size for each buffer.\n\
                \t--complete-type XX\t\t0 for no complete, 1 for serial, 2 for leaf\n\
                \t--total-loaded-leaves XX\tNumber of leaves to load at each fetch\n\
                \t--min-checked-leaves XX\t\tNumber of leaves to check at minimum\n\
                \t--tight-bound XX\tSet for tight bounds.\n\
                \t--aggressive-check XX\t\tSet for aggressive check\n\
                \t--serial\t\t\tSet for serial scan\n\
                \t--queue-number\t\tset the number of priority queues\n\
                \t--in-memory\t\t\tSet for in-memory search\n\
                \t--function-type\t\t\tSet for query answering type\n\
                \t\t\tin memory  only index creation: 0\n\
                \t\t\tParIS-TS: 1\n\
                \t\t\tParIS: 2\n\
                \t\t\t\\MESSI-mq: 3\n\
                \t--cpu-type\t\t\tSet for how many cores you want to used and in 1 or 2 cpu\n\
                \t--help\n\n\
                \tCPU type code:\t\t\t21 : 2 core in 1 CPU\n\
                \t\t\t\t\t22 : 2 core in 2 CPUs\n\
                \t\t\t\t\t41 : 4 core in 1 CPU\n\
                \t\t\t\t\t42 : 4 core in 2 CPUs\n\
                \t\t\t\t\t61 : 6 core in 1 CPU\n\
                \t\t\t\t\t62 : 6 core in 2 CPUs\n\
                \t\t\t\t\t81 : 8 core in 1 CPU\n\
                \t\t\t\t\t82 : 8 core in 2 CPUs\n\
                \t\t\t\t\t101: 10 core in 1 CPU\n\
                \t\t\t\t\t102: 10 core in 2 CPUs\n\
                \t\t\t\t\t121: 12 core in 1 CPU\n\
                \t\t\t\t\t122: 12 core in 2 CPUs\n\
                \t\t\t\t\t181: 18 core in 1 CPU\n\
                \t\t\t\t\t182: 18 core in 2 CPUs\n\
                \t\t\t\t\t242: 24 core in 2 CPUs\n\
                \t\t\t\t\t562: 56 core in 2 CPUs\n\
                \t\t\t\t\tOther: 1 core in 1 CPU\n\
                ");
                return 0;
                break;
            case 'a':
                use_index = 1;
                break;
            case 'v':
                inmemory_flag = 1;
                break;
            default:
                exit(-1);
                break;
        }
    }
    INIT_STATS();
	cpu_set_t mask,get;
    CPU_ZERO(&mask);
    CPU_ZERO(&get);
    if(cpu_control_type==21)
    {
        CPU_SET(0, &mask);
        CPU_SET(2, &mask);
        calculate_thread=2;
        maxquerythread=2;
    }
    else if(cpu_control_type==22)
    {
        CPU_SET(0, &mask);
        CPU_SET(1, &mask);
        calculate_thread=2;
        maxquerythread=2;
    }
    else if(cpu_control_type==41)
    {
        CPU_SET(0, &mask);
        CPU_SET(2, &mask);
        CPU_SET(4, &mask);
        CPU_SET(6, &mask);
	    calculate_thread=4;
        maxquerythread=4;
    }
    else if(cpu_control_type==42)
    {
        CPU_SET(0, &mask);
        CPU_SET(1, &mask);
        CPU_SET(2, &mask);
        CPU_SET(3, &mask);
        calculate_thread=4;
        maxquerythread=4;
    }
    else if(cpu_control_type==61)
    {
        CPU_SET(0, &mask);
        CPU_SET(2, &mask);
        CPU_SET(4, &mask);
        CPU_SET(6, &mask);
        CPU_SET(8, &mask);
        CPU_SET(10, &mask);
        calculate_thread=6;
        maxquerythread=6;
    }
    else if(cpu_control_type==62)
    {
        CPU_SET(0, &mask);
        CPU_SET(1, &mask);
        CPU_SET(2, &mask);
        CPU_SET(3, &mask);
        CPU_SET(4, &mask);
        CPU_SET(5, &mask);
        calculate_thread=6;
        maxquerythread=6;
    }
    else if (cpu_control_type==81)
    {
        CPU_SET(0, &mask);
        CPU_SET(2, &mask);
        CPU_SET(4, &mask);
        CPU_SET(6, &mask);
        CPU_SET(8, &mask);
        CPU_SET(10, &mask);
        CPU_SET(12, &mask);
        CPU_SET(14, &mask);
        calculate_thread=8;
        maxquerythread=8;
    }
    else if (cpu_control_type==82)
    {
        CPU_SET(0, &mask);
        CPU_SET(1, &mask);
        CPU_SET(2, &mask);
        CPU_SET(3, &mask);
        CPU_SET(4, &mask);
        CPU_SET(5, &mask);
        CPU_SET(6, &mask);
        CPU_SET(7, &mask);
        calculate_thread=8;
        maxquerythread=8;
    }
    else if (cpu_control_type==101)
    {
        CPU_SET(0, &mask);
        CPU_SET(2, &mask);
        CPU_SET(4, &mask);
        CPU_SET(6, &mask);
        CPU_SET(8, &mask);
        CPU_SET(10, &mask);
        CPU_SET(12, &mask);
        CPU_SET(14, &mask);
        CPU_SET(16, &mask);
        CPU_SET(18, &mask);
        calculate_thread=10;
        maxquerythread=10;
    }
    else if (cpu_control_type==102)
    {
        CPU_SET(0, &mask);
        CPU_SET(1, &mask);
        CPU_SET(2, &mask);
        CPU_SET(3, &mask);
        CPU_SET(4, &mask);
        CPU_SET(5, &mask);
        CPU_SET(6, &mask);
        CPU_SET(7, &mask);
        CPU_SET(8, &mask);
        CPU_SET(9, &mask);
        calculate_thread=10;
        maxquerythread=10;
    }
    else if (cpu_control_type==121)
    {
        CPU_SET(0, &mask);
        CPU_SET(2, &mask);
        CPU_SET(4, &mask);
        CPU_SET(6, &mask);
        CPU_SET(8, &mask);
        CPU_SET(10, &mask);
        CPU_SET(12, &mask);
        CPU_SET(14, &mask);
        CPU_SET(16, &mask);
        CPU_SET(18, &mask);
        CPU_SET(20, &mask);
        CPU_SET(22, &mask);
        calculate_thread=12;
        maxquerythread=12;
    }
    else if (cpu_control_type==122)
    {
        CPU_SET(0, &mask);
        CPU_SET(1, &mask);
        CPU_SET(2, &mask);
        CPU_SET(3, &mask);
        CPU_SET(4, &mask);
        CPU_SET(5, &mask);
        CPU_SET(6, &mask);
        CPU_SET(7, &mask);
        CPU_SET(8, &mask);
        CPU_SET(9, &mask);
        CPU_SET(10, &mask);
        CPU_SET(11, &mask);
        calculate_thread=12;
        maxquerythread=12;
    }
    else if (cpu_control_type==182)
    {
            CPU_SET(0, &mask);
            CPU_SET(1, &mask);
            CPU_SET(2, &mask);
            CPU_SET(3, &mask);
            CPU_SET(4, &mask);
            CPU_SET(5, &mask);
            CPU_SET(6, &mask);
            CPU_SET(7, &mask);
            CPU_SET(8, &mask);
            CPU_SET(9, &mask);
            CPU_SET(10, &mask);
            CPU_SET(11, &mask);
            CPU_SET(12, &mask);
            CPU_SET(13, &mask);
            CPU_SET(14, &mask);
            CPU_SET(15, &mask);
            CPU_SET(16, &mask);
            CPU_SET(17, &mask);
            calculate_thread=18;
            maxquerythread=18;
    }
    else if (cpu_control_type==242)
    {
            CPU_SET(0, &mask);
            CPU_SET(1, &mask);
            CPU_SET(2, &mask);
            CPU_SET(3, &mask);
            CPU_SET(4, &mask);
            CPU_SET(5, &mask);
            CPU_SET(6, &mask);
            CPU_SET(7, &mask);
            CPU_SET(8, &mask);
            CPU_SET(9, &mask);
            CPU_SET(10, &mask);
            CPU_SET(11, &mask);
            CPU_SET(12, &mask);
            CPU_SET(13, &mask);
            CPU_SET(14, &mask);
            CPU_SET(15, &mask);
            CPU_SET(16, &mask);
            CPU_SET(17, &mask);
            CPU_SET(18, &mask);
            CPU_SET(19, &mask);
            CPU_SET(20, &mask);
            CPU_SET(21, &mask);
            CPU_SET(22, &mask);
            CPU_SET(23, &mask);
            calculate_thread=24;
            maxquerythread=24;
    }
    else if (cpu_control_type==562)
    {
            CPU_SET(0, &mask);
            CPU_SET(1, &mask);
            CPU_SET(2, &mask);
            CPU_SET(3, &mask);
            CPU_SET(4, &mask);
            CPU_SET(5, &mask);
            CPU_SET(6, &mask);
            CPU_SET(7, &mask);
            CPU_SET(8, &mask);
            CPU_SET(9, &mask);
            CPU_SET(10, &mask);
            CPU_SET(11, &mask);
            CPU_SET(12, &mask);
            CPU_SET(13, &mask);
            CPU_SET(14, &mask);
            CPU_SET(15, &mask);
            CPU_SET(16, &mask);
            CPU_SET(17, &mask);
            CPU_SET(18, &mask);
            CPU_SET(19, &mask);
            CPU_SET(20, &mask);
            CPU_SET(21, &mask);
            CPU_SET(22, &mask);
            CPU_SET(23, &mask);
            CPU_SET(24, &mask);
            CPU_SET(25, &mask);
            CPU_SET(26, &mask);
            CPU_SET(27, &mask);
            CPU_SET(28, &mask);
            CPU_SET(29, &mask);
            CPU_SET(30, &mask);
            CPU_SET(31, &mask);
            CPU_SET(32, &mask);
            CPU_SET(33, &mask);
            CPU_SET(34, &mask);
            CPU_SET(35, &mask);
            CPU_SET(36, &mask);
            CPU_SET(37, &mask);
            CPU_SET(38, &mask);
            CPU_SET(39, &mask);
            CPU_SET(40, &mask);
            CPU_SET(41, &mask);
            CPU_SET(42, &mask);
            CPU_SET(43, &mask);
            CPU_SET(44, &mask);
            CPU_SET(45, &mask);
            CPU_SET(46, &mask);
            CPU_SET(47, &mask);
            CPU_SET(48, &mask);
            CPU_SET(49, &mask);
            CPU_SET(50, &mask);
            CPU_SET(51, &mask);
            CPU_SET(52, &mask);
            CPU_SET(53, &mask);
            CPU_SET(54, &mask);
            CPU_SET(55, &mask);

            calculate_thread=56;
            maxquerythread=56;
    }
    else if(cpu_control_type==1)
    {

        calculate_thread=1;
        maxquerythread=1;
    }
    else
    {
        calculate_thread=cpu_control_type;
        maxquerythread=cpu_control_type;

        for (int i = 0; i < cpu_control_type; i++)
        {
            //CPU_SET(i, &mask);
        }
    }


    if (pthread_setaffinity_np(pthread_self(), sizeof(mask), &mask) < 0) 
    {
        fprintf(stderr, "set thread affinity failed\n");
    }

    if (pthread_getaffinity_np(pthread_self(), sizeof(get), &get) < 0) 
    {
        fprintf(stderr, "get thread affinity failed\n");
    }
    if (use_index) 
    {
    	isax_index *idx = index_read(index_path);
    	idx->settings->tight_bound = tight_bound;
    	idx->settings->aggressive_check = aggressive_check;
    	idx->settings->total_loaded_leaves = total_loaded_leaves;
        idx->settings->min_leaf_size = min_leaf_size;
    	print_settings(idx->settings);
        //fprintf(stderr,"total_records: %ld\n", idx->total_records);
        //fprintf(stderr,"loaded_records: %ld\n", idx->loaded_records);

    	//create_wedges(idx, NULL);
 
    	char sanity_test = 0;
    	if(sanity_test)
        {
    		cache_sax_file(idx);
            isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves, &sanity_check_query);
		}
    	else 
        {
    		if(serial_scan) 
            {
    			cache_sax_file(idx);
                if(knnlabel)
                {
                    if(function_type==0)
                    {
                        isax_knn_query_binary_file(queries,labelset, queries_size, idx, minimum_distance, min_checked_leaves,k_size,labelsize, &exact_topk_serial); //ADS+ knn
                    }
                    else if(function_type==1)
                    {
                        isax_knn_query_binary_file(queries,labelset, queries_size, idx, minimum_distance, min_checked_leaves,k_size,labelsize, &exact_topk_serial_ParIS); //ParIS knn
                    }

                }
                else if(topk)
                {
                    if(function_type==0)
                    {
                        isax_topk_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves,k_size, &exact_topk_serial);//ADS+ topk
                    }
                    else if(function_type==1)
                    {
                        isax_topk_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves,k_size, &exact_topk_serial_ParIS);//ParIS topk
                    }
                }
                else{
                    if(function_type==0)
                    {
                        isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_serial);//ADS+ similarity search
                    }
                    else if(function_type==1)
                    {
                        isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_serial_ParIS);//ParIS similarity search
                    }
                else if(function_type==2)
                    {
                        isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_serial_ParIS_nb);//ParIS-nb similarity search
                    }
                    else if(function_type==3)
                    {
                    //isax_query_binary_file_para(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_serial);
                    }
                }

			}
			else
            {
                if(knnlabel)
                {
                    if(function_type==0)
                    {
                        isax_knn_query_binary_file(queries,labelset, queries_size, idx, minimum_distance, min_checked_leaves,k_size,labelsize, &exact_topk);
                    }

                }
                else if(topk)
                {
                    if(function_type==0)
                    {
                        isax_topk_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves,k_size, &exact_topk);
                    }
                }
                else
                {
                if(function_type==0)
                {
                    isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search);
                }
                else
                {
                    isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_m);
                }
                }
			}
    	}
        PRINT_STATS(0.00f)
    	flush_all_leaf_buffers(idx, TMP_ONLY_CLEAN);
    	index_write(idx);

    	//clear_wedges(idx, NULL);
        //printf("this is the total search time is %f\n", total_time);
    	isax_index_destroy(idx, NULL);
    }
    else
    {
        char rm_command[256];

        if(!inmemory_flag)
    	{
            sprintf(rm_command, "rm -rf %s", index_path);
    	   system(rm_command);
           
        }


    	isax_index_settings * index_settings = isax_index_settings_init(index_path,         // INDEX DIRECTORY
    	                                                                    time_series_size,   // TIME SERIES SIZE
    	                                                                    paa_segments,       // PAA SEGMENTS
    	                                                                    sax_cardinality,    // SAX CARDINALITY IN BITS
    	                                                                    leaf_size,          // LEAF SIZE
    	                                                                    min_leaf_size,      // MIN LEAF SIZE
    	                                                                    initial_lbl_size,   // INITIAL LEAF BUFFER SIZE
    	                                                                    flush_limit,        // FLUSH LIMIT
    	                                                                    initial_fbl_size,   // INITIAL FBL BUFFER SIZE
    	                                                                    total_loaded_leaves,// Leaves to load at each fetch
    																		tight_bound,		// Tightness of leaf bounds
    																		aggressive_check,	// aggressive check
    																		1,inmemory_flag);					// new index
    	
        
        if(!inmemory_flag)
        {
            idx = isax_index_init(index_settings);
            print_settings(idx->settings);
        }
        else
        {
            idx = isax_index_init_inmemory(index_settings);
            print_settings(idx->settings);
        }
    	
        #ifdef CLUSTERED
        char s[255];
        sprintf(s, "rm -rf %s.*", dataset);
        system(s);
        #endif





        COUNT_TOTAL_TIME_START

        //// #### COMMENTED OUT: BUGGY CODE!!! ####
        /*int sorted = 0;
        if(sorted == 1) {
        	isax_sorted_index_binary_file(dataset, dataset_size, idx);
        }
        else if(sorted == 2) {
        	isax_merge_sorted_index_binary_file(dataset, dataset_size, idx);
        }
        else {*/
        /// ########################################

        if (inmemory_flag)
        {
            // MESSI: paralllel in memory index creation 
            index_creation_pRecBuf_new(dataset, dataset_size, idx);
            //index_creation_gpu(dataset, dataset_size, idx);
            //dINIT_STATS()
            if(knnlabel)
            {
                if(dtwwindowsize==0)
                {
                    if(function_type==1)
                    {
                        isax_knn_query_binary_file(queries,labelset, queries_size, idx, minimum_distance, min_checked_leaves,k_size,2000, &exact_search_serial_topk_inmemory);
                    }
                    else if(function_type==4)
                    {
                        isax_knn_query_binary_file(queries,labelset, queries_size, idx, minimum_distance, min_checked_leaves,k_size,2000, &exact_topk_serial_ParIS_inmemory);
                    }                
                    else if(function_type==7)
                    {
                       isax_knn_query_binary_file_traditional(queries,labelset, queries_size, idx, minimum_distance, min_checked_leaves,k_size,2000, &exact_topk_MESSImq_inmemory);
                    }
                }
                else
                {
                    if(function_type==5)
                    {
                        isax_DTWknn_query_binary_file(queries,labelset, queries_size, idx, minimum_distance, min_checked_leaves,dtwwindowsize,k_size,2000, &exact_DTWknn_serial_ParIS_inmemory);
                    }                
                    else if(function_type==8)
                    {
                        isax_DTWknn_query_binary_file_traditional(queries,labelset, queries_size, idx, minimum_distance, min_checked_leaves,dtwwindowsize,k_size,2000, &exact_DTWknn_MESSI_inmemory_hybrid);
                    }
                }

            }
            else if(dtwwindowsize==0)
            {
                if(function_type==0)
                {
                    //isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_inmemory);
                    isax_DTWquery_binary_file_traditional(queries, queries_size, idx, minimum_distance, min_checked_leaves,dtwwindowsize);
                }
                else if(function_type==1)
                {
                    isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_parads_inmemory);
                }
                else if(function_type==4)
                {
                    //isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_inmemory_openmp);
                }
                else if(function_type==7)
                {
                    //isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_serial_ParIS_nb_inmemory);
                }
                else if(function_type==2)
                {
                    isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_serial_ParIS_inmemory);
                }
                else if(function_type==5)
                {
                    //bf=1;
                    //isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_serial_ParIS2_inmemory);
                    //isax_query_binary_file_traditional(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_ParISnew_inmemory_hybrid_workstealing);
                }
                else if(function_type==6)
                {
                    isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_serial_ParIS_openmp_inmemory);
                }
                else if(function_type==3)
                {
                    isax_query_binary_file_traditional(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_ParISnew_inmemory_hybrid);
                    //isax_query_binary_file_batch(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_serial_ParIS_nb_batch_inmemory);
                }
                else if(function_type==8)
                {
                    //isax_query_binary_file_traditional(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_ParISnew_inmemory);
                    //isax_query_binary_file_traditional(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_ParISnew_inmemory_workstealing);
                }
                else if(function_type==9)
                {
                   // isax_query_binary_file_traditional(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_serial_ParGIS_openmp_inmemory);
                }
                else if(function_type==10)
                {
                    
                    //index_generate_inmemory_pRecBuf(dataset, dataset_size, idx);
                    //COUNT_QUEUE_TIME_START
                    //flush_pRecBuf_inmemory((parallel_first_buffer_layer*) idx->fbl, idx);
                    //COUNT_QUEUE_TIME_END
                    //PRINT_STATS(0.0f)

                    //printf("this is the queue time: %f\n", total_queue_time);
                    //INIT_STATS()
                    //isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_serial_inmemory);
                }
            }
            else
            {
                if(function_type==0)
                {                    
                    isax_DTWquery_binary_file_traditional(queries, queries_size, idx, minimum_distance, min_checked_leaves,dtwwindowsize);
                }
            }
            
            PRINT_STATS(0.00f)

        }
        else
        {   
            // ParIS/ParIS+ on disk index creation
            if(function_type==0)
            {
                // ADS+ index creation
                isax_index_binary_file(dataset, dataset_size, idx);

                flush_fbl(idx->fbl, idx);
            }
            else if(function_type==1)
            {
                isax_index_binary_file_m(dataset, dataset_size, idx,calculate_thread);// ParIS indexing program in parallel
            }
            else if(function_type==2)
            {
                isax_index_binary_file_m_new(dataset, dataset_size, idx,calculate_thread);//ParIS+ indexing program in parallel
            }
        }
        
        
        float distance = 0;
        FLUSHES++;
        
        COUNT_TOTAL_TIME_END
        //PRINT_STATS(distance)

        COUNT_TOTAL_TIME_START
        if(complete_type == 1) {
        	fprintf(stderr,">>> Completing index.\n");
        	complete_index(idx, dataset_size);
        }
        else if(complete_type == 2) {
        	fprintf(stderr,">>> Completing index.\n");
            complete_index_leafs(idx);
        }
        COUNT_TOTAL_TIME_END


        if(!inmemory_flag)
        {
            index_write(idx); 
            isax_index_destroy(idx, NULL);
            PRINT_STATS(distance)

        }
        else
        {
            free(rawfile);
            //if(function_type==10)
            {
                isax_index_pRecBuf_destroy(idx, NULL,maxquerythread);
            }
            //else
            //{
                //isax_index_destroy(idx, NULL);
                //MESSI2_index_destroy(idx, NULL);
            //}
        }

        
	    
	    printf("\n");

    }

    return 0;
}
        
void  INThandler(int sig)
{
     char  c;

     signal(sig, SIG_IGN);
     fprintf(stderr, "Do you really want to quit? [y/n] ");
     c = getchar();
     if (c == 'y' || c == 'Y') {
    	 c = getchar();
    	 fprintf(stderr, "Do you want to save the index? [y/n] ");
    	 c = getchar();
    	 if (c == 'y' || c == 'Y') {
    		 flush_fbl(idx->fbl, idx);
    		 index_write(idx);
    	 }
    	 exit(0);
     }
     else
    signal(SIGINT, INThandler);
    getchar(); // Get new line character
}
