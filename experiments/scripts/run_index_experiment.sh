#!/bin/bash

function run_index_experiment ()
{
     
    if [ "$#" -eq "13" ]; then
        ALGORITHM=$1
		DATASET=$2;
		DATASET_SIZE=$3;
		BUFFER_SIZE=$4;
	
		DATASERIES_SIZE=$5;
		LEAF_SIZE=$6;
		WTHREADS=$7; 
		THREADS=$8;
		FLUSH_THRESHOLD=$9;
	    DBSIZE=${10};

		RESDIR=${11};
   		LOG_FILE=${12};
	    ERR_FILE=${13};	

        SCRIPTS="scripts";
        EXPERIMENTS="/home/karimae/projects/SISS_full/experiments";
        DATASETS="/home/karimae/data/synthetic/";

		#saving the current directory
    	CURR_DIR=$PWD

		#executing from the bin directory
		#BIN_DIRECTORY=${EXPERIMENTS}/bin/$ALGORITHM
		#cd $BIN_DIRECTORY
   		cd $EXPERIMENTS         

	 	#calculate the flush limit for ADS and ISAX2+ based on buffer size value
        let "FLUSH_LIMIT=$BUFFER_SIZE*1024*1024/4/$DATASERIES_SIZE"
        let "TS_FFT_TOTAL_BYTE_SIZE=$DATASERIES_SIZE*4+16*4"
        let "SAMPLE_SIZE=$BUFFER_SIZE*1024*1024/$TS_FFT_TOTAL_BYTE_SIZE"
        let "SAMPLE_SIZE=$SAMPLE_SIZE/10"



        case ${ALGORITHM} in
          "dstree")
	       CODE="DSTREE";
 		  COMMAND="bin/dstree";
		  ARGS="--dataset ${DATASETS}/$DATASET.bin";
		  ARGS+=" --dataset-size $DATASET_SIZE";
		  ARGS+=" --buffer-size $BUFFER_SIZE";
		  ARGS+=" --leaf-size $LEAF_SIZE";
		  ARGS+=" --index-path $RESDIR/index/";
		  ARGS+=" --ascii-input 0";
		  ARGS+=" --mode 0";
		  ARGS+=" --timeseries-size ${DATASERIES_SIZE}";
                  ;;
           "hercules")
          CODE="HERCULES";
		  COMMAND="bin/hercules";
		  ARGS="--dataset ${DATASETS}/$DATASET.bin";
		  ARGS+=" --dataset-size $DATASET_SIZE";
		  ARGS+=" --buffer-size $BUFFER_SIZE";
		  ARGS+=" --leaf-size $LEAF_SIZE";
		  ARGS+=" --index-path $RESDIR/index/";
		  ARGS+=" --ascii-input 0";
		  ARGS+=" --timeseries-size ${DATASERIES_SIZE}";
		  ARGS+=" --serial 1";
		  ARGS+=" --sims 1";
		  ARGS+=" --num-threads ${THREADS}";
		  ARGS+=" --flush-threshold ${FLUSH_THRESHOLD}";
		  ARGS+=" --num-write-threads ${WTHREADS}";
		  ARGS+=" --initial-db-size $DBSIZE";
          ;;
 
            "vaplus")
	          CODE="VAPLUS";
 		  COMMAND="bin/vaplus";
		  ARGS="--dataset ${DATASETS}/$DATASET.bin";
		  ARGS+=" --dataset-size $DATASET_SIZE";
		  ARGS+=" --buffer-size $BUFFER_SIZE"; #replace this later with formula to calc num of ts out of buff size
		  ARGS+=" --index-path $RESDIR/index/";
		  ARGS+=" --ascii-input 0";
		  ARGS+=" --mode 0";
		  ARGS+=" --timeseries-size ${DATASERIES_SIZE}";
                  ;;

            "paris")
	          CODE="PARIS";
 		  COMMAND="bin/paris";
		  ARGS="--dataset ${DATASETS}/$DATASET.bin";
		  ARGS+=" --dataset-size $DATASET_SIZE";
		  ARGS+=" --flush-limit $FLUSH_LIMIT"; #replace this later with formula to calc num of ts out of buff size
		  ARGS+=" --leaf-size $LEAF_SIZE";
		  ARGS+=" --index-path $RESDIR/index/";
		  #ARGS+=" --serial";
		  #ARGS+=" --sax-cardinality 8";		  
		  ARGS+=" --initial-lbl-size $LEAF_SIZE";		  
		  ARGS+=" --timeseries-size ${DATASERIES_SIZE}";
          ARGS+=" --function-type 1";  
          ARGS+=" --cpu-type 121";  
                  ;;
            *)
                  echo "Unknown algorithm in run_index_experiment.sh"
		  exit -1
                  ;;

        esac

        echo $COMMAND $ARGS; 

        echo $COMMAND $ARGS >> $RESDIR/logs/$ERR_FILE;

     	stdbuf -oL /usr/bin/time --format='TIME_ELAPSED %e' $COMMAND $ARGS >> $RESDIR/logs/$LOG_FILE 2>>$RESDIR/logs/$ERR_FILE;

        echo -n "Index_disk_size_MB = " >> $RESDIR/logs/$LOG_FILE;
        du -ms  $RESDIR/index >> $RESDIR/logs/$LOG_FILE;

        cp  $RESDIR/logs/$LOG_FILE $EXPERIMENTS/logs;
        cp  $RESDIR/logs/$ERR_FILE $EXPERIMENTS/logs;


    else #if wrong args
      echo "Error: wrong args for the run_index_experiment function";
      #exit -1;   	
    fi

    cd $CURR_DIR;
}

run_index_experiment $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} 
