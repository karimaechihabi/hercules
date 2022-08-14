#!/bin/bash


function run_query_experiment ()
{

    if [[ "$#" -eq 18 ]]; then
        ALGORITHM=$1
	QUERIES=$2;
	QUERIES_SIZE=$3;
	DATASET=$4;	
	DATASET_SIZE=$5;
	BUFFER_SIZE=$6;        
	DATASERIES_SIZE=$7;
	LEAF_SIZE=$8;
    DELTA=${9};
    EPSILON=${10};
    K=${11};
    EAPCA_THRESHOLD=${12};
    SAX_THRESHOLD=${13};
    QTHREADS=${14};
    RTHREADS=${15};

	RESDIR=${16};
    LOG_FILE=${17};
    ERR_FILE=${18};
	
	#saving the current directory
        CURR_DIR=$PWD

        SCRIPTS="scripts";
        EXPERIMENTS="/home/karimae/projects/SISS_full/experiments";
        DATASETS="/home/karimae/data/synthetic/";

    cd $EXPERIMENTS

 	#calculate the flush limit for ADS based on buffer size value
        let "FLUSH_LIMIT=$BUFFER_SIZE*1024*1024/4/256"

        let "TS_FFT_TOTAL_BYTE_SIZE=$DATASERIES_SIZE*4+16*4"
        let "SAMPLE_SIZE=$BUFFER_SIZE*1024*1024/$TS_FFT_TOTAL_BYTE_SIZE"
        let "SAMPLE_SIZE=$SAMPLE_SIZE/10"
     
        case ${ALGORITHM} in
          "dstree")
 		  COMMAND="bin/dstree";
		  ARGS=" --queries ${DATASETS}/${QUERIES}.bin";
		  ARGS+=" --queries-size $QUERIES_SIZE";
		  ARGS+=" --buffer-size $BUFFER_SIZE";
		  #ARGS+=" --leaf-size $LEAF_SIZE";
		  ARGS+=" --index-path $RESDIR/index/";
		  ARGS+=" --ascii-input 0";
		  ARGS+=" --mode 1";
		  ARGS+=" --delta $DELTA";
		  ARGS+=" --epsilon $EPSILON";
		  ARGS+=" --k $K";
		  ARGS+=" --timeseries-size ${DATASERIES_SIZE}";
           ;;

          "hercules")
 		  COMMAND="bin/hercules";
		  ARGS=" --queries ${DATASETS}/${QUERIES}.bin";
		  ARGS+=" --queries-size $QUERIES_SIZE";
		  ARGS+=" --buffer-size $BUFFER_SIZE";
		  #ARGS+=" --leaf-size $LEAF_SIZE";
		  ARGS+=" --index-path $RESDIR/index/";
		  ARGS+=" --ascii-input 0";
		  ARGS+=" --mode 1";
		  ARGS+=" --delta $DELTA";
		  ARGS+=" --epsilon $EPSILON";
		  ARGS+=" --k $K";
		  ARGS+=" --serial 1";
		  ARGS+=" --num-threads $QTHREADS";
		  ARGS+=" --num-read-threads $RTHREADS";
		  ARGS+=" --sims 1";
		  ARGS+=" --timeseries-size ${DATASERIES_SIZE}";
		  ARGS+=" --eapca-threshold ${EAPCA_THRESHOLD}";
		  ARGS+=" --sax-threshold ${SAX_THRESHOLD}";
          ;;

          "vaplus")
	      CODE="VAPLUS";
 		  COMMAND="bin/vaplus";
		  ARGS="--dataset ${DATASETS}/$DATASET.bin";
		  ARGS+=" --queries ${DATASETS}/$QUERIES.bin";
		  ARGS+=" --queries-size ${QUERIES_SIZE}";		  
		  ARGS+=" --dataset-size $DATASET_SIZE";
		  ARGS+=" --buffer-size $BUFFER_SIZE";
		  ARGS+=" --index-path $RESDIR/index/";
		  ARGS+=" --mode 1";
		  ARGS+=" --timeseries-size ${DATASERIES_SIZE}";		 
          ARGS+=" --epsilon $EPSILON";
          ARGS+=" --delta $DELTA";  
          ARGS+=" --k $K";  
                  ;;

         "paris")	
 		  CODE="PARIS"    
 		  COMMAND="bin/paris";
		  ARGS="--queries ${DATASETS}/$QUERIES.bin";
		  ARGS+=" --queries-size ${QUERIES_SIZE}";			  
		  #ARGS+=" --flush-limit $FLUSH_LIMIT"; 
		  ARGS+=" --leaf-size $LEAF_SIZE";
		  ARGS+=" --min-leaf-size $LEAF_SIZE";
		  ARGS+=" --index-path $RESDIR/index/";
		  ARGS+=" --serial";
		  ARGS+=" --function-type 1";
		  ARGS+=" --cpu-type 242 ";
		  ARGS+=" --initial-lbl-size $LEAF_SIZE";		  
		  ARGS+=" --timeseries-size ${DATASERIES_SIZE}";
		  ARGS+=" --use-index";		  
          ARGS+=" --k-size $K";  
          ARGS+=" --topk";  
		  ;;

          "pscan")	    
 		  COMMAND="bin/pscan";
		  ARGS=" --queries ${DATASETS}/$QUERIES.bin";
		  ARGS+=" --queries-size ${QUERIES_SIZE}";
		  ARGS+=" --dataset ${DATASETS}/$DATASET.bin";
		  ARGS+=" --dataset-size $DATASET_SIZE";
		  ARGS+=" --timeseries-size ${DATASERIES_SIZE}";
		  ARGS+=" --k $K";
		  ARGS+=" --dbsize $BUFFER_SIZE";
		  ARGS+=" --num-query-threads $QTHREADS";

		  ;;
            *)
            echo "Unknown algorithm in run_query_experiment.sh"
		  exit -1
                  ;;

        esac

        echo $COMMAND $ARGS;

        echo $COMMAND $ARGS >> $RESDIR/logs/$ERR_FILE;

    	stdbuf -oL /usr/bin/time --format='TIME_ELAPSED %e' $COMMAND $ARGS > $RESDIR/logs/$LOG_FILE 2>>$RESDIR/logs/$ERR_FILE;


        cp  $RESDIR/logs/$LOG_FILE $EXPERIMENTS/logs;
        cp  $RESDIR/logs/$ERR_FILE $EXPERIMENTS/logs;

	
    else #if other args, do sthing else
      echo "Error: wrong args for the run_query_experiment function";
      #exit -1;   	
    fi;
    cd $CURR_DIR;
}

run_query_experiment $1 $2 $3 $4 $5  $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14}  ${15} ${16} ${17} ${18} 

