#!/bin/bash

function index_experiment ()
{
    CURR_DIR=$PWD

    ALGORITHM=$1;
    DATASET=$2;
    DATASET_SIZE=$3;
    BUFFER_SIZE=$4;
    DATASERIES_SIZE=$5;
    LEAF_SIZE=$6;
    WTHREADS=$7;
    THREADS=$8;
    FLUSH_THRESHOLD=$9;
    DBSIZE=${10};

    SCRIPTS="scripts";
    EXPERIMENTS="/home/karimae/projects/SISS_full/experiments";

    case ${ALGORITHM} in
            "dstree")
	          CODE="DSTREE"
                  ;;
            "hercules")
	          CODE="HERCULES"
                  ;;
            "paris")
	          CODE="PARIS"
                  ;;
            "pscan")
	          CODE="PSCAN"
                  ;;
            "vaplus")
	          CODE="VAPLUS"
                  ;;
            *)		
                  echo "Unknown algorithm in index_experiment.sh"
		  exit -1
                  ;;
    esac
	RESULTS_PATH=${EXPERIMENTS}/results/${DATASET}/${DATASET_SIZE}/${BUFFER_SIZE}/${DATASERIES_SIZE}/${LEAF_SIZE}/${ALGORITHM}/;
    RESULTS_LOG_FILE=${CODE}_${DATASET}_${DATASET_SIZE}_${BUFFER_SIZE}_${DATASERIES_SIZE}_${LEAF_SIZE}_${THREADS}_${FLUSH_THRESHOLD}_${WTHREADS}_${DBSIZE}_index.log;
    RESULTS_ERR_FILE=${CODE}_${DATASET}_${DATASET_SIZE}_${BUFFER_SIZE}_${DATASERIES_SIZE}_${LEAF_SIZE}_${THREADS}_${FLUSH_THRESHOLD}_${WTHREADS}_${DBSIZE}_index.err;

    mkdir -p $RESULTS_PATH/logs

    sudo sync;
    sudo sh -c "echo 3 >> /proc/sys/vm/drop_caches";

    $EXPERIMENTS/$SCRIPTS/run_index_experiment.sh  $ALGORITHM ${DATASET} $DATASET_SIZE  $BUFFER_SIZE $DATASERIES_SIZE $LEAF_SIZE $WTHREADS $THREADS $FLUSH_THRESHOLD $DBSIZE $RESULTS_PATH $RESULTS_LOG_FILE $RESULTS_ERR_FILE
    
    #cd  $EXPERIMENTS/logs

    #python $EXPERIMENTS/$SCRIPTS/post_process.py $RESULTS_LOG_FILE;

    cd $CURR_PWD    
}

index_experiment $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10}

