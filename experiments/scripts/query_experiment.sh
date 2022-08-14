#!/bin/bash


function query_experiment ()
{
    CURR_DIR=$PWD

    ALGORITHM=$1
    QUERIES=$2; 
    QUERIES_SIZE=$3;
    DATASET=$4;	
    DATASET_SIZE=$5; 
    BUFFER_SIZE=$6;        
    DATASERIES_SIZE=$7;
    LEAF_SIZE=$8;
    DELTA=$9;
    EPSILON=${10};
    K=${11};
    EAPCA_THRESHOLD=${12};
    SAX_THRESHOLD=${13};
    QTHREADS=${14};
    RTHREADS=${15};


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
                  echo "Unknown algorithm in query_experiment.sh"
		  exit -1
                  ;;
    esac

    RESULTS_PATH=${EXPERIMENTS}/results/${DATASET}/${DATASET_SIZE}/${BUFFER_SIZE}/${DATASERIES_SIZE}/${LEAF_SIZE}/${ALGORITHM};
    RESULTS_LOG_FILE=${CODE}_${DATASET}_${DATASET_SIZE}_${BUFFER_SIZE}_${DATASERIES_SIZE}_${LEAF_SIZE}_${QUERIES}_${DELTA}_${EPSILON}_${K}_${EAPCA_THRESHOLD}_${SAX_THRESHOLD}_${QTHREADS}_${RTHREADS}_${QUERIES_SIZE}_search.log;
    RESULTS_ERR_FILE=${CODE}_${DATASET}_${DATASET_SIZE}_${BUFFER_SIZE}_${DATASERIES_SIZE}_${LEAF_SIZE}_${QUERIES}_${DELTA}_${EPSILON}_${K}_${EAPCA_THRESHOLD}_${SAX_THRESHOLD}_${QTHREADS}_${RTHREADS}_${QUERIES_SIZE}_search.err;
 
    mkdir -p $RESULTS_PATH/logs

   
    sudo sync;
    sudo sh -c "echo 3 >> /proc/sys/vm/drop_caches";


    $EXPERIMENTS/scripts/run_query_experiment.sh  $ALGORITHM ${QUERIES} $QUERIES_SIZE $DATASET $DATASET_SIZE $BUFFER_SIZE $DATASERIES_SIZE $LEAF_SIZE $DELTA $EPSILON $K ${EAPCA_THRESHOLD} ${SAX_THRESHOLD} ${QTHREADS} ${RTHREADS} $RESULTS_PATH $RESULTS_LOG_FILE $RESULTS_ERR_FILE


    cd $CURR_PWD    
}

query_experiment $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10}  ${11} ${12} ${13}  ${14} ${15} 
