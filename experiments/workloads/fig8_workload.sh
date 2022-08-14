##################################
# Test workloads with small dataset
#################################

mkdir -p $EXPERIMENTS/logs    
mkdir -p $EXPERIMENTS/results
mkdir -p $EXPERIMENTS/results_post_processed/logs_csv
mkdir -p $EXPERIMENTS/plots/
mkdir -p $HOME/local/lib/R    

#EXPERIMENTS="/home/karimae/projects/SISS_full/experiments";
#SCRIPTS=$EXPERIMENTS/scripts


#########################################
#Figure 8: Scalability with Series Length
#########################################

#Series Length = 128
CURR_QUERIES=queries_size1K_seed14784_len128_znorm
CURR_DATASET=data_size200M_seed1184_len128_znorm
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 200000000 61440  128 200000 12 24 12 120000;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/index_experiment.sh dstree    ${CURR_DATASET} 200000000 61440  128 200000 1 1 1 0;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh vaplus    ${CURR_DATASET} 200000000 61440  128 200000 1 1 1 0;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh paris     ${CURR_DATASET} 200000000 20480  128 2000 1 1 1 0;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  20480  128 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  120000  128 1 1 0 1 1 1 2 1;

#Series Length = 256
CURR_QUERIES=queries_size1K_seed14784_len256_znorm
CURR_DATASET=data_size2B_seed1184_len256_znorm
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 100000000 61440  256 100000 12 24 12 120000;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/index_experiment.sh dstree ${CURR_DATASET} 100000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/query_experiment.sh dstree ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh vaplus    ${CURR_DATASET} 100000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh paris     ${CURR_DATASET} 100000000 20480  256 2000 1 1 1 0;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  20480  256 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  120000  256 1 1 0 1 1 1 2 1;

#Series Length = 512
CURR_QUERIES=queries_size1K_seed14784_len512_znorm
CURR_DATASET=data_size50M_seed1184_len512_znorm
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 50000000 61440  512 50000 12 24 12 60000;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 50000000  61440  512 50000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/index_experiment.sh dstree    ${CURR_DATASET} 50000000 61440  512 50000 1 1 1 0;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 50000000  61440  512 50000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh vaplus    ${CURR_DATASET} 50000000 61440  512 50000 1 1 1 0;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 50000000  61440  512 50000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh paris     ${CURR_DATASET} 50000000 20480  512 2000 1 1 1 0;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 50000000  20480  512 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 50000000  60000  512 1 1 0 1 1 1 2 1;

#Series Length = 1024
CURR_QUERIES=queries_size1K_seed14784_len1024_znorm
CURR_DATASET=data_size25M_seed1184_len1024_znorm
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 25000000 61440  1024 25000 12 24 12 30000;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 25000000  61440  1024 25000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/index_experiment.sh dstree    ${CURR_DATASET} 25000000 61440  1024 25000 1 1 1 0;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 25000000  61440  1024 25000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh vaplus    ${CURR_DATASET} 25000000 61440  1024 25000 1 1 1 0;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 25000000  61440  1024 25000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh paris     ${CURR_DATASET} 25000000 20480  1024 2000 1 1 1 0;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 25000000  20480  1024 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 25000000  30000  1024 1 1 0 1 1 1 2 1;


#Series Length = 2048
CURR_QUERIES=queries_size1K_seed15631_len2048_znorm
CURR_DATASET=data_size12M500K_seed1184_len2048_znorm
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 12500000 61440  2048 12500 12 24 12 15000;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 12500000  61440  2048 12500 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/index_experiment.sh dstree    ${CURR_DATASET} 12500000 61440  2048 12500 1 1 1 0;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 12500000  61440  2048 12500 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh vaplus    ${CURR_DATASET} 12500000 61440  2048 12500 1 1 1 0;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 12500000  61440  2048 12500 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh paris     ${CURR_DATASET} 12500000 20480  2048 2000 1 1 1 0;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 12500000  20480  2048 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 12500000  15000 2048 1 1 0 1 1 1 2 1;


#Series Length = 4096
CURR_QUERIES=queries_size1K_seed14784_len4096_znorm
CURR_DATASET=data_size6M250K_seed1184_len4096_znorm
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 6250000 61440 4096 6250 12 24 12 7500;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 6250000  61440  4096 6250 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/index_experiment.sh dstree    ${CURR_DATASET} 6250000 61440  4096 6250 1 1 1 0;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 6250000  61440  4096 6250 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh vaplus    ${CURR_DATASET} 6250000 61440 4096 6250 1 1 1 0;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 6250000  61440  4096 6250 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh paris     ${CURR_DATASET} 6250000 20480  4096 2000 1 1 1 0;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 6250000  20480  4096 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 6250000  7500 4096 1 1 0 1 1 1 2 1;


#Series Length = 8192
CURR_QUERIES=queries_size1K_seed14784_len8192_znorm
CURR_DATASET=data_size3M125K_seed1184_len8192_znorm
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 3125000 61440 8192 3125 12 24 12 3500;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 3125000  61440  8192 3125 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/index_experiment.sh dstree    ${CURR_DATASET} 3125000 61440  8192 3125 1 1 1 0;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 3125000  61440  8192 3125 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh vaplus    ${CURR_DATASET} 3125000 61440 8192 3125 1 1 1 0;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 3125000  61440  8192 3125 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh paris     ${CURR_DATASET} 3125000 20480  8192 2000 1 1 1 0;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 3125000  20480  8192 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 3125000  3500 8192 1 1 0 1 1 1 2 1;


#Series Length = 16384
CURR_QUERIES=queries_size1K_seed14784_len16384_znorm
CURR_DATASET=data_size1M62K500_seed1184_len16384_znorm
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 1562500 61440 16384 1562 12 24 12 1500;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 1562500  61440  16384 1562 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/index_experiment.sh dstree    ${CURR_DATASET} 1562500 61440  16384 12500 1 1 1 0;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 1562500  61440  16384 1562 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh vaplus    ${CURR_DATASET} 1562500 61440 16384 12500 1 1 1 0;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 1562500  61440  16384 1562 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh paris     ${CURR_DATASET} 1562500 20480  16384 2000 1 1 1 0;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 1562500  20480  16384 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 1562500  1500 16384 1 1 0 1 1 1 2 1;



CURR_DIR=$PWD


#cd $EXPERIMENTS/results_post_processed/logs_csv/
#mkdir -p fig2/$VERSION/
#cp *_25000000*.* indexes/$VERSION/
#mv *.* fig2/$VERSION/
#cd fig2/$VERSION/

#Rscript $EXPERIMENTS/scripts/plot.R fig2;

#cd $EXPERIMENTS/logs/
#mkdir -p fig2
#mv *.* fig2/

cd $CURR_DIR;

