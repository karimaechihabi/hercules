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


################################################
#Figure 11a: Scalability with Increasing k (SALD)
################################################
CURR_DATASET=data_size899M_sald_len128_znorm 
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 200000000 61440  128 200000 12 24 12 120000;
#$SCRIPTS/index_experiment.sh dstree    ${CURR_DATASET} 200000000 61440  128 200000 1 1 1 0;
#$SCRIPTS/index_experiment.sh vaplus    ${CURR_DATASET} 200000000 61440  128 200000 1 1 1 0;
#$SCRIPTS/index_experiment.sh paris     ${CURR_DATASET} 200000000 20480  128 2000 1 1 1 0;
CURR_QUERIES=queries_hard5p_sald_len128_znorm
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 5 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 10 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 25 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 100 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 5 1 1 1 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 10 1 1 1 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 25 1 1 1 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 100 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 5 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 10 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 25 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 100 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  20480  128 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  20480  128 2000 1 0 5 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  20480  128 2000 1 0 10 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  20480  128 2000 1 0 25 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  20480  128 2000 1 0 100 1 1 1 1;

####################################################
#Figure 11b: Scalability with Increasing k (SEISMIC)
####################################################
CURR_DATASET=data_size100M_seismic_len256_znorm
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 100000000 61440  256 100000 12 24 12 120000;
#$SCRIPTS/index_experiment.sh dstree    ${CURR_DATASET} 100000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/index_experiment.sh vaplus    ${CURR_DATASET} 100000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/index_experiment.sh paris     ${CURR_DATASET} 100000000 20480  256 2000 1 1 1 0;
CURR_QUERIES=queries_hard5p_seismic_len256_znorm
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 5 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 10 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 25 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 100 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 5 1 1 1 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 10 1 1 1 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 25 1 1 1 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 100 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 5 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 10 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 25 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 100 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  20480  256 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  20480  256 2000 1 0 5 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  20480  256 2000 1 0 10 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  20480  256 2000 1 0 25 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  20480  256 2000 1 0 100 1 1 1 1;

####################################################
#Figure 11c: Scalability with Increasing k (DEEP)
####################################################
CURR_DATASET=data_size1B_deep1b_len96_znorm
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 266666667 61440  96 266667 12 24 12 120000;
#$SCRIPTS/index_experiment.sh dstree    ${CURR_DATASET} 266666667 61440  96 266667 1 1 1 0;
#$SCRIPTS/index_experiment.sh vaplus    ${CURR_DATASET} 266666667 61440  96 266667 1 1 1 0;
#$SCRIPTS/index_experiment.sh paris     ${CURR_DATASET} 266666667 20480  96 2000 1 1 1 0;
CURR_QUERIES=queries_hard5p_deep1b_len96_znorm
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 5 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 10 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 25 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 100 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 5 1 1 1 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 10 1 1 1 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 25 1 1 1 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 100 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 5 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 10 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 25 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 100 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  20480  96 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  20480  96 2000 1 0 5 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  20480  96 2000 1 0 10 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  20480  96 2000 1 0 25 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  20480  96 2000 1 0 100 1 1 1 1;


################################################
#Figure 11.d : Extrapolate from logs of Fig 11.a 
################################################

################################################
#Figure 11.e : Extrapolate from logs of Fig 11.b 
################################################

################################################
#Figure 11.f : Extrapolate from logs of Fig 11.c 
################################################



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

