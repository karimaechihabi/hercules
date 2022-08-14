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


####################################################
#Figure 9a: Scalability with Query Difficulty (SALD)
####################################################
CURR_DATASET=data_size899M_sald_len128_znorm 
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 200000000 61440  128 200000 12 24 12 120000;
#$SCRIPTS/index_experiment.sh dstree    ${CURR_DATASET} 200000000 61440  128 200000 1 1 1 0;
#$SCRIPTS/index_experiment.sh vaplus    ${CURR_DATASET} 200000000 61440  128 200000 1 1 1 0;
#$SCRIPTS/index_experiment.sh paris     ${CURR_DATASET} 200000000 20480  128 2000 1 1 1 0;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  120000  128 1 1 0 1 1 1 2 1;
CURR_QUERIES=queries_hard1p_sald_len128_znorm
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  20480  128 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  120000  128 1 1 0 1 1 1 2 1;
CURR_QUERIES=queries_hard2p_sald_len128_znorm
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  20480  128 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  120000  128 1 1 0 1 1 1 2 1;
CURR_QUERIES=queries_hard5p_sald_len128_znorm
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  20480  128 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  120000  128 1 1 0 1 1 1 2 1;
CURR_QUERIES=queries_hard10p_sald_len128_znorm
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  20480  128 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  120000  128 1 1 0 1 1 1 2 1;
CURR_QUERIES=queries_size100_sald_len128_znorm
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  61440  128 200000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  20480  128 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 200000000  120000  128 1 1 0 1 1 1 2 1;

##############################################
#Figure 9.b : Extrapolate from logs of Fig 9.a 
##############################################

#######################################################
#Figure 9c: Scalability with Query Difficulty (SEISMIC)
#######################################################
CURR_DATASET=data_size100M_seismic_len256_znorm
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 100000000 61440  256 100000 12 24 12 120000;
#$SCRIPTS/index_experiment.sh dstree    ${CURR_DATASET} 100000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/index_experiment.sh vaplus    ${CURR_DATASET} 100000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/index_experiment.sh paris     ${CURR_DATASET} 100000000 20480  256 2000 1 1 1 0;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  120000  256 1 1 0 1 1 1 2 1;
CURR_QUERIES=queries_hard1p_seismic_len256_znorm
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  20480  256 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  120000  256 1 1 0 1 1 1 2 1;
CURR_QUERIES=queries_hard2p_seismic_len256_znorm
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  20480  256 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  120000  256 1 1 0 1 1 1 2 1;
CURR_QUERIES=queries_hard5p_seismic_len256_znorm
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  120000  256 1 1 0 1 1 1 2 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  20480  256 2000 1 0 1 1 1 1 1;
CURR_QUERIES=queries_hard10p_seismic_len256_znorm
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  20480  256 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  120000  256 1 1 0 1 1 1 2 1;
CURR_QUERIES=queries_size100_seismic_len256_znorm
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  20480  256 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  120000  256 1 1 0 1 1 1 2 1;


##############################################
#Figure 9.d : Extrapolate from logs of Fig 9.c 
##############################################

#######################################################
#Figure 9e: Scalability with Query Difficulty (SEISMIC)
#######################################################
CURR_DATASET=data_size1B_deep1b_len96_znorm
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 266666667 61440  96 266667 12 24 12 120000;
#$SCRIPTS/index_experiment.sh dstree    ${CURR_DATASET} 266666667 61440  96 266667 1 1 1 0;
#$SCRIPTS/index_experiment.sh vaplus    ${CURR_DATASET} 266666667 61440  96 266667 1 1 1 0;
#$SCRIPTS/index_experiment.sh paris     ${CURR_DATASET} 266666667 20480  96 2000 1 1 1 0;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  120000  96 1 1 0 1 1 1 2 1;
CURR_QUERIES=queries_hard1p_deep1b_len96_znorm
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  20480  96 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  120000  96 1 1 0 1 1 1 2 1;
CURR_QUERIES=queries_hard2p_deep1b_len96_znorm
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  20480  96 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  120000  96 1 1 0 1 1 1 2 1;
CURR_QUERIES=queries_hard5p_deep1b_len96_znorm
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  20480  96 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  120000  96 1 1 0 1 1 1 2 1;
CURR_QUERIES=queries_hard10p_deep1b_len96_znorm
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  20480  96 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  120000  96 1 1 0 1 1 1 2 1;
CURR_QUERIES=queries_orig100_deep1b_len96_znorm
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  20480  96 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  120000  96 1 1 0 1 1 1 2 1;

##############################################
#Figure 9.f : Extrapolate from logs of Fig 9.e 
##############################################



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

