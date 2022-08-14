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


CURR_QUERIES=queries_size1K_seed14784_len256_znorm
CURR_DATASET=data_size2B_seed1184_len256_znorm


####################
#PARAMETERS
#################### 
# HERCULES
#index_script_full_path algorithm dataset dataset_size flush_limit timeseries_size leaf_size write_threads build_threads flush_threshold dbsize
#query_script_full_path algorithm queries queries_size dataset dataset_size flush_limit timeseries_size leaf_size delta epsilon k sax_threshold eapca_threshold query_threads read_threads


########################################################
#Figure 6.a : Scalability with Increasing Dataset Sizes 
########################################################

#1GB Dataset (sanity check)
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 1000000 61440  256 1000 12 24 12 120000;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 1000000  61440  256 1000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/index_experiment.sh  dstree   ${CURR_DATASET} 1000000 61440  256 1000 1 1 1 0;
#$SCRIPTS/query_experiment.sh  dstree   ${CURR_QUERIES} 100 ${CURR_DATASET} 1000000  61440  256 1000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh vaplus    ${CURR_DATASET} 1000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 1000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh paris     ${CURR_DATASET} 1000000 20480  256 2000  1 1 1 0;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 1000000  20480  256 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 1000000  120000  256 1 1 0 1 1 1 2 1;

#25GB Dataset
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 25000000 61440  256 100000 12 24 12 120000;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 25000000  61440  256 100000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/index_experiment.sh dstree    ${CURR_DATASET} 25000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 25000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh vaplus    ${CURR_DATASET} 25000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 25000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh paris     ${CURR_DATASET} 25000000 20480  256 2000 1 1 1 0;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 25000000  20480  256 2000 1 0 1 1 1 1 1;


#50GB Dataset
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 50000000 61440  256 100000 12 24 12 120000;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 50000000  61440  256 100000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/index_experiment.sh dstree    ${CURR_DATASET} 50000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 50000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh vaplus    ${CURR_DATASET} 50000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 50000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh paris     ${CURR_DATASET} 50000000 20480  256 2000 1 1 1 0;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 50000000  20480  256 2000 1 0 1 1 1 1 1;


#100GB Dataset 
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 100000000 61440  256 100000 12 24 12 120000;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/index_experiment.sh dstree    ${CURR_DATASET} 100000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh vaplus    ${CURR_DATASET} 100000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh paris     ${CURR_DATASET} 100000000 20480  256 2000 1 1 1 0;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  20480  256 2000 1 0 1 1 1 1 1;


#250GB Dataset
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 250000000 61440  256 100000 12 24 12 120000;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 250000000  61440  256 100000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/index_experiment.sh dstree    ${CURR_DATASET} 250000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 250000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh vaplus    ${CURR_DATASET} 250000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 250000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh paris     ${CURR_DATASET} 250000000 20480  256 2000 1 1 1 0;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 250000000  20480  256 2000 1 0 1 1 1 1 1;

##############################################
#Figure 6.b : Extrapolate from logs of Fig 6.a 
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

