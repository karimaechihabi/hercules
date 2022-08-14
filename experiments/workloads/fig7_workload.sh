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


#########################################
#Figure 7: Scalability with Dataset Size
#########################################

#1GB Dataset (sanity check)
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 1000000 61440  256 1000 12 24 12 120000;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 1000000  61440  256 1000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/index_experiment.sh dstree    ${CURR_DATASET} 1000000 61440  256 1000 1 1 1 0;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 1000000  61440  256 1000 1 0 1 1 1 1 1;
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
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 25000000  120000  256 1 1 0 1 1 1 2 1;


#50GB Dataset
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 50000000 61440  256 100000 12 24 12 120000;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 50000000  61440  256 100000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/index_experiment.sh dstree    ${CURR_DATASET} 50000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 50000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh vaplus    ${CURR_DATASET} 50000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 50000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh paris     ${CURR_DATASET} 50000000 20480  256 2000 1 1 1 0;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 50000000  20480  256 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 50000000  120000  256 1 1 0 1 1 1 2 1;


#100GB Dataset 
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 100000000 61440  256 100000 12 24 12 120000;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/index_experiment.sh dstree    ${CURR_DATASET} 100000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh vaplus    ${CURR_DATASET} 100000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh paris     ${CURR_DATASET} 100000000 20480  256 2000 1 1 1 0;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  20480  256 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 100000000  120000  256 1 1 0 1 1 1 2 1;


#250GB Dataset
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 250000000 61440  256 100000 12 24 12 120000;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 250000000  61440  256 100000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/index_experiment.sh dstree    ${CURR_DATASET} 250000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 250000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh vaplus    ${CURR_DATASET} 250000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 250000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh paris     ${CURR_DATASET} 250000000 20480  256 2000 1 1 1 0;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 250000000  20480  256 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 250000000  120000  256 1 1 0 1 1 1 2 1;


#1TB Dataset 
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 1000000000 61440  256 100000 12 24 12 120000;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 1000000000  61440  256 100000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/index_experiment.sh dstree    ${CURR_DATASET} 1000000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 1000000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh vaplus    ${CURR_DATASET} 1000000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 1000000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh paris     ${CURR_DATASET} 1000000000 20480  256 2000 1 1 1 0;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 1000000000  20480  256 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 1000000000  120000  256 1 1 0 1 1 1 2 1;


#1.5TB Dataset 
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 1500000000 61440  256 100000 12 24 12 120000;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 1500000000  61440  256 100000 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/index_experiment.sh dstree    ${CURR_DATASET} 1500000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 1500000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh vaplus    ${CURR_DATASET} 1500000000 61440  256 100000 1 1 1 0;
#$SCRIPTS/query_experiment.sh vaplus    ${CURR_QUERIES} 100 ${CURR_DATASET} 1500000000  61440  256 100000 1 0 1 1 1 1 1;
#$SCRIPTS/index_experiment.sh paris     ${CURR_DATASET} 1500000000 20480  256 2000 1 1 1 0;
#$SCRIPTS/query_experiment.sh paris     ${CURR_QUERIES} 100 ${CURR_DATASET} 1500000000  20480  256 2000 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh pscan     ${CURR_QUERIES} 100 ${CURR_DATASET} 1500000000  120000  256 1 1 0 1 1 1 2 1;

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

