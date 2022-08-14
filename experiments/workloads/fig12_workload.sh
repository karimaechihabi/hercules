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
#Figure 12: Ablation Study
####################################################
CURR_DATASET=data_size1B_deep1b_len96_znorm
#$SCRIPTS/index_experiment.sh hercules  ${CURR_DATASET} 266666667 61440  96 266667 12 24 12 120000;
CURR_QUERIES=queries_hard1p_deep1b_len96_znorm
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 0.00 0.00 24 1;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 1.00 1.00 24 1;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 0.25 1.00 24 1;
CURR_QUERIES=queries_hard5p_deep1b_len96_znorm
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 0.00 0.00 24 1;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 1.00 1.00 24 1;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 0.25 1.00 24 1;
CURR_QUERIES=queries_orig100_deep1b_len96_znorm
#$SCRIPTS/query_experiment.sh dstree    ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 1 1 1 1;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 0.25 0.9 24 1;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 0.00 0.00 24 1;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 1.00 1.00 24 1;
#$SCRIPTS/query_experiment.sh hercules  ${CURR_QUERIES} 100 ${CURR_DATASET} 266666667  61440  96 266667 1 0 1 0.25 1.00 24 1;


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

