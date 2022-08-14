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


CURR_QUERIES=queries_size100K_seed14784_len256_znorm
CURR_DATASET=data_size2B_seed1184_len256_znorm


#######################################################
#Figure 10a-10b: Extract numbers from logs of Figure 9a
#######################################################

#######################################################
#Figure 10c-10d: Extract numbers from logs of Figure 9c
#######################################################

#######################################################
#Figure 10e-10bf Extract numbers from logs of Figure 9e
#######################################################


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

