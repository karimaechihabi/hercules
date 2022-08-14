mkdir -p $EXPERIMENTS/logs    
mkdir -p $EXPERIMENTS/results
mkdir -p $EXPERIMENTS/results_post_processed/logs_csv/
mkdir -p $EXPERIMENTS/plots

CURR_DIR=$PWD

##################################
#Run experiments and produce plots
##################################
cd $EXPERIMENTS/
workloads/fig6_workload.sh
workloads/fig7_workload.sh
workloads/fig8_workload.sh
workloads/fig9_workload.sh
workloads/fig10_workload.sh
workloads/fig11_workload.sh
workloads/fig12_workload.sh

###################################
#Regenerate paper from source files
##################################
#cd $PROJECTS_ROOT/paper/
#pdflatex paper.tex
#pdflatex repro_paper.tex

cd $CURR_DIR;

