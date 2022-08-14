# 1. Paper Information

## Title: Karima Echihabi, Panagiota Fatourou, Kostas Zoumpatianos, Themis Palpanas, and Houda Benbrahim. Hercules Against Data Series Similarity Search. PVLDB, 15(10), 2022.

## Abstract: 

In this paper, we propose Hercules, a parallel tree-based technique for exact similarity search on massive disk-based data series collections. We present novel index construction and query answering algorithms that leverage different summarization techniques, carefully schedule costly operations, optimize memory and disk accesses, and exploit the multi-threading and SIMD capabilities of modern hardware to perform CPU-intensive calculations. We demonstrate the superiority and robustness of Hercules with an extensive experimental evaluation against the state-of-the-art techniques, using a variety of synthetic and real datasets, and query workloads of varying difficulty. The results show that Hercules performs up to one order of magnitude faster than the best competitor (which is not always the same). Moreover, Hercules is the only index that outperforms the optimized sequential scan on all scenarios, including the hard query workloads on disk-based datasets. 

## Paper Link: https://github.com/karimaechihabi/hercules/blob/main/paper/p1064-echihabi.pdf

# 2. Reproducibilty

## 2.1. Archive
This archive contains detailed information required to reproduce the experimental results of the above paper.

The archive contains the following 4 subdirectories:

### data
The data subdirectory contains the datasets and queries used in the paper. 

All queries are contained in the archive.

A script to download the real datasets can be found at hercules/experiments/data/real/download_real_data.sh. The script also provides the links to the google drive archives containing the datasets. 

A script to generate the synthetic datasets and queries can be found at hercules/experiments/data/synthetic/generate_synthetic_data.sh.

### tools
The tools subdirectory contains the tools used by the archive to generate or download data.

### experiments
Upon publication of the paper, the experiments subdirectory will contain the following:

   <u>bin</u>: executables required to run the experiments (can be used as a back-box). \
   <u>config</u>: configurations to add to the .bashrc file. \
   <u>logs</u>: stores the logs generated once experiments finish. \
   <u>scripts</u>: scripts used to automate experiments. \
   <u>workloads</u>: contains scripts to schedule experiments.

## 2.2. Software Requirements
Linux Ubuntu 16.04.2 \
GCC 6.2.0 \
Python 2.7.13 \
R 3.4.0 \
jemalloc (http://jemalloc.net/). Make sure to add the location of jemalloc to the path. \

## 2.3. Hardware Requirements
No particular requirement other than a machine equipped with an HDD having at least 75GB of RAM. The experiments in the paper were run on a server with two Intel Xeon E5-2650 v4 2.2GHz CPUs, 75GB of RAM, and 10.8TB (6 x 1.8TB) 10K RPM SAS hard drives in RAID0 with a throughput of 1290 MB/sec. Although using a different machine might lead to running times different from the paper, the overall trends should stay the same.

If the RAM is over 75 GB, restrict its size using the grub as follows:

sudo gedit  /etc/default/grub  \
add or update the GRUB_CMDLINE_LINUX_DEFAULT variable as GRUB_CMDLINE_LINUX_DEFAULT="quiet splash mem=75G" \
sudo update-grub \
sudo reboot 

## 2.4. Reproducibility Steps

### Configure the Setup

Run: chmod u+x hercules/experiments/tools/misc/give_x_rights.sh

Then obtain execute rights on the scripts and binaries by running: hercules/experiments/tools/misc/give_x_rights.sh

Configure your environment using the project.config file in hercules/experiments/config/ and modify the DATASETS and PROJECT_ROOT variables to point to the right path in your system.

### Generate Datasets
Generate the synthetic datasets and download the real datasets. Put all datasets and queries in a directory referenced as $DATASETS.

### Run Experiments

Modify the workloads in hercules/experiments/workloads/ to specify which plots to reproduce. If only a subset of a given figure is to be reproduced, modify the workload corresponding to this figure by uncommenting the commands for the specific subfigures. Some experiments can take several days to complete, so we put an estimate for the completion time of each figure and subfigure inside each workload file. Note that these estimates are based on the running times on our server and may vary given your hardware.

Each workload file is standalone, and currently all experiments are commented out. Some workload files run the same experiment, please only uncomment this experiment once. For instance, fig6_workload.sh and fig9_workload.sh both build indexes, so the commands for building these indexes should be uncommented only for one of them (unless the indexes were deleted for space considerations).  

Currently all commands are commented out. 

Once the workloads are modified, run hercules/experiments/repro_paper.sh

 

