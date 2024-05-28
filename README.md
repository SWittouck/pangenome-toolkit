# SCARAP paper

This repository contains code to:

* Analyze benchmarking results of SCARAP, a toolkit for comparative genomics of prokaryotes ([benchmark-scarap repo](https://github.com/swittouck/benchmark-scarap)) 
* Explore a large-scale *Lactobacillales* pangenome created by SCARAP ([LEGEN repo](https://github.com/swittouck/legen))

A manuscript describing these two analyses has been submitted to an open access journal. 

## How to run 

Step 1: clone the repository: 

    git clone https://github.com/SWittouck/pangenome-toolkit.git

Step 2: download the main output files from the benchmark-scarap and LEGEN repositories: 

    cd pangenome-toolkit
    mkdir data
    cd data
    wget https://github.com/SWittouck/benchmark-scarap/releases/download/v3/benchmark_scarap_v3.tar.gz
    wget https://github.com/SWittouck/legen/releases/download/v4/legen_v4.tar.gz 
    tar xzf benchmark_scarap_v3.tar.gz --one-top-level
    tar xzf legen_v4.tar.gz --one-top-level
    cd ..

Step 3: install the most important dependencies: 

* R v4.2.3
* tidyverse v1.3.1

(Some scripts need extra R packages - they are listed at the top of the script.)

Step 4: create a folder for the results: 

    mkdir results 

Step 5: run your script of interest, e.g. 

    Rscript ./src/01_demonstrate_pan.R