#!/bin/bash

## How to install ascp, in a gist.

## Check for latest link: http://downloads.asperasoft.com/en/downloads/8?list
wget -qO- https://download.asperasoft.com/download/sw/connect/3.9.8/ibm-aspera-connect-3.9.8.176272-linux-g2.12-64.tar.gz | tar xvz

## run it
chmod +x ibm-aspera-connect-3.9.8.176272-linux-g2.12-64.sh
./ibm-aspera-connect-3.9.8.176272-linux-g2.12-64.sh

## add it to the path now and in the future
export PATH=$PATH:~/.aspera/connect/bin/
echo 'export PATH=$PATH:~/.aspera/connect/bin/' >> ~/.bash_profile



## Snakemake
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p ./miniconda

eval "$(./miniconda/bin/conda shell.bash hook)"

conda env create --name snakemake-rnaseq --file environment.yaml
conda activate snakemake-rnaseq




qsub pbs_index_hg_38_star
