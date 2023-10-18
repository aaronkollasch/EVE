#!/bin/bash
# run on c5ad.24xlarge or similar
# preferred OS: Ubuntu

sudo apt update && sudo apt install -y awscli
sudo mkfs.ext4 /dev/nvme1n1
sudo mount /dev/nvme1n1 /tmp
sudo chmod a+w /tmp

curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash "Mambaforge-$(uname)-$(uname -m).sh" -b
~/mambaforge/bin/conda init
source ~/.bashrc

mamba create -y -n py310 python=3.10 tqdm pandas numpy numba
conda activate py310

git clone https://github.com/aaronkollasch/EVE
cd EVE/scripts && aws s3 cp s3://markslab-us-east-2/colabfold/output/popeve/missing_genes_with_priority.tsv .
