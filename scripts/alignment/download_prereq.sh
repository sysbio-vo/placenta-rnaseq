#!/bin/sh

mkdir hg38

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.primary_assembly.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz


wget ftp://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz -O "hg38gz/hg38.fa.gz"
wget ftp://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf.gz -O "hg38gz/hg38.gtf.gz"
wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz -O "hg38gz/hg38.gene_info.gz"

cd hg38gz
gunzip hg38.fa.gz
gunzip hg38.gtf.gz
gunzip hg38.gene_info.gz

