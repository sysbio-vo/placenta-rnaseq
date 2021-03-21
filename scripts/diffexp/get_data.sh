#!/bin/sh

gdown https://drive.google.com/uc?id=1sUWoYRhl27Fny_tfD9RUZ6FOv8SQuvPz

mkdir raw_counts
mv china_2019_study.csv raw_counts/

cp ../aligment/counts/transcripts_count_matrix.csv raw_counts/ || (gdown https://drive.google.com/uc?id=1z9ADd-5LFbW_6VCVr95EBZ0dxEqT1kly && mv transcripts_count_matrix.csv raw_counts/)

python merge_counts.py
