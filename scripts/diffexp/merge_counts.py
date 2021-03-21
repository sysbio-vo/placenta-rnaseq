import glob
import pandas as pd


csv_files = glob.glob("./counts/*.csv")
tsv_files = glob.glob("./counts/*.tsv")

files_to_merge = []

for csv_file in csv_files:
    df = pd.read_csv(csv_file)
    df = df.set_index("gene", drop=True)
    files_to_merge.append(df)


for tsv_file in tsv_files:
    df = pd.read_table(tsv_file, sep='\t')
    df["gene"] = df["gene"].str.split(".").str[0]
    df = df.set_index("gene", drop=True)
    files_to_merge.append(df)


start_file = files_to_merge[0]

for f_merge in files_to_merge[1:]:
    start_file = pd.merge(start_file, f_merge, how='outer', left_index=True, right_index=True)

start_file = start_file.fillna(0)
start_file = start_file[~start_file.index.duplicated(keep='first')]
start_file.to_csv("merged_counts.csv", sep=',')
