import pandas as pd

[print(f) for f in snakemake.input]

counts = [pd.read_table(f, index_col=0, usecols=[0, 2], 
          header=None, skiprows=4) 
          for f in snakemake.input]

for t, sample in zip(counts, snakemake.params.samples):
    t.columns = [sample]

matrix = pd.concat(counts, axis=1)
matrix.index.name = "gene"
matrix.to_csv(snakemake.output[0], sep="\t")
