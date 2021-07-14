import pandas as pd

df = pd.read_csv(snakemake.input.bed, sep="\t")
df = df.astype(str)
df = (
    df.apply(lambda col: col.str.split(";").explode())
    .reset_index()
    .reindex(df.columns, axis=1)
)
df["SAMPLE"] = snakemake.wildcards.sm
df.to_csv(snakemake.output.snv, sep="\t", compression="gzip", index=False)

