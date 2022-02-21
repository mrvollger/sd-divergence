import pandas as pd

df = pd.read_csv(snakemake.input.bed, sep="\t", low_memory=False)
print(df.shape)
df = df[((df["TYPE"] == "SNV") | (df["TYPE"] == "SNP")) & (df["#CHROM"] != "chrY")]
print(df.shape)

df = df.astype(str)
df = (
    df.apply(lambda col: col.str.split(";").explode())
    .reset_index()
    .reindex(df.columns, axis=1)
)
# fix the genotypes for expanded
df.loc[(df["HAP"] == "h1") & (df["GT"] == "1|1"), "GT"] = "1|."
df.loc[(df["HAP"] == "h2") & (df["GT"] == "1|1"), "GT"] = ".|1"

df["SAMPLE"] = snakemake.wildcards.sm
print(df.shape)
df.to_csv(snakemake.output.snv, sep="\t", compression="gzip", index=False)
