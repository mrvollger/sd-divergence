import pandas as pd

df = pd.read_csv(snakemake.input.bed, sep="\t", low_memory=False)
print(f"All things in the vcf: {df.shape[0]/1e6}M")
df = df[((df["TYPE"] == "SNV") | (df["TYPE"] == "SNP")) & (df["#CHROM"] != "chrY")]
print(f"SNP and not chrY: {df.shape[0]/1e6}M")

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
print(f"SNP and not chrY exploded: {df.shape[0]/1e6}M")
df.to_csv(snakemake.output.snv, sep="\t", compression="gzip", index=False)
