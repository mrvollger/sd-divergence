import pandas as pd

df = pd.read_csv(snakemake.input.bed, sep="\t", low_memory=False)
print(f"All in the vcf:\t{df.shape[0]/1e6}M")
df = df[((df["TYPE"] == "SNV") | (df["TYPE"] == "SNP")) & (df["#CHROM"] != "chrY")]
print(f"SNP and not chrY:\t{df.shape[0]/1e6}M")

df = df.astype(str)
df["GT"] = df["GT"].str.strip()
df["HAP"] = df["HAP"].str.strip()
df = (
    df.apply(lambda col: col.str.split(";").explode())
    .reset_index()
    .reindex(df.columns, axis=1)
)

#     638  .|.
#     204  .|0
#      46  0|.
# 2383400  0|1
#   61908  .|1
#   68572  1|.
# 2393786  1|0
# 2325920  1|1
# fix the genotypes for expanded
df.loc[(df["HAP"] == "h1") & (df["GT"] == "1|1"), "GT"] = "1|."
df.loc[(df["HAP"] == "h2") & (df["GT"] == "1|1"), "GT"] = ".|1"

is_h1 = (df["HAP"] == "h1") & (
    (df["GT"] == "1|0") | (df["GT"] == "1|1") | (df["GT"] == "1|.")
)
is_h2 = (df["HAP"] == "h2") & (
    (df["GT"] == "0|1") | (df["GT"] == "1|1") | (df["GT"] == ".|1")
)

df = df[is_h1 | is_h2]

df["SAMPLE"] = snakemake.wildcards.sm
print(f"Final output exploded:\t{df.shape[0]/1e6}M")
df.to_csv(snakemake.output.snv, sep="\t", compression="gzip", index=False)
