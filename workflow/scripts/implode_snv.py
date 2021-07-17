import pandas as pd

li = []
for filename in snakemake.input.beds:
    tdf = pd.read_csv(filename, index_col=None, sep="\t", low_memory=False)
    li.append(tdf)

df = pd.concat(li, axis=0, ignore_index=True)
df = df.astype(str)

id_cols = ["ID"]
if "SAMPLE" in df.columns:
    id_cols = ["ID", "SAMPLE"]

implode = (
    df.sort_values(by=["#CHROM", "POS", "ID", "HAP"])
    .groupby(id_cols)
    .agg(lambda x: ";".join(x.unique()))
    .reset_index()
)
implode.loc[(implode["GT"] == "1|.;.|1"), "GT"] = "1|1"

implode = implode.loc[:, df.columns]

implode.to_csv(snakemake.output.snv, sep="\t", compression="gzip", index=False)
