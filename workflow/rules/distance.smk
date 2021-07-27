rule per_file_distance_snv:
    input:
        snv=rules.filter_snv_by_syntenic.output.snv,
        dist=rules.clean_distance_files.output,
    output:
        temp("temp/distance/{dist}/dist_{sm}_{h}.bed"),
    log:
        "logs/distance/{dist}/dist_{sm}_{h}.log",
    conda:
        "../envs/env.yml"
    threads: 1
    shell:
        """
        gunzip -c {input.snv} \
            | cut -f 1-3 \
            | bedtools closest -d -t first \
                -a - \
                -b {input.dist} \
            | cut -f 1-3,7 \
            | bedtools sort -i - \
            > {output}
        """


rule distance_snv:
    input:
        snv=rules.annotate_snv.output,
        distances=expand(
            rules.per_file_distance_snv.output,
            dist=config["distance_files"].keys(),
            allow_missing=True,
        ),
    output:
        temp("temp/snv/{sm}_{h}/dist_{sm}_{h}.bed"),
    log:
        "logs/snv/{sm}_{h}/dist_{sm}_{h}.log",
    conda:
        "../envs/env.yml"
    threads: 1
    params:
        dist="\t".join([f"dist_{key}" for key in config["distance_files"].keys()]),
        columns=",".join(
            [str(4 * (i + 1)) for i in range(len(config["distance_files"].keys()))]
        ),
    shell:
        """
        HEADER=$(cat {input.snv} | head -n 1 || :)
        HEADER="${{HEADER}}\t{params.dist}"
        echo $HEADER

        paste \
            {input.snv} \
            <(paste {input.distances} | cut -f {params.columns} ) \
            | grep -v "^#" | sed "1s/^/${{HEADER}}\\n/" \
            > {output}
        """


rule all_snv:
    input:
        snv=expand(rules.distance_snv.output, sm=tbl.index, h=[1, 2]),
    output:
        "results/all_snv_exploded.bed.gz",
    log:
        "logs/all_snv.log",
    conda:
        "../envs/env.yml"
    threads: 8
    shell:
        """
        HEADER=$(head -n 1 {input.snv[1]} || :)
        echo $HEADER

        sort -m -k 1,1 -k2,2n {input.snv} \
            | grep -v "^#" \
            | sed "1s/^/${{HEADER}}\\n/" \
            | pigz -p {threads} \
            > {output}
        """


rule small_snv:
    input:
        rules.all_snv.output,
    output:
        "results/small_snv_exploded.bed.gz",
    log:
        "logs/small_snv.log",
    threads: 1
    run:
        df = pd.read_csv(input[0], sep="\t")
        dist_col = [
            col for col in df if col.startswith("dist_") or col.startswith("anno_")
        ]
        names = ["#CHROM", "POS", "END", "REF", "ALT", "HAP", "SAMPLE"] + dist_col
        df[names].to_csv(output[0], sep="\t", index=False, compression="gzip")
