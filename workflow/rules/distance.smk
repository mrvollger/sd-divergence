#
# distance cleaning
#
rule clean_distance_files:
    input:
        annotation_file=lambda w: config["distance_files"][w.dist],
    output:
        temp("temp/distance/{dist}.bed.gz"),
    log:
        "logs/clean_distance_files.{dist}.log",
    conda:
        "../envs/env.yml"
    threads: 1
    shell:
        """
        bedtools sort -i {input.annotation_file} \
            | gzip -c > {output}
        """


#
# Distance for windows
#
rule per_file_distance_window:
    input:
        windows=rules.window_regions.output,
        dist=rules.clean_distance_files.output,
    output:
        temp("temp/distance/{dist}/window_dist.bed"),
    log:
        "logs/distance/{dist}/window_dist.log",
    conda:
        "../envs/env.yml"
    threads: 1
    shell:
        """
        NUM_COL=$(gunzip -c {input.dist} | head -n 1 | awk '{{print NF}}' || :)
        LAST_COL=$((NUM_COL+4))
        echo $NUM_COL, $LAST_COL 

        gunzip -c {input.windows} \
            | grep -v "^#" \
            | cut -f 1-3 \
            | bedtools closest -D b -t first \
                -a - \
                -b {input.dist} \
            | cut -f 1,2,3,$LAST_COL \
            > {output}
        """


rule distance_windows:
    input:
        windows=rules.window_regions.output,
        distances=expand(
            rules.per_file_distance_window.output,
            dist=config["distance_files"].keys(),
            allow_missing=True,
        ),
    output:
        temp("temp/distance_windows.bed.gz"),
    log:
        "logs/windows.log",
    conda:
        "../envs/env.yml"
    threads: 8
    params:
        dist="\t".join([f"dist_{key}" for key in config["distance_files"].keys()]),
        columns=",".join(
            [str(4 * (i + 1)) for i in range(len(config["distance_files"].keys()))]
        ),
    shell:
        """
        HEADER=$(gunzip -c {input.windows} | head -n 1 || :)
        HEADER="${{HEADER}}\t{params.dist}"
        echo $HEADER

        paste \
            <(gunzip -c {input.windows} | grep -v "^#") \
            <(paste {input.distances} | cut -f {params.columns} ) \
            | sed "1s/^/${{HEADER}}\\n/" \
            | pigz -p {threads} \
            > {output}
        """


#
# Add SNV counts to the distance windows
#
rule add_snv_to_windows:
    input:
        windows=rules.distance_windows.output,
        #snv=rules.filter_snv_by_syntenic.output.snv,
        snv="temp/snv/{sm}_{h}/dist_{sm}_{h}.bed",
    output:
        "results/windowed_snv/{sm}_{h}_snvs_haplotype_coverage.bed.gz",
    log:
        "logs/snv_count_windows_{sm}_{h}.log",
    conda:
        "../envs/env.yml"
    threads: 4
    shell:
        """
        HEADER=$(gunzip -c {input.windows} | head -n 1 || :)
        HEADER="${{HEADER}}\tnum_snv\thap"
        echo $HEADER

        gunzip -c {input.windows} \
            | grep "{wildcards.sm}_{wildcards.h}" \
            | bedtools coverage \
                -a - \
                -b <(csvtk filter -C "$" -tT -f "anno_TRF<1" {input.snv}) \
                -counts -sorted \
            | sed "s/$/\\t{wildcards.sm}_{wildcards.h}/" \
            | sed "1s/^/${{HEADER}}\\n/" \
            | csvtk cut -C "$" -tT -f -haps \
            | pigz -p {threads} \
            > {output}
        """


rule long_and_filtered_windows:
    input:
        snv=expand(rules.add_snv_to_windows.output, sm=tbl.index, h=[1, 2]),
    output:
        snv="results/long_windows_with_snv_dist_annotation.bed.gz",
    log:
        "logs/long_and_filter.log",
    conda:
        "../envs/env.yml"
    params:
        names=" ".join(expand("{sm}_{h}", sm=tbl.index, h=[1, 2])),
    resources:
        mem=8,
    threads: 40
    shell:
        """
        HEADER=$(zcat {input.snv[1]} | head -n 1 || :)
        echo $HEADER
        zcat {input.snv} \
            | sort -k 1,1 -k2,2n \
                --parallel={threads} -S 80G \
            | grep -v "^#" \
            | sed "1s/^/${{HEADER}}\\n/" \
            | pigz -p {threads} \
            > {output}
        """
        #sort -m -k 1,1 -k2,2n \
        #--batch-size=500 --parallel {threads} \
        # {input.snv} \


#
# DISTANCE for SNVs
#
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
        NUM_COL=$(gunzip -c {input.dist} | head -n 1 | awk '{{print NF}}' || :)
        LAST_COL=$((NUM_COL+4))
        echo $NUM_COL, $LAST_COL 

        gunzip -c {input.snv} \
            | cut -f 1-3 \
            | bedtools closest -D b -t first \
                -a - \
                -b {input.dist} \
            | cut -f 1-3,$LAST_COL \
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
        "results/snv/snv_{sm}_{h}.bed",
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
            <(grep -v "^#" {input.snv}) \
            <(paste {input.distances} | grep -v "^#" | cut -f {params.columns} ) \
            | grep -v "^#" | sed "1s/^/${{HEADER}}\\n/" \
            | csvtk filter -C "$" -tT -f "anno_TRF<1" \
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
    threads: 8
    shell:
        """
        zcat {input} | cut -f 1-3,15,18- | pigz -p {threads} > {output}
        """


"""
df = pd.read_csv(input[0], sep="\t")
dist_col = [
col for col in df if col.startswith("dist_") or col.startswith("anno_")
]
names = ["#CHROM", "POS", "END", "REF", "ALT", "HAP", "SAMPLE"] + dist_col
df[names].to_csv(output[0], sep="\t", index=False, compression="gzip")
"""
