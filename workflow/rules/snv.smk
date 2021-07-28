rule explode_snv:
    input:
        bed=lambda w: tbl.loc[w.sm][f"snv"],
    output:
        snv=temp("temp/syntenic_and_callable/snv_explode_{sm}.bed.gz"),
    log:
        "logs/explode_snv.{sm}.log",
    conda:
        "../envs/env.yml"
    threads: 1
    script:
        "../scripts/explode_snv.py"


rule filter_snv_by_syntenic:
    input:
        snv=rules.explode_snv.output.snv,
        callable=rules.syntenic_and_callable.output,
    output:
        snv=temp("temp/syntenic_and_callable/snv_{sm}_{h}.bed.gz"),
    log:
        "logs/snv_by_syntenic.{sm}_{h}.log",
    conda:
        "../envs/env.yml"
    threads: 1
    shell:
        """
        gunzip -c {input.snv} \
            | awk -F$'\\t' -v pat="h{wildcards.h}|HAP" '$15 ~ pat' \
            | bedtools intersect -u \
                -a - \
                -b <(gunzip -c {input.callable}) \
                -header \
            | bedtools sort -header -i - \
        | gzip -c > {output}
        """


#
# group results
#
rule merge_filtered_snv_by_syntenic:
    input:
        beds=expand(rules.filter_snv_by_syntenic.output, h=[1, 2], allow_missing=True),
    output:
        snv="results/syntenic_and_callable/snv_{sm}.bed.gz",
    log:
        "logs/merge_filtered_snv_by_syntenic.{sm}.log",
    conda:
        "../envs/env.yml"
    threads: 1
    script:
        "../scripts/implode_snv.py"


rule annotate_snv:
    input:
        snv=rules.filter_snv_by_syntenic.output.snv,
        annotation_files=expand(
            rules.clean_annotation_files.output, anno=config["annotation_files"].keys()
        ),
    output:
        temp("temp/snv/anno_{sm}_{h}.bed"),
    log:
        "logs/{sm}_{h}_snv_anno.log",
    conda:
        "../envs/env.yml"
    params:
        annotation_names="\t".join(
            [f"anno_{key}" for key in config["annotation_files"].keys()]
        ),
    threads: 1
    shell:
        """
        HEADER=$(gunzip -c {input.snv} | head -n 1 || :)
        HEADER="${{HEADER}}\t{params.annotation_names}"
        echo $HEADER

        gunzip -c {input.snv} \
            | bedtools annotate -i - -counts \
                -files {input.annotation_files} \
            | bedtools sort -i - \
            | sed "1s/^/${{HEADER}}\\n/" \
            > {output}
        """


#
# adding in the SNVs
#
rule add_snv_to_windows:
    input:
        windows=rules.distance_windows.output,
        snv=rules.filter_snv_by_syntenic.output.snv,
    output:
        temp("temp/add_snv/{sm}_{h}_snvs_haplotype_coverage.bed"),
    log:
        "logs/snv_count_windows_{sm}_{h}.log",
    conda:
        "../envs/env.yml"
    threads: 1
    shell:
        """
        HEADER=$(gunzip -c {input.windows} | head -n 1 || :)
        HEADER="${{HEADER}}\tnum_snv\thap"
        echo $HEADER

        gunzip -c {input.windows} \
            | grep "{wildcards.sm}_{wildcards.h}" \
            | bedtools coverage \
                -a - \
                -b {input.snv} \
                -counts -sorted \
            | sed "s/$/\\t{wildcards.sm}_{wildcards.h}/" \
            | sed "1s/^/${{HEADER}}\\n/" \
            | csvtk cut -C "$" -tT -f -haps \
            > {output}
        """


#
# make a long filtered form of the df
#
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
        mem=8
    threads: 8
    shell:
        """
        HEADER=$(head -n 1 {input.snv[1]} || :)
        echo $HEADER

        sort -m -k 1,1 -k2,2n \
            -S {recources.mem} --parallel {threads} \
            {input.snv} \
            | grep -v "^#" \
            | sed "1s/^/${{HEADER}}\\n/" \
            | pigz -p {threads} \
            > {output}
        """
