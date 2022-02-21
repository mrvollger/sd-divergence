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
    resources:
        load=100,
    script:
        "../scripts/explode_snv.py"


rule filter_snv_by_syntenic:
    input:
        snv=rules.explode_snv.output.snv,
        callable=rules.syntenic_and_callable.output,
        exclude=config["exclude"],
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
            | awk -F$'\\t' -v pat="h{wildcards.h}|HAP" '$9 ~ pat' \
            | bedtools intersect -u \
                -a - \
                -b <(gunzip -c {input.callable}) \
                -header \
            | bedtools intersect -header -v -a - -b {input.exclude} \
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


def get_annotation_files(wc):
    rtn = expand(
        rules.clean_annotation_files.output, anno=config["annotation_files"].keys()
    )
    anno = f"gene-conversion-{wc.sm}_{wc.h}"
    rtn.append(expand(rules.clean_annotation_files.output, anno=anno)[0])
    return rtn


rule annotate_snv:
    input:
        snv=rules.filter_snv_by_syntenic.output.snv,
        annotation_files=get_annotation_files,
    output:
        temp("temp/snv/anno_{sm}_{h}.bed"),
    log:
        "logs/{sm}_{h}_snv_anno.log",
    conda:
        "../envs/env.yml"
    params:
        annotation_names="\t".join(
            [
                f"anno_{key}"
                for key in list(config["annotation_files"].keys()) + ["IGC"]
            ]
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
