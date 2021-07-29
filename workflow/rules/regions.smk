#
# per bed files
#
rule syntenic_and_callable:
    input:
        callable=lambda w: tbl.loc[w.sm][f"h{w.h}_callable"],
        aln=lambda w: tbl.loc[w.sm][f"h{w.h}_aln"],
    output:
        "results/syntenic_and_callable/{sm}_{h}.bed.gz",
    log:
        "logs/syntenic.{sm}_{h}.log",
    conda:
        "../envs/env.yml"
    params:
        min_syntenic_size=config["min_syntenic_size"],
    threads: 1
    shell:
        """
        bedtools intersect -a {input.callable} \
            -b <( bedtools sort -i {input.aln} | awk '$3 - $2 >= {params.min_syntenic_size}' ) \
            | bedtools sort -i - \
            | bedtools merge -i - \
            | sed 's/$/\\t{wildcards.sm}_{wildcards.h}/g' \
            | gzip -c > {output}
        """


rule clean_annotation_files:
    input:
        annotation_file=lambda w: config["annotation_files"][w.anno],
    output:
        temp("temp/anno/{anno}.bed.gz"),
    log:
        "logs/clean_annotation_files.{anno}.log",
    conda:
        "../envs/env.yml"
    threads: 1
    shell:
        """
        bedtools sort -i {input.annotation_file} \
            | bedtools merge -i - \
            | gzip -c > {output}
        """


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
# combinations of annotation files
#
rule combinations_of_annotation_files:
    input:
        bed=rules.syntenic_and_callable.output,
        anno1="temp/anno/{anno1}.bed.gz",
        anno2="temp/anno/{anno2}.bed.gz",
    output:
        "results/annotation/combos/{sm}_{h}_intersection_{anno1}_and_{anno2}.bed.gz",
    log:
        "logs/combo_{sm}_{h}_{anno1}_{anno2}.log",
    conda:
        "../envs/env.yml"
    threads: 1
    shell:
        """
        bedtools intersect \
            -sorted \
            -a {input.bed} \
            -b {input.anno1} \
        | bedtools intersect \
            -sorted \
            -a - \
            -b {input.anno2} \
        | bedtools sort -i - \
        | bedtools merge -i - \
        | gzip -c \
        > {output}
        """


rule size_of_combinations_of_annotation_files:
    input:
        rules.combinations_of_annotation_files.output,
    output:
        tbl=temp("temp/regions/{sm}_{h}_intersection_{anno1}_and_{anno2}.tbl"),
    log:
        "logs/combo_{sm}_{h}_{anno1}_{anno2}.log",
    conda:
        "../envs/env.yml"
    threads: 1
    shell:
        """
        printf "{wildcards.sm}_{wildcards.h}\t" > {output}
        if [ {wildcards.anno1} == {wildcards.anno2} ]; then
            printf "{wildcards.anno1}_size\t" >> {output}
        else
            printf "{wildcards.anno1}_{wildcards.anno2}_size\t" >> {output}
        fi 
        awk '{{sum += $3 - $2}}END{{print sum}}' <(gunzip -c {input}) >> {output}
        """


def allowed_anno_combinations(wc):
    file_fmt = rules.size_of_combinations_of_annotation_files.output.tbl
    anno = config["annotation_files"].keys()
    out = []
    for sm in tbl.index:
        for h in [1, 2]:
            for a1 in anno:
                for a2 in anno:
                    f = (file_fmt).format(h=h, sm=sm, anno1=a1, anno2=a2)
                    if a1 == a2 or a1 in ["SD", "sd"]:
                        out.append(f)
                    else:
                        continue
    return out


rule make_combos:
    input:
        combos=allowed_anno_combinations,
    output:
        "results/annotation/annotation_sizes.tbl",
    log:
        "logs/combo_sizes.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        printf "hap\tanno\tsize\n" > {output}
        cat {input.combos} | sort >> {output}
        """


rule annotation_sizes_wide:
    input:
        rules.make_combos.output,
    output:
        "results/annotation/annotation_sizes_wide.tbl",
    log:
        "logs/combo_wide_sizes.log",
    threads: 1
    run:
        df = pd.read_csv(str(input), sep="\t")
        out = df.pivot(index="hap", columns="anno", values="size")
        out.to_csv(str(output), sep="\t")
