rule make_windows:
    input:
        fai=f"{config['ref']}.fai",
    output:
        temp("temp/windows.bed.gz"),
    log:
        "logs/windows.log",
    conda:
        "../envs/env.yml"
    params:
        window_size=config["window_size"],
        step_size=config["step_size"],
    threads: 1
    shell:
        """
        bedtools makewindows -g {input.fai} -w {params.window_size} -s {params.step_size} \
            | bedtools sort -i - \
            | gzip -c > {output}
        """


#
# add columns
#
rule haplotype_coverage_over_windows:
    input:
        windows=rules.make_windows.output,
        beds=expand(rules.syntenic_and_callable.output, sm=tbl.index, h=[1, 2]),
    output:
        temp("temp/haplotype_coverage.bed.gz"),
    log:
        "logs/haplotype_coverage.log",
    conda:
        "../envs/env.yml"
    threads: 1
    params:
        names=" ".join(expand("{sm}_{h}", sm=tbl.index, h=[1, 2])),
        overlap=1 - config["window_size"],  # makes sure bedtools merge only merges the same windows
    shell:
        """
        NUM_COL=$(gunzip -c {input.windows} | head -n 1 | awk '{{print NF}}' || :)
        NAME_COL=$((NUM_COL+1))
        BOOL_COL=$((NUM_COL+2))
        echo $NAME_COL $BOOL_COL

        bedtools intersect -f 0.95 -C \
            -sorted -sortout \
            -a {input.windows} -b {input.beds} \
            -names {params.names} \
        | awk -v bool=$BOOL_COL '$bool != 0' \
        | bedtools merge -i - \
            -d {params.overlap} \
            -c $NAME_COL,$BOOL_COL -o distinct,sum \
        | sed '1s/^/#chr\\tstart\\tend\\thaps\\thap_count\\n/' \
        | gzip -c > {output}
        """


rule annotate_windows:
    input:
        windows=rules.haplotype_coverage_over_windows.output,
        annotation_files=expand(
            rules.clean_annotation_files.output, anno=config["annotation_files"].keys()
        ),
    output:
        temp("temp/annotated_windows.bed.gz"),
    log:
        "logs/annotate_windows.log",
    conda:
        "../envs/env.yml"
    params:
        annotation_names="\t".join(
            [f"anno_{key}" for key in config["annotation_files"].keys()]
        ),
    threads: 1
    shell:
        """
        HEADER=$(gunzip -c {input.windows} | head -n 1 || :)
        HEADER="${{HEADER}}\t{params.annotation_names}"
        echo $HEADER

        bedtools annotate -i {input.windows} \
                -files {input.annotation_files} \
            | bedtools sort -i - \
            | sed "1s/^/${{HEADER}}\\n/" \
            | gzip -c > {output}
        """


rule window_regions:
    input:
        bed=rules.annotate_windows.output,
    output:
        bed=temp("temp/no_dist_windows.bed.gz"),
    log:
        "logs/window_regions.log",
    threads: 1
    run:
        df = pd.read_csv(str(input.bed), sep="\t")
        anno_cols = [col for col in df.columns if col.startswith("anno_")]
        df["region"] = "Other"
        for anno in anno_cols:
            df.loc[
                (df[anno] >= 0.95) & (df["region"] == "Other"), "region"
            ] = anno.strip("anno_")
        df.to_csv(output.bed, sep="\t", index=False)


rule uncallable_windows:
    input:
        fai=fai,
        windows=rules.window_regions.output,
    output:
        "results/uncallable_windows.bed.gz",
    log:
        "logs/uncallable_windows.log",
    conda:
        "../envs/env.yml"
    params:
        annotation_names="\t".join(config["annotation_files"].keys()),
    threads: 1
    shell:
        """
        bedtools complement \
            -i <(bedtools merge -i {input.windows}) \
             -g <(cat {input.fai} | sort -k 1,1 -k2,2n ) \
            | gzip -c > {output}
        """
