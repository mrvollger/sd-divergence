
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
        bed="results/windows.bed.gz",
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
        u_condition = (
            (df["anno_SD"] < 0.2) & (df["anno_Sat"] < 0.2) & (df["region"] == "Other")
        )
        df.loc[u_condition, "region"] = "Unique"
        df.to_csv(output.bed, sep="\t", index=False)
