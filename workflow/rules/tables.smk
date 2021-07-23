#
# summary table
#
rule make_summary_table:
    input:
        rules.all_snv.output,
        rules.make_combos.output,
    output:
        tbl="results/tables/snv_per_kbp.tbl",
        wide="results/tables/snv_per_kbp_wide.tbl",
        html=report(
            directory("results/tables/"), category="Tables", htmlindex="index.html"
        ),
    log:
        "logs/table_snps_per_kbp.log",
    conda:
        "../envs/R.yml"
    threads: 8
    script:
        "../scripts/per_kbp_stats.R"
