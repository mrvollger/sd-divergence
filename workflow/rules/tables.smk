#
# summary table
#
# rule region_sizes:
#     input:
#         rules.long_and_filtered_windows.output,
#     output:
#         "results/tables/region_sizes.tbl",
#     log:
#         "logs/region_sizes.log",
#     conda:
#         "envs/R.yml"
#     threads: 8
#     script:
#         "scripts/region_sizes.R"
rule make_summary_table:
    input:
        rules.all_snv.output,
        rules.long_and_filtered_windows.output,
    output:
        tbl="results/tables/snv_per_kbp.tbl",
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
