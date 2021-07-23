#
# summary table
#
rule make_summary_table_per_haplotype:
    input:
        #rules.all_snv.output,
        rules.annotate_snv.output,
        rules.make_combos.output,
    output:
        tbl="results/tables/{sm}_{h}/{sm}_{h}_snv_per_kbp.tbl",
        wide="results/tables/{sm}_{h}/{sm}_{h}_snv_per_kbp_wide.tbl",
    log:
        "logs/{sm}_{h}_snps_per_kbp.log",
    conda:
        "../envs/R.yml"
    threads: 8
    script:
        "../scripts/per_hap_kbp_stats.R"


rule make_summary_table:
    input:
        rules.all_snv.output,
        rules.make_combos.output,
    output:
        tbl="results/tables/snv_per_kbp.tbl",
        wide="results/tables/snv_per_kbp_wide.tbl",
        html=report(
            directory("results/html/"), category="Tables", htmlindex="index.html"
        ),
    log:
        "logs/table_snps_per_kbp.log",
    conda:
        "../envs/R.yml"
    threads: 8
    script:
        "../scripts/per_kbp_stats.R"
