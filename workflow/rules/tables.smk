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
        long=expand(
            rules.make_summary_table_per_haplotype.output.tbl, sm=tbl.index, h=[1, 2]
        ),
        wide=expand(
            rules.make_summary_table_per_haplotype.output.wide, sm=tbl.index, h=[1, 2]
        ),
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


rule make_r_data:
    input:
        tbl=rules.make_summary_table.output.tbl,
        snv=rules.small_snv.output,
        long=rules.long_and_filtered_windows.output,
    output:
        windows="results/Rdata/windows.Rdata",
        tbl="results/Rdata/tbl.Rdata",
        snv="results/Rdata/snv.Rdata",
    log:
        "logs/rdata.log",
    conda:
        "../envs/R.yml"
    threads: 8
    script:
        "../scripts/make_r_data.R"
