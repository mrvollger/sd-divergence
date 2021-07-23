#
# plotting
#
rule plot_cumulative_divergence:
    input:
        rules.long_and_filtered_windows.output,
    output:
        report(
            "results/figures/cumulative_divergence.svg",
            category="Figures",
        ),
        report(
            "results/figures/cumulative_divergence_per_hap.svg",
            category="Figures",
        ),
    log:
        "logs/cumulative_divergence.log",
    conda:
        "../envs/R.yml"
    threads: 8
    script:
        "../scripts/divergence_cum.R"


rule plot_callable_space:
    input:
        windows=rules.long_and_filtered_windows.output,
        beds=expand(rules.syntenic_and_callable.output, sm=tbl.index, h=[1, 2]),
    output:
        report(
            "results/figures/callable_space.svg",
            category="Figures",
        ),
    log:
        "logs/callable_space.log",
    conda:
        "../envs/R.yml"
    threads: 8
    script:
        "../scripts/called_regions.R"


rule plot_violin:
    input:
        rules.make_summary_table.output.tbl,
        config["metadata"],
    output:
        report(
            "results/figures/violin.svg",
            category="Figures",
        ),
    log:
        "logs/violin.log",
    conda:
        "../envs/R.yml"
    threads: 8
    script:
        "../scripts/violin_plots.R"


rule plots:
    input:
        rules.plot_violin.output,
        rules.plot_callable_space.output,
        rules.plot_cumulative_divergence.output,
