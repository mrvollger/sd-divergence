track_db_header = """
track snv_density
compositeTrack off
shortLabel SNV desity 
longLabel  SNV density of HPRC year 1 samples in 1 kbp windows 
visibility hide
priority 30
type bigWig 0 {window_size}
maxItems 100000
group varRep

"""

hub = """
hub SNV_density
shortLabel SNV density
longLabel SNV density of HPRC year 1 samples in 1 kbp windows
genomesFile genomes.txt
email mvollger.edu
"""
genomes = """
genome t2t-chm13-v1.1
trackDb trackDb.chm13.txt
"""

track = """
	track {sm}_{h}_snv_density
	parent snv_density
	bigDataUrl SNVdensity/{sm}_{h}.bigWig
	shortLabel {sm}_{h} snv density
	longLabel {sm}_{h}
	visibility full
	autoScale Off
	maxHeightPixels 128:16:10
	graphTypeDefault Bar
	gridDefault OFF
	windowingFunction Mean
	color {strong_color}
	altColor {weak_color}
	viewLimits 0:5
	type bigWig 0 1000

"""


rule make_bigwig:
    input:
        bed=rules.add_snv_to_windows.output,
        fai=fai,
    output:
        bigwig="results/trackHub/SNVdensity/{sm}_{h}.bigWig",
        bg=temp("temp/trackHub/SNVdensity/{sm}_{h}.bg"),
    threads: 1
    conda:
        "../envs/env.yml"
    log:
        "logs/trackHub/SNVdesnisty/{sm}_{h}.bigWig.log",
    shell:
        """
        csvtk cut -C "$" -tT \
            -f "#chr",start,end,num_snv  \
            {input.bed} \
            > {output.bg}

        bedGraphToBigWig {output.bg} {input.fai} {output.bigwig}
        """


rule all_bigwig:
    input:
        bed="results/long_windows_with_snv_dist_annotation.bed.gz",
        fai=fai,
    output:
        bigwig="results/trackHub/SNVdensity/all_1.bigWig",
        bg="results/trackHub/SNVdensity/all_1.bg.bed",
    threads: 1
    conda:
        "../envs/env.yml"
    log:
        "logs/trackHub/SNVdesnisty/all_1.bigWig.log",
    params:
        window_size=config["window_size"],
        step_size=config["step_size"],
    shell:
        """
        zcat {input.bed} \
            | csvtk -tT -C "$" cut -f "#chr",start,end,hap_count,num_snv \
            | bedtools merge -i - -d -{params.step_size} -c 4,5 -o distinct,sum \
            | awk -v OFS=$'\t' '$3-$2 >= {params.step_size} {{print $1,$2,$3,$5/$4}}' \
        > {output.bg}

        bedGraphToBigWig {output.bg} {input.fai} {output.bigwig}
        """


rule make_trackdb:
    input:
        allbg=rules.all_bigwig.output.bg,
        bigwig=expand(rules.make_bigwig.output, sm=tbl.index, h=[1, 2]),
    output:
        track="results/trackHub/trackDb.chm13.txt",
        hub="results/trackHub/hub.txt",
        genomes="results/trackHub/genomes.txt",
    threads: 1
    log:
        "logs/trackHub/trackHub.log",
    run:
        strong_color = "175,4,4"
        weak_color = "47,79,79"
        out = open(output.track, "w")
        out.write(track_db_header.format(window_size=config["window_size"]))
        for sm in list(tbl.index) + ["all"]:
            for h in [1, 2]:
                if sm == "all" and h == 2:
                    continue
                out.write(
                    track.format(
                        sm=sm, h=h, strong_color=strong_color, weak_color=weak_color
                    )
                )
        out.close()
        open(output.hub, "w").write(hub)
        open(output.genomes, "w").write(genomes)
