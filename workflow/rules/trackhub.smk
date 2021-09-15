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
    longLabel {sm}_{h} snv density in 1 kbp windows
    visibility full
	autoScale Off
	maxHeightPixels 128:36:16
	graphTypeDefault Bar
	gridDefault OFF
	windowingFunction Mean
	color {strong_color}
	altColor {weak_color}
	viewLimits 0:20
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
        cat {input.bed} \
            | csvtk cut -C "$" -tT -f "#chr",start,end,num_snv  \
            > {output.bg}

        bedGraphToBigWig {output.bg} {input.fai} {output.bigwig}
        """


rule make_trackdb:
    input:
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
        for sm in tbl.index:
            for h in [1, 2]:
                out.write(
                    track.format(
                        sm=sm, h=h, strong_color=strong_color, weak_color=weak_color
                    )
                )
        out.close()
        open(output.hub, "w").write(hub)
        open(output.genomes, "w").write(genomes)
