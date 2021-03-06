include: "snv.smk"


VCF = "/net/eichler/vol27/projects/hprc/nobackups/data_table/preqc/chm13/vcf/variants_freeze1_snv_snv_alt.vcf.gz"
REF = "/net/eichler/vol26/projects/chm13_t2t/nobackups/assemblies/chm13_v1.1_plus38Y.fasta"
ANCESTRAL = "/net/eichler/vol28/projects/long_read_archive/nobackups/nhp/Clint_PTR/assemblies/hifiasm/0.15.2/Clint_PTR.hifiasm.bp.hap1.p_ctg.fasta"
SAM = "/net/eichler/vol26/projects/primate_sv/nobackups/nhp_sd_pav/temp/Clint_PTR/align/pre-cut/aligned_tig_h1.sam.gz"
PAIRS = "temp/mutypes.txt"

if not os.path.exists(PAIRS):
    shell(
        f"""
samtools view  \
	<(zcat {SAM}) \
	| cut -f 1,3 \
	| sort | uniq \
	> {PAIRS}
"""
    )

pairs = [line.strip().split() for line in open(PAIRS)]


rule make_psl:
    input:
        sam=SAM,
    output:
        psl=temp("temp/mutyper/psl/clint.psl"),
    log:
        "logs/mutyper/psl.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        module load ucsc
        samtools view -b \
            <(zcat {input.sam}) \
            | bamToPsl /dev/stdin \
            {output}  
        """


rule make_chain:
    input:
        psl=rules.make_psl.output.psl,
    output:
        chain_out_to_ref=temp("temp/mutyper/chain/out-to-ref/{out}-{ref}.chain"),
        chain=temp("temp/mutyper/chain/ref-to-out/{ref}-{out}.chain"),
    log:
        "logs/mutyper/chain/chain_{ref}-{out}.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        module load ucsc
        grep -w {wildcards.ref} {input.psl} \
            | grep -w {wildcards.out} \
            | pslToChain /dev/stdin {output.chain_out_to_ref}

        chainSwap {output.chain_out_to_ref} /dev/stdout \
            | chainSort /dev/stdin /dev/stdout \
            | sed '/^$/d' \
            > {output.chain}
        """


rule setup_vcf:
    input:
        vcf=VCF,
    output:
        bcf=temp("temp/mutyper/all.bcf"),
        csi=temp("temp/mutyper/all.bcf.csi"),
    log:
        "logs/mutyper/vcf.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        bcftools sort \
            -O b -m 8G \
            {input.vcf} \
            > {output.bcf}

        bcftools index -f {output.bcf}
        """


rule subset_vcf:
    input:
        bcf=rules.setup_vcf.output.bcf,
        chain=rules.make_chain.output.chain_out_to_ref,
    output:
        rgn=temp("temp/mutyper/subset/rgn/{ref}-{out}.rgn"),
        bcf=temp("temp/mutyper/subset/vcf/{ref}-{out}.bcf"),
    log:
        "logs/mutyper/vcf/{ref}-{out}.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        grep "^chain" {input.chain} \
             | grep -w {wildcards.ref} \
             | grep -w {wildcards.out} \
             | cut -d " " -f 3,6,7 \
             | awk '{{print $1"\t"$2"\t"$3 }}' \
             | bedtools sort -i - \
             | bedtools merge -i - \
             > {output.rgn}

        bcftools view {input.bcf} \
             --regions-file {output.rgn} \
            | bcftools sort -m 8G - \
            | bcftools +fill-tags \
            -Ob -o {output.bcf}
        """


rule prep_ancestor:
    input:
        reference=REF,
        outgroup=ANCESTRAL,
    output:
        ref=temp("temp/mutyper/input_fasta/ref_{ref}-{out}.fa"),
        fai_ref=temp("temp/mutyper/input_fasta/ref_{ref}-{out}.fa.fai"),
        out=temp("temp/mutyper/input_fasta/out_{out}-{ref}.fa"),
        fai_out=temp("temp/mutyper/input_fasta/out_{out}-{ref}.fa.fai"),
    log:
        "logs/mutyper/ancestral_fasta/{ref}-{out}.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        samtools faidx {input.reference} {wildcards.ref} \
            | seqtk seq -l 60 > {output.ref}
        sleep 5s; samtools faidx {output.ref}

        samtools faidx {input.outgroup} {wildcards.out} \
            | seqtk seq -l 60 > {output.out}
        sleep 5s; samtools faidx {output.out}
        """


rule make_ancestor:
    input:
        bcf=rules.subset_vcf.output.bcf,
        rgn=rules.subset_vcf.output.rgn,
        chain=rules.make_chain.output.chain,
        ref=rules.prep_ancestor.output.ref,
        out=rules.prep_ancestor.output.out,
    output:
        fasta=temp("temp/mutyper/ancestral_fasta/ref-to-out/{ref}-{out}.fa"),
        fai=temp("temp/mutyper/ancestral_fasta/ref-to-out/{ref}-{out}.fa.fai"),
    log:
        "logs/mutyper/ancestral_fasta/ref-to-out/{ref}-{out}.log",
    conda:
        "../envs/mutyper.yml"
    shell:
        """
        rm {output}
        mutyper ancestor \
            --bed {input.rgn} \
            {input.bcf} \
            {input.ref} \
            {input.out} \
            {input.chain} \
        {output.fasta}
        samtools faidx {output.fasta}
        """


rule annotate_vcf:
    input:
        bcf=rules.subset_vcf.output.bcf,
        fasta=rules.make_ancestor.output.fasta,
    output:
        bcf=temp("temp/mutyper/vcf/ref-to-out/{ref}-{out}.bcf"),
        csi=temp("temp/mutyper/vcf/ref-to-out/{ref}-{out}.bcf.csi"),
    log:
        "logs/mutyper/annotate_vcf/ref-to-out/{ref}-{out}.log",
    conda:
        "../envs/mutyper.yml"
    shell:
        """
        mutyper variants {input.fasta} {input.bcf} \
            | bcftools sort -Ob -m 8G - \
            > {output.bcf}
        bcftools index -f {output.bcf}
        """


###############################################################################
###############################################################################
###############################################################################
rule make_ancestor_out_to_ref:
    input:
        bcf=rules.subset_vcf.output.bcf,
        rgn=rules.subset_vcf.output.rgn,
        chain=rules.make_chain.output.chain_out_to_ref,
        ref=rules.prep_ancestor.output.ref,
        out=rules.prep_ancestor.output.out,
    output:
        fasta=temp("temp/mutyper/ancestral_fasta/out-to-ref/{out}-{ref}.fa"),
        fai=temp("temp/mutyper/ancestral_fasta/out-to-ref/{out}-{ref}.fa.fai"),
    log:
        "logs/mutyper/ancestral_fasta/out-to-ref/{ref}-{out}.log",
    conda:
        "../envs/mutyper.yml"
    shell:
        """
        rm {output}
        mutyper ancestor \
            --bed {input.rgn} \
            {input.bcf} \
            {input.ref} \
            {input.out} \
            {input.chain} \
        {output.fasta}
        samtools faidx {output.fasta}
        """


rule annotate_vcf_out_to_ref:
    input:
        bcf=rules.subset_vcf.output.bcf,
        fasta=rules.make_ancestor_out_to_ref.output.fasta,
    output:
        bcf=temp("temp/mutyper/vcf/out-to-ref/{out}-{ref}.bcf"),
        csi=temp("temp/mutyper/vcf/out-to-ref/{out}-{ref}.bcf.csi"),
    log:
        "logs/mutyper/annotate_vcf/out-to-ref/{ref}-{out}.log",
    conda:
        "../envs/mutyper.yml"
    shell:
        """
        mutyper variants {input.fasta} {input.bcf} \
            | bcftools sort -Ob -m 8G - \
            > {output.bcf}
        bcftools index -f {output.bcf}
        """


###############################################################################
###############################################################################
###############################################################################


def get_mutyper_rtn(wc):
    for out, ref in pairs:
        if ref == "*" or out == "*" or out == "h1tg000047l":
            continue
        yield (rules.annotate_vcf.output.bcf).format(ref=ref, out=out)


rule annotated_vcf:
    input:
        bcf=get_mutyper_rtn,
    output:
        bcf="results/mutyper/anno_vcf.bcf",
        csi="results/mutyper/anno_vcf.bcf.csi",
    log:
        "logs/mutyper/annotated_vcf.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        bcftools concat -Ob -a \
            {input.bcf} > {output.bcf}
        bcftools index -f {output.bcf}
        """


rule mutyper_spectra:
    input:
        bcf=rules.annotated_vcf.output.bcf,
    output:
        spectra="results/mutyper/spectra.txt",
    log:
        "logs/mutyper/spectra.log",
    conda:
        "../envs/mutyper.yml"
    shell:
        """
        mutyper spectra {input.bcf} > {output.spectra}
        """


rule mutyper:
    input:
        rules.annotated_vcf.output,
        rules.mutyper_spectra.output,
