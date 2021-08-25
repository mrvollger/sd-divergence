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
        chain_out_to_ref=temp(
            "temp/mutyper/chain/temp/outgroup-to-reference_{rn}-{an}.chain"
        ),
        chain=temp("temp/mutyper/chain/reference-to-outgroup_{rn}-{an}.chain"),
    log:
        "logs/mutyper/chain/chain_{rn}-{an}.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        module load ucsc
        grep -w {wildcards.rn} {input.psl} \
            | grep -w {wildcards.an} \
            | pslToChain /dev/stdin {output.chain_out_to_ref}

        chainSwap {output.chain_out_to_ref} /dev/stdout \
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
        rgn=temp("temp/mutyper/rgn/{rn}-{an}.rgn"),
        bcf=temp("temp/mutyper/vcf/{rn}-{an}.bcf"),
    log:
        "logs/mutyper/vcf/{rn}-{an}.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        grep "^chain" {input.chain} \
             | grep -w {wildcards.rn} \
             | grep -w {wildcards.an} \
             | cut -d " " -f 3,6,7 \
             | awk '{{print $1"\t"$2"\t"$3 }}' \
             | bedtools sort -i - \
             | bedtools merge -i - \
             > {output.rgn}
         cat {output.rgn}

        bcftools view {input.bcf} \
             --regions-file {output.rgn} \
            | bcftools sort -m 8G - \
            | bcftools +fill-tags \
            -Ob -o {output.bcf}
        """


rule make_ancestor:
    input:
        bcf=rules.subset_vcf.output.bcf,
        chain=rules.make_chain.output.chain,
        reference=REF,
        ancestor=ANCESTRAL,
    output:
        an=temp("temp/mutyper/ancestral_fasta/rn/an_{rn}-{an}.fa"),
        fai_an=temp("temp/mutyper/ancestral_fasta/rn/an_{rn}-{an}.fa.fai"),
        rn=temp("temp/mutyper/ancestral_fasta/an/rn_{rn}-{an}.fa"),
        fai_rn=temp("temp/mutyper/ancestral_fasta/an/rn_{rn}-{an}.fa.fai"),
        fasta=temp("temp/mutyper/ancestral_fasta/{rn}-{an}.fa"),
        fai=temp("temp/mutyper/ancestral_fasta/{rn}-{an}.fa.fai"),
    log:
        "logs/mutyper/ancestral_fasta/{rn}-{an}.log",
    conda:
        "../envs/mutyper.yml"
    shell:
        """
        samtools faidx {input.reference} {wildcards.rn} | seqtk seq -l 60 > {output.rn}
        samtools faidx {input.ancestor} {wildcards.an} | seqtk seq -l 60 > {output.an}
        samtools faidx {output.rn}
        samtools faidx {output.an}

        mutyper ancestor \
            {input.bcf} \
            {output.rn} \
            {output.an} \
            {input.chain} \
         {output.fasta}
         samtools faidx {output.fasta}
        """


rule annotate_vcf:
    input:
        bcf=rules.subset_vcf.output.bcf,
        fasta=rules.make_ancestor.output.fasta,
    output:
        bcf=temp("temp/mutyper/anno_vcf/{rn}-{an}.bcf"),
        csi=temp("temp/mutyper/anno_vcf/{rn}-{an}.bcf.csi"),
    log:
        "logs/mutyper/annotate_vcf/{rn}-{an}.log",
    conda:
        "../envs/mutyper.yml"
    shell:
        """
        mutyper variants {input.fasta} {input.bcf} \
            | bcftools sort -Ob -m 8G - \
            > {output.bcf}
        bcftools index -f {output.bcf}
        """


def get_mutyper_rtn(wc):
    for an, rn in pairs:
        if rn == "*" or an == "*" or an == "h1tg000047l":
            continue
        yield (rules.annotate_vcf.output.bcf).format(rn=rn, an=an)


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
