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


rule make_chain:
    input:
        sam=SAM,
    output:
        chain="temp/mutyper/clint.chain",
    log:
        "logs/mutyper/chain.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        module load ucsc
        samtools view -b \
            <(zcat {input.sam}) \
            | bamToPsl /dev/stdin /dev/stdout \
            | pslToChain /dev/stdin \
            {output}  
        """


rule make_ancestor:
    input:
        vcf=VCF,
        sam=SAM,
        chain=rules.make_chain.output.chain,
        reference=REF,
        ancestor=ANCESTRAL,
    output:
        an="temp/mutyper/ancestral_fasta/{an}.fa",
        rn="temp/mutyper/ancestral_fasta/{rn}.fa",
        fasta="temp/mutyper/ancestral_fasta/{rn}-{an}.fa",
    log:
        "logs/mutyper/ancestral_fasta/{rn}-{an}.log",
    conda:
        "../envs/mutyper.yml"
    shell:
        """
        samtools faidx {input.reference} {wildcards.rn} > {output.rn}
        samtools faidx {input.ancestor} {wildcards.an} > {output.an}

        mutyper ancestor \
            {input.vcf} \
            {output.rn} \
            {output.an} \
             {input.chain} \
         {output.fasta}
        """


def get_mutyper_fastas(wc):
    for an, rn in pairs:
        if rn == "*" or an == "*":
            continue
        yield (rules.make_ancestor.output.fasta).format(rn=rn, an=an)


rule mutyper_setup:
    input:
        fastas=get_mutyper_fastas,
