#!/usr/bin/env bash
set -euo pipefail

mkdir -p data 

ls -d /net/eichler/vol27/projects/hprc/nobackups/variants/pav/202103/pav_mm2_chm13/results/* \
  | head -n 4 \
  | parallel -n 1 \
  'mkdir -p data/{/.} \
  && cp {}/align/pre-cut/aligned_tig_h1.bed.gz data/{/.}/. \
  && cp {}/callable/callable_regions_h1_500.bed.gz data/{/.}/. \
  && cp {}/align/pre-cut/aligned_tig_h2.bed.gz data/{/.}/. \
  && cp {}/callable/callable_regions_h2_500.bed.gz data/{/.}/. \
  && zcat {}/bed/snv_snv.bed.gz | head -n 1000 | gzip > data/{/.}/snv_snv.bed.gz
'


printf "sample\th1_aln\th2_aln\th1_callable\th2_callable\tsnv\n"

ls -d data/* \
  | parallel -n 1 --will-cite \
  'printf "{/.}\t" \
  && readlink -f {}/* | tr "\n" "\t" \
  && printf "\n"'
