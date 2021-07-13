#!/usr/bin/env bash
set -euo pipefail


printf "sample\th1_aln\th1_callable\th2_aln\th2_callable\tsnv\n"
ls -d /net/eichler/vol27/projects/hprc/nobackups/variants/pav/202103/pav_mm2_chm13/results/* | \
  parallel -n 1 \
  'printf "\
{/.}\
\t{}/align/pre-cut/aligned_tig_h1.bed.gz\t{}/callable/callable_regions_h1_500.bed.gz\
\t{}/align/pre-cut/aligned_tig_h2.bed.gz\t{}/callable/callable_regions_h2_500.bed.gz\
\t{}/bed/snv_snv.bed.gz\n\
"' \
  | sort 
