#!/usr/bin/env bash
set -euo pipefail
snakemake --configfile config/real_config.yaml --use-conda --cores 200 \
  -p --notemp \
  $@
exit
--resources load=1000 \
