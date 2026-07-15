#!/usr/bin/env bash
# Example call for the future final pipeline.
# This is not yet the final executable pipeline; it defines the intended interface.

./suz_pipeline.bash \
  --longrna_samples metadata/longRNA_samples.tsv \
  --smallrna_samples metadata/smallRNA_samples.tsv \
  --fastq_manifest metadata/fastq_manifest.tsv \
  --genome references/Rattus_norvegicus.GRCr8.dna.toplevel.fa \
  --gtf references/Rattus_norvegicus.GRCr8.115.gtf \
  --tflink references/TFLink_All.tsv \
  --threads 16 \
  --out results
