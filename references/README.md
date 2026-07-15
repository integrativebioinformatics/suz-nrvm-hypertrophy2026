# References Directory

This directory contains genome FASTA and GTF annotation files required by the pipeline.

## Default References

The pipeline will automatically download these files if they don't already exist:
- `Rattus_norvegicus.GRCr8.dna.toplevel.fa` — Ensembl rat genome FASTA
- `Rattus_norvegicus.GRCr8.115.gtf` — Ensembl rat GTF annotation (release 115)

These are downloaded by `scripts/bash/00_download_references.sh` via `curl` from Ensembl FTP.

## Custom References

To use different genome/GTF files, pass them as CLI arguments to the pipeline:

```bash
./suz_pipeline.bash \
  --genome-fasta /path/to/custom_genome.fasta \
  --gtf /path/to/custom_annotation.gtf \
  ...
```

The defaults are used only if these flags are omitted.

## Files Not in Git

These reference files are **NOT committed to git** (see `reference_manifest.tsv`):
- Genome FASTA (too large)
- GTF annotation (large, regenerable from pipeline)
- TFLink All interactions (too large, provided separately by coauthor)

All are managed by:
1. Automatic download script (`00_download_references.sh`)
2. User-provided paths (via CLI flags)
3. External manifest tracking (in `reference_manifest.tsv`)
