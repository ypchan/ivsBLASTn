# ivsBLASTn

`ivsBLASTn` detects intervening sequences in 16S/SSU rRNA sequences from BLASTN HSP geometry.

The core signal is a split alignment pattern:

```text
query:    HSP1 ---- candidate IVS ---- HSP2
subject:  HSP1 ---- nearly continuous - HSP2
```

If the query contains an IVS and the reference subject does not, BLASTN often reports two query HSPs separated by a large query gap while the subject coordinates remain nearly continuous.

## Installation

Install from this repository:

```bash
python -m pip install .
```

For development:

```bash
python -m pip install -e .
```

External BLAST+ tools are required for modes that run BLAST internally:

- `blastn`
- `makeblastdb`

Python dependencies are declared in `pyproject.toml`.

## Repository Layout

```text
src/ivsblastn/
  cli.py          command-line interface and workflow orchestration
  models.py       shared dataclasses
  fasta.py        FASTA, SILVA header, and species-name helpers
  reference.py    reference selection and optional self-cleaning
  blast.py        BLAST+ execution and outfmt 6 parsing
  algorithm.py    HSP-gap IVS detection logic
  outputs.py      TSV, FASTA, BED, report, and terminal summary writers
  paths.py        output path construction
  taxonomy.py     taxonomy parsing and rank helpers
```

The installed command is `ivsBLASTn`. The package can also be run as a module during development:

```bash
PYTHONPATH=src python -m ivsblastn --help
```

## Quick Start

Run with a SILVA-style reference FASTA:

```bash
ivsBLASTn \
  --query query_16s.fa \
  --ref-fasta SILVA_NR99.fa.gz \
  --outdir ivs_detection \
  --threads 8
```

Run with an existing BLAST database:

```bash
ivsBLASTn \
  --query query_16s.fa \
  --db reference_db_prefix \
  --taxonomy reference.tax.tsv \
  --outdir ivs_detection \
  --threads 8
```

Run from an existing BLAST outfmt 6 table:

```bash
ivsBLASTn \
  --query query_16s.fa \
  --blast query_vs_ref.blastn.tsv \
  --taxonomy reference.tax.tsv \
  --outdir ivs_detection
```

The BLAST table must contain these columns in order:

```text
qseqid sseqid pident length qstart qend sstart send evalue bitscore
```

## Reference Selection

In `--ref-fasta` mode, reference records are selected as follows:

1. Keep only allowed domains, default `Archaea,Bacteria`.
2. Normalize sequence text and keep only non-empty `A/T/G/C` sequences.
3. Require a clear species name in the taxonomy field.
4. Skip records where the species field has only a genus name.
5. Skip records where the species epithet, the second word of the species name, contains digits.
6. Keep the longest `--ref-per-species` records per clear species, default `1`.

Example accepted species field:

```text
Vibrio halioticoli
```

Examples skipped as unclear:

```text
Vibrio
Vibrio 1234
Vibrio sp001
```

## Main Outputs

For a query FASTA named `sample.fa`, outputs are written under:

```text
ivs_detection/
  reference/
  blast/
  results/
    sample.summary.tsv
    sample.supporting_hsps.tsv
    sample.intron_free.fa
    sample.introns.fa
    sample.exons.bed
    sample.introns.bed
    sample.report.md
```

The summary table reports per-query classification, IVS coordinates, confidence, support counts, subject gap statistics, and supporting taxonomy.

## Important Parameters

- `--ref-per-species`: maximum selected references per clear species, default `1`.
- `--min-pident`: minimum HSP percent identity, default `75.0`.
- `--min-hsp-len`: minimum HSP length, default `100`.
- `--min-intron-len`: minimum query gap accepted as IVS, default `25`.
- `--max-intron-len`: maximum query gap accepted as IVS, default `2000`.
- `--max-ref-gap`: maximum absolute subject gap or overlap, default `30`.
- `--top-subjects`: number of top subjects analyzed per query, default `100`.
- `--threads`: Python worker threads and BLASTN `-num_threads`, default `4`.

## Large Datasets

For million-sequence query sets, do not rely on a single `blastn -num_threads` process scaling linearly. The algorithm is query-independent, so the preferred HPC strategy is:

1. Build or prepare the reference DB once.
2. Split the query FASTA into chunks.
3. Run each chunk against the same reference DB as a Slurm array job.
4. Merge per-chunk `summary.tsv`, `supporting_hsps.tsv`, FASTA, BED, and reports.

This repository currently contains the single-job core program. Slurm chunking and merge helpers should be added as a separate workflow layer.

## Documentation

See [ivsBLASTn_technical_doc.md](ivsBLASTn_technical_doc.md) for the full algorithm description, parameter reference, and output semantics.

## Tests

Run the current unit tests from a checkout:

```bash
PYTHONPATH=src python -m unittest discover -s tests
```
