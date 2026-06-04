# ivsBLASTn

`ivsBLASTn` detects intervening sequences (IVSs) in 16S/SSU rRNA gene sequences using BLASTN HSP geometry.

The main evidence pattern is:

```text
query:    exon-left  [candidate IVS]  exon-right
subject:  exon-left  nearly continuous exon-right
```

In BLASTN output this appears as two HSPs on the query separated by a large query gap, while the subject-side coordinates remain adjacent or nearly adjacent.

## Install

Clone or download this repository, then install:

```bash
python -m pip install .
```

For development, install editable:

```bash
python -m pip install -e .
```

Check the command:

```bash
ivsBLASTn --help
ivsBLASTn run --help
```

## Update

From a git checkout:

```bash
git pull
python -m pip install -e .
```

If installed non-editably:

```bash
git pull
python -m pip install --upgrade .
```

## Requirements

Python dependencies are installed from `pyproject.toml`:

- `rich`
- `rich-argparse`

External tools are required when `ivsBLASTn` runs BLAST internally:

- `blastn`
- `makeblastdb`

Check them before a large run:

```bash
blastn -version
makeblastdb -version
```

If BLAST+ is not on `PATH`, pass explicit paths:

```bash
ivsBLASTn run ... \
  --blastn-bin /path/to/blastn \
  --makeblastdb-bin /path/to/makeblastdb
```

## Commands

```text
ivsBLASTn run           Run IVS detection on one query FASTA or one chunk
ivsBLASTn split         Split a large query FASTA into chunks
ivsBLASTn submit-slurm  Generate and optionally submit a Slurm array job
ivsBLASTn merge         Merge per-chunk outputs
```

Old single-command usage is still accepted:

```bash
ivsBLASTn --query query.fa --db ref_db --outdir out
```

Internally this is treated as:

```bash
ivsBLASTn run --query query.fa --db ref_db --outdir out
```

## Reference Data

### Starting From SILVA-Style FASTA

Use `--ref-fasta` when the reference FASTA header contains taxonomy after the sequence ID:

```text
>AB000393.1.1510 Bacteria;Pseudomonadota;...;Vibrio;Vibrio halioticoli
```

Example:

```bash
ivsBLASTn run \
  --query query_16s.fa \
  --ref-fasta SILVA_NR99.fa.gz \
  --outdir ivs_run \
  --threads 8
```

Reference preprocessing rules:

1. Keep only domains in `--ref-domains`, default `Archaea,Bacteria`.
2. Normalize sequence text, convert `U` to `T`, and keep only non-empty `A/T/G/C` sequences.
3. Require a clear species name in the seventh taxonomy field.
4. Skip records where the species field has only a genus name.
5. Skip records where the species epithet, the second word of the species name, contains digits.
6. Keep the longest `--ref-per-species` records per clear species, default `1`.

Accepted:

```text
Vibrio halioticoli
```

Skipped as unclear:

```text
Vibrio
Vibrio 1234
Vibrio sp001
```

Outputs under `outdir/reference/` include:

```text
raw_reference.fa
raw_reference.tax.tsv
raw_reference_db.*
```

### Starting From An Existing BLAST DB

Use this when you already have a curated reference DB:

```bash
ivsBLASTn run \
  --query query_16s.fa \
  --db reference_db_prefix \
  --taxonomy reference.tax.tsv \
  --outdir ivs_run \
  --threads 8
```

The taxonomy TSV should have two columns:

```text
subject_id<TAB>taxonomy
```

### Starting From An Existing BLAST Table

Use this when BLASTN was run externally:

```bash
ivsBLASTn run \
  --query query_16s.fa \
  --blast query_vs_reference.blastn.tsv \
  --taxonomy reference.tax.tsv \
  --outdir ivs_run
```

The BLAST table must be tab-delimited outfmt 6 with columns:

```text
qseqid sseqid pident length qstart qend sstart send evalue bitscore
```

## Key Parameters

For 16S/SSU IVS screening, these are usually the most important:

```bash
--top-subjects 100
--blast-max-hsps 5
--min-pident 75
--min-hsp-len 100
--min-intron-len 25
--max-intron-len 2000
--max-ref-gap 30
```

Query BLASTN uses:

```text
-max_target_seqs = --top-subjects
-max_hsps        = --blast-max-hsps
```

So the maximum reported HSP count per query is bounded by:

```text
--top-subjects * --blast-max-hsps
```

For IVS detection in 16S rRNA genes, IVSs are expected to be sparse, so the default `--blast-max-hsps 5` keeps BLAST output smaller than the earlier conservative value of 20.

## Single-Node Run

Recommended when the query FASTA is modest:

```bash
ivsBLASTn run \
  --query query_16s.fa \
  --db reference_db_prefix \
  --taxonomy reference.tax.tsv \
  --outdir ivs_run \
  --threads 8 \
  --top-subjects 100 \
  --blast-max-hsps 5
```

Output layout:

```text
ivs_run/
  blast/
    query_16s.vs_reference.blastn.tsv
  results/
    query_16s.summary.tsv
    query_16s.supporting_hsps.tsv
    query_16s.intron_free.fa
    query_16s.introns.fa
    query_16s.exons.bed
    query_16s.introns.bed
    query_16s.report.md
```

## Large Datasets And Slurm

For million-sequence query FASTA files, do not rely on one large `blastn -num_threads 96` process. BLASTN often does not scale linearly with CPU count.

The better strategy is query splitting:

```text
large query FASTA
  -> many chunk FASTA files
  -> Slurm array, one chunk per task
  -> merge outputs
```

The algorithm is query-independent, so splitting query records is safe. Each query is classified from its own BLAST HSPs, reference taxonomy, and parameters.

### 1. Split Query FASTA

```bash
ivsBLASTn split \
  --query all_16s.fa \
  --chunks-dir batch01/chunks \
  --chunk-size 5000
```

This writes:

```text
batch01/chunks/
  chunks.tsv
  query.000001.fa
  query.000002.fa
  ...
```

### 2. Generate Slurm Array Script

Dry run first. This writes a script but does not submit:

```bash
ivsBLASTn submit-slurm \
  --chunks-dir batch01/chunks \
  --db reference_db_prefix \
  --taxonomy reference.tax.tsv \
  --outdir batch01 \
  --threads 8 \
  --cpus-per-task 8 \
  --mem 16G \
  --time 12:00:00 \
  --array-concurrency 40
```

Review:

```text
batch01/slurm/ivsBLASTn_array.sbatch
```

Then submit manually:

```bash
sbatch batch01/slurm/ivsBLASTn_array.sbatch
```

Or submit directly:

```bash
ivsBLASTn submit-slurm \
  --chunks-dir batch01/chunks \
  --db reference_db_prefix \
  --taxonomy reference.tax.tsv \
  --outdir batch01 \
  --threads 8 \
  --cpus-per-task 8 \
  --mem 16G \
  --time 12:00:00 \
  --array-concurrency 40 \
  --submit
```

Typical starting points:

```text
chunk-size:          2,000-20,000 query records
threads/task:        8-16
cpus-per-task:       same as threads
array-concurrency:   20-100, depending on cluster policy
```

Prefer many moderate BLAST jobs over one 96-CPU BLAST job.

### 3. Merge Chunk Outputs

After all Slurm array tasks finish:

```bash
ivsBLASTn merge \
  --chunk-results-dir batch01/chunks \
  --outdir batch01/final \
  --label all_16s
```

Merged files are written under:

```text
batch01/final/results/
  all_16s.summary.tsv
  all_16s.supporting_hsps.tsv
  all_16s.intron_free.fa
  all_16s.introns.fa
  all_16s.exons.bed
  all_16s.introns.bed
  all_16s.merge_report.md
```

Do not use `--gzip-fasta-output` for chunk runs if you want direct FASTA merging.

## Interpreting Results

Start with `*.summary.tsv`.

Important columns:

```text
query_id
classification
confidence
intron_start
intron_end
intron_len
support_subjects
support_taxa_at_genus
median_subject_gap
median_pident
best_subject
best_subject_taxonomy
reasons
```

Confidence tiers:

```text
HIGH    many supporting subjects and enough taxonomic spread
MEDIUM  multiple supporting subjects and taxa
LOW     at least one subject supports the HSP-gap pattern
NONE    no supported IVS signal
```

Review candidates by checking:

1. `support_subjects`: more independent subjects is stronger.
2. `support_taxa_at_genus`: support across genera is stronger than one narrow group.
3. `median_subject_gap`: values near 0 are best.
4. `median_pident`: should be reasonable for the reference distance.
5. `intron_len`: very short or near `--max-intron-len` needs manual inspection.
6. `*.supporting_hsps.tsv`: supporting HSP pairs should agree on breakpoint coordinates.

Output FASTA files:

```text
*.intron_free.fa   all query sequences; IVS removed only for candidates passing --min-output-confidence
*.introns.fa       candidate IVS sequences passing --min-output-confidence
```

BED files:

```text
*.introns.bed      candidate IVS intervals
*.exons.bed        exon intervals after IVS removal
```

## Troubleshooting

Missing BLAST+:

```text
FileNotFoundError: blastn
```

Install BLAST+ or pass `--blastn-bin`.

Too slow on large input:

```text
Use split + submit-slurm + merge.
Lower --chunk-size if tasks run too long.
Use 8-16 CPUs per task before trying 96 CPUs.
```

Huge BLAST output:

```text
Lower --top-subjects or --blast-max-hsps.
Defaults are --top-subjects 100 and --blast-max-hsps 5.
```

Few MEDIUM/HIGH calls:

```text
Provide --taxonomy.
Increase --top-subjects if reference may contain IVS-bearing near neighbors.
Inspect *.supporting_hsps.tsv before changing confidence thresholds.
```

## Development

Run tests:

```bash
python -m unittest discover -s tests
```

Run the package without installing:

```bash
python -m ivsblastn --help
```

Repository layout:

```text
ivsblastn/
  cli.py
  algorithm.py
  blast.py
  fasta.py
  reference.py
  split.py
  slurm.py
  merge.py
  outputs.py
  taxonomy.py
  models.py
tests/
figures/
```

See [ivsBLASTn_technical_doc.md](ivsBLASTn_technical_doc.md) for the detailed algorithm reference.
