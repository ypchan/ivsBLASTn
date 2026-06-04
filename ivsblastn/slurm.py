from __future__ import annotations

import re
import shlex
import subprocess
from pathlib import Path
from typing import Optional

from .logging import LOG


def render_slurm_array_script(
    *,
    chunks_dir: Path,
    outdir: Path,
    db: Path,
    taxonomy: Optional[Path],
    threads: int,
    top_subjects: int,
    blast_max_hsps: int,
    cpus_per_task: int,
    mem: str,
    time: str,
    partition: Optional[str],
    array_concurrency: Optional[int],
    extra_run_args: str,
) -> str:
    """Render a Slurm array script that runs one ivsBLASTn chunk per task."""

    manifest = chunks_dir / "chunks.tsv"
    with manifest.open("rt", encoding="utf-8") as handle:
        chunk_count = max(0, sum(1 for _line in handle) - 1)
    if chunk_count <= 0:
        raise ValueError(f"No chunks found in manifest: {manifest}")
    concurrency = f"%{array_concurrency}" if array_concurrency else ""
    partition_line = f"#SBATCH --partition={partition}\n" if partition else ""
    db_arg = shlex.quote(str(db))
    taxonomy_line = f"  --taxonomy {shlex.quote(str(taxonomy))} \\\n" if taxonomy else ""
    blast_hsps_suffix = " \\" if extra_run_args else ""
    extra_line = f"  {extra_run_args}\n" if extra_run_args else ""
    return f"""#!/usr/bin/env bash
#SBATCH --job-name=ivsBLASTn
#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --mem={mem}
#SBATCH --time={time}
{partition_line}#SBATCH --output={outdir}/slurm/logs/%A_%a.out
#SBATCH --error={outdir}/slurm/logs/%A_%a.err
#SBATCH --array=1-{chunk_count}{concurrency}

set -euo pipefail

MANIFEST=\"{manifest}\"
CHUNK_FASTA=$(awk -v task_id=\"${{SLURM_ARRAY_TASK_ID}}\" 'NR == task_id + 1 {{print $2}}' \"$MANIFEST\")
CHUNK_NAME=$(basename \"$CHUNK_FASTA\" .fa)
CHUNK_OUTDIR=\"{outdir}/chunks/${{CHUNK_NAME}}\"

ivsBLASTn run \\
  --query \"$CHUNK_FASTA\" \\
  --db {db_arg} \\
{taxonomy_line}\
  --outdir \"$CHUNK_OUTDIR\" \\
  --threads {threads} \\
  --top-subjects {top_subjects} \\
  --blast-max-hsps {blast_max_hsps}{blast_hsps_suffix}
{extra_line}\
"""


def write_slurm_array_script(script: str, outdir: Path) -> Path:
    """Write a Slurm script under outdir/slurm."""

    slurm_dir = outdir / "slurm"
    (slurm_dir / "logs").mkdir(parents=True, exist_ok=True)
    script_path = slurm_dir / "ivsBLASTn_array.sbatch"
    script_path.write_text(script, encoding="utf-8")
    LOG.info("Wrote Slurm script: %s", script_path)
    return script_path


def submit_sbatch(script_path: Path) -> str:
    """Submit a Slurm script and return the parsed job id."""

    completed = subprocess.run(["sbatch", str(script_path)], check=True, text=True, capture_output=True)
    output = completed.stdout.strip()
    match = re.search(r"Submitted batch job\s+(\S+)", output)
    job_id = match.group(1) if match else output
    LOG.info("Submitted Slurm job: %s", job_id)
    return job_id
