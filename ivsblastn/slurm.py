from __future__ import annotations

import re
import shlex
import subprocess
from pathlib import Path
from typing import List, Optional, Sequence

from .logging import LOG


def render_sbatch_header(
    *,
    job_name: str,
    cpus_per_task: int,
    mem: str,
    time: str,
    output: Path,
    error: Path,
    array: str,
    partition: Optional[str] = None,
    account: Optional[str] = None,
    qos: Optional[str] = None,
    nodes: Optional[int] = None,
    ntasks: Optional[int] = None,
    constraint: Optional[str] = None,
    gres: Optional[str] = None,
    exclude: Optional[str] = None,
    nodelist: Optional[str] = None,
    extra_sbatch_options: Optional[Sequence[str]] = None,
) -> str:
    """Render SBATCH resource directives."""

    lines = [
        "#!/usr/bin/env bash",
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --cpus-per-task={cpus_per_task}",
        f"#SBATCH --mem={mem}",
        f"#SBATCH --time={time}",
    ]
    optional_options = [
        ("partition", partition),
        ("account", account),
        ("qos", qos),
        ("nodes", nodes),
        ("ntasks", ntasks),
        ("constraint", constraint),
        ("gres", gres),
        ("exclude", exclude),
        ("nodelist", nodelist),
    ]
    for name, value in optional_options:
        if value is not None:
            lines.append(f"#SBATCH --{name}={value}")
    lines.extend(
        [
            f"#SBATCH --output={output}",
            f"#SBATCH --error={error}",
            f"#SBATCH --array={array}",
        ]
    )
    for option in extra_sbatch_options or []:
        stripped = option.strip()
        if not stripped:
            continue
        if stripped.startswith("#SBATCH"):
            lines.append(stripped)
        else:
            lines.append(f"#SBATCH {stripped}")
    return "\n".join(lines)


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
    account: Optional[str] = None,
    qos: Optional[str] = None,
    nodes: Optional[int] = None,
    ntasks: Optional[int] = None,
    constraint: Optional[str] = None,
    gres: Optional[str] = None,
    exclude: Optional[str] = None,
    nodelist: Optional[str] = None,
    extra_sbatch_options: Optional[Sequence[str]] = None,
    array_concurrency: Optional[int] = None,
    extra_run_args: str = "",
    run_args: Optional[Sequence[str]] = None,
) -> str:
    """Render a Slurm array script that runs one ivsBLASTn chunk per task."""

    manifest = chunks_dir / "chunks.tsv"
    with manifest.open("rt", encoding="utf-8") as handle:
        chunk_count = max(0, sum(1 for _line in handle) - 1)
    if chunk_count <= 0:
        raise ValueError(f"No chunks found in manifest: {manifest}")
    concurrency = f"%{array_concurrency}" if array_concurrency else ""
    array_spec = f"1-{chunk_count}{concurrency}"
    manifest_arg = shlex.quote(str(manifest.resolve()))
    outdir_abs = outdir.resolve()
    outdir_arg = shlex.quote(str(outdir_abs))
    db_arg = shlex.quote(str(db.resolve()))
    header = render_sbatch_header(
        job_name="ivsBLASTn",
        cpus_per_task=cpus_per_task,
        mem=mem,
        time=time,
        output=outdir_abs / "slurm" / "logs" / "%A_%a.out",
        error=outdir_abs / "slurm" / "logs" / "%A_%a.err",
        array=array_spec,
        partition=partition,
        account=account,
        qos=qos,
        nodes=nodes,
        ntasks=ntasks,
        constraint=constraint,
        gres=gres,
        exclude=exclude,
        nodelist=nodelist,
        extra_sbatch_options=extra_sbatch_options,
    )
    command_lines: List[str] = [
        "ivsBLASTn run",
        '  --query "$CHUNK_FASTA"',
        f"  --db {db_arg}",
    ]
    if taxonomy:
        command_lines.append(f"  --taxonomy {shlex.quote(str(taxonomy.resolve()))}")
    command_lines.extend(
        [
            '  --outdir "$CHUNK_OUTDIR"',
            f"  --threads {threads}",
            f"  --top-subjects {top_subjects}",
            f"  --blast-max-hsps {blast_max_hsps}",
        ]
    )
    for item in run_args or []:
        command_lines.append(f"  {item}")

    rendered_command = ""
    if extra_run_args:
        rendered_command = " \\\n".join(command_lines) + " \\\n" + f"  {extra_run_args}\n"
    else:
        rendered_command = " \\\n".join(command_lines) + "\n"

    return f"""{header}

set -euo pipefail

MANIFEST={manifest_arg}
CHUNK_FASTA=$(awk -F '\\t' -v task_id=\"${{SLURM_ARRAY_TASK_ID}}\" 'NR == task_id + 1 {{print $2}}' \"$MANIFEST\")
if [[ -z \"$CHUNK_FASTA\" ]]; then
  echo \"No chunk FASTA found for task ${{SLURM_ARRAY_TASK_ID}} in $MANIFEST\" >&2
  exit 2
fi

CHUNK_BASE=$(basename \"$CHUNK_FASTA\")
CHUNK_NAME=\"${{CHUNK_BASE%.gz}}\"
CHUNK_NAME=\"${{CHUNK_NAME%.fasta}}\"
CHUNK_NAME=\"${{CHUNK_NAME%.fa}}\"
CHUNK_NAME=\"${{CHUNK_NAME%.fna}}\"
CHUNK_OUTDIR={outdir_arg}/chunks/${{CHUNK_NAME}}

{rendered_command}\
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
