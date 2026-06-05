from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Optional

from .logging import LOG


MERGE_SUFFIXES = [
    ".summary.tsv",
    ".supporting_hsps.tsv",
    ".ivs_free.fa",
    ".ivs_free.fa.gz",
    ".ivs.fa",
    ".ivs.fa.gz",
    ".ivs.bed",
    ".intron_free.fa",
    ".intron_free.fa.gz",
    ".introns.fa",
    ".introns.fa.gz",
    ".exons.bed",
    ".introns.bed",
]


def sorted_chunk_result_dirs(chunk_results_dir: Path, exclude: Optional[Path] = None) -> List[Path]:
    """Return chunk result directories in stable lexical order."""

    excluded = exclude.resolve() if exclude is not None else None
    chunk_dirs: List[Path] = []
    for path in chunk_results_dir.iterdir():
        if not path.is_dir():
            continue
        if excluded is not None and path.resolve() == excluded:
            continue
        chunk_dirs.append(path)
    return sorted(chunk_dirs)


def iter_files_by_suffix(chunk_dirs: Iterable[Path], suffix: str) -> Iterable[Path]:
    for chunk_dir in chunk_dirs:
        results_dir = chunk_dir / "results"
        if not results_dir.exists():
            continue
        yield from sorted(results_dir.glob(f"*{suffix}"))


def concatenate_text_files(files: Iterable[Path], output: Path, keep_one_header: bool) -> int:
    """Concatenate text files, optionally keeping only the first header line."""

    output.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    wrote_header = False
    with output.open("wt", encoding="utf-8") as out:
        for path in files:
            with path.open("rt", encoding="utf-8") as handle:
                for line_index, line in enumerate(handle):
                    if keep_one_header and line_index == 0:
                        if wrote_header:
                            continue
                        wrote_header = True
                    out.write(line)
                    written += 1
    return written


def concatenate_binary_files(files: Iterable[Path], output: Path) -> int:
    """Concatenate binary files and return the number of files written."""

    output.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with output.open("wb") as out:
        for path in files:
            with path.open("rb") as handle:
                while True:
                    chunk = handle.read(1024 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)
            written += 1
    return written


def merge_chunk_outputs(chunk_results_dir: Path, outdir: Path, label: str = "merged") -> List[Path]:
    """Merge standard ivsBLASTn per-chunk outputs into a final output directory."""

    chunk_dirs = sorted_chunk_result_dirs(chunk_results_dir, exclude=outdir)
    final_results = outdir / "results"
    final_results.mkdir(parents=True, exist_ok=True)
    outputs: List[Path] = []

    for suffix in MERGE_SUFFIXES:
        files = list(iter_files_by_suffix(chunk_dirs, suffix))
        if not files:
            continue
        output = final_results / f"{label}{suffix}"
        if suffix.endswith(".gz"):
            concatenate_binary_files(files, output)
            outputs.append(output)
            continue
        keep_one_header = suffix in {".summary.tsv", ".supporting_hsps.tsv"}
        concatenate_text_files(files, output, keep_one_header=keep_one_header)
        outputs.append(output)

    report = final_results / f"{label}.merge_report.md"
    with report.open("wt", encoding="utf-8") as handle:
        print("# ivsBLASTn merge report", file=handle)
        print("", file=handle)
        print(f"- Chunk result directory: `{chunk_results_dir}`", file=handle)
        print(f"- Chunks discovered: `{len(chunk_dirs)}`", file=handle)
        print(f"- Output directory: `{final_results}`", file=handle)
        print("", file=handle)
        print("## Merged files", file=handle)
        print("", file=handle)
        for output in outputs:
            print(f"- `{output}`", file=handle)
    outputs.append(report)
    LOG.info("Merged %s chunk directories into %s", len(chunk_dirs), final_results)
    return outputs
