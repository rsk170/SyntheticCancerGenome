#!/usr/bin/env python3
"""Build cumulative per-clone somatic mutation lists for tumor simulation."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build final per-clone SCG mutation files from converted hg38 clone "
            "mutation files. Descendant clones receive ancestor/founding mutations "
            "plus their own clone-specific mutations."
        )
    )
    parser.add_argument("manifest", type=Path, help="patient_manifest.csv from Step 2/4")
    parser.add_argument(
        "--out-dir",
        type=Path,
        help="Output directory; defaults to MANIFEST_DIR/final_clone_mutations",
    )
    parser.add_argument(
        "--founding-clone",
        help="Founding clone ID; auto-detected from clone_type=founding when omitted",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing final mutation files and summaries",
    )
    return parser.parse_args()


def relpath(path: Path, start: Path) -> str:
    try:
        return str(path.resolve().relative_to(start.resolve()))
    except ValueError:
        return str(path.resolve())


def resolve_path(value: str | Path, repo_root: Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else repo_root / path


def read_manifest(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise SystemExit(f"Manifest does not exist: {path}")
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    if not rows:
        raise SystemExit(f"Manifest has no rows: {path}")

    required = {"patient_id", "clone_id", "clone_fraction"}
    missing = required - set(rows[0])
    if missing:
        raise SystemExit(f"Manifest missing columns: {', '.join(sorted(missing))}")
    return rows


def detect_founding_clone(rows: list[dict[str, str]], requested: str | None) -> str:
    if requested:
        return requested

    founding = [row["clone_id"] for row in rows if row.get("clone_type") == "founding"]
    if len(founding) == 1:
        return founding[0]
    if len(founding) > 1:
        raise SystemExit(f"Multiple founding clones in manifest: {', '.join(sorted(founding))}")

    parentless = [row["clone_id"] for row in rows if not row.get("parent_clone_id")]
    if len(parentless) == 1:
        return parentless[0]
    raise SystemExit("Could not auto-detect founding clone. Pass --founding-clone.")


def default_mutation_path(manifest_dir: Path, clone_id: str) -> Path:
    return manifest_dir / "hg38_mutations" / f"{clone_id}.hg38.scg_mutations"


def mutation_path_for_clone(
    clone_id: str,
    row_by_clone: dict[str, dict[str, str]],
    manifest_dir: Path,
    repo_root: Path,
) -> Path:
    row = row_by_clone.get(clone_id)
    if row and row.get("clone_hg38_mutations_path"):
        return resolve_path(row["clone_hg38_mutations_path"], repo_root)
    return default_mutation_path(manifest_dir, clone_id)


def lineage_for_clone(
    clone_id: str,
    *,
    row_by_clone: dict[str, dict[str, str]],
    founding_clone: str,
    manifest_dir: Path,
    seen: tuple[str, ...] = (),
) -> list[str]:
    if clone_id in seen:
        chain = " -> ".join(seen + (clone_id,))
        raise SystemExit(f"Cycle in clone lineage: {chain}")

    row = row_by_clone.get(clone_id)
    if row is None:
        if clone_id == founding_clone and default_mutation_path(manifest_dir, clone_id).exists():
            return [clone_id]
        raise SystemExit(f"Clone {clone_id} is needed as an ancestor but is not in the manifest")

    parent = row.get("parent_clone_id", "").strip()
    if not parent:
        return [clone_id]
    return lineage_for_clone(
        parent,
        row_by_clone=row_by_clone,
        founding_clone=founding_clone,
        manifest_dir=manifest_dir,
        seen=seen + (clone_id,),
    ) + [clone_id]


def read_mutations(path: Path) -> list[str]:
    if not path.exists():
        raise SystemExit(f"Missing clone mutation file: {path}")
    with path.open() as handle:
        return [line.strip() for line in handle if line.strip()]


def write_mutations(path: Path, mutations: list[str], overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise SystemExit(f"Output already exists: {path}. Use --overwrite to replace it.")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(mutations) + ("\n" if mutations else ""))


def main() -> int:
    args = parse_args()
    repo_root = Path.cwd()
    manifest = args.manifest.resolve()
    manifest_dir = manifest.parent
    rows = read_manifest(manifest)
    row_by_clone = {row["clone_id"]: row for row in rows}
    if len(row_by_clone) != len(rows):
        raise SystemExit("Manifest contains duplicate clone_id values")

    patient_ids = {row["patient_id"] for row in rows}
    if len(patient_ids) != 1:
        raise SystemExit(f"Manifest contains multiple patient IDs: {', '.join(sorted(patient_ids))}")
    patient_id = next(iter(patient_ids))

    founding_clone = detect_founding_clone(rows, args.founding_clone)
    out_dir = (args.out_dir or manifest_dir / "final_clone_mutations").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_path = out_dir / "final_clone_mutation_summary.tsv"
    final_manifest_path = out_dir / "patient_manifest.final_clone_mutations.csv"
    if (summary_path.exists() or final_manifest_path.exists()) and not args.overwrite:
        raise SystemExit(f"Summary outputs already exist in {out_dir}. Use --overwrite to replace them.")

    summary_rows: list[dict[str, str]] = []
    final_manifest_rows: list[dict[str, str]] = []

    for row in rows:
        clone_id = row["clone_id"]
        lineage = lineage_for_clone(
            clone_id,
            row_by_clone=row_by_clone,
            founding_clone=founding_clone,
            manifest_dir=manifest_dir,
        )
        if founding_clone not in lineage:
            lineage = [founding_clone] + lineage

        seen: set[str] = set()
        final_mutations: list[str] = []
        source_paths: list[Path] = []
        source_count = 0
        duplicate_count = 0

        for lineage_clone in lineage:
            mutation_path = mutation_path_for_clone(lineage_clone, row_by_clone, manifest_dir, repo_root)
            source_paths.append(mutation_path)
            for mutation in read_mutations(mutation_path):
                source_count += 1
                if mutation in seen:
                    duplicate_count += 1
                    continue
                seen.add(mutation)
                final_mutations.append(mutation)

        output_path = out_dir / f"{clone_id}.final.hg38.scg_mutations"
        write_mutations(output_path, final_mutations, args.overwrite)

        summary_rows.append(
            {
                "patient_id": patient_id,
                "clone_id": clone_id,
                "clone_type": row.get("clone_type", ""),
                "parent_clone_id": row.get("parent_clone_id", ""),
                "clone_fraction": row.get("clone_fraction", ""),
                "lineage_clone_ids": ";".join(lineage),
                "source_mutation_files": ";".join(relpath(path, repo_root) for path in source_paths),
                "source_mutations": str(source_count),
                "final_mutations": str(len(final_mutations)),
                "duplicates_removed": str(duplicate_count),
                "final_clone_hg38_mutations_path": relpath(output_path, repo_root),
            }
        )

        final_row = dict(row)
        final_row["final_clone_hg38_mutations_path"] = relpath(output_path, repo_root)
        final_row["final_clone_hg38_mutation_count"] = str(len(final_mutations))
        final_row["final_clone_lineage"] = ";".join(lineage)
        final_manifest_rows.append(final_row)

    with summary_path.open("w", newline="") as handle:
        fieldnames = [
            "patient_id",
            "clone_id",
            "clone_type",
            "parent_clone_id",
            "clone_fraction",
            "lineage_clone_ids",
            "source_mutation_files",
            "source_mutations",
            "final_mutations",
            "duplicates_removed",
            "final_clone_hg38_mutations_path",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    with final_manifest_path.open("w", newline="") as handle:
        fieldnames = list(final_manifest_rows[0].keys())
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(final_manifest_rows)

    print(f"Wrote final clone mutation files: {out_dir}")
    print(f"Wrote summary: {summary_path}")
    print(f"Wrote final clone manifest: {final_manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
