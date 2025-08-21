#!/usr/bin/env python3

import argparse
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Dict, Optional, Tuple, List, Iterable

from ete3 import Tree
import pandas as pd
import sys


# =====================
# Layout & resolution
# =====================
DEFAULT_MARKERS = ("TerL", "MCP", "portal")


@dataclass
class DefaultLayout:
    markers: Tuple[str, ...] = DEFAULT_MARKERS

    def resolve_for_marker(self, base: Path, marker: str) -> Tuple[Path, Path]:
        tree = base / "5_phylogenies" / "3_iToL" / f"{marker}_iToL.nwk"
        out = tree.with_name(f"{marker}_iToL_renamed.nwk")
        return tree, out

    def resolve_misc(self, base: Path) -> Tuple[Path, Path]:
        rename_table = base / "1_rename" / "prophage.table"
        func_in = base / "4_ORF" / "3_functional_annot_tables"
        func_out = base / "4_ORF" / "3_functional_annot_table_renamed.tsv"
        return rename_table, func_in, func_out
    
    def resolve_marker_fasta(self, base: Path, marker: str) -> Tuple[Path, Path]:
        """
        Input:  5_phylogenies/0_marker_genes/1_final/{marker}.faa
        Output: 5_phylogenies/0_marker_genes/1_final/{marker}_renamed.faa
        """
        in_faa = base / "5_phylogenies" / "0_marker_genes" / "1_final" / f"{marker}.faa"
        out_faa = in_faa.with_name(f"{marker}_renamed.faa")
        return in_faa, out_faa


# =====================
# Argument parsing
# =====================

def parse_arguments() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Rename leaves of phylogenetic tree(s) using a rename table. "
            "Supports CrassUS layout via --crassus_output and runs automatically for markers: TerL, MCP, portal. "
            "You can override any path or limit to a single explicit tree with --tree_file/--output_tree."
        )
    )

    # Convenience base path
    p.add_argument(
        "--crassus_output",
        type=Path,
        help=(
            "Base directory of CrassUS results (…/results/prophage_analysis). "
            "If provided, the script will rename all markers (TerL, MCP, portal) by default."
        ),
    )

    # Markers control (only used with --crassus_output)
    p.add_argument(
        "--markers",
        type=str,
        default=",".join(DEFAULT_MARKERS),
        help=(
            "Comma‑separated list of markers to process when --crassus_output is used. "
            f"Default: {','.join(DEFAULT_MARKERS)}"
        ),
    )

    p.add_argument(
        "--rename_all_codings",
        action="store_true",
        help=(
            "Also rename files in 4_ORF/0_all_codings → 4_ORF/0_all_codings_renamed "
            "(filenames to old_id; also rewrite .gff seqid & seqhdr)."
        ),
    )
    p.add_argument(
        "--rewrite_fasta_headers",
        action="store_true",
        help=(
            "When used with --rename_all_codings, also rewrite FASTA headers in .faa/.fna "
            "from new_id to old_id (only affects leading identifier or leading 'new_id|…')."
        ),
    )

    p.add_argument(
        "--rename_best_coding",
        action="store_true",
        help=(
            "Also rename files in 4_ORF/1_best_coding → 4_ORF/1_best_coding_renamed "
            "(filenames to old_id; also rewrite .gff seqid & seqhdr)."
        ),
    )

    # Individual paths (override or single‑file mode)
    p.add_argument("--tree_file", type=Path, help="Input tree file (Newick). If set, processes only this file.")
    p.add_argument("--output_tree", type=Path, help="Path to save the renamed tree (required if --tree_file).")
    p.add_argument("--rename_table", type=Path, help="TSV with columns: new_id and old_id.")
    p.add_argument("--func_in", type=Path, help="Directory with .table files to rename (optional).")
    p.add_argument("--func_out", type=Path, help="Path to save concatenated & renamed .table file (optional).")

    # Misc
    p.add_argument("--dry_run", action="store_true", help="Print resolved paths and exit.")
    p.add_argument("--strict", action="store_true", help="Error if a required path is missing.")
    return p.parse_args()


# =====================
# Utilities
# =====================

def ensure_parent_exists(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def error_or_warn(msg: str, strict: bool) -> None:
    if strict:
        raise FileNotFoundError(msg)
    print(msg, file=sys.stderr)


def require_file_exists(path: Path, desc: str, strict: bool) -> None:
    if not path or not path.is_file():
        error_or_warn(f"[ERROR] {desc} not found: {path}", strict)


def require_dir_exists(path: Path, desc: str, strict: bool) -> None:
    if not path or not path.is_dir():
        error_or_warn(f"[ERROR] {desc} directory not found: {path}", strict)


# =====================
# Core renaming logic
# =====================

def load_rename_dict(rename_table: Optional[Path]) -> Dict[str, str]:
    if rename_table is None:
        return {}
    df = pd.read_csv(rename_table, sep="\t", dtype=str)
    if not {"new_id", "old_id"}.issubset(df.columns):
        raise ValueError("Rename table must contain 'new_id' and 'old_id' columns.")
    return dict(zip(df["new_id"], df["old_id"]))


def rename_tree_leaves(tree: Tree, rename_dict: Dict[str, str]) -> Tuple[int, int, List[str]]:
    renamed = 0
    not_renamed = 0
    warnings: List[str] = []
    for leaf in tree.iter_leaves():
        parts = leaf.name.split("|")
        if parts:
            new_id = parts[0]
            if new_id in rename_dict:
                parts[0] = rename_dict[new_id]
                leaf.name = "|".join(parts)
                renamed += 1
            else:
                not_renamed += 1
                warnings.append(f"[WARNING] ID '{new_id}' not found in rename table for leaf: {leaf.name}")
        else:
            not_renamed += 1
            warnings.append(f"[WARNING] Unexpected leaf name format: {leaf.name}")
    return renamed, not_renamed, warnings


def save_warnings(warnings: List[str], output_tree: Path) -> None:
    log_file = output_tree.with_suffix(".warnings.log")
    ensure_parent_exists(log_file)
    with open(log_file, "w") as f:
        for line in warnings:
            f.write(line + "\n")


def rename_functional_tables(func_in: Path, func_out: Path, rename_dict: Dict[str, str]) -> None:
    table_files = sorted(func_in.glob("*.table"))
    if not table_files:
        print(f"[WARNING] No .table files found in {func_in}")
        return

    all_rows: List[pd.DataFrame] = []
    header = None

    for i, file in enumerate(table_files):
        df = pd.read_csv(file, sep="\t", dtype=str)
        if i == 0:
            header = df.columns.tolist()
        else:
            df.columns = header  # enforce consistent columns

        # Rename genome
        if "genome" in df.columns:
            df["genome"] = df["genome"].apply(lambda x: rename_dict.get(x, x))

        # Rename protein_id (only part before first "|")
        if "protein_id" in df.columns:
            def rename_protein_id(pid: str) -> str:
                if isinstance(pid, str) and "|" in pid:
                    parts = pid.split("|")
                    parts[0] = rename_dict.get(parts[0], parts[0])
                    return "|".join(parts)
                else:
                    return rename_dict.get(str(pid), str(pid))

            df["protein_id"] = df["protein_id"].apply(rename_protein_id)

        all_rows.append(df)

    result_df = pd.concat(all_rows, ignore_index=True)
    ensure_parent_exists(func_out)
    result_df.to_csv(func_out, sep="\t", index=False)
    print(f"[INFO] Concatenated functional table written to: {func_out}")


# =====================
# Orchestration
# =====================

def process_single_tree(tree_path: Path, out_path: Path, rename_dict: Dict[str, str]) -> None:
    print(f"[INFO] Processing tree: {tree_path}")
    tree = Tree(str(tree_path), format=1)
    renamed, not_renamed, warnings = rename_tree_leaves(tree, rename_dict)
    ensure_parent_exists(out_path)
    tree.write(outfile=str(out_path), format=1)
    save_warnings(warnings, out_path)
    print(f"[INFO] → Saved: {out_path} | renamed={renamed}, not_renamed={not_renamed}")


def _parse_all_codings_filename(path: Path) -> Optional[re.Match]:
    """
    Match files like {new_id}_tbl-{code}.{ext}
    where code ∈ {11, TAG, TGA} and ext ∈ {gff, faa, fna}.
    """
    return re.match(r"^(?P<new_id>.+?)_tbl-(?P<code>11|TAG|TGA)\.(?P<ext>gff|faa|fna)$", path.name)


def _rewrite_gff_text(text: str, new_id: str, old_id: str) -> str:
    """
    - In '# Sequence Data:' line, replace seqhdr="new_id" → seqhdr="old_id"
    - In feature lines (non-#), replace first column new_id → old_id
    """
    out_lines: List[str] = []
    seqhdr_re = re.compile(r'(seqhdr=")([^"]+)(")')
    for line in text.splitlines():
        if line.startswith("# Sequence Data:"):
            line = seqhdr_re.sub(
                lambda m: m.group(1) + (old_id if m.group(2) == new_id else m.group(2)) + m.group(3),
                line,
            )
            out_lines.append(line)
            continue
        if line.startswith("#") or not line.strip():
            out_lines.append(line)
            continue
        parts = line.split("\t")
        if parts and parts[0] == new_id:
            parts[0] = old_id
            line = "\t".join(parts)
        out_lines.append(line)
    return "\n".join(out_lines) + ("\n" if text.endswith("\n") else "")


def _rewrite_fasta_headers_text(text: str, new_id: str, old_id: str) -> str:
    """
    Change leading identifier in FASTA headers:
      >new_id
      >new_id|something
      >new_id_1
    → old_id[…], old_id|something, old_id_1  (preserve suffix and rest of line)
    """
    out_lines: List[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            m = re.match(r"^(>)(\S+)(.*)$", line)
            if m:
                prefix, ident, rest = m.groups()
                if ident == new_id:
                    ident = old_id
                elif ident.startswith(new_id + "|"):
                    ident = old_id + ident[len(new_id):]
                elif ident.startswith(new_id + "_"):
                    ident = old_id + ident[len(new_id):]
                line = f"{prefix}{ident}{rest}"
        out_lines.append(line)
    return "\n".join(out_lines) + ("\n" if text.endswith("\n") else "")

def _rewrite_marker_faa_headers_text(text: str, rename_dict: Dict[str, str]) -> str:
    """
    Headers look like:
      >new_id|{protein_length}|{gene_ordinal_number}
    Replace only the first field (new_id) using rename_dict, preserve the rest.
    """
    out_lines: List[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            # Split header at first whitespace; we only rename the identifier chunk
            m = re.match(r"^(>)(\S+)(.*)$", line)
            if m:
                prefix, ident, rest = m.groups()
                # ident is something like new_id|len|ord
                parts = ident.split("|")
                if parts:
                    parts[0] = rename_dict.get(parts[0], parts[0])
                    ident = "|".join(parts)
                line = f"{prefix}{ident}{rest}"
        out_lines.append(line)
    return "\n".join(out_lines) + ("\n" if text.endswith("\n") else "")


def rename_marker_proteins(base: Path, markers: Iterable[str], rename_dict: Dict[str, str]) -> None:
    """
    For each marker in markers:
      Input : 5_phylogenies/0_marker_genes/1_final/{marker}.faa
      Output: 5_phylogenies/0_marker_genes/1_final/{marker}_renamed.faa
    Rewrite headers: {new_id}|len|ord  ->  {old_id}|len|ord
    """
    for marker in markers:
        in_faa = base / "5_phylogenies" / "0_marker_genes" / "1_final" / f"{marker}.faa"
        out_faa = in_faa.with_name(f"{marker}_renamed.faa")

        if not in_faa.is_file():
            print(f"[WARNING] Marker FASTA not found for '{marker}': {in_faa}")
            continue

        try:
            text = in_faa.read_text()
        except Exception as e:
            print(f"[ERROR] Failed to read {in_faa}: {e}")
            continue

        new_text = _rewrite_marker_faa_headers_text(text, rename_dict)
        ensure_parent_exists(out_faa)
        out_faa.write_text(new_text)
        print(f"[INFO] Renamed marker proteins → {out_faa}")

def rename_all_codings(base: Path, rename_dict: Dict[str, str], rewrite_fasta_headers: bool = False) -> None:
    """
    Read from 4_ORF/0_all_codings, write to 4_ORF/0_all_codings_renamed:
    - Filenames: {new_id}_tbl-{code}.{ext} → {old_id}_tbl-{code}.{ext}
    - .gff content: change seqid (col 1) and seqhdr="..."
    - .faa/.fna content: optionally rewrite headers if rewrite_fasta_headers=True
    """
    in_dir = base / "4_ORF" / "0_all_codings"
    out_dir = base / "4_ORF" / "0_all_codings_renamed"
    out_dir.mkdir(parents=True, exist_ok=True)

    if not in_dir.exists() or not in_dir.is_dir():
        print(f"[WARNING] Skipping 0_all_codings: missing dir {in_dir}")
        return

    processed = skipped = missing = 0

    for path in sorted(in_dir.iterdir()):
        if not path.is_file():
            continue
        m = _parse_all_codings_filename(path)
        if not m:
            skipped += 1
            continue

        new_id = m.group("new_id")
        code = m.group("code")
        ext  = m.group("ext")

        old_id = rename_dict.get(new_id, new_id)
        if old_id == new_id:
            missing += 1
            print(f"[WARNING] No mapping for {new_id} (filename: {path.name}); keeping new_id in name/content.")

        out_name = f"{old_id}_tbl-{code}.{ext}"
        out_path = out_dir / out_name

        text = path.read_text()
        if ext == "gff":
            text = _rewrite_gff_text(text, new_id, old_id)
        elif rewrite_fasta_headers and ext in {"faa", "fna"}:
            text = _rewrite_fasta_headers_text(text, new_id, old_id)

        out_path.write_text(text)
        processed += 1
        print(f"[INFO] Wrote {out_path}")

    print(f"[INFO] 0_all_codings: processed={processed}, skipped(non-matching)={skipped}, missing_mappings={missing}")

def rename_best_coding(base: Path, rename_dict: Dict[str, str], rewrite_fasta_headers: bool = False) -> None:
    """
    Read from 4_ORF/1_best_coding, write to 4_ORF/1_best_coding_renamed:
    - Filenames: {new_id}_tbl-{code}.{ext} → {old_id}_tbl-{code}.{ext}
    - .gff content: change seqid (col 1) and seqhdr="..."
    - .faa/.fna content: optionally rewrite headers if rewrite_fasta_headers=True
    """
    in_dir = base / "4_ORF" / "1_best_coding"
    out_dir = base / "4_ORF" / "1_best_coding_renamed"
    out_dir.mkdir(parents=True, exist_ok=True)

    if not in_dir.exists() or not in_dir.is_dir():
        print(f"[WARNING] Skipping 1_best_coding: missing dir {in_dir}")
        return

    processed = skipped = missing = 0

    for path in sorted(in_dir.iterdir()):
        if not path.is_file():
            continue
        m = _parse_all_codings_filename(path)
        if not m:
            skipped += 1
            continue

        new_id = m.group("new_id")
        code = m.group("code")
        ext  = m.group("ext")

        old_id = rename_dict.get(new_id, new_id)
        if old_id == new_id:
            missing += 1
            print(f"[WARNING] No mapping for {new_id} (filename: {path.name}); keeping new_id in name/content.")

        out_name = f"{old_id}_tbl-{code}.{ext}"
        out_path = out_dir / out_name

        text = path.read_text()
        if ext == "gff":
            text = _rewrite_gff_text(text, new_id, old_id)
        elif rewrite_fasta_headers and ext in {"faa", "fna"}:
            text = _rewrite_fasta_headers_text(text, new_id, old_id)

        out_path.write_text(text)
        processed += 1
        print(f"[INFO] Wrote {out_path}")

    print(f"[INFO] 1_best_coding: processed={processed}, skipped(non-matching)={skipped}, missing_mappings={missing}")



def main() -> None:
    args = parse_arguments()

    # Mode A: explicit single tree
    if args.tree_file:
        if not args.output_tree:
            print("[ERROR] --output_tree is required when --tree_file is provided.", file=sys.stderr)
            sys.exit(2)

        rename_table = args.rename_table
        if not rename_table and not args.crassus_output:
            print("[ERROR] Provide --rename_table or --crassus_output to locate the rename table.", file=sys.stderr)
            sys.exit(2)

        # If base is given, fill missing misc paths from it
        if args.crassus_output and not rename_table:
            layout = DefaultLayout()
            rename_table, _, _ = layout.resolve_misc(args.crassus_output)

        # Report paths
        print("[INFO] Resolved paths (single tree mode):")
        print(f"  tree_file   : {args.tree_file}")
        print(f"  output_tree : {args.output_tree}")
        print(f"  rename_table: {rename_table}")
        print(f"  func_in     : {args.func_in}")
        print(f"  func_out    : {args.func_out}")

        if args.dry_run:
            print("[INFO] Dry run complete. Exiting.")
            return

        require_file_exists(args.tree_file, "Tree file", strict=args.strict)
        require_file_exists(rename_table, "Rename table", strict=args.strict)

        rename_dict = load_rename_dict(rename_table)
        process_single_tree(args.tree_file, args.output_tree, rename_dict)

        # Optional functional tables
        if args.func_in and args.func_out:
            if args.func_in.is_dir():
                rename_functional_tables(args.func_in, args.func_out, rename_dict)
            else:
                print(f"[WARNING] Skipping functional tables: func_in not a directory: {args.func_in}", file=sys.stderr)
        return

    # Mode B: multi‑marker mode via --crassus_output
    if not args.crassus_output:
        print(
            "[ERROR] Either provide --tree_file/--output_tree (single mode) or --crassus_output (multi‑marker mode).",
            file=sys.stderr,
        )
        sys.exit(2)

    base = args.crassus_output
    markers: Iterable[str] = [m.strip() for m in args.markers.split(",") if m.strip()]
    layout = DefaultLayout(markers=tuple(markers))

    # Resolve shared resources
    rename_table, func_in, func_out = layout.resolve_misc(base)

    # Allow overrides
    if args.rename_table:
        rename_table = args.rename_table
    if args.func_in:
        func_in = args.func_in
    if args.func_out:
        func_out = args.func_out

    # Resolve trees per marker
    trees: List[Tuple[str, Path, Path]] = []  # (marker, tree_path, out_path)
    for m in layout.markers:
        t, o = layout.resolve_for_marker(base, m)
        trees.append((m, t, o))

    # Report paths
    print("[INFO] Resolved paths (multi‑marker mode):")
    print(f"  base        : {base}")
    print(f"  markers     : {', '.join(layout.markers)}")
    print(f"  rename_table: {rename_table}")
    print(f"  func_in     : {func_in}")
    print(f"  func_out    : {func_out}")
    for m, t, o in trees:
        print(f"  [{m}] tree_file -> {t}")
        print(f"  [{m}] output   -> {o}")

    if args.dry_run:
        print("[INFO] Dry run complete. Exiting.")
        return

    # Validate inputs
    require_file_exists(rename_table, "Rename table", strict=args.strict)
    for m, t, _ in trees:
        if not t.is_file():
            error_or_warn(f"[ERROR] Tree file for marker '{m}' not found: {t}", args.strict)

    # Load once
    rename_dict = load_rename_dict(rename_table)

    # Process each marker that has an existing tree
    for m, t, o in trees:
        if t.is_file():
            process_single_tree(t, o, rename_dict)
        else:
            print(f"[WARN] Skipping '{m}' (tree file missing): {t}")

    # Optional functional tables once per run
    if func_in and func_out and func_in.is_dir():
        rename_functional_tables(func_in, func_out, rename_dict)
    else:
        print(f"[WARNING] Skipping functional tables (missing or not a dir): {func_in}")

    # Optionally rename 4_ORF/0_all_codings files
    if args.rename_all_codings:
        rename_all_codings(base, rename_dict, rewrite_fasta_headers=args.rewrite_fasta_headers)
    
    if args.rename_best_coding:
        rename_best_coding(base, rename_dict, rewrite_fasta_headers=args.rewrite_fasta_headers)

    # Rename marker protein FASTAs in 5_phylogenies/0_marker_genes/1_final
    rename_marker_proteins(base, layout.markers, rename_dict)


if __name__ == "__main__":
    main()
