#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path
from tqdm import tqdm
from typing import Dict, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Get CrassUS genomes and convert headers from new_id to old_id based on mapping table."
    )
    parser.add_argument(
        "--input_dir", type=Path, required=True,
        help="Path to CrassUS tool output directory"
    )
    parser.add_argument(
        "--output_directory", type=Path, required=True,
        help="Where to save the converted .fna file and log"
    )
    parser.add_argument(
        "--write_warnings", action="store_true",
        help="If set, log warnings for missing old_id mappings"
    )
    return parser.parse_args()


def setup_logging(log_path: Path) -> None:
    logging.basicConfig(
        filename=log_path,
        filemode='w',
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=logging.INFO
    )


def load_id_mapping(table_file: Path) -> Dict[str, str]:
    mapping = {}
    with table_file.open() as f:
        next(f)  # skip header
        for line_num, line in enumerate(f, start=2):
            stripped = line.strip()
            if not stripped:
                continue
            parts = stripped.split("\t")
            if len(parts) != 2:
                continue  # skip malformed
            old, new = parts
            if new:
                mapping[new] = old
    return mapping


def convert_fasta_headers(
    fasta_files: list[Path],
    id_map: Dict[str, str],
    output_fasta: Path,
    write_warnings: bool
) -> Tuple[int, int]:
    converted_count = 0
    skipped_count = 0

    with output_fasta.open("w") as out_f:
        for fasta_file in tqdm(fasta_files, desc="Renaming genomes", unit="genome"):
            new_id = fasta_file.stem
            if new_id not in id_map:
                if write_warnings:
                    logging.warning(f"Missing old_id for {new_id} (no mapping found)")
                skipped_count += 1
                continue

            old_id = id_map[new_id]
            header = f">{old_id}"
            seq_lines = []
            with fasta_file.open() as f:
                for line in f:
                    line = line.strip()
                    if not line.startswith(">"):  # skip original header
                        seq_lines.append(line)
            out_f.write(f"{header}\n" + "\n".join(seq_lines) + "\n")
            converted_count += 1

    return converted_count, skipped_count


def main(input_dir: Path, output_dir: Path, write_warnings: bool) -> None:
    rename_table = input_dir / "1_rename"
    table_files = list(rename_table.glob("*.table"))
    if not table_files:
        raise FileNotFoundError(f"No .table file found in {rename_table}")
    table_file = table_files[0]

    fasta_dir = input_dir / "3_crass_contigs"
    all_fasta_files = list(fasta_dir.glob("*.fasta"))
    if not all_fasta_files:
        raise FileNotFoundError(f"No .fasta files found in {fasta_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)
    log_file = output_dir / "conversion.log"
    setup_logging(log_file)

    id_map = load_id_mapping(table_file)

    # Only keep FASTA files whose stem (new_id) is in the mapping
    valid_fasta_files = [f for f in all_fasta_files if f.stem in id_map]

    output_fasta = output_dir / f"{input_dir.name}_crassus_genomes.fna"
    converted_count, skipped_count = convert_fasta_headers(valid_fasta_files, id_map, output_fasta, write_warnings)

    summary_msg = (
        f"Conversion complete. Output: {output_fasta}\n"
        f"Converted sequences: {converted_count}\n"
        f"Skipped sequences (no matching new_id in mapping or not renamed): {skipped_count}\n"
        f"Log file: {log_file}"
    )

    print(summary_msg)
    logging.info(summary_msg)


if __name__ == "__main__":
    args = parse_args()
    main(args.input_dir, args.output_directory, args.write_warnings)
