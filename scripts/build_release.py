"""Build an annual curated pydipole data release."""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path
from typing import Optional

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pydipole.__main__ import to_tex
from pydipole.utils import SEP

from validate_data import REFERENCES_DIR, TABLE_FIELDS, iter_submission_files, load_record, validate_all

RELEASES_DIR = ROOT / "data" / "releases"
PACKAGE_DATA_DIR = ROOT / "src" / "pydipole" / "data"


def record_to_row(record: dict) -> dict:
    return {target: record[source] for target, source in TABLE_FIELDS.items()}


def record_release_year(record: dict, default_year: int) -> int:
    return int(record.get("release_year", default_year))


def collect_rows(year: int, after_year: Optional[int] = None) -> list[dict]:
    rows = []
    for path in iter_submission_files():
        record = load_record(path)
        release_year = record_release_year(record, year)
        if release_year > year or (after_year is not None and release_year <= after_year):
            continue
        source_row = record.get("source_row")
        sort_key = int(source_row) if source_row is not None else 1_000_000 + release_year
        rows.append((sort_key, str(path), record_to_row(record)))
    return [row for _, _, row in sorted(rows)]


def latest_previous_release(year: int) -> Optional[Path]:
    candidates = []
    for path in RELEASES_DIR.iterdir():
        if path.is_dir() and path.name.isdigit() and int(path.name) < year:
            candidates.append(path)
    return max(candidates, key=lambda path: int(path.name), default=None)


def load_base_data(year: int):
    previous_release = latest_previous_release(year)
    if previous_release is not None:
        return (
            pd.read_csv(previous_release / "database.csv", sep=SEP),
            int(previous_release.name),
        )
    return pd.DataFrame(columns=TABLE_FIELDS.keys()), None


def collect_reference_paths() -> list[Path]:
    return sorted(path for path in REFERENCES_DIR.glob("*.bib") if path.is_file())


def write_references(reference_paths: list[Path], output_path: Path) -> None:
    seen = set()
    contents = []
    for path in reference_paths:
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        contents.append(path.read_text().strip())
        contents.append("")
    output_path.write_text("\n".join(contents).rstrip() + "\n")


def build_release(year: int, sync_package_data: bool = False) -> Path:
    errors = validate_all()
    if errors:
        raise ValueError("\n".join(errors))

    release_dir = RELEASES_DIR / str(year)
    release_dir.mkdir(parents=True, exist_ok=True)

    base_df, base_year = load_base_data(year)
    rows = collect_rows(year, after_year=base_year)
    if rows:
        base_df = pd.concat([base_df, pd.DataFrame(rows)], ignore_index=True)

    database_path = release_dir / "database.csv"
    references_path = release_dir / "references.bib"
    table_path = release_dir / "table.tex"

    base_df.to_csv(database_path, sep=SEP, index=False)
    write_references(collect_reference_paths(), references_path)
    table_df = base_df.copy()
    table_df["Atom"] = table_df["Atom"].mask(table_df["Atom"].duplicated(), "")
    table_df["Z"] = table_df["Z"].mask(table_df["Z"].duplicated(), "")
    to_tex(table_df, table_path, label=f"tab:table_{year}")

    if sync_package_data:
        shutil.copy(database_path, PACKAGE_DATA_DIR / "database.csv")
        shutil.copy(references_path, PACKAGE_DATA_DIR / "references.bib")

    return release_dir


def main(args: Optional[list] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--year", type=int, required=True, help="Release year to build.")
    parser.add_argument(
        "--sync-package-data",
        action="store_true",
        help="Also update src/pydipole/data with the generated release data.",
    )
    parsed_args = parser.parse_args(args)

    try:
        release_dir = build_release(parsed_args.year, parsed_args.sync_package_data)
    except ValueError as exc:
        print(exc)
        return 1

    print(f"Built release data in {release_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
