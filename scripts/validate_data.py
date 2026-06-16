"""Validate accepted community data records."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Optional

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pydipole.periodic import Element

DATA_DIR = ROOT / "data"
SUBMISSIONS_DIR = DATA_DIR / "submissions"
REFERENCES_DIR = DATA_DIR / "references"
CATEGORIES = {"theoretical", "experimental"}
REQUIRED_FIELDS = {
    "z",
    "atom",
    "category",
    "refs",
    "state",
    "alpha",
    "year",
    "comments",
}
TABLE_FIELDS = {
    "Z": "z",
    "Atom": "atom",
    "Refs.": "refs",
    "State": "state",
    "Alpha": "alpha",
    "Year": "year",
    "Comments": "comments",
}


def iter_submission_files(root: Path = SUBMISSIONS_DIR):
    for category in sorted(CATEGORIES):
        yield from sorted((root / category).glob("*.json"))


def load_record(path: Path) -> dict:
    with path.open() as f:
        return json.load(f)


def validate_record(path: Path, record: dict) -> list[str]:
    errors = []
    missing = sorted(REQUIRED_FIELDS.difference(record))
    if missing:
        errors.append(f"{path}: missing required fields: {', '.join(missing)}")

    category = record.get("category")
    if category not in CATEGORIES:
        errors.append(f"{path}: category must be one of {sorted(CATEGORIES)}")
    elif path.parent.name != category:
        errors.append(f"{path}: category does not match containing directory")

    try:
        z = int(record.get("z"))
        expected_atom = Element(z).symbol
    except (TypeError, ValueError):
        errors.append(f"{path}: z must be a valid atomic number")
    else:
        if record.get("atom") != expected_atom:
            errors.append(f"{path}: atom must be {expected_atom!r} for z={z}")

    for field in ("refs", "state", "alpha", "comments"):
        value = record.get(field)
        if value is None or str(value).strip() == "":
            errors.append(f"{path}: {field} must not be empty")
        elif "|" in str(value):
            errors.append(f"{path}: {field} must not contain the pipe separator")

    for field in ("year", "release_year"):
        if field in record:
            try:
                int(record[field])
            except (TypeError, ValueError):
                errors.append(f"{path}: {field} must be an integer year")

    reference_files = record.get("reference_files", [])
    if reference_files is None:
        reference_files = []
    if not isinstance(reference_files, list):
        errors.append(f"{path}: reference_files must be a list")
    else:
        for rel_path in reference_files:
            ref_path = REFERENCES_DIR / rel_path
            if ref_path.suffix != ".bib":
                errors.append(f"{path}: reference file must end with .bib: {rel_path}")
            if not ref_path.is_file():
                errors.append(f"{path}: missing reference file: {ref_path}")

    return errors


def validate_all(root: Path = SUBMISSIONS_DIR) -> list[str]:
    errors = []
    seen_rows = {}
    for path in iter_submission_files(root):
        try:
            record = load_record(path)
        except json.JSONDecodeError as exc:
            errors.append(f"{path}: invalid JSON: {exc}")
            continue
        errors.extend(validate_record(path, record))
        row_key = tuple(str(record.get(source, "")).strip() for source in TABLE_FIELDS.values())
        if row_key in seen_rows:
            errors.append(f"{path}: duplicate accepted record; first seen in {seen_rows[row_key]}")
        else:
            seen_rows[row_key] = path
    return errors


def main(args: Optional[list] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--submissions-dir",
        type=Path,
        default=SUBMISSIONS_DIR,
        help="Directory containing theoretical/ and experimental/ accepted records.",
    )
    parsed_args = parser.parse_args(args)

    errors = validate_all(parsed_args.submissions_dir)
    if errors:
        for error in errors:
            print(error)
        return 1

    print("Accepted data records are valid.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
