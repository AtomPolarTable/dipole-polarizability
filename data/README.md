# Curated Data Workflow

Curated records are stored as one JSON file per record under
`data/submissions/theoretical/` or `data/submissions/experimental/`. GitHub
issues are the intake channel; files in this directory are the reviewed,
maintainer-approved source of truth for table data.

Rows imported from the historical package database include
`source: "legacy_database"` and `source_row`. To correct old data, edit the
matching JSON file, validate, then rebuild the release snapshot.
Duplicate historical rows were removed during import, and validation rejects
duplicate accepted records.

Do not edit `database.csv` by hand for routine data updates. Release CSV files
are generated snapshots, and `src/pydipole/data/database.csv` is only refreshed
when maintainers intentionally update the package default data.

## Accepted Record Schema

Each JSON file must contain:

```json
{
  "z": 18,
  "atom": "Ar",
  "category": "theoretical",
  "refs": "[\\citenum{Example2026}]",
  "state": "$^1S_0$",
  "alpha": "$11.083$",
  "year": 2026,
  "comments": "CCSD(T), static dipole polarizability.",
  "reference_files": ["example-2026.bib"]
}
```

- `category` must match the containing folder: `theoretical` or `experimental`.
- `z` and `atom` must identify the same element.
- `alpha` should be in atomic units unless the comments explicitly document a
  conversion.
- `reference_files` are BibTeX files stored in `data/references/`.

## Annual Releases

Maintainers publish yearly snapshots with:

```bash
python scripts/build_release.py --year 2026
```

The command validates records, starts from the latest earlier annual release
when one exists, and adds JSON records whose `release_year` belongs in the
target snapshot. If no earlier release exists, the snapshot is built directly
from JSON records.

Use `--sync-package-data` only when preparing a package release that should
update `src/pydipole/data/database.csv` and
`src/pydipole/data/references.bib`.
