# Repository Guidelines

## Project Structure & Module Organization

This repository contains the `pydipole` Python package for atomic dipole polarizability tables. Core code lives in `src/pydipole/`, with the CLI in `src/pydipole/__main__.py`. Bundled package data is in `src/pydipole/data/`. Curated records live in `data/submissions/`, BibTeX files in `data/references/`, and generated annual snapshots in `data/releases/`. Tests are in `tests/`; docs and assets are under `docs/`; table artifacts are in `tables/2023/`.

## Build, Test, and Development Commands

- `pip install -e .[dev,tests]`: install the package in editable mode with pre-commit and pytest extras.
- `pytest`: run the test suite configured in `pyproject.toml`; tests are discovered from `tests/`.
- `pre-commit run --all-files`: run formatting, linting, and basic file checks.
- `write-table output.tex references.bib`: run the installed CLI and generate a LaTeX table plus references file.
- `write-table --release 2026 output.tex --bib references.bib`: generate a table from an annual release in `data/releases/`.
- `python scripts/validate_data.py`: validate accepted JSON submissions and referenced BibTeX files.
- `python scripts/build_release.py --year 2026`: build generated release files for a yearly snapshot.
- `cd docs && ./gen_api.sh && make html`: build the Sphinx documentation locally.
- `python -m build`: build release artifacts after installing `build`.

## Coding Style & Naming Conventions

Use Python 3.8+ compatible code. Format Python with Black at 100 columns and isort's Black profile. Ruff enforces `E` and `F`, with `E501` and `E741` ignored. Use 4-space indentation, lowercase module names, `snake_case` functions and variables, and clear test names such as `test_get_element_data`.

## Testing Guidelines

Tests use pytest. Add tests in `tests/` for any behavior change in `src/pydipole/`, especially parsing, table generation, and periodic-table lookup logic. Name new files `test_<module>.py` and new tests `test_<behavior>`. Run `pytest` before submitting changes; use focused runs such as `pytest tests/test_database.py` while developing.

## Commit & Pull Request Guidelines

The history uses short, imperative commit messages such as `fix bug`, `update docs`, and `Prepare for release 0.0.1`. Keep messages concise and action-oriented. Pull requests should describe the change, list test results, link issues, and note regenerated docs or LaTeX/PDF table artifacts.

## Security & Configuration Tips

Do not commit virtual environments, build outputs, credentials, or private data. Treat JSON files in `data/submissions/` as the editable source for accepted updates. `database.csv` files are generated snapshots; update bundled package data only through `scripts/build_release.py --sync-package-data`.
