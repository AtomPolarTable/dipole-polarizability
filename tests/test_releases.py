from pathlib import Path

from pydipole.__main__ import main
from pydipole.utils import get_data_paths, load_db


def write_release(root: Path, year: int, alpha: str) -> Path:
    release_dir = root / "data" / "releases" / str(year)
    release_dir.mkdir(parents=True)
    (release_dir / "database.csv").write_text(
        "Z|Atom|Refs.|State|Alpha|Year|Comments\n"
        f"1|H|[\\citenum{{Example{year}}}]|$^2S$|${alpha}$|{year}|test release\n"
    )
    (release_dir / "references.bib").write_text(f"@article{{Example{year}, year = {{{year}}}}}\n")
    return release_dir


def test_load_db_from_release(tmp_path, monkeypatch):
    write_release(tmp_path, 2026, "4.5")
    monkeypatch.chdir(tmp_path)

    database, references = get_data_paths("2026")
    assert database.name == "database.csv"
    assert references.name == "references.bib"

    df = load_db(release="2026")
    assert df.iloc[0]["Atom"] == "H"
    assert df.iloc[0]["Year"] == 2026


def test_latest_release_uses_newest_numeric_directory(tmp_path, monkeypatch):
    write_release(tmp_path, 2026, "4.5")
    write_release(tmp_path, 2027, "4.6")
    monkeypatch.chdir(tmp_path)

    df = load_db(release="latest")
    assert df.iloc[0]["Alpha"] == "$4.6$"


def test_write_table_from_release(tmp_path, monkeypatch):
    write_release(tmp_path, 2026, "4.5")
    monkeypatch.chdir(tmp_path)
    table_path = tmp_path / "table.tex"
    bib_path = tmp_path / "references-copy.bib"

    assert main([str(table_path), "--release", "2026", "--bib", str(bib_path)]) == 0

    assert table_path.is_file()
    assert "test release" in table_path.read_text()
    assert bib_path.read_text().startswith("@article{Example2026")
