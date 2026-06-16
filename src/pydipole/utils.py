# PYDIPOLE: molecular density partition schemes based on HORTON package.
# Copyright (C) 2023-2026 The PYDIPOLE Development Team
#
# This file is part of PYDIPOLE
#
# PYDIPOLE is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# PYDIPOLE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --
import os
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from importlib_resources import files

__all__ = [
    "DATA_PATH",
    "FNAME_DB",
    "FNAME_BIB",
    "SEP",
    "find_project_data_dir",
    "get_data_paths",
    "load_db",
]

# Constants
DATA_PATH = files("pydipole.data")
FNAME_DB = DATA_PATH.joinpath("database.csv")
FNAME_BIB = DATA_PATH.joinpath("references.bib")
SEP = "|"


def find_project_data_dir(start: Optional[Path] = None) -> Optional[Path]:
    """Find the repository data directory from the current working tree."""
    env_data_dir = os.environ.get("PYDIPOLE_DATA_DIR")
    if env_data_dir:
        data_dir = Path(env_data_dir)
        if (data_dir / "releases").is_dir():
            return data_dir

    current = Path.cwd() if start is None else Path(start)
    for parent in (current, *current.parents):
        data_dir = parent / "data"
        if (data_dir / "releases").is_dir():
            return data_dir
    return None


def get_data_paths(release: Optional[Union[str, int]] = None):
    if release is None:
        return FNAME_DB, FNAME_BIB

    data_dir = find_project_data_dir()
    if data_dir is None:
        raise FileNotFoundError(
            "Could not find data/releases. Run from the repository or set PYDIPOLE_DATA_DIR."
        )

    releases_dir = data_dir / "releases"
    if str(release) == "latest":
        years = sorted(path.name for path in releases_dir.iterdir() if path.is_dir() and path.name.isdigit())
        if not years:
            raise FileNotFoundError(f"No annual releases found in {releases_dir}")
        release = years[-1]

    release_dir = releases_dir / str(release)
    database = release_dir / "database.csv"
    references = release_dir / "references.bib"
    if not database.is_file() or not references.is_file():
        raise FileNotFoundError(f"Release {release!r} is missing database.csv or references.bib")
    return database, references


def load_db(release: Optional[Union[str, int]] = None):
    database, _ = get_data_paths(release)
    return pd.read_csv(database, sep=SEP)
