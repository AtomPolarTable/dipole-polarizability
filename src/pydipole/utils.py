# PYDIPOLE: molecular density partition schemes based on HORTON package.
# Copyright (C) 2023-2024 The PYDIPOLE Development Team
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
import pandas as pd
from importlib_resources import files

__all__ = ["DATA_PATH", "FNAME_DB", "FNAME_BIB", "SEP", "load_db"]

# Constants
DATA_PATH = files("pydipole.data")
FNAME_DB = DATA_PATH.joinpath("database.csv")
FNAME_BIB = DATA_PATH.joinpath("references.bib")
SEP = "|"


def load_db():
    return pd.read_csv(FNAME_DB, sep=SEP)
