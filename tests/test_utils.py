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

from pydipole.periodic import Element
from pydipole.utils import load_db


def test_load_db():
    df = load_db()
    assert len(df.columns) == 7
    assert df[df["Atom"] == "He"].size > 0
    assert df[df["Atom"] == Element(119).symbol].size > 0
    # print(df[df["Atom"] == Element(120).symbol])
    assert df[df["Atom"] == Element(120).symbol].size > 0
