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
import argparse
import logging
import shutil
import sys

from pydipole import FNAME_BIB, load_db

__all__ = ["main"]

logger = logging.getLogger(__name__)


def to_tex(df, fn_out, **kwargs):
    df.rename(columns={"Alpha": r"$\alpha$", "Refs": "Refs."}, inplace=True)
    label = "tab:table_2023"
    latex_code = df.to_latex(
        longtable=True, escape=False, na_rep="$--$", index=False, label=label, **kwargs
    )
    with open(fn_out, "w") as f:
        f.write(latex_code)


def main(args=None) -> int:
    """
    Main program for the pydipole application.

    This program accepts a configuration file as input and executes
    the 'part-gen' and 'part-dens' commands with the given configuration.

    Parameters
    ----------
    args : list
        List of arguments passed from command line.

    Returns
    -------
    int :
        Exit status (0 for success, non-zero for errors).
    """

    # Program description
    description = "pydipole main program."

    # Argument parsing
    parser = argparse.ArgumentParser(prog="write-table", description=description)
    parser.add_argument("filename", type=str, help="Output tex filename.")
    parser.add_argument(
        "-b",
        "--bib",
        type=str,
        default="references.bib",
        help="The filename for bibtex.",
    )
    parsed_args = parser.parse_args(args)

    df = load_db()
    df["Atom"] = df["Atom"].mask(df["Atom"].duplicated(), "")
    df["Z"] = df["Z"].mask(df["Z"].duplicated(), "")
    to_tex(df, parsed_args.filename)
    shutil.copy(FNAME_BIB, parsed_args.bib)
    return 0


if __name__ == "__main__":
    sys.exit(main())
