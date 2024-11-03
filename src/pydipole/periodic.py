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
import numpy as np

__all__ = ["Element"]
# fmt: off
_PERIODIC_TABLE = [
    'dummy', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S',
    'Cl', 'Ar','K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge',
    'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
    'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm',
    'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',
    'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'
]
# fmt: on

_SPECIAL_ELEMENTS = {
    24: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^1 3d^5",
    29: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^1 3d^10",
    41: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^4",
    42: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^5",
    44: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^7",
    45: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^8",
    46: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 4s^10",
    47: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^10",
    57: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^10 5p^6 6s^2 5d^1",
    58: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^10 5p^6 6s^2 4f^1 5d^1",
    64: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^10 5p^6 6s^2 4f^7 5d^1",
    78: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^1 4f^14 5d^9",
    79: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^1 4f^14 5d^10",
    89: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^2 4f^14 5d^10 6p^6 7s^2 6d^1",
    90: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^2 4f^14 5d^10 6p^6 7s^2 6d^2",
    91: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^2 4f^14 5d^10 6p^6 7s^2 5f^2 6d^1",
    92: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^2 4f^14 5d^10 6p^6 7s^2 5f^3 6d^1",
    93: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^2 4f^14 5d^10 6p^6 7s^2 5f^4 6d^1",
    96: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^2 4f^14 5d^10 6p^6 7s^2 5f^7 6d^1",
    103: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^2 4f^14 5d^10 6p^6 7s^2 5f^14 7p^1",
    110: "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^2 4d^10 5p^6 6s^2 4f^14 5d^10 6p^6 7s^1 5f^14 6d^9",
}


class Element:
    """
    Represents a chemical element in the periodic table.

    Attributes
    ----------
    atomic_number : int
        The atomic number of the element.
    symbol : str
        The symbol of the element.
    """

    ROW_SIZES = (2, 8, 8, 18, 18, 32, 32)

    def __init__(self, element):
        """
        Initialize the Element object.

        Parameters
        ----------
        element : str or int
            The symbol or atomic number of the element.
        """
        if isinstance(element, (str, int, np.integer)):
            idx, symbol = self._get_element(element)
        else:
            raise ValueError(f"Expected a <str> or <int>, got: {type(element):s}")

        self.atomic_number = idx
        self.symbol = symbol

    @staticmethod
    def _get_element(element):
        """
        Get the element by the given id (either symbol or atomic number).

        Parameters
        ----------
        element : str or int
            The symbol or atomic number of the element.

        Returns
        -------
        tuple
            A tuple containing the atomic number and symbol of the element.
        """

        def get_symbol(Z):
            """
            Get element by given index.
            """
            if isinstance(Z, (int, np.integer)):
                if Z < 0 or Z > 118:
                    raise ValueError(f"Wrong atomic number: {Z}")
                return _PERIODIC_TABLE[Z]
            else:
                raise RuntimeError("Type id should be int")

        def get_id(symbol):
            """
            Get element number by given element type.
            """

            if symbol not in _PERIODIC_TABLE:
                raise ValueError(f"Wrong atomic symbol {symbol}")
            idx = _PERIODIC_TABLE.index(symbol)
            if idx >= 0:
                return idx
            else:
                raise RuntimeError("Wrong element type!")

        if isinstance(element, str):
            return get_id(element), element
        else:
            return element, get_symbol(element)

    @property
    def block(self):
        """
        Get the block character (s, p, d, or f) of the element.

        Returns
        -------
        str
            The block character of the element.
        """
        if (self.is_actinoid or self.is_lanthanoid) and self.atomic_number not in [
            71,
            103,
        ]:
            return "f"
        if self.is_actinoid or self.is_lanthanoid:
            return "d"
        if self.group in [1, 2]:
            return "s"
        if self.group in range(13, 19):
            return "p"
        if self.group in range(3, 13):
            return "d"
        raise ValueError("unable to determine block")

    @property
    def is_lanthanoid(self):
        """
        Check if the element is a lanthanoid.

        Returns
        -------
        bool
            True if the element is a lanthanoid, False otherwise.
        """
        return 56 < self.atomic_number < 72

    @property
    def is_actinoid(self):
        """
        Check if the element is an actinoid.

        Returns
        -------
        bool
            True if the element is an actinoid, False otherwise.
        """
        return 88 < self.atomic_number < 104

    @property
    def row(self):
        """
        Get the periodic table row of the element.

        Returns
        -------
        int
            The row of the element in the periodic table.
        """
        z = self.Z
        total = 0
        if 57 <= z <= 71:
            return 8
        if 89 <= z <= 103:
            return 9
        for i, size in enumerate(Element.ROW_SIZES):
            total += size
            if total >= z:
                return i + 1
        return 8

    @property
    def group(self):
        """
        Get the periodic table group of the element.

        Returns
        -------
        int
            The group of the element in the periodic table.
        """
        z = self.atomic_number
        if (self.is_actinoid or self.is_lanthanoid) and z not in [71, 103]:
            return None

        if z == 1:
            return 1
        if z == 2:
            return 18
        if 3 <= z <= 18:
            if (z - 2) % 8 == 0:
                return 18
            if (z - 2) % 8 <= 2:
                return (z - 2) % 8
            return 10 + (z - 2) % 8

        if 19 <= z <= 54:
            if (z - 18) % 18 == 0:
                return 18
            return (z - 18) % 18

        assert z >= 54
        rest = (z - 54) % 32
        # group 18
        if rest == 0:
            return 18
        # group 1 and 2
        elif 1 <= rest <= 2:
            return rest
        #  group 3-12
        elif 2 < rest < 16:
            return None
        # group 13-17
        elif rest >= 16:
            return rest - 14
        else:
            raise ValueError("You found a bug!")

    @property
    def group_symbol(self):
        """
        Get the group symbol of the element.

        Returns
        -------
        str
            The group symbol of the element.
        """
        return f"group-{self.group}"

    @property
    def Z(self):
        """
        Get the atomic number of the element.

        Returns
        -------
        int
            The atomic number of the element.
        """
        return self.atomic_number

    @property
    def period(self):
        """
        Get the period of the element in the periodic table.

        Returns
        -------
        int
            The period of the element in the periodic table.
        """
        return self.row

    def get_elec_config(self):
        """
        Get the electron configuration of the element.

        Returns
        -------
        str
            The electron configuration of the element.
        """
        n = self.Z
        if n in _SPECIAL_ELEMENTS.keys():
            return _SPECIAL_ELEMENTS[n]

        rule = "1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s 5f 6d 7p"
        nb_dict = {"s": 2, "p": 6, "d": 10, "f": 14}
        orbitals = [(i, nb_dict[i[-1]]) for i in rule.split()]
        output = []

        for orbital, size in orbitals:
            k = min(size, n)
            output.append("%s^%d" % (orbital, k))
            n -= k
            if n <= 0:
                break

        orbital_info = " ".join(output)
        return orbital_info
