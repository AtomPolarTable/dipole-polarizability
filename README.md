# Pydipole
[![Python 3.8](https://img.shields.io/badge/python-3.8-blue.svg)](https://docs.python.org/3.8/)
[![Python 3.9](https://img.shields.io/badge/python-3.9-blue.svg)](https://docs.python.org/3.9/)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://docs.python.org/3.10/)
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://docs.python.org/3.11/)
[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://docs.python.org/3.12/)

<div align="center">
  <img src="./docs/pydipole.svg"  width="150px" />
</div>

[`Pydipole`](https://github.com/AtomPolarTable/dipole-polarizability) is a Python package to maintain Table of atomic dipole polarizability for neutural atoms in the Periodic Table.

This version contains contributions from:
Peter Schwerdtfeger (1)
Jeffrey K. Nagle (2)
YingXing Cheng (3),

- (1) Centre for Theoretical Chemistry and Physics, Massey University, Auckland, New Zealand.
- (2) Department of Chemistry, Bowdoin College, Brunswick, ME, USA.
- (3) Numerical Mathematics for High-Performance Computing (NMH), University of Stuttgart, Stuttgart, Germany.

The `pydipole` source code is hosted on GitHub and is released under the GNU General Public License v3.0.
Please report any issues you encounter while using the `pydipole` library on [GitHub Issues](https://github.com/AtomPolarTable/dipole-polarizability/issues/new).
For further information and inquiries, please contact us at yxcheng2buaa@gmail.com.

## Tables

- The 2023 Table is available [here](https://github.com/AtomPolarTable/dipole-polarizability/blob/main/tables/2023/main.pdf).

## License

![GPLv3 License](https://img.shields.io/badge/license-GPLv3-blue.svg)


`pydipole` is distributed under GPL License version 3 (GPLv3).

## Dependencies

The following dependencies will be necessary for `pydipole` to build properly,

* pandas : https://github.com/pandas-dev/pandas


## Installation

To install `pydipole` with version `0.0.x`:

```bash
pip install pydipole
```

To install latest `pydipole`:

```bash
git clone http://github.com/AtomPolarTable/dipole-polarizability
cd dipole-polarizability
pip install .
```

To run test, one needs to add tests dependencies for `tests`:

```bash
pip install .[tests]
```

For developers, one could need all dependencies:
```bash
pip install -e .[dev,tests]
```

## Citations

Please use the following citations in any publication using `pydipole` library:

[1] Schwerdtfeger, P.; Nagle, J. K.
**2018 Table of static dipole polarizabilities of the neutral elements in the periodic table.**
Molecular Physics 2019, 117 (9-12), 1200–1225.
[![DOI: 10.1080/00268976.2018.1535143](https://img.shields.io/badge/DOI-10.1080/00268976.2018.1535143-blue)](https://doi.org/10.1080/00268976.2018.1535143)
