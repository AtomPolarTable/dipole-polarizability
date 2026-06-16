Quick Start
###########

Load the Curated Database
=========================

Install ``pydipole`` and load the bundled table as a pandas DataFrame:

.. code-block:: python

   from pydipole import Element, load_db

   df = load_db()
   print(df.head(20))

Filter rows with ordinary pandas expressions:

.. code-block:: python

   hydrogen = df[df["Atom"] == "H"]
   element_119 = df[df["Atom"] == Element(119).symbol]

Use Annual Release Data
=======================

When working from a repository checkout, generate tables from a curated annual
release instead of the packaged default:

.. code-block:: bash

   write-table --release 2026 table.tex --bib references.bib
   write-table --latest table.tex --bib references.bib

The ``--latest`` option selects the newest numeric directory under
``data/releases/``.

Generate a LaTeX Table
======================

The default command uses the bundled package data:

.. code-block:: bash

   write-table table.tex --bib references.bib

Compile the generated ``table.tex`` with LaTeX. A complete example is available
in ``tables/2023/``.

Table Notes
===========

Static scalar dipole polarizabilities are listed in atomic units for neutral
atoms. Unless otherwise indicated by the state symmetry, values are
``M_L``/``M_J`` averaged. Common abbreviations include ``exp.`` for
experimental values, ``NR`` for nonrelativistic calculations, ``R`` for
relativistic calculations, ``MBPT`` for many-body perturbation theory, ``CI``
for configuration interaction, and ``CCSD(T)`` for coupled cluster singles and
doubles with perturbative triples.
