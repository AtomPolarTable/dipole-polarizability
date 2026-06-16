Curated Data Workflow
#####################

Source of Truth
===============

Curated records are stored as one JSON file per table row in
``data/submissions/theoretical/`` or ``data/submissions/experimental/``. The
historical packaged table was imported into these files with
``source: "legacy_database"`` and ``source_row`` metadata.

Do not edit ``database.csv`` directly for routine updates. CSV files under
``data/releases/YEAR/`` are generated snapshots, and
``src/pydipole/data/database.csv`` is refreshed only when preparing a package
release.

Accept a Submitted Record
=========================

New data should arrive through the GitHub issue templates. After review, add an
accepted record as JSON:

.. code-block:: json

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

Place referenced BibTeX files in ``data/references/``. Duplicate accepted
records are rejected by validation.

Validate and Publish
====================

Validate accepted records before building a release:

.. code-block:: bash

   python scripts/validate_data.py

Publish an annual snapshot:

.. code-block:: bash

   python scripts/build_release.py --year 2026

When preparing a new package release, also refresh the bundled default data:

.. code-block:: bash

   python scripts/build_release.py --year 2026 --sync-package-data
