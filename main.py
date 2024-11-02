#!/usr/bin/env python

import os

import pandas as pd


def to_tex(df, fn_out):
    df.rename(columns={"Alpha": r"$\alpha$", "Refs": "Refs."}, inplace=True)

    # caption = r"""
    # Static scalar dipole polarizabilities (in atomic units) for neutral atoms.
    # If not otherwise indicated by the state symmetry, $M_L (M_J)$-averaged polarizabilities are
    # listed; $M_L (M_J)$ respectively denotes that the polarizability for each $M_L (M_J)$ state
    # can be found in the reference given. Abbreviations used (uncertainties given here consistently
    # as $\pm$ values): exp.: experimentally determined value; NR: nonrelativistic; R: Relativistic,
    # DK: Scalar relativistic Douglas-Kroll; MVD: mass-velocity-Darwin; SO: Spin-orbit coupled;
    # SF: Dyall's spin-free formalism (scalar relativistic); PP: relativistic pseudopotential;
    # LDA: local (spin) density approximation; PW91: Perdew-Wang 91 functional;
    # RPA: Random phase approximation; PolPot: Polarization potential; MBPT: many-body
    # perturbation theory; CI: configuration interaction; CCSD(T): coupled cluster singles
    # doubles (SD) with perturbative triples; FS Fock-space; CEPA: coupled electron pair
    # approximation; MR: multi-reference; CAS: complete active space; VPA: variational perturbation
    # approach. For all other abbreviations see text or references. If the symmetry of the state
    # is not clearly specified as in Doolen's calculations, the calculation was most likely set
    # at a specific configuration (orbital occupancy) as listed in the Desclaux tables
    # \citenum{Müller1984}, reflecting the ground state symmetry of the specific atom.
    # NB: 1 a.u. $= 0.1481847113 \, \text{\AA}^3 = 1.6487773 \times 10^{-41} \,
    # \text{C m}^2/\text{V}$.
    # """

    label = "tab:table_2023"

    # Generate LaTeX code
    latex_code = df.to_latex(
        longtable=True,
        escape=False,
        na_rep="$--$",
        index=False,
        label=label,
    )
    # # Insert the caption before the label
    # caption_insert_position = latex_code.find("\\label{")
    # latex_code_with_caption = (
    #     latex_code[:caption_insert_position]
    #     + f"\\caption{{{caption}}}\n"
    #     + latex_code[caption_insert_position:]
    # )

    os.makedirs(os.path.dirname(fn_out), exist_ok=True)

    # Save to a file
    with open(fn_out, "w") as f:
        # f.write(latex_code_with_caption)
        f.write(latex_code)


def main():
    dir_all = "pdf"
    os.makedirs(dir_all, exist_ok=True)

    fn_in = "data/database.csv"
    fn_out = f"{dir_all}/table_2023.tex"

    df = pd.read_csv(fn_in, sep=";")
    df["Atom"] = df["Atom"].mask(df["Atom"].duplicated(), "")
    df["Z"] = df["Z"].mask(df["Z"].duplicated(), "")
    to_tex(df, fn_out)


if __name__ == "__main__":
    main()
