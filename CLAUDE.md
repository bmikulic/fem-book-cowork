# FEM Book — Project Notes

## Goal

Boris is writing **Introduction to the Finite Element Method — Theory, Implementation, and MATLAB Examples**, a senior-undergraduate / early-graduate textbook for structural and mechanical engineers. The pedagogy emphasises connecting theory to a runnable MATLAB example in every chapter.

The author is a Croatian engineering professor; uploaded reference materials may be in Croatian, but the book itself is in English.

## Build

LaTeX project rooted at `main.tex`. Build with:

```
latexmk -pdf main.tex
```

Bibliography is plain BibTeX (`bibliography.bib`, `\bibliographystyle{plain}`). The preamble lives in `preamble.tex`.

## File Layout

```
main.tex                  # book root: \include's the chapters and appendix
preamble.tex              # all packages, custom commands, theorem envs, listings style
bibliography.bib          # BibTeX entries
chapters/
  ch01_intro.tex          # Ch.1 Introduction
  ch02_weighted_residuals.tex   # Ch.2 Weighted residual methods (collocation, subdomain, Galerkin)
  ch02_1d_elements.tex    # Ch.3 1D bar/rod elements   (filename starts ch02_ — see Quirks)
  ch04_truss.tex          # Ch.4 2D truss element
  ch05_beam.tex           # Ch.5 Euler-Bernoulli beam and frame elements
  ch04_2d_continuum.tex   # Ch.5 2D continuum (CST, Q4 isoparametric)
  ch05_dynamics.tex       # Ch.6 Structural dynamics, Newmark-β
appendix/
  matlab_listings.tex     # All MATLAB code listings, lst:bar1d / lst:stepped_bar / lst:pvw_verification / lst:beam2d / lst:cst2d / lst:dynamics
figures/                  # PDFs (used by LaTeX) and PNGs (for previewing) — produced by scripts in matlab/
matlab/                   # Self-contained .m scripts; each saves PDF+PNG into ../figures/
```

## Notation Conventions (from preamble.tex)

Always use these; the printed book is consistent on them.

| Macro | Meaning |
|---|---|
| `\vect{u}` | column vector (bold italic) |
| `\matr{K}` | matrix (bold italic) |
| `\stiff` | global stiffness `\matr{K}` |
| `\mass`, `\damp` | global mass / damping |
| `\transpose` | `^{\mathsf{T}}` |
| `\inv` | `^{-1}` |
| `\dd` | upright differential `\,\mathrm{d}` |
| `\pd{f}{x}` | partial derivative |

The `importantbox` tcolorbox environment (blue, "Key Result") is used to highlight defining equations like the bar stiffness matrix — keep this for major results only.

## Symbol Choices

- **`p(x)` and `p_0`** for distributed *axial* load throughout (not `q`). Earlier drafts used `q`; switched to `p` per author's preference. Don't reintroduce `q` for axial loads.
- **`q`** is reserved for the *transverse* uniform beam load in Ch.4 (truss/beam) — that's the standard convention there, leave it alone.
- **`r, s`** are local DOF indices in the assembly notation `K_{ij} += k^e_{rs}` (Ch.3). Avoid `p, q` for indices since `p` is the load.

## Pedagogical Order in Ch.3 (Bar Elements)

The chapter deliberately presents the bar element in this order — please respect it:

1. **Direct Formulation** — derive `k^e = (EA/L)[[1,−1],[−1,1]]` from physical equilibrium + Hooke's law, no integrals.
2. **Displacement Field Approximation** — introduce shape functions `N_1, N_2`, the B-matrix, then *enrich the trial* with a particular solution: `u(x) = N_1 u_1 + N_2 u_2 + u_p(x)` where `u_p` solves `−EA u_p'' = p(x)` with homogeneous endpoint BCs. This is the chapter's signature.
3. **Element Load Vector** — derive `f^e = (p_0 L / 2)[1; 1]` from the enriched trial, note it agrees with the variational `∫ N^T p dx`.
4. **Global Assembly**, **Boundary Conditions**.
5. **Worked Example: Three-Element Bar** — simplest illustration.
6. **Worked Example: Stepped Bar with Combined Loading** — full workflow including the *graphical* `K_global` assembly figure (color-coded 3×3 matrices with overlap shading); MATLAB script `matlab/stepped_bar_2elem.m`.
7. **Principle of Virtual Work** — placed *after* the worked examples on purpose. Statement → derivation from equilibrium → recovery of `k^e` from PVW → verification example using the already-solved stepped bar.
8. **Natural Coordinates** — closes the chapter, bridges to Ch.5.

The PVW worked example is a *verification* (compute `δW_int` and `δW_ext` separately for an admissible `δu`, show they're equal numerically). Don't replace it with a "use PVW to solve a new problem" example — that role is served by Ch.2's weighted-residual treatment.

## MATLAB Conventions

- One `.m` file per worked example, named to match the example: `weighted_residual_comparison.m`, `stepped_bar_2elem.m`, `pvw_linear_load.m`, `pvw_verification_stepped_bar.m`.
- Top of each script: `clear; clc; close all;` and a header comment naming the chapter and example.
- Each script saves both `.pdf` (vector, for LaTeX `\includegraphics`) and `.png` (300 dpi, for previewing) into `../figures/`. Use `exportgraphics` with `'ContentType','vector'` for PDFs.
- All scripts are self-contained and use only core MATLAB — no toolboxes required. Keep that constraint.
- Plot labels in MATLAB code should match the math symbols used in the book (e.g. `uEA/(p_0 L^2)`, not `uEA/(q_0 L^2)`).
- When listings are embedded in `appendix/matlab_listings.tex`, *replace em-dashes (—) with `--` inside the listing block* — `lstlisting` chokes on UTF-8.

## Figure Conventions

- Geometry/free-body TikZ figures use this color code:
  - **Blue** (`blue!70!black`) for *displacements* (kinematics)
  - **Red** (`red!75!black`, ultra thick) for *forces* (statics)
  - Distributed loads in lighter blue with multiple short arrows
  - Hatched walls (`pattern=north east lines`) for fixed supports
- TikZ style aliases: don't name a style just `z` — it collides with TikZ's built-in z-key. Use `zc` or similar.
- The graphical `K_global` assembly diagram (Ch.3 stepped bar) uses three 3×3 colored matrices side-by-side with `+` and `=`: blue for element 1's contribution, red for element 2's, violet for the overlap cell.

## Chapter Numbering Quirk

The filenames `chapters/ch02_*` don't match the *printed* chapter numbers (Ch.3 in the printed book is `ch02_1d_elements.tex`, etc.). This came from an early reordering (Ch.2 weighted residuals was inserted later). Don't rename the files — the cross-references work via `\label` keys (`ch:1d`, `ch:weighted_residuals`, `ch:truss`, `ch:beam`, `ch:2d_continuum`, `ch:dynamics`), and renaming would propagate through the includes for no reader-visible benefit.

## Known Editorial TODOs

- The `fem_bar_1d.m` script in `appendix/matlab_listings.tex` (lst:bar1d, three-element bar) has a suspicious stress-print line near the end (`sigma * Le / (EA/Le / (EA/L))`) — likely a leftover from an earlier draft. Worth simplifying when you next touch that listing.
- `\listoflistings` is not enabled. If a list of MATLAB listings is wanted in the front matter, that's a small addition (load `\usepackage[chapter]{minted}` alternative or use `lstlistoflistings`).
- The fancyhdr "headheight too small" warning in the build log is cosmetic — fixable with `\setlength{\headheight}{14.5pt}` in `preamble.tex` if it bothers you.

## Author Preferences (from prior conversations)

- Replace existing constructs rather than accumulating versions when the new is clearly better (e.g. PVW placement, symbol changes).
- When showing assembly or matrix-building, prefer **graphical** illustrations (color-coded matrices, side-by-side with `+`/`=`) over plain `bmatrix` walls.
- Worked examples should be *complete*: discretization tables, every element matrix written numerically, assembly visual, BC step, solution, reactions, postprocessing fields, and a MATLAB listing in the appendix that reproduces all of it.
- Sanity-check numerics before presenting: every claim like "u_2 = −0.1058 mm" should be confirmed by a MATLAB run, not just hand algebra.

## When in Doubt

Read the relevant chapter `.tex` file end-to-end before editing — the prose, equations, figures, and MATLAB are tightly cross-referenced (`\cref`, equation labels, listing labels). Breaking one usually breaks several.
