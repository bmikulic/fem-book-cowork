---
name: fem-tikz-figures
description: Use when authoring or fixing TikZ figures, MATLAB result plots, or LaTeX figure environments in this FEM textbook repo. Triggers include drawing element diagrams, fixing label/support/arrow collisions, adding K-assembly visuals, regenerating MATLAB plots, or styling figure floats.
---

# FEM TikZ Figures

## Overview

In FEM diagrams, geometry, supports, DOF arrows, force arrows, dimensions, and node labels all compete for the crowded margin around each node. Collisions are the *default* outcome unless every label position is verified against neighboring symbols **before** rendering.

**Core principle:** render-then-fix is more expensive than think-then-render. Bounding-box reasoning before `pdflatex` saves a round trip — and avoids the locked-PDF dance with the author's viewer.

## When to Use

- Authoring or editing any TikZ environment in `chapters/*.tex`
- Writing or modifying any MATLAB script in `matlab/` that produces a figure
- Triaging a figure the author flagged (overlap, wrong-direction arrow, wrong color, drifted float)
- Adding a new K-assembly or graphical-equation visualization

Do NOT use for: prose edits, equation typesetting, BibTeX work.

## Pre-Render Checklist

Verify all of these mentally (or by sketching coords) **before** running `pdflatex`:

1. **Sign convention.** Every DOF arrow points along `+x` (or `+y`, etc.). Don't flip signs to make arrows point "outward" — the convention is assumed-positive directions, regardless of which side of the geometry the node sits on.
2. **Label vs. support.** For each node with a support, the label sits in the free quadrant *opposite* the support. Pin below → label above-left or above-right; roller on left → label on right; roller above → label below.
3. **Label vs. arrows.** Same check, opposite direction from any DOF/force arrow extending from the node.
4. **Pin "ground line."** The horizontal line below a pin's triangle goes from `(-w, y)` to `(+w, y)` — both endpoints at the *same* `y`. A slanted line is a typo.
5. **Dimension clearance.** Each dim line has ≥0.05 cm clear of every label edge and support edge along its full extent. If a label widens, push the dim line out — don't accept the collision.
6. **Float placement.** Results figures inside `\begin{example}…\end{example}` use `[H]` (from package `float`), never `[htbp]`. Anything else lets LaTeX drift the figure past `\end{example}`.
7. **Color convention.** Blue (`blue!70!black`, thick) for displacements/kinematics; red (`red!75!black`, ultra thick) for forces/statics; lighter blue with multiple short arrows for distributed loads; hatched walls (`pattern=north east lines`) for fixed supports. No new accent colors — they're load-bearing.
8. **MATLAB typography.** Every plot script calls `book_style(fig)` immediately before the first `exportgraphics`. Default 10 pt is too small for textbook print.

## Conventions (from this codebase)

- **Style aliases.** `disp/.style={->,thick,blue!70!black}`, `force/.style={->,ultra thick,red!75!black}`. Never name a TikZ style just `z` — collides with the built-in `z` key. Use `zc` or similar.
- **Symbols.** `p` and `p_0` for axial loads (don't reintroduce `q` for axial); `q` reserved for transverse beam load; `r,s` for local DOF indices in assembly.
- **Deformed-shape plots.** `axis equal off`, no grid, modest scale (×100–×500 typical). ×800+ exaggerates so much the deformed shape reads as a second mesh.
- **Element badges.** `\textcircled{}` has tight margins. When the numeral must breathe, use `rectangle 'Position',[c-r,c-r,2r,2r] 'Curvature',[1 1]` plus a centered `text`. Radius around 0.04–0.05 of the geometry unit.
- **K-assembly visualization.** Three same-size matrix grids side by side with `+` and `=`; blue for element 1's contribution, red for element 2's, violet for the overlap cell. `\footnotesize` cell font for six-digit entries; schematic-only (no numbers) when entries grow past five digits.

## Recurring Failure Modes

| Symptom | Root cause | Fix |
|---|---|---|
| DOF arrow flipped | Author drew "outward" by negating | Always `+x`/`+y`; place the *label* opposite, not the arrow |
| Node label overlaps support | Label and support in same quadrant | Move label to opposite quadrant; if all quadrants are crowded, move the dim line out and add label padding |
| Pin looks tilted | "Ground line" drawn diagonally | Both endpoints same `y` |
| Figure drifts past `\end{example}` | `[htbp]` is too permissive | `[H]` |
| Page-width overflow on a TikZ plot | Raw model coords (mm) with `x=0.52mm` → meters wide | Pick `x=1cm` plus a sensible per-unit scale |
| Slope diagram jumps at element boundaries | Hermite derivative algebra wrong | Use clean polynomial form: `H1p = -6ξ/L² + 6ξ²/L³`, etc. Don't expand by hand |
| Orange/blue or other custom accents on a beam fig | Author drift from convention | Strip custom colors; only black geometry + blue/red conventions |
| K-assembly cells too cramped | `\scriptsize` (~7 pt) for six-digit numbers | Bump to `\footnotesize` |
| `! I can't write on file 'main.pdf'` | Author's PDF viewer holds the file | Ask author to close, then retry — do NOT delete the PDF |
| `latexmk: All targets up-to-date` after a real edit | Stale aux state | Run `pdflatex` directly, or `latexmk -C` then rebuild |

## Verification Workflow

1. Edit the `.tex` or `.m` file.
2. If MATLAB: regenerate via the matlab MCP `run_matlab_file` (don't rely on cached PDFs).
3. Build with `pdflatex -interaction=nonstopmode main.tex` — direct invocation avoids latexmk's "nothing to do" trap when only included files changed.
4. If lock error: stop, ask author to close the PDF viewer, retry. Never `rm` the PDF.
5. Open the affected page; visually verify each checklist item against the render.
6. Commit with a focused message naming the figure and the specific change. Don't bundle multiple unrelated figure fixes in one commit.

## Things to Resist

- Fixing one collision by creating another (label moved into the dim line's path).
- Reintroducing `q` as the axial-load symbol.
- Adding accent colors "for clarity."
- Claiming success on a figure change without rebuilding the PDF and visually checking the rendered page.
- Renaming the `chapters/ch02_*` files to match printed chapter numbers — see CLAUDE.md.

## Cross-References

- Static conventions and chapter-level notation: `CLAUDE.md` at repo root.
- MATLAB typography helper: `matlab/book_style.m`.
- Existing K-assembly examples: `fig:stepped_bar_K_assembly`, `fig:beam_simple_K_assembly`, `fig:truss_K_assembly`.
