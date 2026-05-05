# Design: Expand Chapter 6 — 2D Frame Element

**Date:** 2026-05-05  
**File:** `chapters/ch06_frame.tex`  
**Current state:** 27-line stub (T₆ matrix + one-sentence description)  
**Goal:** Full chapter matching the depth of Ch. 5 (beam), with lean theory and one comprehensive worked example.

---

## Decisions

| Question | Decision |
|---|---|
| Worked example type | Portal frame (two columns + beam) |
| Loading | Combined: lateral point load at top-left joint + uniform vertical load on beam |
| Theory depth | Lean — derive 6×6 local stiffness and transformation; cross-reference Ch. 3 and Ch. 5 for PVW/PMPE derivations |
| Structure | Approach A: theory sections → one comprehensive worked example |

---

## Section Outline

### §6.1 Frame Element *(replaces current stub)*

- Prose: what a frame element is (beam-column combining axial + bending stiffness), when frames are used vs. trusses (moment-carrying joints)
- TikZ figure: two-node frame element in local frame, showing all 6 DOFs (ū₁, v̄₁, θ₁, ū₂, v̄₂, θ₂) with conjugate forces (N̄₁, V̄₁, M̄₁, N̄₂, V̄₂, M̄₂)
- Follow book figure conventions: displacements in blue, forces in red, hatched supports

### §6.2 Local Stiffness Matrix

- Derive k̄^e (6×6) as block-diagonal: bar (2×2) in DOFs {1,4} + beam (4×4) in DOFs {2,3,5,6}
- Write the full 6×6 matrix explicitly with EA/L and EI/L³ terms
- `importantbox` highlighting the final result
- Remark on null space: 3 rigid-body modes (2 translations + 1 rotation)

### §6.3 Element Load Vector

- Combined 6×1 f̄^e:
  - Axial distributed load p₀ → DOFs 1, 4: (p₀L/2)[1; 1] from Ch. 3
  - Transverse distributed load q₀ → DOFs 2, 3, 5, 6: (q₀L/12)[6; L; 6; −L] from Ch. 5
- Cross-reference Ch. 3 §3.3 and Ch. 5 §5.3 for derivations

### §6.4 Coordinate Transformation *(expands current content)*

- TikZ figure: frame element inclined at angle α in global frame, showing both local (x̄, ȳ) and global (x, y) axes
- Explain T₆ as two stacked copies of the 3×3 rotation block (angle α, CCW from positive x-axis)
- Global stiffness: k^e = T₆ᵀ k̄^e T₆
- Angle convention box: α measured CCW from positive x-axis to element axis (node 1 → node 2)
- Remark: columns are typically vertical (α = 90°, so c = 0, s = 1); T₆ simplifies

### §6.5 Worked Example: Portal Frame under Combined Loading

**Geometry:**
- 4 nodes, 3 elements (left column, beam, right column)
- Column height H = 3 m; beam span B = 4 m
- Both column bases fixed (nodes 1 and 4)
- All members: E = 210,000 N/mm², same I and A (chosen for clean numbers)

**Loading:**
- Horizontal point force F = 10 kN at node 2 (top-left joint, rightward)
- Uniform vertical load q = 5 kN/m downward on beam (element 2)

**DOF map:**
- Node 1 (bottom-left, fixed): DOFs 1,2,3 — all zero
- Node 2 (top-left, free): DOFs 4,5,6
- Node 3 (top-right, free): DOFs 7,8,9
- Node 4 (bottom-right, fixed): DOFs 10,11,12 — all zero
- 12 global DOFs, 6 free after BCs

**Element data:**
- Element 1: node 1→2, left column, α = 90°, L = H
- Element 2: node 2→3, beam, α = 0°, L = B
- Element 3: node 3→4, right column, α = −90° (pointing downward), L = H

**Steps:**
1. Element data table (L, α, c, s, EA/L, EI/L³ for each element)
2. Local stiffness matrices k̄^e (write numerically for each element)
3. Transformation matrices T₆ for each element (exploit c=0,s=±1 for columns)
4. Global element stiffnesses k^e = T₆ᵀ k̄^e T₆
5. Assembly into 12×12 global K (describe DOF mapping)
6. BC reduction → 6×6 reduced system
7. Solution (MATLAB-verified)
8. Reactions at fixed bases (equilibrium check: ΣFx = 0, ΣFy = 0, ΣM = 0)
9. Element internal forces: recover local DOFs via ū^e = T₆ u^e; compute N, V, M along each element
10. Deflected shape, bending moment, shear force, axial force diagrams via MATLAB script
11. Reference to MATLAB listing in Appendix

### §6.6 Internal Force Recovery

- Half-page standalone section explaining the general procedure:
  1. Extract global nodal DOFs for each element
  2. Transform to local: ū^e = T₆ u^e
  3. Compute local forces: f̄^e = k̄^e ū^e − f̄^e_load
  4. Distribute along element using bar/beam shape functions
- Cross-reference MATLAB script
- Remark: for vertical columns, N is axial (vertical), V is shear (horizontal); orientation matters for sign convention

---

## MATLAB Script

New file: `matlab/fem_frame_portal.m`

- Header: chapter 6, portal frame example
- `clear; clc; close all;`
- Self-contained, no toolboxes
- Outputs: deflected shape, bending moment diagram, shear force diagram, axial force diagram
- Saves PDF + PNG to `../figures/` using `exportgraphics` with `'ContentType','vector'`
- Listing label: `lst:frame_portal`

---

## Notation

- Local quantities use overbar: k̄^e, f̄^e, ū^e (consistent with Ch. 4 truss)
- Angle α for frame inclination (θ already used for rotation DOF in beam chapter)
- Use `\vect`, `\matr`, `\stiff`, `\transpose`, `\dd`, `\pd` macros from preamble.tex
- No `q` for axial load (reserved); use `p₀` for axial distributed load

---

## Cross-References

- §6.2 cross-refs: `eq:bar_stiffness` (Ch. 3), `eq:beam_stiffness` (Ch. 5)
- §6.3 cross-refs: Ch. 3 §3.3 (bar load vector), Ch. 5 §5.3 (beam load vector)
- §6.4 cross-refs: `eq:truss_T` (Ch. 4 truss transformation, same 2×2 block pattern)
- §6.6 cross-refs: MATLAB listing `lst:frame_portal`

---

## Constraints

- All numerics must be MATLAB-verified before committing to the chapter
- No toolboxes in MATLAB script
- Listings in appendix: replace em-dashes with `--` inside `lstlisting` blocks
- TikZ style alias: use `zc` not `z` (collides with TikZ z-key)
