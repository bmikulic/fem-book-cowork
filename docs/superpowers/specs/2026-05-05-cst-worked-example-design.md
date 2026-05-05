# Design: CST Full Worked Example (§7.2.4)

**Date:** 2026-05-05
**File:** `chapters/ch04_2d_continuum.tex`
**Goal:** Replace the 3-line §7.4 stub with a complete two-element CST worked example inserted as §7.2.4 inside the CST section.

---

## Decisions

| Question | Decision |
|---|---|
| Mesh | Two CST triangles forming a 1m × 1m square, diagonal (0,0)→(1,1) |
| Loading | Uniaxial tension in x: σx = 100 MPa, applied as nodal forces at nodes 2 and 3 |
| BCs | Pin at node 1 (u1=v1=0), horizontal roller at node 4 (u4=0) |
| Material | E = 200 GPa, ν = 0.3, plane stress, t = 0.01 m |
| Structure | §7.2.4 subsection within §7.2; old §7.4 removed |
| Approach | Full step-by-step (Approach A): every matrix written numerically |

---

## Problem Definition

- **Plate:** 1 m × 1 m × 0.01 m, plane stress
- **Material:** E = 200×10⁹ Pa, ν = 0.3
- **Nodes:** 1(0,0), 2(1,0), 3(1,1), 4(0,1) — coordinates in metres
- **Elements:** T1 = nodes {1,2,3} (lower-right), T2 = nodes {1,3,4} (upper-left)
- **Global DOFs:** [u1,v1,u2,v2,u3,v3,u4,v4] = DOFs [1…8]
- **BCs:** u1=v1=0 (pin, node 1); u4=0 (roller, node 4) → 3 fixed DOFs, 5 free
- **Load:** σx = 100 MPa → F = σx × t × height = 100×10⁶ × 0.01 × 1 = 1×10⁶ N → F₂x = F₃x = 500,000 N (DOFs 3 and 5)

---

## Section Structure

Remove: `\section{Worked Example: CST Plate in Tension}` (current §7.4, lines 185–200)

Add inside `\section{Constant Strain Triangle (CST)}` (§7.2):

```
\subsection{Worked Example: Two-Element CST Mesh in Uniaxial Tension}
\begin{example}[Two-element CST plate under uniaxial tension]
\label{ex:cst_tension}
...steps 1–10...
\end{example}
```

The label `ex:cst_tension` reuses the old label (no cross-ref chain to break).

---

## Steps Inside the Example

### Step 1 — Problem description and TikZ figure

TikZ figure (`fig:cst_example`) showing:
- 1×1 square, diagonal from (0,0) to (1,1)
- Node labels 1–4 at correct corners
- Element labels ①, ② in triangle centroids
- Pin support at node 1 (triangle symbol + hatching)
- Horizontal roller at node 4 (triangle symbol, free vertical)
- Force arrows at nodes 2 and 3 pointing rightward, labeled F/2 = 500 kN
- Coordinate axes

Node/element connectivity table:

| Elem | Nodes | DOFs |
|---|---|---|
| 1 | 1, 2, 3 | 1,2,3,4,5,6 |
| 2 | 1, 3, 4 | 1,2,5,6,7,8 |

### Step 2 — Element geometry: areas and B-matrices

**Element 1** (nodes 1=(0,0), 2=(1,0), 3=(1,1)):
- b1 = y2−y3 = −1, b2 = y3−y1 = 1, b3 = y1−y2 = 0
- c1 = x3−x2 = 0, c2 = x1−x3 = −1, c3 = x2−x1 = 1
- 2Ae = 1 → Ae = 0.5 m²

B¹ = [−1, 0, 1, 0, 0, 0; 0, 0, 0, −1, 0, 1; 0, −1, −1, 1, 1, 0]

**Element 2** (local: 1→global 1=(0,0), 2→global 3=(1,1), 3→global 4=(0,1)):
- b1 = y2−y3 = 0, b2 = y3−y1 = 1, b3 = y1−y2 = −1
- c1 = x3−x2 = −1, c2 = x1−x3 = 0, c3 = x2−x1 = 1
- 2Ae = 1 → Ae = 0.5 m²

B² = [0, 0, 1, 0, −1, 0; 0, −1, 0, 0, 0, 1; −1, 0, 0, 1, 1, −1]

Both matrices written as explicit 3×6 LaTeX arrays.

### Step 3 — D-matrix

Plane stress:
D = E/(1−ν²) × [1, ν, 0; ν, 1, 0; 0, 0, (1−ν)/2]
  = 219,780.2 MPa × [1, 0.3, 0; 0.3, 1, 0; 0, 0, 0.35]

Written numerically (entries in Pa for consistency with load vector).

### Step 4 — Element stiffness matrices

ke = t × Ae × Bᵀ D B for each element.

Both 6×6 matrices written numerically (MATLAB-verified, entries in N/m).

### Step 5 — Global assembly

8×8 global K assembled from:
- Element 1 → DOFs {1,2,3,4,5,6}
- Element 2 → DOFs {1,2,5,6,7,8}

DOFs 1,2,5,6 (nodes 1 and 3) receive contributions from both elements. Describe overlap pattern in prose; full 8×8 matrix computed by MATLAB (referenced via `\Cref{lst:cst2d}`).

### Step 6 — Load vector and BC reduction

Global load vector: F3 = F5 = 500,000 N, all others zero.

Apply BCs: fix DOFs {1,2,7} (u1=v1=u4=0). Free DOFs: {3,4,5,6,8} = {u2,v2,u3,v3,v4}.

Write the reduced 5×5 system Kff uf = Ff.

### Step 7 — Solution (MATLAB-verified)

Nodal displacements (from `matlab/fem_cst_2d.m`):
- u2 = σx/E = 100×10⁶/200×10⁹ = 5×10⁻⁴ m (exact)
- v2 = 0 (by symmetry / BCs)
- u3 = 5×10⁻⁴ m
- v3 = −ν σx/E = −1.5×10⁻⁴ m
- v4 = −1.5×10⁻⁴ m

All values to be MATLAB-verified before committing to LaTeX.

### Step 8 — Stress recovery

For each element:
- εe = Be ue → constant strain vector
- σe = D εe → constant stress vector

Both elements should give:
- σx = 100 MPa
- σy = 0
- τxy = 0

Written out numerically for element 1 (element 2 symmetric).

### Step 9 — Analytical verification

Exact solution for uniform uniaxial tension:
- σx = F/(t×L) = 1×10⁶/(0.01×1) = 100 MPa ✓
- σy = τxy = 0 ✓
- u(x,y) = (σx/E)×x → u2 = u3 = 5×10⁻⁴ m ✓
- v(x,y) = −(ν σx/E)×y → v3 = v4 = −1.5×10⁻⁴ m ✓

Key remark: CST recovers the exact solution because the true displacement field is linear — the same polynomial space used by the CST shape functions.

---

## MATLAB Script

**File:** `matlab/fem_cst_2d.m` (full rewrite of existing stub)

Contents:
- Problem parameters (E, ν, t, node coordinates, connectivity)
- Loop over elements: compute b, c, Ae, B, ke
- Assemble global K
- Apply BCs, solve
- Print nodal displacements
- Compute and print element stresses (both elements)
- Plot: mesh + deformed shape (magnified), with node labels
- Save PDF + PNG to `../figures/cst_example_results.pdf/.png`

Listing label: `lst:cst2d` (reuse existing label — appendix section already exists, just update its content)

---

## Changes to Existing Files

| File | Change |
|---|---|
| `chapters/ch04_2d_continuum.tex` | Remove §7.4 (lines 185–200); add §7.2.4 subsection inside §7.2 |
| `matlab/fem_cst_2d.m` | Full rewrite (new file — old content was only in lstlisting block) |
| `appendix/matlab_listings.tex` | Replace embedded `lstlisting` block with `\lstinputlisting{matlab/fem_cst_2d.m}` (consistent with portal-frame listing style) |

---

## Constraints

- All numerics MATLAB-verified before committing to LaTeX
- No toolboxes — core MATLAB only
- Units: SI throughout (Pa, m, N)
- Figure saved as PDF (vector) + PNG (300 dpi) to `../figures/`
- `lstlisting` block: no UTF-8 em-dashes inside the listing
- The label `ex:cst_tension` is reused — check no other chapter references it (currently only referenced by `lst:cst2d` caption, which will be updated)
