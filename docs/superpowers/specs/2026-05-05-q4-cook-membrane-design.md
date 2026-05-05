# Design: Q4 Cook Membrane Worked Example (§7.3.4)

**Date:** 2026-05-05
**File:** `chapters/ch04_2d_continuum.tex`
**Goal:** Add §7.3.4 — a full Q4 worked example (Cook membrane, one element) inside the existing §7.3 Q4 section.

---

## Decisions

| Question | Decision |
|---|---|
| Location | §7.3.4 subsection inside §7.3 Q4 section |
| Problem | Cook membrane, single Q4 element |
| Material | E = 70,000 N/mm², ν = 0.3, t = 1 mm, plane stress |
| Loading | P = 100 N total vertical (shear) force on right edge → 50 N at each of nodes 2 and 3 |
| Jacobian depth | Approach C: full symbolic J(ξ,η) derived analytically; Gauss-point table numerically |
| B-matrix depth | One representative Gauss point shown in full; others via MATLAB |
| k^e | Gauss sum formula stated; 8×8 result MATLAB-verified |

---

## Problem Definition

**Cook membrane geometry (mm):**

| Node | (x, y) | (ξ, η) | Global DOFs |
|---|---|---|---|
| 1 | (0, 0) | (−1, −1) | 1 (u₁), 2 (v₁) |
| 2 | (48, 44) | (+1, −1) | 3 (u₂), 4 (v₂) |
| 3 | (48, 60) | (+1, +1) | 5 (u₃), 6 (v₃) |
| 4 | (0, 44) | (−1, +1) | 7 (u₄), 8 (v₄) |

**BCs:** Left edge clamped — u₁=v₁=u₄=v₄=0 (DOFs 1,2,7,8 fixed). Free DOFs: {3,4,5,6}.

**Load:** Right edge (nodes 2→3, length = 60−44 = 16 mm). Uniform traction → consistent nodal forces F₄ = F₆ = 50 N. All other load components zero.

---

## Key Symbolic Results (Analytically Derived)

### Shape function natural derivatives

| Node | ∂Nᵢ/∂ξ | ∂Nᵢ/∂η |
|---|---|---|
| 1 | −(1−η)/4 | −(1−ξ)/4 |
| 2 | +(1−η)/4 | −(1+ξ)/4 |
| 3 | +(1+η)/4 | +(1+ξ)/4 |
| 4 | −(1+η)/4 | +(1−ξ)/4 |

### Jacobian J(ξ,η)

J₁₁ = ∂x/∂ξ = 24 (constant)
J₁₂ = ∂y/∂ξ = 15 − 7η (linear in η)
J₂₁ = ∂x/∂η = 0 (exact zero — left and right edges both vertical)
J₂₂ = ∂y/∂η = 15 − 7ξ (linear in ξ)

**det|J| = 24(15 − 7ξ)** — varies with ξ only, not η.

**J⁻¹ = [1/24, −(15−7η)/(24(15−7ξ)); 0, 1/(15−7ξ)]**

Key remark: J₂₁=0 because x-coordinates are identical on each vertical edge (x₁=x₄=0, x₂=x₃=48), so ∂x/∂η = 0 everywhere.

Key remark: det|J| depends only on ξ — so Gauss points at the same ξ but different η have the same det|J|.

### Gauss-point numerical values (g = 1/√3 ≈ 0.5774)

| GP | (ξ, η) | J₁₂ = 15−7η | J₂₂ = 15−7ξ | det\|J\| | J⁻¹[1,2] | J⁻¹[2,2] |
|---|---|---|---|---|---|---|
| 1 | (−g, −g) | 19.04 | 19.04 | 457.0 | −1/24 | 0.05251 |
| 2 | (+g, −g) | 19.04 | 10.96 | 263.0 | −0.07238 | 0.09125 |
| 3 | (+g, +g) | 10.96 | 10.96 | 263.0 | −1/24 | 0.09125 |
| 4 | (−g, +g) | 10.96 | 19.04 | 457.0 | −0.02398 | 0.05251 |

(J⁻¹[1,1] = 1/24 = 0.04167 and J⁻¹[2,1] = 0 at all Gauss points.)

Note: GPs 1 and 4 share det|J|=457.0 (both at ξ=−g); GPs 2 and 3 share det|J|=263.0 (both at ξ=+g).

### B-matrix at GP1 (ξ=η=−g)

Physical derivatives via J⁻¹ ∂N/∂ξ:

∂N/∂x = [0, 1/48, 0, −1/48] mm⁻¹ = [0, 0.02083, 0, −0.02083] mm⁻¹
∂N/∂y = [−0.02071, −0.005549, 0.005549, 0.02071] mm⁻¹

B(GP1) [3×8, mm⁻¹]:
```
[  0,       0,      0.02083,  0,      0,        0,     -0.02083,  0      ]
[  0,      -0.02071, 0,      -0.005549, 0,      0.005549, 0,      0.02071]
[-0.02071,  0,      -0.005549, 0.02083, 0.005549, 0,      0.02071,-0.02083]
```

---

## Section Structure

### §7.3.4 Worked Example: Single Q4 Element (Cook Membrane)

Added as `\subsection` immediately after the Gauss quadrature table in §7.3 (after line 501 in current file).

**Step 1:** Problem description + TikZ figure
- Tapered quadrilateral with node labels 1–4, element number, coordinates
- Hatched clamped left edge (nodes 1 and 4)
- Force arrows on right edge labeled F/2 = 50 N at nodes 2 and 3
- Dimension annotations

**Step 2:** Shape function natural derivatives — table of ∂Nᵢ/∂ξ and ∂Nᵢ/∂η

**Step 3:** Jacobian derivation
- Compute J₁₁, J₁₂, J₂₁, J₂₂ symbolically
- State det|J| and J⁻¹ analytically
- Remark on J₂₁=0 and det|J| varying with ξ only

**Step 4:** Gauss-point table — compact 4-row table: (ξᵢ,ηⱼ), J matrix values, det|J|, J⁻¹ entries (MATLAB-verified)

**Step 5:** B-matrix at GP1
- Show ∂Nᵢ/∂x and ∂Nᵢ/∂y computed via J⁻¹
- Write out full 3×8 B(GP1) numerically

**Step 6:** D-matrix (plane stress, E=70000, ν=0.3 — numerical values)

**Step 7:** Element stiffness
- State formula: k^e = t Σᵢ Σⱼ wᵢwⱼ Bᵀ D B |J| at (ξᵢ,ηⱼ), all weights=1
- Remark: B non-constant → no closed form; must use numerical integration
- Present MATLAB-verified 8×8 k^e (×10³ N/mm scaling)

**Step 8:** BCs and reduced system
- Fix DOFs {1,2,7,8}; free DOFs {3,4,5,6}
- Load vector: F₄=F₆=50 N
- Show 4×4 Kff and RHS (MATLAB-verified)

**Step 9:** Solution + stress recovery
- MATLAB-verified displacements: u₂, v₂, u₃, v₃
- Stresses σ = D B(ξᵢ,ηⱼ) u^e at all 4 Gauss points in a table
- Remark: stresses vary across element (unlike CST with constant stress)
- Figure: deformed shape + σx or σy distribution from MATLAB

---

## MATLAB Script

**File:** `matlab/fem_q4_cook.m` (new file)

Contents:
- Problem setup (coords, D, BCs, load)
- Loop: compute B at each Gauss point, accumulate k^e via Gauss sum
- Print: J at each GP, det|J|, B at GP1, full k^e, Kff, u, stresses
- Figure: undeformed mesh + deformed shape (scaled); color plot of σy
- Save: `../figures/q4_cook_results.pdf` and `.png`

Listing label: `lst:q4cook` — new section in appendix using `\lstinputlisting`

---

## Changes to Existing Files

| File | Change |
|---|---|
| `chapters/ch04_2d_continuum.tex` | Append §7.3.4 subsection after Gauss quadrature table |
| `matlab/fem_q4_cook.m` | Create new script |
| `appendix/matlab_listings.tex` | Add new section with `\lstinputlisting{matlab/fem_q4_cook.m}` |

---

## Constraints

- All numerics MATLAB-verified before committing to LaTeX
- No toolboxes — core MATLAB only
- Units: N and mm throughout (E in N/mm², coordinates in mm, forces in N)
- TikZ: use `zc` not `z` for any custom style; no `\alph` redefinition
- lstlisting: no UTF-8 em-dashes in MATLAB script
- New labels: `sec:q4_example`, `ex:q4_cook`, `fig:q4_cook`, `fig:q4_cook_results`, `eq:jacobian_cook`, `lst:q4cook`
- Cross-refs in example: `eq:q4_shape`, `eq:jacobian`, `eq:q4_stiffness`, `eq:D_planestress`
