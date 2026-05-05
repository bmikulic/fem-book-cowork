# Q4 Cook Membrane Worked Example — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add §7.3.4 — a full Q4 Cook membrane worked example (one element, symbolic Jacobian, 2×2 Gauss integration, MATLAB-verified solution and stresses) — as a subsection inside §7.3 in `chapters/ch04_2d_continuum.tex`.

**Architecture:** Create the MATLAB script first (verified numbers); then write the LaTeX in three passes (setup+Jacobian → Gauss table+B+k^e → BCs+solution+stresses); update the appendix listing last; build and verify.

**Tech Stack:** LaTeX (`latexmk -pdf`), MATLAB (R2016b+), TikZ, BibTeX.

---

## Pre-Verified Numbers (MATLAB-Confirmed)

All values below are exact MATLAB output for E=70000 N/mm², ν=0.3, t=1 mm, nodes (0,0)-(48,44)-(48,60)-(0,44).

**g = 1/√3 ≈ 0.5774**

**D-matrix (N/mm²):**
```
[76923,  23077,     0]
[23077,  76923,     0]
[    0,      0, 26923]
```

**Gauss-point data** (J₁₁=24 and J₂₁=0 at all GPs; J⁻¹[1,1]=0.04167, J⁻¹[2,1]=0 at all GPs):
```
GP1 (−g,−g): J₁₂=19.04, J₂₂=19.04, det|J|=457.0, J⁻¹[1,2]=−0.04167, J⁻¹[2,2]=0.05252
GP2 (+g,−g): J₁₂=19.04, J₂₂=10.96, det|J|=263.0, J⁻¹[1,2]=−0.07240, J⁻¹[2,2]=0.09125
GP3 (+g,+g): J₁₂=10.96, J₂₂=10.96, det|J|=263.0, J⁻¹[1,2]=−0.04167, J⁻¹[2,2]=0.09125
GP4 (−g,+g): J₁₂=10.96, J₂₂=19.04, det|J|=457.0, J⁻¹[1,2]=−0.02398, J⁻¹[2,2]=0.05252
```

**B at GP1 (×10⁻³ mm⁻¹, DOFs u₁,v₁,u₂,v₂,u₃,v₃,u₄,v₄):**
```
Row 1: [  0,       0,     20.83,   0,      0,      0,   −20.83,   0    ]
Row 2: [  0,     −20.71,   0,    −5.549,   0,    5.549,   0,     20.71 ]
Row 3: [−20.71,   0,      −5.549, 20.83,  5.549,   0,    20.71, −20.83 ]
```

**k^e (×10³ N/mm, symmetric, DOFs u₁,v₁,u₂,v₂,u₃,v₃,u₄,v₄):**
```
[ 14.33, -1.28,   0.98,  -9.95,  -0.98,  -3.51, -14.33,  14.74]
[ -1.28, 34.78,  -8.02,  19.75,  -3.51, -19.75,  12.82, -34.78]
[  0.98, -8.02,  85.76, -34.66, -37.68,  21.20, -49.06,  21.49]
[ -9.95, 19.75, -34.66,  77.89,  23.13, -61.06,  21.49, -36.58]
[ -0.98, -3.51, -37.68,  23.13,  37.68,  -9.66,   0.98,  -9.95]
[ -3.51,-19.75,  21.20, -61.06,  -9.66,  61.06,  -8.02,  19.75]
[-14.33, 12.82, -49.06,  21.49,   0.98,  -8.02,  62.40, -26.28]
[ 14.74,-34.78,  21.49, -36.58,  -9.95,  19.75, -26.28,  51.60]
```

**Kff (×10³ N/mm, free DOFs u₂,v₂,u₃,v₃ = DOFs 3,4,5,6):**
```
[ 85.76, -34.66, -37.68,  21.20]
[-34.66,  77.89,  23.13, -61.06]
[-37.68,  23.13,  37.68,  -9.66]
[ 21.20, -61.06,  -9.66,  61.06]
```

**Load:** F₄ = F₆ = 50 N (v₂ and v₃ DOFs)

**Solution (mm):**
```
u₂ = −9.3×10⁻⁵,   v₂ = +8.210×10⁻³
u₃ = −2.926×10⁻³, v₃ = +8.597×10⁻³
```

**Stresses at Gauss points (N/mm²):**
```
GP1: σx = −0.100,  σy = +0.121,  τxy = +4.181
GP2: σx = +5.436,  σy = +2.608,  τxy = +1.608
GP3: σx = +0.173,  σy = +1.029,  τxy = +1.860
GP4: σx = −3.128,  σy = −0.788,  τxy = +4.327
```

---

## File Map

| File | Action | Purpose |
|---|---|---|
| `matlab/fem_q4_cook.m` | Create | Cook membrane Q4 script; generates figure |
| `chapters/ch04_2d_continuum.tex` | Modify (append) | Add §7.3.4 after line 501 |
| `appendix/matlab_listings.tex` | Modify (append) | Add `lst:q4cook` section |

---

## Task 1: Write and Run the MATLAB Script

**Files:**
- Create: `matlab/fem_q4_cook.m`

- [ ] **Step 1: Create `matlab/fem_q4_cook.m` with this exact content**

```matlab
% Chapter 7 -- Q4 Cook membrane under in-plane shear
% One Q4 element: nodes 1(0,0) 2(48,44) 3(48,60) 4(0,44) [mm]
% E=70000 N/mm^2, nu=0.3, t=1 mm, plane stress
% BCs: left edge clamped (nodes 1,4); load: 50 N vertical at nodes 2,3
clear; clc; close all;

%% Parameters
E  = 70000;  nu = 0.3;  t = 1;
xy = [0,0; 48,44; 48,60; 0,44];

%% Constitutive matrix (plane stress)
D = E/(1-nu^2) * [1,nu,0; nu,1,0; 0,0,(1-nu)/2];

%% 2x2 Gauss points and weights
g    = 1/sqrt(3);
gpts = [-g,-g; g,-g; g,g; -g,g];
w    = [1,1,1,1];

%% Assemble element stiffness
ke   = zeros(8);
B_all = cell(4,1);
J_all = cell(4,1);

for gp = 1:4
    xi  = gpts(gp,1);  eta = gpts(gp,2);

    % Natural derivatives of shape functions
    dN_dxi  = [-(1-eta), (1-eta), (1+eta), -(1+eta)] / 4;
    dN_deta = [-(1-xi), -(1+xi),  (1+xi),  (1-xi) ] / 4;

    % Jacobian and its inverse
    J    = [dN_dxi;  dN_deta] * xy;
    detJ = det(J);
    Jinv = inv(J);

    % Physical derivatives
    dN_dx = Jinv(1,:) * [dN_dxi; dN_deta];
    dN_dy = Jinv(2,:) * [dN_dxi; dN_deta];

    % Strain-displacement matrix B (3x8)
    B = zeros(3,8);
    B(1,1:2:7) = dN_dx;
    B(2,2:2:8) = dN_dy;
    B(3,1:2:7) = dN_dy;
    B(3,2:2:8) = dN_dx;

    ke = ke + t * w(gp) * detJ * B' * D * B;
    B_all{gp} = B;
    J_all{gp} = J;
end

%% Print Gauss-point data
fprintf('=== JACOBIAN AT GAUSS POINTS ===\n')
for gp = 1:4
    J = J_all{gp};
    fprintf('GP%d: J=[%.4f,%.4f;%.4f,%.4f]  detJ=%.4f\n', ...
            gp, J(1,1),J(1,2),J(2,1),J(2,2), det(J_all{gp}))
end
fprintf('\n=== B MATRIX AT GP1 (x1e3) ===\n')
disp(B_all{1}*1e3)
fprintf('\n=== ke (x1e-3 N/mm) ===\n')
disp(ke/1e3)

%% BCs and load
fixed = [1,2,7,8];   % u1,v1,u4,v4
free  = setdiff(1:8, fixed);
F = zeros(8,1);
F(4) = 50;   % v2
F(6) = 50;   % v3

u = zeros(8,1);
u(free) = ke(free,free) \ F(free);

fprintf('\n=== SOLUTION (mm) ===\n')
fprintf('u2=%+.6f  v2=%+.6f\n', u(3), u(4))
fprintf('u3=%+.6f  v3=%+.6f\n', u(5), u(6))

%% Stresses at Gauss points
fprintf('\n=== STRESSES AT GAUSS POINTS (N/mm^2) ===\n')
for gp = 1:4
    sig = D * B_all{gp} * u;
    fprintf('GP%d: sx=%+.4f  sy=%+.4f  txy=%+.4f\n', gp, sig(1), sig(2), sig(3))
end

%% Figure: mesh and deformed shape
scale = 5000;
fig = figure('Units','centimeters','Position',[2,2,14,12]);
hold on; axis equal; grid on;

% Undeformed
patch('XData', xy([1,2,3,4,1],1), 'YData', xy([1,2,3,4,1],2), ...
      'FaceColor','none', 'EdgeColor','k', 'LineWidth',1.2)

% Deformed
xy_def = xy + scale * reshape(u,2,[])';
patch('XData', xy_def([1,2,3,4,1],1), 'YData', xy_def([1,2,3,4,1],2), ...
      'FaceColor','none', 'EdgeColor','b', 'LineStyle','--', 'LineWidth',1.5)

% Nodes
for n = 1:4
    plot(xy(n,1), xy(n,2), 'ko', 'MarkerSize', 5, 'MarkerFaceColor','k')
    text(xy(n,1)+1.5, xy(n,2)+1, sprintf('%d',n), 'FontSize',10)
end

xlabel('x (mm)'); ylabel('y (mm)')
title(sprintf('Q4 Cook membrane: undeformed (solid) and deformed (dashed, \\times%d)', scale))
legend({'Undeformed','Deformed'}, 'Location','northwest')

exportgraphics(fig, '../figures/q4_cook_results.pdf', 'ContentType','vector')
exportgraphics(fig, '../figures/q4_cook_results.png', 'Resolution', 300)
fprintf('\nFigures saved to ../figures/\n')
```

- [ ] **Step 2: Run in MATLAB**

Navigate to `matlab/` and run:
```
fem_q4_cook
```

- [ ] **Step 3: Confirm output matches pre-verified values**

Expected (summary):
```
GP1: J=[24.0000,19.0415;0.0000,19.0415]  detJ=456.9948
GP2: J=[24.0000,19.0415;0.0000,10.9585]  detJ=263.0052
...
u2=-0.000093  v2=+0.008210
u3=-0.002926  v3=+0.008597
GP1: sx=-0.0995  sy=+0.1208  txy=+4.1814
```

- [ ] **Step 4: Commit**

```
git add matlab/fem_q4_cook.m
git commit -m "Ch.7: add Q4 Cook membrane MATLAB script"
```

---

## Task 2: Write §7.3.4 Part A — Problem Setup, Jacobian Derivation, and Gauss-Point Table

**Files:**
- Modify: `chapters/ch04_2d_continuum.tex`

The current §7.3 ends after the Gauss quadrature table (around line 501). Append the following before `\end{document}` or after the Gauss table closing `\end{center}`.

- [ ] **Step 1: Append §7.3.4 opening through Step 3 (symbolic Jacobian)**

```latex

% =============================================================================
\subsection{Worked Example: Single Q4 Element (Cook Membrane)}
\label{sec:q4_example}
% =============================================================================

\begin{example}[Cook membrane --- single Q4 element under in-plane shear]
\label{ex:q4_cook}

The Cook membrane is a tapered quadrilateral plate used as a benchmark for
2D elements. A single Q4 element discretises the plate, as shown in
\cref{fig:q4_cook}. Material properties: $E = 70{,}000\ \text{N/mm}^2$,
$\nu = 0.3$, $t = 1\ \text{mm}$, plane stress. The left edge is clamped
(nodes~1 and~4). A total vertical force $P = 100\ \text{N}$ acts uniformly
on the right edge; the consistent nodal forces are $F/2 = 50\ \text{N}$
upward at each of nodes~2 and~3.

\begin{figure}[htbp]
\centering
\begin{tikzpicture}[>=Stealth, x=0.07cm, y=0.07cm,
    force/.style={->,ultra thick,red!75!black}]

  %% Plate outline
  \fill[gray!12]  (0,0) -- (48,44) -- (48,60) -- (0,44) -- cycle;
  \draw[thick]    (0,0) -- (48,44) -- (48,60) -- (0,44) -- cycle;

  %% Clamped left edge
  \fill[pattern=north east lines,pattern color=gray!70]
       (-5,0) rectangle (0,44);
  \draw[thick] (0,0) -- (0,44);

  %% Nodes
  \filldraw (0, 0) circle (3pt) node[left=4pt] {\small\bfseries 1 $(0,0)$};
  \filldraw (48,44) circle (3pt) node[right=4pt] {\small\bfseries 2 $(48,44)$};
  \filldraw (48,60) circle (3pt) node[right=4pt] {\small\bfseries 3 $(48,60)$};
  \filldraw (0, 44) circle (3pt) node[left=4pt] {\small\bfseries 4 $(0,44)$};

  %% Forces at nodes 2 and 3
  \draw[force] (48,44) -- ++(0,14) node[right=2pt,midway] {$50\ \text{N}$};
  \draw[force] (48,60) -- ++(0,14) node[right=2pt,midway] {$50\ \text{N}$};

  %% Dimensions
  \draw[<->,thin] (0,-8) -- (48,-8)
       node[fill=white,inner sep=1pt,midway] {$48\ \text{mm}$};
  \draw[<->,thin] (-10,0) -- (-10,44)
       node[fill=white,inner sep=1pt,midway,rotate=90] {$44\ \text{mm}$};
  \draw[<->,thin] (54,44) -- (54,60)
       node[fill=white,inner sep=1pt,midway,rotate=90] {$16\ \text{mm}$};

\end{tikzpicture}
\caption{Single Q4 element for the Cook membrane benchmark.
Node coordinates in mm. Left edge clamped; $50\ \text{N}$ vertical
nodal forces at nodes~2 and~3.}
\label{fig:q4_cook}
\end{figure}

\medskip
\noindent\textbf{DOF ordering:}
$[u_1,v_1,u_2,v_2,u_3,v_3,u_4,v_4] = \text{DOFs}\ [1\ldots8]$.\\
Fixed: $\{1,2,7,8\}$. Free: $\{3,4,5,6\}$ ($u_2,v_2,u_3,v_3$).

\medskip
\noindent\textbf{Step 1: Shape function natural derivatives.}

From \cref{eq:q4_shape}, the derivatives are:

\begin{table}[htbp]
\centering
\small
\caption{Natural derivatives of Q4 shape functions.}
\label{tab:q4_dN}
\begin{tabular}{@{}ccc@{}}
\toprule
Node $i$ & $\partial N_i/\partial\xi$ & $\partial N_i/\partial\eta$ \\
\midrule
1 & $-(1-\eta)/4$ & $-(1-\xi)/4$ \\
2 & $+(1-\eta)/4$ & $-(1+\xi)/4$ \\
3 & $+(1+\eta)/4$ & $+(1+\xi)/4$ \\
4 & $-(1+\eta)/4$ & $+(1-\xi)/4$ \\
\bottomrule
\end{tabular}
\end{table}

\medskip
\noindent\textbf{Step 2: Symbolic Jacobian.}

Substituting the node coordinates $(x_1,y_1)=(0,0)$,
$(x_2,y_2)=(48,44)$, $(x_3,y_3)=(48,60)$, $(x_4,y_4)=(0,44)$
into \cref{eq:jacobian}:
\begin{align}
    J_{11} &= \frac{\partial x}{\partial\xi}
    = \sum_i \frac{\partial N_i}{\partial\xi} x_i
    = \frac{-(1-\eta)}{4}\cdot0 + \frac{(1-\eta)}{4}\cdot48
      + \frac{(1+\eta)}{4}\cdot48 + \frac{-(1+\eta)}{4}\cdot0
    = 24 \notag\\
    J_{12} &= \frac{\partial y}{\partial\xi}
    = \frac{(1-\eta)}{4}\cdot44 + \frac{(1+\eta)}{4}\cdot60
      - \frac{(1+\eta)}{4}\cdot44
    = 15 - 7\eta \notag\\
    J_{21} &= \frac{\partial x}{\partial\eta}
    = \frac{-(1+\xi)}{4}\cdot48 + \frac{(1+\xi)}{4}\cdot48
    = 0 \notag\\
    J_{22} &= \frac{\partial y}{\partial\eta}
    = \frac{-(1+\xi)}{4}\cdot44 + \frac{(1+\xi)}{4}\cdot60
      + \frac{(1-\xi)}{4}\cdot44
    = 15 - 7\xi \notag
\end{align}

\begin{equation}
    \matr{J}(\xi,\eta) =
    \begin{bmatrix} 24 & 15-7\eta \\ 0 & 15-7\xi \end{bmatrix},
    \qquad
    \det|\matr{J}| = 24(15-7\xi)
    \label{eq:jacobian_cook}
\end{equation}

\begin{equation}
    \matr{J}^{-1} =
    \frac{1}{24(15-7\xi)}
    \begin{bmatrix} 15-7\xi & -(15-7\eta) \\ 0 & 24 \end{bmatrix}
    =
    \begin{bmatrix}
        \dfrac{1}{24} & \dfrac{-(15-7\eta)}{24(15-7\xi)} \\[8pt]
        0 & \dfrac{1}{15-7\xi}
    \end{bmatrix}
    \label{eq:Jinv_cook}
\end{equation}

\begin{remark}
$J_{21} = 0$ because the left and right edges are both vertical
($x_1 = x_4 = 0$, $x_2 = x_3 = 48\ \text{mm}$), so $x$ depends only
on $\xi$, not $\eta$. Consequently $\det|\matr{J}| = 24(15-7\xi)$
varies with $\xi$ only --- Gauss points at the same $\xi$ share the
same determinant regardless of $\eta$.
\end{remark}

\medskip
\noindent\textbf{Step 3: Gauss-point table.}

For $2\times2$ quadrature with $g = 1/\sqrt{3}$, $w_i = 1$:

\begin{table}[htbp]
\centering
\small
\caption{Jacobian data at the four Gauss points ($g=1/\sqrt{3}\approx0.5774$).
$J_{11}=24$ and $J_{21}=0$ at all points; $J^{-1}_{11}=1/24\approx0.04167$
and $J^{-1}_{21}=0$ at all points.}
\label{tab:q4_gauss}
\begin{tabular}{@{}cccccccc@{}}
\toprule
GP & $(\xi,\eta)$ & $J_{12}=15{-}7\eta$ & $J_{22}=15{-}7\xi$
   & $\det|\matr{J}|$ & $J^{-1}_{12}$ & $J^{-1}_{22}$ \\
\midrule
1 & $(-g,-g)$ & 19.04 & 19.04 & 457.0 & $-1/24$ & 0.05252 \\
2 & $(+g,-g)$ & 19.04 & 10.96 & 263.0 & $-0.07240$ & 0.09125 \\
3 & $(+g,+g)$ & 10.96 & 10.96 & 263.0 & $-1/24$ & 0.09125 \\
4 & $(-g,+g)$ & 10.96 & 19.04 & 457.0 & $-0.02398$ & 0.05252 \\
\bottomrule
\end{tabular}
\end{table}

\end{example}
```

Note: `\end{example}` temporarily closes the example so the file compiles. Tasks 3 and 4 will reopen it.

- [ ] **Step 2: Compile**

```
cd "C:\Users\Korisnik\Documents\FEM_Book Cowork"
latexmk -pdf main.tex
```

Expected: compiles; §7.3.4 appears with TikZ figure and two tables. Warning for `lst:q4cook` is expected.

- [ ] **Step 3: Commit**

```
git add chapters/ch04_2d_continuum.tex
git commit -m "Ch.7: §7.3.4 Part A -- Cook membrane setup, Jacobian, Gauss table"
```

---

## Task 3: Write §7.3.4 Part B — B-matrix, D-matrix, k^e

**Files:**
- Modify: `chapters/ch04_2d_continuum.tex`

Find the temporary `\end{example}` just added and replace it with the following content + a new `\end{example}`.

- [ ] **Step 1: Remove temporary `\end{example}` and append Steps 4–6**

```latex

\medskip
\noindent\textbf{Step 4: Strain-displacement matrix at GP1.}

The physical derivatives $\partial N_i/\partial x$ and $\partial N_i/\partial y$
follow from $[\partial N/\partial x,\; \partial N/\partial y]\transpose =
\matr{J}^{-1}[\partial N/\partial\xi,\; \partial N/\partial\eta]\transpose$.
At GP1 $(\xi\!=\!\eta\!=\!-g)$, substituting $J^{-1}_{12} = -1/24$
and $J^{-1}_{22} = 1/19.04$:

\begin{align*}
    \frac{\partial N_i}{\partial x}
    &= \frac{1}{24}\frac{\partial N_i}{\partial\xi}
      -\frac{1}{24}\frac{\partial N_i}{\partial\eta},\\[4pt]
    \frac{\partial N_i}{\partial y}
    &= \frac{1}{19.04}\frac{\partial N_i}{\partial\eta}.
\end{align*}

Evaluating for all four nodes yields (in $10^{-3}\ \text{mm}^{-1}$):
\[
    \frac{\partial\vect{N}}{\partial x}
    = [0,\; 20.83,\; 0,\; -20.83],\qquad
    \frac{\partial\vect{N}}{\partial y}
    = [-20.71,\; -5.549,\; 5.549,\; 20.71]
\]

Assembling $\matr{B} = \matr{B}(\xi,\eta)$ at GP1:
\begin{equation}
    \matr{B}(\text{GP1}) = 10^{-3}\times
    \begin{bmatrix}
         0     & 0      & 20.83  & 0      & 0     & 0     & -20.83 & 0     \\
         0     & -20.71 & 0      & -5.549 & 0     & 5.549 & 0      & 20.71 \\
       -20.71  & 0      & -5.549 & 20.83  & 5.549 & 0     & 20.71  & -20.83
    \end{bmatrix}\ \text{mm}^{-1}
    \label{eq:B_q4_gp1}
\end{equation}

\begin{remark}
Unlike the CST where $\matr{B}$ is constant, the Q4 $\matr{B}(\xi,\eta)$
changes at every Gauss point because the physical derivatives depend on
$\matr{J}^{-1}(\xi,\eta)$, which varies over the element. Gauss quadrature
is therefore \emph{essential} --- no closed-form integral exists.
\end{remark}

\medskip
\noindent\textbf{Step 5: Constitutive matrix.}

\begin{equation}
    \matr{D} = \frac{70{,}000}{1-0.09}
    \begin{bmatrix}1&0.3&0\\0.3&1&0\\0&0&0.35\end{bmatrix}
    =
    \begin{bmatrix}
        76{,}923 & 23{,}077 & 0 \\
        23{,}077 & 76{,}923 & 0 \\
        0 & 0 & 26{,}923
    \end{bmatrix}\ \text{N/mm}^2
    \label{eq:D_q4_ex}
\end{equation}

\medskip
\noindent\textbf{Step 6: Element stiffness matrix.}

The $8\times8$ element stiffness is accumulated over the four Gauss points
(\cref{eq:q4_stiffness}), with all weights $w_i = 1$:
\[
    \matr{k}^e = t\sum_{p=1}^{4} \matr{B}_p\transpose\,\matr{D}\,\matr{B}_p\,
    \det|\matr{J}_p|
\]
Since $\matr{B}$ varies at each Gauss point, this sum must be computed
numerically (\Cref{lst:q4cook}). The MATLAB-verified result is:

\begin{equation}
    \matr{k}^e = 10^3\times
    \begin{bmatrix*}[r]
         14.33 &  -1.28 &   0.98 &  -9.95 &  -0.98 &  -3.51 & -14.33 &  14.74 \\
         -1.28 &  34.78 &  -8.02 &  19.75 &  -3.51 & -19.75 &  12.82 & -34.78 \\
          0.98 &  -8.02 &  85.76 & -34.66 & -37.68 &  21.20 & -49.06 &  21.49 \\
         -9.95 &  19.75 & -34.66 &  77.89 &  23.13 & -61.06 &  21.49 & -36.58 \\
         -0.98 &  -3.51 & -37.68 &  23.13 &  37.68 &  -9.66 &   0.98 &  -9.95 \\
         -3.51 & -19.75 &  21.20 & -61.06 &  -9.66 &  61.06 &  -8.02 &  19.75 \\
        -14.33 &  12.82 & -49.06 &  21.49 &   0.98 &  -8.02 &  62.40 & -26.28 \\
         14.74 & -34.78 &  21.49 & -36.58 &  -9.95 &  19.75 & -26.28 &  51.60
    \end{bmatrix*}\ \tfrac{\text{N}}{\text{mm}}
    \label{eq:ke_q4_cook}
\end{equation}

\end{example}
```

- [ ] **Step 2: Compile**

```
latexmk -pdf main.tex
```

- [ ] **Step 3: Commit**

```
git add chapters/ch04_2d_continuum.tex
git commit -m "Ch.7: §7.3.4 Part B -- B-matrix at GP1, D-matrix, 8x8 ke"
```

---

## Task 4: Write §7.3.4 Part C — BCs, Solution, Stress Recovery, and Results Figure

**Files:**
- Modify: `chapters/ch04_2d_continuum.tex`

Remove the temporary `\end{example}` and append Steps 7–9.

- [ ] **Step 1: Remove `\end{example}` and append Steps 7–9**

```latex

\medskip
\noindent\textbf{Step 7: Boundary conditions and reduced system.}

Clamping the left edge sets DOFs $\{1,2,7,8\} = 0$. The four free DOFs
$\{3,4,5,6\}$ ($u_2,v_2,u_3,v_3$) satisfy the $4\times4$ reduced system:

\begin{equation}
    10^3\times
    \begin{bmatrix}
         85.76 & -34.66 & -37.68 &  21.20 \\
        -34.66 &  77.89 &  23.13 & -61.06 \\
        -37.68 &  23.13 &  37.68 &  -9.66 \\
         21.20 & -61.06 &  -9.66 &  61.06
    \end{bmatrix}
    \begin{bmatrix} u_2 \\ v_2 \\ u_3 \\ v_3 \end{bmatrix}
    =
    \begin{bmatrix} 0 \\ 50 \\ 0 \\ 50 \end{bmatrix}\ \text{N}
    \label{eq:Kff_q4}
\end{equation}

\medskip
\noindent\textbf{Step 8: Solution (MATLAB-verified, \Cref{lst:q4cook}).}

\begin{align*}
    u_2 &= -9.3\times10^{-5}\ \text{mm}, &
    v_2 &= +8.210\times10^{-3}\ \text{mm}, \\
    u_3 &= -2.926\times10^{-3}\ \text{mm}, &
    v_3 &= +8.597\times10^{-3}\ \text{mm}.
\end{align*}

The dominant response is vertical (in the direction of loading): both right-edge
nodes deflect upward by approximately $8.4\ \mu\text{m}$. The small horizontal
displacements ($\sim\!3\ \mu\text{m}$) reflect the shear-induced distortion of
the tapered geometry. The deformed shape is shown in \cref{fig:q4_cook_results}.

\begin{figure}[htbp]
\centering
\includegraphics[width=0.62\textwidth]{figures/q4_cook_results}
\caption{Cook membrane: undeformed (solid) and deformed shape (dashed,
magnified $\times5000$). Both right-edge nodes deflect vertically;
the tapered geometry induces small lateral movement.}
\label{fig:q4_cook_results}
\end{figure}

\medskip
\noindent\textbf{Step 9: Stress recovery at Gauss points.}

The stresses $\vect{\sigma}^e = \matr{D}\,\matr{B}(\xi_p,\eta_p)\,\vect{u}^e$
are evaluated at each Gauss point:

\begin{table}[htbp]
\centering
\small
\caption{Stresses at the four Gauss points (N/mm$^2$).}
\label{tab:q4_stresses}
\begin{tabular}{@{}ccccc@{}}
\toprule
GP & $(\xi,\eta)$ & $\sigma_x$ & $\sigma_y$ & $\tau_{xy}$ \\
\midrule
1 & $(-g,-g)$ & $-0.100$ & $+0.121$ & $+4.181$ \\
2 & $(+g,-g)$ & $+5.436$ & $+2.608$ & $+1.608$ \\
3 & $(+g,+g)$ & $+0.173$ & $+1.029$ & $+1.860$ \\
4 & $(-g,+g)$ & $-3.128$ & $-0.788$ & $+4.327$ \\
\bottomrule
\end{tabular}
\end{table}

\begin{remark}
The stresses vary significantly across the element --- a key difference
from the CST, which yields constant stress. The large $\tau_{xy}$ values
at GPs~1 and~4 (left side of the element) reflect the shear dominance
near the clamped edge, where bending is most constrained. This stress
variation within the element is what the Q4 brings over the CST: it can
represent linearly varying strains, which better captures the true
behaviour of the Cook membrane under in-plane shear.
\end{remark}

\end{example}
```

- [ ] **Step 2: Compile**

```
latexmk -pdf main.tex
```

Expected: all 9 steps present; figure `fig:q4_cook_results` loads from `figures/q4_cook_results.pdf`.

- [ ] **Step 3: Commit**

```
git add chapters/ch04_2d_continuum.tex
git commit -m "Ch.7: §7.3.4 Part C -- BCs, solution, stress recovery, complete example"
```

---

## Task 5: Add Appendix Listing and Final Build

**Files:**
- Modify: `appendix/matlab_listings.tex`

The appendix currently ends with `lst:q4cook` absent. Append after the last listing entry (after the `lst:cst2d` section for Chapter 7, around line 500).

- [ ] **Step 1: Read the end of `appendix/matlab_listings.tex` to find insertion point**

Use the Read tool to read from line 490 onward and identify the last entry.

- [ ] **Step 2: Append the Q4 Cook membrane listing**

After the last existing listing block, append:

```latex

% ---------------------------------------------------------------
\section{Chapter~7 --- Q4 Cook Membrane}
\label{lst:q4cook}

\lstinputlisting[
    caption={Single Q4 element for the Cook membrane benchmark
             (\Cref{ex:q4_cook}).},
    label={code:q4cook}
]{matlab/fem_q4_cook.m}
```

**Before appending:** verify `matlab/fem_q4_cook.m` contains no UTF-8 em-dashes. The file uses only `%--` double-hyphen separators and standard ASCII — it is clean.

- [ ] **Step 3: Final build**

```
cd "C:\Users\Korisnik\Documents\FEM_Book Cowork"
latexmk -pdf main.tex
```

- [ ] **Step 4: Verify the PDF**

Check:
- [ ] §7.3.4 appears in TOC as subsection of §7.3 Q4
- [ ] TikZ figure `fig:q4_cook` renders (tapered quad, nodes 1–4, clamped left, force arrows)
- [ ] Shape-derivative table `tab:q4_dN` present
- [ ] Symbolic Jacobian `eq:jacobian_cook` and inverse `eq:Jinv_cook` present
- [ ] Gauss-point table `tab:q4_gauss` with all 4 GPs and det|J| values 457.0 / 263.0
- [ ] B-matrix at GP1 `eq:B_q4_gp1` (3×8 matrix)
- [ ] D-matrix `eq:D_q4_ex` numerical values (76923, 23077, 26923)
- [ ] 8×8 k^e `eq:ke_q4_cook` (×10³ N/mm)
- [ ] 4×4 Kff `eq:Kff_q4` and solution values u₂,v₂,u₃,v₃
- [ ] Stress table `tab:q4_stresses` with all 4 GPs
- [ ] Figure `fig:q4_cook_results` loads from `figures/q4_cook_results.pdf`
- [ ] `lst:q4cook` in appendix resolves (no `??` in Step 6 and Step 8 references)
- [ ] Zero new undefined-reference warnings

- [ ] **Step 5: Commit and push**

```
git add chapters/ch04_2d_continuum.tex appendix/matlab_listings.tex
git commit -m "Ch.7: §7.3.4 Q4 Cook membrane -- appendix listing, chapter complete"
git push
```

---

## Self-Review

### Spec Coverage

| Spec requirement | Task |
|---|---|
| §7.3.4 subsection inside §7.3 | Task 2 |
| TikZ figure: tapered quad, clamped left, force arrows | Task 2 |
| Shape function natural derivatives table | Task 2 |
| Symbolic Jacobian J(ξ,η) derived | Task 2 |
| J₂₁=0 remark | Task 2 |
| det|J|=24(15−7ξ) depends on ξ only remark | Task 2 |
| J⁻¹ symbolic | Task 2 |
| Gauss-point table (4 GPs, J values, det|J|, J⁻¹) | Task 2 |
| B-matrix at GP1 (3×8, numerical) | Task 3 |
| Remark: B non-constant → Gauss integration essential | Task 3 |
| D-matrix (numerical) | Task 3 |
| k^e (8×8, Gauss sum, MATLAB-verified) | Task 3 |
| 4×4 Kff + RHS (F₄=F₆=50 N) | Task 4 |
| MATLAB-verified solution (u₂,v₂,u₃,v₃) | Task 4 |
| Deformed-shape figure `fig:q4_cook_results` | Task 4 |
| Stress table at 4 Gauss points | Task 4 |
| Remark: stress varies over element (vs. CST) | Task 4 |
| MATLAB script `fem_q4_cook.m` | Task 1 |
| Appendix listing `lst:q4cook` | Task 5 |
| All numerics MATLAB-verified | Task 1 → Tasks 2–4 |

### Placeholder Scan

All numeric values in Tasks 2–4 are drawn from the Pre-Verified Numbers section (MATLAB-confirmed). No TBD or TODO markers.

### Type Consistency

- `eq:jacobian_cook` defined Task 2, referenced in remark prose — consistent.
- `eq:Jinv_cook` defined Task 2, referenced in Task 3 Step 4 prose — consistent.
- `eq:B_q4_gp1` defined Task 3, label standalone — no forward-reference issues.
- `eq:ke_q4_cook` defined Task 3, referenced in Task 4 Step 7 (via Kff prose) — consistent.
- `eq:Kff_q4` defined Task 4 — no other references needed.
- `fig:q4_cook` defined Task 2, referenced in example opening — consistent.
- `fig:q4_cook_results` defined Task 4, filename matches `figures/q4_cook_results.pdf` saved by MATLAB in Task 1 — consistent.
- `lst:q4cook` defined Task 5, referenced by `\Cref{lst:q4cook}` in Tasks 3 and 4 — will produce `??` until Task 5 is complete; this is expected.
- Cross-refs `eq:q4_shape`, `eq:jacobian`, `eq:q4_stiffness`, `eq:D_planestress` all exist in the current file — verified before writing plan.
