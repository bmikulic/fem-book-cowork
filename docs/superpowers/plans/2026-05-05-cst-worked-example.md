# CST Two-Element Worked Example — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the 3-line §7.4 CST stub with a complete 9-step two-element worked example inserted as §7.2.4, plus a new self-contained MATLAB script and updated appendix listing.

**Architecture:** Run MATLAB first to generate the figure; write LaTeX sections in order; update the appendix listing last; build and verify. All numeric values in the LaTeX are taken from the MATLAB-verified output below.

**Tech Stack:** LaTeX (`latexmk -pdf`), MATLAB (R2016b+, local functions in scripts), TikZ, BibTeX.

---

## Pre-verified Numbers (MATLAB-confirmed)

These values are exact (the CST recovers the linear displacement field exactly for uniform tension):

**B-matrices** (`1/(2Ae) × [...]`, units m⁻¹):

B¹ (elem 1, nodes 1-2-3):
```
[-1,  0,  1,  0,  0,  0]
[ 0,  0,  0, -1,  0,  1]
[ 0, -1, -1,  1,  1,  0]
```

B² (elem 2, local nodes 1→global 1, 2→global 3, 3→global 4):
```
[ 0,  0,  1,  0, -1,  0]
[ 0, -1,  0,  0,  0,  1]
[-1,  0,  0,  1,  1, -1]
```

**D matrix** (plane stress, ×10¹¹ Pa):
```
[2.1978  0.6593  0     ]
[0.6593  2.1978  0     ]
[0       0       0.7692]
```

**k^e1** (×10⁸ N/m, DOFs u1,v1,u2,v2,u3,v3):
```
[ 10.9890   0       -10.9890   3.2967    0       -3.2967]
[  0        3.8462    3.8462  -3.8462   -3.8462   0     ]
[-10.9890   3.8462   14.8352  -7.1429   -3.8462   3.2967]
[  3.2967  -3.8462   -7.1429  14.8352    3.8462 -10.9890]
[  0       -3.8462   -3.8462   3.8462    3.8462   0     ]
[ -3.2967   0         3.2967 -10.9890    0       10.9890]
```

**k^e2** (×10⁸ N/m, DOFs u1,v1,u3,v3,u4,v4):
```
[  3.8462   0         0       -3.8462   -3.8462   3.8462]
[  0       10.9890   -3.2967   0         3.2967  -10.9890]
[  0       -3.2967   10.9890   0        -10.9890   3.2967]
[ -3.8462   0         0        3.8462    3.8462   -3.8462]
[ -3.8462   3.2967  -10.9890   3.8462   14.8352   -7.1429]
[  3.8462 -10.9890    3.2967  -3.8462   -7.1429   14.8352]
```

**K_ff** (×10⁸ N/m, free DOFs u2,v2,u3,v3,v4 = DOFs 3,4,5,6,8):
```
[ 14.8352  -7.1429  -3.8462   3.2967   0     ]
[ -7.1429  14.8352   3.8462 -10.9890   0     ]
[ -3.8462   3.8462  14.8352   0        3.2967]
[  3.2967 -10.9890   0       14.8352  -3.8462]
[  0        0        3.2967  -3.8462  14.8352]
```

**Solution** (u in metres, DOF order [u1,v1,u2,v2,u3,v3,u4,v4]):
```
u = [0, 0, 5.000e-4, 0, 5.000e-4, -1.500e-4, 0, -1.500e-4]
```

**Stresses** (both elements, in Pa):
```
σx = 1.000e+08 Pa = 100 MPa
σy = 0 Pa
τxy = 0 Pa
```

---

## File Map

| File | Action | Purpose |
|---|---|---|
| `matlab/fem_cst_2d.m` | Create (replaces old appendix listing) | Two-element CST script; generates figure |
| `chapters/ch04_2d_continuum.tex` | Modify | Remove §7.4 stub; add §7.2.4 subsection |
| `appendix/matlab_listings.tex` | Modify | Replace embedded `lstlisting` with `\lstinputlisting` |

---

## Task 1: Write and Run the MATLAB Script

**Files:**
- Create: `matlab/fem_cst_2d.m`

- [ ] **Step 1: Create `matlab/fem_cst_2d.m` with this exact content**

```matlab
% Chapter 7 -- Two-element CST mesh under uniaxial tension
% Plate: 1 m x 1 m x 0.01 m, plane stress, E=200 GPa, nu=0.3
% Nodes: 1(0,0) 2(1,0) 3(1,1) 4(0,1)
% Elements: 1=[1,2,3], 2=[1,3,4]
% BCs: pin at node 1 (u1=v1=0), roller at node 4 (u4=0)
% Load: F_x=500 kN at nodes 2 and 3
clear; clc; close all;

%% Parameters
E  = 200e9;   % Pa
nu = 0.3;
t  = 0.01;    % m

%% Nodes and connectivity
xy   = [0,0; 1,0; 1,1; 0,1];
conn = [1,2,3; 1,3,4];

%% Constitutive matrix (plane stress)
D = E/(1-nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];

%% Assemble global K
n_dof = 8;
K = zeros(n_dof);
B_all = cell(2,1);

for e = 1:2
    nd = conn(e,:);
    x  = xy(nd,1);  y = xy(nd,2);
    Ae = 0.5*abs((x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)));
    b  = [y(2)-y(3); y(3)-y(1); y(1)-y(2)];
    c  = [x(3)-x(2); x(1)-x(3); x(2)-x(1)];
    B  = 1/(2*Ae) * [b(1),0,b(2),0,b(3),0; ...
                     0,c(1),0,c(2),0,c(3); ...
                     c(1),b(1),c(2),b(2),c(3),b(3)];
    ke = t * Ae * B' * D * B;
    B_all{e} = B;
    dofs = [2*nd(1)-1, 2*nd(1), 2*nd(2)-1, 2*nd(2), 2*nd(3)-1, 2*nd(3)];
    K(dofs,dofs) = K(dofs,dofs) + ke;
end

%% Load vector and BCs
F = zeros(n_dof,1);
F(3) = 500e3;   % u2
F(5) = 500e3;   % u3

fixed = [1, 2, 7];   % u1, v1, u4
free  = setdiff(1:n_dof, fixed);
u = zeros(n_dof,1);
u(free) = K(free,free) \ F(free);

%% Print results
fprintf('=== NODAL DISPLACEMENTS ===\n')
for n = 1:4
    fprintf('Node %d: u=%.4e m  v=%.4e m\n', n, u(2*n-1), u(2*n))
end

fprintf('\n=== ELEMENT STRESSES ===\n')
for e = 1:2
    nd   = conn(e,:);
    dofs = [2*nd(1)-1,2*nd(1),2*nd(2)-1,2*nd(2),2*nd(3)-1,2*nd(3)];
    sig  = D * B_all{e} * u(dofs);
    fprintf('Elem %d: sx=%.4e Pa  sy=%.4e Pa  txy=%.4e Pa\n', ...
            e, sig(1), sig(2), sig(3))
end

%% Figure: mesh + deformed shape
scale = 800;
fig = figure('Units','centimeters','Position',[2,2,14,10]);
hold on; axis equal; grid on;

tri_x = xy([1,2,3,1],1); tri_y = xy([1,2,3,1],2);
plot(tri_x, tri_y, 'k-', 'LineWidth', 1.2)
tri_x = xy([1,3,4,1],1); tri_y = xy([1,3,4,1],2);
plot(tri_x, tri_y, 'k-', 'LineWidth', 1.2)

% Deformed nodes
xy_def = xy + scale * reshape(u,2,[])';
tri_x = xy_def([1,2,3,1],1); tri_y = xy_def([1,2,3,1],2);
plot(tri_x, tri_y, 'b--', 'LineWidth', 1.5)
tri_x = xy_def([1,3,4,1],1); tri_y = xy_def([1,3,4,1],2);
plot(tri_x, tri_y, 'b--', 'LineWidth', 1.5)

% Node labels
for n = 1:4
    text(xy(n,1)-0.07, xy(n,2)+0.04, sprintf('  %d',n), 'FontSize',10)
end

% Element labels
c1 = mean(xy([1,2,3],:)); text(c1(1),c1(2),'\textcircled{1}','Interpreter','latex','FontSize',12,'HorizontalAlignment','center')
c2 = mean(xy([1,3,4],:)); text(c2(1),c2(2),'\textcircled{2}','Interpreter','latex','FontSize',12,'HorizontalAlignment','center')

xlabel('x (m)'); ylabel('y (m)')
title(sprintf('CST mesh (solid) and deformed shape \\times%d (dashed)', scale))
legend({'Undeformed','','Deformed (\\times800)',''}, 'Location','northwest')

exportgraphics(fig,'../figures/cst_example_results.pdf','ContentType','vector')
exportgraphics(fig,'../figures/cst_example_results.png','Resolution',300)
fprintf('\nFigures saved to ../figures/\n')
```

- [ ] **Step 2: Run the script in MATLAB**

Open MATLAB, navigate to `matlab/`, run:
```
fem_cst_2d
```

- [ ] **Step 3: Confirm output matches these values**

Expected terminal output:
```
=== NODAL DISPLACEMENTS ===
Node 1: u= 0.0000e+00 m  v= 0.0000e+00 m
Node 2: u= 5.0000e-04 m  v= 0.0000e+00 m
Node 3: u= 5.0000e-04 m  v=-1.5000e-04 m
Node 4: u= 0.0000e+00 m  v=-1.5000e-04 m

=== ELEMENT STRESSES ===
Elem 1: sx= 1.0000e+08 Pa  sy~0  txy~0
Elem 2: sx= 1.0000e+08 Pa  sy~0  txy~0
```

- [ ] **Step 4: Commit**

```
git add matlab/fem_cst_2d.m
git commit -m "Ch.7: add two-element CST MATLAB script"
```

---

## Task 2: Remove §7.4 Stub and Add §7.2.4 — Problem Setup and Figure

**Files:**
- Modify: `chapters/ch04_2d_continuum.tex`

The file currently ends at line 200 with `\end{example}` closing the §7.4 stub. The stub spans lines 185–200.

- [ ] **Step 1: Delete the §7.4 stub section (lines 185–200)**

Remove from `chapters/ch04_2d_continuum.tex`:
```latex
\section{Worked Example: CST Plate in Tension}

\begin{example}[Single CST element under tension]
\label{ex:cst_tension}

A single CST element: nodes at $(0,0)$, $(1,0)$, $(0,1)$, thickness $t=0.01$ m,
$E = 200\times10^9$ Pa, $\nu = 0.3$, plane stress. Apply $\sigma_x = 100$ MPa
(force $F = 0.5$ MN at nodes 2 and 3).

\medskip
\noindent $A_e = 0.5$ m$^2$,\quad
$\matr{B}$: constant from \cref{eq:B_cst},\quad
$\matr{D}$: from \cref{eq:D_planestress}.

\noindent See \Cref{lst:cst2d} in the Appendix for the MATLAB implementation.
\end{example}
```

- [ ] **Step 2: Append §7.2.4 inside the `\section{Constant Strain Triangle (CST)}` block**

Find the line `\section{Isoparametric Quadrilateral Element (Q4)}` (currently line 123). Insert the following block immediately before it:

```latex

% =============================================================================
\subsection{Worked Example: Two-Element CST Mesh in Uniaxial Tension}
\label{sec:cst_example}
% =============================================================================

\begin{example}[Two-element CST plate under uniaxial tension]
\label{ex:cst_tension}

A $1\ \text{m} \times 1\ \text{m}$ square plate of thickness $t = 0.01\ \text{m}$
is modelled with two CST elements formed by splitting the square along its diagonal,
as shown in \cref{fig:cst_example}. Material properties: $E = 200\ \text{GPa}$,
$\nu = 0.3$, plane stress. A uniform tensile stress $\sigma_x = 100\ \text{MPa}$
is applied to the right edge, equivalent to a total force
$F = \sigma_x \cdot t \cdot 1\ \text{m} = 1\ \text{MN}$
split equally between nodes~2 and~3 as $F/2 = 500\ \text{kN}$ each.
The left-edge nodes are restrained: pin at node~1 ($u_1 = v_1 = 0$) and
horizontal roller at node~4 ($u_4 = 0$).

\begin{figure}[htbp]
\centering
\begin{tikzpicture}[>=Stealth, x=3.2cm, y=3.2cm,
    force/.style={->, ultra thick, red!75!black},
    disp/.style={->, thick, blue!70!black}]

  %% Square outline and diagonal
  \fill[gray!12] (0,0) -- (1,0) -- (1,1) -- (0,1) -- cycle;
  \draw[thick] (0,0) -- (1,0) -- (1,1) -- (0,1) -- cycle;
  \draw[thick] (0,0) -- (1,1);

  %% Element labels
  \node at (0.62,0.25) {\small\textcircled{\small 1}};
  \node at (0.25,0.70) {\small\textcircled{\small 2}};

  %% Nodes
  \filldraw (0,0) circle (2.5pt) node[below left=2pt] {\small\bfseries 1 $(0,0)$};
  \filldraw (1,0) circle (2.5pt) node[below right=2pt] {\small\bfseries 2 $(1,0)$};
  \filldraw (1,1) circle (2.5pt) node[above right=2pt] {\small\bfseries 3 $(1,1)$};
  \filldraw (0,1) circle (2.5pt) node[above left=2pt] {\small\bfseries 4 $(0,1)$};

  %% Pin support at node 1
  \draw[fill=gray!30,thick] (0,0) -- (-0.12,-0.18) -- (0.12,-0.18) -- cycle;
  \fill[pattern=north east lines,pattern color=gray!70]
       (-0.14,-0.18) rectangle (0.14,-0.25);
  \draw[thick] (-0.14,-0.18) -- (0.14,-0.18);

  %% Roller at node 4 (horizontal roller: free v, fixed u)
  \draw[fill=gray!30,thick] (0,1) -- (-0.12,0.82) -- (-0.12,1.18) -- cycle;
  \fill[pattern=north east lines,pattern color=gray!70]
       (-0.25,0.82) rectangle (-0.18,1.18);
  \draw[thick] (-0.18,0.82) -- (-0.18,1.18);

  %% Applied forces at nodes 2 and 3
  \draw[force] (1.0,0) -- (1.45,0)
       node[right=2pt] {$F/2 = 500\ \text{kN}$};
  \draw[force] (1.0,1) -- (1.45,1)
       node[right=2pt] {$F/2 = 500\ \text{kN}$};

  %% Dimensions
  \draw[<->,thin] (0,-0.35) -- (1,-0.35)
       node[fill=white,inner sep=1pt,midway] {$1\ \text{m}$};
  \draw[<->,thin] (-0.35,0) -- (-0.35,1)
       node[fill=white,inner sep=1pt,midway,rotate=90] {$1\ \text{m}$};

\end{tikzpicture}
\caption{Two-element CST mesh: 4 nodes, 2 triangular elements formed by the
diagonal from node~1 to node~3. Pin support at node~1, horizontal roller at
node~4, and $500\ \text{kN}$ horizontal force at nodes~2 and~3.}
\label{fig:cst_example}
\end{figure}

\medskip
\noindent\textbf{DOF ordering:}
$[u_1,v_1,u_2,v_2,u_3,v_3,u_4,v_4] = \text{DOFs}\ [1\ldots8]$.

\begin{table}[htbp]
\centering
\small
\caption{Node coordinates and element connectivity.}
\label{tab:cst_elements}
\begin{tabular}{@{}cccccc@{}}
\toprule
Elem & Nodes & Node coords (m) & $A_e$ (m$^2$) & Global DOFs \\
\midrule
1 & $1\to2\to3$ & $(0,0),(1,0),(1,1)$ & $0.5$ & 1,2,3,4,5,6 \\
2 & $1\to3\to4$ & $(0,0),(1,1),(0,1)$ & $0.5$ & 1,2,5,6,7,8 \\
\bottomrule
\end{tabular}
\end{table}
```

- [ ] **Step 3: Compile and check**

```
cd "C:\Users\Korisnik\Documents\FEM_Book Cowork"
latexmk -pdf main.tex
```

Expected: compiles; §7.2.4 subsection appears with TikZ figure and table. The old §7.4 heading is gone. Some forward-reference warnings (`lst:cst2d`) are expected.

- [ ] **Step 4: Commit**

```
git add chapters/ch04_2d_continuum.tex
git commit -m "Ch.7: remove §7.4 stub, add §7.2.4 problem setup and TikZ figure"
```

---

## Task 3: Write Steps 2–4 — B-matrices, D-matrix, Element Stiffness Matrices

**Files:**
- Modify: `chapters/ch04_2d_continuum.tex` (append inside the example, before `\end{example}`)

Find the `\end{example}` that closes the example added in Task 2 and insert the following before it:

- [ ] **Step 1: Append Steps 2–4 inside the example**

```latex

\medskip
\noindent\textbf{Step 1: Element geometry and B-matrices.}

For a CST with nodes at $(x_1,y_1)$, $(x_2,y_2)$, $(x_3,y_3)$, the
coefficients are $b_i = y_j - y_k$ and $c_i = x_k - x_j$
(cyclic, \cref{eq:cst_shape}), and $\matr{B}$ is constant (\cref{eq:B_cst}).

\smallskip
\noindent\textit{Element~1} (nodes $1\to2\to3$, coordinates $(0,0),(1,0),(1,1)$):
\[
    b_1 = -1,\quad b_2 = 1,\quad b_3 = 0,\qquad
    c_1 = 0,\quad c_2 = -1,\quad c_3 = 1,\qquad 2A_e = 1\ \text{m}^2
\]
\begin{equation}
    \matr{B}^{(1)} =
    \begin{bmatrix}
        -1 &  0 &  1 &  0 &  0 &  0 \\
         0 &  0 &  0 & -1 &  0 &  1 \\
         0 & -1 & -1 &  1 &  1 &  0
    \end{bmatrix}\ \text{m}^{-1}
    \label{eq:B1_cst_ex}
\end{equation}

\noindent\textit{Element~2} (nodes $1\to3\to4$, coordinates $(0,0),(1,1),(0,1)$):
\[
    b_1 = 0,\quad b_2 = 1,\quad b_3 = -1,\qquad
    c_1 = -1,\quad c_2 = 0,\quad c_3 = 1,\qquad 2A_e = 1\ \text{m}^2
\]
\begin{equation}
    \matr{B}^{(2)} =
    \begin{bmatrix}
         0 &  0 &  1 &  0 & -1 &  0 \\
         0 & -1 &  0 &  0 &  0 &  1 \\
        -1 &  0 &  0 &  1 &  1 & -1
    \end{bmatrix}\ \text{m}^{-1}
    \label{eq:B2_cst_ex}
\end{equation}

The two B-matrices differ because the diagonal splits the square asymmetrically:
element~1 has its right-angle corner at node~2, element~2 at node~4.

\medskip
\noindent\textbf{Step 2: Constitutive matrix.}

Plane stress with $E = 200\ \text{GPa}$, $\nu = 0.3$:
\begin{equation}
    \matr{D} = \frac{200\times10^9}{1-0.09}
    \begin{bmatrix}1 & 0.3 & 0 \\ 0.3 & 1 & 0 \\ 0 & 0 & 0.35\end{bmatrix}
    =
    \begin{bmatrix}
        2.198\times10^{11} & 6.593\times10^{10} & 0 \\
        6.593\times10^{10} & 2.198\times10^{11} & 0 \\
        0 & 0 & 7.692\times10^{10}
    \end{bmatrix}\ \text{Pa}
    \label{eq:D_cst_ex}
\end{equation}

\medskip
\noindent\textbf{Step 3: Element stiffness matrices.}

Using $\matr{k}^e = t\,A_e\,\matr{B}\transpose\matr{D}\matr{B}$ (\cref{eq:cst_stiffness})
with $t\,A_e = 0.01 \times 0.5 = 0.005\ \text{m}^3$:

\begin{equation}
    \matr{k}^{(1)} = 10^8 \times
    \begin{bmatrix}
        10.99 &  0     & -10.99 &  3.30  &  0     & -3.30  \\
         0    &  3.85  &   3.85 & -3.85  & -3.85  &  0     \\
       -10.99 &  3.85  &  14.84 & -7.14  & -3.85  &  3.30  \\
         3.30 & -3.85  &  -7.14 & 14.84  &  3.85  & -10.99 \\
         0    & -3.85  &  -3.85 &  3.85  &  3.85  &  0     \\
        -3.30 &  0     &   3.30 & -10.99 &  0     & 10.99
    \end{bmatrix}\ \tfrac{\text{N}}{\text{m}}
    \quad \text{(DOFs }u_1,v_1,u_2,v_2,u_3,v_3\text{)}
    \label{eq:ke1_cst}
\end{equation}

\begin{equation}
    \matr{k}^{(2)} = 10^8 \times
    \begin{bmatrix}
         3.85 &  0     &  0     & -3.85  & -3.85  &  3.85  \\
         0    & 10.99  & -3.30  &  0     &  3.30  & -10.99 \\
         0    & -3.30  & 10.99  &  0     & -10.99 &  3.30  \\
        -3.85 &  0     &  0     &  3.85  &  3.85  & -3.85  \\
        -3.85 &  3.30  & -10.99 &  3.85  & 14.84  & -7.14  \\
         3.85 & -10.99 &  3.30  & -3.85  & -7.14  & 14.84
    \end{bmatrix}\ \tfrac{\text{N}}{\text{m}}
    \quad \text{(DOFs }u_1,v_1,u_3,v_3,u_4,v_4\text{)}
    \label{eq:ke2_cst}
\end{equation}

\begin{remark}
$\matr{k}^{(1)} \neq \matr{k}^{(2)}$ even though both elements have the same area
and material. The two triangles are geometrically distinct (different orientations of
the hypotenuse) so their B-matrices differ, leading to different stiffness matrices.
\end{remark}
```

- [ ] **Step 2: Compile and check**

```
latexmk -pdf main.tex
```

Expected: Steps 1–3 appear in the example. Equations `eq:B1_cst_ex`, `eq:B2_cst_ex`, `eq:D_cst_ex`, `eq:ke1_cst`, `eq:ke2_cst` defined.

- [ ] **Step 3: Commit**

```
git add chapters/ch04_2d_continuum.tex
git commit -m "Ch.7: §7.2.4 Steps 1-3 -- B-matrices, D-matrix, element stiffness"
```

---

## Task 4: Write Steps 4–7 — Global Assembly, BCs, and MATLAB-Verified Solution

**Files:**
- Modify: `chapters/ch04_2d_continuum.tex` (append inside the example, before `\end{example}`)

- [ ] **Step 1: Append Steps 4–7 inside the example**

```latex

\medskip
\noindent\textbf{Step 4: Global assembly.}

The $8\times8$ global stiffness $\stiff$ is assembled from the two element
matrices. Element~1 contributes to DOFs $\{1,2,3,4,5,6\}$; element~2 to
DOFs $\{1,2,5,6,7,8\}$. DOFs $\{1,2,5,6\}$ (nodes~1 and~3) receive
contributions from both elements:

\begin{equation}
    \stiff = 10^8 \times
    \begin{bmatrix*}[r]
        14.84 &  0     & -10.99 &  3.30  &  0     & -7.14  & -3.85  &  3.85  \\
         0    & 14.84  &  3.85  & -3.85  & -7.14  &  0     &  3.30  & -10.99 \\
       -10.99 &  3.85  & 14.84  & -7.14  & -3.85  &  3.30  &  0     &  0     \\
         3.30 & -3.85  & -7.14  & 14.84  &  3.85  & -10.99 &  0     &  0     \\
         0    & -7.14  & -3.85  &  3.85  & 14.84  &  0     & -10.99 &  3.30  \\
        -7.14 &  0     &  3.30  & -10.99 &  0     & 14.84  &  3.85  & -3.85  \\
        -3.85 &  3.30  &  0     &  0     & -10.99 &  3.85  & 14.84  & -7.14  \\
         3.85 & -10.99 &  0     &  0     &  3.30  & -3.85  & -7.14  & 14.84
    \end{bmatrix*}\ \tfrac{\text{N}}{\text{m}}
    \label{eq:K_cst_global}
\end{equation}

Note that all diagonal entries equal $14.84 \times 10^8\ \text{N/m}$ --- a
consequence of the mesh symmetry.

\medskip
\noindent\textbf{Step 5: Load vector and boundary conditions.}

The global load vector has two non-zero entries:
$F_3 = F_5 = 500\,000\ \text{N}$ (forces at DOFs~3 and~5, i.e.\ $u_2$ and $u_3$).

Essential BCs: $u_1 = 0$ (DOF~1), $v_1 = 0$ (DOF~2), $u_4 = 0$ (DOF~7).
The five free DOFs are $\{3,4,5,6,8\}$ ($u_2, v_2, u_3, v_3, v_4$). The
reduced $5\times5$ system is:

\begin{equation}
    10^8 \times
    \begin{bmatrix}
        14.84 & -7.14  & -3.85  &  3.30  &  0     \\
        -7.14 & 14.84  &  3.85  & -10.99 &  0     \\
        -3.85 &  3.85  & 14.84  &  0     &  3.30  \\
         3.30 & -10.99 &  0     & 14.84  & -3.85  \\
         0    &  0     &  3.30  & -3.85  & 14.84
    \end{bmatrix}
    \begin{bmatrix} u_2 \\ v_2 \\ u_3 \\ v_3 \\ v_4 \end{bmatrix}
    =
    \begin{bmatrix} 500\,000 \\ 0 \\ 500\,000 \\ 0 \\ 0 \end{bmatrix}\ \text{N}
    \label{eq:Kff_cst}
\end{equation}

\medskip
\noindent\textbf{Step 6: Solution (MATLAB-verified, \Cref{lst:cst2d}).}

\begin{align*}
    u_2 &= +5.000\times10^{-4}\ \text{m}, &
    v_2 &= 0\ \text{m}, \\
    u_3 &= +5.000\times10^{-4}\ \text{m}, &
    v_3 &= -1.500\times10^{-4}\ \text{m}, \\
    u_4 &= 0\ \text{m}, &
    v_4 &= -1.500\times10^{-4}\ \text{m}.
\end{align*}

The right edge moves uniformly rightward ($u_2 = u_3 = 0.5\ \text{mm}$) with
no shear distortion. The top-right and top-left nodes contract vertically by
$0.15\ \text{mm}$ due to the Poisson effect.
```

- [ ] **Step 2: Compile and check**

```
latexmk -pdf main.tex
```

Expected: Steps 4–6 appear; the 8×8 and 5×5 matrices render; all numeric values present.

- [ ] **Step 3: Commit**

```
git add chapters/ch04_2d_continuum.tex
git commit -m "Ch.7: §7.2.4 Steps 4-6 -- global assembly, BCs, and solution"
```

---

## Task 5: Write Steps 7–9 — Stress Recovery, Verification, and Results Figure

**Files:**
- Modify: `chapters/ch04_2d_continuum.tex` (append inside the example, before `\end{example}`)

- [ ] **Step 1: Append Steps 7–9 and close the example**

```latex

\medskip
\noindent\textbf{Step 7: Stress recovery.}

For each element, the constant strain and stress vectors are:
\[
    \vect{\varepsilon}^e = \matr{B}^e\,\vect{u}^e, \qquad
    \vect{\sigma}^e = \matr{D}\,\vect{\varepsilon}^e .
\]

\noindent\textit{Element~1} ($\vect{u}^{(1)} = [0,0,5\!\times\!10^{-4},0,5\!\times\!10^{-4},-1.5\!\times\!10^{-4}]\transpose$ m):
\begin{align*}
    \vect{\varepsilon}^{(1)} &= \matr{B}^{(1)}\vect{u}^{(1)}
    = \begin{bmatrix}5\times10^{-4} \\ -1.5\times10^{-4} \\ 0\end{bmatrix}
    \begin{bmatrix}\text{m/m}\\\text{m/m}\\\text{rad}\end{bmatrix} \\[4pt]
    \vect{\sigma}^{(1)} &= \matr{D}\,\vect{\varepsilon}^{(1)}
    = \begin{bmatrix}100 \\ 0 \\ 0\end{bmatrix}\ \text{MPa}
\end{align*}

\noindent\textit{Element~2} ($\vect{u}^{(2)} = [0,0,5\!\times\!10^{-4},-1.5\!\times\!10^{-4},0,-1.5\!\times\!10^{-4}]\transpose$ m):
\[
    \vect{\varepsilon}^{(2)} = \vect{\varepsilon}^{(1)}, \qquad
    \vect{\sigma}^{(2)} = \vect{\sigma}^{(1)} =
    \begin{bmatrix}100 \\ 0 \\ 0\end{bmatrix}\ \text{MPa}
\]

Both elements give identical stresses: the CST mesh is in a state of uniform
$\sigma_x = 100\ \text{MPa}$, $\sigma_y = 0$, $\tau_{xy} = 0$.

\medskip
\noindent\textbf{Step 8: Analytical verification.}

The exact solution for a uniformly loaded elastic plate is:
\begin{align*}
    \sigma_x &= \frac{F}{t \cdot L} = \frac{1\times10^6}{0.01 \times 1}
              = 100\ \text{MPa}\ \checkmark \\
    \sigma_y &= 0\ \checkmark, \qquad \tau_{xy} = 0\ \checkmark \\[4pt]
    u(x,y) &= \frac{\sigma_x}{E}\,x = 5\times10^{-4}\,x\ \text{m}
              \;\Rightarrow\; u_2 = u_3 = 5\times10^{-4}\ \text{m}\ \checkmark \\
    v(x,y) &= -\frac{\nu\sigma_x}{E}\,y = -1.5\times10^{-4}\,y\ \text{m}
              \;\Rightarrow\; v_3 = v_4 = -1.5\times10^{-4}\ \text{m}\ \checkmark
\end{align*}

\begin{importantbox}
The CST recovers the \emph{exact} analytical solution for uniform stress
states. This is because the true displacement field ($u$ linear in $x$,
$v$ linear in $y$) lies within the CST polynomial space — a polynomial
that the element can represent exactly. The FEM error is zero for any
problem whose exact solution is linear.
\end{importantbox}

The deformed mesh, magnified $\times800$, is shown in \cref{fig:cst_results}.

\begin{figure}[htbp]
\centering
\includegraphics[width=0.62\textwidth]{figures/cst_example_results}
\caption{CST mesh (solid) and deformed shape (dashed, magnified $\times800$).
The right edge moves uniformly rightward; the top nodes contract downward by
the Poisson effect. No shear distortion occurs.}
\label{fig:cst_results}
\end{figure}

\end{example}
```

- [ ] **Step 2: Compile and check**

```
latexmk -pdf main.tex
```

Expected: Steps 7–9 appear with stress recovery, verification equations with checkmarks, importantbox, and the MATLAB figure. Warning about `lst:cst2d` may still appear — that resolves in Task 6.

- [ ] **Step 3: Commit**

```
git add chapters/ch04_2d_continuum.tex
git commit -m "Ch.7: §7.2.4 Steps 7-9 -- stress recovery, verification, exact-solution importantbox"
```

---

## Task 6: Update Appendix Listing, Final Build and Push

**Files:**
- Modify: `appendix/matlab_listings.tex`

The current appendix section at lines 446–494 has:
```
\section{CST Element for Plane Stress}
\label{lst:cst2d}
\begin{lstlisting}[caption={...}, label={code:cst2d}]
... old embedded code ...
\end{lstlisting}
```

- [ ] **Step 1: Replace the embedded lstlisting with lstinputlisting**

In `appendix/matlab_listings.tex`, find the block from `\section{CST Element for Plane Stress}` through the closing `\end{lstlisting}` (~lines 446–494). Replace the entire block with:

```latex
% ---------------------------------------------------------------
\section{Chapter~7 --- Two-Element CST Mesh}
\label{lst:cst2d}

\lstinputlisting[
    caption={Two-element CST plate in uniaxial tension
             (\Cref{ex:cst_tension}).},
    label={code:cst2d}
]{matlab/fem_cst_2d.m}
```

- [ ] **Step 2: Final build**

```
cd "C:\Users\Korisnik\Documents\FEM_Book Cowork"
latexmk -pdf main.tex
```

- [ ] **Step 3: Verify the PDF**

Check the compiled PDF:
- [ ] §7.2.4 subsection appears inside §7.2 in the TOC (before §7.3 Q4)
- [ ] §7.4 heading "Worked Example: CST Plate in Tension" is gone
- [ ] TikZ figure `fig:cst_example` renders (4 nodes, diagonal, pin, roller, force arrows)
- [ ] All three B-matrices written as 3×6 arrays
- [ ] Both 6×6 k^e matrices present (labeled `eq:ke1_cst`, `eq:ke2_cst`)
- [ ] 8×8 global K and 5×5 Kff appear
- [ ] Solution: u2=u3=5e-4 m, v3=v4=-1.5e-4 m
- [ ] `importantbox` with "CST recovers exact solution" present
- [ ] MATLAB figure `fig:cst_results` loads from `figures/cst_example_results.pdf`
- [ ] Appendix section title updated to "Chapter 7 -- Two-Element CST Mesh"
- [ ] `lst:cst2d` resolves (no `??` for this reference in the chapter)
- [ ] No new undefined-reference warnings in Chapter 7

- [ ] **Step 4: Commit and push**

```
git add chapters/ch04_2d_continuum.tex appendix/matlab_listings.tex
git commit -m "Ch.7: §7.2.4 complete -- appendix listing updated, chapter done"
git push
```

---

## Self-Review

### Spec Coverage

| Spec requirement | Task |
|---|---|
| Remove §7.4 stub | Task 2 Step 1 |
| Add §7.2.4 subsection within §7.2 | Task 2 Step 2 |
| TikZ figure: nodes, elements, BCs, forces | Task 2 Step 2 |
| Node/element connectivity table | Task 2 Step 2 |
| B¹ and B² matrices (3×6, numerically) | Task 3 |
| D-matrix (plane stress, numerically) | Task 3 |
| k^e1 and k^e2 (6×6, numerically) | Task 3 |
| Remark: k^e1 ≠ k^e2 explanation | Task 3 |
| 8×8 global K assembly | Task 4 |
| Load vector and 5×5 Kff | Task 4 |
| MATLAB-verified solution (all 8 DOFs) | Task 4 |
| Physical interpretation of displacements | Task 4 |
| Stress recovery: ε = B u, σ = D ε | Task 5 |
| Stress recovery for both elements | Task 5 |
| Analytical verification with checkmarks | Task 5 |
| `importantbox` on exact CST for linear fields | Task 5 |
| MATLAB figure `fig:cst_results` | Task 5 |
| MATLAB script `fem_cst_2d.m` | Task 1 |
| Appendix updated to lstinputlisting | Task 6 |
| label `ex:cst_tension` reused | Task 2 |
| label `lst:cst2d` preserved | Task 6 |

### Placeholder Scan

All numeric values in Tasks 3–5 are drawn from the Pre-verified Numbers section (MATLAB-confirmed). No TBD or TODO markers.

### Type Consistency

- `fig:cst_example` defined Task 2, referenced Task 5 — consistent.
- `fig:cst_results` defined Task 5, matches `figures/cst_example_results.pdf` saved by MATLAB in Task 1 — consistent.
- `eq:ke1_cst`, `eq:ke2_cst` defined Task 3 — not cross-referenced elsewhere (standalone example equations).
- `\Cref{lst:cst2d}` in chapter text references section label `lst:cst2d` defined in appendix Task 6 — consistent.
- `ex:cst_tension` reuses old label — verify no other file cross-references it before Task 2 deletes the old definition (a quick grep shows the only reference is in the appendix caption, which is updated in Task 6).
