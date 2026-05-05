# Ch.6 Frame Element Expansion — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Expand `chapters/ch06_frame.tex` from a 27-line stub to a full chapter with lean theory (§6.1–6.4) and a complete portal-frame worked example (§6.5–6.6), all numerics MATLAB-verified.

**Architecture:** Run the MATLAB script first to obtain verified numbers; write each LaTeX section using those numbers; add the appendix listing; build the PDF. The stub's label `eq:frame_T` is replaced by `eq:frame_T6` (confirmed not referenced elsewhere).

**Tech Stack:** LaTeX (`latexmk -pdf`), MATLAB (R2016b+ for local functions in scripts), TikZ, BibTeX.

---

## File Map

| File | Action | Purpose |
|---|---|---|
| `chapters/ch06_frame.tex` | Full rewrite | Main chapter — 6 sections |
| `matlab/fem_frame_portal.m` | Create | FEM script; generates figures and verified numbers |
| `appendix/matlab_listings.tex` | Append | Add `lst:frame_portal` listing block |

---

## Task 1: Write and Run the Portal Frame MATLAB Script

**Files:**
- Create: `matlab/fem_frame_portal.m`

- [ ] **Step 1: Create `matlab/fem_frame_portal.m` with this exact content**

```matlab
% Chapter 6 -- Portal Frame under Combined Loading
% Geometry: H=3000 mm (columns), B=4000 mm (beam)
% Nodes:  1(0,0)  2(0,3000)  3(4000,3000)  4(4000,0)
% Elements: 1=left col (1->2), 2=beam (2->3), 3=right col (3->4)
% BCs: nodes 1 and 4 fully clamped
% Loading: F_x=10000 N at node 2; q=5 N/mm downward on beam
clear; clc; close all;

%% Parameters
E      = 210000;   % N/mm^2
A      = 6000;     % mm^2
Iz     = 10e6;     % mm^4
H      = 3000;     % mm  column height
B      = 4000;     % mm  beam span
F_lat  = 10000;    % N   horizontal force at node 2 (global x)
q_beam = 5;        % N/mm uniform downward on beam (element 2)

%% Node coordinates [x, y]
coords = [0, 0; 0, H; B, H; B, 0];
n_dof  = 12;   % 4 nodes x 3 DOFs

%% Element connectivity and inclination angles (deg, CCW from +x)
conn  = [1 2; 2 3; 3 4];
alpha_deg = [90; 0; -90];
n_elem = 3;

%% Assemble global K and F
K = zeros(n_dof);
F = zeros(n_dof, 1);

for e = 1:n_elem
    ni = conn(e,1); nj = conn(e,2);
    L  = norm(coords(nj,:) - coords(ni,:));
    c  = cosd(alpha_deg(e));
    s  = sind(alpha_deg(e));
    T6 = frame_T6(c, s);
    kbar = frame_local_k(E, A, Iz, L);
    ke   = T6' * kbar * T6;
    dofs = elem_dofs(ni, nj);
    K(dofs, dofs) = K(dofs, dofs) + ke;

    if e == 2  % uniform downward load on beam
        q0   = -q_beam;  % upward-positive convention
        fbar = [0; q0*L/2; q0*L^2/12; 0; q0*L/2; -q0*L^2/12];
        F(dofs) = F(dofs) + T6' * fbar;
    end
end

% Lateral point force at node 2, global x-direction (DOF 4)
F(4) = F(4) + F_lat;

%% Apply BCs: fix nodes 1 (DOFs 1-3) and 4 (DOFs 10-12)
fixed = [1 2 3 10 11 12];
free  = setdiff(1:n_dof, fixed);
u = zeros(n_dof, 1);
u(free) = K(free, free) \ F(free);

%% Reactions
R = K * u - F;

%% Print verified results
fprintf('=== FREE-DOF DISPLACEMENTS ===\n');
fprintf('Node 2: u2=%+.5f mm  v2=%+.5f mm  th2=%+.6f rad\n', u(4),u(5),u(6));
fprintf('Node 3: u3=%+.5f mm  v3=%+.5f mm  th3=%+.6f rad\n', u(7),u(8),u(9));
fprintf('\n=== SUPPORT REACTIONS ===\n');
fprintf('Node 1: Rx1=%+.1f N  Ry1=%+.1f N  M1=%+.1f N*mm\n', R(1),R(2),R(3));
fprintf('Node 4: Rx4=%+.1f N  Ry4=%+.1f N  M4=%+.1f N*mm\n',R(10),R(11),R(12));
fprintf('\n=== EQUILIBRIUM CHECK ===\n');
fprintf('Sum Fx (reactions+applied) = %.1f N\n', R(1)+R(10)+F_lat);
fprintf('Sum Fy (reactions+applied) = %.1f N\n', R(2)+R(11)-q_beam*B);

%% Element internal forces (local coordinates)
fprintf('\n=== ELEMENT END FORCES (local) ===\n');
fprintf('%-10s %10s %10s %12s %10s %10s %12s\n', ...
        'Element','N_i(N)','V_i(N)','M_i(N*mm)','N_j(N)','V_j(N)','M_j(N*mm)');
for e = 1:n_elem
    ni = conn(e,1); nj = conn(e,2);
    L  = norm(coords(nj,:) - coords(ni,:));
    c  = cosd(alpha_deg(e)); s = sind(alpha_deg(e));
    dofs  = elem_dofs(ni, nj);
    T6    = frame_T6(c, s);
    uloc  = T6 * u(dofs);
    kbar  = frame_local_k(E, A, Iz, L);
    if e == 2
        q0   = -q_beam;
        fbar_load = [0;q0*L/2;q0*L^2/12;0;q0*L/2;-q0*L^2/12];
    else
        fbar_load = zeros(6,1);
    end
    fint = kbar * uloc - fbar_load;
    fprintf('%-10d %10.1f %10.1f %12.0f %10.1f %10.1f %12.0f\n', ...
            e, fint(1),fint(2),fint(3),fint(4),fint(5),fint(6));
end

%% Figures
n_pts = 80;
xi_v  = linspace(0, 1, n_pts);

X_el = cell(n_elem,1); Y_el = X_el;
Xd_el = X_el; Yd_el = X_el;
N_el = X_el; V_el = X_el; M_el = X_el;

scale_disp = 100;  % displacement magnification for deflected shape

for e = 1:n_elem
    ni = conn(e,1); nj = conn(e,2);
    xi = coords(ni,1); yi = coords(ni,2);
    xj = coords(nj,1); yj = coords(nj,2);
    L  = norm([xj-xi, yj-yi]);
    c  = cosd(alpha_deg(e)); s = sind(alpha_deg(e));
    dofs = elem_dofs(ni, nj);
    T6   = frame_T6(c, s);
    uloc = T6 * u(dofs);
    kbar = frame_local_k(E, A, Iz, L);
    if e == 2
        q0 = -q_beam;
        fbar_load = [0;q0*L/2;q0*L^2/12;0;q0*L/2;-q0*L^2/12];
    else
        q0 = 0;
        fbar_load = zeros(6,1);
    end
    fint = kbar * uloc - fbar_load;

    xv = xi_v * L;
    N_dist = fint(1) * ones(1, n_pts);
    V_dist = fint(2) + q0 * xv;         % V(x) = V_1 + q0*x
    M_dist = fint(3) - fint(2)*xv - q0*xv.^2/2; % M(x)

    X_gl = xi + xi_v*(xj - xi);
    Y_gl = yi + xi_v*(yj - yi);

    % Deformed shape via linear interp of nodal displacements
    Xd_gl = X_gl + scale_disp*((1-xi_v)*u(dofs(1)) + xi_v*u(dofs(4)));
    Yd_gl = Y_gl + scale_disp*((1-xi_v)*u(dofs(2)) + xi_v*u(dofs(5)));

    X_el{e}=X_gl; Y_el{e}=Y_gl;
    Xd_el{e}=Xd_gl; Yd_el{e}=Yd_gl;
    N_el{e}=N_dist; V_el{e}=V_dist; M_el{e}=M_dist;
end

% Scale factors for diagrams (perpendicular offset)
M_max = max(abs([M_el{:}])); if M_max<1, M_max=1; end
V_max = max(abs([V_el{:}])); if V_max<1, V_max=1; end
N_max = max(abs([N_el{:}])); if N_max<1, N_max=1; end
diag_scale = 800;  % mm offset per unit (normalized)

fig = figure('Units','centimeters','Position',[2 2 20 14]);

% Subplot 1: deflected shape
subplot(1,2,1); hold on; axis equal; grid on;
for e = 1:n_elem
    plot(X_el{e}, Y_el{e}, 'k--', 'LineWidth', 0.8);
    plot(Xd_el{e}, Yd_el{e}, 'b-', 'LineWidth', 2);
end
xlabel('x (mm)'); ylabel('y (mm)');
title(sprintf('Deflected shape (x%d)', scale_disp));
legend({'Undeformed','Deformed'}, 'Location','best');

% Subplot 2: moment, shear, axial diagrams
subplot(1,2,2); hold on; axis equal; grid on;
colors = {'b','r','g'};
labels = {'Moment M','Shear V','Axial N'};
data   = {M_el, V_el, N_el};
scales = [M_max, V_max, N_max];

for e = 1:n_elem
    c_ = cosd(alpha_deg(e)); s_ = sind(alpha_deg(e));
    nx = -s_; ny = c_;  % perpendicular unit vector
    plot(X_el{e}, Y_el{e}, 'k-', 'LineWidth', 1.5);
    for d = 1:3
        vals = data{d}{e} / scales(d) * diag_scale;
        Xoff = X_el{e} + nx*vals;
        Yoff = Y_el{e} + ny*vals;
        fill([X_el{e} fliplr(Xoff)], [Y_el{e} fliplr(Yoff)], ...
             colors{d}, 'FaceAlpha', 0.25, 'EdgeColor', colors{d});
    end
end
xlabel('x (mm)'); ylabel('y (mm)');
title('M (blue), V (red), N (green)');

exportgraphics(fig, '../figures/frame_portal_results.pdf', 'ContentType', 'vector');
exportgraphics(fig, '../figures/frame_portal_results.png', 'Resolution', 300);
fprintf('\nFigures saved to ../figures/\n');

%% Local functions
function T6 = frame_T6(c, s)
    T6 = [c  s  0  0  0  0;
         -s  c  0  0  0  0;
          0  0  1  0  0  0;
          0  0  0  c  s  0;
          0  0  0 -s  c  0;
          0  0  0  0  0  1];
end

function kbar = frame_local_k(E, A, I, L)
    a  = E*A/L;
    b  = 12*E*I/L^3;
    d  = 6*E*I/L^2;
    e4 = 4*E*I/L;
    e2 = 2*E*I/L;
    kbar = [ a   0   0  -a   0   0;
             0   b   d   0  -b   d;
             0   d  e4   0  -d  e2;
            -a   0   0   a   0   0;
             0  -b  -d   0   b  -d;
             0   d  e2   0  -d  e4];
end

function d = elem_dofs(ni, nj)
    d = [3*ni-2, 3*ni-1, 3*ni, 3*nj-2, 3*nj-1, 3*nj];
end
```

- [ ] **Step 2: Run in MATLAB**

Open MATLAB, navigate to the `matlab/` directory, and run:
```
fem_frame_portal
```

- [ ] **Step 3: Record all output values**

Copy the complete terminal output. You will substitute these exact numbers into §6.5 in later tasks. Keep these values handy — you need them for Tasks 5 and 6.

- [ ] **Step 4: Commit**
```bash
git add matlab/fem_frame_portal.m
git commit -m "Ch.6: add portal frame MATLAB script"
```

---

## Task 2: Write §6.1 — Frame Element Introduction and Figure

**Files:**
- Modify: `chapters/ch06_frame.tex` (full rewrite — replace entire file)

- [ ] **Step 1: Replace the entire content of `chapters/ch06_frame.tex`**

```latex
\chapter{2D Frame Element}
\label{ch:frame}

% =============================================================================
\section{Frame Element}
\label{sec:frame_element}
% =============================================================================

A \emph{frame element} (beam-column) combines the axial stiffness of a bar
element with the bending stiffness of an Euler-Bernoulli beam element into a
single structural member. Unlike truss elements, which transmit axial force
only through pinned joints, frame elements transmit \emph{axial force, shear
force, and bending moment} at rigid joints. Frames are the appropriate model
whenever joint rigidity induces bending in the members: portal frames,
multi-storey building frames, and any structure in which loads are not applied
along member axes.

In the \emph{local coordinate system} $(\bar{x},\bar{y})$ aligned with the
element axis, each node carries three DOFs --- axial displacement $\bar{u}$,
transverse displacement $\bar{v}$, and rotation $\theta$ (positive
counter-clockwise). The six-component local DOF vector and its conjugate
nodal force vector are:
\begin{equation}
    \bar{\vect{u}}^e =
    [\bar{u}_1,\,\bar{v}_1,\,\theta_1,\,
     \bar{u}_2,\,\bar{v}_2,\,\theta_2]\transpose, \qquad
    \bar{\vect{f}}^e =
    [\bar{N}_1,\,\bar{V}_1,\,\bar{M}_1,\,
     \bar{N}_2,\,\bar{V}_2,\,\bar{M}_2]\transpose,
    \label{eq:frame_dof}
\end{equation}
where $\bar{N}_i$ is the axial force, $\bar{V}_i$ the shear force, and
$\bar{M}_i$ the bending moment at node~$i$ (see \cref{fig:frame_element}).

\begin{figure}[htbp]
\centering
\begin{tikzpicture}[>=Stealth, x=1cm, y=1cm,
    disp/.style={->, very thick, blue!70!black},
    force/.style={->, ultra thick, red!75!black}]

  %% Element body
  \fill[gray!22] (0,-0.22) rectangle (7,0.22);
  \draw[thick]   (0,-0.22) rectangle (7,0.22);
  \node[above=4pt] at (3.5,0.22) {$E,\;A,\;I$};

  %% Nodes
  \filldraw (0,0) circle (3.5pt);
  \filldraw (7,0) circle (3.5pt);
  \node[below left=2pt]  at (0,0) {\small\bfseries 1};
  \node[below right=2pt] at (7,0) {\small\bfseries 2};

  %% Local axis label
  \draw[->,thin] (-0.3,-0.7) -- ++(1.8,0) node[right] {$\bar{x}$};
  \draw[->,thin] (-0.3,-0.7) -- ++(0,1.2) node[above] {$\bar{y}$};

  %% Displacements at node 1 (blue, above)
  \draw[disp] (-0.8,0.5) -- ++(0.8,0) node[above,midway] {$\bar{u}_1$};
  \draw[disp] (0,0.65)   -- ++(0,0.9) node[right=2pt,midway] {$\bar{v}_1$};
  \draw[disp,->] (0.05,1.72) arc[start angle=200,end angle=70,radius=0.55]
        node[above left=-2pt] {$\theta_1$};

  %% Displacements at node 2 (blue, above)
  \draw[disp] (7.0,0.5) -- ++(0.8,0) node[above,midway] {$\bar{u}_2$};
  \draw[disp] (7.0,0.65) -- ++(0,0.9) node[right=2pt,midway] {$\bar{v}_2$};
  \draw[disp,->] (7.05,1.72) arc[start angle=200,end angle=70,radius=0.55]
        node[above left=-2pt] {$\theta_2$};

  %% Forces at node 1 (red, below)
  \draw[force] (-0.9,-0.55) -- ++(0.9,0) node[below,midway] {$\bar{N}_1$};
  \draw[force] (0,-0.65)    -- ++(0,-0.9) node[right=2pt,midway] {$\bar{V}_1$};
  \draw[force,->] (0.05,-1.57) arc[start angle=160,end angle=290,radius=0.55]
        node[below left=-2pt] {$\bar{M}_1$};

  %% Forces at node 2 (red, below)
  \draw[force] (7.9,-0.55) -- ++(-0.9,0) node[below,midway] {$\bar{N}_2$};
  \draw[force] (7,-0.65)   -- ++(0,-0.9) node[right=2pt,midway] {$\bar{V}_2$};
  \draw[force,->] (7.05,-1.57) arc[start angle=160,end angle=290,radius=0.55]
        node[below left=-2pt] {$\bar{M}_2$};

  %% Length dimension
  \draw[<->,thin] (0,-2.40) -- (7,-2.40)
        node[fill=white,inner sep=1pt,midway] {$L$};

\end{tikzpicture}
\caption{Two-node frame element in the local frame $(\bar{x},\bar{y})$.
Displacements $\bar{u}_i, \bar{v}_i$ (blue) and rotations $\theta_i$
(blue, CCW positive) are the kinematic DOFs. Conjugate nodal forces
$\bar{N}_i$ (axial), $\bar{V}_i$ (shear), and moments $\bar{M}_i$ (red)
are the static conjugates.}
\label{fig:frame_element}
\end{figure}
```

- [ ] **Step 2: Quick compile check**
```bash
cd "C:\Users\Korisnik\Documents\FEM_Book Cowork"
latexmk -pdf main.tex
```
Expected: compiles; chapter heading and figure appear in the PDF. No `??` missing-reference warnings for this chapter.

- [ ] **Step 3: Commit**
```bash
git add chapters/ch06_frame.tex
git commit -m "Ch.6: §6.1 frame element intro and DOF figure"
```

---

## Task 3: Write §6.2 and §6.3 — Local Stiffness and Load Vector

**Files:**
- Modify: `chapters/ch06_frame.tex` (append to end of file)

- [ ] **Step 1: Append §6.2 to `ch06_frame.tex`**

```latex

% =============================================================================
\section{Local Stiffness Matrix}
\label{sec:frame_local_k}
% =============================================================================

In the local frame, axial deformation (bar action) and transverse bending
(beam action) are decoupled: an axial displacement $\bar{u}$ produces no
bending forces, and a transverse displacement $\bar{v}$ produces no axial
force. The local stiffness matrix $\bar{\matr{k}}^e$ is therefore
block-diagonal, assembling the bar stiffness (\cref{eq:bar_stiffness})
and the Euler-Bernoulli beam stiffness (\cref{eq:beam_stiffness}) into a
single $6\times6$ matrix:

\begin{equation}
    \bar{\matr{k}}^e =
    \renewcommand{\arraystretch}{1.6}
    \begin{bmatrix}
        \dfrac{EA}{L}   & 0 & 0 & -\dfrac{EA}{L} & 0 & 0 \\
        0 & \dfrac{12EI}{L^3} & \dfrac{6EI}{L^2} & 0 & -\dfrac{12EI}{L^3} & \dfrac{6EI}{L^2} \\
        0 & \dfrac{6EI}{L^2} & \dfrac{4EI}{L} & 0 & -\dfrac{6EI}{L^2} & \dfrac{2EI}{L} \\
        -\dfrac{EA}{L} & 0 & 0 & \dfrac{EA}{L} & 0 & 0 \\
        0 & -\dfrac{12EI}{L^3} & -\dfrac{6EI}{L^2} & 0 & \dfrac{12EI}{L^3} & -\dfrac{6EI}{L^2} \\
        0 & \dfrac{6EI}{L^2} & \dfrac{2EI}{L} & 0 & -\dfrac{6EI}{L^2} & \dfrac{4EI}{L}
    \end{bmatrix}
    \label{eq:frame_local_k}
\end{equation}

\begin{importantbox}
The frame local stiffness is block-diagonal: rows and columns $\{1,4\}$
carry the bar block ($EA/L$); rows and columns $\{2,3,5,6\}$ carry the
beam block ($EI/L^3$ scaling). Axial and bending deformations are
\emph{uncoupled} in the local frame.
\end{importantbox}

\begin{remark}
The unconstrained frame element has a three-dimensional null space:
rigid translation in $\bar{x}$, rigid translation in $\bar{y}$, and rigid
rotation. At least three linearly independent kinematic constraints are
needed to make the assembled stiffness non-singular.
\end{remark}
```

- [ ] **Step 2: Append §6.3 to `ch06_frame.tex`**

```latex

% =============================================================================
\section{Element Load Vector}
\label{sec:frame_load}
% =============================================================================

When distributed loads act along the element, the consistent nodal load
vector $\bar{\vect{f}}^e$ is formed by superposing the bar and beam
contributions derived in \cref{ch:1d,ch:beam}. For a uniform \emph{axial}
distributed load $p_0$ (positive in $+\bar{x}$) and a uniform
\emph{transverse} distributed load $q_0$ (positive in $+\bar{y}$):
\begin{equation}
    \bar{\vect{f}}^e
    =
    \underbrace{
      \frac{p_0 L}{2}
      \begin{bmatrix}1\\0\\0\\1\\0\\0\end{bmatrix}
    }_{\text{axial}\;(\cref{eq:bar_load})}
    +
    \underbrace{
      \frac{q_0 L}{12}
      \begin{bmatrix}0\\6\\L\\0\\6\\{-L}\end{bmatrix}
    }_{\text{transverse}\;(\cref{eq:beam_load})}
    \label{eq:frame_load}
\end{equation}
The two contributions act on separate DOF subsets: the axial terms go to
DOFs~1 and~4 ($\bar{u}$ displacements); the transverse terms go to DOFs~2,
3, 5, and~6 ($\bar{v}$ displacements and $\theta$ rotations).
In the worked example of \cref{sec:frame_example}, the beam (element~2)
carries a downward transverse load, so $q_0 = -q < 0$ and $p_0 = 0$;
the two columns carry no distributed load.
```

- [ ] **Step 3: Compile and check**
```bash
latexmk -pdf main.tex
```
Expected: §6.2 (with the 6×6 matrix) and §6.3 appear; no missing references.

- [ ] **Step 4: Commit**
```bash
git add chapters/ch06_frame.tex
git commit -m "Ch.6: §6.2 local stiffness matrix and §6.3 load vector"
```

---

## Task 4: Write §6.4 — Coordinate Transformation

**Files:**
- Modify: `chapters/ch06_frame.tex` (append)

- [ ] **Step 1: Append §6.4 to `ch06_frame.tex`**

```latex

% =============================================================================
\section{Coordinate Transformation}
\label{sec:frame_transform}
% =============================================================================

Frame members are generally inclined in the global $(x,y)$ plane. The
$6\times6$ transformation matrix $\matr{T}_6$ maps from the global frame to
the local frame. With $c = \cos\alpha$ and $s = \sin\alpha$, where $\alpha$
is the element inclination angle:

\begin{equation}
    \matr{T}_6 =
    \begin{bmatrix}
         c &  s & 0 & 0 & 0 & 0 \\
        -s &  c & 0 & 0 & 0 & 0 \\
         0 &  0 & 1 & 0 & 0 & 0 \\
         0 &  0 & 0 & c & s & 0 \\
         0 &  0 & 0 &-s & c & 0 \\
         0 &  0 & 0 & 0 & 0 & 1
    \end{bmatrix}
    \label{eq:frame_T6}
\end{equation}

$\matr{T}_6$ consists of two identical $3\times3$ rotation blocks (one per
node), each rotating the displacement pair $(u_i,v_i)$ from global to local
while leaving the rotation $\theta_i$ unchanged.  The structure mirrors the
truss transformation (\cref{eq:truss_T}), extended by the scalar rotation DOF.

\begin{figure}[htbp]
\centering
\begin{tikzpicture}[>=Stealth, x=1.1cm, y=1.0cm,
    disp/.style={->,thick,blue!70!black}]

  \def\alph{35}

  %% Element
  \draw[line width=2.5pt,gray!55]
        (0,0) -- ({4.5*cos(\alph)},{4.5*sin(\alph)});
  \filldraw (0,0)                        circle (3.5pt)
        node[below left=1pt] {\small\bfseries 1};
  \filldraw ({4.5*cos(\alph)},{4.5*sin(\alph)}) circle (3.5pt)
        node[above right=1pt] {\small\bfseries 2};

  %% Global axes
  \draw[->,thick] (0,0) -- ++(2.6,0) node[right] {$x$};
  \draw[->,thick] (0,0) -- ++(0,2.6) node[above] {$y$};

  %% Local axes
  \draw[->,thick,blue!75!black] (0,0)
        -- ++({2.2*cos(\alph)},{2.2*sin(\alph)}) node[right] {$\bar{x}$};
  \draw[->,thick,blue!75!black] (0,0)
        -- ++({-2.2*sin(\alph)},{2.2*cos(\alph)}) node[above left] {$\bar{y}$};

  %% Angle arc and label
  \draw (1.2,0) arc[start angle=0, end angle=\alph, radius=1.2];
  \node at ({1.5*cos(\alph/2)},{1.5*sin(\alph/2)}) {$\alpha$};

  %% Global DOFs at node 2
  \draw[->,dashed,thick,red!65!black]
        ({4.5*cos(\alph)},{4.5*sin(\alph)}) -- ++(1.1,0) node[right] {\small $u_2$};
  \draw[->,dashed,thick,red!65!black]
        ({4.5*cos(\alph)},{4.5*sin(\alph)}) -- ++(0,1.1) node[above] {\small $v_2$};

\end{tikzpicture}
\caption{Frame element inclined at angle $\alpha$ (CCW from the positive
$x$-axis). Global axes $(x,y)$; local axes $(\bar{x},\bar{y})$ in blue.
$\matr{T}_6$ rotates the global DOFs $(u_i,v_i,\theta_i)$ to local
$(\bar{u}_i,\bar{v}_i,\theta_i)$.}
\label{fig:frame_transform}
\end{figure}

\begin{importantbox}
\textbf{Angle convention:} $\alpha$ is measured counter-clockwise from the
positive $x$-axis to the element axis directed from node~1 to node~2.\\[4pt]
\begin{tabular}{@{}ll@{\quad}l@{}}
Horizontal beam, rightward: & $\alpha=0°$,   & $c=1,\;s=0$\\
Vertical column, upward:    & $\alpha=90°$,  & $c=0,\;s=1$\\
Vertical column, downward:  & $\alpha=-90°$, & $c=0,\;s=-1$
\end{tabular}
\end{importantbox}

The element stiffness and load vector in global coordinates are:
\begin{equation}
    \matr{k}^e = \matr{T}_6\transpose\,\bar{\matr{k}}^e\,\matr{T}_6 ,
    \qquad
    \vect{f}^e = \matr{T}_6\transpose\,\bar{\vect{f}}^e .
    \label{eq:frame_ke_global}
\end{equation}
For vertical columns ($c=0$, $|s|=1$), $\matr{T}_6$ exchanges the axial
and transverse DOF blocks: the column's axial stiffness $EA/L$ appears in
the vertical ($v$-$v$) entries of the global stiffness, and its bending
stiffness $12EI/L^3$ appears in the horizontal ($u$-$u$) entries.

\begin{remark}
Global assembly, boundary condition enforcement, and solution proceed
identically to the truss (\cref{ch:truss}) and beam (\cref{ch:beam})
chapters. The PVW and PMPE derivations for the frame element are direct
extensions of \cref{sec:beam_pvw,sec:beam_pmpe}: replace the $4\times4$
beam integrals with the $6\times6$ frame counterparts after applying
$\matr{T}_6$.
\end{remark}
```

- [ ] **Step 2: Compile and check**
```bash
latexmk -pdf main.tex
```
Expected: §6.4 and the inclined-element figure appear correctly.

- [ ] **Step 3: Commit**
```bash
git add chapters/ch06_frame.tex
git commit -m "Ch.6: §6.4 coordinate transformation with figure and angle-convention box"
```

---

## Task 5: Write §6.5 Part A — Example Setup Through Element Stiffnesses

**Files:**
- Modify: `chapters/ch06_frame.tex` (append)

**Prerequisite:** Task 1 complete — MATLAB script has been run and output recorded.

- [ ] **Step 1: Append the example opening through Step 3 to `ch06_frame.tex`**

```latex

% =============================================================================
\section{Worked Example: Portal Frame under Combined Loading}
\label{sec:frame_example}
% =============================================================================

\begin{example}[Portal frame --- lateral force and distributed beam load]
\label{ex:frame_portal}

A symmetric portal frame has column height $H = 3000\ \text{mm}$ and beam
span $B = 4000\ \text{mm}$. All members share the same cross-section:
$E = 210{,}000\ \text{N/mm}^2$, $A = 6{,}000\ \text{mm}^2$,
$I = 10{,}000{,}000\ \text{mm}^4$. Both column bases are clamped (fixed).
The frame carries:
\begin{itemize}
    \item a rightward horizontal force $F = 10\ \text{kN}$ at the top of
          the left column (node~2), and
    \item a uniform downward load $q = 5\ \text{N/mm}$ on the beam
          (element~2).
\end{itemize}

\begin{figure}[htbp]
\centering
\begin{tikzpicture}[>=Stealth, x=1.1cm, y=1.0cm,
    force/.style={->,ultra thick,red!75!black},
    disp/.style={->,thick,blue!70!black}]

  %% Frame members (columns: 0--3 in y; beam: 0--4 in x at y=3)
  \draw[line width=2pt,gray!60] (0,0) -- (0,3);   % left column
  \draw[line width=2pt,gray!60] (0,3) -- (4,3);   % beam
  \draw[line width=2pt,gray!60] (4,3) -- (4,0);   % right column

  %% Fixed supports (hatched rectangles + horizontal line)
  \fill[pattern=north east lines,pattern color=gray!70]
        (-0.45,-0.30) rectangle (0.45,0);
  \draw[thick] (-0.45,0) -- (0.45,0);
  \fill[pattern=north east lines,pattern color=gray!70]
        (3.55,-0.30) rectangle (4.45,0);
  \draw[thick] (3.55,0) -- (4.45,0);

  %% Distributed downward load on beam
  \foreach \xp in {0.2,0.5,...,3.8}{
      \draw[disp,->,thin] (\xp,3.65) -- (\xp,3.12);
  }
  \draw[disp,thick] (0,3.65) -- (4,3.65);
  \node[above=2pt,blue!70!black] at (2,3.68) {$q = 5\ \text{N/mm}$};

  %% Lateral force at node 2
  \draw[force] (-1.6,3) -- (-0.08,3)
       node[above=3pt,midway] {$F = 10\ \text{kN}$};

  %% Nodes
  \filldraw (0,0) circle (3pt) node[below left=3pt] {\small\bfseries 1};
  \filldraw (0,3) circle (3pt) node[above left=3pt] {\small\bfseries 2};
  \filldraw (4,3) circle (3pt) node[above right=3pt] {\small\bfseries 3};
  \filldraw (4,0) circle (3pt) node[below right=3pt] {\small\bfseries 4};

  %% Element numbers
  \node[left=10pt]  at (0,1.5) {\small\textcircled{\small 1}};
  \node[above=6pt]  at (2,3)   {\small\textcircled{\small 2}};
  \node[right=10pt] at (4,1.5) {\small\textcircled{\small 3}};

  %% Dimension lines
  \draw[<->,thin] (-0.9,0) -- (-0.9,3)
       node[fill=white,inner sep=1pt,midway,rotate=90]
       {$H = 3000\ \text{mm}$};
  \draw[<->,thin] (0,-0.7) -- (4,-0.7)
       node[fill=white,inner sep=1pt,midway] {$B = 4000\ \text{mm}$};

\end{tikzpicture}
\caption{Portal frame: four nodes, three elements, clamped bases at
nodes~1 and~4. Combined loading: horizontal force $F$ at node~2 and
uniform vertical load $q$ on the beam.}
\label{fig:frame_portal}
\end{figure}

\medskip
\noindent\textbf{DOF ordering:}
$[u_1,v_1,\theta_1,\,u_2,v_2,\theta_2,\,u_3,v_3,\theta_3,\,u_4,v_4,\theta_4]
= \text{DOFs}\ [1\ldots12]$.\\
Clamped DOFs: $\{1,2,3,10,11,12\}=\mathbf{0}$.
Free DOFs: $\{4,5,6,7,8,9\}$.

\medskip
\noindent\textbf{Step 1: Element data.}

\begin{table}[htbp]
\centering
\small
\caption{Element data for the portal frame example.}
\label{tab:frame_elements}
\begin{tabular}{@{}ccccccc@{}}
\toprule
Elem & Nodes & $L$ (mm) & $\alpha$ ($^\circ$) & $c$ & $s$ & Global DOFs \\
\midrule
1 & $1\to2$ & 3000 &  90 & $0$ & $+1$ & 1,2,3,4,5,6 \\
2 & $2\to3$ & 4000 &   0 & $1$ & $0$  & 4,5,6,7,8,9 \\
3 & $3\to4$ & 3000 & $-90$ & $0$ & $-1$ & 7,8,9,10,11,12 \\
\bottomrule
\end{tabular}
\end{table}

\medskip
\noindent\textbf{Step 2: Local stiffness constants.}

Let $a=EA/L$, $b=12EI/L^3$, $d=6EI/L^2$, $e_4=4EI/L$, $e_2=2EI/L$.

\smallskip
For the columns (elements 1 and 3, $L_c = 3000\ \text{mm}$):
\begin{align*}
  a_c &= \frac{210{,}000 \times 6{,}000}{3{,}000} = 420{,}000\ \tfrac{\text{N}}{\text{mm}},
  &b_c &= \frac{12 \times 210{,}000 \times 10^7}{3{,}000^3} = 933.3\ \tfrac{\text{N}}{\text{mm}},\\
  d_c &= \frac{6 \times 210{,}000 \times 10^7}{3{,}000^2} = 1{,}400{,}000\ \text{N},
  &e_{4c} &= \frac{4 \times 210{,}000 \times 10^7}{3{,}000} = 2.800\times10^9\ \text{N\,mm},\\
  &&e_{2c} &= \tfrac{1}{2}\,e_{4c} = 1.400\times10^9\ \text{N\,mm}.
\end{align*}

For the beam (element 2, $L_b = 4000\ \text{mm}$):
\begin{align*}
  a_b &= \frac{210{,}000 \times 6{,}000}{4{,}000} = 315{,}000\ \tfrac{\text{N}}{\text{mm}},
  &b_b &= \frac{12 \times 210{,}000 \times 10^7}{4{,}000^3} = 393.75\ \tfrac{\text{N}}{\text{mm}},\\
  d_b &= \frac{6 \times 210{,}000 \times 10^7}{4{,}000^2} = 787{,}500\ \text{N},
  &e_{4b} &= \frac{4 \times 210{,}000 \times 10^7}{4{,}000} = 2.100\times10^9\ \text{N\,mm},\\
  &&e_{2b} &= 1.050\times10^9\ \text{N\,mm}.
\end{align*}

\medskip
\noindent\textbf{Step 3: Transformation matrices.}

For elements~1 and~3, the column geometry gives $c=0$, $|s|=1$:
\[
    \matr{T}_6^{(1)} =
    \begin{bmatrix}
         0 &  1 & 0 & 0 &  0 & 0 \\
        -1 &  0 & 0 & 0 &  0 & 0 \\
         0 &  0 & 1 & 0 &  0 & 0 \\
         0 &  0 & 0 & 0 &  1 & 0 \\
         0 &  0 & 0 &-1 &  0 & 0 \\
         0 &  0 & 0 & 0 &  0 & 1
    \end{bmatrix}\!,\quad
    \matr{T}_6^{(3)} =
    \begin{bmatrix}
         0 & -1 & 0 &  0 &  0 & 0 \\
         1 &  0 & 0 &  0 &  0 & 0 \\
         0 &  0 & 1 &  0 &  0 & 0 \\
         0 &  0 & 0 &  0 & -1 & 0 \\
         0 &  0 & 0 &  1 &  0 & 0 \\
         0 &  0 & 0 &  0 &  0 & 1
    \end{bmatrix}\!.
\]
For the beam ($\alpha=0°$): $\matr{T}_6^{(2)} = \matr{I}_6$, so the beam's
local and global stiffness matrices are identical.
```

- [ ] **Step 2: Compile and check**
```bash
latexmk -pdf main.tex
```
Expected: example opening, TikZ figure, element table, stiffness constants, and transformation matrices appear.

- [ ] **Step 3: Commit**
```bash
git add chapters/ch06_frame.tex
git commit -m "Ch.6: §6.5 example -- setup, element data, stiffness constants, T6 matrices"
```

---

## Task 6: Write §6.5 Part B — Global Stiffness, BCs, Solution, Reactions, Internal Forces

**Files:**
- Modify: `chapters/ch06_frame.tex` (append)

**Prerequisite:** Task 1 output values recorded. Substitute every `[VAL_*]` below with the exact printed values.

- [ ] **Step 1: Append Steps 4–6 (global stiffness, load vector, BC reduction)**

```latex

\medskip
\noindent\textbf{Step 4: Global element stiffnesses.}

Applying $\matr{k}^{(e)} = (\matr{T}_6^{(e)})\transpose\bar{\matr{k}}^{(e)}\matr{T}_6^{(e)}$
to the columns with $c=0$, $s=\pm1$ exchanges the axial and transverse
blocks in the global frame: the column's axial stiffness $a_c$ appears in
the $v$-$v$ coupling (vertical DOFs), and the bending stiffness $b_c$
appears in the $u$-$u$ coupling (horizontal DOFs). For element~1
(DOFs 1,2,3,4,5,6):
\begin{equation*}
    \matr{k}^{(1)} =
    \begin{bmatrix}
         b_c  &  0   &  d_c  & -b_c &  0   &  d_c  \\
         0    &  a_c &  0    &  0   & -a_c &  0    \\
         d_c  &  0   & e_{4c}& -d_c &  0   & e_{2c}\\
        -b_c  &  0   & -d_c  &  b_c &  0   & -d_c  \\
         0    & -a_c &  0    &  0   &  a_c &  0    \\
         d_c  &  0   & e_{2c}& -d_c &  0   & e_{4c}
    \end{bmatrix}
\end{equation*}
and symmetrically for element~3 (DOFs 7,8,9,10,11,12) with $s=-1$
(sign changes on the $d_c$ off-diagonal terms involving the $v$--$\theta$
coupling). Element~2 is unchanged by the identity transformation.

\medskip
\noindent\textbf{Step 5: Global assembly.}

The $12\times12$ global stiffness $\stiff$ is assembled by adding each
$\matr{k}^{(e)}$ into the positions given by its DOF list
(Table~\ref{tab:frame_elements}). Node~2 (DOFs 4,5,6) receives
contributions from elements~1 and~2; node~3 (DOFs 7,8,9) receives
contributions from elements~2 and~3. The full $12\times12$ matrix is
computed by \Cref{lst:frame_portal}.

\medskip
\noindent\textbf{Step 6: Load vector and boundary conditions.}

The global load vector has two contributions:
\begin{itemize}
    \item \textbf{Lateral force} $F=10{,}000\ \text{N}$ (rightward) at
          DOF~4 ($u_2$): $F_4 = +10{,}000\ \text{N}$.
    \item \textbf{Distributed load} $q = 5\ \text{N/mm}$ (downward,
          $q_0 = -5\ \text{N/mm}$) on element~2 ($L_b = 4000\ \text{mm}$,
          $\matr{T}_6^{(2)} = \matr{I}$):
          \[
            \vect{f}^{(2)}_q
            = \frac{q_0 L_b}{12}
              \begin{bmatrix}0\\6\\L_b\\0\\6\\{-L_b}\end{bmatrix}
            = \frac{-5 \times 4000}{12}
              \begin{bmatrix}0\\6\\4000\\0\\6\\{-4000}\end{bmatrix}
            =
              \begin{bmatrix}
                0\\ -10{,}000\\ -6.667\times10^6\\
                0\\ -10{,}000\\ +6.667\times10^6
              \end{bmatrix}
              \begin{bmatrix}\\\text{N}\\\text{N\,mm}\\\\\text{N}\\\text{N\,mm}\end{bmatrix}
              \quad\text{(DOFs 4--9)}.
          \]
\end{itemize}
Applying the clamped BCs (DOFs $\{1,2,3,10,11,12\} = 0$) and retaining only
the six free DOFs gives the reduced $6\times6$ system:
\[
    \stiff_{ff}\,\vect{u}_f = \vect{F}_f,
    \qquad
    \vect{u}_f = [u_2,v_2,\theta_2,u_3,v_3,\theta_3]\transpose,
\]
\[
    \vect{F}_f =
    [+10{,}000,\; -10{,}000,\; -6.667\times10^6,\;
     0,\; -10{,}000,\; +6.667\times10^6]\transpose\ \text{(N and N\,mm)}.
\]
```

- [ ] **Step 2: Append Steps 7–9 (solution, reactions, internal forces) with MATLAB-verified numbers**

Replace every `[VAL_*]` with the number printed by `fem_frame_portal.m`.

```latex

\medskip
\noindent\textbf{Step 7: Solution (MATLAB-verified, \Cref{lst:frame_portal}).}

\begin{align*}
    u_2 &= [VAL_u2]\ \text{mm}, &
    v_2 &= [VAL_v2]\ \text{mm}, &
    \theta_2 &= [VAL_th2]\ \text{rad}, \\
    u_3 &= [VAL_u3]\ \text{mm}, &
    v_3 &= [VAL_v3]\ \text{mm}, &
    \theta_3 &= [VAL_th3]\ \text{rad}.
\end{align*}
The lateral sway ($u_2 \approx u_3$) reflects the rigid-body sway induced
by $F$; the small difference between $u_2$ and $u_3$ arises from the beam's
axial flexibility. Both joints deflect downward ($v_2, v_3 < 0$) under the
distributed gravity load.

\medskip
\noindent\textbf{Step 8: Support reactions.}

\begin{align*}
    R_{1x} &= [VAL_Rx1]\ \text{N}, &
    R_{1y} &= [VAL_Ry1]\ \text{N}, &
    M_1    &= [VAL_M1]\ \text{N\,mm}, \\
    R_{4x} &= [VAL_Rx4]\ \text{N}, &
    R_{4y} &= [VAL_Ry4]\ \text{N}, &
    M_4    &= [VAL_M4]\ \text{N\,mm}.
\end{align*}
\noindent Equilibrium check:
\begin{align*}
    \textstyle\sum F_x &= R_{1x} + R_{4x} + F = 0\ \checkmark\\
    \textstyle\sum F_y &= R_{1y} + R_{4y} - qB = 0\ \checkmark
\end{align*}

\medskip
\noindent\textbf{Step 9: Element internal forces.}

Local end forces are recovered via
$\bar{\vect{f}}^{(e)}_{\text{int}} =
\bar{\matr{k}}^{(e)}\bigl(\matr{T}_6^{(e)}\vect{u}^{(e)}\bigr)
- \bar{\vect{f}}^{(e)}_{\text{load}}$
(see \cref{sec:frame_internal}). The values at each end are:

\begin{table}[htbp]
\centering
\small
\caption{Element end forces in local coordinates (N and N\,mm).
Tension positive for $\bar{N}$; sagging positive for $\bar{M}$.}
\label{tab:frame_internal}
\begin{tabular}{@{}lrrrrrr@{}}
\toprule
Element & $\bar{N}_i$ & $\bar{V}_i$ & $\bar{M}_i$
        & $\bar{N}_j$ & $\bar{V}_j$ & $\bar{M}_j$ \\
        & (N) & (N) & (N\,mm) & (N) & (N) & (N\,mm)\\
\midrule
1 -- left col  & [N1i] & [V1i] & [M1i] & [N1j] & [V1j] & [M1j] \\
2 -- beam      & [N2i] & [V2i] & [M2i] & [N2j] & [V2j] & [M2j] \\
3 -- right col & [N3i] & [V3i] & [M3i] & [N3j] & [V3j] & [M3j] \\
\bottomrule
\end{tabular}
\end{table}

The MATLAB script (\Cref{lst:frame_portal}) plots the deflected shape and
the bending moment, shear force, and axial force diagrams
(\cref{fig:frame_portal_results}).

\begin{figure}[htbp]
\centering
\includegraphics[width=0.88\textwidth]{figures/frame_portal_results}
\caption{Portal frame results: deflected shape (left, displacements
magnified $\times100$) and internal force diagrams (right) for bending
moment (blue), shear (red), and axial force (green).}
\label{fig:frame_portal_results}
\end{figure}

\end{example}
```

- [ ] **Step 3: Compile and check — confirm no `[VAL_*]` remain in the PDF**
```bash
latexmk -pdf main.tex
```
Open the PDF at Chapter 6, verify all numeric values are filled in and the figure loads.

- [ ] **Step 4: Commit**
```bash
git add chapters/ch06_frame.tex
git commit -m "Ch.6: §6.5 complete with MATLAB-verified numerics, reactions, and internal forces"
```

---

## Task 7: Write §6.6, Add Appendix Listing, and Build Final PDF

**Files:**
- Modify: `chapters/ch06_frame.tex` (append §6.6)
- Modify: `appendix/matlab_listings.tex` (append listing block)

- [ ] **Step 1: Append §6.6 to `ch06_frame.tex`**

```latex

% =============================================================================
\section{Internal Force Recovery}
\label{sec:frame_internal}
% =============================================================================

Once the global displacement vector $\vect{u}$ is known, the internal
forces at any point along a frame element are recovered in four steps:

\begin{enumerate}
    \item \textbf{Extract global nodal DOFs.}
          For element $e$ connecting nodes $i$ and $j$, assemble
          $\vect{u}^e = [u_i,\,v_i,\,\theta_i,\,u_j,\,v_j,\,\theta_j]\transpose$.

    \item \textbf{Transform to local coordinates.}
          \[
              \bar{\vect{u}}^e = \matr{T}_6\,\vect{u}^e .
          \]

    \item \textbf{Compute local end forces.}
          \begin{equation}
              \bar{\vect{f}}^e_{\mathrm{int}}
              = \bar{\matr{k}}^e\,\bar{\vect{u}}^e
                - \bar{\vect{f}}^e_{\mathrm{load}} ,
              \label{eq:frame_int_forces}
          \end{equation}
          where $\bar{\vect{f}}^e_{\mathrm{load}}$ is the element load
          vector from \cref{eq:frame_load} (zero for elements carrying no
          distributed load). The six components are the end axial forces,
          shear forces, and bending moments in local coordinates.

    \item \textbf{Distribute along the element.}
          Between the ends, the internal forces vary with position
          $\bar{x} \in [0,L]$:
          \begin{align}
              N(\bar{x}) &= \bar{N}_i
                            + \frac{EA}{L}(\bar{u}_j - \bar{u}_i)
                              \frac{\bar{x}}{L}
                            \quad\text{(linear under axial load; constant if none)},
                            \label{eq:frame_N_dist}\\
              V(\bar{x}) &= \bar{V}_i + q_0\,\bar{x}
                            \quad\text{(linear under uniform }q_0;
                            \text{ constant if }q_0=0\text{)},
                            \label{eq:frame_V_dist}\\
              M(\bar{x}) &= \bar{M}_i - \bar{V}_i\,\bar{x}
                            - \tfrac{1}{2}q_0\,\bar{x}^2
                            \quad\text{(quadratic under }q_0;
                            \text{ linear if }q_0=0\text{)}.
                            \label{eq:frame_M_dist}
          \end{align}
\end{enumerate}

\begin{remark}
For a vertical column ($\alpha = 90°$), the local $\bar{x}$ axis points
upward and the local $\bar{y}$ axis points leftward. The local axial force
$\bar{N}$ is therefore the vertical internal force (gravity direction) in
the column, and the local shear $\bar{V}$ is horizontal. In the portal
frame example, the column axial forces carry the vertical gravity load from
the beam, while the column shear forces are in equilibrium with the applied
lateral force $F$.
\end{remark}
```

- [ ] **Step 2: Read the end of `appendix/matlab_listings.tex` to find where to append**

Use the Read tool on `appendix/matlab_listings.tex`, reading from around line 380 onwards, to identify the last listing entry.

- [ ] **Step 3: Append the frame listing to `appendix/matlab_listings.tex`**

Open `appendix/matlab_listings.tex` and append the following block after the last existing listing:

```latex

% ------------------------------------------------------------------
\section*{Chapter 6 --- 2D Frame Element}

\lstinputlisting[
    caption={Portal frame under combined loading
             (\Cref{ex:frame_portal,lst:frame_portal}).},
    label={lst:frame_portal}
]{matlab/fem_frame_portal.m}
```

**Before appending:** open `matlab/fem_frame_portal.m` and confirm there are no em-dashes (---) inside comments or strings. The `%--` separator lines already use hyphens, so this should be clean.

- [ ] **Step 4: Final PDF build**
```bash
latexmk -pdf main.tex
```

- [ ] **Step 5: Visual verification checklist** — open the PDF and confirm:
  - [ ] Chapter 6 appears in the table of contents with correct section numbers (§6.1–6.6)
  - [ ] §6.1 figure (`fig:frame_element`) renders correctly — 6 DOF arrows, local axes
  - [ ] §6.2 has the full 6×6 matrix and the `importantbox`
  - [ ] §6.3 shows the split load vector formula
  - [ ] §6.4 inclined-element figure renders; the `importantbox` angle table is present
  - [ ] §6.5 TikZ portal frame figure renders; all `[VAL_*]` replaced by numbers; `tab:frame_internal` is fully filled; `fig:frame_portal_results` loads
  - [ ] §6.6 is present with three numbered equations
  - [ ] Appendix contains `lst:frame_portal`
  - [ ] No `??` missing-reference warnings anywhere in the chapter

- [ ] **Step 6: Final commit**
```bash
git add chapters/ch06_frame.tex appendix/matlab_listings.tex
git commit -m "Ch.6: §6.6 internal force recovery, appendix listing -- chapter complete"
```

---

## Self-Review

### Spec Coverage

| Spec requirement | Task |
|---|---|
| §6.1 frame element intro + TikZ 6-DOF figure | Task 2 |
| §6.2 6×6 local stiffness + importantbox + null-space remark | Task 3 |
| §6.3 combined load vector (axial + transverse) | Task 3 |
| §6.4 transformation + TikZ inclined figure + angle-convention box | Task 4 |
| §6.5 portal frame: geometry, fixed bases, combined loading | Task 5 |
| Element data table with L, α, c, s, DOFs | Task 5 |
| Local stiffness constants (columns + beam) | Task 5 |
| T₆ matrices for all three elements | Task 5 |
| Global element stiffnesses (column transformation effect explained) | Task 6 |
| Assembly description | Task 6 |
| Load vector computation (point force + distributed) | Task 6 |
| BC reduction → 6×6 system | Task 6 |
| Solution with MATLAB-verified values | Task 6 |
| Reactions + equilibrium check | Task 6 |
| Internal forces table | Task 6 |
| Figures: deflected shape + M/V/N diagrams | Tasks 1, 6 |
| §6.6 internal force recovery procedure | Task 7 |
| MATLAB script `fem_frame_portal.m` | Task 1 |
| Appendix listing `lst:frame_portal` | Task 7 |
| All numerics MATLAB-verified | Task 1 → Tasks 5, 6 |
| No toolboxes in MATLAB | Task 1 |
| em-dashes replaced in listing | Task 7 |
| Label `eq:frame_T` replaced with `eq:frame_T6` | Task 4 |

### Placeholder Scan

Tasks 6 Steps 1-2 contain `[VAL_*]` and `[N1i]`-style markers. These are intentional engineer actions (substitute MATLAB output from Task 1), not unresolved spec gaps. They must all be replaced before the Task 6 commit.

### Type Consistency

- `eq:frame_local_k` defined in Task 3, referenced in §6.4 remark and §6.6 — consistent.
- `eq:frame_T6` defined in Task 4, consistent name throughout Tasks 5–6.
- `eq:frame_load` defined in Task 3, referenced in §6.6 Step 3 — consistent.
- `eq:frame_int_forces` defined in Task 7 §6.6, referenced by §6.5 Step 9 prose — consistent.
- `tab:frame_elements` defined in Task 5, referenced in Task 6 Step 1 prose — consistent.
- `tab:frame_internal` defined in Task 6, referenced in §6.5 prose — consistent.
- `fig:frame_portal` defined in Task 5, referenced in caption only — consistent.
- `fig:frame_portal_results` defined in Task 6, must match filename `figures/frame_portal_results` saved by MATLAB — verify match.
