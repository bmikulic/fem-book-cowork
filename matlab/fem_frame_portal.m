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
    V_dist = fint(2) + q0 * xv;
    M_dist = fint(3) - fint(2)*xv - q0*xv.^2/2;

    X_gl = xi + xi_v*(xj - xi);
    Y_gl = yi + xi_v*(yj - yi);

    Xd_gl = X_gl + scale_disp*((1-xi_v)*u(dofs(1)) + xi_v*u(dofs(4)));
    Yd_gl = Y_gl + scale_disp*((1-xi_v)*u(dofs(2)) + xi_v*u(dofs(5)));

    X_el{e}=X_gl; Y_el{e}=Y_gl;
    Xd_el{e}=Xd_gl; Yd_el{e}=Yd_gl;
    N_el{e}=N_dist; V_el{e}=V_dist; M_el{e}=M_dist;
end

M_max = max(abs([M_el{:}])); if M_max<1, M_max=1; end
V_max = max(abs([V_el{:}])); if V_max<1, V_max=1; end
N_max = max(abs([N_el{:}])); if N_max<1, N_max=1; end
diag_scale = 800;

fig = figure('Units','centimeters','Position',[2 2 20 14]);

subplot(1,2,1); hold on; axis equal; grid on;
for e = 1:n_elem
    plot(X_el{e}, Y_el{e}, 'k--', 'LineWidth', 0.8);
    plot(Xd_el{e}, Yd_el{e}, 'b-', 'LineWidth', 2);
end
xlabel('x (mm)'); ylabel('y (mm)');
title(sprintf('Deflected shape (x%d)', scale_disp));
legend({'Undeformed','Deformed'}, 'Location','best');

subplot(1,2,2); hold on; axis equal; grid on;
colors = {'b','r','g'};
data   = {M_el, V_el, N_el};
scales = [M_max, V_max, N_max];

for e = 1:n_elem
    c_ = cosd(alpha_deg(e)); s_ = sind(alpha_deg(e));
    nx = -s_; ny = c_;
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

book_style(fig);
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
