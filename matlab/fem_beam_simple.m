% fem_beam_simple.m -- Chapter 5: Worked Example
% Simply supported Euler-Bernoulli beam, 2 equal elements
% L = 2 m, EI = 1e4 N*m^2, P = 1000 N at midspan
% Sign convention: upward-positive displacement, CCW-positive rotation
% Units: N, m
clear; clc; close all;

%% --- Parameters ---
L  = 2;      % total beam length [m]
Le = L/2;    % element length [m]
EI = 1e4;    % bending stiffness [N*m^2]
P  = 1000;   % applied load [N] (downward, enters as -P)

nNodes = 3;  nEl = 2;  nDOF = 2*nNodes;   % 6 DOFs total
% DOF map: node i -> (2i-1) = v_i, (2i) = theta_i
% DOFs: [v1 th1 v2 th2 v3 th3] = [1 2 3 4 5 6]

%% --- Element stiffness matrix ---
K0 = EI / Le^3;
fprintf('K0 = EI/Le^3 = %.4g N/m\n\n', K0);

Ce = [12   6*Le   -12   6*Le;
      6*Le 4*Le^2 -6*Le 2*Le^2;
     -12  -6*Le    12  -6*Le;
      6*Le 2*Le^2 -6*Le 4*Le^2];

ke = K0 * Ce;
fprintf('=== Element stiffness matrix k_e [N/m, N, N*m] ===\n');
fprintf('(K0 = %.4g N/m; Le = %.4g m)\n', K0, Le);
disp(ke);

%% --- Global assembly ---
K = zeros(nDOF);
for e = 1:nEl
    n1 = e; n2 = e+1;
    dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2];
    K(dofs, dofs) = K(dofs, dofs) + ke;
end
fprintf('=== Assembled 6x6 Global K [N/m or N or N*m] ===\n');
disp(K/K0);
fprintf('(all entries above multiplied by 1/K0)\n\n');

%% --- Load vector and BCs ---
F = zeros(nDOF, 1);
F(3) = -P;   % downward load at node 2 (DOF 3 = v2)

fixedDOFs = [1, 5];   % v1=0, v3=0 (pin and roller)
freeDOFs  = [2, 3, 4, 6];

Kff = K(freeDOFs, freeDOFs);
Ff  = F(freeDOFs);

fprintf('=== Reduced 4x4 System (free DOFs: th1 v2 th2 th3) ===\n');
fprintf('Kff =\n'); disp(Kff);
fprintf('Ff = [%.1f; %.1f; %.1f; %.1f] N or N*m\n\n', Ff(1),Ff(2),Ff(3),Ff(4));

%% --- Solve ---
u_free = Kff \ Ff;
u = zeros(nDOF, 1);
u(freeDOFs) = u_free;

fprintf('=== Nodal Displacements ===\n');
fprintf('v1    = %+.6f m  (fixed)\n', u(1));
fprintf('th1   = %+.6f rad\n', u(2));
fprintf('v2    = %+.6f m  (= %.4f mm)\n', u(3), u(3)*1000);
fprintf('th2   = %+.6f rad\n', u(4));
fprintf('v3    = %+.6f m  (fixed)\n', u(5));
fprintf('th3   = %+.6f rad\n', u(6));

% Analytical check
v2_analytic = P*L^3/(48*EI);
fprintf('\nAnalytical v2 = PL^3/(48EI) = %.6f m = %.4f mm\n', v2_analytic, v2_analytic*1000);
fprintf('FEM error = %.2e %%\n\n', abs(u(3)-v2_analytic)/v2_analytic*100);

%% --- Reactions ---
R = K*u - F;
fprintf('=== Support Reactions ===\n');
fprintf('R1 (v1): %+.4f N\n', R(1));
fprintf('R3 (v3): %+.4f N\n', R(5));
fprintf('Sum F:   %+.6f N (should be 0)\n\n', R(1)+R(5)-P);

%% --- Bending moment along beam ---
n_pts = 51;
M_all = []; x_all = [];
for e = 1:nEl
    n1 = e; n2 = e+1;
    dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2];
    ue = u(dofs);
    xi = linspace(0, Le, n_pts);
    % Second derivative of Hermite functions (curvature-displacement matrix)
    H1pp = (12*xi/Le^3 - 6/Le^2);
    H2pp = (6*xi/Le^2 - 4/Le);
    H3pp = (-12*xi/Le^3 + 6/Le^2);
    H4pp = (6*xi/Le^2 - 2/Le);
    % Bending moment M = EI * v''  (sagging positive)
    My = EI * (H1pp*ue(1) + H2pp*ue(2) + H3pp*ue(3) + H4pp*ue(4));
    x_el = (e-1)*Le + xi;
    M_all = [M_all, My];
    x_all = [x_all, x_el];
end
M_max = max(abs(M_all));
fprintf('=== Bending Moment ===\n');
fprintf('M_max = %.4f N*m  (analytical PL/4 = %.4f N*m)\n\n', M_max, P*L/4);

%% --- Figure ---
figure('Color','white','Position',[100 100 800 500]);
subplot(2,1,1);
x_plot = [0, Le, 2*Le]; v_plot = [u(1), u(3), u(5)]*1000;
xi_fine = linspace(0, L, 200);
v_fine = zeros(size(xi_fine));
for i = 1:length(xi_fine)
    xi = xi_fine(i);
    e = min(floor(xi/Le)+1, nEl);
    x_loc = xi - (e-1)*Le;
    n1=e; n2=e+1; dofs=[2*n1-1,2*n1,2*n2-1,2*n2]; ue=u(dofs);
    h1 = 1 - 3*(x_loc/Le)^2 + 2*(x_loc/Le)^3;
    h2 = x_loc*(1-x_loc/Le)^2;
    h3 = 3*(x_loc/Le)^2 - 2*(x_loc/Le)^3;
    h4 = x_loc*((x_loc/Le)^2 - x_loc/Le);
    v_fine(i) = (h1*ue(1)+h2*ue(2)+h3*ue(3)+h4*ue(4))*1000;
end
plot(xi_fine, v_fine, 'b-', 'LineWidth', 2);
xlabel('x [m]'); ylabel('v [mm]');
title('Deflection'); grid on; set(gca,'YDir','reverse');

subplot(2,1,2);
plot(x_all, M_all, 'r-', 'LineWidth', 2);
xlabel('x [m]'); ylabel('M [N\cdotm]');
title('Bending moment'); grid on;
yline(0,'k--');

book_style(gcf);
exportgraphics(gcf,'../figures/beam_simple_results.pdf','ContentType','vector');
exportgraphics(gcf,'../figures/beam_simple_results.png','Resolution',300);
fprintf('Figures saved.\n');
