% fem_beam_3elem.m -- Chapter 5: Worked Example
% Three-span continuous beam (from Analiza linijskog nosaca.pdf)
% Nodes: 1=0mm, 2=3000mm, 3=5000mm, 4=6000mm
% Cross-section: I-profile, I = 5481140 mm^4, z_max = 70 mm
% E = 210000 N/mm^2  -> EI in N*mm^2
% Loading: q = 10 N/mm (downward) on element 1 only
%          P = 30000 N (downward) at node 3
% BCs: pin at node 1 (w1=0), roller at node 2 (w2=0), roller at node 4 (w4=0)
%      rotations free at all nodes
% Sign convention: UPWARD-positive displacement and CCW-positive rotation
% (deflections downward will appear as negative values)
clear; clc; close all;

%% --- Parameters ---
E  = 210000;      % N/mm^2
Iy = 5481140;     % mm^4
EI = E * Iy;      % N*mm^2
z_max = 70;       % mm (half-height of I-section)

xn = [0; 3000; 5000; 6000];   % node coordinates [mm]
Le = diff(xn);                % element lengths: [3000, 2000, 1000] mm
nEl = 3;  nNodes = 4;  nDOF = 2*nNodes;
% DOF map: node i -> 2i-1 = w_i, 2i = theta_i
% DOFs: [w1 th1 w2 th2 w3 th3 w4 th4] = [1 2 3 4 5 6 7 8]

q1 = 10;      % N/mm, distributed load on element 1 (downward)
P  = 30000;   % N, concentrated load at node 3 (downward)

%% --- Element stiffness matrices ---
fprintf('=== Element Properties ===\n');
for e = 1:nEl
    L = Le(e);
    EIL3 = EI / L^3;
    ke{e} = EIL3 * [12,   6*L,  -12,   6*L;
                     6*L,  4*L^2, -6*L,  2*L^2;
                    -12,  -6*L,   12,  -6*L;
                     6*L,  2*L^2, -6*L,  4*L^2];
    fprintf('Element %d: L = %g mm, EI/L^3 = %.6g N/mm^3\n', e, L, EIL3);
    fprintf('  k_e(×10^-9) =\n'); disp(ke{e}*1e9);
end

%% --- Element load vectors ---
% Uniform load q on element 1 (downward -> negative in upward-positive convention)
L1 = Le(1);
fq1 = -q1 * L1/12 * [6; L1; 6; -L1];   % [N, N*mm, N, N*mm]
fprintf('Element 1 load vector (distributed q=%g N/mm):\n', q1);
fprintf('  fq1 = [%.4g; %.4g; %.4g; %.4g] [N or N*mm]\n\n', fq1(1),fq1(2),fq1(3),fq1(4));

%% --- Global assembly ---
K = zeros(nDOF);
for e = 1:nEl
    n1 = e; n2 = e+1;
    dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2];
    K(dofs, dofs) = K(dofs, dofs) + ke{e};
end

%% --- Global load vector ---
F = zeros(nDOF, 1);
% Distributed load on element 1 (DOFs 1,2,3,4)
F([1,2,3,4]) = F([1,2,3,4]) + fq1;
% Concentrated load at node 3, DOF 5 (w3), downward
F(5) = F(5) - P;

fprintf('=== Global Load Vector FG [N or N*mm] ===\n');
fprintf('  FG = ['); fprintf('%.4g ', F'); fprintf(']\n\n');

%% --- BCs: fixed DOFs = {1,3,7} (w1=0, w2=0, w4=0) ---
fixedDOFs = [1, 3, 7];
freeDOFs  = [2, 4, 5, 6, 8];

%% --- Reduced system ---
Kff = K(freeDOFs, freeDOFs);
Ff  = F(freeDOFs);

fprintf('=== Reduced 5x5 Stiffness Matrix Kff (×10^-9) ===\n');
disp(Kff * 1e9);
fprintf('Ff = ['); fprintf('%.4g ', Ff'); fprintf('] [N or N*mm]\n\n');

%% --- Solve ---
u_free = Kff \ Ff;
u = zeros(nDOF, 1);
u(freeDOFs) = u_free;

fprintf('=== Nodal DOFs ===\n');
labels = {'w1','th1','w2','th2','w3','th3','w4','th4'};
for i = 1:nDOF
    if any(fixedDOFs == i)
        fprintf('  %s = 0 (fixed)\n', labels{i});
    else
        fprintf('  %s = %+.6g', labels{i}, u(i));
        if mod(i,2)==1, fprintf(' mm\n'); else, fprintf(' rad\n'); end
    end
end

%% --- Reactions ---
R = K*u - F;
fprintf('\n=== Support Reactions ===\n');
fprintf('  R1 (w1): %+.4g N\n',  R(1));
fprintf('  R2 (w2): %+.4g N\n',  R(3));
fprintf('  R4 (w4): %+.4g N\n',  R(7));
fprintf('  Sum F_vertical: %.6g N (should be 0)\n\n', R(1)+R(3)+R(7)-P-q1*Le(1));

%% --- Internal fields ---
n_pts = 101;
w_all = []; th_all = []; My_all = []; Qz_all = []; x_all = [];

for e = 1:nEl
    n1 = e; n2 = e+1;
    dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2];
    ue  = u(dofs);
    L   = Le(e);
    xi  = linspace(0, L, n_pts);   % local coordinate
    q_e = 0;
    if e == 1, q_e = q1; end

    % Hermite shape functions and derivatives
    H1  = 1 - 3*(xi/L).^2 + 2*(xi/L).^3;
    H2  = xi.*(1 - xi/L).^2;
    H3  = 3*(xi/L).^2 - 2*(xi/L).^3;
    H4  = xi.*((xi/L).^2 - xi/L);
    H1p = (-6*xi/L^2 + 6*xi.^2/L^3);
    H2p = (1 - xi/L).^2 + xi.*2.*(1-xi/L).*(-1/L);
    H3p = (6*xi/L^2 - 6*xi.^2/L^3);
    H4p = (xi/L).^2 + xi.*2.*(xi/L).*(1/L) - xi./L.^2.*L;  % simplified below
    % Cleaner:
    H1pp = 12*xi/L^3 - 6/L^2;
    H2pp = 6*xi/L^2 - 4/L;
    H3pp = -12*xi/L^3 + 6/L^2;
    H4pp = 6*xi/L^2 - 2/L;
    H1ppp = 12/L^3 * ones(size(xi));
    H2ppp = 6/L^2  * ones(size(xi));
    H3ppp = -12/L^3 * ones(size(xi));
    H4ppp = 6/L^2  * ones(size(xi));

    % FEM displacement and slope
    w_el  = H1*ue(1) + H2*ue(2) + H3*ue(3) + H4*ue(4);
    th_el = H1p.*ue(1) + H2p.*ue(2) + H3p.*ue(3) + H4p.*ue(4);

    % Add particular solution for element 1 under uniform load
    % v_p(xi) = -q/(24EI) * xi^2*(L-xi)^2  [upward-positive -> downward load -> negative]
    if e == 1
        vp    = -q_e/(24*EI) * xi.^2 .* (L - xi).^2;
        w_el  = w_el + vp;
        thp   = -q_e/(24*EI) * (2*xi.*(L-xi).^2 - 2*xi.^2.*(L-xi));
        th_el = th_el + thp;
    end

    % Bending moment: M = +EI * v''  (sagging positive convention)
    v2nd = H1pp*ue(1) + H2pp*ue(2) + H3pp*ue(3) + H4pp*ue(4);
    if e == 1
        % v_p'' = -q/(24EI) * (2L^2 - 12Lx + 12x^2)
        vp2nd = -q_e/(24*EI) * (2*L^2 - 12*L*xi + 12*xi.^2);
        v2nd  = v2nd + vp2nd;
    end
    My_el = EI * v2nd;

    % Shear force: Q = +EI * v'''  (positive = upward on left face)
    v3rd = H1ppp*ue(1) + H2ppp*ue(2) + H3ppp*ue(3) + H4ppp*ue(4);
    if e == 1
        % v_p''' = q/(2EI) * (L - 2x)  [derived from v_p = -q/(24EI)*x^2*(L-x)^2]
        vp3rd = q_e/(2*EI) * (L - 2*xi);
        v3rd  = v3rd + vp3rd;
    end
    Qz_el = EI * v3rd;

    x_el   = xn(e) + xi;
    w_all  = [w_all,  w_el];
    th_all = [th_all, th_el];
    My_all = [My_all, My_el];
    Qz_all = [Qz_all, Qz_el];
    x_all  = [x_all,  x_el];
end

%% --- Max values ---
fprintf('=== Internal Force Extremes ===\n');
[My_max, idx] = max(My_all);
fprintf('  M_max = %.4g N*mm at x = %.4g mm\n', My_max, x_all(idx));
[My_min, idx] = min(My_all);
fprintf('  M_min = %.4g N*mm at x = %.4g mm\n', My_min, x_all(idx));

sig_max = max(abs([My_max, My_min])) * z_max / Iy;
fprintf('\n=== Maximum Normal Stress ===\n');
fprintf('  sigma_max = M_max * z_max / Iy = %.4g N/mm^2 = %.4g MPa\n', sig_max, sig_max);

%% --- Figure ---
figure('Color','white','Position',[100 100 900 900]);

subplot(4,1,1);
plot(x_all, -w_all, 'b-', 'LineWidth', 2);  % plot downward deflection as positive
xlabel('x [mm]'); ylabel('w [mm] (downward +)');
title('Deflection'); grid on;
yline(0,'k--'); set(gca,'YDir','reverse');

subplot(4,1,2);
plot(x_all, th_all*1e3, 'b-', 'LineWidth', 2);
xlabel('x [mm]'); ylabel('\theta [mrad]');
title('Slope (rotation)'); grid on; yline(0,'k--');

subplot(4,1,3);
plot(x_all, My_all*1e-6, 'r-', 'LineWidth', 2);
xlabel('x [mm]'); ylabel('M [kN\cdotm]');
title('Bending moment'); grid on; yline(0,'k--');

subplot(4,1,4);
plot(x_all, Qz_all*1e-3, 'm-', 'LineWidth', 2);
xlabel('x [mm]'); ylabel('Q [kN]');
title('Shear force'); grid on; yline(0,'k--');

book_style(gcf);
exportgraphics(gcf,'../figures/beam_continuous_results.pdf','ContentType','vector');
exportgraphics(gcf,'../figures/beam_continuous_results.png','Resolution',300);
fprintf('\nFigures saved.\n');
