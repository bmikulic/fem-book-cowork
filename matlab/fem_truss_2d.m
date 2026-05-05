% fem_truss_2d.m -- Chapter 4: Worked Example, 3-element 2D truss
% Nodes: 1=(0,0), 2=(1000,800), 3=(0,800) [mm]
% BCs:   node 1 pinned (U1=V1=0), node 3 pinned (U3=V3=0)
% Load:  Fy = -10000 N (downward) at node 2
% E = 210000 N/mm^2, A = 225 mm^2
clear; clc; close all;

%% --- Geometry and material ---
coords = [   0    0;    % node 1
           1000  800;   % node 2
              0  800];  % node 3
conn = [1 2;            % element 1: nodes 1->2 (diagonal)
        1 3;            % element 2: nodes 1->3 (vertical)
        3 2];           % element 3: nodes 3->2 (horizontal)
E = 210000;  % N/mm^2
A = 225;     % mm^2
nEl = 3;  nNodes = 3;  nDOF = 2*nNodes;

%% --- Element properties ---
fprintf('=== Element Properties ===\n');
fprintf('%-4s %-8s %-8s %-8s %-8s %-12s\n', ...
        'Elem', 'L [mm]', 'theta', 'c', 's', 'EA/L [N/mm]');
for e = 1:nEl
    n1 = conn(e,1); n2 = conn(e,2);
    dx = coords(n2,1) - coords(n1,1);
    dy = coords(n2,2) - coords(n1,2);
    Le(e)  = sqrt(dx^2 + dy^2);
    ce(e)  = dx / Le(e);
    se(e)  = dy / Le(e);
    EAL(e) = E * A / Le(e);
    theta  = atan2d(dy, dx);
    fprintf('%-4d %-8.4f %-8.2f %-8.4f %-8.4f %-12.4f\n', ...
            e, Le(e), theta, ce(e), se(e), EAL(e));
end

%% --- Element global stiffness matrices ---
K = zeros(nDOF);
for e = 1:nEl
    c = ce(e); s = se(e); k = EAL(e);
    ke = k * [ c^2   c*s  -c^2  -c*s;
               c*s   s^2  -c*s  -s^2;
              -c^2  -c*s   c^2   c*s;
              -c*s  -s^2   c*s   s^2];
    n1 = conn(e,1); n2 = conn(e,2);
    dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2];
    K(dofs, dofs) = K(dofs, dofs) + ke;
    fprintf('\nElement %d global stiffness [DOFs %d %d %d %d]:\n', ...
            e, dofs(1), dofs(2), dofs(3), dofs(4));
    disp(round(ke));
end

fprintf('\n=== Assembled 6x6 Global Stiffness Matrix K [N/mm] ===\n');
disp(round(K));

%% --- Loads and BCs ---
F = zeros(nDOF, 1);
F(4) = -10000;            % Fy at node 2 (DOF 4)

fixedDOFs = [1 2 5 6];   % node 1 and 3 pinned
freeDOFs  = [3 4];        % node 2 free (U2, V2)

%% --- Solve reduced system ---
Kff = K(freeDOFs, freeDOFs);
Ff  = F(freeDOFs);

fprintf('\n=== Reduced 2x2 Stiffness Matrix Kff [N/mm] ===\n');
disp(round(Kff));
fprintf('Load vector Ff [N]: [%.1f; %.1f]\n', Ff(1), Ff(2));

u_free = Kff \ Ff;
u = zeros(nDOF, 1);
u(freeDOFs) = u_free;

fprintf('\n=== Nodal Displacements ===\n');
fprintf('U2 = %+.6f mm\n', u(3));
fprintf('V2 = %+.6f mm\n', u(4));

%% --- Element axial forces ---
fprintf('\n=== Element Axial Forces ===\n');
N = zeros(1, nEl);
for e = 1:nEl
    c = ce(e); s = se(e);
    n1 = conn(e,1); n2 = conn(e,2);
    dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2];
    delta = [-c -s c s] * u(dofs);
    N(e) = EAL(e) * delta;
    if N(e) > 1e-6
        type = 'tension';
    elseif N(e) < -1e-6
        type = 'compression';
    else
        type = 'zero';
    end
    fprintf('N%d = %+.4f N  (%s)\n', e, N(e), type);
end

%% --- Reactions ---
R = K * u - F;
fprintf('\n=== Support Reactions ===\n');
fprintf('R1x = %+.4f N\n', R(1));
fprintf('R1y = %+.4f N\n', R(2));
fprintf('R3x = %+.4f N\n', R(5));
fprintf('R3y = %+.4f N\n', R(6));

fprintf('\n=== Global Equilibrium Check ===\n');
Fext = F;
fprintf('SumFx = %+.8f N\n', sum(R(1:2:end)) + sum(Fext(1:2:end)));
fprintf('SumFy = %+.8f N\n', sum(R(2:2:end)) + sum(Fext(2:2:end)));

%% --- Deformed shape figure ---
scale = 500;
figure('Color','white');
hold on; axis equal; box on;
colors = {'b','r','g'};
for e = 1:nEl
    n1 = conn(e,1); n2 = conn(e,2);
    xo = [coords(n1,1) coords(n2,1)];
    yo = [coords(n1,2) coords(n2,2)];
    xd = xo + scale * [u(2*n1-1) u(2*n2-1)]';
    yd = yo + scale * [u(2*n1  ) u(2*n2  )]';
    plot(xo, yo, 'k--', 'LineWidth', 1.2);
    plot(xd, yd, '-', 'Color', colors{e}, 'LineWidth', 2.5);
end
% Nodes undeformed
scatter(coords(:,1), coords(:,2), 60, 'ko', 'filled');
% Nodes deformed
xu = coords(:,1) + scale * u(1:2:end);
yu = coords(:,2) + scale * u(2:2:end);
scatter(xu, yu, 60, 'b^', 'filled');
xlabel('x [mm]', 'FontSize', 12);
ylabel('y [mm]', 'FontSize', 12);
title(sprintf('Deformed shape (magnification \\times %d)', scale), 'FontSize', 12);
legend('Undeformed', 'Elem 1 def.', 'Elem 2 def.', 'Elem 3 def.', ...
       'Nodes orig.', 'Nodes def.', 'Location', 'best');
book_style(gcf);
exportgraphics(gcf, '../figures/truss2d_deformed.pdf', 'ContentType', 'vector');
exportgraphics(gcf, '../figures/truss2d_deformed.png', 'Resolution', 300);
fprintf('\nFigures saved to ../figures/truss2d_deformed.pdf/.png\n');
