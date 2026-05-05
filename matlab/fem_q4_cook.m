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
