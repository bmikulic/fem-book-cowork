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

book_style(fig);
exportgraphics(fig,'../figures/cst_example_results.pdf','ContentType','vector')
exportgraphics(fig,'../figures/cst_example_results.png','Resolution',300)
fprintf('\nFigures saved to ../figures/\n')
