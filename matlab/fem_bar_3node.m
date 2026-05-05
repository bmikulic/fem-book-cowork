% Chapter 3 -- Three-node quadratic rod element
% Worked example (sec:ex_bar_quad): cantilever bar of length L, axial
% stiffness EA, uniform distributed load p0, single 3-node element.
% Demonstrates the 2-point Gauss-Legendre loop that produces ke and fe.

clear; clc; close all;

%% Parameters
L  = 1.0;       % length [m]
EA = 1.0;       % axial stiffness [N]
p0 = 1.0;       % uniform load per unit length [N/m]

%% Quadratic Lagrange shape functions on the parent domain xi in [-1, 1]
% Each handle returns a 3-by-1 column vector [N1; N2; N3] or its derivative
N    = @(xi) [0.5*xi*(xi - 1); 1 - xi*xi; 0.5*xi*(xi + 1)];
dNdx = @(xi) [xi - 0.5;        -2*xi;     xi + 0.5      ];

%% Two-point Gauss-Legendre rule on [-1, 1]
g  = 1/sqrt(3);
gp = [-g, +g];
w  = [1,  1 ];

%% Element matrices via Gauss quadrature
% Three nodes are equally spaced at x = 0, L/2, L, so the Jacobian
% J = dx/dxi = L/2 is constant; the loop is written generally so that the
% same code applies to a distorted interior node where J(xi) varies.
ke = zeros(3);
fe = zeros(3, 1);
for q = 1:2
    xi   = gp(q);
    Nv   = N(xi);
    dNv  = dNdx(xi);
    Jq   = dNv' * [0; L/2; L];      % geometry interpolation
    Binv = 1/Jq;
    B    = Binv * dNv';             % 1-by-3 strain-displacement row
    ke   = ke + w(q) * EA * (B' * B) * Jq;
    fe   = fe + w(q) * Nv * p0 * Jq;
end

fprintf('=== ELEMENT STIFFNESS (multiplied by 3L/EA) ===\n');
disp(ke * 3*L/EA);
fprintf('=== ELEMENT LOAD (multiplied by 6/(p0 L)) ===\n');
disp(fe * 6/(p0*L));

%% Apply boundary condition u1 = 0 and solve
free    = [2, 3];
u       = zeros(3, 1);
u(free) = ke(free, free) \ fe(free);

fprintf('\n=== NODAL DISPLACEMENTS (in units p0 L^2 / EA) ===\n');
fprintf('u1 = %+.4f   (exact:  0.0000)\n', u(1) * EA / (p0*L^2));
fprintf('u2 = %+.4f   (exact:  0.3750 = 3/8)\n', u(2) * EA / (p0*L^2));
fprintf('u3 = %+.4f   (exact:  0.5000 = 1/2)\n', u(3) * EA / (p0*L^2));

%% Strain at the three nodes (xi = -1, 0, +1)
xi_nodes = [-1, 0, +1];
strain   = zeros(3, 1);
for k = 1:3
    Bk = (2/L) * dNdx(xi_nodes(k))';
    strain(k) = Bk * u;
end
fprintf('\n=== NODAL STRAIN (in units p0 L / EA) ===\n');
fprintf('eps1 = %+.4f   (exact:  1.0000)\n', strain(1) * EA / (p0*L));
fprintf('eps2 = %+.4f   (exact:  0.5000)\n', strain(2) * EA / (p0*L));
fprintf('eps3 = %+.4f   (exact:  0.0000)\n', strain(3) * EA / (p0*L));

%% Reaction at the fixed support
R1 = ke(1, :) * u - fe(1);
fprintf('\n=== REACTION AT FIXED END ===\n');
fprintf('R1 = %+.4f p0 L   (exact: -1.0000)\n', R1 / (p0*L));

%% Plot FE vs analytical displacement field
nx       = 201;
xx       = linspace(0, L, nx);
xi_grid  = 2*xx/L - 1;
u_fe     = zeros(1, nx);
for k = 1:nx
    u_fe(k) = N(xi_grid(k))' * u;
end
u_exact = (p0/EA) * (L*xx - 0.5*xx.^2);

fig = figure('Units','centimeters','Position',[2, 2, 14, 9]);
plot(xx/L, u_exact * EA/(p0*L^2), 'k-',  'LineWidth', 1.8); hold on;
plot(xx/L, u_fe    * EA/(p0*L^2), 'b--', 'LineWidth', 1.8);
plot([0, 0.5, 1.0], u * EA/(p0*L^2), 'ro', ...
     'MarkerSize', 8, 'MarkerFaceColor', 'r');
grid on; box on;
xlabel('x / L');
ylabel('u EA / (p_0 L^2)');
legend({'Exact','FE field','FE nodes'}, 'Location','northwest');
title('Cantilever bar -- single 3-node quadratic element');

if ~exist('../figures', 'dir'), mkdir('../figures'); end
exportgraphics(fig, '../figures/bar_3node_results.pdf', 'ContentType','vector');
exportgraphics(fig, '../figures/bar_3node_results.png', 'Resolution', 300);
fprintf('\nFigure saved to ../figures/bar_3node_results.pdf and .png\n');
