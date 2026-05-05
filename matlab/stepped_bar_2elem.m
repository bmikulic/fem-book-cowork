% stepped_bar_2elem.m
% Worked example: stepped axial bar with two bar elements
% (Chapter 3 — One-Dimensional Bar and Rod Elements)
%
% Geometry: cantilever bar with two segments
%   Element 1: length 500 mm, cross-section area  A
%   Element 2: length 400 mm, cross-section area 2A
% Loading:
%   Distributed axial load p1 = -4 N/mm on element 1 only
%   Concentrated forces F2 = +6000 N at node 2, F3 = -15000 N at node 3
% Boundary condition: node 1 fully fixed (u1 = 0)

clear; clc; close all;

%% Material and geometry
E  = 210000;        % Young's modulus [MPa = N/mm^2]
A  = 225;           % cross-section area, element 1 [mm^2]

%% Distributed load (acts only on element 1)
p_e1 = -4;          % N/mm  (negative = directed toward fixed end)
p_e2 =  0;          % N/mm

%% Concentrated forces
P_node2 =  6000;    % +x at node 2 [N]
P_node3 = -15000;   % -x at node 3 [N]

%% Discretization
nodes_x   = [0, 500, 900];    % node coordinates [mm]
nElem     = 2;
nNode     = numel(nodes_x);
nDof      = nNode;            % 1D bar: 1 DOF per node
elem_conn = [1 2; 2 3];       % element connectivity
elem_area = [A, 2*A];         % cross-section per element
elem_p    = [p_e1, p_e2];     % distributed load per element
elem_L    = diff(nodes_x);    % element lengths

%% Initialize global arrays
KG    = zeros(nDof, nDof);    % global stiffness matrix
FG_p  = zeros(nDof, 1);       % consistent (distributed) load vector
FG_c  = zeros(nDof, 1);       % concentrated nodal forces

%% Element-by-element assembly
for e = 1:nElem
    L    = elem_L(e);
    Ae   = elem_area(e);
    pe   = elem_p(e);
    dofs = elem_conn(e,:);

    % Element stiffness matrix
    ke = (E*Ae/L) * [1 -1; -1 1];

    % Consistent distributed-load vector for uniform p
    fe_p = (pe*L/2) * [1; 1];

    % Add to global arrays at the element DOFs
    KG(dofs, dofs) = KG(dofs, dofs) + ke;
    FG_p(dofs)     = FG_p(dofs)     + fe_p;

    fprintf('Element %d:  L = %g mm, A = %g mm^2, p = %g N/mm\n', ...
            e, L, Ae, pe);
    fprintf('  k^e = \n'); disp(ke);
    fprintf('  f_p^e = \n'); disp(fe_p);
end

%% Concentrated nodal forces
FG_c(2) = P_node2;
FG_c(3) = P_node3;

%% Total global force vector
FG = FG_p + FG_c;

fprintf('\nGlobal stiffness matrix KG [N/mm]:\n'); disp(KG);
fprintf('Distributed-load vector FG_p [N]:\n');    disp(FG_p);
fprintf('Concentrated-force vector FG_c [N]:\n');  disp(FG_c);
fprintf('Total force vector FG [N]:\n');           disp(FG);

%% Apply boundary conditions
fixed_dofs = 1;                                % node 1 fixed
free_dofs  = setdiff(1:nDof, fixed_dofs);

%% Solve reduced system
UG             = zeros(nDof, 1);
UG(free_dofs)  = KG(free_dofs, free_dofs) \ FG(free_dofs);

%% Compute reactions
RG = KG*UG - FG;

fprintf('Nodal displacements UG [mm]:\n'); disp(UG);
fprintf('Reactions RG [N]:\n');           disp(RG);

%% Postprocess: displacement and internal axial force inside each element
%  Using the enriched trial: u(x) = N1*u1 + N2*u2 + p*x*(L-x)/(2EA)
%  Internal force:           N(x) = EA*(u2-u1)/L + p*(L-2x)/2

%% Element 1
L1 = elem_L(1);  Ae1 = elem_area(1);  pe1 = elem_p(1);
u1 = UG(1);  u2 = UG(2);
x1 = linspace(0, L1, 200);
xi1 = x1 / L1;
u_e1 = (1-xi1)*u1 + xi1*u2 + pe1*x1.*(L1 - x1)/(2*E*Ae1);
N_e1 = E*Ae1*(u2-u1)/L1 + pe1*(L1 - 2*x1)/2;
xg1  = nodes_x(1) + x1;

%% Element 2
L2 = elem_L(2);  Ae2 = elem_area(2);  pe2 = elem_p(2);
u1 = UG(2);  u2 = UG(3);
x2 = linspace(0, L2, 200);
xi2 = x2 / L2;
u_e2 = (1-xi2)*u1 + xi2*u2 + pe2*x2.*(L2 - x2)/(2*E*Ae2);
N_e2 = E*Ae2*(u2-u1)/L2 + pe2*(L2 - 2*x2)/2;
xg2  = nodes_x(2) + x2;

%% Plot displacement field
fig1 = figure('Name','Stepped Bar — Displacement','Color','w', ...
              'Units','inches','Position',[1 1 7 4]);
plot(xg1, u_e1, 'b-', 'LineWidth', 1.7); hold on;
plot(xg2, u_e2, 'r-', 'LineWidth', 1.7);
plot(nodes_x, UG, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'y');
grid on; box on;
xlabel('Axial coordinate, x [mm]');
ylabel('Axial displacement, u(x) [mm]');
legend('Element 1', 'Element 2', 'Nodal values', 'Location','southwest');
title('Displacement field along the stepped bar');

%% Plot internal axial force
fig2 = figure('Name','Stepped Bar — Internal Force','Color','w', ...
              'Units','inches','Position',[1 1 7 4]);
plot(xg1, N_e1, 'b-', 'LineWidth', 1.7); hold on;
plot(xg2, N_e2, 'r-', 'LineWidth', 1.7);
% Visualize the jump at node 2 due to concentrated force
plot([nodes_x(2), nodes_x(2)], [N_e1(end), N_e2(1)], 'k--', 'LineWidth', 1.0);
grid on; box on;
xlabel('Axial coordinate, x [mm]');
ylabel('Internal axial force, N(x) [N]');
legend('Element 1', 'Element 2', 'Jump at node 2', 'Location','best');
title('Internal axial force along the stepped bar');

%% Save the figures into ../figures/
thisFile = mfilename('fullpath');
[here,~,~] = fileparts(thisFile);
figDir = fullfile(here, '..', 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end
book_style(fig1);
exportgraphics(fig1, fullfile(figDir, 'stepped_bar_displacement.pdf'), ...
               'ContentType','vector');
exportgraphics(fig1, fullfile(figDir, 'stepped_bar_displacement.png'), ...
               'Resolution', 300);
book_style(fig2);
exportgraphics(fig2, fullfile(figDir, 'stepped_bar_internal_force.pdf'), ...
               'ContentType','vector');
exportgraphics(fig2, fullfile(figDir, 'stepped_bar_internal_force.png'), ...
               'Resolution', 300);

fprintf('\nFigures saved to %s\n', figDir);
