% pvw_linear_load.m
% Worked example for the Principle of Virtual Work (Chapter 3)
%
% Cantilever bar of length L, axial stiffness EA, fixed at x = 0,
% loaded by a linearly increasing distributed force p(x) = p0*x/L.
%
% Apply PVW with the linear one-parameter trial   u(x) = a*x
% and compare to the exact solution
%   u_exact(x) = p0 * x * (3L^2 - x^2) / (6 * EA * L)

clear; clc; close all;

%% Parameters (consistent units; values chosen for visualization)
L  = 1.0;          % bar length
EA = 1.0;          % axial stiffness  (so u is in units of p0)
p0 = 1.0;          % peak distributed load at x = L

%% PVW one-term solution
%   delta W_int = EA * L * a * delta_a
%   delta W_ext = (p0 * L^2 / 3) * delta_a
%   ==>  a = p0*L / (3*EA)
a_pvw  = p0 * L / (3*EA);
u_pvw  = @(x) a_pvw .* x;

%% Exact solution
u_exact = @(x) p0 .* x .* (3*L^2 - x.^2) ./ (6*EA*L);

%% Tip values
fprintf('PVW tip displacement   u(L) = %.6f\n', u_pvw(L));
fprintf('Exact tip displacement u(L) = %.6f\n', u_exact(L));
fprintf('Ratio (should be 1)         = %.6f\n', u_pvw(L)/u_exact(L));

%% Plot
xx = linspace(0, L, 400);
fig = figure('Name','PVW vs Exact', 'Color','w', ...
             'Units','inches','Position',[1 1 7 4]);
plot(xx, u_exact(xx), 'k-',  'LineWidth', 1.8); hold on;
plot(xx, u_pvw(xx),   'b--', 'LineWidth', 1.5);
plot([0, L], [u_pvw(0), u_pvw(L)], 'ro', ...
     'MarkerSize', 8, 'MarkerFaceColor','y');
grid on; box on;
xlabel('Normalized coordinate, x/L');
ylabel('Normalized displacement, u \cdot EA / (p_0 L^2)');
legend({'Exact', 'PVW one-term linear trial', 'Endpoints'}, ...
       'Location', 'northwest');
title('PVW with a linear trial vs. exact solution');

%% Save
thisFile = mfilename('fullpath');
[here,~,~] = fileparts(thisFile);
figDir = fullfile(here, '..', 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end
exportgraphics(fig, fullfile(figDir, 'pvw_linear_load.pdf'), ...
               'ContentType','vector');
exportgraphics(fig, fullfile(figDir, 'pvw_linear_load.png'), ...
               'Resolution', 300);
fprintf('Figure saved to %s\n', figDir);
