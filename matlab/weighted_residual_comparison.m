% weighted_residual_comparison.m
% One-term weighted residual approximation for an axially loaded bar
% (Chapter 2 — Weighted Residual Methods)

clear; clc; close all;

%% Problem parameters
L  = 1.0;        % bar length [m]
E  = 210e9;      % Young's modulus [Pa]
A  = 1.0e-4;     % cross-sectional area [m^2]
p0 = 1.0e4;      % load intensity scale [N/m]
EA = E*A;        % axial stiffness [N]

%% Method coefficients for u_tilde = a*(xi^2 - 2*xi)
a_col = -p0*L^2/(8*EA);       % collocation at xi = 1/2
a_sub = -p0*L^2/(6*EA);       % one subdomain over [0,L]
a_gal = -9*p0*L^2/(40*EA);    % Galerkin weighting

methodNames = {'Collocation'; 'Subdomain'; 'Galerkin'};
aValues = [a_col; a_sub; a_gal];
tipValues = -aValues;         % xi=1 gives u_tilde(L) = -a

%% Exact tip displacement and percent error
uTipExact = p0*L^2/(4*EA);
percentError = abs(tipValues - uTipExact)/uTipExact * 100;

fprintf('Weighted residual comparison for axial bar\n');
fprintf('L = %.3f m, EA = %.4e N, p0 = %.4e N/m\n\n', L, EA, p0);
fprintf('%-14s %14s %18s %12s\n', ...
        'Method', 'a [m]', 'u_tilde(L) [m]', 'Error [%]');
fprintf('%s\n', repmat('-', 1, 62));

for i = 1:numel(methodNames)
    fprintf('%-14s %14.6e %18.6e %12.3f\n', ...
        methodNames{i}, aValues(i), tipValues(i), percentError(i));
end
fprintf('%-14s %14s %18.6e %12.3f\n', 'Exact', '--', uTipExact, 0);

%% Normalized displacement curves
xi = linspace(0, 1, 400);
phi = xi.^2 - 2*xi;

u_col_norm = -(1/8)  * phi;
u_sub_norm = -(1/6)  * phi;
u_gal_norm = -(9/40) * phi;
u_exact_norm = xi/3 - xi.^4/12;

figFile = figure('Name','Weighted Residual Comparison','Color','w', ...
                 'Units','inches','Position',[1 1 7 4.5]);
plot(xi, u_exact_norm, 'k-',  'LineWidth', 1.8); hold on;
plot(xi, u_col_norm,   'b--', 'LineWidth', 1.3);
plot(xi, u_sub_norm,   'r-.', 'LineWidth', 1.3);
plot(xi, u_gal_norm,   'g:',  'LineWidth', 1.9);
grid on; box on;
xlabel('Normalized coordinate, \xi = x/L');
ylabel('Normalized displacement, uEA/(p_0L^2)');
legend('Exact', 'Collocation', 'Subdomain', 'Galerkin', ...
       'Location', 'northwest');

%% Save the figure into ../figures/  (PDF for LaTeX, PNG for previewing)
thisFile = mfilename('fullpath');
[here, ~, ~] = fileparts(thisFile);
figDir = fullfile(here, '..', 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

set(figFile, 'PaperPositionMode', 'auto');
exportgraphics(figFile, fullfile(figDir, 'weighted_residual_comparison.pdf'), ...
               'ContentType', 'vector');
exportgraphics(figFile, fullfile(figDir, 'weighted_residual_comparison.png'), ...
               'Resolution', 300);

fprintf('\nFigure saved to %s\n', figDir);
