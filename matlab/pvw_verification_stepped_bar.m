% pvw_verification_stepped_bar.m
% Numerically verify the principle of virtual work for the SOLVED
% stepped bar (the worked example of Section "Stepped Bar").
% For three different kinematically admissible virtual displacements,
% the script computes the internal and external virtual work separately
% and shows they are equal.

clear; clc;

%% Stepped bar parameters
E  = 210000;
A  = 225;
EA1 = E*A;        % element 1 axial stiffness [N]
EA2 = E*2*A;      % element 2 axial stiffness [N]
L1  = 500;        % element 1 length [mm]
L2  = 400;        % element 2 length [mm]
p1  = -4;         % element 1 distributed load [N/mm]
F2  =  6000;      % concentrated force at node 2 [N]
F3  = -15000;     % concentrated force at node 3 [N]

%% Exact (FEM) nodal displacements from the Stepped Bar example
u1 = 0;            % fixed
u2 = -20/189;      % mm   (= -0.10582)
u3 = -32/189;      % mm   (= -0.16931)

%% Internal axial force in each element using the enriched trial:
%   N(x) = EA*(u_b - u_a)/L  +  p*(L - 2x_local)/2
N1 = @(x) EA1*(u2-u1)/L1 + p1*(L1 - 2*x)/2;        % x measured from node 1
N2 = @(x) EA2*(u3-u2)/L2 + 0*x;                    % x measured from node 2

%% Three admissible virtual displacements (each must satisfy du(x=0)=0)
%   Each entry stores du(x), du'(x), and the nodal values du(node_k)
testCases = struct( ...
  'name',     {}, ...
  'du',       {}, ...
  'du_dx',    {}, ...
  'du_n2',    {}, ...
  'du_n3',    {});

% Case A — tent: du_2 = 1, du_3 = 0
testCases(end+1) = struct( ...
  'name',  'Case A:  (du_1, du_2, du_3) = (0, 1, 0)', ...
  'du',    @(x_local, elem) (elem==1).*(x_local/L1) + (elem==2).*(1 - x_local/L2), ...
  'du_dx', @(x_local, elem) (elem==1)*(1/L1)         + (elem==2)*(-1/L2), ...
  'du_n2', 1, ...
  'du_n3', 0);

% Case B — ramp + plateau: du_2 = 1, du_3 = 1
testCases(end+1) = struct( ...
  'name',  'Case B:  (du_1, du_2, du_3) = (0, 1, 1)', ...
  'du',    @(x_local, elem) (elem==1).*(x_local/L1) + (elem==2).*1, ...
  'du_dx', @(x_local, elem) (elem==1)*(1/L1)         + (elem==2)*0, ...
  'du_n2', 1, ...
  'du_n3', 1);

% Case C — global quadratic du(x_glob) = (x_glob / Ltot)^2
Ltot = L1 + L2;
testCases(end+1) = struct( ...
  'name',  'Case C:  du(x) = (x / 900)^2  (smooth, non-piecewise-linear)', ...
  'du',    @(x_local, elem) ((elem==1).*x_local + (elem==2).*(L1 + x_local)).^2 / Ltot^2, ...
  'du_dx', @(x_local, elem) 2*((elem==1).*x_local + (elem==2).*(L1 + x_local)) / Ltot^2, ...
  'du_n2', (L1/Ltot)^2, ...
  'du_n3', 1);

%% Verify PVW for each case
fprintf('Principle of Virtual Work verification — stepped bar\n');
fprintf('Solved nodal displacements:  u1 = 0,  u2 = %.5f,  u3 = %.5f mm\n\n', u2, u3);
fprintf('%-58s %14s %14s %14s\n', 'Admissible virtual displacement', ...
        'dW_int [Nmm]', 'dW_ext [Nmm]', '|diff|');
fprintf('%s\n', repmat('-', 1, 102));

for k = 1:length(testCases)
    tc = testCases(k);

    % Internal virtual work — sum element contributions
    dWint_e1 = integral(@(x) N1(x).*tc.du_dx(x,1), 0, L1);
    dWint_e2 = integral(@(x) N2(x).*tc.du_dx(x,2), 0, L2);
    dWint    = dWint_e1 + dWint_e2;

    % External virtual work — distributed load (element 1 only) and
    % concentrated forces; reaction at node 1 contributes 0 because du(0)=0.
    dWext_d  = integral(@(x) p1.*tc.du(x,1), 0, L1);
    dWext_c  = F2*tc.du_n2 + F3*tc.du_n3;
    dWext    = dWext_d + dWext_c;

    fprintf('%-58s %14.5f %14.5f %14.2e\n', tc.name, dWint, dWext, abs(dWint-dWext));
end

fprintf(['\nFor every kinematically admissible virtual displacement,\n', ...
         'dW_int = dW_ext to numerical precision: PVW holds.\n']);
