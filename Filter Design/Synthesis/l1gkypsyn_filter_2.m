function [sys_tf, gam_gkyp, gam_l1] = l1gkypsyn_filter_2(order, w, gam_gkyp, gam_l1)

n = order;

if gam_gkyp==-1
    gam_gkyp = sdpvar(1); % Optimization problem.
else
    gam_gkyp = gam_gkyp.^2; % Sub-optimal (feasibility) problem.
end
if gam_l1==-1
    gam_l1 = sdpvar(1); % Optimization problem.
else
    gam_l1 = gam_l1.^2; % Sub-optimal (feasibility) problem.
end

if w(1)==0 % Low frequency case, Psi_d = [-2*cos(w(1)) 1; 1 0]
    w = [-w(2) w(2)];
elseif w(2)==pi % High  frequency case, Psi_d = [2*cos(w(2)) -1; -1 0]
    w = [w(1) 2*pi-w(1)];
end

wc = (w(2) + w(1))/2;
w0 = (w(2) - w(1))/2;
Psi_d = [0 exp(1i*wc); exp(-1i*wc) -2*cos(w0)];
Phi_d = [1 0; 0 -1];

b = sdpvar(n, 1, 'full'); % Transfer function numerator coefficients.
a = sdpvar(n, 1, 'full'); % Transfer function denominator coefficients.

% LMI variables for H-Inf (GKYP) minimization.
P{1} = sdpvar(n, n, 'symmetric');
Q_script = sdpvar(n, n, 'symmetric');
aa = sdpvar(n, n, 'symmetric'); % Convex relaxation, aa is constrained to force it to equal a*a'.

% LMI variables for l1/(*) minimization.
alpha = sdpvar(1);
P{2} = sdpvar(n, n, 'symmetric');
mu = sdpvar(1);
nu = sdpvar(1);

A_c = [zeros(n-1, 1) eye(n-1); zeros(1, n)];
B = [zeros(n-1, 1); 1];
D = 1;

outerFactor = [A_c B; eye(n) zeros(n, 1)];
Theta = @(Q, P, a) kron((Phi_d+[0 0; 0 a]), P) + kron(Psi_d, Q);

ymLmi{1} = [outerFactor'*Theta(Q_script, P{1}, 0)*outerFactor - [aa a; a' 1]...
    [b; D];...
    b'...
    D' ...
    -gam_gkyp];
ymLmi{2} = [alpha*P{2}...
    zeros(n, 1)...
    b - a;...
    zeros(1, n)...
    (mu - 1)...
    D';...
    (b - a)'...
    D...
    nu];
ymLmi{3} = outerFactor'*Theta(zeros(n), P{2}, 1-alpha)*outerFactor - [aa a; a' 1];
ymLmi{4} = [gam_l1 mu nu; mu 1 0; nu 0 1];
ymLmi{5} = [aa a; a' 1]; % Shor relaxation to try to force aa==(a*a').
ymConstraint = [ymLmi{1}<=0, Q_script>=0]; % H-Inf (GKYP) constraints.
ymConstraint = [ymConstraint, ymLmi{2}>=0, ymLmi{3}<=0, ymLmi{4}>=0, P{2}>=0, nu>=0, mu>=0]; % l1/(*) constraints.
ymConstraint = [ymConstraint, ymLmi{5}>=0]; % Convex relaxation constraint.

ymSolver = 'sdpt3';%'SeDuMi';
ymOptions = sdpsettings('verbose', true, 'solver', ymSolver, 'sdpt3.inftol', 1e-12, 'sdpt3.gaptol', 1e-12, 'sdpt3.maxit', 1000, 'debug', 1);
%ymObjective = gam_l1 + gam_gkyp + trace(aa);
ymObjective = gam_l1 + gam_gkyp% + norm(ymLmi{5}, 'nuclear');

ymLineSearch = optimizer(ymConstraint, ymObjective, ymOptions, alpha, ymObjective);
%               % Produce plot of alpha vs. optimized value.
%               figure;
%               plot(0:0.001:1, ymLineSearch(0:0.001:1));
lineSearchOptions = optimset('Display', 'iter', 'FunValCheck', 'on');
[alphaOpt, ymLineSearchOpt] = fminbnd(@(x) ymLineSearch(x), 0, 1, lineSearchOptions);
ymConstraint = replace(ymConstraint, alpha, alphaOpt);

ymDiagnostics = optimize(ymConstraint, ymObjective, ymOptions);
if ymDiagnostics.problem==1
    error('dtsyn:solver_infeasible', ymDiagnostics.info);
elseif ymDiagnostics.problem~=0
    warning('on');
    warning(ymDiagnostics.info);
end

aa = value(aa);
a = value(a);
b = value(b);
gam_gkyp = sqrt(value(gam_gkyp));
mag2db(gam_gkyp)
gam_l1 = sqrt(value(gam_l1));
gam_l1

A = [zeros(n-1, 1) eye(n-1); -a'];
C = (b - a)';
D = 1;
sys_ss = ss(A, B, C, D, 1);
sys_tf = tf([1 fliplr(b')], [1 fliplr(a')], 1);

end