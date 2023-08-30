n = 4;
w = [0 pi/32];
gam_gkyp = -1;
gam_hinf = 1.5;

if gam_gkyp==-1
    gam_gkyp = sdpvar(1); % Optimization problem.
else
    gam_gkyp = gam_gkyp.^2; % Sub-optimal (feasibility) problem.
end
if gam_hinf==-1
    gam_hinf = sdpvar(1); % Optimization problem.
else
    gam_hinf = gam_hinf.^2; % Sub-optimal (feasibility) problem.
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
adot = sdpvar(n, 1, 'full');
P{1} = sdpvar(n, n, 'symmetric');
P{2} = sdpvar(n, n, 'symmetric');
Q_script = sdpvar(n, n, 'symmetric');

A_c = [zeros(n-1, 1) eye(n-1); zeros(1, n)];
B = [zeros(n-1, 1); 1];

outerFactor = [A_c B; eye(n) zeros(n, 1)];
Theta = @(Q, P) kron(Phi_d, P) + kron(Psi_d, Q);

ymLmi{1} = [outerFactor'*Theta(Q_script, P{1})*outerFactor - [adot*a' + a*adot' a; a' 1]...
    [b adot; 1 0];...
    [b' 1; adot' 0]...
    [-gam_gkyp 0; 0 -1]];
ymLmi{2} = [outerFactor'*Theta(zeros(n), P{2})*outerFactor - [adot*a' + a*adot' a; a' 1]...
    [b adot; 1 0];...
    [b' 1; adot' 0]...
    [-gam_hinf 0; 0 -1]];
ymConstraint = [ymLmi{1}<=0, ymLmi{2}<=0, Q_script>=0, P{2}>=0];

ymSolver = 'mosek'; %'SeDuMi';
%ymOptions = sdpsettings('verbose', true, 'solver', ymSolver, 'lmilab.reltol', 1e-6, 'lmilab.maxiter', 500, 'lmilab.feasradius', 1e9, 'lmilab.L', 100);
ymOptions = sdpsettings('verbose', true, 'solver', ymSolver, 'sdpt3.inftol', 1e-9, 'sdpt3.gaptol', 1e-9, 'sdpt3.maxit', 500);
ymObjective = gam_gkyp + gam_hinf;

ymOptimizeA = optimizer(ymConstraint, ymObjective, ymOptions, adot, a);
ymOptimizeAdot = optimizer(ymConstraint, ymObjective, ymOptions, a, [adot; gam_gkyp; gam_hinf]);

epsilon = 1e-9;
maxIter = 10000;
%ntf = tf(synthesizeNTF(4, 32, 3));
%a_val = fliplr(ntf.den{1}(2:end))';
a_val = [0 0 0 0]';
%a_val = [0.2401;-1.36374876879259;2.91645060693392;-2.78316075263793];
%a_val = [0.444444297096043;-2.13578877231391;3.89195992817283;-3.19447366909881]; % synthesizeNTF(4, 32, 0)
objective_val(1:2, :) = 1e6;
relTol = 1e6;
k = 1;

while any(relTol>epsilon) && k<maxIter
    ymResult = ymOptimizeAdot(ymOptimizeA(a_val(:, k)));
    k = k + 1;
    a_val(:, k) = ymResult(1:n);
    objective_val(:, k) = ymResult(n+1:n+2);
    relTol = abs(objective_val(:, k) - objective_val(:, k-1))./abs(objective_val(:, k-1));
end

ymSolver = 'lmilab'; %'SeDuMi';
ymOptions = sdpsettings('verbose', true, 'solver', ymSolver, 'lmilab.reltol', 1e-6, 'lmilab.maxiter', 500, 'lmilab.feasradius', 1e9, 'lmilab.L', 100);

ymConstraint = replace(ymConstraint, adot, a_val(:, end));
ymDiagnostics = optimize(ymConstraint, ymObjective, ymOptions);
if ymDiagnostics.problem==1
    error('dtsyn:solver_infeasible', ymDiagnostics.info);
elseif ymDiagnostics.problem~=0
    warning('on');
    warning(ymDiagnostics.info);
end

a = value(a);
b = value(b);
gam_gkyp = sqrt(value(gam_gkyp));
gam_hinf = sqrt(value(gam_hinf));

A = [zeros(n-1, 1) eye(n-1); -a'];
C = (b - a)';
D = 1;
sys_ss = ss(A, B, C, D, 1);
sys_tf = tf([1 fliplr(b')], [1 fliplr(a')], 1);