n = 2;
w = [0 pi/32];
gam_l1 = -1;

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

ntf = tf(zpk([0.95 0.95], [0.5+1i*0.1 0.5-1i*0.1], 1, 1)); % Test TF.

b = fliplr(ntf.num{1}(2:end))'; % Transfer function numerator coefficients.
a = fliplr(ntf.den{1}(2:end))';% Transfer function denominator coefficients.

% LMI variables for l1/(*) minimization.
P = sdpvar(n, n, 'symmetric');
alpha = sdpvar(1);
mu = sdpvar(1);
nu = sdpvar(1);
aa = sdpvar(n, n, 'symmetric'); % Convex relaxation, aa is constrained to force it to equal a*a'.

A = [zeros(n-1, 1) eye(n-1); -a'];
A_c = [zeros(n-1, 1) eye(n-1); zeros(1, n)];
B = [zeros(n-1, 1); 1];
D = 1;

outerFactor = [A_c B; eye(n) zeros(n, 1)];
Theta = @(Q, P, a) kron((Phi_d+[0 0; 0 a]), P) + kron(Psi_d, Q);

ymLmi{1} = [alpha*P...
    zeros(n, 1)...
    b - a;...
    zeros(1, n)...
    (mu - 1)...
    D';...
    (b - a)'...
    D...
    nu];
ymLmi{2} = outerFactor'*Theta(zeros(n), P, alpha)*outerFactor - [a*a' a; a' 1];
ymLmi{3} = [gam_l1 mu nu; mu 1 0; nu 0 1];
ymConstraint = [ymLmi{1}>=0, ymLmi{2}<=0, ymLmi{3}>=0, P>=0, nu>=0, mu>=0]; % l1/(*) constraints.

ymSolver = 'sdpt3'; %'SeDuMi';
ymOptions = sdpsettings('verbose', true, 'solver', ymSolver, 'sdpt3.inftol', 1e-9, 'sdpt3.gaptol', 1e-9, 'sdpt3.maxit', 200);
ymObjective = gam_l1;

ymLineSearch = optimizer(ymConstraint, ymObjective, ymOptions, alpha, ymObjective);
lineSearchOptions = optimset('Display', 'iter', 'FunValCheck', 'on');
[alphaOpt, ymLineSearchOpt] = fminbnd(@(x) ymLineSearch(x), 0, 1-(max(abs(eig(A))))^2, lineSearchOptions);
ymConstraint = replace(ymConstraint, alpha, alphaOpt);
            
ymDiagnostics = optimize(ymConstraint, ymObjective, ymOptions);
if ymDiagnostics.problem==1
    error('dtsyn:solver_infeasible', ymDiagnostics.info);
elseif ymDiagnostics.problem~=0
    warning('on');
    warning(ymDiagnostics.info);
end

gam_l1 = sqrt(value(gam_l1))

A = [zeros(n-1, 1) eye(n-1); -a'];
C = (b - a)';
D = 1;
sys_ss = ss(A, B, C, D, 1);
sys_tf = tf([1 fliplr(b')], [1 fliplr(a')], 1);