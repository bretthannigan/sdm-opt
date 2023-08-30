n = 5;
w = [0 pi/64];
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
Psi_d = [-1 1i*wc; -1i*wc -w(1)*w(2)];
Phi_d = [0 1; 1 0];

b = sdpvar(n, 1, 'full'); % Transfer function numerator coefficients.
a = sdpvar(n, 1, 'full'); % Transfer function denominator coefficients.
test = sdpvar(1);
P{1} = sdpvar(n, n, 'symmetric');
P{2} = sdpvar(n, n, 'symmetric');
Q_script = sdpvar(n, n, 'symmetric');
aa = sdpvar(n, n, 'symmetric'); % Convex relaxation, aa is constrained to force it to equal a*a'.

A_c = [zeros(n-1, 1) eye(n-1); zeros(1, n)];
B = [zeros(n-1, 1); 1];

outerFactor = [A_c B; eye(n) zeros(n, 1)];
Theta = @(Q, P) kron(Phi_d, P) + kron(Psi_d, Q);

ymLmi{1} = [outerFactor'*Theta(Q_script, P{1})*outerFactor - [aa a; a' 1]...
    [b; 1];...
    b'...
    1 ...
    -gam_gkyp];
ymLmi{2} = [outerFactor'*Theta(zeros(n), P{2})*outerFactor - [aa a; a' 1]...
    [b; 1];...
    b'...
    1 ...
    -gam_hinf];
ymLmi{3} = [aa a; a' 1]; % Shor relaxation to try to force aa==(a*a').
ymConstraint = [ymLmi{1}<=0, ymLmi{2}<=0, Q_script>=0, P{2}>=0, ymLmi{3}>=0];

ymSolver = 'sdpt3'; %'SeDuMi';
ymOptions = sdpsettings('verbose', true, 'solver', ymSolver, 'sdpt3.inftol', 1e-9, 'sdpt3.gaptol', 1e-9, 'sdpt3.maxit', 200);
ymObjective = gam_gkyp + gam_hinf + trace(aa);

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
gam_hinf = sqrt(value(gam_hinf));

A = [zeros(n-1, 1) eye(n-1); -a'];
C = (b - a)';
D = 1;
sys_ss = ss(A, B, C, D);
sys_tf = tf([1 fliplr(b')], [1 fliplr(a')]);