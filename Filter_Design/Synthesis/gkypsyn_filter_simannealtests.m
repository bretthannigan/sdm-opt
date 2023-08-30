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

%b = [0.9767;-3.9205;5.9106;-3.9668];
b = sdpvar(n, 1);
a0 = [0.4442;-2.1086;3.8462;-3.17];
a = sdpvar(n, 1);
%b = [-0.7737809375;4.07253125;-8.57375;9.025;-4.75];
test = sdpvar(1);
P{1} = sdpvar(n, n, 'symmetric');
P{2} = sdpvar(n, n, 'symmetric');
Q_script = sdpvar(n, n, 'symmetric');
aa = sdpvar(n, n, 'symmetric');

A_c = [zeros(n-1, 1) eye(n-1); zeros(1, n)];
B = [zeros(n-1, 1); 1];

outerFactor = [A_c B; eye(n) zeros(n, 1)];
Theta = @(Q, P) kron(Phi_d, P) + kron(Psi_d, Q);

ymLmi{1} = [outerFactor'*Theta(Q_script, P{1})*outerFactor - [a*a' a; a' 1]...
    [b; 1];...
    b'...
    1 ...
    -gam_gkyp];
ymLmi{2} = [outerFactor'*Theta(zeros(n), P{2})*outerFactor - [a*a' a; a' 1]...
    [b; 1];...
    b'...
    1 ...
    -gam_hinf];
ymLmi{3} = [diag(a) a; a' 1]; % Shor relaxation to try to force aa==(a*a').
ymConstraint = [ymLmi{1}<=0, ymLmi{2}<=0, Q_script>=0, P{2}>=0];

ymSolver = 'sdpt3'; %'SeDuMi';
ymOptions = sdpsettings('verbose', true, 'solver', ymSolver);%, 'sdpt3.inftol', 1e-12, 'sdpt3.gaptol', 1e-12, 'sdpt3.steptol', 1e-12, 'sdpt3.maxit', 2000);
%ymOptions = sdpsettings('verbose', true, 'solver', ymSolver, 'lmirank.solver', 'sedumi', 'sedumi.eps', 0, 'lmirank.maxiter', 10000);
ymObjective = gam_gkyp;

ymOptimizer = optimizer(ymConstraint, ymObjective, ymOptions, a, ymObjective);

anew = simulannealbnd(@(x) siman(ymOptimizer, x), a0);
        
ymConstraint = replace(ymConstraint, a, anew)
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
sys_ss = ss(A, B, C, D, 1);
sys_tf = tf([1 fliplr(b')], [1 fliplr(a')], 1);

function out = siman(funhandle, x)
    [a, b] = funhandle(x);
    if b~=0
        out = 10000;
    else
        out = a;
    end
end