function [sys_tf, gam_gkyp, gam_hinf] = hinfgkypsyn_filter(order, w, gam_gkyp, gam_hinf)

n = order;

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

%b = sdpvar(n, 1, 'full'); % Transfer function numerator coefficients.
%a = sdpvar(n, 1, 'full'); % Transfer function denominator coefficients.
b = [-0.9039207968;4.60196838300844;-9.38171911425817;9.57318276965119;-4.88951073001943];
a = [-0.442070331395993;2.56974468437732;-6.00870809167249;7.06873226323112;-4.18702296473404];

% LMI variables for H-Inf (GKYP) minimization.
P{1} = sdpvar(n, n, 'symmetric');
Q_script = sdpvar(n, n, 'symmetric');
aa = sdpvar(n, n, 'symmetric'); % Convex relaxation, aa is constrained to force it to equal a*a'.

% LMI variables for H-Inf (non-GKYP) minimization.
P{2} = sdpvar(n, n, 'symmetric');

A_c = [zeros(n-1, 1) eye(n-1); zeros(1, n)];
B = [zeros(n-1, 1); 1];
D = 1;

outerFactor = [A_c B; eye(n) zeros(n, 1)];
Theta = @(Q, P, a) kron((Phi_d+[0 0; 0 a]), P) + kron(Psi_d, Q);

ymLmi{1} = [outerFactor'*Theta(Q_script, P{1}, 0)*outerFactor - [a*a' a; a' 1]...
    [b; D];...
    b'...
    D' ...
    -gam_gkyp];
ymLmi{2} = [outerFactor'*Theta(zeros(n), P{2}, 0)*outerFactor - [a*a' a; a' 1]...
    [b; D];...
    b'...
    D'...
    -gam_hinf];
%ymLmi{5} = [aa a; a' 1]; % Shor relaxation to try to force aa==(a*a').
ymConstraint = [ymLmi{1}<=0, Q_script>=0]; % H-Inf (GKYP) constraints.
ymConstraint = [ymConstraint, ymLmi{2}<=0, P{2}>=0, gam_gkyp>=0, gam_hinf>=0];
%ymConstraint = [ymConstraint, ymLmi{5}>=0]; % Convex relaxation constraint.

ymSolver = 'lmilab';%'SeDuMi';
ymOptions = sdpsettings('verbose', true, 'solver', ymSolver, 'lmilab.reltol', 1e-100, 'lmilab.maxiter', 5000, 'lmilab.feasradius', -1, 'lmilab.L', 1000);
%ymOptions = sdpsettings('verbose', true, 'solver', ymSolver, 'sdpt3.inftol', 1e-9, 'sdpt3.gaptol', 1e-9, 'sdpt3.maxit', 1000);
ymObjective = gam_hinf + gam_gkyp% + trace(aa);

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
gam_hinf = sqrt(value(gam_hinf));
gam_hinf

A = [zeros(n-1, 1) eye(n-1); -a'];
C = (b - a)';
D = 1;
sys_ss = ss(A, B, C, D, 1);
sys_tf = tf([1 fliplr(b')], [1 fliplr(a')], 1);

end