n = 2;
w = [0 pi/32];
gam_gkyp = -1;
gam_hinf = db2mag(5);

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

ntf = tf(zpk([0.95 0.95], [0.5+1i*0.1 0.5-1i*0.1], 1, 1)); % Test TF.

b = fliplr(ntf.num{1}(2:end))'; % Transfer function numerator coefficients.
a = fliplr(ntf.den{1}(2:end))';% Transfer function denominator coefficients.

P_1 = sdpvar(n, n, 'symmetric');
P_2 = sdpvar(n, n, 'symmetric');
Q_script = sdpvar(n, n, 'symmetric');

A_c = [zeros(n-1, 1) eye(n-1); zeros(1, n)];
B = [zeros(n-1, 1); 1];

outerFactor = [A_c B; eye(n) zeros(n, 1)];
Theta = @(Q, P) kron(Phi_d, P) + kron(Psi_d, Q);

ymLmi{1} = [A_c'*P_1*A_c - P_1 + exp(1i*wc)*A_c'*Q_script + exp(-1i*wc)*Q_script*A_c - 2*cos(w0)*Q_script - a*a'...
    A_c'*P_1*B + exp(-1i*wc)*Q_script*B - a...
    b;...
    B'*P_1*A_c + exp(1i*wc)*B'*Q_script - a'...
    B'*P_1*B - 1 ...
    1;...
    b'...
    1 ...
    -gam_gkyp];

ymLmi{2} = [A_c'*P_2*A_c - P_2 - a*a'...
    A_c'*P_2*B - a...
    b;...
    B'*P_2*A_c - a'...
    B'*P_2*B - 1 ...
    1;...
    b'...
    1 ...
    -gam_hinf];

ymConstraint = [ymLmi{1}<=0, ymLmi{2}<=0, Q_script>=0, P_2>=0];

ymSolver = 'sdpt3'; %'SeDuMi';
ymOptions = sdpsettings('verbose', true, 'solver', ymSolver, 'sdpt3.inftol', 1e-12, 'sdpt3.gaptol', 1e-12, 'sdpt3.maxit', 200);
ymObjective = gam_gkyp;

ymDiagnostics = optimize(ymConstraint, ymObjective, ymOptions);
if ymDiagnostics.problem==1
    error('dtsyn:solver_infeasible', ymDiagnostics.info);
elseif ymDiagnostics.problem~=0
    warning('on');
    warning(ymDiagnostics.info);
end

gam_gkyp = sqrt(value(gam_gkyp));
mag2db(gam_gkyp)
gam_hinf = sqrt(value(gam_hinf));

A = [zeros(n-1, 1) eye(n-1); -a'];
C = (b - a)';
D = 1;
sys_ss = ss(A, B, C, D, 1);
sys_tf = tf([1 fliplr(b')], [1 fliplr(a')], 1);