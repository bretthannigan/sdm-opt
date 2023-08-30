n = 4;
w = [0 pi/32];
gam_gkyp = -1;
gam_hinf = 1.5;

if gam_gkyp==-1
    gam_gkyp = sdpvar(1); % Optimization problem.
    gam_gkypdot = sdpvar(1);
else
    gam_gkyp = gam_gkyp.^2; % Sub-optimal (feasibility) problem.
    gam_gkypdot = 0;
end
if gam_hinf==-1
    gam_hinf = sdpvar(1); % Optimization problem.
    gam_hinfdot = sdpvar(1);
else
    gam_hinf = gam_hinf.^2; % Sub-optimal (feasibility) problem.
    gam_hinfdot = 0;
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
bdot = sdpvar(n, 1, 'full');
a = sdpvar(n, 1, 'full'); % Transfer function denominator coefficients.
adot = sdpvar(n, 1, 'full');
P{1} = sdpvar(n, n, 'symmetric');
Pdot{1} = sdpvar(n, n, 'symmetric');
P{2} = sdpvar(n, n, 'symmetric');
Pdot{2} = sdpvar(n, n, 'symmetric');
Q_script = sdpvar(n, n, 'symmetric');
Q_scriptdot = sdpvar(n, n, 'symmetric');

A_c = [zeros(n-1, 1) eye(n-1); zeros(1, n)];
B = [zeros(n-1, 1); 1];
D = 1;

outerFactor = [A_c B; eye(n) zeros(n, 1)];
Theta = @(Q, P) kron(Phi_d, P) + kron(Psi_d, Q);

ymLmi{1} = -[outerFactor'*Theta(Q_script + Q_scriptdot, P{1} + Pdot{1})*outerFactor - [zeros(n) (a + adot); (a + adot)' 1]...
    [(b + bdot); D];...
    (b + bdot)'...
    D' ...
    -(gam_gkyp)] + [a*adot' + adot*a' + adot*adot' zeros(n, 2); zeros(2, n) zeros(2, 2)];
ymLmi{2} = -[outerFactor'*Theta(zeros(n), P{2} + Pdot{2})*outerFactor - [zeros(n) (a + adot); (a + adot)' 1]...
    [(b + bdot); D];...
    (b + bdot)'...
    D'...
    -(gam_hinf)] + [a*adot' + adot*a' + adot*adot' zeros(n, 2); zeros(2, n) zeros(2, 2)];
ymConstraint = [ymLmi{1}>=0, ymLmi{2}>=0, (Q_script + Q_scriptdot)>=0, (P{2} + Pdot{2})>=0];

ymSolver = 'mosek';
%ymOptions = sdpsettings('verbose', true, 'solver', ymSolver, 'lmilab.reltol', 1e-6, 'lmilab.maxiter', 500, 'lmilab.feasradius', 1e9, 'lmilab.L', 100);
ymOptions = sdpsettings('verbose', true, 'solver', ymSolver);%, 'sdpt3.inftol', 1e-9, 'sdpt3.gaptol', 1e-9, 'sdpt3.maxit', 200);

epsilon = 1e-6;
maxIter = 5000;
%ntf = tf(synthesizeNTF(4, 32, 3));
%a_val = fliplr(ntf.den{1}(2:end))';
%a_val = [0 0 0 0]';
%a_val = [0.2401;-1.36374876879259;2.91645060693392;-2.78316075263793];
%a_val = [0.444444297096043;-2.13578877231391;3.89195992817283;-3.19447366909881]; % synthesizeNTF(4, 32, 0)
% Output from gkypsyn_filter2.m:
a_val = [0.0527046286435801;0.0145037390691871;-0.0023769083058262;-0.00559069650500656];
b_val = [0.0977288797365547;-0.0244535905981645;-0.305386650896797;-0.745326196719701];
Q_val = [732.177213758847 -305.08396811617 -298.107189949917 -114.180760674288;-305.08396811617 937.981023307973 -275.554542338686 -339.293958609526;-298.107189949917 -275.554542338686 931.19676003202 -340.37221128484;-114.180760674288 -339.293958609526 -340.37221128484 812.841045707383];
P_val_1 = [-613.654710262479 119.807384474716 321.008291037806 171.872622910439;119.807384474716 -574.950702097776 82.7306317983261 359.05409200347;321.008291037806 82.7306317983261 -623.634292104916 196.23601012712;171.872622910439 359.05409200347 196.23601012712 -758.348125573109];
P_val_2 = [0.0230025284079773 0.0363593196304074 0.0282241031450044 0.00190314113558974;0.0363593196304074 0.0978001531322096 0.122499477364535 0.106103631581994;0.0282241031450044 0.122499477364535 0.231369891581577 0.301595724217668;0.00190314113558974 0.106103631581994 0.301595724217668 0.544533237050318];
gam_gkyp_val = 0.1;
gam_hinf_val = 1.5;

objective_val(1:2, :) = 1e6;
relTol = 1;
k = 1;
kappa = 0.1;

ymObjective = gam_gkyp;

while k<maxIter %any(relTol(k)>epsilon) && k<maxIter
    ymConstraintA = replace(ymConstraint, adot, a_val(:, k));
    ymConstraintA = replace(ymConstraintA, bdot, b_val(:, k));
    ymConstraintA = replace(ymConstraintA, Q_scriptdot, Q_val(:, :, k));
    ymConstraintA = replace(ymConstraintA, Pdot{1}, P_val_1(:, :, k));
    ymConstraintA = replace(ymConstraintA, Pdot{2}, P_val_2(:, :, k));
%     ymConstraintA = replace(ymConstraintA, gam_gkypdot, gam_gkyp_val(k));
%     ymConstraintA = replace(ymConstraintA, gam_hinf, gam_hinf_val(k));
    
    ymDiagnostics = optimize(ymConstraintA, ymObjective, ymOptions);
    k = k + 1;
    
    a_val(:, k) = a_val(:, k-1) + value(a);
    b_val(:, k) = b_val(:, k-1) + value(b);
    Q_val(:, :, k) = Q_val(:, :, k-1) + value(Q_script);
    P_val_1(:, :, k) = P_val_1(:, :, k-1) + value(P{1});
    P_val_2(:, :, k) = P_val_2(:, :, k-1) + value(P{2});
    gam_gkyp_val(k) = value(gam_gkyp);
    %gam_hinf_val(k) = gam_hinf_val(k-1) + value(gam_hinf);
    
    objective_val(:, k) = [value(gam_gkyp); value(gam_hinf)];
    relTol(k) = norm(a_val(:, k)- a_val(:, k-1));
end

% ymSolver = 'lmilab'; %'SeDuMi';
% ymOptions = sdpsettings('verbose', true, 'solver', ymSolver, 'lmilab.reltol', 1e-6, 'lmilab.maxiter', 500, 'lmilab.feasradius', 1e9, 'lmilab.L', 100);
% 
% ymConstraint = replace(ymConstraint, adot, a_val(:, end));
% ymDiagnostics = optimize(ymConstraint, ymObjective, ymOptions);
% if ymDiagnostics.problem==1
%     error('dtsyn:solver_infeasible', ymDiagnostics.info);
% elseif ymDiagnostics.problem~=0
%     warning('on');
%     warning(ymDiagnostics.info);
% end

a = a_val(:, end);%value(a);
b = b_val(:, end);%value(b);
gam_gkyp = sqrt(value(gam_gkyp));
gam_hinf = sqrt(value(gam_hinf));

A = [zeros(n-1, 1) eye(n-1); -a'];
C = (b - a)';
D = 1;
sys_ss = ss(A, B, C, D, 1);
sys_tf = tf([1 fliplr(b')], [1 fliplr(a')], 1);