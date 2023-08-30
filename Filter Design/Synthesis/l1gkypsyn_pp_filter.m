function l1gkypsyn_pp_filter
n = 4;
osr = 32;
w = [0 pi/osr];
gam_gkyp = -1;
gam_l1 = 8;

% Initial condition.
ntf = tf(synthesizeNTF(n, osr)); % Test TF.
a_init = fliplr(ntf.den{1}(2:end))';% Transfer function denominator coefficients.

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

% LMI variables for regional pole placement.
outerAngle = pi/4;
P{3} = sdpvar(n, n, 'symmetric');
Y = sdpvar(1, n, 'full');
L = [0 0; 0 0];
M = [-sin(outerAngle) cos(outerAngle); -cos(outerAngle) -sin(outerAngle)];

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
ymLmi{3} = outerFactor'*Theta(zeros(n), P{2}, alpha)*outerFactor - [aa a; a' 1];
ymLmi{4} = [gam_l1 mu nu; mu 1 0; nu 0 1];
ymLmi{5} = [aa a; a' 1]; % Shor relaxation to try to force aa==(a*a').
ymLmi{6} = kron(L, P{3}) + kron(M, P{3}*A_c - Y'*B') + kron(M', A_c*P{3} - B*Y);
ymConstraint = [ymLmi{1}<=0, Q_script>=0]; % H-Inf (GKYP) constraints.
ymConstraint = [ymConstraint, ymLmi{2}>=0, ymLmi{3}<=0, ymLmi{4}>=0, P{2}>=0, nu>=0, mu>=0]; % l1/(*) constraints.
ymConstraintCvx = [ymConstraint, ymLmi{5}>=0]; % Convex relaxation constraint.

E = [zeros(n-1, 1) eye(n-1); -a_init'];
ymLmi{7} = kron(L, P{3}) + kron(M, P{3}*E') + kron(M', E*P{3});
ymConstraint1 = replace(ymConstraint, a, a_init);
ymConstraint1 = replace(ymConstraint1, aa, a_init*a_init');
ymConstraint1 = [ymConstraint1, ymLmi{7}<=0, P{3}>=0];

evaluateLmis(ymConstraint1, false);
P3val = value(P{3})
gg1 = value(gam_gkyp)

%ymConstraint2 = replace(ymConstraintCvx, P{3}, P3val);
ymLmi{8} = kron(L, P3val) + kron(M, P3val*A_c' - Y'*B') + kron(M', A_c*P3val - B*Y);
ymConstraint2 = [ymConstraintCvx, ymLmi{8}<=0];
ymConstraint2 = replace(ymConstraint2, a, (Y/P3val)');
evaluateLmis(ymConstraint2, true);
YVal = value(Y)
gg2 = value(gam_gkyp)

a_init = (YVal/P3val)';
E = [zeros(n-1, 1) eye(n-1); -a_init'];
ymLmi{7} = kron(L, P{3}) + kron(M, P{3}*E') + kron(M', E*P{3});
ymConstraint3 = replace(ymConstraint, a, a_init);
ymConstraint3 = replace(ymConstraint3, aa, a_init*a_init');
ymConstraint3 = [ymConstraint3, ymLmi{7}<=0, P{3}>=0];

evaluateLmis(ymConstraint3, false);
P3val = value(P{3})
gg3 = value(gam_gkyp)

ymLmi{8} = kron(L, P3val) + kron(M, P3val*A_c' - Y'*B') + kron(M', A_c*P3val - B*Y);
ymConstraint4 = [ymConstraintCvx, ymLmi{8}<=0];
ymConstraint4 = replace(ymConstraint4, a, (Y/P3val)');
evaluateLmis(ymConstraint4, true);
YVal = value(Y)
gg4 = value(gam_gkyp)

a_init = (YVal/P3val)';
E = [zeros(n-1, 1) eye(n-1); -a_init'];
ymLmi{7} = kron(L, P{3}) + kron(M, P{3}*E') + kron(M', E*P{3});
ymConstraint5 = replace(ymConstraint, a, a_init);
ymConstraint5 = replace(ymConstraint5, aa, a_init*a_init');
ymConstraint5 = [ymConstraint5, ymLmi{7}<=0, P{3}>=0];

evaluateLmis(ymConstraint5, false);
P3val = value(P{3})
gg5 = value(gam_gkyp)

ymLmi{8} = kron(L, P3val) + kron(M, P3val*A_c' - Y'*B') + kron(M', A_c*P3val - B*Y);
ymConstraint6 = [ymConstraintCvx, ymLmi{8}<=0];
ymConstraint6 = replace(ymConstraint6, a, (Y/P3val)');
evaluateLmis(ymConstraint6, true);
YVal = value(Y)
gg6 = value(gam_gkyp)

a_init = (YVal/P3val)';
E = [zeros(n-1, 1) eye(n-1); -a_init'];
ymLmi{7} = kron(L, P{3}) + kron(M, P{3}*E') + kron(M', E*P{3});
ymConstraint7 = replace(ymConstraint, a, a_init);
ymConstraint7 = replace(ymConstraint7, aa, a_init*a_init');
ymConstraint7 = [ymConstraint7, ymLmi{7}<=0, P{3}>=0];

evaluateLmis(ymConstraint7, false);
P3val = value(P{3})
gg7 = value(gam_gkyp)

ymLmi{8} = kron(L, P3val) + kron(M, P3val*A_c' - Y'*B') + kron(M', A_c*P3val - B*Y);
ymConstraint8 = [ymConstraintCvx, ymLmi{8}<=0];
ymConstraint8 = replace(ymConstraint8, a, (Y/P3val)');
evaluateLmis(ymConstraint8, true);
YVal = value(Y)
gg8 = value(gam_gkyp)

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

function evaluateLmis(ymConstraint, cvxRelaxation)
    ymSolver = 'sdpt3'; %'SeDuMi';
    ymOptions = sdpsettings('verbose', true, 'solver', ymSolver, 'sdpt3.inftol', 1e-9, 'sdpt3.gaptol', 1e-9, 'sdpt3.maxit', 200);
    ymObjective = gam_l1 + gam_gkyp;
    if cvxRelaxation
        ymObjective = ymObjective + trace(aa);
    end

    ymLineSearch = optimizer(ymConstraint, ymObjective, ymOptions, alpha, ymObjective);
    lineSearchOptions = optimset('Display', 'iter', 'FunValCheck', 'on');
    [alphaOpt, ymLineSearchOpt] = fminbnd(@(x) ymLineSearch(x), 0, 1, lineSearchOptions);
    ymConstraintLin = replace(ymConstraint, alpha, alphaOpt);

    ymDiagnostics = optimize(ymConstraintLin, ymObjective, ymOptions);
    if ymDiagnostics.problem==1
        error('dtsyn:solver_infeasible', ymDiagnostics.info);
    elseif ymDiagnostics.problem~=0
        warning('on');
        warning(ymDiagnostics.info);
    end
end

function feasLmis(ymConstraint, cvxRelaxation)
    ymSolver = 'sdpt3'; %'SeDuMi';
    ymOptions = sdpsettings('verbose', true, 'solver', ymSolver, 'sdpt3.inftol', 1e-9, 'sdpt3.gaptol', 1e-9, 'sdpt3.maxit', 200);
    ymObjective = [];

%     ymLineSearch = optimizer(ymConstraint, ymObjective, ymOptions, alpha, ymObjective);
%     lineSearchOptions = optimset('Display', 'iter', 'FunValCheck', 'on');
%     [alphaOpt, ymLineSearchOpt] = fminbnd(@(x) ymLineSearch(x), 0, 1, lineSearchOptions);
    alphaOpt = (1-(max(abs(eig(E))))^2)/2;
    ymConstraintLin = replace(ymConstraint, alpha, alphaOpt);

    ymDiagnostics = optimize(ymConstraintLin, ymObjective, ymOptions);
    if ymDiagnostics.problem==1
        error('dtsyn:solver_infeasible', ymDiagnostics.info);
    elseif ymDiagnostics.problem~=0
        warning('on');
        warning(ymDiagnostics.info);
    end
end
end