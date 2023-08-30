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

k_q = 1;%sdpvar(1);
b = sdpvar(n, 1, 'full'); % Transfer function numerator coefficients.
bdot = sdpvar(n, 1, 'full');
a = sdpvar(n, 1, 'full'); % Transfer function denominator coefficients.
adot = sdpvar(n, 1, 'full');
P{1} = sdpvar(n, n, 'symmetric');

P{2} = sdpvar(n, n, 'symmetric');

Q_script = sdpvar(n, n, 'symmetric');

A_c = [zeros(n-1, 1) eye(n-1); zeros(1, n)];
B = [zeros(n-1, 1); 1];
D = 1;

outerFactor = [A_c B; eye(n) zeros(n, 1)];
Theta = @(Q, P) kron(Phi_d, P) + kron(Psi_d, Q);

ymLmi{1} = -[outerFactor'*Theta(Q_script, P{1})*outerFactor - [zeros(n) (a + k_q*b + adot + k_q*bdot); (a + k_q*b + adot + k_q*bdot)' 1]...
    [(a + adot); D];...
    (a + adot)'...
    D' ...
    (gam_gkyp)] + [(a + k_q*b)*(adot + k_q*bdot)' + (adot + k_q*bdot)*(a + k_q*b)' + (adot + k_q*bdot)*(adot + k_q*bdot)' zeros(n, 2); zeros(2, n) zeros(2, 2)];
ymLmi{2} = -[outerFactor'*Theta(zeros(n), P{2})*outerFactor - [zeros(n) (a + k_q*b + adot + k_q*bdot); (a + k_q*b + adot + k_q*bdot)' 1]...
    [(a + adot); D];...
    (a + adot)'...
    D'...
    -(gam_hinf)] + [(a + k_q*b)*(adot + k_q*bdot)' + (adot + k_q*bdot)*(a + k_q*b)' + (adot + k_q*bdot)*(adot + k_q*bdot)' zeros(n, 2); zeros(2, n) zeros(2, 2)];
ymConstraint = [ymLmi{1}>=0, ymLmi{2}>=0, (Q_script)>=0, (P{2})>=0];
% ymConstraint = [ymLmi{1}>=0, (Q_script + Q_scriptdot)>=0];
%ymConstraint = [ymConstraint, 0.999 <= k_q <= 1.001, uncertain(k_q)];

ymSolver = 'lmilab';
ymOptions = sdpsettings('verbose', true, 'solver', ymSolver, 'lmilab.reltol', 1e-6, 'lmilab.maxiter', 500, 'lmilab.feasradius', 1e9, 'lmilab.L', 100);
%ymOptions = sdpsettings('verbose', true, 'solver', ymSolver);%, 'sdpt3.inftol', 1e-9, 'sdpt3.gaptol', 1e-9, 'sdpt3.maxit', 200);

epsilon = 1e-6;
maxIter = 500; %4995;
%ntf = tf(synthesizeNTF(4, 32, 3));
%a_val = fliplr(ntf.den{1}(2:end))';
%a_val = [0 0 0 0]';
%a_val = [0.2401;-1.36374876879259;2.91645060693392;-2.78316075263793];
%a_val = [0.444444297096043;-2.13578877231391;3.89195992817283;-3.19447366909881]; % synthesizeNTF(4, 32, 0)
% Output from gkypsyn_filter2.m:
% a_val = [0.0527044035466642;0.0145048559447814;-0.00237602832108952;-0.00559032375529501];
% b_val = [0.0977261790556061;-0.024451995254671;-0.305384923269644;-0.745327501431005];
Q_val = [27402915.6342544 -81719998.6098839 81431432.9911447 -27112960.4308305;-81719998.6098839 245611375.709121 -246648332.99835 82761959.1169966;81431432.9911447 -246648332.99835 249579092.012427 -84376291.4254008;-27112960.4308305 82761959.1169966 -84376291.4254008 28735065.2485522];
P_val_1 = [-25773346.2673357 76982566.3002704 -76828146.2177914 25618207.8413842;76982566.3002704 -231825628.080862 233248582.468432 -78412367.4872749;-76828146.2177914 233248582.468432 -236559402.736401 80154880.7489142;25618207.8413842 -78412367.4872749 80154880.7489142 -27369136.3234336];
P_val_2 = [0.0219007792548243 0.0348835439119598 0.0267488144597126 0.000805202607148859;0.0348835439119598 0.0958193794048255 0.120520873180772 0.104629630248212;0.0267488144597126 0.120520873180772 0.229391311555292 0.30012351112534;0.000805202607148859 0.104629630248212 0.30012351112534 0.543437458103619];
a_val = [0.0977261790556061;-0.024451995254671;-0.305384923269644;-0.745327501431005];
b_val = ([0.0527044035466642;0.0145048559447814;-0.00237602832108952;-0.00559032375529501] - [0.0977261790556061;-0.024451995254671;-0.305384923269644;-0.745327501431005])/k_q;

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
%     ymConstraintA = replace(ymConstraintA, gam_gkypdot, gam_gkyp_val(k));
%     ymConstraintA = replace(ymConstraintA, gam_hinf, gam_hinf_val(k));
    
    ymDiagnostics = optimize(ymConstraintA, ymObjective, ymOptions);
    k = k + 1;
    
    a_val(:, k) = a_val(:, k-1) + value(a);
    b_val(:, k) = b_val(:, k-1) + value(b);
    %Q_val(:, :, k) = Q_val(:, :, k-1) + value(Q_script);
    %P_val_1(:, :, k) = P_val_1(:, :, k-1) + value(P{1});
    %P_val_2(:, :, k) = P_val_2(:, :, k-1) + value(P{2});
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

% A = [zeros(n-1, 1) eye(n-1); -(a' + ];
% C = (b - a)';
% D = 1;
% sys_ss = ss(A, B, C, D, 1);
sys_tf = tf([1 fliplr(a')], [1 fliplr(a') + 1*fliplr(b')], 1);