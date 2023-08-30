% The following based on Corollary 1 in Section 4.2 and "Digital filter 
% design" in Section 7.1 of [Iwasaki2003a].

n = 4;
w_s = [0 pi/32];

if w_s(1)==0 % Low frequency case, Psi_d = [-2*cos(w(1)) 1; 1 0]
    w_s = [-w_s(2) w_s(2)];
elseif w_s(2)==pi % High  frequency case, Psi_d = [2*cos(w(2)) -1; -1 0]
    w_s = [w_s(1) 2*pi-w_s(1)];
end

wc_s = (w_s(2) + w_s(1))/2;
w0_s = (w_s(2) - w_s(1))/2;
Psi = [0 exp(1i*wc_s); exp(-1i*wc_s) -2*cos(w0_s)];

B = sdpvar(ones(1, n+1), ones(1, n+1));
A = [sdpvar(ones(1, n), ones(1, n)) 1];

gam = 0.001;

P = sdpvar(n, n, 'hermitian', 'complex');
Q = sdpvar(n, n, 'hermitian', 'complex');
G = sdpvar(1, n+1);

X_B = sdpvar(n+1, n+1, 'symmetric');
Q_B = sdpvar(n, n, 'symmetric');

X_A = sdpvar(n+1, n+1, 'symmetric');
Q_A = sdpvar(n, n, 'symmetric');

X_G = sdpvar(n+1, n+1, 'symmetric');
Q_G = sdpvar(n, n, 'symmetric');

J = [1 zeros(1, n)];
U = [zeros(n, 1) eye(n)];
V = [eye(n) zeros(n, 1)];
F = [U; V];

Phi = [1 0; 0 -1];
Phi_test = [-1 0; 0 -1];
Psi_B = [0 1; 1 -2];
Psi_A = Psi_B;

lmi{1} = X_B + F'*kron(Phi, Q_B)*F;
lmi{2} = X_A + F'*kron(Phi, Q_A)*F;
lmi{3} = X_G + F'*kron(Phi, Q_G)*F;

ymConstraint = [];
for iDiag=0:n
    ymConstraint = [ymConstraint, B{iDiag + 1} == sum(diag(lmi{1}, iDiag))]; %#ok<AGROW>
end
ymConstraint = [ymConstraint, X_B>=0, Q_B>=0];

for iDiag=0:n
    ymConstraint = [ymConstraint, A{iDiag + 1} == sum(diag(lmi{2}, iDiag))]; %#ok<AGROW>
end
ymConstraint = [ymConstraint, X_A>=0, Q_A>=0];

ymConstraint = [ymConstraint, (B{1} - gam^2*A{1})/2 == trace(lmi{3})];
for iDiag=1:n
    ymConstraint = [ymConstraint, (B{iDiag + 1} - gam^2*A{iDiag + 1}) == sum(diag(lmi{3}, iDiag))]; %#ok<AGROW>
end
ymConstraint = [ymConstraint, X_G<=0, Q_G<=0];

% G(1) = (B{1} - gam^2*A{1})/2;
% for iG=2:(n+1)
%     G(iG) = B{iG} - gam^2*A{iG}; %#ok<SAGROW>
% end
% lmi{3} = F'*(kron(Phi, P) + kron(Phi_test, Q))*F - (G'*J + J'*G);
% ymConstraint = [ymConstraint, (lmi{3} + lmi{3}')<=0, (Q + Q')>=0];

ymSolver = 'sedumi'; %'SeDuMi';
ymOptions = sdpsettings('verbose', true, 'solver', ymSolver, 'sdpt3.inftol', 1e-9, 'sdpt3.gaptol', 1e-9, 'sdpt3.maxit', 200);
ymObjective = [];

ymDiagnostics = optimize(ymConstraint, ymObjective, ymOptions);
if ymDiagnostics.problem==1
    error('dtsyn:solver_infeasible', ymDiagnostics.info);
elseif ymDiagnostics.problem~=0
    warning('on');
    warning(ymDiagnostics.info);
end

B = fliplr(cellfun(@value, B));
A = fliplr(cellfun(@value, A));
sys_tf = tf(B, A, 1);

B_sys = [B(1:end-1) 2*B(end) B(end-1:-1:1)];
A_sys = [A(1:end-1) 2*A(end) A(end-1:-1:1)];
H_sys = tf(B_sys, A_sys, 1);
%H_sys = tf([B B(end-1:-1:1)], [A A(end-1:-1:1)], 1);
%H_sys = sys_tf + sys_tf';
[G_sys, S_sys] = spectralfact(H_sys);