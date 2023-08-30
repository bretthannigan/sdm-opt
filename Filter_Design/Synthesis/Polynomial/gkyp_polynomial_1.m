% The following based on Corollary 1 in Section 4.2 and "Digital filter 
% design" in Section 7.1 of [Iwasaki2003a].

% This version enforces the spectral factorization scalar term to be 1 so
% that the optimization is done without any shift in magnitude.

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

% B = [sdpvar(ones(1, n), ones(1, n)) 1];
C = [sdpvar(ones(1, n), ones(1, n)) 0];
A = [sdpvar(ones(1, n+1), ones(1, n+1))];

gam = db2mag(-50);
gam_pb = 1.5;

% P = sdpvar(n, n, 'hermitian', 'complex');
% P = sdpvar(n, n, 'hermitian', 'complex');
P = sdpvar(n, n, 'symmetric');
Q = sdpvar(n, n, 'symmetric');
P2 = sdpvar(n, n, 'symmetric');
Q2 = sdpvar(n, n, 'symmetric');
G = sdpvar(1, n+1);
G2 = sdpvar(1, n+1);

X_B = sdpvar(n+1, n+1, 'symmetric');
Q_B = sdpvar(n, n, 'symmetric');
X_B(n+1, n+1) = 1;

X_A = sdpvar(n+1, n+1, 'symmetric');
Q_A = sdpvar(n, n, 'symmetric');

J = [1 zeros(1, n)];
U = [zeros(n, 1) eye(n)];
V = [eye(n) zeros(n, 1)];
F = [U; V];

Phi = [1 0; 0 -1];
Phi2 = [-1 0; 0 100^2];

lmi{1} = X_B + F'*kron(Phi, Q_B)*F;
lmi{2} = X_A + F'*kron(Phi, Q_A)*F;

ymConstraint = [];
ymConstraint = [ymConstraint, C{1} == sum(diag(lmi{1}, 0))]; %#ok<AGROW>
for iDiag=1:n
    ymConstraint = [ymConstraint, C{iDiag + 1} == sum(diag(lmi{1}, iDiag))]; %#ok<AGROW>
end
ymConstraint = [ymConstraint, X_B>=0, Q_B>=0];

ymConstraint = [ymConstraint, A{1} == sum(diag(lmi{2}, 0))]; %#ok<AGROW>
for iDiag=1:n
    ymConstraint = [ymConstraint, A{iDiag + 1} == sum(diag(lmi{2}, iDiag))]; %#ok<AGROW>
end
ymConstraint = [ymConstraint, X_A>=0, Q_A>=0];

G(1) = (C{1} - (gam - 1)^2*A{1})/2;
for iG=2:(n+1)
    G(iG) = (C{iG} - (gam - 1)^2*A{iG}); %#ok<SAGROW>
end
lmi{3} = F'*(kron(-Phi, P) + kron(Psi, Q))*F - (G'*J + J'*G);
ymConstraint = [ymConstraint, lmi{3}<=0, Q>=0];

G2(1) = -(C{1})/2;
for iG=2:(n+1)
    G2(iG) = -(C{iG}); %#ok<SAGROW>
end
lmi{4} = F'*(kron(-Phi, P2) + kron(Psi, Q2))*F - (G2'*J + J'*G2);
ymConstraint = [ymConstraint, lmi{4}>=0, Q2>=0];

% G_pb(1) = -(B{1} - gam_pb^2*A{1})/2;
% for iG_pb=2:(n+1)
%     G_pb(iG_pb) = -(B{iG_pb} - gam_pb^2*A{iG_pb}); %#ok<SAGROW>
% end
% lmi{4} = F'*(kron(-Phi, P2) + kron(Phi, Q2))*F - (G_pb'*J + J'*G_pb);
% ymConstraint = [ymConstraint, lmi{4}<=0, Q2>=0];

% G2(1) = (B{1} - A{1})/2;
% for iG=2:(n+1)
%     G2(iG) = (B{iG} - A{iG}); %#ok<SAGROW>
% end
% lmi{4} = F'*(kron(Phi2, P2) + kron(Phi2, Q2))*F - (G2'*J + J'*G2);
% ymConstraint = [ymConstraint, lmi{4}<=0, Q2>=0];

ymSolver = 'mosek'; %'SeDuMi';
ymOptions = sdpsettings('verbose', true, 'solver', ymSolver, 'sdpt3.steptol', 1e-12, 'sdpt3.inftol', 1e-12, 'sdpt3.gaptol', 1e-12, 'sdpt3.maxit', 200);
ymObjective = [];

ymDiagnostics = optimize(ymConstraint, ymObjective, ymOptions);
if ymDiagnostics.problem==1
    error('dtsyn:solver_infeasible', ymDiagnostics.info);
elseif ymDiagnostics.problem~=0
    warning('on');
    warning(ymDiagnostics.info);
end

C = fliplr(cellfun(@value, C));
A = fliplr(cellfun(@value, A));
sys_tf = tf(C, A, 1);
C_sys = [C(2:end-1) C(end) C(end-1:-1:2)];
A_sys = [A(1:end-1) A(end) A(end-1:-1:1)];
H_sys = tf(C_sys, A_sys, 1);
% H_sys = tf([B B(end-1:-1:1)], [A A(end-1:-1:1)], 1);
%H_sys = sys_tf + sys_tf';
%[G_sys, S_sys] = spectralfact(H_sys);
G_sys = tf([0 -sfact(C_sys)]+sfact(A_sys), sfact(A_sys), 1);
bode(G_sys)