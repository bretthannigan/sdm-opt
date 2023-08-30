% Equation 
% Check if matrix sizes are compatible.

ny = 1;
nu = 2;
nz = 3;
nw = 5;
np = 7;
nc = 11;
n = np + nc;
P_script = sym('P_script', [n n]);
Q_script = sym('Q_script', [n n]);
X = sym('X', [np np]);
L = sym('L', [nu np]);
Y = sym('Y', [np np]);
F = sym('F', [np ny]);
Q = sym('Q', [np np]);
R = sym('R', [nu ny]);
S = sym('S', [np np]);
A = sym('A', [np np]);
B1 = sym('B1', [np nw]);
B2 = sym('B2', [np nu]);
C1 = sym('C1', [nz np]);
C2 = sym('C2', [ny np]);
D11 = sym('D11', [nz nw]);
D12 = sym('D12', [nz nu]);
D21 = sym('D21', [ny nw]);
D22 = sym('D22', [ny nu]);
gam = sym('gam');
Pi = [-gam*eye(nz) zeros(nz, nw); zeros(nw, nz) eye(nw)];
w1 = 0;
w2 = pi/4;
wc = (w2 + w1)/2;
w0 = (w2 - w1)/2;
Phi_d = [-1 0; 0 1];
Psi_d = [-2*cos(w0) exp(1i*wc); exp(-1i*wc) 0];
R_script = [eye(2*np) exp(1i*wc)*eye(2*np) zeros(2*np, nz)];
A_script = [A*X + B2*L A + B2*R*C2;...
    Q Y*A + F*C2;...
    -X -eye(np);...
    -S -Y;...
    C1*X + D12*L C1 + D12*R*C2];
B_script = [B1 + B2*R*D21;...
    Y*B1 + F*D21;...
    zeros(2*np, nw);...
    D11 + D12*R*D21];
multiplier = blkdiag(eye(4*np, 2*n), eye(nz));
multiplier = [multiplier B_script];
lhs = multiplier*...
    [kron(Phi_d, P_script) + kron(Psi_d, Q_script) zeros(2*n, nz + nw);...
    zeros(nz + nw, 2*n) Pi]...
    *multiplier';
rhs = A_script*R_script + (A_script*R_script)';
lhs + rhs;