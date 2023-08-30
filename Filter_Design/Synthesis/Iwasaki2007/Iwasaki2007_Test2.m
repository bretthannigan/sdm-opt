% Equation 
% Check if matrix sizes are compatible.

ny = 1;
nu = 1;
nz = 2;
nw = 2;
np = 3;
nc = 3;
n = np + nc;
P_script = sdpvar(n, n, 'symmetric', 'complex');%sym('P_script', [n n]);
Q_script = sdpvar(n, n, 'symmetric', 'complex');%sym('Q_script', [n n]);
X = sdpvar(np, np, 'full');%sym('X', [np np]);
L = sdpvar(nu, np, 'full');%sym('L', [nu np]);
Y = sdpvar(np, np, 'full');%sym('Y', [np np]);
F = sdpvar(np, ny, 'full');%sym('F', [np ny]);
Q = sdpvar(np, np, 'full');%sym('Q', [np np]);
R = sdpvar(nu, ny, 'full');%sym('R', [nu ny]);
S = sdpvar(np, np, 'full');%sym('S', [np np]);
sys = drss(np, nz+ny, nw+nu);
A = sys.A;
B1 = sys.B(:, 1:nw);
B2 = sys.B(:, nw+1:end);
C1 = sys.C(1:nz, :);
C2 = sys.C(nz+1:end, :);
D11 = sys.D(1:nz, 1:nw);
D12 = sys.D(1:nz, nw+1:end);
D21 = sys.D(nz+1:end, 1:nw);
D22 = sys.D(nz+1:end, nw+1:end);

gam = sdpvar(1);%sym('gam');
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