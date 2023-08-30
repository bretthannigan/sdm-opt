isExtended = false;
ny = 1;
nu = 1;
nz = 1;
nw = 1;
np = 2;
nc = 2;
n = np + nc;
P_script = sdpvar(n, n, 'symmetric', 'complex');
Q_script = sdpvar(n, n, 'symmetric', 'complex');
X = sdpvar(np, np, 'full');
L = sdpvar(nu, np, 'full');
Y = sdpvar(np, np, 'full');
F = sdpvar(np, ny, 'full');
Q = sdpvar(np, np, 'full');
R = sdpvar(nu, ny, 'full');
S = sdpvar(np, np, 'full');

A_script = [A*X + B{2}*L A + B{2}*R*C{2};...
    Q Y*A + F*C{2};...
    -X -eye(np);...
    -S -Y;...
    C{1}*X + D{1,2}*L C{1} + D{1,2}*R*C{2}];

B_script = [B{1} + B{2}*R*D{2,1};...
    Y*B{1} + F*D{2,1};...
    zeros(2*np, nw);...
    D{1,1} + D{1,2}*R*D{2,1}];
multiplier = [eye(2*n) zeros(2*n, nz); zeros(2*n, 2*n + nz); zeros(nz, 2*n) eye(nz)];
multiplier = [multiplier B_script];

%% Frequency properties.
if isPositiveReal
    Pi = [zeros(nz) -eye(nz, nw); -eye(nw, nz) zeros(nw)];
else % Bounded real (small gain theorem) - H infinity norm.
    gam = sdpvar(1);
    Pi = [-gam*eye(nz) zeros(nz, nw); zeros(nw, nz) eye(nw)];
end

%% Finite frequency interval.
w1 = 0;
w2 = pi/4;
wc = (w2 + w1)/2;
w0 = (w2 - w1)/2;
Phi_d = [-1 0; 0 1];
Psi_d = [-2*cos(w0) exp(1i*wc); exp(-1i*wc) 0];

ymLmi{1} = multiplier*...
    [kron(Phi_d, P_script) + kron(Psi_d, Q_script) zeros(nz + nw);...
    zeros(nz + nw) Pi]...
    *multiplier;
