%SLIDING_MODE Simple Sliding Mode Control
%   Example of sliding mode control performed to control the position and
%   velocity of a unit mass subject to a disturbance force in continuous 
%   time.
%
%   See also:
%       Shtessel Y., Edwards C., Fridman L., Levant A., Sliding Mode
%       Control and Observation. 2014. pp.1-9.

%% Simulation.
t_s = 0.0001; % s.
t_f = 8; % s.
t = 0:t_s:t_f;
u = zeros(length(t), 1);
y = zeros(length(t), 1);

%% Sliding mode.
omega = zeros(length(t), 1);
V = zeros(length(t), 1);
rho = 2;
c = 1.5;

%% State space.
A = [0 1; 0 0]; % State transition matrix.
B = [0; 1]; % Input matrix. 
C = [1 1]; % Output matrix.
D = 0; % Feedforward matrix.
sys = ss(A, B, C, D);
sys_d = c2d(sys, t_s);
x_0 = [1; -2]; % Initial conditions.
x = zeros(length(x_0), length(t));
x(:, 1) = x_0;
f = @(t) [0; sin(2*t)]; % Input disturbance.
state_f = @(t, x, u) sys.a*x + sys.b*u + [0; sin(2*t)];
output_f = @(t, x, u) sys.c*x + sys.d*u;
for iStep=1:(length(t)-1)
    %% Simulation.
    [~, x_step] = ode45(@(t, x) state_f(t, x, u(iStep)), [t(iStep) t(iStep+1)], x(:, iStep));
    x(:, iStep+1) = x_step(end, :)';
    y(iStep, :) = sys_d.c*x(:, iStep) + sys_d.d*u(iStep, :);
    
    %% Control.
    omega(iStep) = [c 1]*x(:, iStep); % Sliding variable.
    V(iStep) = 0.5*omega(iStep)^2; % Lyapunov function.
    u(iStep+1) = -c*x(2, iStep) - rho*sign(omega(iStep));
end
