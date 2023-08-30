%OBEROI2004_42 Example from Section 2 of [1].
%   Script used to compare the derivation of the generalized plant shown in
%   state-space form in the thesis to my derivation.
%
%   See also:
%       [1] A. Oberoi, “A Convex Optimization Approach to the Design of 
%           Multiobjective Discrete Time Systems,” Rochester Institute of 
%           Technology, 2004. p. 54.

T_s = 1e-6;
OSR = 10;

%% State-space generalized plant from p. 56-57
A_42 = [0.5335 0 0; 0.2214 0.7304 0; 0 0 0.0019];
B_42 = [0.4989 0 0; 0.2070 -0.4438 -0.4438; 0 0 0.6995];
C_42 = [0.1345 0.4438 0; 0 0 -0.6995; 0.4989 0 0];
D_42 = [0.1258 -0.2696 -0.2696; 0 0 0.5217; 0.4665 -1 0];
G_42 = ss(A_42, B_42, C_42, D_42, T_s);

%% My derivation of the generalized plant.
W_r(1,1) = tf([0.46651 0], [1 -0.5335], T_s);
W_r(2,2) = tf(1, 1, T_s);
W_r(3,3) = tf(1, 1, T_s);
W_e = tf([0.2696 0], [1 -0.7304], T_s);
W_u(1) = tf(0, 1, T_s);
W_u(2) = tf([0.5217 -0.9397*0.5217], [1 -0.001867], T_s);
P = [tf(1, 1, T_s) tf(1, 1, T_s)];
P.InputName{1} = 'd';
P.InputName{2} = 'u';
P.OutputName{1} = 'y';
G = generalizeplant(ss(P), ss(W_e), ss(W_u), []);
G = series(ss(W_r), G);
% The results are the same with the exception of the feedthrough term D_eu.
% [Oberoi2004] omits the feedthrough but no "loop shifting" seems to have
% been done as claimed in the thesis.

%% Controllers from p. 67
K1 = zpk([0.4498 0.001763 0.0002342], [0.5947 0.7304 0.3676], 0.47794, T_s); % H2-HInf with lambda = 0.9.
K2 = zpk([0.4664 0.002027 -0.0007528], [0.4716 0.7304 0.02155], 0.7910, T_s); % H2-HInf with lambda = 0.1.
K3 = zpk([0.4566 0.001831 -0.0001911], [0.5077 0.7304 0.2655], 0.6251, T_s); % H2-HInf with lambda = 0.5.
K4 = zpk([-0.9816 0.472 0.001867], [0.08013 0.7304 0.4962], 0.25669, T_s); % l1(*).

G_42_1 = lft(G_42, K1);
G_42_2 = lft(G_42, K2);
G_42_3 = lft(G_42, K3);
G_42_4 = lft(G_42, K4);