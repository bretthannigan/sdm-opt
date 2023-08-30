u_max = 0.5;
alpha = ureal('alpha', 4, 'range', [1/u_max 10]);
P = c2d(ss(tf(1)), 1);
P = alpha*P;
W_S = makeweight(0.01, 1e-3, 30, 1); %ss(makeweight(0.5, 1e-7, 10, 1));
W_R = tf(u_max);
W_T = [];%1/makeweight(0.9, 1e-3, 5, 1);
P_gen = generalizeplant(P, W_S, W_R, W_T, 'UncertaintySz', 0);

[K, CL, gam] = my_hinfsyn(P_gen, 1, 1, 'DISPLAY', 'on'); % P_mu
bodemag(CL)
H = K*P;
bodemag(1/(1 + H), 'r') % Sensitivity function.
hold on
if ~isempty(W_S)
    bodemag(W_S, '--r') % Sensitivity specification.
end
bodemag(K/(1 + H), 'g') % Control transfer function.
if ~isempty(W_R)
    bodemag(W_R, '--g') % Control specification.
end
bodemag(H/(1 + H), 'b') % Complementary sensitivity function.
if ~isempty(W_T)
    bodemag(W_T, '--b') % Complementary sensitivity specification.
end
legend('Sensitivity', 'Desired sensitivity', 'Control', 'Desired control', 'Complementary sensitivity', 'Desired complementary sensitivity')
isstable(feedback(P, K))