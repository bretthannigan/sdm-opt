%% Parameters.
order = 5;
osr = 32; % Oversampling ratio.
fs = 200; % Signal Nyquist frequeuncy (Hz).
DESIGN_TYPE = 1;
INITIAL_GUESS_TYPE = 2; % Choose the method used to generate an initial condition.

%% Define the optimzation goals.
switch DESIGN_TYPE
    case 1 % H-infinity design like in Thesis Section 5.1.
        termParam = struct('maxIter', 500, 'kappa', 0.01, 'epsilon', 1e-4);
        lee_criterion = 1.5; % Maximum NTF out-of-band gain.
        synOpt = defineSynOpt(-1, Inf, [0 pi/osr], 2, 2,... % Performance target.
            lee_criterion, Inf, [0 pi], 2, 2); % H-inf stability criterion (Lee's rule).
        k = 1; % Quantizer gain.
    case 2 % Root locus design like in Thesis Section 5.2.
        termParam = struct('maxIter', 2000, 'kappa', 0.005, 'epsilon', 1e-4);
        synOpt = defineSynOpt(-1, Inf, [0 pi/osr], 2, 2,... % Performance target.
            1, Inf, [0 pi], 1, 1); % Root locus stability criterion.
        k = ureal('k', 1, 'Range', [0.1 2]); % The upper range usually doesn't matter as long as it well above 1.
    case 3 % H-2 design like in Thesis Section 5.3.
        termParam = struct('maxIter', 150, 'kappa', 0, 'epsilon', 1e-6);
        h2_criterion = sqrt(1.79);
        synOpt = defineSynOpt(-1, Inf, [0 pi/osr], 2, 2,... % Performance target.
            h2_criterion, 2, [0 pi], 2, 2); % H-2 stability criterion.
        k = 1; % Quantizer gain.
    case 4 % l-1/star norm design like in Thesis Section 5.4.
        termParam = struct('maxIter', 100, 'kappa', 0.001, 'epsilon', 1e-6);
        star_criterion = 4;
        synOpt = defineSynOpt(-1, Inf, [0 pi/osr], 2, 2,... % Performance target.
            star_criterion, 1, [0 pi], 2, 2); % l-1 stability criterion.
        k = 1; % Quantizer gain.
end

%% Obtain initial condition.
switch INITIAL_GUESS_TYPE
    case 1
        % Let's use a simple initial condition transfer function:
        h_0 = zpk(zeros(order - 1, 1), zeros(order, 1), 1, 1/fs);
    case 2
        % Alternately, can use an initial guess from Schreier's DSToolbox, 
        % with a couple small changes.
        S_0 = synthesizeNTF(order, osr, 1);
        S_0.z{1} = 0.98*S_0.z{1}; % Slightly contract the zeros as they are located on the unit circle by default.
        h_0 = minreal((1 - S_0)/S_0); % Must be provided as a loop filter rather than NTF.
    case 3
        % The third option is to use a convex relaxation to generate an
        % initial transfer function, where the h_0 input to SDSYN is a
        % numeric scalar specifying the desired system order.
        [h_0, ~, ~, ~, ~] = sdsyn(5, synOpt, k, 'Display', 'off', 'TermParam', termParam, 'InitialGuess', true);
end

%% Run the optimization.
tic
[H, S, CL, normz, diagn, iterProgress] = sdsyn(h_0, synOpt, k, 'Display', 'on', 'TermParam', termParam);
toc