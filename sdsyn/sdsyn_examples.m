%SDSYN_EXAMPLES Script to generate similar sigma delta filters to the
%thesis and paper.
%
%   See also:
%       [1] B. C. Hannigan, "On the Design of Stable, High Performance 
%           Sigma Delta Modulators", MASc. thesis, University of British 
%           Columbia, 2018.
%       [2] B. C. Hannigan, C. L. Petersen, A. M. Mallinson, and G. A. 
%           Dumont, "An Optimization Framework for the Design of Noise 
%           Shaping Loop Filters with Improved Stability Properties", 
%           Circuits, Systems, and Signal Processing, vol. 39, no. 12, pp. 
%           6276–6298, 2020.

%% H-Inf Stability Criterion 
% Reproduces a design similar to that from Section 5.1 of [1] and Section
% 4.1 of [2].

synOpt_51 = defineSynOpt(-1, Inf, [0 pi/32], 2, 2, 1.5, Inf, [0 pi], 2, 2);
S0 = synthesizeNTF(5, 32, 1, 0.5, 0); % Start with an initial guess NTF from the Delta-Sigma Toolbox.
H0 = zpk(minreal((1 - S0)/S0)); % Convert the NTF to loop filter form using Equation (1.3) from [1].
H0.p{1} = H0.p{1}*0.995; % Slightly scale the poles so that the loop filter is strictly stable.
[H, S, CL, normz, diagn, iterProgress] = sdsyn(H0, synOpt_51, 'display', 'on');
[snr, amp] = predictSNR(S, 32);
figure; 
plot(amp, snr, '.-k')
legend('location', 'NorthWest')
title('SQNR vs. Amplitude for H-Infinity Design')
xlabel('Amplitude (dB)')
ylabel('SQNR (dB)')

%% Root Locus Stability Criterion
% Reproduces a design similar to that from Section 5.2 of [1] and 4.2 of
% [2].


