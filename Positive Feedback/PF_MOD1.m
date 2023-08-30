%PF_MOD1 Positive Feedback Sigma-Delta Modulator
%   
%   See also:
%       [1] A. M. Mallinson and S. J. Damphousse, “System and Method for 
%           Compensating for Error in a Sigma Delta Circuit,” 7259704 B2, 
%           2007.
%       [2] https://youtu.be/8_VKUVrVErA

%% Conventional sigma-delta modulator circuit.
% Implementing the sigma-delta modulator in electrical circuit form with an
% integrating op-amp. See [2] at 00:17:43 for the circuit.
R_1 = 1e3; % 1 kohm
R_2 = 1e3; % 1 kohm
C = 1e-9; % 1 nF
%stf = tf(-1, [R_1*C R_1/R_2]); % Inverting.
stf = tf(-R_2/R_1, [R_2*C, 1]);
ntf = tf([R_2*C 0], [R_2*C 1]);

%% Frequency response plot.
% Compared to [2], the peak magnitude of both transfer functions is lower.
% This is because [2] plots the output voltage when the input voltage is a
% summation of sinusoids with an RMS >1 whereas this script plots the
% transfer function ratio.
w = 2*pi*logspace(3, 6, 30);
[mag_stf, ~, wout] = bode(stf, w);
[mag_ntf, ~, ~] = bode(ntf, w);
figure;
subplot(2, 1, 1)
semilogx(wout/(2*pi), 20*log10(squeeze(mag_stf)))
hold on
semilogx(wout/(2*pi), 20*log10(squeeze(mag_ntf)))
legend('STF', 'NTF')
grid on
legend('STF', 'NTF', 'Location', 'SouthEast')
title('Frequency Response of Conventional Sigma Delta Modulator Circuit')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

%% Positive feedback sigma-delta modulator circuit.
% These are the transfer functions for the positive feedback modulator
% presented in [2].
R_1 = 1e3; % 1 kohm, input resistor.
R_2 = 1e3; % 1 kohm, negative feedback resistor.
R_3 = 999; % 999 ohm positive feedback resistor.
C = 1e-9; % 1 nF.
stf_pf = tf(-R_2, [(R_1-R_3)*R_2*C R_1]);
ntf_pf = tf([(R_1-R_3)*R_2*C 0], [(R_1-R_3)*R_2*C R_1]);
[mag_stf_pf, ~, ~] = bode(stf_pf, w);
[mag_ntf_pf, ~, ~] = bode(ntf_pf, w);
subplot(2, 1, 2)
semilogx(wout/(2*pi), 20*log10(squeeze(mag_stf_pf)))
hold on
semilogx(wout/(2*pi), 20*log10(squeeze(mag_ntf_pf)))
legend('STF', 'NTF')
grid on
legend('STF', 'NTF', 'Location', 'SouthEast')
title('Frequency Response of Positive Feedback Sigma Delta Modulator Circuit')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')