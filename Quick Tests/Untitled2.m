ntf_1 = synthesizeNTF(3, 8, 0);
r_2 = 0.9;
theta_2 = 15/180*pi;
lambda_2 = 1.05;
b_2 = 0.5;
r_3 = 0.9;
theta_3 = 15/180*pi;
lambda_3 = 1.05;
b_3 = 1;
h_1 = minreal((1 - ntf_1)/ntf_1);
h_2 = tf([0, b_2], [1 -lambda_2], 1, 'Variable', 'z^-1') + ...
    tf([0 2*r_2*cos(theta_2) -r_2^2],  [1 -2*r_2*cos(theta_2) r_2^2], 1, 'Variable', 'z^-1');
h_3 = tf([0, b_3], [1 -lambda_3], 1, 'Variable', 'z^-1') + ...
    tf([0 2*r_3*cos(theta_3) -r_3^2],  [1 -2*r_3*cos(theta_3) r_3^2], 1, 'Variable', 'z^-1');
ntf_2 = zpk(minreal(1/(1 + h_2)));
ntf_3 = zpk(minreal(1/(1 + h_3)));

figure
bode(h_1);
hold on
bode(h_2);
bode(h_3);
title('Bode Plot of Loop Filter')
legend('DSToolbox', 'Mladenov 1', 'Mladenov 2')

figure
bode(ntf_1);
hold on
bode(ntf_2);
bode(ntf_3);
title('Bode Plot of NTF')
legend('DSToolbox', 'Mladenov 1', 'Mladenov 2')

figure
pzplot(ntf_1);
hold on
pzplot(ntf_2);
pzplot(ntf_3);
title('Pole-Zero Plot of NTF')
legend('DSToolbox', 'Mladenov 1', 'Mladenov 2')

[snr_1, amp_1] = simulateSNR(ntf_1, 8);
[snr_2, amp_2] = simulateSNR(ntf_2, 8);
[snr_3, amp_3] = simulateSNR(ntf_3, 8);
figure
plot(amp_1, snr_1, '-.');
hold on
plot(amp_2, snr_2, '-.');
plot(amp_3, snr_3, '-.');
title('SNR Curve')
legend('DSToolbox', 'Mladenov 1', 'Mladenov 2')

u = 0.5*sin(2*pi*85/8192*(0:8191));
[v_1, xn_1] = simulateDSM(u, ntf_1);
[v_2, xn_2] = simulateDSM(u, ntf_2);
[v_3, xn_3] = simulateDSM(u, ntf_3);
figure
plot(u)
hold on
stairs(v_1)
stairs(v_2)
stairs(v_3)
title('Time Domain Output')
xlabel('Sample')
ylabel('Amplitude')
legend('Input', 'DSToolbox', 'Mladenov 1', 'Mladenov 2')
xlim([0 85])

figure
subplot(3,1,1)
plot(xn_1')
title('Internal States - DSToolbox')
xlabel('Sample')
ylabel('Amplitude')
ylim([-5 5])
subplot(3,1,2)
plot(xn_2')
title('Internal States - Mladenov 1')
xlabel('Sample')
ylabel('Amplitude')
ylim([-5 5])
subplot(3,1,3)
plot(xn_3')
title('Internal States - Mladenov 2')
xlabel('Sample')
ylabel('Amplitude')
ylim([-5 5])

b = 1/8*ones(1, 8);
o_1 = filter(b, 1, filter(b, 1, filter(b, 1, v_1)));
o_2 = filter(b, 1, filter(b, 1, filter(b, 1, v_2)));
o_3 = filter(b, 1, filter(b, 1, filter(b, 1, v_3)));
figure
plot(u)
hold on
stairs(o_1)
stairs(o_2)
stairs(o_3)
title('Output')
xlabel('Sample')
ylabel('Amplitude')
legend('Input', 'DSToolbox', 'Mladenov 1', 'Mladenov 2')