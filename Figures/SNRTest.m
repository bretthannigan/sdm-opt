osr = 64;

H = zpk([], 1, 1, 1);
H_p = zpk(0.5469, [1 0.9942], 2.5185, 1);

NTF = 1/(1 + H);
NTF_2 = synthesizeNTF(2, 64, 1);
NTF_p = 1/(1 + H_p);

f = logspace(-4, log10(1/(osr*2.1)), 40);
a = -logspace(2, 0, 60);

snr = zeros(length(a), length(f));
snr_2 = zeros(length(a), length(f));
snr_p = zeros(length(a), length(f));

for iFreq=1:length(f)
    [snr(:, iFreq), ~] = simulateSNR(NTF, 64, a, 0, 2, f(iFreq));
    [snr_2(:, iFreq), ~] = simulateSNR(NTF_2, 64, a, 0, 2, f(iFreq));
    [snr_p(:, iFreq), ~] = simulateSNR(NTF_p, 64, a, 0, 2, f(iFreq));
%     snr(iFreq) = max(snr(:, iFreq));
%     snr_p(iFreq) = max(s1_p);
end

figure
semilogx(f*(2*pi*512000), max(snr), '-s');
hold on
semilogx(f*(2*pi*512000), max(snr_p), '-s');
semilogx(f*(2*pi*512000), max(snr_2), '-s');
legend('MOD1', 'US7259704');