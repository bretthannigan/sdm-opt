order = 4;
osr = 32;

oobGain = (1:0.1:10)';
n = length(oobGain);
hinfNormOob = zeros(n, 1);
maxStabAmp = zeros(n, 1);

for ii=1:n
    
    [sys_tf, gam_gkyp, gam_l1] = hinfgkypsyn_filter(order, [0 pi/osr], -1, oobGain(ii));
    hinfNorm(ii) = gam_gkyp;
    l1NormOob(ii) = norm(sys_tf, Inf);
    [snr, amp] = simulateSNR(zpk(sys_tf), osr, [-100:0.1:0]);
    test = amp(find(snr>=(0.8*max(snr)),1 , 'last'));
    if ~isempty(test)
        maxStabAmp(ii) = test;
    else
        maxStabAmp(ii) = NaN;
    end
    
end