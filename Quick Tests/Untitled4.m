order = 4;
osr = 32;

starNorm = (2.2:0.1:10)';
n = length(starNorm);
l1Norm = zeros(n, 1);
l1NormOpt = zeros(n, 1);
hinfNorm = zeros(n, 1);
maxStabAmp = zeros(n, 1);

for ii=1:n
    
    [sys_tf, gam_gkyp, gam_l1] = l1gkypsyn_filter_2(order, [0 pi/osr], -1, starNorm(ii));
    hinfNorm(ii) = gam_gkyp;
    l1Norm(ii) = norm1(sys_tf);
    H = minreal((1 - sys_tf)/sys_tf);
    [test, l1NormOpt(ii)] = fminbnd(@(x) norm1(1/(1 + x*H)), 0, 10);
    [snr, amp] = simulateSNR(zpk(sys_tf), osr, [-100:0.1:0]);
    test = amp(find(snr>=(0.8*max(snr)),1 , 'last'));
    if ~isempty(test)
        maxStabAmp(ii) = test;
    else
        maxStabAmp(ii) = NaN;
    end
    
end