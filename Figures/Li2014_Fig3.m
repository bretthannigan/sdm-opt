synOpt = defineSynOpt(-1, Inf, [0 pi/32], 2, 2, 1.5, Inf, [0 pi], 2, 2);
gam_gkyp = {};

for ii=3:40
    validNTF = false;
    while ~validNTF
        h = drss(4);
        hs = 1/(1 + h);
        if hs.D==1 && isstable(hs)
            validNTF = true;
        end
    end
    
    [~, S{ii}, ~, ~, ~, iterProgress] = sdsyn(h, synOpt, 1, 'Display', 'off');
    iterProgress1 = vertcat(iterProgress(:).gam);
    gam_gkyp{ii} = iterProgress1(:, 1);
end
    
for ii=1:length(gam_gkyp)
    best_gam(ii) = gam_gkyp{ii}(end);
end
best_gam = mag2db(sqrt(best_gam));

for ii=1:40
    if ~isBest(ii)
        pzplot(S{ii})
        hold on
    end
end