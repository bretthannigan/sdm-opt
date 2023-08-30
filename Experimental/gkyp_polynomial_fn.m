function gkyp_polynomial_fn
% The following based on Corollary 1 in Section 4.2 and "Digital filter 
% design" in Section 7.1 of [Iwasaki2003a].

    n = 4;
    gam_s = db2mag(-20);
    gam_p = 1;
    w_s = [0 pi/32];
    w_p = [w_s(2)*10^(-mag2db(gam_s)/(20*n)) pi];

    if w_s(1)==0 % Low frequency case, Psi_d = [-2*cos(w(1)) 1; 1 0]
        w_s = [-w_s(2) w_s(2)];
    elseif w_s(2)==pi % High  frequency case, Psi_d = [2*cos(w(2)) -1; -1 0]
        w_s = [w_s(1) 2*pi-w_s(1)];
    end
    if w_p(1)==0 % Low frequency case, Psi_d = [-2*cos(w(1)) 1; 1 0]
        w_p = [-w_p(2) w_p(2)];
    elseif w_p(2)==pi % High  frequency case, Psi_d = [2*cos(w(2)) -1; -1 0]
        w_p = [w_p(1) 2*pi-w_p(1)];
    end

    wc_s = (w_s(2) + w_s(1))/2;
    w0_s = (w_s(2) - w_s(1))/2;
    wc_p = (w_p(2) + w_p(1))/2;
    w0_p = (w_p(2) - w_p(1))/2;
    Psi_s = [0 exp(1i*wc_s); exp(-1i*wc_s) -2*cos(w0_s)];
    Psi_p = [0 exp(1i*wc_p); exp(-1i*wc_p) -2*cos(w0_p)];

    B = sdpvar(ones(1, n+1), ones(1, n+1));
    A = [sdpvar(ones(1, n), ones(1, n)) 1];

    X_B = sdpvar(n+1, n+1, 'symmetric', 'complex');
    Q_B = sdpvar(n, n, 'symmetric', 'complex');

    X_A = sdpvar(n+1, n+1, 'symmetric', 'complex');
    Q_A = sdpvar(n, n, 'symmetric', 'complex');

    J = [1 zeros(1, n)];
    U = [zeros(n, 1) eye(n)];
    V = [eye(n) zeros(n, 1)];
    F = [U; V];

    Phi = [1 0; 0 -1];
    Psi_B = [0 1; 1 -2];
    Psi_A = Psi_B;

    lmi{1} = X_B + F'*kron(Psi_B, Q_B)*F;
    lmi{2} = X_A + F'*kron(Psi_A, Q_A)*F;

    ymConstraint = [];
    for iDiag=0:n
        ymConstraint = [ymConstraint, B{iDiag + 1} == sum(diag(lmi{1}, iDiag))]; %#ok<AGROW>
    end
    ymConstraint = [ymConstraint, X_B>=0, Q_B>=0];

    for iDiag=0:n
        ymConstraint = [ymConstraint, A{iDiag + 1} == sum(diag(lmi{2}, iDiag))]; %#ok<AGROW>
    end
    ymConstraint = [ymConstraint, X_A>=0, Q_A>=0];

    [lmi{3}, Q_p] = polynomial_lmi(B, A, gam_s, Psi_s);
    ymConstraint = [ymConstraint, lmi{3}<=0, Q_p>=0];

    [lmi{4}, Q_s] = polynomial_lmi(B, A, gam_p, Psi_p);
    ymConstraint = [ymConstraint, lmi{4}>=0, Q_s>=0];
    
    [lmi{5}, Q] = polynomial_lmi(B, A, 1.5, Psi_p);
    ymConstraint = [ymConstraint, lmi{5}<=0, Q>=0];
    
    ymSolver = 'sedumi'; %'SeDuMi';
    ymOptions = sdpsettings('verbose', true, 'solver', ymSolver);
    ymObjective = [];

    ymDiagnostics = optimize(ymConstraint, ymObjective, ymOptions);
    if ymDiagnostics.problem==1
        error('dtsyn:solver_infeasible', ymDiagnostics.info);
    elseif ymDiagnostics.problem~=0
        warning('on');
        warning(ymDiagnostics.info);
    end

    B = fliplr(cellfun(@value, B));
    A = fliplr(cellfun(@value, A));
    sys_tf = tf(B, A, 1);

    H_sys = tf([B B(end-1:-1:1)], [A A(end-1:-1:1)], 1);
    %H_sys = sys_tf + sys_tf';
    [G_sys, S_sys] = spectralfact(H_sys);

    function [lmi, Q] = polynomial_lmi(B, A, gam, Psi)
        P = sdpvar(n, n, 'symmetric', 'complex');
        Q = sdpvar(n, n, 'symmetric', 'complex');
        G = sdpvar(1, n+1);
        G(1) = (B{1} - gam^2*A{1})/2;
        for iG=2:(n+1)
            G(iG) = B{iG} - gam^2*A{iG}; %#ok<SAGROW>
        end
        lmi = F'*(kron(Phi, P) + kron(Psi, Q))*F - (G'*J + J'*G);
    end
end