function [K, CL, normz] = my_h2hinfsyn2(varargin)
%MY_H2HINFSYN2 Discrete-Time Mixed H-2/H-Infinity Multiobjective Control
%using Extended LMI Approach
%   Generates optimal or sub-optimal H-infinity controller using the
%   extended LMI technique from [1]. Produces discrete-time controllers for
%   discrete-time plants.
%
%   Inputs:
%       P: stabilizable and detectable real rational proper LTI plant.
%       NMEAS: number of measured variables or controller y inputs. 
%       (optional, default: 1 less than total plant outputs.)
%       NCON: number of controlled variables or controller u outputs. 
%       (optional, default: 1 less than total plant inputs).
%       N2: number of variables under H-2 norm constraint. The number of 
%       variables under H-Inf norm constraint is taken as the first 
%       Nin - N2(1) - Ncon input and first Nout - N2(2) - Nmeas output 
%       channels. If size(N2)==1, the input channel(s) are shared between 
%       constraints. (optional, default=1).
%       H2MAX: if this name-value parameter is specified, a sub-optimal 
%       H-2 controller is designed to the value given. If unspecified or 
%       negative, an optimal H-2 controller is designed to minimize the H-2 
%       norm of the performance bound H2MAX. If H2MAX == Inf, constraint is
%       ignored.
%       HINFMAX: if this name-value parameter is specified, a sub-optimal 
%       H-Inf controller is designed to the value given. If unspecified or 
%       negative, an optimal H-Inf controller is designed to minimize the 
%       H-Inf norm of the robustness bound HINFMAX. If HINFMAX == Inf,
%       constraint is ignored.
%       DISPLAY: 'on' | 'off' (default), name-value parameter to specify if
%       the LMI solver should output progress information.
%       EXTENDED: logical true (default) | false, name-value parameter to
%       specify if DeOliveira Extended method is used or the results shall
%       be reduced to the particular case of the Scherer (1997)
%       parameterization.
%       
%   Outputs:
%       K: H-infinity controller state-space model.
%       CL: closed loop system, equivalent to lft(P, K).
%       normz: H-infinity bound from w to z, minimized if GAM input was
%       unspecified.
%   
%   See also: H2HINFSYN, MY_H2HINFSYN
%       [1] M. C. De Oliveira, J. C. Geromel, and J. Bernussou, “Extended 
%           H2 and H? norm characterizations and controller 
%           parametrizations for discrete-time systems,” Int. J. Control, 
%           vol. 75, no. 9, pp. 666–679, 2002.
%       [3] EECE 508 Notes, 12. Multiobjective Control.
%
%   $Author: BH $    $Date: 2017-07-29 $    $Revision: 1 $
%
% REVISION 1:
%   2017-03-06 by BH: Fixed order of M, N matrices to match [2] and correct
%   issue where optimized normz was correct but the controller produced did
%   not satisfy the norm bound.

%% Deal with variable function inputs.
    p = inputParser;
    addRequired(p, 'P', @(x) isa(x, 'lti'));
    parse(p, varargin{1});
    G = p.Results.P;
    if ~isa(G, 'ss') % Convert other LTI models to state space form.
        G = ss(G);
    end
    addOptional(p, 'NMEAS', max(size(G.C, 1) - 2, 1), @(x) isnumeric(x) && isscalar(x));
    addOptional(p, 'NCON', max(size(G.B, 2) - 2, 1), @(x) isnumeric(x) && isscalar(x));
    addOptional(p, 'N2', 1, @(x) isnumeric(x));
    addParameter(p, 'HINFMAX', -1, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'H2MAX', -1, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'DISPLAY', 'off', @(x) ismember(x, {'on', 'off'}));
    addParameter(p, 'EXTENDED', true, @(x) isscalar(x) && islogical(x));
    parse(p, varargin{:});
    n = length(G.A);
    ny = p.Results.NMEAS;
    nz_2 = p.Results.N2(1);
    nz_inf = size(G.C, 1) - ny - nz_2;  
    nu = p.Results.NCON;
    gam_in = p.Results.HINFMAX;
    nu_in = p.Results.H2MAX;
    traceOutput = strcmp(p.Results.DISPLAY, 'on');
    isExtended = p.Results.EXTENDED;
    if length(p.Results.N2)==1 % input channels w are shared between constraints.
        nw_2 = size(G.B, 2) - nu;
        nw_inf = nw_2; 
    else
        nw_2 = p.Results.N2(2);
        nw_inf = size(G.B, 2) - nu - nw_2;
    end
    if nw_2==0 || nz_2==0
        nu_in = Inf;
    end
    if nw_inf==0 || nz_inf==0
        gam_in = Inf;
    end
    
    if G.Ts == 0
        %error('my_h2hinfsyn2:ctplant', 'Continuous-time plants are not supported.');
        [K, CL, normz] = my_h2hinfsyn(varargin{:}); % Continuous-time case.
        return
    end
    
    %% Plant partitioning:
    %       In:    x    w_inf   w_2     u
    %              |      |__/ __|      |
    % Out:         V      V      V      V   Size:
    %          .---------------------------.
    % dx/dt <- |   A  |  B_1 |  B_2 |  B_2 | n
    %          |------|------|------|------|
    % z_inf <- |  C_1 | D_11 | D_12 | D_13 | nz_inf
    %          |------|------|------|------|
    %   z_2 <- |  C_2 | D_21 | D_22 | D_23 | nz_2
    %          |------|------|------|------|
    %     y <- |  C_3 | D_31 | D_32 | D_33 | ny
    %          '---------------------------'
    %      Size:   n   nw_inf  nw_2    nu
    A = G.A;
    C = mat2cell(G.C, [nz_inf nz_2 ny], size(G.C, 2));
    if length(p.Results.N2)==1
        B = mat2cell(G.B, size(G.B, 1), [nw_inf nu]);
        B(3) = B(2);
        B(2) = B(1);
        D = mat2cell(G.D, [nz_inf nz_2 ny], [nw_inf nu]);
        D(:,3) = D(:,2);
        D(:,2) = D(:,1);
    else
        B = mat2cell(G.B, size(G.B, 1), [nw_inf nw_2 nu]);
        D = mat2cell(G.D, [nz_inf nz_2 ny], [nw_inf nw_2 nu]);
    end
    
    %% Check Assumptions.
    P_2 = ss(A, B{3}, C{3}, D{3,3}, G.Ts); % State-space model for testing stability.
    [~, P_us] = stabsep(P_2);
    if rank(ctrb(P_us.A, P_us.B)) < size(P_us.A, 1) % Unstabilizable.
        error('my_hinfsyn:unstab', 'Stabilizability assumption is violated.');
    end
    if rank(obsv(P_us.A, P_us.C)) < size(P_us.A, 2) % Undetectable.
        error('my_hinfsyn:undet', 'Detectability assumption is violated.');
    end

    setlmis([]);
    if gam_in == -1
        gam_inf = lmivar(1, [1 1]);
    end
    P = lmivar(1, [n 1]);
    H = lmivar(1, [n 1]);
    W = lmivar(1, [nz_2 1]);
    if isExtended
        X = lmivar(2, [n n]);
        S = lmivar(2, [n n]);
        J = lmivar(2, [n n]);
        Y = lmivar(2, [n n]);
    end
    L = lmivar(2, [nu n]);
    F = lmivar(2, [n ny]);
    Q = lmivar(2, [n n]);
    if nu_in ~= Inf
        R = -D{2,3}\D{2,2}/D{3,2}; % R is a linear matrix equality.
    else
        R = lmivar(2, [nu ny]);
    end
    
    %% Form H-2 output feedback LMIs.
    % Comments use notation from [1], for H-2 output feedback:
    %   B{2} = B_w, B{3} = B_u, C{2} = C_z, C{3} = C_y, D{2,2} = D_zw, 
    %   D{2,3} = D_zu, D{3,2} = D_yw, D{3,3} = D_yu.
    if nu_in ~= Inf
%       % LMI #1:
%       lmiterm([1 1 1 W], 1, 1);
%       if nu_in == -1 % Optimize nu.
%           lmiterm([-1 1 1 nu_2], 1, 1); % mu
%       else % Sub-optimal nu.
%           lmiterm([-1 1 1 0], nu_in); % mu
%       end
        % LMI #1:
        if nu_in ~= -1 % Sub-optimal nu.
            if nz_2>0
                diagElement = zeros(nz_2);
                diagElement(1) = 1;
            end
            for inz_2=1:nz_2
                lmiterm([1 1 1 W], circshift(diagElement, inz_2 - 1), circshift(diagElement, inz_2 - 1)'); % trace(W)
            end
            lmiterm([-1 1 1 0], nu_in^2); % mu
        end
        % LMI #2:
        lmiterm([-2 1 1 W], 1, 1); % W
        if isExtended
            lmiterm([-2 1 2 X], C{2}, 1); % C_z*X + ...
        else
            lmiterm([-2 1 2 P], C{2}, 1); % C_z*P + ...
        end
        lmiterm([-2 1 2 L], D{2,3}, 1); % ... + D_zu*L
        lmiterm([-2 1 3 0], C{2}); % C_z + ...
        %lmiterm([-2 1 3 R], D{2,3}, C{3}); % ... + D_zu*R*C_y
        lmiterm([-2 1 3 0], D{2,3}*R*C{3}); % ... + D_zu*R*C_y
        if isExtended
            lmiterm([-2 2 2 X], 1, 1, 's'); % X + X' + ...
            lmiterm([-2 2 2 P], -1, 1); % ... - P
            lmiterm([-2 2 3 0], 1); % I_n + ...
            lmiterm([-2 2 3 -S], 1, 1); % ... + S' + ...
            lmiterm([-2 2 3 J], -1, 1); % ... - J
            lmiterm([-2 3 3 Y], 1, 1, 's'); % Y + Y' + ...
            lmiterm([-2 3 3 H], -1, 1); % ... - H
        else
            lmiterm([-2 2 2 P], 1, 1); % P
            lmiterm([-2 2 3 0], 1); % I_n
            lmiterm([-2 3 3 H], 1, 1); % H
        end

        % LMI #3:
        lmiterm([-3 1 1 P], 1, 1); % P
        if isExtended
            lmiterm([-3 1 2 J], 1, 1); % J
            lmiterm([-3 1 3 X], A, 1); % A*X + ...
        else
            lmiterm([-3 1 2 0], 1); % I_n
            lmiterm([-3 1 3 P], A, 1); % A*P + ...
        end
        lmiterm([-3 1 3 L], B{3}, 1); % ... + B_u*L
        %lmiterm([-3 1 4 0], A); % A + ...
        %lmiterm([-3 1 4 R], B{3}, C{3}); % ... + B_u*R*C_y
        lmiterm([-3 1 4 0], A + B{3}*R*C{3}); % A + B_u*R*C_y
        %lmiterm([-3 1 5 0], B{2}); % B_w + ...
        %lmiterm([-3 1 5 R], B{3}, D{3,2}); % ... + B_u*R*D_yw
        lmiterm([-3 1 5 0], B{2} + B{3}*R*D{3,2}); % B_w + B_u*R*D_yw
        lmiterm([-3 2 2 H], 1, 1); % H
        lmiterm([-3 2 3 Q], 1, 1); % Q
        if isExtended
            lmiterm([-3 2 4 Y], 1, A); % Y*A + ...
        else
            lmiterm([-3 2 4 H], 1, A); % H*A + ...
        end
        lmiterm([-3 2 4 F], 1, C{3}); % ... + F*C_y
        if isExtended
            lmiterm([-3 2 5 Y], 1, B{2}); % Y*B_w + ...
        else
            lmiterm([-3 2 5 H], 1, B{2}); % H*B_w + ...
        end
        lmiterm([-3 2 5 F], 1, D{3,2}); % ... + F*D_yw
        if isExtended
            lmiterm([-3 3 3 X], 1, 1, 's'); % X + X' + ...
            lmiterm([-3 3 3 P], -1, 1); % ... - P
            lmiterm([-3 3 4 0], 1); % I_n + ...
            lmiterm([-3 3 4 -S], 1, 1); % + S' + ...
            lmiterm([-3 3 4 J], -1, 1); % ... - J
            lmiterm([-3 4 4 Y], 1, 1, 's'); % Y + Y' + ...
            lmiterm([-3 4 4 H], -1, 1); % ... - H
        else
            lmiterm([-3 3 3 P], 1, 1); % P
            lmiterm([-3 3 4 0], 1); % I_n
            lmiterm([-3 4 4 H], 1, 1); % H
        end
        lmiterm([-3 5 5 0], 1); % I_nw_2
    end
    
    %% Form H-inf output feedback LMIs.
    % Comments use notation from [1], for H-inf output feedback:
    %   B{1} = B_w, B{3} = B_u, C{1} = C_z, C{3} = C_y, D{1,1} = D_zw, 
    %   D{1,3} = D_zu, D{3,1} = D_yw, D{3,3} = D_yu.
    
    if gam_in ~= Inf
        % LMI #4:
        lmiterm([-4 1 1 P], 1, 1); % P
        if isExtended
            lmiterm([-4 1 2 J], 1, 1); % J
            lmiterm([-4 1 3 X], A, 1); % A*X + ...
        else
            lmiterm([-4 1 2 0], 1); % I_n
            lmiterm([-4 1 3 P], A, 1); % A*P + ...
        end
        lmiterm([-4 1 3 L], B{3}, 1); % ... + B_u*L
        if nu_in ~= Inf
            lmiterm([-4 1 4 0], A + B{3}*R*C{3}); % A + B_u*R*C_y
            lmiterm([-4 1 5 0], B{1} + B{3}*R*D{3,1}); % B_w + B_u*R*D_yw
        else
            lmiterm([-4 1 4 0], A); % A + ...
            lmiterm([-4 1 4 R], B{3}, C{3}); % ... + B_u*R*C_y
            lmiterm([-4 1 5 0], B{1}); % B_w + ...
            lmiterm([-4 1 5 R], B{3}, D{3,1}); % ... + B_u*R*D_yw
        end
        lmiterm([-4 2 2 H], 1, 1); % H
        lmiterm([-4 2 3 Q], 1, 1); % Q
        if isExtended
            lmiterm([-4 2 4 Y], 1, A); % Y*A + ...
        else
            lmiterm([-4 2 4 H], 1, A); % H*A + ...
        end
        lmiterm([-4 2 4 F], 1, C{3}); % ... + F*C_y
        if isExtended
            lmiterm([-4 2 5 Y], 1, B{1}); % Y*B_w + ...
        else
            lmiterm([-4 2 5 H], 1, B{1}); % H*B_w + ...
        end
        lmiterm([-4 2 5 F], 1, D{3,1}); % ... + F*D_yw
        if isExtended
            lmiterm([-4 3 3 X], 1, 1, 's'); % X + X' + ...
            lmiterm([-4 3 3 P], -1, 1); % ... - P
            lmiterm([-4 3 4 0], 1); % I_n + ...
            lmiterm([-4 3 4 -S], 1, 1); % ... + S' + ...
            lmiterm([-4 3 4 J], -1, 1); % ... - J
            lmiterm([-4 3 6 -X], 1, C{1}'); % X'*C_z' + ...
        else
            lmiterm([-4 3 3 P], 1, 1); % P
            lmiterm([-4 3 4 0], 1); % I_n
            lmiterm([-4 3 6 P], 1, C{1}'); % P*C_z' + ...
        end
        lmiterm([-4 3 6 -L], 1, D{1,3}'); % ... + L'*D_zu'
        if isExtended
            lmiterm([-4 4 4 Y], 1, 1, 's'); % Y + Y' + ...
            lmiterm([-4 4 4 H], -1, 1); % ... - H
        else
            lmiterm([-4 4 4 H], 1, 1); % H
        end
        if nu_in ~= Inf
            lmiterm([-4 4 6 0], C{1}' + C{3}'*R'*D{1,3}'); % C_z' + C_y'*R'*D_zu'
            lmiterm([-4 5 5 0], 1); % I_nw_inf
            lmiterm([-4 5 6 0], D{1,1}' + D{3,1}'*R'*D{1,3}'); % D_zw' + D_yw'*R'*D_zu'
        else
            lmiterm([-4 4 6 0], C{1}'); % C_z' + ...
            lmiterm([-4 4 6 -R], C{3}', D{1,3}'); % ... + C_y'*R'*D_zu'
            lmiterm([-4 5 5 0], 1); % I_nw_inf
            lmiterm([-4 5 6 0], D{1,1}'); % D_zw' + ...
            lmiterm([-4 5 6 -R], D{3,1}', D{1,3}'); % ... + D_yw'*R'*D_zu'
        end
        if gam_in == -1 % Optimize gamma.
            lmiterm([-4 6 6 gam_inf], 1, eye(nz_inf)); % mu*I_nz_inf
        else % Sub-optimal gamma.
            lmiterm([-4 6 6 0], gam_in^2*eye(nz_inf)); % mu*I_nz_inf
        end
    end
    
    %% Solve LMI Minimization Problem.
    lmisys = getlmis;
    opt = [1e-5 0 0 0 ~traceOutput];
    if nu_in == Inf 
        isOptimizedVariable = num2cell(zeros(1,12));
    else
        isOptimizedVariable = num2cell(zeros(1,11));
    end
    if nu_in == -1
        isOptimizedVariable(4) = {eye(nz_2)};
    end
    if ~isExtended
        isOptimizedVariable(5:8) = [];
    end
    if gam_in == -1
        isOptimizedVariable{1} = 1;
    else
        isOptimizedVariable(1) = [];
    end
    c = mat2dec(lmisys, isOptimizedVariable{:});
    if ~any(c)
        [~, xfeas] = feasp(lmisys, opt);
    else
        [~, xfeas] = mincx(lmisys, c, opt);
    end
    if gam_in == -1
        gam_inf = dec2mat(lmisys, xfeas, gam_inf);
    else
        gam_inf = gam_in.^2;
    end
    if nu_in == -1
        W = dec2mat(lmisys, xfeas, W);
        nu_2 = trace(W);
        %nu_2 = dec2mat(lmisys, xfeas, nu_2);
    else
        nu_2 = nu_in.^2;
    end
    P = dec2mat(lmisys, xfeas, P);
    H = dec2mat(lmisys, xfeas, H);
    if isExtended
        X = dec2mat(lmisys, xfeas, X);
        S = dec2mat(lmisys, xfeas, S);
        Y = dec2mat(lmisys, xfeas, Y);
    else
        X = P;
        S = eye(n);
        Y = H;
    end
    L = dec2mat(lmisys, xfeas, L);
    F = dec2mat(lmisys, xfeas, F);
    Q = dec2mat(lmisys, xfeas, Q);
    if nu_in == Inf
        R = dec2mat(lmisys, xfeas, R);
    end
    
    %% Synthesize Controller.
    [U, Sigma, V] = svd(S - Y*X);
    N = U*sqrt(Sigma);
    M = sqrt(Sigma)*V';
    % Reconstruct controller state-space variables:
    K = [inv(N) -N\Y*B{3}; zeros(nu, n) eye(nu)]*...
        [Q - Y*A*X F; L R]*[inv(M) zeros(n, ny); -C{3}*X/M eye(ny)];
    K = ss(K(1:n, 1:n), K(1:n, n+1:n+ny), K(n+1:n+nu, 1:n), K(n+1:n+nu, n+1:n+ny), G.Ts);
    
    normz = sqrt([gam_inf nu_2]);
    CL = lft(G, K); % Close the loop (emulate output from H2HINFSYN.m).
end