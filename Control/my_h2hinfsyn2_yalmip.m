function [K, CL, normz] = my_h2hinfsyn2_yalmip(varargin)
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
%       N2: 2 element numeric vector of number of variables under H-2 norm 
%       constraint. The number of variables under H-Inf norm constraint is 
%       taken as the first Nin - N2(1) - Ncon input and first 
%       Nout - N2(2) - Nmeas output channels. (optional, default: 
%       [Nin - Ncon, 1]).
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
%       [2] EECE 508 Notes, 12. Multiobjective Control.
%
%   $Author: BH $    $Date: 2017-07-29 $    $Revision: 1 $
%
% REVISION 1:
%   2017-03-06 by BH: Fixed order of M, N matrices to match [2] and correct
%   issue where optimized normz was correct but the controller produced did
%   not satisfy the norm bound.

%% Deal with variable function inputs.
    p = inputParser;
    addRequired(p, 'P', @(x) isa(x, 'ss') || isa(x, 'tf'));
    parse(p, varargin{1});
    G = p.Results.P;
    addOptional(p, 'NMEAS', max(size(G.C, 1) - 2, 1), @(x) isnumeric(x) && isscalar(x));
    addOptional(p, 'NCON', max(size(G.B, 2) - 2, 1), @(x) isnumeric(x) && isscalar(x));
    addOptional(p, 'N2', [0 1], @(x) isnumeric(x) && length(x)==2);
    addParameter(p, 'HINFMAX', -1, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'H2MAX', -1, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'DISPLAY', 'off', @(x) ismember(x, {'on', 'off'}));
    addParameter(p, 'EXTENDED', true, @(x) isscalar(x) && islogical(x));
    addParameter(p, 'COMBINED', [false false], @(x) length(x)==2 && all(islogical(x)));
    parse(p, varargin{:});
    n = length(G.A);
    ny = p.Results.NMEAS;
    nu = p.Results.NCON;
    gam_in = p.Results.HINFMAX;
    nu_in = p.Results.H2MAX;
    traceOutput = strcmp(p.Results.DISPLAY, 'on');
    isExtended = p.Results.EXTENDED;
    isCombined = p.Results.COMBINED;
    nw_2 = p.Results.N2(1)*~isCombined(1); % Will equal 0 if combined.
    nz_2 = p.Results.N2(2)*~isCombined(2); % Will equal 0 if combined.
    nw_inf = size(G.B, 2) - nu - nw_2;
    nz_inf = size(G.C, 1) - ny - nz_2;
    
    if ~isa(G, 'ss') % Convert other LTI models to state space form.
        G = ss(G);
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
    B = mat2cell(G.B, size(G.B, 1), [nw_inf nw_2 nu]);
    C = mat2cell(G.C, [nz_inf nz_2 ny], size(G.C, 2));
    D = mat2cell(G.D, [nz_inf nz_2 ny], [nw_inf nw_2 nu]);
    if isCombined(1) % Input channels are shared between constraints.
        B(2) = B(1);
        D(:,2) = D(:,1);
        nw_2 = nw_inf;
    end
    if isCombined(2) % Output channels are shared between constraints.
        C(2) = C(1);
        D(2,:) = D(1,:);
        nz_2 = nz_inf;
    end
    
    %% Check assumptions.
    P_2 = ss(A, B{3}, C{3}, D{3,3}, G.Ts); % State-space model for testing stability.
    [~, P_us] = stabsep(P_2);
    if rank(ctrb(P_us.A, P_us.B)) < size(P_us.A, 1) % Unstabilizable.
        error('my_hinfsyn:unstab', 'Stabilizability assumption is violated.');
    end
    if rank(obsv(P_us.A, P_us.C)) < size(P_us.A, 2) % Undetectable.
        error('my_hinfsyn:undet', 'Detectability assumption is violated.');
    end
    if nw_inf==0 || nz_inf==0
        gam_in = Inf; % Do not optimize H-Inf norm.
    end
    if nw_2==0 || nz_2==0
        nu_in = Inf; % Do not optimize H-2 norm.
    end

    %% Define optimization variables.
    if gam_in == -1
        gam_inf = sdpvar(1); % Optimization problem.
    else
        gam_inf = gam_in^2; % Sub-optimal (feasibility) problem.
    end
    if nu_in == -1
        nu_2 = sdpvar(1); % Optimization problem.
    else
        nu_2 = nu_in^2; % Sub-optimal (feasibility) problem.
    end
    P = sdpvar(n);
    H = sdpvar(n);
    W = sdpvar(nz_2);
    if isExtended
        X = sdpvar(n, n, 'full');
        S = sdpvar(n, n, 'full');
        J = sdpvar(n, n, 'full');
        Y = sdpvar(n, n, 'full');
    else
        X = P;
        S = eye(n);
        J = eye(n);
        Y = H;
    end
    L = sdpvar(nu, n, 'full');
    F = sdpvar(n, ny, 'full');
    Q = sdpvar(n, n, 'full');
    R = sdpvar(nu, ny, 'full');
    
    %% Form H2 output feedback LMIs.
    % Comments use notation from [1], for H-2 output feedback:
    %   B{2} = B_w, B{3} = B_u, C{2} = C_z, C{3} = C_y, D{2,2} = D_zw, 
    %   D{2,3} = D_zu, D{3,2} = D_yw, D{3,3} = D_yu.
    ymLmi = {};
    ymConstraint = [];
    if nu_in ~= Inf
        % LMI #1:
        ymLmi{1} = trace(W);
        ymConstraint = [ymConstraint, ymLmi{1} <= nu_2];
        % LMI #2:
%         ymLmi{2} = [W C{2}*X+D{2,3}*L C{2}+D{2,3}*R*C{3};...
%             X'*C{2}'+L'*D{2,3}' X+X'-P eye(n)+S'-J;...
%             C{2}'+C{3}'
        ymLmi{2} = blkvar;
        ymLmi{2}(1,1) = W; % W
        ymLmi{2}(1,2) = C{2}*X + D{2,3}*L; % C_z*X + D_zu*L
        ymLmi{2}(1,3) = C{2} + D{2,3}*R*C{4}; % C_z + D_zu*R*C_y
        ymLmi{2}(2,2) = X + X' - P; % X + X' - P
        ymLmi{2}(2,3) = eye(n) + S' - J; % I_n + S' - J
        ymLmi{2}(3,3) = Y + Y' - H; % Y + Y' - H
        ymLmi{2} = sdpvar(ymLmi{2});
        ymConstraint = [ymConstraint, ymLmi{2} >= 0];
        % LMI #3:
        ymLmi{3} = blkvar;
        ymLmi{3}(1,1) = P; % P
        ymLmi{3}(1,2) = J; % J
        ymLmi{3}(1,3) = A*X + B{3}*L; % A*X + B_u*L
        ymLmi{3}(1,4) = A + B{3}*R*C{3}; % A + B_u*R*C_y
        ymLmi{3}(1,5) = B{2} + B{3}*R*D{3,2}; % B_w + B_u*R*D_yw
        ymLmi{3}(2,2) = H; % H
        ymLmi{3}(2,3) = Q; % Q
        ymLmi{3}(2,4) = Y*A + F*C{3}; % Y*A + F*C_y
        ymLmi{3}(2,5) = Y*B{2} + F*D{3,2}; % Y*B_w + F*D_yw
        ymLmi{3}(3,3) = X + X' - P; % X + X' - P
        ymLmi{3}(3,4) = eye(n) + S' - J; % I_n + S' - J
        ymLmi{3}(4,4) = Y + Y' - H; % Y + Y' - H
        ymLmi{3}(5,5) = eye(nw_2); % I_nw_2
        ymLmi{3} = sdpvar(ymLmi{3});
        ymConstraint = [ymConstraint, ymLmi{3} >= 0];
        % LMI #4:
        ymLmi{4} = D{2,2} + D{2,3}*R*D{3,2}; % D_zw + D_zu*R*D_yw
        ymConstraint = [ymConstraint, ymLmi{4} == zeros(nz_2, nw_2)];
    end
    
    %% Form HInf output feedback LMIs.
    % Comments use notation from [1], for HInf output feedback:
    %   B{1} = B_w, B{3} = B_u, C{1} = C_z, C{3} = C_y, D{1,1} = D_zw, 
    %   D{1,3} = D_zu, D{3,1} = D_yw, D{3,3} = D_yu.
    if gam_in ~= Inf
        % LMI #4:
        ymLmi{5} = blkvar;
        ymLmi{5}(1,1) = P; % P
        ymLmi{5}(1,2) = J; % J
        ymLmi{5}(1,3) = A*X + B{3}*L; % A*X + B_u*L
        ymLmi{5}(1,4) = A + B{3}*R*C{3}; % A + B_u*R*C_y
        ymLmi{5}(1,5) = B{1} + B{3}*R*D{3,1}; % B_w + B_u*R*D_yw
        ymLmi{5}(2,2) = H; % H
        ymLmi{5}(2,3) = Q; % Q
        ymLmi{5}(2,4) = Y*A + F*C{3}; % Y*A + F*C_y
        ymLmi{5}(2,5) = Y*B{1} + F*D{3,1}; % Y*B_w + F*D_yw
        ymLmi{5}(3,3) = X + X' - P; % X + X' - P
        ymLmi{5}(3,4) = eye(n) + S' - J; % I_n + S' - J
        ymLmi{5}(3,6) = X'*C{1}' + L'*D{1,3}'; % X'*C_z' + L'*D_zu'
        ymLmi{5}(4,4) = Y + Y' - H; % Y + Y' - H;
        ymLmi{5}(4,6) = C{1}' + C{3}'*R'*D{1,3}'; % C_z' + C_y'*R'*D_zu'
        ymLmi{5}(5,5) = eye(nw_inf); % I_nw_inf
        ymLmi{5}(5,6) = D{1,1}' + D{3,1}'*R'*D{1,3}'; % D_zw' + D_yw'*R'*D_zu'
        ymLmi{5}(6,6) = gam_inf*eye(nz_inf); % mu*I_nz_inf
        ymLmi{5} = sdpvar(ymLmi{5});
        ymConstraint = [ymConstraint, ymLmi{5} >= 0];
    end
    
    %% Solve LMI Minimization Problem.
    if gam_in==-1 && nu_in==-1
        ymObjective = gam_inf + nu_2;
    elseif gam_in==-1
        ymObjective = gam_inf;
    elseif nu_in==-1
        ymObjective = nu_2;
    else
        ymObjective = [];
    end
    ymOptions = sdpsettings('verbose', traceOutput, 'solver', 'sdpt3');
    [model, ~] = export(ymConstraint, ymObjective, ymOptions);
    xfeas = optimize(ymConstraint, ymObjective, ymOptions);
    if xfeas.problem==0
        
    else
        xfeas.info
    end

    X = value(X);
    S = value(S);
    Y = value(Y);
    L = value(L);
    F = value(F);
    Q = value(Q);
    R = value(R);
    
    %% Synthesize Controller.
    [U, Sigma, V] = svd(S - Y*X);
    N = U*sqrt(Sigma);
    M = sqrt(Sigma)*V';
    % Reconstruct controller state-space variables:
    K = [inv(N) -N\Y*B{3}; zeros(nu, n) eye(nu)]*...
        [Q - Y*A*X F; L R]*[inv(M) zeros(n, ny); -C{3}*X/M eye(ny)];
    K = ss(K(1:n, 1:n), K(1:n, n+1:n+ny), K(n+1:n+nu, 1:n), K(n+1:n+nu, n+1:n+ny), G.Ts);
    
    normz = sqrt([value(gam_inf) value(nu_2)]);
    CL = lft(G, K); % Close the loop (emulate output from H2HINFSYN.m).
end