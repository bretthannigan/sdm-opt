function [K, CL, normz] = my_l1syn_ext(varargin)
%MY_L1SYN Discrete-Time l1 Optimal Control Using Extended Simplified LMI 
% Approach
%   Generates optimal or sub-optimal controller in the l1 norm sense using 
%   the simplified LMI approach from [1]. Produces discrete-time 
%   controllers for discrete-time plants with the extended 
%
%   Inputs:
%       P: stabilizable and detectable real rational proper LTI plant.
%       NMEAS: number of measured variables or controller y inputs. 
%       (optional, default: 1 less than total plant outputs.)
%       NCON: number of controlled variables or controller u outputs. 
%       (optional, default: 1 less than total plant inputs).
%       L1MAX: if this name-value parameter is specified, a sub-optimal 
%       H-2 controller is designed to the value given. If unspecified or 
%       negative, an optimal H-2 controller is designed to minimize the H-2 
%       norm of the performance bound H2MAX. If H2MAX == Inf, constraint is
%       ignored.
%       DISPLAY: 'on' | 'off' (default), name-value parameter to specify if
%       the LMI solver should output progress information.
%       
%   Outputs:
%       K: H-infinity controller state-space model.
%       CL: closed loop system, equivalent to lft(P, K).
%       normz: H-infinity bound from w to z, minimized if GAM input was
%       unspecified.
%   
%   See also: H2HINFSYN, MY_H2HINFSYN
%       [1] A. Oberoi and J. C. Cockburn, “A simplified LMI approach to l1 
%           Controller Design,” Proc. 2005 Am. Control Conf., no. July, pp. 
%           1788–1792, 2005.
%       [2] M. C. De Oliveira, J. C. Geromel, and J. Bernussou, “Extended 
%           H2 and H? norm characterizations and controller 
%           parametrizations for discrete-time systems,” Int. J. Control, 
%           vol. 75, no. 9, pp. 666–679, 2002.
%
%   $Author: BH $    $Date: 2017-12-04 $    $Revision: 1 $

%% Deal with variable function inputs.
    p = inputParser;
    addRequired(p, 'P', @(x) isa(x, 'ss') || isa(x, 'tf'));
    parse(p, varargin{1});
    G = p.Results.P;
    if ~isa(G, 'ss') % Convert other LTI models to state space form.
        G = ss(G);
    end
    addOptional(p, 'NMEAS', max(size(G.C, 1) - 2, 1), @(x) isnumeric(x) && isscalar(x));
    addOptional(p, 'NCON', max(size(G.B, 2) - 2, 1), @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'L1MAX', -1, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'DISPLAY', 'off', @(x) ismember(x, {'on', 'off'}));
    addParameter(p, 'EXTENDED', true, @(x) isscalar(x) && islogical(x));
    parse(p, varargin{:});
    n = length(G.A);
    ny = p.Results.NMEAS;
    nz = size(G.C, 1) - ny;  
    nu = p.Results.NCON;
    nw = size(G.B, 1) - nu;
    psi_in = p.Results.L1MAX;
    traceOutput = strcmp(p.Results.DISPLAY, 'on');
    isExtended = p.Results.EXTENDED;

    if G.Ts == 0
        %error('my_h2hinfsyn2:ctplant', 'Continuous-time plants are not supported.');
    end
    
    %% Plant partitioning:
    %       In:    x      w      u
    %              |      |      |
    % Out:         V      V      V
    %          .--------------------.
    % dx/dt <- |   A  |  B_1 |  B_2 | n
    %          |------|------|------|
    %     z <- |  C_1 | D_11 | D_12 | nz
    %          |------|------|------|
    %     y <- |  C_2 | D_21 | D_22 | ny
    %          '--------------------'
    %      Size:   n     nw     nu
    A = G.A;
    C = mat2cell(G.C, [nz ny], size(G.C, 2));
    B = mat2cell(G.B, size(G.B, 1), [nw nu]);
    D = mat2cell(G.D, [nz ny], [nw nu]);
    
    %% Check Assumptions.
    P_2 = ss(A, B{2}, C{2}, D{2,2}, G.Ts); % State-space model for testing stability.
    [~, P_us] = stabsep(P_2);
    if rank(ctrb(P_us.A, P_us.B)) < size(P_us.A, 1) % Unstabilizable.
        error('my_hinfsyn:unstab', 'Stabilizability assumption is violated.');
    end
    if rank(obsv(P_us.A, P_us.C)) < size(P_us.A, 2) % Undetectable.
        error('my_hinfsyn:undet', 'Detectability assumption is violated.');
    end
    
    if psi_in == -1
        psi = sdpvar(1); % Optimization problem.
    else
        psi = psi_in.^2; % Sub-optimal (feasibility) problem.
    end
    
    function [linObjective, xfeas] = linOptimize(linValue, isVerbose)

        % Scalars:
        alpha = linValue;
        mu = sdpvar(1);
        nu_1 = sdpvar(1);
        % Matrices:
        P = sdpvar(n);
        H = sdpvar(n);
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

        ymLmi = {};
        ymConstraint = [];%[norm(X)<=1e3, norm(Y)<=1e3]; % Added to fix numerical problems.
        ymLmi{1} = blkvar;
        ymLmi{1}(1,1) = alpha*(X + X' - P); % alpha*(X + X' - P)
        ymLmi{1}(1,2) = alpha*(eye(n) + S' - J); % alpha*(I_n + S' - J)
        ymLmi{1}(1,3) = zeros(n, nw);
        ymLmi{1}(1,4) = X'*C{1}' + L'*D{1,2}'; % X'*C_z' + L'*D_zu'
        ymLmi{1}(2,2) = alpha*(Y + Y' - H); % alpha*(Y + Y' - H)
        ymLmi{1}(2,3) = zeros(n, nw);
        ymLmi{1}(2,4) = C{1}' + C{2}'*R'*D{1,2}'; % C_z' + C_y'*R'*D_zu'
        ymLmi{1}(3,3) = (mu - 1)*eye(nw); % (mu - 1)*I_nw
        ymLmi{1}(3,4) = D{1,1}' + D{2,1}'*R'*D{1,2}'; % D_zw' + D_yw'*R'*D_zu'
        ymLmi{1}(4,4) = nu_1*eye(nz); % nu*I_nz
        ymLmi{1} = sdpvar(ymLmi{1});
        ymConstraint = [ymConstraint, ymLmi{1}>=0];
        ymLmi{2} = blkvar;
        ymLmi{2}(1,1) = P; % P
        ymLmi{2}(1,2) = J; % J
        ymLmi{2}(1,3) = A*X + B{2}*L; % A*X + B_u*L
        ymLmi{2}(1,4) = A + B{2}*R*C{2}; % A + B_u*R*C_y
        ymLmi{2}(1,5) = B{1} + B{2}*R*D{2,1}; % B_w + B_u*R*D_yw
        ymLmi{2}(2,2) = H; % H
        ymLmi{2}(2,3) = Q; % Q
        ymLmi{2}(2,4) = Y*A + F*C{2}; % Y*A + F*C_y
        ymLmi{2}(2,5) = Y*B{1} + F*D{2,1}; % Y*B_w + F*D_yw
        ymLmi{2}(3,3) = (1 - alpha)*(X + X' - P); % (1 - alpha)*(X + X' - P)
        ymLmi{2}(3,4) = (1 - alpha)*(eye(n) + S' - J); % (1 - alpha)*(I_n + S' - J)
        ymLmi{2}(3,5) = zeros(n, nw);
        ymLmi{2}(4,4) = (1 - alpha)*(Y + Y' - H); % (1 - alpha)*(Y + Y' - H)
        ymLmi{2}(4,5) = zeros(n, nw);
        ymLmi{2}(5,5) = eye(nw); % I_nw
        ymLmi{2} = sdpvar(ymLmi{2});
        ymConstraint = [ymConstraint, ymLmi{2}>=0];
        ymLmi{3} = blkvar;
        ymLmi{3}(1,1) = psi;
        ymLmi{3}(1,2) = mu;
        ymLmi{3}(1,3) = nu_1;
        ymLmi{3}(2,2) = 1;
        ymLmi{3}(3,3) = 1;
        ymLmi{3} = sdpvar(ymLmi{3});
        ymConstraint = [ymConstraint, ymLmi{3}>=0];   

        %% Solve LMI Minimization Problem.
        if psi_in==-1
            ymObjective = psi;
        else
            ymObjective = [];
        end
        ymOptions = sdpsettings('verbose', isVerbose, 'solver', 'sedumi');
        xfeas = optimize(ymConstraint, ymObjective, ymOptions);
        linObjective = value(psi);
    end

    %% Line search over alpha.
    if traceOutput
        linOptions = optimset('Display', 'iter', 'TolX', 1e-2);
    else
        linOptions = optimset('Display', 'off', 'TolX', 1e-2);
    end
    linValue = fminbnd(@(x) linOptimize(x, false), 0, 1, linOptions);
    % Repeat last iteration to produce output.
    [~, xfeas] = linOptimize(linValue, traceOutput);
    if xfeas.problem==0
        error('my_l1syn:infeasible', xfeas.info);
    end
    
    psi = value(psi);
    X = value(X);
    S = value(S);
    Y = value(Y);
    L = value(L);
    F = value(F);
    Q = value(Q);
    R = value(R);
    
    %% Synthesize Controller.
    [U, Sigma, V] = svd(S - X*Y); % Unsure what to do if (eye(n) - X*Y) is not full rank.
    M = U*sqrt(Sigma);
  	M = M'; % Transposing M and N seems to generate stable CL system.
    N = sqrt(Sigma)*V';
    N = N';
    % Reconstruct controller state-space variables:
    K = [inv(N) -N\Y*B{2}; zeros(nu, n) eye(nu)]*...
        [Q - Y*A*X F; L R]*[inv(M) zeros(n, ny); -C{2}*X/M eye(ny)];
    K = ss(K(1:n, 1:n), K(1:n, n+1:n+ny), K(n+1:n+nu, 1:n), K(n+1:n+nu, n+1:n+ny), G.Ts);
    
    normz = sqrt(psi);
    CL = lft(G, K); % Close the loop (emulate output from H2HINFSYN.m).
end