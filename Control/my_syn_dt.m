function [K, CL, normz, ymDiagnostics] = my_syn_dt(varargin)
%MY_SYN Discrete-Time Multiobjective Optimal Control Using Extended LMI 
% Approach
%   Generates optimal or sub-optimal controller in the l1, H2, or Hinf 
%   sense using the simplified LMI approach from [1], [2]. Produces 
%   discrete-time controllers for discrete-time plants.
%
%   Inputs:
%       P: stabilizable and detectable real rational proper LTI plant.
%       NMEAS: number of measured variables or controller y inputs. 
%       (optional, default: 1 less than total plant outputs.)
%       NCON: number of controlled variables or controller u outputs. 
%       (optional, default: 1 less than total plant inputs).
%       W1: vector indicating which generalized plant disturbance input
%       channels are under l1 constraint. (optional, default: [])
%       W2: vector indicating which generalized plant disturbance input
%       channels are under H2 constraint. (optional, default: [])
%       WINF: vector indicating which generalized plant disturbance input
%       channels are under HInf constraint. (optional, default: [])
%       Z1: vector indicating which generalized plant performance output
%       channels are under l1 constraint. (optional, default: [])
%       Z2: vector indicating which generalized plant performance output
%       channels are under H2 constraint. (optional, default: [])
%       ZINF: vector indicating which generalized plant performance output
%       channels are under HInf constraint. (optional, default: [])
%       L1MAX: if this name-value parameter is specified, a sub-optimal 
%       l1 controller is designed to the value given. If unspecified or 
%       negative, an optimal l1 controller is designed to minimize the l1 
%       norm of the performance bound L1MAX. If L1MAX == Inf, constraint is
%       ignored.  
%       H2MAX: if this name-value parameter is specified, a sub-optimal 
%       H2 controller is designed to the value given. If unspecified or 
%       negative, an optimal H2 controller is designed to minimize the H2 
%       norm of the performance bound H2MAX. If H2MAX == Inf, constraint is
%       ignored.
%       HINFMAX: if this name-value parameter is specified, a sub-optimal 
%       HInf controller is designed to the value given. If unspecified or 
%       negative, an optimal HInf controller is designed to minimize the 
%       HInf norm of the performance bound HINFMAX. If HINFMAX == Inf, 
%       constraint is ignored.
%       DISPLAY: 'on' | 'off' (default), name-value parameter to specify if
%       the LMI solver should output progress information.
%       
%   Outputs:
%       K: H-infinity controller state-space model.
%       CL: closed loop system, equivalent to lft(P, K).
%       normz: H-infinity bound from w to z, minimized if GAM input was
%       unspecified.
%       xfeas: Sedumi optimizer structure output.
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
%   $Author: BH $    $Date: 2017-12-20 $    $Revision: 2 $
%
% REVISION 1:
%   2017-12-22 by BH: Changed YALMIP LMI format from that using BLKVAR
%   function to traditional MATLAB syntax (BLKVAR does not work with empty
%   matrix blocks and is not recommended by YALMIP documentation).
%
% REVISION 2:
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
    addOptional(p, 'W1', 1, @(x) isnumeric(x));
    addOptional(p, 'W2', 2, @(x) isnumeric(x));
    addOptional(p, 'WINF', 3, @(x) isnumeric(x));
    addOptional(p, 'Z1', 1, @(x) isnumeric(x));
    addOptional(p, 'Z2', 2, @(x) isnumeric(x));
    addOptional(p, 'ZINF', 3, @(x) isnumeric(x));
    addParameter(p, 'L1MAX', -1, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'H2MAX', -1, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'HINFMAX', -1, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'DISPLAY', 'off', @(x) ismember(x, {'on', 'off'}));
    addParameter(p, 'EXTENDED', true, @(x) isscalar(x) && islogical(x));
    parse(p, varargin{:});
    n = length(G.A);
    ny = p.Results.NMEAS; 
    nu = p.Results.NCON;
    psi_in = p.Results.L1MAX;
    gam_in = p.Results.HINFMAX;
    nu_in = p.Results.H2MAX;
    traceOutput = strcmp(p.Results.DISPLAY, 'on');
    isExtended = p.Results.EXTENDED;

    if G.Ts == 0
        error('my_syn_dt:ctplant', 'Continuous-time plants are not supported.');
    end
    
    [A, B, C, D] = partitionplant(G, ny, nu, p.Results.W1, p.Results.W2, p.Results.WINF, p.Results.Z1, p.Results.Z2, p.Results.ZINF);
    nw_1 = size(B{1}, 2); nw_2 = size(B{2}, 2); nw_inf = size(B{3}, 2);
    nz_1 = size(C{1}, 1); nz_2 = size(C{2}, 1); nz_inf = size(C{3}, 1);
    if nw_1==0 || nz_1==0
        psi_in = Inf;
    end
    if nw_2==0 || nz_2==0
        nu_in = Inf;
    end
    if nw_inf==0 || nz_inf==0
        gam_in = Inf;
    end
    
    %% Check Assumptions.
    P_2 = ss(A, B{4}, C{4}, D{4,4}, G.Ts); % State-space model for testing stability.
    [~, P_us] = stabsep(P_2);
    if rank(ctrb(P_us.A, P_us.B)) < size(P_us.A, 1) % Unstabilizable.
        error('my_syn_dt:unstab', 'Stabilizability assumption is violated.');
    end
    if rank(obsv(P_us.A, P_us.C)) < size(P_us.A, 2) % Undetectable.
        error('my_syn_dt:undet', 'Detectability assumption is violated.');
    end
    if all(D{4,4}==0)
        isFeedThrough = false;
    else
        isFeedThrough = true;
    end
    
    %% Define optimization variables.
    if psi_in == -1
        psi_1 = sdpvar(1); % Optimization problem.
    else
        psi_1 = psi_in.^2; % Sub-optimal (feasibility) problem.
    end
    if nu_in == -1
        nu_2 = sdpvar(1); % Optimization problem.
    else
        nu_2 = nu_in^2; % Sub-optimal (feasibility) problem.
    end
    if gam_in == -1
        gam_inf = sdpvar(1); % Optimization problem.
    else
        gam_inf = gam_in^2; % Sub-optimal (feasibility) problem.
    end
    
    function [linObjective, xfeas] = ymOptimize(linValue, isVerbose)

        % Scalars:
        alpha = linValue;
        mu_1 = sdpvar(1);
        nu_1 = sdpvar(1);
        % Matrices:
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
        R = sdpvar(nu, ny, 'full'); % Setting R as an equality constraint seems to make H-2 optimization problem infeasible.
%         if ~isinf(nu_in)
%             R = -D{2,4}\D{2,2}/D{4,2};
%         else
%             R = sdpvar(nu, ny, 'full');
%         end
        
        ymLmi = {};
        ymConstraint = [];%[norm(X)<=1e3, norm(Y)<=1e3]; % Added to fix numerical problems.
        
        %% Form l1 output feedback LMIs
        % Comments use notation from [2], for l1 output feedback:
        %   B{1} = B_w, B{4} = B_u, C{1} = C_z, C{4} = C_y, D{1,1} = D_zw, 
        %   D{1,4} = D_zu, D{4,1} = D_yw, D{4,4} = D_yu.
        if ~isinf(psi_in)
            % LMI 1.1:
            ymLmi{1} = [alpha*(X + X' - P)...	% Block (1,1): alpha*(X + X' - P)
                alpha*(eye(n) + S' - J)...      % Block (1,2): alpha*(I_n + S' - J)
                zeros(n, nw_1)...               % Block (1,3): 0_n,nw
                X'*C{1}' + L'*D{1,4}';...       % Block (1,4): X'*C_z' + L'*D_zu'
                alpha*(eye(n) + S - J')...      % Block (2,1): (Block (1,2))'
                alpha*(Y + Y' - H)...           % Block (2,2): alpha*(Y + Y' - H)
                zeros(n, nw_1)...               % Block (2,3): 0_n,nw
                C{1}' + C{4}'*R'*D{1,4}';...    % Block (2,4): C_z' + C_y'*R'*D_zu'
                zeros(nw_1, n)...               % Block (3,1): (Block (1,3))'
                zeros(nw_1, n)...               % Block (3,2): (Block (2,3))'
                (mu_1 - 1)*eye(nw_1)...         % Block (3,3): (mu - 1)*I_nw
                D{1,1}' + D{4,1}'*R'*D{1,4}';...% Block (3,4): D_zw' + D_yw'*R'*D_zu'
                C{1}*X + D{1,4}*L...            % Block (4,1): (Block (4,1))'
                C{1} + D{1,4}*R*C{4}...         % Block (4,2): (Block (2,4))'
                D{1,1} + D{1,4}*R*D{4,1}...     % Block (4,3): (Block (3,4))'
                nu_1*eye(nz_1)];                % Block (4,4): nu*I_nz
            ymConstraint = [ymConstraint, ymLmi{1}>=0];
            % LMI 1.2:
            ymLmi{2} = [P...                    % Block (1,1): P
                J...                            % Block (1,2): J
                A*X + B{4}*L...                 % Block (1,3): A*X + B_u*L
                A + B{4}*R*C{4}...              % Block (1,4): A + B_u*R*C_y
                B{1} + B{4}*R*D{4,1};...        % Block (1,5): B_w + B_u*R*D_yw
                J'...                           % Block (2,1): (Block (1,2))'
                H...                            % Block (2,2): H
                Q...                            % Block (2,3): Q
                Y*A + F*C{4}...                 % Block (2,4): Y*A + F*C_y
                Y*B{1} + F*D{4,1};...           % Block (2,5): Y*B_w + F*D_yw
                X'*A' + L'*B{4}'...             % Block (3,1): (Block (1,3))'
                Q'...                           % Block (3,2): (Block (2,3))'
                (1 - alpha)*(X + X' - P)...     % Block (3,3): (1 - alpha)*(X + X' - P)
                (1 - alpha)*(eye(n) + S' - J)...% Block (3,4): (1 - alpha)*(eye(n) + S' - J)
                zeros(n, nw_1);...              % Block (3,5): 0_n,nw
                A' + C{4}'*R'*B{4}'...          % Block (4,1): (Block (1,4))'
                A'*Y' + C{4}'*F'...             % Block (4,2): (Block (2,4))'
                (1 - alpha)*(eye(n) + S - J')...% Block (4,3): (Block (3,4))'
                (1 - alpha)*(Y + Y' - H)...     % Block (4,4): (1 - alpha)*(Y + Y' - H)
                zeros(n, nw_1);...              % Block (4,5): 0_n,nw
                B{1}' + D{4,1}'*R'*B{4}'...     % Block (5,1): (Block (1,5))'
                B{1}'*Y' + D{4,1}'*F'...        % Block (5,2): (Block (2,5))'
                zeros(nw_1, n)...               % Block (5,3): (Block (3,5))'
                zeros(nw_1, n)...               % Block (5,4): (Block (4,5))'
                eye(nw_1)];                     % Block (5,5): I_nw
            ymConstraint = [ymConstraint, ymLmi{2}>=0];
            % LMI 1.3:
            ymLmi{3} = [psi_1 mu_1 nu_1; mu_1 1 0; nu_1 0 1];
            ymConstraint = [ymConstraint, ymLmi{3}>=0];
        end
        
        %% Form H2 output feedback LMIs.
        % Comments use notation from [2], for H2 output feedback:
        %   B{2} = B_w, B{4} = B_u, C{2} = C_z, C{4} = C_y, D{2,2} = D_zw, 
        %   D{2,4} = D_zu, D{4,2} = D_yw, D{4,4} = D_yu.
        if ~isinf(nu_in)
            % LMI 2.1:
            ymLmi{4} = trace(W);
            ymConstraint = [ymConstraint, ymLmi{4} <= nu_2];
            % LMI 2.2:
            ymLmi{5} = [W...                % Block (1,1): W
                C{2}*X + D{2,4}*L...        % Block (1,2): % C_z*X + D_zu*L
                C{2} + D{2,4}*R*C{4};...	% Block (1,3): C_z + D_zu*R*C_y
                X'*C{2}' + L'*D{2,4}'...	% Block (2,1): (Block (1,2))'
                X + X' - P...               % Block (2,2): X + X' - P
                eye(n) + S' - J;...         % Block (2,3): I_n + S' - J
                C{2}' + C{4}'*R'*D{2,4}'...	% Block (3,1): (Block (1,3))'
                eye(n) + S - J'...          % Block (3,2): (Block (2,3))'
                Y + Y' - H];                % Block (3,3): Y + Y' - H
            ymConstraint = [ymConstraint, ymLmi{5} >= 0];
            % LMI 2.3:
            ymLmi{6} = [P...                % Block (1,1): P
                J...                        % Block (1,2): J
                A*X + B{4}*L...             % Block (1,3): A*X + B_u*L
                A + B{4}*R*C{4}...          % Block (1,4): A + B_u*R*C_y
                B{2} + B{4}*R*D{4,2};...    % Block (1,5): B_w + B_u*R*D_yw
                J'...                       % Block (2,1): (Block (1,2))'
                H...                        % Block (2,2): H
                Q...                        % Block (2,3): Q
                Y*A + F*C{4}...             % Block (2,4): Y*A + F*C_y
                Y*B{2} + F*D{4,2};...       % Block (2,5): Y*B_w + F*D_yw
                X'*A' + L'*B{4}'...         % Block (3,1): (Block (1,3))'
                Q'...                       % Block (3,2): (Block (2,3))'
                X + X' - P...               % Block (3,3): X + X' - P
                eye(n) + S' - J...          % Block (3,4): I_n + S' - J
                zeros(n, nw_2);...          % Block (3,5): 0_n,nw
                A' + C{4}'*R'*B{4}'...      % Block (4,1): (Block (1,4))'
                A'*Y' + C{4}'*F'...         % Block (4,2): (Block (2,4))'
                eye(n) + S - J'...          % Block (4,3): (Block (3,4))'
                Y + Y' - H...               % Block (4,4): Y + Y' - H
                zeros(n, nw_2);...          % Block (4,5): 0_n,nw
                B{2}' + D{4,2}'*R'*B{4}'... % Block (5,1): (Block (1,5))'
                B{2}'*Y' + D{4,2}'*F'...    % Block (5,2): (Block (2,5))'
                zeros(nw_2, n)...           % Block (5,3): (Block (3,5))'
                zeros(nw_2, n)...           % Block (5,4): (Block (4,5))'
                eye(nw_2)];                 % Block (5,5): I_nw_2
            ymConstraint = [ymConstraint, ymLmi{6} >= 0];
            % The below equality contraint (44) in [2] seems to make the
            % H-2 optimization problem infeasible. Instead, solve directly
            % for R if H-2 optimization is to be done.
            % LMI 2.4:
            %ymLmi{7} = D{2,2} + D{2,4}*R*D{4,2}; % D_zw + D_zu*R*D_yw
            %ymConstraint = [ymConstraint, ymLmi{7} == zeros(nz_2, nw_2)];
            % 2018-01-24: LMI 2.4 should not be required because
            % discrete-time biproper systems have a valid H2 norm.
        end
        
        %% Form HInf output feedback LMIs.
        % Comments use notation from [2], for HInf output feedback:
        %   B{3} = B_w, B{4} = B_u, C{3} = C_z, C{4} = C_y, D{3,3} = D_zw, 
        %   D{3,4} = D_zu, D{4,3} = D_yw, D{4,4} = D_yu.
        if ~isinf(gam_in)
            % LMI 3.1:
            ymLmi{8} = [P...                    % Block (1,1): P
                J...                            % Block (1,2): J
                A*X + B{4}*L...                 % Block (1,3): A*X + B_u*L
                A + B{4}*R*C{4}...              % Block (1,4): A + B_u*R*C_y
                B{3} + B{4}*R*D{4,3}...         % Block (1,5): B_w + B_u*R*D_yw
                zeros(n, nz_inf);...            % Block (1,6): 0_n,nz
                J'...                           % Block (2,1): (Block (1,2))'
                H...                            % Block (2,2): H
                Q...                            % Block (2,3): Q
                Y*A + F*C{4}...                 % Block (2,4): Y*A + F*C_y
                Y*B{3} + F*D{4,3}...            % Block (2,5): Y*B_w + F*D_yw
                zeros(n, nz_inf);...            % Block (2,6): 0_n,nz
                X'*A' + L'*B{4}'...             % Block (3,1): (Block (1,3))'
                Q'...                           % Block (3,2): (Block (2,3))'
                X + X' - P...                   % Block (3,3): X + X' - P
                eye(n) + S' - J...              % Block (3,4): I_n + S' - J
                zeros(n, nw_inf)...             % Block (3,5): 0_n,nw
                X'*C{3}' + L'*D{3,4}';...       % Block (3,6): X'*C_z' + L'*D_zu'
                A' + C{4}'*R'*B{4}'...          % Block (4,1): (Block (1,4))'
                A'*Y' + C{4}'*F'...             % Block (4,2): (Block (2,4))'
                eye(n) + S - J'...              % Block (4,3): (Block (3,4))'
                Y + Y' - H...                   % Block (4,4): Y + Y' - H
                zeros(n, nw_inf)...             % Block (4,5): 0_n,nw
                C{3}' + C{4}'*R'*D{3,4}';...    % Block (4,6): C_z' + C_y'*R'*D_zu'
                B{3}' + D{4,3}'*R'*B{4}'...     % Block (5,1): (Block (1,5))'
                B{3}'*Y' + D{4,3}'*F'...        % Block (5,2): (Block (2,5))'
                zeros(nw_inf, n)...             % Block (5,3): (Block (3,5))'
                zeros(nw_inf, n)...             % Block (5,4): (Block (4,5))'
                eye(nw_inf)...                  % Block (5,5): I_nw_inf
                D{3,3}' + D{4,3}'*R'*D{3,4}';...% Block (5,6): D_zw' + D_yw'*R'*D_zu'
                zeros(nz_inf, n)...             % Block (6,1): (Block (1,6))'
                zeros(nz_inf, n)...             % Block (6,2): (Block (2,6))'
                C{3}*X + D{3,4}*L...            % Block (6,3): (Block (3,6))'
                C{3} + D{3,4}*R*C{4}...         % Block (6,4): (Block (4,6))'
                D{3,3} + D{3,4}*R*D{4,3}...     % Block (6,5): (Block (5,6))'
                gam_inf*eye(nz_inf)];           % Block (6,6): mu*I_nz_inf
            ymConstraint = [ymConstraint, ymLmi{8} >= 0];
        end

        %% Solve LMI Minimization Problem.
        if psi_in==-1 || nu_in==-1 || gam_in==-1 % Optimization problem.
            ymObjective = 0;
            if psi_in==-1 % l1 optimization.
                ymObjective = ymObjective + psi_1;
            end
            if nu_in==-1 % H2 optimization.
                ymObjective = ymObjective + nu_2;
            end
            if gam_in==-1 % HInf optimization.
                ymObjective = ymObjective + gam_inf;
            end
        else % Feasibility problem.
            ymObjective = [];
        end
        ymOptions = sdpsettings('verbose', isVerbose, 'solver', 'sdpt3', 'sdpt3.inftol', 1e-9, 'sdpt3.gaptol', 1e-9);
        xfeas = optimize(ymConstraint, ymObjective, ymOptions);
        linObjective = value(psi_1);
    end

    %% Solve optimization/feasibility problem.
    % Line search over alpha (if l1 optimization).
    if psi_in==-1
        if traceOutput
            linOptions = optimset('Display', 'iter', 'TolX', 1e-2);
        else
            linOptions = optimset('Display', 'off', 'TolX', 1e-2);
        end
        linValue = fminbnd(@(x) ymOptimize(x, false), 0, 1, linOptions)
    else
        linValue = 0;
    end
    % Repeat last iteration to produce output.
    [~, ymDiagnostics] = ymOptimize(linValue, traceOutput);
    if ymDiagnostics.problem==1
        error('my_syn_dt:solver_infeasible', ymDiagnostics.info);
    elseif ymDiagnostics.problem~=0
        warning(ymDiagnostics.info);
    end
    
    %% Synthesize controller.
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
    Theta = [inv(N) -N\Y*B{4}; zeros(nu, n) eye(nu)]*...
        [Q - Y*A*X F; L R]*[inv(M) zeros(n, ny); -C{4}*X/M eye(ny)];
    
    %% Reduce Nonzero D_22 Case (Optional).
    if isFeedThrough
        % Below variables have hat in notation from [2].
        A_K = Theta(1:n, 1:n);
        B_K = Theta(1:n, n+1:n+ny);
        C_K = Theta(n+1:n+nu, 1:n);
        D_K = Theta(n+1:n+nu, n+1:n+ny);
        Theta = [A_K B_K; zeros(size(Theta, 1)-n, size(Theta, 2))]...
            + [-B_K*D{4,4}; eye(size(B_K*D{4,4}, 2))]...
            /(eye(size(B_K*D{4,4}, 2)) + D_K*D{4,4})...
            *[C_K D_K];
    end
    K = ss(Theta(1:n, 1:n), Theta(1:n, n+1:n+ny), Theta(n+1:n+nu, 1:n), Theta(n+1:n+nu, n+1:n+ny), G.Ts);
    
    normz = sqrt([value(psi_1) value(nu_2) value(gam_inf)]);
    
    CL = lft(G, K); % Close the loop (emulate output from H2HINFSYN.m).
end

function [A, B, C, D] = partitionplant(P, ny, nu, w_1, w_2, w_inf, z_1, z_2, z_inf)
%PARTITIONPLANT Modify plant for multiobjective synthesis.
%   Permutes the inputs and outputs of the plant P so that each
%   input/output group is separate for specifying performance criteria.
%   
%   Inputs:
%       P: state-space generalized plant.
%       ny: numeric scalar number of output channels (controller 
%       measurement channels). Outputs are the last ny output channels.
%       nu: numeric scalar number of input channels (controller controlled 
%       channels). Inputs are the last nu input channels.
%       w_1: numeric vector of disturbance input channel indices under l1 
%       norm constraint.
%       w_2: numeric vector of disturbance input channel indices under H2 
%       norm constraint.
%       w_inf: numeric vector of disturbance input channel indices under 
%       HInf norm constraint.
%       z_1: numeric vector of performance output channel indices under l1 
%       norm constraint.
%       z_2: numeric vector of performance output channel indices under H2 
%       norm constraint.
%       z_inf: numeric vector of performance output channel indices under 
%       HInf norm constraint.
%       
%   Outputs:
%       A: state-space A matrix (equal to P.A).
%       B: 1x4 cell array of partitioned state-space B matrices.
%       C: 4x1 cell array of partitioned state-space C matrices.
%       D: 4x4 cell array of partitioned state-space D matrices.
%
%   $Author: BH $    $Date: 2017-12-20 $    $Revision: 1 $

    %% Plant partitioning:
    %       In:    x     w_1    w_2   w_inf    u
    %              |      |      |      |      |
    % Out:         V      V      V      V      V
    %          .----------------------------------.
    % dx/dt <- |   A  |  B_1 |  B_2 |  B_3 |  B_4 | n
    %          |------|------|------|------|------|
    %   z_1 <- |  C_1 | D_11 |(D_12)|(D_13)| D_14 | nz_1
    %          |------|------|------|------|------|
    %   z_2 <- |  C_2 |(D_21)| D_22 |(D_23)| D_24 | nz_2
    %          |------|------|------|------|------|
    % z_inf <- |  C_3 |(D_31)|(D_32)| D_33 | D_34 | nz_inf
    %          |------|------|------|------|------|
    %     y <- |  C_4 | D_41 | D_42 | D_43 | D_44 | ny
    %          '----------------------------------'
    %      Size:   n    nw_1   nw_2  nw_inf   nu
    
    %% Permute rows and columns.
    ABCD = [P.A P.B; P.C P.D]; % Doyle representation.
    nx = length(P.A);
    ABCD = ABCD(:,[1:nx nx+w_1 nx+w_2 nx+w_inf end-nu+1:end]);
    ABCD = ABCD([1:nx nx+z_1 nx+z_2 nx+z_inf end-ny+1:end],:);
    
    %% Number of channels for each input/output group.
    nw_1 = numel(w_1); nw_2 = numel(w_2); nw_inf = numel(w_inf);
    nz_1 = numel(z_1); nz_2 = numel(z_2); nz_inf = numel(z_inf);
    
    %% Split into cell array.
    A = P.A;
    B = mat2cell(ABCD(1:nx,nx+1:end), nx, [nw_1 nw_2 nw_inf nu]);
    C = mat2cell(ABCD(nx+1:end,1:nx), [nz_1 nz_2 nz_inf ny], nx);
    D = mat2cell(ABCD(nx+1:end,nx+1:end), [nz_1 nz_2 nz_inf ny], [nw_1 nw_2 nw_inf nu]);
end