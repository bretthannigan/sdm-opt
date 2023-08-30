function [H, S, CL, normz, diagn, iterProgress] = sdsyn(varargin) 
%SDSYN Sigma Delta Loop Filter Synthesis
%   Design a sigma delta modulator loop filter according to performance and
%   stability criteria defined by H-Infinity (GKYP), H-Infinity, H-2, and 
%   l-1 system norms. Optionally permits an uncertain quantizer gain in the
%   loop around which stability can be enforced.
%
%   Inputs:
%       H_0: lti loop filter prototype model initial guess which also 
%           defines the system order and sample time (required).
%       synOpt: strucure of optimization targets defined with 
%           DEFINESYNOPT.M (required).
%       UncertainGain: scalar ureal of multiplicative uncertainty modelling
%           the quantizer gain OR scalar fixed gain OR 2-by-2 matrix of 
%           upper LFT transform of the uncertain gain 
%           (optional, default: 1).
%       Display: 'on' | 'off' (default) name-value string parameter to
%           specify if the LMI solver should output progress information
%           (optional).
%
%   Outputs:
%       H: loop filter, in format specified by input H_0 OR numeric scalar
%           indicating the order of the filter to be synthesized (default:
%           discrete-time system with H zeros and H poles at the origin).
%       S: sensitivity function (NTF), in format specified by input H_0.
%       CL: closed-loop system with uncertainty block (if present).
%       normz: achieved system norms.
%       feas: optimizer strucure output.
%
%       See also: DEFINESYNOPT, DTSYN.
%           [1] T. Iwasaki and S. Hara, “MATHEMATICAL ENGINEERING 
%               Generalized KYP Lemma : Unified Characterization of 
%               Frequency Domain Inequalities with Applications to System 
%               Design,” no. August, 2003.
%           [2] X. Li, C. Yu, and H. Gao, “Design of delta-sigma modulators 
%               via generalized Kalman-Yakubovich-Popov lemma,” Automatica, 
%               vol. 50, no. 10, pp. 2700–2708, 2014.
%           [3] S. L. Shishkin, “Optimization under non-convex Quadratic 
%               Matrix Inequality constraints with application to design of
%               optimal sparse controller,” IFAC-PapersOnLine, vol. 50, no.
%               1, pp. 10754–10759, 2017.
%           [4] I. Masubuchi, A. Ohara, and N. Suda, “LMI-based controller 
%               synthesis: a unified formulation and solution,” Robust 
%               Nonlinear Control, vol. 8, no. 9, pp. 669–686, 1998.
%           [5] A. Oberoi and J. C. Cockburn, “A simplified LMI approach to 
%               l1 Controller Design,” Proc. 2005 Am. Control Conf., no. 
%               July, pp. 1788–1792, 2005.
%           [6] J. Bu and M. Sznaier, “Linear matrix inequality approach to
%               synthesizing low-order suboptimal mixed l1/Hp controllers,”
%               Automatica, vol. 36, no. 7, pp. 957–963, 2000.
%           [7] J. Löfberg, “YALMIP: A Toolbox for Modeling and 
%               Optimization in MATLAB,” In Proceedings of the CACSD 
%               Conference. Taipei, Taiwan, 2004.
%           [8] S. Boyd and L. Vandenberghe, “Semidefinite Programming 
%               Relaxations of Non-Convex Problems in Control and 
%               Combinatorial Optimization,” Commun. Comput. Control Signal
%               Process. A Tribut. to Thomas Kailath, pp. 279–288, 1997.
%
%   $Author: BH $    $Date: 2018-08-15 $    $Revision: 2 $
%
% REVISION 1:
%   2018-10-17 by BH: Simplified the derivation of b_s.
%
% REVISION 2:
%   2018-11-16 by BH: Added functionality to make initial guess using
%   Shor's relaxation.

    %% Deal with variable function inputs.
    p = inputParser;
    addRequired(p, 'H_0', @(x) (isnumeric(x) && isscalar(x)) || isa(x, 'lti') || isempty(x));
    addRequired(p, 'SynOpt', @isstruct);
    addOptional(p, 'UncertainGain', 1, @(x) isscalar(x) || (isnumeric(x) && all(size(x)==[2 2])));
    addParameter(p, 'Display', 'off', @(x) ismember(x, {'on', 'off'}));
    addParameter(p, 'Solver', 'lmilab', @(x) ischar(x)); % LMILAB is by far the slowest solver but seems to give the best results.
    addParameter(p, 'TermParam', struct('maxIter', 500, 'kappa', 0, 'epsilon', 1e-6), @(x) isstruct(x) && isfield(x, 'maxIter') && isfield(x, 'kappa') && isfield(x, 'epsilon'));
    addParameter(p, 'InitialGuess', false, @(x) isscalar(x) && islogical(x));
    parse(p, varargin{:});
    if isa(p.Results.H_0, 'lti')
        H_0 = p.Results.H_0;
    else
        H_0 = zpk(zeros(1, p.Results.H_0), zeros(1, p.Results.H_0), 1, 1); % Default n-th order DT system.
    end
    synOpt = p.Results.SynOpt;
    K = p.Results.UncertainGain;
    isVerbose = strcmp(p.Results.Display, 'on');
    ymSolver = p.Results.Solver;
    termParam = p.Results.TermParam;
    isInitialGuess = p.Results.InitialGuess;
    ymCellValue = @(x) cellfun(@value, x, 'UniformOutput', false);
    
    %% Take LFT to separate the uncertain quantizer gain.
    if isa(K, 'ureal')
        [M, delta] = lftdata(K); % delta is the norm-bounded uncertainty.
    elseif isnumeric(K) && isscalar(K)
        M = [0 0; 0 K];
        delta = 0; % There is no uncertainty.
    elseif isnumeric(K) && all(size(K) == [2, 2])
        M = K;
        delta = 0;
    else
        error('sdsyn:invalidUncertainGain', 'Uncertain gain must be a numeric scalar, ureal scalar, or 2-by-2 matrix.');
    end
    
    n = order(H_0);
    Ts = H_0.Ts;
    
    b = sdpvar(n, 1, 'full'); % Loop filter numerator coefficients.
    b_op = sdpvar(n, 1, 'full');
    a = sdpvar(n, 1, 'full'); % Loop filter denominator coefficients.
    a_op = sdpvar(n, 1, 'full');
    
    P = sdpvar(n*ones(1, length(synOpt)), n*ones(1, length(synOpt)), 'symmetric');
    Q = sdpvar(n*ones(1, length(synOpt)), n*ones(1, length(synOpt)), 'symmetric');
    
    alpha = sdpvar(ones(1, length(synOpt)), ones(1, length(synOpt))); % l-1/(*) norm bilinear matrix inequality parameter.
    gam = sdpvar(ones(1, length(synOpt)), ones(1, length(synOpt))); % Optimization problem.
    
    if length(synOpt)==1
        P = {P};
        Q = {Q};
        alpha = {alpha};
        gam = {gam};
    end
    
    if isInitialGuess
        aa = sdpvar(n, n);
    else
        aa = zeros(n, n);
    end
    
    %% Partition plant.
    [synOpt, CL_in] = sdsyn_buildplant(synOpt, b, a, M);
    
    %% Loop through constraints and build LMIs.
    ymConstraint = [];
    for iConstraint=1:length(synOpt)
        
        %% Load in generalized plant matrices.
        
        A = synOpt(iConstraint).A;
        B = synOpt(iConstraint).B;
        C = synOpt(iConstraint).C;
        D = synOpt(iConstraint).D; % Closed-loop feedthrough matrix.
        %a_s = -synOpt(iConstraint).A(end, :)'; % Sensitivity function denominator coefficients (equivalent to vector d from [2]).
        %b_s = synOpt(iConstraint).C' + a_s; % Sensitivity function numerator coefficients (equivalent to vector c from [2]).
        if ~all(abs(synOpt(iConstraint).A(1:n-1,1:n)-[zeros(n-1, 1) eye(n-1)])<1e-6)
        	error('sdsyn:AMatrix', ['Invalid A matrix, constraint ' num2str(iConstraint) ' is not in controllable canonical form.']);
        end
        if all(abs(synOpt(iConstraint).B(1:n-1))<1e-6)
            T = eye(n)./synOpt(iConstraint).B(n);
        else
            error('sdsyn:BMatrix', ['Invalid B matrix, constraint ' num2str(iConstraint) ' is not in controllable canonical form.']);
        end
        
        %% Transform to obtain CCF.
        % SDSYN_BUILDPLANT does not always produce complete controllable
        % canonical form representations. Often, the B matrix has one
        % non-unity term.
        
        b_s = T'\C'; % Previously: b_s = T'\C' - (T'\A'*T')*T*B*D';
        b_s_op = replace(b_s, [b a], [b_op a_op]);
        a_s = -(T'\A'*T')*T*B;
        a_s_op = replace(a_s, [b a], [b_op a_op]);
        B = T*B; % Closed-loop input-state matrix.
        C = C/T; % Closed-loop state-output matrix.
        C_op = replace(C, [b a], [b_op a_op]);
        
        %% Define optimization or feasibility targets.
        
        if isinf(synOpt(iConstraint).constraint)
            break
        elseif synOpt(iConstraint).constraint~=-1
            gam{iConstraint} = synOpt(iConstraint).constraint.^2; % Sub-optimal (feasibility) problem.
        else
            gam{iConstraint} = sdpvar(1); % Optimization problem.
        end
        
        switch synOpt(iConstraint).norm
            case 1
                %% Form l-1/(*) output feedback LMI.
                partialConstraint = sdsyn_l1_lmi(gam{iConstraint}, P{iConstraint}, alpha{iConstraint}, synOpt(iConstraint).ffi, Ts);
            case 2
                %% Form H-2 output feedback LMI.
                partialConstraint = sdsyn_h2_lmi(gam{iConstraint}, P{iConstraint}, synOpt(iConstraint).ffi, Ts);
                alpha{iConstraint} = 0; % Alpha only used for l-1/(*) norm LMI.
            case Inf
                %% Form H-Inf output feedback LMI.
                partialConstraint = sdsyn_gkyp_lmi(gam{iConstraint}, P{iConstraint}, Q{iConstraint}, synOpt(iConstraint).ffi, Ts);
                alpha{iConstraint} = 0; % Alpha only used for l-1/(*) norm LMI.
        end
        ymConstraint = [ymConstraint, partialConstraint]; %#ok<AGROW>
    end

    ymObjective = sum([gam{:}]); 
    ymOptions = sdpsettings('verbose', isVerbose, 'solver', ymSolver, 'lmilab.reltol', 1e-9, 'lmilab.maxiter', 500, 'lmilab.feasradius', 1e9, 'lmilab.L', 100); % Settings from [2].
    if isInitialGuess % Add the constraint and objective term for Shor relaxation [8].
        ymConstraint = [ymConstraint, [aa a_op; a_op' 1]>=0];
        ymObjective = ymObjective + trace(aa);
        [~, diagn] = sdsyn_initialguess(zeros(n, 1), zeros(n, 1));
        [H, S, CL] = sdsyn_tolti(value(b), value(a), CL_in, delta);
        normz = sqrt(cell2mat(ymCellValue(gam)));
        iterProgress = struct('a', value(a), 'b', value(b), 'gam', cell2mat(ymCellValue(gam)));
    else
        if ~isproper(H_0)
            error('sdsyn_initialguess:improperFilter', 'The initial guess system must be proper.');
        end
        if ~isstable(H_0)
            error('sdsyn_initialguess:unstableFilter', 'The initial guess system must be stable.');
        end
        H_0_tf = tf(H_0);
        b_0 = fliplr(H_0_tf.num{1}(2:end))';
        a_0 = fliplr(H_0_tf.den{1}(2:end))';
        sdsyn_initialguess(b_0, a_0);
        [iterProgress, diagn] = sdsyn_nlinhandler(termParam.maxIter, termParam.kappa, termParam.epsilon);
        [H, S, CL] = sdsyn_tolti(iterProgress(end).b, iterProgress(end).a, CL_in, delta);
        normz = sqrt(iterProgress(end).gam);
    end
    
    function ymConstraint = sdsyn_l1_lmi(gam, P, alpha, ffi, Ts)
    %FILTERSYNTHESIS_L1_LMI l-1 Linear Matrix Inequality
    %   Establishes the linear matrix inequality for the l-1/(*) norm
    %   optimization or feasibility program. Matrix variables that are 
    %   shared among LMIs are passed globally.
    %
    %   Input:
    %       gam: the optimized variable (constaint if sub-optimal
    %       feasibility test is to be done).
    %   Outputs:
    %       ymConstraint: a YALMIP constraint object containing the LMIs.
    %
    %   See also: SDSYN
    %       [1] T. Iwasaki and S. Hara, “MATHEMATICAL ENGINEERING 
    %           Generalized KYP Lemma : Unified Characterization of 
    %           Frequency Domain Inequalities with Applications to System 
    %           Design,” no. August, 2003.
    %       [2] A. Oberoi and J. C. Cockburn, “A simplified LMI approach to 
    %           l1 Controller Design,” Proc. 2005 Am. Control Conf., no. 
    %           July, pp. 1788–1792, 2005.
    %       [3] J. Bu and M. Sznaier, “Linear matrix inequality approach to
    %           synthesizing low-order suboptimal mixed l1/Hp controllers,”
    %           Automatica, vol. 36, no. 7, pp. 957–963, 2000.
    %
    %   $Author: BH $    $Date: 2018-08-27 $    $Revision: 0 $
    
        %% Define local l-1/(*) optimization variables.
        % These are local to the l-1/(*) LMIs and not shared among 
        % constraints.
        
        if Ts==0 % Continuous-time.
            Phi = [0 1; 1 0];
            if ~all(ffi==[0 Inf])
                warning('Finite frequency interval ignored for l-1/(*) constraint.')
            end
        else % Discrete-time.
            Phi = [1 0; 0 -1];
            if ~all(ffi==[0 pi])
                warning('Finite frequency interval ignored for l-1/(*) constraint.')
            end
        end
        
        A_c = [zeros(n-1, 1) eye(n-1); zeros(1, n)];
        outerFactor = [A_c B; eye(n) zeros(n, 1)];
        Theta = @(P, a) kron((Phi + [0 0; 0 a]), P);
        
        mu = sdpvar(1);
        nu = sdpvar(1);
        
        %% Form l-1/(*) output feedback LMI.
        
        ymLmi{1} = [alpha*P zeros(n, 1) (C + C_op)';...
            zeros(1, n) (mu - 1) D';...
            (C + C_op) D nu];
        ymLmi{2} = -outerFactor'*Theta(P, alpha)*outerFactor...
            + [a_s*a_s_op' + a_s_op*a_s' + a_s_op*a_s_op' + aa (a_s + a_s_op); (a_s + a_s_op)' 1];
        ymLmi{3} = [gam mu nu; mu 1 0; nu 0 1];
        ymConstraint = [ymLmi{1}>=0, ymLmi{2}>=0, ymLmi{3}>=0, P>=0, nu>=0, mu>=0];
    end
    
    function ymConstraint = sdsyn_h2_lmi(gam, P, ffi, Ts)
    %FILTERSYNTHESIS_H2_LMI H-2 Linear Matrix Inequality
    %   Establishes the linear matrix inequality for the H-2 norm
    %   optimization or feasibility program. Matrix variables that are 
    %   shared among LMIs are passed globally.
    %
    %   Input:
    %       gam: the optimized variable (constaint if sub-optimal
    %       feasibility test is to be done).
    %   Outputs:
    %       ymConstraint: a YALMIP constraint object containing the LMIs.
    %
    %   See also: SDSYN
    %       [1] T. Iwasaki and S. Hara, “MATHEMATICAL ENGINEERING 
    %           Generalized KYP Lemma : Unified Characterization of 
    %           Frequency Domain Inequalities with Applications to System 
    %           Design,” no. August, 2003.
    %       [2] I. Masubuchi, A. Ohara, and N. Suda, “LMI-based controller 
    %           synthesis: a unified formulation and solution,” Robust 
    %           Nonlinear Control, vol. 8, no. 9, pp. 669–686, 1998.
    %
    %   $Author: BH $    $Date: 2018-08-27 $    $Revision: 0 $
    
        %% Define local H-2 optimization variables.
        % These are local to the H-2 LMIs and not shared among constraints.
        
        if Ts==0 % Continuous-time.
            Phi = [0 1; 1 0];
            if ~all(ffi==[0 Inf])
                warning('Finite frequency interval ignored for H-2 constraint.')
            end
        else % Discrete-time.
            Phi = [1 0; 0 -1];
            if ~all(ffi==[0 pi])
                warning('Finite frequency interval ignored for H-2 constraint.')
            end
        end
        
        A_c = [zeros(n-1, 1) eye(n-1); zeros(1, n)];
        outerFactor = [A_c B; eye(n) zeros(n, 1)];
        Theta = @(P) kron(Phi, P);
        
        %% Form H-2 output feedback LMI.
        % Notation is the same as [2], except gam = R.
        
        ymLmi{1} = [gam (C + C_op) D;...
            (C + C_op)' P zeros(n, 1);...
            D' zeros(1, n) 1];
        ymLmi{2} = -outerFactor'*Theta(P)*outerFactor...
            + [a_s*a_s_op' + a_s_op*a_s' + a_s_op*a_s_op' + aa (a_s + a_s_op); (a_s + a_s_op)' 1];
        ymConstraint = [ymLmi{1}>=0, ymLmi{2}>=0];
    end

    function ymConstraint = sdsyn_gkyp_lmi(gam, P, Q, ffi, Ts)
    %FILTERSYNTHESIS_GKYP_LMI H-Infinity Linear Matrix Inequality Using the
    % Generalized KYP Lemma
    %   Establishes the linear matrix inequality for the H-Inf norm
    %   optimization or feasibility program using the Generalized 
    %   Kalman-Yakubovich-Popov Lemma to restrict the LMI to a finite 
    %   frequency interval. Matrix variables that are shared among LMIs are
    %   passed globally.
    %
    %   Input:
    %       gam: the optimized variable (constaint if sub-optimal
    %       feasibility test is to be done).
    %       ffi: a 2-element vector defining the finite frequency interval 
    %           (for H-Inf GKYP condition), otherwise [0 2*pi];
    %   Outputs:
    %       ymConstraint: a YALMIP constraint object containing the LMIs.
    %
    %   See also: SDSYN
    %       [1] T. Iwasaki and S. Hara, “MATHEMATICAL ENGINEERING 
    %           Generalized KYP Lemma : Unified Characterization of 
    %           Frequency Domain Inequalities with Applications to System 
    %           Design,” no. August, 2003.
    %
    %   $Author: BH $    $Date: 2018-08-15 $    $Revision: 1 $
    %
    % REVISION 1:
    %   2018-09-12 by BH: Fixed CT frequency calculations. 
    
        %% Define local H-Inf optimization variables.
        % These are local to the H-Inf LMIs and not shared among 
        % constraints.
        
        % Frequency range
        % Reflect across real axis.
        if Ts==0 % Continuous-time.
            if isinf(ffi(2))
                w = [0 0];
            elseif ffi(1)==0
                w = [-ffi(2) ffi(2)];
            else
                w = [ffi(1) ffi(2)];
            end
            wc = (w(2) + w(1))/2;
            Phi = [0 1; 1 0];
            Psi = [-1 1i*wc; -1i*wc -w(1)*w(2)];
        else % Discrete-time.
            if ffi(1)==0 % Low frequency case, Psi_d = [-2*cos(w(1)) 1; 1 0]
                w = [-ffi(2) ffi(2)];
            elseif ffi(2)==pi % High  frequency case, Psi_d = [2*cos(w(2)) -1; -1 0]
                w = [ffi(1) 2*pi-ffi(1)];
            else
                w = [ffi(1) ffi(2)];
            end
            wc = (w(2) + w(1))/2;
            Phi = [1 0; 0 -1];
            w0 = (w(2) - w(1))/2;
            Psi = [0 exp(1i*wc); exp(-1i*wc) -2*cos(w0)];
        end
        
        A_c = [zeros(n-1, 1) eye(n-1); zeros(1, n)];
        outerFactor = [A_c B; eye(n) zeros(n, 1)];
        Theta = @(Q, P) kron(Phi, P) + kron(Psi, Q);
        
        if all(ffi==[0 pi]) || all(ffi==[0 Inf]) % All frequencies.
            Q = zeros(n);
            ymConstraint = P>=0;
        else
            ymConstraint = [];
        end
        
        %% Form H-Inf GKYP output feedback LMI.
        % Notation is the same as [1], except a_s = d, b_s = c, gam =
        % gam^2, A_c = A, A = A^, C = C^, and the presence of D.
%         ymLmi{1} = -[outerFactor'*Theta(Q, P)*outerFactor - [zeros(n) (a_s + a_s_op); (a_s + a_s_op)' 1]...
%             [b_s + b_s_op; D];...
%             (b_s + b_s_op)'...
%             D' ...
%             -gam] + [a_s*a_s_op' + a_s_op*a_s' + a_s_op*a_s_op' zeros(n, 2); zeros(2, n) zeros(2, 2)];
        ymLmi{1} = -[outerFactor'*Theta(Q, P)*outerFactor - [zeros(n) (a_s + a_s_op); (a_s + a_s_op)' 1]...
            [(C + C_op)' + (a_s + a_s_op)*D'; D'];...
            (C + C_op) + (a_s + a_s_op)'*D...
            D...
            -gam] + [a_s*a_s_op' + a_s_op*a_s' + a_s_op*a_s_op' + aa zeros(n, 2); zeros(2, n) zeros(2, 2)];
        ymConstraint = [ymConstraint, ymLmi{1}>=0, Q>=0, gam>=0];
    end

    function [ymConstraintIG, ymDiagnostics] = sdsyn_initialguess(b_0, a_0)
    %SDSYN_INITIALGUESS Generate Initial Guess for LMI Sovler
    %   Uses an initial guess to produce the auxiliary variables for the
    %   next steps.
    %
    %   Inputs:
    %       a_0: vector of loop filter denominator coefficients that are
    %           feasible for the initial guess LMI problem.
    %       b_0: vector of loop filter numerator coefficients that are
    %           feasible for the initial guess LMI problem.
    %
    %   Outputs:
    %       ymConstraintIG: a YALMIP constraint object containing the 
    %           initial guess LMIs.
    %       ymDiagnostics: YALMIP diagnostics structure.
    %
    %   See also: SDSYN
    %       [1] X. Li, C. Yu, and H. Gao, “Design of delta-sigma modulators 
    %               via generalized Kalman-Yakubovich-Popov lemma,” 
    %               Automatica, vol. 50, no. 10, pp. 2700–2708, 2014.
    %
    %   $Author: BH $    $Date: 2018-08-16 $    $Revision: 3 $
    %
    % REVISION 1:
    %   2018-08-21 by BH: Removed P_0, Q_0 outputs, added some error
    %   checking.
    %
    % REVISION 2:
    %   2018-10-25 by BH: Corrected docstrings only.
    %
    % REVISION 3:
    %   2018-11-16 by BH: Passing a_0, b_0 as vectors rather than a
    %   transfer function to allow zero denominator.

        %ymConstraintIG = replace(ymConstraint, [a b], [zeros(size(a)) zeros(size(b))]);
        %ymConstraintIG = replace(ymConstraintIG, [a_op b_op], [a_0 b_0]);
        %ymConstraintIG = replace(ymConstraint, a, zeros(size(a)));
        %ymConstraintIG = replace(ymConstraintIG, b, zeros(size(b)));
        ymConstraintIG = replace(ymConstraint, [a_op b_op], [a_0 b_0]);
        ymDiagnostics = sdsyn_blinhandler(ymConstraintIG, ymObjective, ymOptions);
        if ymDiagnostics.problem==1
            error('sdsyn_initialguess:infeasIG', 'Initial guess is infeasible.');
        elseif ymDiagnostics.problem~=0
            warning(ymDiagnostics.info);
        end
    end

    function [iterProgress, ymDiagnostics] = sdsyn_nlinhandler(maxIter, kappa, epsilon)
    %SDSYN_NLINHANDLER Handler for SDP Problems with Non-Convex QMI
    % Constraints
    %   Handler for an iterative method that allows solving of quadratic 
    %   matrix inequality constraints with linear objectives. Non-convex
    %   constrants are replaced by a sum of 2 terms, the first of which is
    %   the linearized/convexified operating point and the second of which
    %   is an incremental change (see [1] for background).
    %
    %   Inputs:
    %       maxIter: maximum number of interations.
    %       kappa: tuning parameter for termination.
    %       epsilon: tolerance for change in output.
    %
    %   Output:
    %       iterProgress: structure containing iteration by interation
    %           variables and objective values.
    %
    %   See also: SDSYN, SDSYN_BLINHANDLER
    %       [1] S. L. Shishkin, “Optimization under non-convex Quadratic 
    %           Matrix Inequality constraints with application to design of
    %           optimal sparse controller,” IFAC-PapersOnLine, vol. 50, no.
    %           1, pp. 10754–10759, 2017.
    %       [2] J. Löfberg, “YALMIP: A Toolbox for Modeling and 
    %           Optimization in MATLAB,” In Proceedings of the CACSD 
    %           Conference. Taipei, Taiwan, 2004.
    
%         ymCellValue = @(x) cellfun(@value, x, 'UniformOutput', false);
%         iterProgress = struct('a', a_0, 'b', b_0, 'gam', ones(size(gam))*1e6);
        iterProgress = struct('a', a_0, 'b', b_0, 'gam', cell2mat(ymCellValue(gam)));
        k = 1;
        iterChange = 1;
        while k<maxIter && any(iterChange>=epsilon)
            % Convexification of problem (Step (2) from [1]).
            ymConstraintNL = replace(ymConstraint, [a_op b_op], [iterProgress(k).a iterProgress(k).b]);
            % Optimization (Equations (10)-(12) from [1]).
            ymDiagnostics = sdsyn_blinhandler(ymConstraintNL, ymObjective + kappa*norm(a)^2, ymOptions);
            if ymDiagnostics.problem==1
                warning('sdsyn_nlinhandler:infeas', ['Iteration ' num2str(k) ' is infeasible.']);
            elseif ymDiagnostics.problem~=0
                warning(ymDiagnostics.info);
            end
            k = k + 1; % Update k.
            % Assign solution for next convexification point (Equation (13)
            % from [1])
            iterProgress(k).a = iterProgress(k-1).a + value(a);
            iterProgress(k).b = iterProgress(k-1).b + value(b);
            %ymConstraintNL2 = replace(ymConstraint, [a b], [zeros(size(a)) zeros(size(b))]);
            %ymConstraintNL2 = replace(ymConstraintNL2, [a_op b_op], [iterProgress(k).a iterProgress(k).b]);
            %sdsyn_blinhandler(ymConstraintNL2, ymObjective, ymOptions);
            iterProgress(k).gam = cell2mat(ymCellValue(gam));
            iterChange = norm(iterProgress(k).a - iterProgress(k-1).a);
        end
    end

    function diagn = sdsyn_blinhandler(constraint, objective, options)
    %SDSYN_BLINHANDLER Bilinear Optimization Handler using Bisection Method
    %   Handles calls to YALMIP optimize() to allow the bisection method to
    %   be used with different solvers. Bisection method is required for
    %   the line search over alpha if an l-1 norm constraint is in effect.
    %
    %   Inputs:
    %       constraint: YALMIP constraint object.
    %       objective: YALMIP objective function.
    %       options: YALMIP options structure defined with sdpsettings().
    %
    %   Outputs:
    %       diag: YALMIP diagnostics structure.
    %
    %   See also: SDSYN, SDSYN_BLINHANDLER
    %       [1] J. Löfberg, “YALMIP: A Toolbox for Modeling and 
    %           Optimization in MATLAB,” In Proceedings of the CACSD 
    %           Conference. Taipei, Taiwan, 2004.
    %       
    %   $Author: BH $    $Date: 2018-08-16 $    $Revision: 0 $
        iAlpha = find(cellfun(@(x) isa(x, 'sdpvar'), alpha));
        if any(iAlpha)
            if isVerbose
                lineSearchOptions = optimset('Display', 'iter', 'FunValCheck', 'off');
                options.verbose = false;
            else
                lineSearchOptions = optimset('Display', 'off', 'FunValCheck', 'off');
            end
            if strcmpi(options.solver, 'lmilab') % Must use YALMIP optimize with variable replacement (slower).
                alpha_opt = fminsearch(@(x) sdsyn_blinhandler_lmilab(x), zeros(length(alpha), 1), lineSearchOptions);
                options.verbose = isVerbose; % Restore verbosity.
                [~, diagn] = sdsyn_blinhandler_lmilab(alpha_opt);
            else % Supports YALMIP optimizer pre-compiled object (much faster).
                lineSearch = optimizer(constraint, objective, options, alpha, ymObjective);
                alpha_opt = fminbnd(@(x) lineSearch(x), 0, 1, lineSearchOptions);
                options.verbose = isVerbose; % Restore verbosity.
                constraint = replace(constraint, alpha, alpha_opt);
                diagn = optimize(constraint, objective, options);
            end
        else 
            diagn = optimize(constraint, objective, options);
        end
        
        function [objective_val, diagn] = sdsyn_blinhandler_lmilab(alpha_val)
        %SDSYN_BLINHANDLER_LMILAB Bilinear SDP handler for LMILAB solver
        %   This function handles the case where LMILAB is used, because
        %   the YALMIP optimizer() function is not compatible with LMILAB.
        %   Instead, the bilinear variable is manually replaced and the SDP
        %   solved. This is very slow!
        %
        %   Input:
        %       alpha_val: vector of alpha value query points, one value
        %           for each l-1/(*) constraint.
        %
        %   Outputs:
        %       objective_val: the value of the YALMIP objective with the
        %           bilinear SDP evaluated at the given alpha_val.
        %       diag: YALMIP diagnostics object.
        %
        %   See also: SDSYN_BLINHANDLER, OPTIMIZER
 
            for ia=iAlpha
                linearConstraint = replace(constraint, alpha{ia}, alpha_val(ia));
            end
            diagn = optimize(linearConstraint, objective, options);
            objective_val = value(objective);
        end
    end

    function [H, S, CL] = sdsyn_tolti(b_val, a_val, CL, delta)
    %SDSYN_TOLTI Convert optimization variables to LTI object
    %   Forms a MATLAB LTI object from the optimization variables, the LTI
    %   object is in the form of the prototype model H_0 passed into SDSYN.
    %
    %   Inputs:
    %       b_val: the optimal loop filter transfer function numerator 
    %           coefficients.
    %       a_val: the optimal loop filter transfer function denominator 
    %           coefficients.
    %       CL: the closed-loop state space form with uncertainty and
    %           reference input channels and uncertainty and feedback error
    %           output channels. Uncertain block is not present.
    %       delta: the norm-bounded uncertainty extracted by LFT.
    %
    %   Outputs: 
    %       H: the loop filter in the form specified by H_0 (tf, zpk, or ss
    %           model). Does not include the uncertain gain in series.
    %       S: the sensitivity function in the form specified by H_0 (tf,
    %           zpk, or ss model). Does not include the uncertain gain in
    %           series.
    %       CL: the closed-loop state space form with uncertainty included
    %           by LFT.
    %
    %   See also: SDSYN
    %       
    %   $Author: BH $    $Date: 2018-08-16 $    $Revision: 0 $
        ymCellReplace = @(x) cellfun(@(y) replace(y, [a b], [a_val b_val]), x, 'UniformOutput', false);
        A_cl = replace(CL.A, [a b], [a_val b_val]);
        B_cl = cell2mat(ymCellReplace(CL.B));
        C_cl = cell2mat(ymCellReplace(CL.C));
        D_cl = cell2mat(ymCellReplace(CL.D));
        CL = ss(A_cl, B_cl, C_cl, D_cl, Ts);
        CL = lft(delta, CL);
        if isa(H_0, 'tf') || isa(H_0, 'zpk')
            H = tf(fliplr(b_val'), [1 fliplr(a_val')], Ts);
            S = tf([1 fliplr(a_val')], [1 fliplr(a_val') + fliplr(b_val')], Ts);
            if isa(H_0, 'zpk')
                H = zpk(H);
                S = zpk(S);
            end
        elseif isa(H_0, 'ss') % Controllable canonical form representation.
            H = ss([zeros(n-1, 1) eye(n-1); -a_val'], [zeros(n-1, 1); 1], b_val', 0, Ts);
            S = ss([zeros(n-1, 1) eye(n-1); -(a_val + b_val)'], [zeros(n-1, 1); 1], -b_val', 1, Ts);
        end
    end

end

function [synOpt, CL] = sdsyn_buildplant(synOpt, b, a, K)
%PARTITIONPLANT Build Controllable Canonical Form Generalized Plant
%   Generates a state space model of a closed-loop plant with uncertain 
%   quantizer gain. Inputs to the plant are the reference input and
%   norm-bounded uncertainty channel. Outputs of the plant are the
%   feedback error (sensitivity function) and uncertainty output.
%   
%   Inputs:
%       synOpt: optimization target structure generated with
%           DEFINESYNOPT.M (required).
%       b: vector of loop filter numerator coefficients, ordered from
%           lowest to highest power of z/s.
%       a: vector of loop filter denominator coefficients, ordered from
%           lowest to highest power of z/s.
%       K: 2-by-2 matrix of certain gains obtained by lower LFT of the
%           uncertain quantizer gain.
%       
%   Outputs:
%       synOpt: optimization target structure augmented with A, B, C, D
%           matrices for each optimization target.
%       CL: structure containing SDPVAR state space data for the closed
%           loop augmented system.
%
%   $Author: BH $    $Date: 2018-08-15 $    $Revision: 2 $
%
% REVISION 0:
%   2018-08-15 by BH: forked from PARTITIONPLANT in file dtsyn_func.m with
%       major modifications.
%
% REVISION 1:
%   2018-08-21 by BH: added CL output.
%
% REVISION 2:
%   2018-09-19 by BH: added complementary sensitivity function (STF)
%   channel r -> y in an attempt to fix the slow roll-off issue for
%   direct continuous-time design of sigma delta modulators. Also,
%   simplified logic.

    %% Plant partitioning:
    %       In:        x          w      r
    %                  |          |      |
    % Out:             V          V      V  Sz:
    %          .---------------------------.
    % dx/dt <- | -A-K_22*B*C | -K_21*B | B | n
    %          |-------------|---------|---|
    %     z <- |    K_12*C   |   K_11  | 0 | 1
    %          |-------------|---------|---|
    %     e <- |   -K_22*C   |  -K_21  | 1 | 1
    %          |-------------|---------|---|
    %     y <- |    K_22*C   |   K_21  | 0 | 1
    %          '---------------------------'
    %      Size:       n          1      1
    
    n = length(a);
    A = [zeros(n-1, 1) eye(n-1); -a' - K(2,2)*b'];
    B{1,1} = -K(2,1)*[zeros(n-1, 1); 1];
    B{1,2} = [zeros(n-1, 1); 1];
    C{1,1} = K(1,2)*b';
    C{2,1} = -K(2,2)*b';
    C{3,1} = K(2,2)*b';
    D{1,1} = K(1,1);
    D{1,2} = 0;
    D{2,1} = -K(2,1);
    D{2,2} = 1;
    D{3,1} = K(2,1);
    D{3,2} = 0;
    CL = struct('A', A, 'B', {B}, 'C', {C}, 'D', {D});
    
    for iStruct=1:length(synOpt)
        if ~any(synOpt(iStruct).inputChannel == [1 2])
            error('filtersynthesis_buildplant:invalidInputChannel', ['Input channel ' num2str(synOpt(iStruct).inputChannel) ' is invalid (must be 1 or 2).']);
        end
        if ~any(synOpt(iStruct).outputChannel == [1 2 3])
            error('filtersynthesis_buildplant:invalidOutputChannel', ['Output channel ' num2str(synOpt(iStruct).outputChannel) ' is invalid (must be 1, 2, or 3).']);
        end
        synOpt(iStruct).A = A; % A is global and redundant, but kept in structure array for completeness.
        synOpt(iStruct).B = B{synOpt(iStruct).inputChannel};
        synOpt(iStruct).C = C{synOpt(iStruct).outputChannel};
        synOpt(iStruct).D = D{synOpt(iStruct).outputChannel, synOpt(iStruct).inputChannel};
    end
end

function sdsyn_checkconvergence(iterProgress)
    a = [iterProgress(:).a]';
    b = [iterProgress(:).b]';
    gam = reshape([iterProgress(:).gam]', [size(iterProgress(1).gam, 2) length(iterProgress)])';
end