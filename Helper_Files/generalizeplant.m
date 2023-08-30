function P_gen = generalizeplant(varargin)
%GENERALIZEPLANT Create Generalized Plant.
%   Creates a generalized plant formulation for control problems by
%   augmenting the plant with weighting functions on the sensitivity
%   function, controller output, and complementary sensitivity function.
%
%   Inputs:
%       P: LTI plant model. Type of P model (ss, tf, etc.) determines type 
%           of P_gen model. *Function has only been tested with P type ss*
%       W_S: (Optional) desired sensitivity function LTI model. If no 
%           weight is desired, pass the empty matrix [].
%       W_R: (Optional) desired control signal weight LTI model. If no 
%           weight is desired, pass the empty matrix [].
%       W_T: (Optional) desired complementary sensitivity function LTI 
%           model. If no weight is desired, pass the empty matrix [].
%       DOF: (Optional) name-value parameter to specify the degrees of
%           freedom of the controller. DOF must be 1 (default) or 2.
%       UncertaintySz: (Optional) name-value parameter to specify the 
%           number of uncertainty block outputs, inputs [nz nw]. If P is 
%           type USS, UncertaintySz is the size of the LFT extracted 
%           uncertainty block. Default: [0 0].
%
%   Outputs:
%       P_gen: generalized plant with inputs:
%                 [w; r; u]
%                 uncertain:
%                       w: uncertainty channel(s).
%                 reference:
%                       r: reference signal.
%                 control:
%                       u: control input signal.
%              and outputs:
%                   [z; e_s; e_r; e_t; e] (1 DOF), or
%                   [z; e_s; e_r; e_t; r; y] (2 DOF)
%                   uncertain:
%                       z: uncertainty channel(s).
%                   performance:
%                       e_s, e_r, e_t: performance channel(s).
%                   error:
%                       e: error signal (1 DOF), or
%                       r: reference signal, and
%                       y: output feedback signal.                       
%
%              and (if state space) states: 
%                   [x; x_S; x_R; x_T]
%                   x: system states.
%                   x_S, x_R, x_T: states for weighting functions W_S, W_R,
%                   W_T.
%   
%   See also: AUGW
%       [1] R. Nagamune, "Derivations of the Generalized Plant", 
%           MECH528/EECE508 Supplementary Note. 2017. 
%
%   $Author: BH $    $Date: 2017-10-13 $    $Revision: 2 $
%
% REVISION 1:
%   2018-03-02 by BH: Changed behaviour so the weighting functions aren't
%   inverted to match built-in function AUGW, removed stability checks from
%   inputs (as inverse of weighting functions aren't always stable).
%
% REVISION 2:
%   2018-04-11 by BH: Fixed bug where sampling time of sensitivity function
%   was being compared to the plant for all weighting function inputs.

    %% Parse inputs.
    p = inputParser;
    addRequired(p, 'P', @(x) isa(x, 'lti'));
    addOptional(p, 'W_S', [], @(x) isa(x, 'lti') || isempty(x));
    addOptional(p, 'W_R', [], @(x) isa(x, 'lti') || isempty(x));
    addOptional(p, 'W_T', [], @(x) isa(x, 'lti') || isempty(x));
    addParameter(p, 'DOF', 1, @(x) x==1 || x==2);
    addParameter(p, 'UncertaintySz', [0 0], @(x) isnumeric(x) && isrow(x) && length(x)<=2);
    parse(p, varargin{:});
    P = p.Results.P;
    % If empty, modify size. Absence of weighting function should be an LTI
    % system with 1 input, 0 outputs.
    if isempty(p.Results.W_S)
        W_S = zeros(0, 1);
    else
        W_S = p.Results.W_S;
        if W_S.Ts~=P.Ts
            error('generalizeplant:ts', 'Sampling times must agree.')
        end
    end
    if isempty(p.Results.W_R)
        W_R = zeros(0, 1);
    else
        W_R = p.Results.W_R;
        if W_R.Ts~=P.Ts
            error('generalizeplant:ts', 'Sampling times must agree.')
        end
    end
    if isempty(p.Results.W_T)
        W_T = zeros(0, 1);
    else
        W_T = p.Results.W_T;
        if W_T.Ts~=P.Ts
            error('generalizeplant:ts', 'Sampling times must agree.')
        end
    end
    UncertaintySz = p.Results.UncertaintySz;
    if length(UncertaintySz)==1
        UncertaintySz = [UncertaintySz UncertaintySz];
    end
    DOF = p.Results.DOF;

    %% Take upper LFT of system.
    if isa(P, 'uss')
        [M, Delta] = lftdata(P);
        UncertaintySz(1) = size(Delta, 1);
        UncertaintySz(2) = size(Delta, 2);
    else
        M = P; Delta = [];
    end

    %% Convert models to the same type.
    if isa(M, 'ss')
        if ~isa(W_S, 'ss') 
            W_S = ss(W_S);
        end
        if ~isa(W_R, 'ss')
            W_R = ss(W_R);
        end
        if ~isa(W_T, 'ss')
            W_T = ss(W_T);
        end
    else
        if isa(W_S, 'ss')
            W_S = tf(W_S);
        end
        if isa(W_R, 'ss')
            W_R = tf(W_R);
        end
        if isa(W_T, 'ss')
            W_T = tf(W_T);
        end
    end

    %% Determine sizes of signals.
    if isa(M, 'ss')
        nx = size(M.A, 1);
        nS = size(W_S.A, 1);
        nR = size(W_R.A, 1);
        nT = size(W_T.A, 1);
    end
    nw = UncertaintySz(2);
    nz = UncertaintySz(1);
    if isa(M, 'ss')
        nu = size(M.B, 2) - nw;
        ny = size(M.C, 1) - nz;
        ne_S = size(W_S.C, 1);
        ne_R = size(W_R.C, 1);
        ne_T = size(W_T.C, 1);
    else
        nu = size(M, 2) - nw;
        ny = size(M, 1) - nz;
        ne_S = size(W_S, 1);
        ne_R = size(W_R, 1);
        ne_T = size(W_T, 1);
    end
    if DOF==1
        ne = ny;
    end
    nr = ny;

    %% Compute generalized plant.
    if isa(M, 'ss') % State space format.
        A = M.A;
        B{1} = M.B(:,1:nw);
        B{2} = M.B(:,nw+1:end);
        C{1} = M.C(1:nz,:);
        C{2} = M.C(nz+1:end,:);
        D{1,1} = M.D(1:nz,1:nw);
        D{1,2} = M.D(1:nz,nw+1:end);
        D{2,1} = M.D(nz+1:end,1:nw);
        D{2,2} = M.D(nz+1:end,nw+1:end);

        P_gen = ss;
        P_gen.Ts = P.Ts;
        P_gen.A = [A zeros(nx, nS + nR + nT);...
            -W_S.B*C{2} W_S.A zeros(nS, nR + nT);...
            zeros(nR, nx + nS) W_R.A zeros(nR, nT);
            W_T.B*C{2} zeros(nT, nS + nR) W_T.A];
        P_gen.B = [B{1} zeros(nx, nr) B{2};...
            -W_S.B*D{2,1} W_S.B -W_S.B*D{2,2};...
            zeros(nR, nw + nr) W_R.B;...
            W_T.B*D{2,1} zeros(nT, nr) W_T.B*D{2,2}];
        if DOF==1
            P_gen.C = [C{1} zeros(nz, nS + nR + nT);...
                -W_S.D*C{2} W_S.C zeros(ne_S, nR + nT);...
                zeros(ne_R, nx + nS) W_R.C zeros(ne_R, nT);...
                W_T.D*C{2} zeros(ne_T, nS + nR) W_T.C;...
                -C{2} zeros(ne, nS + nR + nT)];
        else
            P_gen.C = [C{1} zeros(nz, nS + nR + nT);...
                -W_S.D*C{2} W_S.C zeros(ne_S, nR + nT);...
                zeros(ne_R, nx + nS) W_R.C zeros(ne_R, nT);...
                W_T.D*C{2} zeros(ne_T, nS + nR) W_T.C;...
                zeros(nr, nx + nS + nR + nT);...
                C{2} zeros(ny, nS + nR + nT)];
        end
        if DOF==1
            P_gen.D = [D{1,1} zeros(nz, nr) D{1,2};...
                -W_S.D*D{2,1} W_S.D -W_S.D*D{2,2};...
                zeros(ne_R, nw + nr) W_R.D;...
                W_T.D*D{2,1} zeros(ne_T, nr) W_T.D*D{2,2};...
                -D{2,1} eye(ne, nr) -D{2,2}];
        else
            P_gen.D = [D{1,1} zeros(nz, nr) D{1,2};...
                -W_S.D*D{2,1} W_S.D -W_S.D*D{2,2};...
                zeros(ne_R, nw + nr) W_R.D;...
                W_T.D*D{2,1} zeros(ne_T, nr) W_T.D*D{2,2};...
                zeros(nr, nw) eye(nr) zeros(nr, nu);...
                D{2,1} zeros(ny, nr) D{2,2}];
        end
        % Assign state names: 
        %   x: system states.
        %   x_S: states for sensitivty weighting function.
        %   x_R: states for control signal weighting function.
        %   x_T: states for complementary sensitivity weighting function.
      	nameState = nameSignal(nx, 'x', nS, 'x_S', nR, 'x_R', nT, 'x_T');
        P_gen.StateName = nameState;
    else % Transfer function format.
        P_gen = [M(1:nz, 1:nw) tf(zeros(nz, nr)) M(1:nz, nw+1:end);...
            -W_S*M(nz+1:end, 1:nw) W_S -W_S*M(nz+1:end, nw+1:end);...
            tf(zeros(ne_R, nw + nr)) W_R;...
            W_T*M(nz+1:end, 1:nw) tf(zeros(ne_T, nr)) W_T*M(nz+1:end, nw+1:end);];
        if DOF==1
            P_gen = [P_gen;...
                -M(nz+1:end, 1:nw) tf(eye(ne, nr)) -M(nz+1:end, nw+1:end)];
        else
            P_gen = [P_gen;...
                tf(zeros(nr, nw)) tf(eye(nr)) tf(zeros(nr, nu));...
                M(nz+1:end, 1:nw) tf(zeros(ny, nr)) M(nz+1:end, nw+1:end)];
        end
    end

    %% Assign signal names.
    % Inputs:
    %   uncertain:
    %       w: uncertainty channel(s).
    %   reference:
    %       r: reference signal.
    %   control:
    %       u: control input signal.
    % Outputs:
    %   uncertain:
    %       z: uncertainty channel(s).
    %   performance:
    %       e_s: sensitivity function performance channel.
    %       e_r: control signal weighting performance channel.
    %       e_t: complementary sensitivity function performance channel.
    %   error:
    %       e: error signal.
    
    % Signal names.
    if DOF==1
        nameOutput = nameSignal(nz, 'z', ne_S, 'e_S', ne_R, 'e_R', ne_T, 'e_T', ne, 'e');
    else
        nameOutput = nameSignal(nz, 'z', ne_S, 'e_S', ne_R, 'e_R', ne_T, 'e_T', nr, 'r', ny, 'y');
    end
    P_gen.OutputName = nameOutput;
    nameInput = nameSignal(nw, 'w', nr, 'r', nu, 'u');
    P_gen.InputName = nameInput;
    % Input/output groups.
    P_gen.InputGroup.uncertain = 1:nw;
    P_gen.InputGroup.reference = nw+1:nw+nr;
    P_gen.InputGroup.control = nw+nr+1:nw+nr+nu;
    P_gen.OutputGroup.uncertain = 1:nz;
    P_gen.OutputGroup.performance = nz+1:nz+ne_S+ne_R+ne_T;
    if DOF==1
        P_gen.OutputGroup.error = nz+ne_S+ne_R+ne_T+1:nz+ne_S+ne_R+ne_T+ne;
    else
        P_gen.OutputGroup.error = nz+ne_S+ne_R+ne_T+1:nz+ne_S+ne_R+ne_T+nr+ny;
    end
    P_gen.Notes = 'Generalized plant produced by GENERALIZEPLANT.M.';
    
    if ~isempty(Delta)
        P_gen = lft(Delta, P_gen);
    end

        %% Helper functions.
        function name = nameSignal(varargin)
        %NAMESIGNAL Compile signal names.
        %   Returns a cell array of signal names from an ordered list of 
        %   prefixes and channel sizes (of any length).
        %
        %   EXAMPLE 
        %       name = nameSignal(count1, prefix1, count2, prefix2, ....)
        %
        %   Inputs:
        %       count: scalar number of signals in given channel.
        %       prefix: string of prefix for channel names.
        %   Output:
        %       name: cell array of signal names.
        %
        %   $Author: BH $    $Date: 2017-10-13 $    $Revision: 1 $
            name = cell(1, sum(horzcat(varargin{1:2:end})));
            iName = 1;
            for ii=1:2:length(varargin)
                if varargin{ii}~=0
                    name(iName:iName+varargin{ii}-1) = cellfun(@(x) sprintf([varargin{ii+1} '(%u)'], x), num2cell(1:varargin{ii}), 'UniformOutput', false);
                    iName = iName+varargin{ii};
                end
            end
        end

end