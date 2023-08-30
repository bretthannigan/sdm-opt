classdef Sdm
    properties (Access = public)
        %% Loop filter and quantizer.
        H @lti
        H_pre = tf(1);
        q @SdmQuantizer
        %% Sampling.
        Ts @double
        osr @double
    end
    properties (Dependent = true)
        R % Control senstivity function / error transfer function (ETF).
        S % Sensitivity function / noise transfer function (NTF).
        T % Complementary sensitivity function / signal transfer function (STF).
    end
    properties (Dependent = true, Hidden = true)
        NTF % Sensitivity function / noise transfer function (NTF).
        STF % Complementary sensitivity function / signal transfer function (STF).
    end
    methods
        function obj = Sdm(varargin)
        end
        function obj = set.T(obj, T)
            obj.H = minreal([T/(1 - T) -T/(1 - T)]); % Assume 1-DOF.
        end
        function T = get.T(obj)
            T = minreal(obj.H(1)/(1 - obj.H(2)));
        end
        function obj = set.S(obj, S)
            obj.H = [minreal((1 - S)/S) minreal(-(1 - S)/S)]; % Assume 1-DOF.
        end
        function S = get.S(obj)
            S = 1/(1 - obj.H(2));
        end
        function obj = set.H(obj, H)
            if size(H, 1)~=1
                error('Sdm:setH_nOutput', 'Loop filter H must have 1 output.');
            end
            if size(H, 2)==1 % 1-DOF loop filter.
                obj.H = [H -H];
            elseif size(H, 2)==2 % 2-DOF loop filter.
                obj.H = H;
            else
                error('Sdm:setH_nInput', 'Loop filter H must have 1 or 2 inputs.');
            end
        end
        function H = get.H(obj)
            H = obj.H;
        end
        function obj = setST(S, T)
            obj.H(1) = minreal(T/S); 
            obj.H(2) = minreal(1 - 1/S);
        end
        function obj = set.STF(obj, STF)
            obj.T = STF;
        end
        function STF = get.STF(obj)
            STF = obj.T;
        end
        function obj = set.NTF(obj, NTF)
            obj.S = NTF;
        end
        function NTF = get.NTF(obj)
            NTF = obj.S;
        end
        function sim_out = simulate(obj, varargin)
            sim_in = struct;
            if nargin==2
                sim_in(1).inputSignal = varargin{1};
                sim_in(1).isLinearModel = false;
            elseif nargin==3
                sim_in(1).inputSignal = timeseries(varargin{2}, varargin{1}, 'Name', 'r');
                sim_in(1).isLinearModel = false;
            elseif nargin==4
                sim_in(1).inputSignal = timeseries(varargin{2}, varargin{1}, varargin{3});
                sim_in(1).isLinearModel = false;
            elseif nargin==5
                sim_in(1).inputSignal = timeseries(varargin{2}, varargin{1}, varargin{3});
                sim_in(1).isLinearModel = varargin{4};
            else
            end
            sim_out = obj.simHandler(sim_in, true);
        end
        function [mag, phase, w] = bode(obj, varargin)
            if nargin==2
                w = varargin{1};
            else
                if ~isempty(obj.Ts)
                    w = {(2*pi/obj.Ts)/(obj.osr^1.25) pi/obj.Ts};
                else
                   	[~, ~, Tw] = bode(obj.T);
                    [~, ~, Sw] = bode(obj.S);
                    w = {min(Tw(1), Sw(1)), max(Tw(end), Sw(end))};
                end
            end
            if nargout==0
                f = figure;
                bode(obj.T, w, '-b')
                hold on;
                bode(obj.S, w, '-r')
                ax = findall(f, 'type', 'axes');
                yLimits = ax(2).YLim;
                if ~isempty(obj.osr)
                    line(ax(2), [2*pi/(obj.Ts*obj.osr) 2*pi/(obj.Ts*obj.osr)], [-1e10 1e10], 'LineStyle', '--', 'Color', 'k');
                    ylim(ax(2), yLimits);
                    legend(ax(2), 'STF', 'NTF', 'OSR');
                else
                    legend(ax(2), 'STF', 'NTF');
                end
                yLimits = ax(3).YLim;
                if ~isempty(obj.osr)
                    line(ax(3), [2*pi/(obj.Ts*obj.osr) 2*pi/(obj.Ts*obj.osr)], [-1e10 1e10], 'LineStyle', '--', 'Color', 'k');
                    ylim(ax(3), yLimits);
                    legend(ax(3), 'STF', 'NTF', 'OSR');
                else
                    legend(ax(3), 'STF', 'NTF');
                end
            else
                [Tmag, Tphase, ~] = bode(obj.T, w);
                [Smag, Sphase, ~] = bode(obj.S, w);
                mag = [squeeze(Tmag) squeeze(Smag)];
                phase = [squeeze(Tphase) squeeze(Sphase)];
            end
        end
    end
    methods (Access = private)
        function sim_out = simHandler(obj, sim_in, isVerbose)
            if isVerbose
                fprintf('Opening Simulink model... ')
                open_system('SdmTestbench');
            else
                load_system('SdmTestbench');
            end
            simWs = get_param('SdmTestbench', 'modelworkspace');
            simWs.clear;
            assignin('base', 'sdm_q_sim', obj.q); % Workaround: class cannot be read in from model workspace.
            simWs.assignin('H', obj.H);
            simWs.assignin('H_pre', obj.H_pre);
            simWs.assignin('Ts', obj.Ts);
            if isVerbose
                fprintf(' Done.\n');
            end
            for iSim=1:length(sim_in)
                if isVerbose
                    fprintf('\tSimulating case %d of %d... ', iSim, length(sim_in));
                end
                simWs.assignin('x', sim_in(iSim).inputSignal);
                simWs.assignin('isLinearModel', sim_in(iSim).isLinearModel);
                if strcmp(sim_in(iSim).inputSignal.name, 'r')
                    simWs.assignin('isReferenceInput', true);
                elseif strcmp(sim_in(iSim).inputSignal.name, 'd')
                    simWs.assignin('isReferenceInput', false);
                else
                    error('simHandler:InputSignalName', ['No input signal of name ' sim_in(iSim).InputSignal.name ' exists.']);
                end
                simWs.assignin('isReferenceInput', strcmp(sim_in(iSim).inputSignal.name, 'r'));
                if obj.H.Ts==0
                    sim_out = sim('SdmTestbench', 'StopTime', num2str(sim_in(iSim).inputSignal.Time(end)), 'Solver', 'ode45');
                else
                    sim_out = sim('SdmTestbench', 'StopTime', num2str(sim_in(iSim).inputSignal.Time(end)), 'Sovler', 'VariableStepDiscrete');
                end
                sim_in(iSim).outputSignal(1) = sim_out.y;
                sim_in(iSim).outputSignal(2) = sim_out.u;
                sim_in(iSim).outputSignal(3) = sim_out.e;
                if isVerbose
                    fprintf('Done.\n');
                end
            end
            sim_out = sim_in;
            clear('sdm_q_sim')
            if isVerbose
                sprintf('Completed %d simulations successfully.\n', length(sim_in));
            end
        end
        function warn = checkProperties(obj, isSilent)
            warn = 0;
            % Initialized properties:
            if isempty(obj.H)
                if ~isSilent
                    warning('Loop filter H is not specified.');
                end
                warn = warn + 1;
            end
            if isempty(obj.Ts)
                if ~isSilent
                    warning('Sample time Ts is not specified.');
                end
                warn = warn + 1;
            end
            if ~isempty(obj.H) && ~isempty(obj.Ts)
                if obj.H.Ts~=0 && obj.H.Ts~=obj.Ts
                    if ~isSilent
                        warning('Sample times of discrete-time loop filter and system do not match.');
                    end
                    warn = warn + 1;
                end
            end
            if isempty(obj.osr)
                if ~isSilent
                    warning('Oversampling rate OSR is not specified.');
                end
                warn = warn + 1;
            else
                if obj.osr<1
                    if ~isSilent
                        warning('Oversampling ratio is less than 1.');
                    end
                    warn = warn + 1;
                end
            end
        end
    end
    methods (Static = true)
%         function sim_out = simulationHandler(sim_in)
%             
%         end
    end
end