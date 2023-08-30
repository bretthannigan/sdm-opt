classdef SdmQuantizer
   	%SDMQUANTIZER Unit quantizer superclass.
    %   The simplest version of a quantizer. Has quantization spacing of 1
    %   centred on (0,0) and has infinite range.
    % SdmQuantizer Family Class Inheritance Hierarchy:
    %    
    % Base class:                           SdmQuantizer(this)
    %                      ,---------------------|-----------------.
    % Infinite:    SdmUInfQuantizer       SdmLInfQuantizer         |
    %                      |                     |                 |
    % Saturating:    SdmUQuantizer         SdmLQuantizer      SdmAQuantizer
    %                      
    % Spacing:        ^Uniform^            ^Logarithmic^       ^Arbirtary^
    %
    %   $Author: BH $    $Date: 2017-11-16 13:52:00 $    $Revision: 1 $
    properties (Access = public)
        roundHalfTo@double scalar % A numeric scalar defining the type of rounding.
        inputName@char % A character vector for the input channel name.
        outputName@char % A character vector for the output channel name.
        inputUnit@char % A character vector for the input channel unit.
        outputUnit@char % A character vector for the output channel unit.
    end
    methods
        function obj = SdmQuantizer(varargin)
            %SDMQUANTIZER SdmQuantizer object constructor.
            %   Constructs an SdmQuantizer object.
            %
            %   Inputs:
            %       RoundHalfTo: name-value parameter for numeric scalar 
            %       decribing rounding behaviour when the input is 0.5*n 
            %       where n is an integer.
            %           Inf: round toward infinity.
            %           -Inf: round toward negative infinity.
            %           0: round toward zero, if zero, round toward
            %           positive infinity.
            %           1 (default): round away from zero, if zero, round 
            %           toward negative infinity.
            %           2: round toward even (TODO).
            %           -2: round toward odd (TODO).
            %           NaN: do not round.
            %           otherwise: round up or down randomly.
            %       InputName: name-value parameter for character vector
            %       describing the input channel name (default: '').
            %       InputUnit: name-value parameter for character vector
            %       describing the input channel unit (default: '').
            %       OutputName: name-value parameter for character vector
            %       describing the output channel name (default: '').
            %       OutputUnit: name-value parameter for character vector
            %       describing the outptu channel unit (default: '').
            %   Output:
            %       obj: An SdmQuantizer object.
            %
            %   $Author: BH $    $Date: 2017-11-15 $    $Revision: 1 $
            p = inputParser;
            addParameter(p, 'RoundHalfTo', 1, @isscalar);
            addParameter(p, 'InputName', '', @ischar);
            addParameter(p, 'OutputName', '', @ischar);
            addParameter(p, 'InputUnit', '', @ischar);
            addParameter(p, 'OutputUnit', '', @ischar);
            parse(p, varargin{:});
            obj.roundHalfTo = p.Results.RoundHalfTo;
            obj.inputName = p.Results.InputName;
            obj.outputName = p.Results.OutputName;
            obj.inputUnit = p.Results.InputUnit;
            obj.outputUnit = p.Results.OutputUnit;
        end
        function y = quantize(obj, u)
            %QUANTIZE Quantize the input.
            %
            %   Rounds the input to the nearest integer based on the 
            %   rounding behaviour defined by 'RoundHalfTo'.
            %
            %   Input:
            %       u: numeric vector input.
            %   Output:
            %       y: quantized numeric vector output.
            %
            %   See also: ROUND
            %
            %   $Author: BH $    $Date: 2017-11-14 $    $Revision: 2 $
            %
            % REVISION 2:
            %   2017-11-15 by BH: Modified for the base class quantizer
            %   with spacing=1.
            rem = abs(mod(u, 1));
            y = u - rem;
            isRoundedUp = false(size(u));
            if obj.roundHalfTo==Inf % Toward Inf.
                isRoundedUp = rem >= 0.5;
            elseif obj.roundHalfTo==-Inf % Toward -Inf.
                isRoundedUp = rem > 0.5;
            elseif obj.roundHalfTo==0 % Toward 0, if zero then toward Inf.
                isRoundedUp(u<0) = rem(u<0) > 0.5;
                isRoundedUp(u>=0) = rem(u>=0) >= 0.5;
            elseif obj.roundHalfTo==1 % Away from 0, if zero then toward Inf.
                isRoundedUp(u<0) = rem(u<0) > 0.5;
                isRoundedUp(u>=0) = rem(u>=0) >= 0.5;
            elseif obj.roundHalfTo==2 % Toward even (TODO).
            elseif obj.roundHalfTo==-2 % Toward odd (TODO).
            elseif isnan(obj.roundHalfTo) % Do not round (adds new quantization level).
                isRoundedUp = rem > 0.5;
                y(rem == 0.5) = y(rem == 0.5) + 0.5;
            else % Randomly round.
                isRoundedUp = rem > 0.5;
                isRoundedUp(rem == 0.5) = randi([0 1], [1 sum(rem == 0.5)]);
            end
            y(isRoundedUp) = y(isRoundedUp) + 1;
        end
    end
    methods (Abstract = true)
        h = plot(obj)
    end
end