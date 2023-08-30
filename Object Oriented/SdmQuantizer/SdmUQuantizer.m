classdef SdmUQuantizer < SdmUInfQuantizer
    %SDMUQUANTIZER Uniform saturating quantizer subclass.
    %   Implements a quantizer with uniform spacing and saturation by
    %   limiting the output of the SdmUInfQuantizer superclass.
    % SdmQuantizer Family Class Inheritance Hierarchy:
    %    
    % Base class:                           SdmQuantizer
    %                      ,---------------------|-----------------.
    % Infinite:    SdmUInfQuantizer       SdmLInfQuantizer         |
    %                      |                     |                 |
    % Saturating:    SdmUQuantizer(this)   SdmLQuantizer      SdmAQuantizer
    % 
    % Spacing:        ^Uniform^            ^Logarithmic^       ^Arbirtary^
    %
    %   $Author: BH $    $Date: 2017-11-16 14:22:00 $    $Revision: 1 $
    properties (SetAccess = protected, Hidden = true)
        uMin % Numeric scalar input range lower bound.
        uMax % Numeric scalar input range upper bound.
    end
    properties (Dependent = true)
        uRange % 2-element numeric vector [uMin uMax].
        nLev % Numeric scalar number of quantization levels.
        levels % Numeric vector of quantization level values.
        bits % Numeric scalar number of quantization bits.
    end
    properties (Access = private, Constant = true)
        intTol = 1e-8; % Numeric scalar tolerance for integer check.
    end
    methods
        function obj = SdmUQuantizer(varargin)
            %SDMUQUANTIZER SdmUQuantizer object constructor.
            %   Constructs an SdmUQuantizer object.
            %
            %   Inputs:
            %       nLev: optional scalar positive integer number of 
            %       quantization levels (default:2).
            %       uMin: optional scalar numeric lower saturation level
            %       (default:-1).
            %       uMax: optional scalar numeric upper saturation level
            %       (default:1).
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
            %   Output:
            %       obj: An SdmUQuantizer object.
            %
            %   See also: SDMUINFQUANTIZER
            %
            %   $Author: BH $    $Date: 2017-11-15 $    $Revision: 1 $
            p = inputParser;
            p.KeepUnmatched = true;
            addOptional(p, 'nLev', 2, @isscalar);
            addOptional(p, 'uMin', -1, @isscalar);
            addOptional(p, 'uMax', 1, @isscalar);
            addParameter(p, 'roundHalfTo', 1, @isscalar);
            parse(p, varargin{:});
            if p.Results.uMax<p.Results.uMin
                error('SdmUQuantizer:set_uRange:inputValues', 'uRange(1) must be less than or equal to uRange(2).');
            end
            % Order of assignment in next 3 lines is important.
            obj.uMin = p.Results.uMin;
            obj.uMax = p.Results.uMax;
            obj.nLev = p.Results.nLev; 
            obj.offset = mod(p.Results.uMax, obj.delta);
            obj.roundHalfTo = p.Results.roundHalfTo;
        end
        function obj = set.uRange(obj, uRange)
            % SET.URANGE Quantizer input range setter.
            %
            %   Input:
            %       uRange: 2-element numeric vector [uMin uMax].
            %
            %   See also: GET.URANGE
            %
            %   $Author: BH $    $Date: 2017-11-22 $    $Revision: 1 $
            if uRange(1)>uRange(2)
                error(['Minimum quantizer input (uMin = ' num2str(uRange(1)) ') is greater than maximum quantizer input (uMax = ' num2str(uRange(2)) ').']);
            end
            obj.uRangeCheck(uRange, obj.delta, false);
            obj.uMin = uRange(1);
            obj.uMax = uRange(2);
            obj.offset = mod(uRange(2), obj.delta);
        end
        function uRange = get.uRange(obj)
            % GET.URANGE Quantizer input range getter.
            %
            %   Output:
            %       uRange: 2-element numeric vector [uMin uMax].
            %
            %   See also: SET.URANGE
            %
            %   $Author: BH $    $Date: 2017-11-22 $    $Revision: 1 $
            uRange = [obj.uMin obj.uMax];
            obj.uRangeCheck(uRange, obj.delta, false);
        end
        function obj = set.nLev(obj, nLev)
            % SET.NLEV Number of quantization levels setter.
            %   Sets the superclass quantization spacing delta to divide
            %   the input range into nLev uniform quantization levels.
            %
            %   Input:
            %       nLev: numeric scalar number of quantization levels.
            %
            %   See also: GET.URANGE
            %
            %   $Author: BH $    $Date: 2017-11-22 $    $Revision: 1 $
            obj.nLevCheck(nLev, false);
            obj.delta = (obj.uMax - obj.uMin)/(nLev - 1);
        end
        function nLev = get.nLev(obj)
            % GET.NLEV Number of quantization levels getter.
            %   Gets the number of quantization levels from the superclass 
            %   quantization spacing delta.
            %
            %   Output:
            %       nLev: numeric scalar number of quantization levels.
            %
            %   See also: SET.URANGE
            %
            %   $Author: BH $    $Date: 2017-11-22 $    $Revision: 1 $
            nLev = (obj.uMax - obj.uMin)/obj.delta + 1;
            obj.nLevCheck(nLev, false);
        end
        function obj = set.bits(obj, bits)
            % SET.BITS Quantizer bits getter.
            %   Gets the number of bits from the number of quantization 
            %   levels. 
            %   
            %   Inputs:
            %       bits: numeric scalar number of quantizer bits.
            %
            %   See also: GET.BITS
            %
            %   $Author: BH $    $Date: 2017-11-22 $    $Revision: 1 $
            obj.nLev = bits + 1;
        end
        function bits = get.bits(obj)
            % SET.BITS Quantizer bits setter.
            %   Sets the number of quantization levels from the number of
            %   quantizer bits. 
            %   
            %   Output:
            %       bits: numeric scalar number of quantizer bits.
            %
            %   See also: SET.BITS
            %
            %   $Author: BH $    $Date: 2017-11-22 $    $Revision: 1 $
            bits = obj.nLev - 1;
        end
        function levels = get.levels(obj)
            % GET.BITS Quantization levels getter.
            %   
            %   Output:
            %       levels: numeric vector of quantization output levels.
            %
            %   $Author: BH $    $Date: 2017-11-22 $    $Revision: 1 $
            levels = obj.uMin:obj.delta:obj.uMax;
        end
        function [y, sat] = quantize(obj, u)
            %QUANTIZE Quantize the input.
            %
            %   Calls the superclass quantization function and saturates
            %   the result.
            %
            %   Input:
            %       u: numeric vector input.
            %   Output:
            %       y: quantized numeric vector output.
            %       sat: numeric vector indicating if no (sat=0), positive 
            %       (sat=1), or negative (sat=-1) saturation has occurred. 
            %
            %   See also: SDMUINFQUANTZER.QUANTIZE
            %
            %   $Author: BH $    $Date: 2017-11-16 $    $Revision: 1 $
            sat = zeros(size(u));
            y = quantize@SdmUInfQuantizer(obj, u);
            % Saturation.
            sat(u > obj.uMax) = 1;
            y(u > obj.uMax) = obj.uMax;
            sat(u < obj.uMin) = -1;
            y(u < obj.uMin) = obj.uMin;
        end
        function h = plot(obj, varargin)
            %PLOT Plot function overload for SdmUQuantizer class.
            %
            %   Input:
            %       varargin: additional arguments to native plot()
            %       function.
            %   Output:
            %       h: vector of handles to the plot object.
            %
            %   See also: SDMUINFQUANTIZER.PLOT, PLOT
            %
            %   $Author: BH $    $Date: 2017-11-15 $    $Revision: 1 $
            h = plot@SdmUInfQuantizer(obj, obj.uRange, varargin{:});
        end
    end
    methods (Access = protected)
        function isInt = nLevCheck(obj, nLev, isWarningHidden)
            %NLEVCHECK Check that number of quantization levels is an
            % integer.
            %
            %   Input:
            %       nLev: numeric scalar number of quantization levels.
            %       isWarningHidden: logical scalar true if warning is to
            %       be suppressed if nLev is not an integer.
            %   Output:
            %       isInt: logical scalar true if nLev is within obj.intTol
            %       of an integer value.
            %
            %   See also: URANGECHECK
            %
            %   $Author: BH $    $Date: 2017-11-15 $    $Revision: 1 $
            isInt = abs(mod(nLev, 1)) < obj.intTol;
            if ~isInt && ~isWarningHidden
                warning(['Number of quantization levels is not an integer (nLev = ' num2str(nLev) ').'])
            end
        end
        function isInt = uRangeCheck(obj, range, delta, isWarningHidden)
            %URANGECHECK Check that quantization range results in an
            % integer number of quantization levels given quantization
            % spacing delta.
            %
            %   Input:
            %       range: 2-element numeric vector [min max]. 
            %       isWarningHidden: logical scalar true if warning is to
            %       be suppressed if nLev is not an integer.
            %   Output:
            %       isInt: logical scalar true if quantization interval
            %       divides by delta within obj.intTol of an integer value. 
            %
            %   See also: NLEVCHECK
            %
            %   $Author: BH $    $Date: 2017-11-15 $    $Revision: 1 $
            isInt = abs(mod((range(2) - range(1))/delta, 1)) < obj.intTol;
            if ~isInt && ~isWarningHidden
                warning(['Quantizer range (uRange = ' num2str(range(2) - range(1)) ') is not evenly divisible by quantization spacing (delta = ' num2str(delta) ').'])
            end
        end
    end
end