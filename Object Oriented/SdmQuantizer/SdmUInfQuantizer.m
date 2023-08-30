classdef SdmUInfQuantizer < SdmQuantizerSystem
    %SDMUINFQUANTIZER Uniform infinite quantizer subclass.
    %   Implements a quantizer with uniform spacing and infinite range by
    %   applying an affine transformation on the input to and output from
    %   the SdmQuantizer base class.
    % SdmQuantizer Family Class Inheritance Hierarchy:
    %    
    % Base class:                           SdmQuantizer
    %                      ,---------------------|-----------------.
    % Infinite:    SdmUInfQuantizer(this) SdmLInfQuantizer         |
    %                      |                     |                 |
    % Saturating:    SdmUQuantizer         SdmLQuantizer      SdmAQuantizer
    %                      
    % Spacing:        ^Uniform^            ^Logarithmic^       ^Arbirtary^
    %
    %   $Author: BH $    $Date: 2017-11-16 14:20:00 $    $Revision: 1 $
    properties (Access = protected)
        T = eye(3); % Affine transformation matrix.
    end
    properties (Dependent = true)
        offset % Numeric scalar indicating translation of base quantizer offset units right and offset units down.
        delta % Numeric scalar quantization spacing.
    end
    methods
        function obj = SdmUInfQuantizer(varargin)
            %SDMUINFQUANTIZER SdmUInfQuantizer object constructor.
            %   Constructs an SdmUInfQuantizer object.
            %
            %   Inputs:
            %       delta: numeric scalar quantization level spacing
            %       (optional, default:1).
            %       offset: numeric scalar horizontal and vertical offset
            %       (optional, default:0).
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
            %       obj: An SdmUInfQuantizer object.
            %
            %   See also: SDMQUANTIZER
            %
            %   $Author: BH $    $Date: 2017-11-15 $    $Revision: 1 $
            p = inputParser;
            addOptional(p, 'delta', 1, @isscalar);
            addOptional(p, 'offset', 0, @isscalar);
            addParameter(p, 'RoundHalfTo', 1, @isscalar);
            parse(p, varargin{:});
            obj.delta = p.Results.delta;
            obj.offset = p.Results.offset;
            obj.roundHalfTo = p.Results.RoundHalfTo;
        end
    end
    methods
        function y = quantize(obj, u)
            %QUANTIZE Quantize the input.
            %
            %   Applies the affine transformation to the unit quantizer,
            %   calls the base class quantization function, and applies the
            %   inverse transformation on the result.
            %
            %   Input:
            %       u: numeric vector input.
            %   Output:
            %       y: quantized numeric vector output.
            %
            %   See also: SDMQUANTZER.QUANTIZE
            %
            %   $Author: BH $    $Date: 2017-11-14 $    $Revision: 1 $
            if ~isrow(u)
                u = u';
                isTransposed = true;
            else
                isTransposed = false;
            end
            [u, ~] = obj.itransform(u, zeros(size(u)));
            y = quantize@SdmQuantizerSystem(obj, u);
            [~, y] = obj.transform(zeros(size(y)), y);
            if isTransposed, y = y'; end
        end
        function obj = plus(obj, b)
            obj.T(2,3) = obj.T(2,3) + b;
        end
        function obj = minus(obj, b)
            obj.T(2,3) = obj.T(2,3) - b;
        end
        function obj = set.delta(obj, delta)
            % SET.DELTA Quantization spacing setter.
            %
            %   Input:
            %       delta: numeric scalar quantization spacing.
            %
            %   See also: GET.DELTA
            %
            %   $Author: BH $    $Date: 2017-11-20 $    $Revision: 1 $
            obj.T(1,1) = delta;
            obj.T(2,2) = delta;
        end
        function delta = get.delta(obj)
            % SET.DELTA Quantization spacing getter.
            %
            %   Output:
            %       delta: numeric scalar quantization spacing.
            %
            %   See also: SET.DELTA
            %
            %   $Author: BH $    $Date: 2017-11-20 $    $Revision: 1 $
            delta = obj.T(2,2);
        end
        function obj = set.offset(obj, offset)
            % SET.OFFSET Quantizer offset setter.
            %
            %   Input:
            %       offset: quantizer offset from zero.
            %
            %   See also: GET.OFFSET
            %
            %   $Author: BH $    $Date: 2017-11-20 $    $Revision: 1 $
            obj.T(1:2,3) = [-offset; -offset];
        end
        function offset = get.offset(obj)
            % GET.OFFSET Quantizer offset getter.
            %
            %   Output:
            %       offset: quantizer offset from zero.
            %
            %   See also: SET.OFFSET
            %
            %   $Author: BH $    $Date: 2017-11-20 $    $Revision: 1 $
            offset = obj.T(1,3);
        end
        function h = plot(obj, uLim, varargin)
            %PLOT Plot function overload for SdmUInfQuantizer class
            %
            %   Input:
            %       uLim: x-axis limits.
            %       varargin: additional arguments to native plot()
            %       function.
            %   Output:
            %       h: vector of handles to the plot objects.
            %
            %   See also: PLOT
            %
            %   $Author: BH $    $Date: 2017-11-15 $    $Revision: 2 $
            %
            % REVISION 2:
            %   2017-11-16 by BH: Modified for operation with
            %   transformation matrices.
            %% Generate points.
            [uLim, ~] = obj.itransform(uLim);
            midPts = floor(uLim(1)):ceil(uLim(2)); % Midpoints of quantization intervals.
            splitPts = (midPts(1)+0.5):(midPts(end)-0.5); % Splits between quantization intervals.
            closedPts = [obj.transform(splitPts)' obj.quantize(obj.transform(splitPts))']; % Included points.
            allPts = [obj.transform([midPts+0.5 midPts-0.5])' repmat(obj.quantize(obj.transform(midPts))', 2, 1)];
            openPts = setdiff(allPts, closedPts, 'rows'); % Excluded points.
            %% Produce plots.
            h1 = plot(closedPts(:,1), closedPts(:,2), varargin{:});
            hold on
            h2 = plot([allPts(1:end/2,1) allPts(end/2+1:end,1)]', [allPts(1:end/2,2) allPts(end/2+1:end,2)]', varargin{:});
            h3 = plot(openPts(:,1), openPts(:,2), varargin{:});
            %% Customize plots.
            set(h1, 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', get(h1, 'Color'));
            set(h2, 'LineStyle', '-', 'Color', get(h1, 'Color'));
            set(h3, 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'None', 'MarkerEdgeColor', get(h1, 'Color'));
            title('Quantizer Transfer Curve')
            xLabelString = 'From:';
            if ~isempty(obj.inputName)
                xLabelString = [xLabelString ' ' obj.inputName];
            else
                xLabelString = [xLabelString ' In'];
            end
            if ~isempty(obj.inputUnit)
                xLabelString = [xLabelString ' (' obj.inputUnit ')'];
            end
            xlabel(xLabelString)
            yLabelString = 'To:';
            if ~isempty(obj.outputName)
                yLabelString = [yLabelString ' ' obj.outputName];
            else
                yLabelString = [yLabelString ' Out'];
            end
            if ~isempty(obj.outputUnit)
                yLabelString = [yLabelString ' (' obj.outputUnit ')'];
            end
            ylabel(yLabelString)
            xlim(obj.transform(uLim))
            ylim(obj.quantize(obj.transform(uLim)))
            h = [h1; h2; h3];
        end
    end
    methods (Access = protected)
        function [x_t, y_t] = transform(obj, varargin)
            %TRANSFORM Apply affine transformation.
            %
            %   Uses the affine transformation matrix to convert the
            %   input/output of the unit quantizer in the SdmQuantizer base
            %   class to the input/output of the uniform quantizer.
            %
            %   Input:
            %       x: numeric vector x-coordinate input for unit
            %       quantizer.
            %       y: (optional) numeric vector y-coordinate output for
            %       unit quantizer.
            %   Output:
            %       x_t: transformed numeric vector x-coordinate input for
            %       uniform quantizer.
            %       y_t: transformed numeric vector y-coordinate output for
            %       uniform quantizer.
            %
            %   See also: ITRANSFORM
            %
            %   $Author: BH $    $Date: 2017-11-14 $    $Revision: 1 $
            p = inputParser;
            addRequired(p, 'x', @(x) isnumeric(x) && isvector(x));
            parse(p, varargin{1});
            x = p.Results.x;
            addOptional(p, 'y', zeros(size(x)), @(x) isnumeric(x) && isvector(x));
            parse(p, varargin{:});
            y = p.Results.y;
            if ~isrow(x), x = x'; end
            if ~isrow(y), y = y'; end
            if length(x)<length(y)
                x = [x zeros(1, length(y) - length(x))];
            elseif length(y)<length(x)
                y = [y zeros(1, length(x) - length(y))];
            end
            xy_t = obj.T*[x; y; ones(size(x))];
            x_t = xy_t(1,:);
            y_t = xy_t(2,:);
        end
        function [x, y] = itransform(obj, varargin)
            %ITRANSFORM Apply inverse affine transformation.
            %
            %   Uses the affine transformation matrix to convert the
            %   input/output of the uniform quantizer to the input/output
            %   of the unit quantizer in the SdmQuantizer base class.
            %
            %   Input:
            %       x: transformed numeric vector x-coordinate input for
            %       uniform quantizer.
            %       y: (optional) transformed numeric vector y-coordinate
            %       output for uniform quantizer.
            %   Output:
            %       x_t: numeric vector x-coordinate input for unit 
            %       quantizer.
            %       y_t: numeric vector y-coordinate output for unit 
            %       quantizer.
            %
            %   See also: TRANSFORM
            %
            %   $Author: BH $    $Date: 2017-11-14 $    $Revision: 1 $
            p = inputParser;
            addRequired(p, 'x_t', @(x) isnumeric(x) && isvector(x));
            parse(p, varargin{1});
            x_t = p.Results.x_t;
            addOptional(p, 'y_t', zeros(size(x_t)), @(x) isnumeric(x) && isvector(x));
            parse(p, varargin{:});
            y_t = p.Results.y_t;
            if ~isrow(x_t), x_t = x_t'; end
            if ~isrow(y_t), y_t = y_t'; end
            if length(x_t)<length(y_t)
                x_t = [x_t zeros(1, length(y_t) - length(x_t))];
            elseif length(y_t)<length(x_t)
                y_t = [y_t zeros(1, length(x_t) - length(y_t))];
            end
            xy = obj.T\[x_t; y_t; ones(size(x_t))];
            x = xy(1,:);
            y = xy(2,:);
        end
    end
end