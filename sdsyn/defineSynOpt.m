function synOpt = defineSynOpt(varargin)
%DEFINESYNOPT Generate Controller Synthesis Options Structure
%   Produces a structure array with synthesis options readable by the 
%   function DTSYN. Specify optimization goals one at a time in the order
%   shown below.
%
%   Inputs:
%       constraint: a numeric scalar defining the target sub-optimal value 
%           of the specified norm (L1MAX, H2MAX, or HINFMAX). If 
%           constraint==Inf, the specified norm is not optimized. If 
%           constranit==-1, the optimal value of the specified norm is 
%           computed.
%       norm: a numeric scalar equal to 1 for l-1/(*) optimization, 2 for 
%           H-2 optimization, or Inf for H-Inf optimization.
%       ffi: a 2-element vector defining the finite frequency interval (for
%           H-Inf GKYP condition), otherwise [0 pi];
%       inputChannel: a numeric vector indicating which generalized plant
%           disturbance input channel(s) are under the specified 
%           constraint (W1, W2, or WINF).
%       outputChannel: a numeric vector indicating which generalized plant
%           performance output channel(s) are under the specified 
%           constraint (Z1, Z2, or ZINF).
%       
%   Outputs:
%       synOpt: structure array to feed into synthesis function.
%   
%   See also: SYNOPT2STR, DTSYN
%
%   $Author: BH $    $Date: 2018-02-19 $    $Revision: 3 $
%
% REVISION 1:
%   2018-02-20 by BH: Removed console output code and placed into
%   SYNOPT2STR so it can be shared with the DTSYN function.
%
% REVISION 2:
%   2018-03-26 by BH: Added 'ffi' input to define the finite frequency
%   interval for H-Inf GKYP optimization.
%
% REVISION 3:
%   2018-04-04 by BH: Corrected error where the unrestricted frequency
%   range was specified as [0, 2*pi] rather than [0, pi].

    %% Validate function inputs.
    if rem(nargin, 5)~=0
        error('defineSynOpt:nInputGroups', '5 arguments per constraint must be provided.')
    end
    nInputGroup = nargin/4;
    fInput1Val = @(x) isnumeric(x) && isscalar(x);
    fInput2Val = @(x) x==1 || x==2 || x==Inf;
    fInput3Val = @(x) isnumeric(x) && isvector(x) && length(x)==2;
    fInput45Val = @(x) isnumeric(x) && isvector(x);
    synOpt = struct;
    
    %% Loop through input groups and assign to structure.
    for iInputGroup=1:nInputGroup
        input1 = varargin{5*iInputGroup-4};
        if fInput1Val(input1)
            synOpt(iInputGroup).constraint = input1;
        else
            error('defineSynOpt:input1', ['constraint(' num2str(5*iInputGroup-4) ') must be a numeric scalar.']);
        end
        input2 = varargin{5*iInputGroup-3};
        if fInput2Val(input2)
            synOpt(iInputGroup).norm = input2;
        else
            error('defineSynOpt:input2', ['norm(' num2str(5*iInputGroup-3) ') must be one of {1, 2, Inf}.']);
        end
        input3 = varargin{5*iInputGroup-2};
        if fInput3Val(input3)
            if (synOpt(iInputGroup).norm ~= Inf) && ~(all(input3 == [0 pi]) || all(input3 == [0 Inf]))
                error('defineSynOpt:input3', ['norm(' num2str(5*iInputGroup-3) ') must be Inf if a finite frequency interval is specified.']);
            end
            synOpt(iInputGroup).ffi = input3;
        else
            synOpt(iInputGroup).ffi = [0 pi];
        end
        input4 = varargin{5*iInputGroup-1};
        if fInput45Val(input4)
            synOpt(iInputGroup).inputChannel = input4;
        else
            error('defineSynOpt:input4', ['inputChannel(' num2str(5*iInputGroup-1) ') must be a numeric vector.'])
        end
        input5 = varargin{5*iInputGroup};
        if fInput45Val(input5)
            synOpt(iInputGroup).outputChannel = input5;
        else
            error('defineSynOpt:input5', ['outputChannel(' num2str(5*iInputGroup) ') must be a numeric vector.']);
        end
    end
    fprintf(1, synOpt2Str(synOpt, true));
end