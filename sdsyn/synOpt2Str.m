function str = synOpt2Str(synOpt, isHeaderOn)
%SYNOPT2STR Generate String from SynOpt Structure.
%   Produces a string which can be printed to the console output describing
%   the contents of a synOpt object.
%
%   Inputs:
%       synOpt: optimization target structure generated with
%           DEFINESYNOPT.M (required).
%       
%   Outputs:
%       str: a character array of a table of the contents of synOpt.
%   
%   See also: DEFINESYNOPT, DTSYN
%
%   $Author: BH $    $Date: 2018-02-20 $    $Revision: 0 $

    if isstruct(synOpt)
        str = '\n';
        if isHeaderOn
            str = [str '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n| synOpt: Summary of Optimization Targets                                   |\n'];
        end
        str = [str '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n| Target | Constraint | Norm    | Input ch.(s) | Output ch.(s) | Opt. Value |\n-----------------------------------------------------------------------------\n'];
        for iSynOpt=1:length(synOpt)
            str = [str sprintf('| %-6s ', num2str(iSynOpt))]; %#ok<*AGROW>
            str = [str sprintf('| %+10s ', num2str(synOpt(iSynOpt).constraint, 10))];
            switch synOpt(iSynOpt).norm
                case 1
                    str = [str '| l-1/(*) '];
                case 2
                    str = [str '| H-2     '];
                case Inf
                    if ~all(synOpt(iSynOpt).ffi==[0 pi])
                        str = [str '| GKYP    '];
                    else
                        str = [str '| H-Inf   '];
                    end
            end
            str = [str sprintf('| %+12s ', mat2str(synOpt(iSynOpt).inputChannel, 1))];
            str = [str sprintf('| %+13s ', mat2str(synOpt(iSynOpt).outputChannel, 1))];
            if isfield(synOpt, 'normz')
                normz = synOpt(iSynOpt).normz;
            else
                normz = NaN;
            end
            str = [str sprintf('| %+10s |\n', mat2str(normz, 4))];
        end
        str = [str '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'];
    end
end