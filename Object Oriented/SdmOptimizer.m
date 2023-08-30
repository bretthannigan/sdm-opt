classdef SdmOptimizer
    properties (Access = public)
        %% Weighting functions.
        W_r @lti % input weighting filter.
        W_s @lti % sensitivity weighting function.
        W_u @lti % control signal weighting function.
        W_t @lti % complementary sensitivity weighting function.
        %% Modulator.
        sdm @Sdm
    end
    properties (Dependent = true)
        
    end
    properties (Access = private)
    end
    methods
    end
end