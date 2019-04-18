classdef params_PF < handle

    properties
        Np
        model_stdev
        meas_stdev
        init_stdev
        bound_stdev
        meas_pt
    end
    methods
        function a = getMeasStdev(obj)
            a = obj.meas_stdev;
        end
    end
end

