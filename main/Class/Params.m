classdef Params < handle
    %PARAMS Parameters for formin calculations
    %   construct with: Params(k_cap,k_del,r_cap,r_del,k_rel,c_PA)

    properties
        k_cap double % rate constant for PRM + profilin-actin binding (capture) | μM^(-1)s^(-1)
        k_del double % rate constant for barbed end + PRM-profilin-actin binding (delivery) | μM^(-1)s^(-1)
        r_cap double % rate constant for PRM + profilin-actin dissociation (reverse capture) | s^(-1)
        r_del double % rate constant for barbed end + PRM-profilin-actin dissociation (reverse delivery) | s^(-1)
        k_rel double % rate constant for PRM + profilin-actin-barbed end dissociation (release) | s^(-1)
        c_PA double % concentration of profilin-actin | μM
    end

    methods
        function obj = Params(k_cap,k_del,r_cap,r_del,k_rel,c_PA)
            %PARAMS Construct an instance of Params
            %   Assigns input values to all properties
            obj.k_cap = k_cap;
            obj.k_del = k_del;
            obj.r_cap = r_cap;
            obj.r_del = r_del;
            obj.k_rel = k_rel;
            obj.c_PA = c_PA;
        end

    end
end