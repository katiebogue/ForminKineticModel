classdef Lookuptable <handle
    %LOOKUPTABLE Summary of this class goes here
    %   Detailed explanation goes here

    properties
        max_length Submodel % Max filament length before extrapolation
        poccs dictionary
        prs dictionary
    end

    properties(Dependent)
        pocc_0s dictionary
    end

    methods

        function obj = Lookuptable(max_single, max_double,max_dimer)
            %LOOKUPTABLE Construct an instance of Lookuptable
            %   
            if nargin>0
                obj.max_length = Submodel(max_single,max_double,max_dimer);
                obj.poccs
            end
        end


        function value = get.pocc_0s(obj)
           
           
        end
    end
end