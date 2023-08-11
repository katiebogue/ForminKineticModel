classdef Formin < handle
    %FORMIN Summary of this class goes here
    %   Detailed explanation goes here

    properties
        name string
        length double
        PRMCount double
        PRMList (1,:) PRM 
        kcap double
        kdel double
        kpoly double
        params Params
        lookup Lookuptable
    end

    properties(Dependent)
        pocc_0 Submodel
        extrapolation Submodel
    end

    methods
        function obj = Formin(inputArg1,inputArg2)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end

        function value = get.extrapolation(obj)
            exsingle=isextrapolated('single');
            exdouble=isextrapolated('double');
            exdimer=isextrapolated('dimer');
            value=Submodel(exsingle,exdouble,exdimer);

            function out = isextrapolated(type)
                out = obj.length > obj.lookup.max_length.(type);
            end
        end
    end
end