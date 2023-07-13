classdef PRM
    %PRM Summary of this class goes here
    %   Detailed explanation goes here

    properties
        formin Formin
        kcap Submodel
        kdel Submodel
        kpoly Submodel
        dist_FH2 double
        dist_NT double
        size double
        pocc Submodel
        pocc_0 Submodel
        pr Submodel
    end

    methods
        function obj = untitled6(inputArg1,inputArg2)
            %UNTITLED6 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end