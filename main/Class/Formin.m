classdef Formin
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
    end
end