classdef PRM < handle
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
    end

    properties(Dependent)
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

        function value=get.pocc(obj)
            lookup=obj.formin.lookup;
            extrapolation=obj.formin.extrapolation;

            if extrapolation.single
                pocc_single=lookup.pocc_extrapolation.single(obj.formin.length,obj.dist_FH2);
            else
                pocc_single=lookup.poccs.single{obj.formin.length}(obj.dist_FH2);
            end
            
            function output=pocc(filament)
                if filament=="single"
                    if extrapolation.single
                        output=lookup.pocc_extrapolation.single(obj.formin.length,obj.dist_FH2);
                    else
                        output=lookup.poccs.single{obj.formin.length}(obj.dist_FH2);
                    end
                else
                    type={'a','b'};
                    if extrapolation.(filament)
                    else
                    end
                    for x=1:2
                        t=type{x};
                        
                    end
                end
            end
        end

        function value=get.pocc_0(obj)
            value=obj.formin.pocc_0;
        end
    end
end