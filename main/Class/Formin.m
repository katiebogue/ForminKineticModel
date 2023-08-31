classdef Formin < handle
    %FORMIN Summary of this class goes here
    %   Detailed explanation goes here

    properties
        name string
        length double
        PRMCount double
        PRMList (1,:) PRM
        params Params
        lookup Lookuptable
        opts Options
    end

    properties(Dependent)
        kcap Submodel
        kdel Submodel
        kpoly Submodel
        pocc_0 Submodel
        extrapolation Submodel
    end

    methods
        function obj = Formin(name,length)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
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

        function value= get.pocc_0(obj)
            p_occs=obj.lookup.poccs;
            p_occs_ext = obj.lookup.pocc_extrapolation;
            if obj.extrapolation.single
                psingle=p_occs_ext.single(obj.length,1);
            else
                psingle=p_occs.single{obj.length}(1);
            end
            
            if obj.extrapolation.double
                pdoublea=p_occs_ext.double.a(obj.length,1);
                pdoubleb=p_occs_ext.double.b(obj.length,1);
            else
                pdoublea=p_occs.double.a{obj.length}(1);
                pdoubleb=p_occs.double.b{obj.length}(1);
            end

            if obj.extrapolation.dimer
                pdimera=p_occs_ext.dimer.a(obj.length,1);
                pdimerb=p_occs_ext.dimer.b(obj.length,1);
            else
                pdimera=p_occs.dimer.a{obj.length}(1);
                pdimerb=p_occs.dimer.b{obj.length}(1);
            end
            
            value=Submodel(psingle,pdoublea,pdimera,pdoubleb,pdimerb);
        end

        function value=get.kcap(obj)
            value=get_rate(obj,"kcap");
        end

        function value=get.kdel(obj)
            value=get_rate(obj,"kdel");
        end

        function value=get.kpoly(obj)
            value=get_rate(obj,"kpoly");
        end

        function r=get_rate(obj,rate)
            rate_sum=Submodel;
            for i=1:obj.PRMCount
                rate_sum=rate_sum+obj.PRMList(1).(rate);
            end
            r=rate_sum.add_fils;
        end
    end
end