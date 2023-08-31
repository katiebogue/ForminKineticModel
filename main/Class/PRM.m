classdef PRM < handle
    %PRM Summary of this class goes here
    %   Detailed explanation goes here

    properties
        formin Formin
        dist_FH2 double
        dist_NT double
        size double % number of prolines
    end

    properties(Dependent)
        kpoly Submodel % rate of polymerization by Submodel
        pocc Submodel
        pocc_0 Submodel
        pr Submodel
        % Parameters all by filament and submodel:
        kcap Submodel % rate of PRM + profilin-actin binding (capture)| s^(-1)
        kdel Submodel % rate of barbed end + PRM-profilin-actin binding (delivery) | s^(-1)
        rcap Submodel % rate of PRM + profilin-actin dissociation (reverse capture) | s^(-1)
        rdel Submodel % rate of barbed end + PRM-profilin-actin dissociation (reverse delivery) | s^(-1)
        krel Submodel % rate of PRM + profilin-actin-barbed end dissociation (release) | s^(-1)
    end

    methods
        function obj = PRM(inputArg1,inputArg2)
            %UNTITLED6 Construct an instance of this class
            %   Detailed explanation goes here
            if nargin>0
            obj.Property1 = inputArg1 + inputArg2;
                end
        end

        function value=get.pocc(obj)
            lookup=obj.formin.lookup;
            extrapolation=obj.formin.extrapolation;

            pocc_single=p_occ("single");
            pocc_double=p_occ("double");
            pocc_dimer=p_occ("dimer");

            value=Submodel(pocc_single,pocc_double,pocc_dimer);
            
            function output=p_occ(filament)
                if filament=="single"
                    if extrapolation.single
                        output=lookup.pocc_extrapolation.single(obj.formin.length,obj.dist_FH2);
                    else
                        output=lookup.poccs.single{obj.formin.length}(obj.dist_FH2);
                    end
                else
                    type={'a','b'};
                    output=Filament();
                    for x=1:2
                        t=type{x};
                        if extrapolation.(filament)
                            output.(t)=lookup.pocc_extrapolation.(filament).(t)(obj.formin.length,obj.dist_FH2);
                        else
                            output.(t)=lookup.poccs.(filament).(t)(obj.formin.length,obj.dist_FH2);
                        end
                    end
                end
            end
        end

        function value=get.pr(obj)
            lookup=obj.formin.lookup;
            extrapolation=obj.formin.extrapolation;

            pr_single=pr_fun("single");
            pr_double=pr_fun("double");
            pr_dimer=pr_fun("dimer");

            value=Submodel(pr_single,pr_double,pr_dimer);
            
            function output=pr_fun(filament)
                if filament=="single"
                    if extrapolation.single
                        output=lookup.pr_extrapolation.single(obj.formin.length,obj.dist_FH2);
                    else
                        output=lookup.prs.single{obj.formin.length}(obj.dist_FH2);
                    end
                else
                    type={'a','b'};
                    output=Filament();
                    for x=1:2
                        t=type{x};
                        if extrapolation.(filament)
                            output.(t)=lookup.pr_extrapolation.(filament).(t)(obj.formin.length,obj.dist_FH2);
                        else
                            output.(t)=lookup.prs.(filament).(t)(obj.formin.length,obj.dist_FH2);
                        end
                    end
                end
            end
        end

        function value=get.pocc_0(obj)
            value=obj.formin.pocc_0;
        end

        function value=get.kcap(obj)
            if isKey(obj.formin.opts.scaling, "kcap") 
                kcap_scale=obj.formin.params.k_cap*obj.formin.opts.scaling{"kcap"}(obj);
            else
                kcap_scale=obj.formin.params.k_cap;
            end
            kcap_single=kcapture(kcap_scale,obj.pocc.single,obj.formin.params.c_PA);
            kcap_doublea=kcapture(kcap_scale,obj.pocc.double.a,obj.formin.params.c_PA);
            kcap_doubleb=kcapture(kcap_scale,obj.pocc.double.b,obj.formin.params.c_PA);
            kcap_dimera=kcapture(kcap_scale,obj.pocc.dimer.a,obj.formin.params.c_PA);
            kcap_dimerb=kcapture(kcap_scale,obj.pocc.dimer.b,obj.formin.params.c_PA);
            value=Submodel(kcap_single,kcap_doublea,kcap_dimera,kcap_doubleb,kcap_dimerb);
        end

        function value=get.kdel(obj)
            if isKey(obj.formin.opts.scaling, "kdel") 
                kdel_scale=obj.formin.params.k_del*obj.formin.opts.scaling{"kdel"}(obj);
            else
                kdel_scale=obj.formin.params.k_del;
            end
            kdel_single=kdelivery(kdel_scale,obj.pocc_0.single,obj.pr.single);
            kdel_doublea=kdelivery(kdel_scale,obj.pocc_0.double.a,obj.pr.double.a);
            kdel_doubleb=kdelivery(kdel_scale,obj.pocc_0.double.b,obj.pr.double.b);
            kdel_dimera=kdelivery(kdel_scale,obj.pocc_0.dimer.a,obj.pr.dimer.a);
            kdel_dimerb=kdelivery(kdel_scale,obj.pocc_0.dimer.b,obj.pr.dimer.b);
            value=Submodel(kdel_single,kdel_doublea,kdel_dimera,kdel_doubleb,kdel_dimerb);
        end

        function value=get.rcap(obj)
            if isKey(obj.formin.opts.scaling, "rcap") 
                rcap_scale=obj.formin.params.r_cap*obj.formin.opts.scaling{"rcap"}(obj);
            else
                rcap_scale=obj.formin.params.r_cap;
            end
            value=Submodel(rcap_scale);
        end

        function value=get.rdel(obj)
            if isKey(obj.formin.opts.scaling, "rdel") 
                rdel_scale=obj.formin.params.r_del*obj.formin.opts.scaling{"rdel"}(obj);
            else
                rdel_scale=obj.formin.params.r_del;
            end
            value=Submodel(rdel_scale);
        end

        function value=get.krel(obj)
            if isKey(obj.formin.opts.scaling, "krel") 
                krel_scale=obj.formin.params.k_rel*obj.formin.opts.scaling{"krel"}(obj);
            else
                krel_scale=obj.formin.params.k_rel;
            end
            value=Submodel(krel_scale);
        end

        function value=get.kpoly(obj)
            kpoly_type = obj.formin.opts.kpoly_type;
            value=kpolymerization(kpoly_type,obj.kcap,obj.kdel,obj.rcap,obj.rdel,obj.krel);
        end

        function r=calculate_kpoly(obj,kpoly_type)
            arguments
                obj PRM
                kpoly_type string = obj.formin.opts.kpoly_type
            end
            r=kpolymerization(kpoly_type,obj.kcap,obj.kdel,obj.rcap,obj.rdel,obj.krel);
        end
    end
end