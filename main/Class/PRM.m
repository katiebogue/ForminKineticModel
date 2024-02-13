classdef PRM < handle & dynamicprops
    %PRM Summary of this class goes here
    %   Detailed explanation goes here

    properties(SetAccess={?Formin})
        formin Formin
        dist_FH2 double % distance from center of PRM to FH2
        dist_NT double % distance from center of PRM to NT
        size double % number of prolines
        dist_FH2_start double % distance from CT P of PRM to FH2
        stat_props struct=struct % contains meta.DynamicProperty objects
    end

    properties(Dependent)
        gating double % gating factor
        c_PA double % concentration of profilin-actin | μM
        kpoly FilType % rate of polymerization by FilType
        fh1length double
        FH2dist_frac double % fractional distance from center of PRM to FH2
        
        % Parameters all by filament and FilType:
        kcap FilType % rate of PRM + profilin-actin binding (capture)| s^(-1)
        kdel FilType % rate of barbed end + PRM-profilin-actin binding (delivery) | s^(-1)
        rcap FilType % rate of PRM + profilin-actin dissociation (reverse capture) | s^(-1)
        rdel FilType % rate of barbed end + PRM-profilin-actin dissociation (reverse delivery) | s^(-1)
        krel FilType % rate of PRM + profilin-actin-barbed end dissociation (release) | s^(-1)
    end

    methods
        function obj = PRM(formin, dist_FH2, dist_NT, size, dist_FH2_start)
            %UNTITLED6 Construct an instance of this class
            %   Detailed explanation goes here
            if nargin>0
                obj.formin=formin;
                obj.dist_FH2=dist_FH2;
                obj.dist_NT=dist_NT;
                obj.size=size;
                obj.dist_FH2_start=dist_FH2_start;

                statlist=obj.formin.lookup.StatNames;
                for i=1:length(statlist)
                    obj.addStat(statlist(i));
                end
                addlistener(obj.formin.lookup,'StatNames','PostSet',@obj.updateStats);
            end
        end

        function value=get.gating(obj)
            value=obj.formin.gating;
        end

        function value=get.c_PA(obj)
            value=obj.formin.c_PA;
        end

        function value=get.kcap(obj)
            value=obj.formin.opts.k_cap*obj.formin.opts.equations.kcap(obj);
        end

        function value=get.kdel(obj)
            value=obj.formin.opts.k_del*obj.formin.opts.equations.kdel(obj);
        end

        function value=get.rcap(obj)
            value=obj.formin.opts.r_cap*obj.formin.opts.equations.rcap(obj);
        end

        function value=get.rdel(obj)
            value=obj.formin.opts.r_del*obj.formin.opts.equations.rdel(obj);
        end

        function value=get.krel(obj)
            value=obj.formin.opts.k_rel*obj.formin.opts.equations.krel(obj);
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

        function value=get.fh1length(obj)
            value=obj.formin.length;
        end

        function value=get.FH2dist_frac(obj)
            value=obj.dist_FH2./obj.formin.length;
        end

        function obj = addStat(obj,stat)
            if ~isprop(obj,stat)
                prop=obj.addprop(stat);
                prop.Dependent=true;
                prop.GetMethod = @(obj) obj.formin.lookup.getstat(iSite=obj.dist_FH2,N=obj.formin.length,Stat=stat,NName=obj.formin.NName);
                obj.stat_props.(stat)=prop;
            end
            obj.addBaseStat(stat);
        end

        function obj = addBaseStat(obj,stat)
            if ~isprop(obj,strcat(stat,"_Base"))
                prop=obj.addprop(strcat(stat,"_Base"));
                prop.Dependent=true;
                prop.GetMethod = @(obj) obj.formin.lookup.getstat(iSite=1,N=obj.formin.length,Stat=stat,NName=obj.formin.NName);
                obj.stat_props.(strcat(stat,"_Base"))=prop;
            end
        end

        function obj = updateStats(obj,src,evnt)
           fields=fieldnames(obj.stat_props);
           for i=1:length(fields)
               delete(obj.stat_props.(fields{i}))
           end
           obj.stat_props=struct;
           statlist=obj.formin.lookup.StatNames;
           for i=1:length(statlist)
               obj.addStat(statlist(i));
           end
        end

        function names=getfields(obj)
            if length(obj)==1
                names=fieldnames(obj);
            else
                fields1=fieldnames(obj);
                fields2=fieldnames([obj.stat_props]);
                names=[fields1;fields2];
            end
        end

        function out=getprop(obj,prop)
            if length(obj)==1
                out=obj.(prop);
            else
                out=[];
                for i=1:length(obj)
                    out=[out,obj(i).(prop)];
                end
            end
        end

    end

    methods (Static)
        function obj = loadobj(s)
         obj=s;
         obj.updateStats;
         addlistener(obj.formin.lookup,'StatNames','PostSet',@obj.updateStats);
      end
    end
    

end