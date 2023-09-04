
    properties
        max_length Submodel % Max filament length before extrapolation
        poccs Submodel
        prs Submodel
        pocc_extrapolation Submodel
        pr_extrapolation Submodel
    end

    properties(Dependent)
        pocc_0s Submodel
    end

    methods

        function obj = Lookuptable(max_single,max_double,max_dimer,file_single,file_double,file_dimer)
            %LOOKUPTABLE Construct an instance of Lookuptable
            %   
            if nargin>0
                obj.max_length = Submodel(max_single,max_double,max_dimer);
                obj.poccs=Submodel();
                obj.prs=Submodel();
                obj.pocc_extrapolation=Submodel();
                obj.pr_extrapolation=Submodel();
                set_fil(obj,file_single,"single");
                set_fil(obj,file_double,"double");
                set_fil(obj,file_dimer,"dimer");
            end
        end


        function value = get.pocc_0s(obj)
            x=obj.poccs;
            for ikey=keys(x.single)
                x.single(ikey)=x.single{ikey}(1);
            end
            fils={"double","dimer"};
            for i=1:2
                type=fils{i};
                for ikey=keys(x.(type))
                    x.(type).a(ikey)=x.(type).a{ikey}(1);
                    x.(type).b(ikey)=x.(type).a{ikey}(1);
                end
            end
            value=x;
        end

        function set_fil(obj,file,filament,maxL)
            if filament=="single"
                if nargin>3
                    obj.max_length.(filament)=maxL;
                end
                [obj.poccs.(filament), obj.prs.(filament),obj.pocc_extrapolation.(filament),obj.pr_extrapolation.(filament)]=readlookup(1,0,file,obj.max_length.(filament));
            else
                obj.poccs.(filament)=Filament();
                type={'a','b'};
                for x=1:2
                    t=type{x};
                    if nargin>3
                        obj.max_length.(filament).(t)=maxL;
                    end
                    [obj.poccs.(filament).(t), obj.prs.(filament).(t),obj.pocc_extrapolation.(filament).(t),obj.pr_extrapolation.(filament).(t)]=readlookup(2,(x-1),file,obj.max_length.(filament).(t));
                end
            end
        end
    end
end