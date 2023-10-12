classdef Experiment 
    %EXPERIMENT Summary of this class goes here
    %   Detailed explanation goes here

    properties
        ForminList (:,:) Formin
        opts Options
        params Params
        data (:,:) struct
    end

    methods
        function obj = Experiment(options,params,file,type,cPA)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                options Options
                params Params
                file string % .txt file, comma separated
                type string {mustBeMember(type,{'seq','uniprot'})}
                cPA double =2.5 % [profilin-actin] to assign to all formins
            end
            if nargin>0
                obj.opts=options;
                obj.params=params;
                forminlist = char(importdata(file)); 
                forminlist = strsplit(forminlist);
                obj.ForminList=Formin.empty((length(forminlist)/2),0);
                for i = 1:length(forminlist)/2
                    forminname = convertCharsToStrings(forminlist(2*i -1));   %takes the names (every other string)
                    input = convertCharsToStrings(forminlist(2*i));  
                    if type=="seq"
                        tempFormin=Formin(forminname,obj.params,obj.opts,c_PA=cPA,sequence=input);
                    elseif type=="uniprot"
                        tempFormin=Formin(forminname,obj.params,obj.opts,c_PA=cPA,uniprotID=input);
                    end
                    obj.ForminList(i)=tempFormin;
                end
            end
        end

        function obj=add_formin(obj,formin)
            arguments
                obj Experiment
                formin Formin
            end
            obj.ForminList(1,end+1)=formin;
        end

        function obj= set.opts(obj,input)
            obj.opts=input;
            obj.syncforminsettings;
        end

        function obj= set.params(obj,input)
            obj.params=input;
            obj.syncforminsettings;
        end

        function syncforminsettings(obj)
            if isfield(obj,"ForminList")
                for i=1:length(obj.ForminList)
                    if obj.ForminList(i).opts~=obj.opts
                        obj.ForminList(i).opts=obj.opts;
                    end
                    if obj.ForminList(i).params~=obj.params
                        obj.ForminList(i).params=obj.params;
                    end
                end
            end
        end

        function T=formintable(obj,save)
            arguments
                obj Experiment
                save logical=false
            end
            fh1names=[obj.ForminList.name];
            T= table('RowNames',fh1names);
            fields=fieldnames([obj.ForminList]);
            for i=1:length(fields)
                varname=fields{i};
                val=[obj.ForminList.(varname)];
                if class(val)=="FilType"
                    T.(strcat("Single ",varname))=[val.single]';
                    T.(strcat("Double ",varname))=[val.double]';
                    T.(strcat("Dimer ",varname))=[val.dimer]';
                    T.(strcat("Dimer/Double ",varname))=[val.ratio]';
                elseif class(val)=="double" && ~isempty(val) && length(val)==length(fh1names)
                    T.(varname)=[val]';
                elseif class(val)=="string"
                    T.(varname)=[val]';
                end
            end
            
            
            if save
                writetable(T,'Per Formin Data.csv')
            end

        end

        function T=PRMtable(obj,save)
            arguments
                obj Experiment
                save logical=false
            end
            PRMs=[obj.ForminList.PRMList];
            formins=[PRMs.formin];
            fh1names=[formins.name];
            T= table;
            T.("formin")=[fh1names]';
            fields=PRMs.getfields;
            for i=1:length(fields)
                varname=fields{i};
                val=[PRMs.getprop((varname))];
                if class(val)=="FilType" 
                    T.(strcat("Single ",varname))=[val.single]';
                    double=[val.double];
                    dimer=[val.dimer];
                    ratio=[val.ratio];
                    if class(double)=="Filament" && length([double.a])==length(fh1names)
                        T.(strcat("Double a ",varname))=[double.a]';
                        T.(strcat("Double a ",varname))=[double.b]';
                    else
                        T.(strcat("Double ",varname))=[double]';
                    end
                    if class(dimer)=="Filament" && length([dimer.a])==length(fh1names)
                        T.(strcat("Dimer a ",varname))=[dimer.a]';
                        T.(strcat("Dimer a ",varname))=[dimer.b]';
                    else
                        T.(strcat("Dimer ",varname))=[dimer]';
                    end
                    if class(ratio)=="Filament" && length([ratio.a])==length(fh1names)
                        T.(strcat("Dimer/Double a ",varname))=[ratio.a]';
                        T.(strcat("Dimer/Double a ",varname))=[ratio.b]';
                    else
                        T.(strcat("Dimer/Double ",varname))=[ratio]';
                    end
                elseif class(val)=="double" && ~isempty(val) && length(val)==length(fh1names)
                    T.(varname)=[val]';
                elseif class(val)=="string"
                    T.(varname)=[val]';
                end
            end
            
            
            if save
                writetable(T,'Per PRM Data.csv')
            end

        end

        function obj = add_data(obj,forminname,value,type,NameValueArgs)
            arguments
                obj Experiment
                forminname string 
                value double
                type string {mustBeMember(type,{'single','double','dimer','ratio'})}
                NameValueArgs.errplus double
                NameValueArgs.errminus double
                NameValueArgs.errtop double
                NameValueArgs.errbot double
                NameValueArgs.errperc double {mustBeLessThanOrEqual(NameValueArgs.errperc,1)}
            end
            formin=-1;
            for i=1:length(obj.ForminList)
                if obj.ForminList(i).name==forminname
                    if formin==-1
                        formin=obj.ForminList(i);
                    else
                        error("multiple formins in this Experiment have the name %s",forminname)
                    end
                end
            end
            if formin==-1
                error("no formin in this Experiment has the name %s",forminname)
            end

            datastruct.type=type;
            datastruct.formin=formin;
            datastruct.value=value;

            if length(fieldnames(NameValueArgs))>2
                error("can only enter one type of error")
            elseif isfield(NameValueArgs,"errplus") && isfield(NameValueArgs,"errtop")
                 error("cannot enter two types of plus error")
            elseif isfield(NameValueArgs,"errminus") && isfield(NameValueArgs,"errbot")
                error("cannot enter two types of minus error")
            end


            errtop=0;
            errbot=0;
            if isfield(NameValueArgs,"errperc")
                if length(fieldnames(NameValueArgs))>1
                    error("can only enter one type of error")
                end
                errtop=value+(value*NameValueArgs.errperc);
                errbot=value-(value*NameValueArgs.errperc);
            else
                if isfield(NameValueArgs,"errplus")
                    errtop=value+NameValueArgs.errplus;
                end
                if isfield(NameValueArgs,"errtop")
                    if NameValueArgs.errtop<value
                        error("top error cannot be greater than the value")
                    end
                    errtop=NameValueArgs.errtop;
                end
                if isfield(NameValueArgs,"errminus")
                    errbot=value-NameValueArgs.errminus;
                end
                if isfield(NameValueArgs,"errbot")
                    if NameValueArgs.errbot>value
                        error("bottom error cannot be less than the value")
                    end
                    errbot=NameValueArgs.errbot;
                end
            end

            datastruct.errtop=errtop;
            datastruct.errbot=errbot;

            if isempty(obj.data)
                obj.data=[datastruct];
            else
                obj.data(1,1+end)=datastruct;
            end

        end

        function tab=optionstable(obj)
            
            tabstruct=struct;

            addel(obj.params);
            addel(obj.opts);
            addel(obj.opts.equationstext);

            Rowtitles=fieldnames(tabstruct);
            Datavars=struct2cell(tabstruct);
            tab=cell2table([Datavars]);
            tab.Properties.RowNames=Rowtitles;
            tab.Properties.VariableNames={' '};

            function addel(object)
                optfields=fieldnames(object);
                for i=1:length(optfields)
                    if class(object.(optfields{i}))=="double" 
                        tabstruct.(optfields{i})=num2str(object.(optfields{i}));
                    elseif class(object.(optfields{i}))=="string" && length(object.(optfields{i}))==1
                        tabstruct.(optfields{i})=object.(optfields{i});
                    end
                end
            end
        end
    end
end