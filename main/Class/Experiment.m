classdef Experiment 
    %EXPERIMENT Summary of this class goes here
    %   Detailed explanation goes here

    properties
        ForminList (:,:) Formin
        opts Options
        data (:,:) struct
    end

    methods
        function obj = Experiment(options,file,type,cPA)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                options Options
                file string % .txt file, comma separated
                type string {mustBeMember(type,{'seq','uniprot'})}
                cPA double =2.5 % [profilin-actin] to assign to all formins
            end
            if nargin>0
                obj.opts=options;
                forminlist = char(importdata(file)); 
                forminlist = strsplit(forminlist);
                obj.ForminList=Formin.empty((length(forminlist)/2),0);
                for i = 1:length(forminlist)/2
                    forminname = convertCharsToStrings(forminlist(2*i -1));   %takes the names (every other string)
                    input = convertCharsToStrings(forminlist(2*i));  
                    if type=="seq"
                        tempFormin=Formin(forminname,obj.opts,c_PA=cPA,sequence=input);
                    elseif type=="uniprot"
                        tempFormin=Formin(forminname,obj.opts,c_PA=cPA,uniprotID=input);
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

        function syncforminsettings(obj)
            if isfield(obj,"ForminList")
                for i=1:length(obj.ForminList)
                    if obj.ForminList(i).opts~=obj.opts
                        obj.ForminList(i).opts=obj.opts;
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
                if ~exist(fullfile(obj.opts.resultsdir,obj.opts.resultsfolder),'dir')
                    mkdir (obj.opts.resultsdir,obj.opts.resultsfolder)
                end
                writetable(T,fullfile(obj.opts.resultsdir,obj.opts.resultsfolder,"Per Formin Data.csv"))
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
                if ~exist(fullfile(obj.opts.resultsdir,obj.opts.resultsfolder),'dir')
                    mkdir (obj.opts.resultsdir,obj.opts.resultsfolder)
                end
                writetable(T,fullfile(obj.opts.resultsdir,obj.opts.resultsfolder,"Per PRM Data.csv"))
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
                NameValueArgs.groups string=""
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
            datastruct.groups=NameValueArgs.groups;
            NameValueArgs=rmfield(NameValueArgs,"groups");

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
                errtop=(value*NameValueArgs.errperc);
                errbot=(value*NameValueArgs.errperc);
            else
                if isfield(NameValueArgs,"errplus")
                    errtop=NameValueArgs.errplus;
                end
                if isfield(NameValueArgs,"errtop")
                    if NameValueArgs.errtop<value
                        error("top error cannot be greater than the value")
                    end
                    errtop=NameValueArgs.errtop-value;
                end
                if isfield(NameValueArgs,"errminus")
                    errbot=NameValueArgs.errminus;
                end
                if isfield(NameValueArgs,"errbot")
                    if NameValueArgs.errbot>value
                        error("bottom error cannot be less than the value")
                    end
                    errbot=value-NameValueArgs.errbot;
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

        function fig=kpolyplot(obj,type,parameter,xlab,lab_limit,scale,save)
            arguments
                obj Experiment
                type string {mustBeMember(type,{'formin','PRM'})}
                parameter string
                xlab string
                lab_limit double=1
                scale string="none"
                save logical=false
            end
            if type=="formin"
                fig=kpolyplot(obj.ForminList, parameter,xlab,lab_limit,scale,obj.opts,save);
            elseif type=="PRM"
                fig=PRMplot(obj.ForminList,parameter,xlab,scale,obj.opts,save);
            end
        end

        function fig=NTDplot(obj,type,parameter,xlab,scale,save)
            arguments
                obj Experiment
                type string {mustBeMember(type,{'formin','PRM'})}
                parameter string
                xlab string
                scale string="log2"
                save logical=false
            end
            if type=="PRM"
                fig = PRM_NTD_plot(obj.ForminList,parameter,xlab,scale,obj.opts,save);
            elseif type=="formin"
                fig=NTDplot(obj.ForminList,parameter,xlab,scale,obj.opts,save);
            end
        end

        function fig=forminbar(obj,save,kpolyscale,ratioscale)
            arguments
                obj Experiment
                save logical=false
                kpolyscale string {mustBeMember(kpolyscale,{'none','log2','log10','ln'})}="none"
                ratioscale string {mustBeMember(ratioscale,{'none','log2','log10','ln'})}="log2"
            end
            fig=forminbar(obj.ForminList,obj.opts,save,kpolyscale=kpolyscale,ratioscale=ratioscale);
        end

        function fig=polymerstat_change_plot(obj,parameter,stat,xlab,ylab,scale,minus,save)
            arguments
                obj Experiment
                parameter string
                stat string
                xlab string
                ylab string
                scale string="log2"
                minus logical=false
                save logical=false
            end
            fig=polymerstat_change_plot(obj.ForminList,parameter,stat,xlab,ylab,scale,obj.opts,minus,save);
        end

        function fig=formingraphic(obj,save)
            arguments
                obj Experiment
                save logical=false
            end
            nformins=length(obj.ForminList);
            fig=figure().empty(0,nformins);
            for i=1:nformins
                fig(i)=obj.ForminList(i).formingraphic(save);
            end
        end

        function fig=expdatabar(obj,scale,save,NameValueArgs)
            arguments
                obj Experiment
                scale string {mustBeMember(scale,{'none','log2','log10','ln'})}="none"
                save logical=false
                NameValueArgs.group string
            end
            if isfield(NameValueArgs,"group")
                fig=expdatabar(obj.data,obj.opts,scale,save,group=NameValueArgs.group);
            else
                fig=expdatabar(obj.data,obj.opts,scale,save);
            end
        end

        function makeresults(obj)
            arguments
                obj Experiment
            end

            set(groot,'defaultfigureposition',[400 250 900 750]) % helps prevent cut offs in figs
            
            %set colors and points
            obj.opts.update_points(max([obj.ForminList.PRMCount,length(obj.ForminList)]));
    
            % formin schematics
            obj.formingraphic(true);
    
            % Overview figures
            obj.forminbar(true);
            
            % Correlation plots
            
            % Change in Polymerization Rates vs Number of PRMs
            obj.NTDplot("formin","PRMCount",'Number of PRMs',"log2",true);
            
            % Polymerization Rates vs Number of PRMs
            obj.kpolyplot("formin","PRMCount",'Number of PRMs',20,"none",true);
            
            % Change in Polymerization Rates vs Length of FH1 Domain
            obj.NTDplot("formin","length",'Length of FH1 domain (1st PRM to FH2)',"log2",true);
            
            % Polymerization Rates vs Length of FH1 Domain
            obj.kpolyplot("formin","length",'Length of FH1 domain (1st PRM to FH2)',400,"none",true);
            
            % Change in Polymerization Rates vs Mean PRM size
            obj.NTDplot("formin","meanPRMsize",'Mean PRM size',"log2",true);
    
            % Change in Polymerization Rates vs Mean PRM size x Number of PRMs
            obj.NTDplot("formin","numPs",'Mean PRM size x #PRM',"log2",true);
            
            % Vs. PRM length
            obj.kpolyplot("PRM","size",'Length of PRMs',20,"none",true);
            
            % Distance from PRM to FH2
            obj.kpolyplot("PRM","dist_FH2",'Distance from PRM to FH2',20,"none",true);
            
            % Distance from PRM to end
            obj.kpolyplot("PRM","dist_NT",'Distance from PRM to N-term',20,"none",true);
            
            % Change in Polymerization Rates vs. PP length per individual PRM
            obj.NTDplot("PRM","size",'Length of PRM',"log2",true);
            
            % Change in Polymerization Rates vs. PP dist to FH2 per individual PRM
            obj.NTDplot("PRM","dist_FH2",'Distance from PRM to FH2',"log2",true);
            
            % Change in Polymerization Rates vs. PP dist to NT per individual PRM
            obj.NTDplot("PRM","dist_NT",'Distance from PRM to N-term',"log2",true);
            
            % Polymer stats plots
            
            % Change in polymer stats vs. FH1 length per individual PRM
            obj.polymerstat_change_plot("fh1length","POcclude",'FH1 length','Occlusion Probability',"log2",true,true);
            obj.polymerstat_change_plot("fh1length","POcclude_Base",'FH1 length','Occlusion Probability at the Barbed End',"log2",true,true);
            obj.polymerstat_change_plot("fh1length","Prvec0",'FH1 length','Concentration at the Barbed End',"log2",false,true);
            
            % Change in polymer stats vs. PP dist to NT per individual PRM
            obj.polymerstat_change_plot("dist_NT","POcclude",'Distance from PRM to N-term','Occlusion Probability',"log2",true,true);
            obj.polymerstat_change_plot("dist_NT","POcclude_Base",'Distance from PRM to N-term','Occlusion Probability at the Barbed End',"log2",true,true);
            obj.polymerstat_change_plot("dist_NT","Prvec0",'Distance from PRM to N-term','Concentration at the Barbed End',"log2",false,true);
            
            % Change in polymer stats vs. Fractional distance from FH2 per individual PRM
            obj.polymerstat_change_plot("FH2dist_frac","POcclude",'Fractional Distance from PRM to FH2','Occlusion Probability',"log2",true,true);
            obj.polymerstat_change_plot("FH2dist_frac","POcclude_Base",'Fractional Distance from PRM to FH2','Occlusion Probability at the Barbed End',"log2",true,true);
            obj.polymerstat_change_plot("FH2dist_frac","Prvec0",'Fractional Distance from PRM to FH2','Concentration at the Barbed End',"log2",false,true);
            
            % fit to experimental data 
            obj.expdatabar("none",true);

            % Make overview tables
            obj.formintable(true);
            obj.PRMtable(true);
        end


        function set_gating(obj,forminname,gating)
            arguments
                obj Experiment
                forminname string 
                gating double
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
            formin.gating=gating;
        end
    end
end