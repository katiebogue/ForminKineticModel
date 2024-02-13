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


        function out=SOSlist(obj,scaled,NameValueArgs)
            arguments
                obj Experiment
                scaled logical=false %only applies if group is specified
                NameValueArgs.formin string
                NameValueArgs.group string
                NameValueArgs.recalc logical=true
            end
            if isfield(NameValueArgs,"formin")
                for i=1:length(obj.data)
                    if obj.data(i).formin.name==NameValueArgs.formin
                        if NameValueArgs.recalc
                            kpolys=obj.data(i).formin.kpoly;
                        end
                        out=getSOS(obj.data(i),1);
                        return
                    end
                end
            end
            out=-1;

            if NameValueArgs.recalc
                kpolys=obj.ForminList.kpoly;
            end
            if isfield(NameValueArgs,"group")
                groupdata=obj.data;
                for i=length(groupdata):-1:1
                    if ~ismember(NameValueArgs.group,groupdata(i).groups)
                        groupdata(i)=[];
                    end
                end
                if scaled
                    if groupdata(1).type=="ratio"
                        scaler=1;
                    else
                        expdata=[groupdata.value];
                        minloc=expdata==min(expdata);
                        scaler=groupdata(minloc).formin.lastkpoly.(groupdata(minloc).type)/min(expdata);
                    end
                else
                    scaler=1;
                end
                for i=1:length(groupdata)
                    tstruct.name=groupdata(i).formin.name;
                    tstruct.value=getSOS(groupdata(i),scaler);
                    if class(out)~="struct"
                        out=tstruct;
                    else
                        out(1,end+1)=tstruct;
                    end
                end
            else
                for i=1:length(obj.data)
                    tstruct.name=obj.data(i).formin.name;
                    tstruct.value=getSOS(obj.data(i),1);
                    if class(out)~="struct"
                        out=tstruct;
                    else
                        out(1,end+1)=tstruct;
                    end
                end
            end

            function val=getSOS(datastruct,scaler)
                if scaler==0
                    scaler=1;
                end
                simvalue=datastruct.formin.lastkpoly.(datastruct.type)./scaler;
                max=datastruct.value+datastruct.errtop;
                min=datastruct.value-datastruct.errbot;
                if (simvalue<=max) && (simvalue>=min)
                    val=0;
                elseif isnan(simvalue)
                    val=10000;
                else
                    val=abs((simvalue-datastruct.value)^2);
                end
            end
        end

        function out=SOS(obj,scaled,NameValueArgs)
            arguments
                obj Experiment
                scaled logical=false %only applies if group is specified
                NameValueArgs.formin string
                NameValueArgs.group string
                NameValueArgs.recalc logical=true
            end
            if isfield(NameValueArgs,"formin")
                out=obj.SOSlist(formin=NameValueArgs.formin,recalc=NameValueArgs.recalc);
                return
            end
            if isfield(NameValueArgs,"group")
                list=obj.SOSlist(scaled,group=NameValueArgs.group,recalc=NameValueArgs.recalc);
            else
                list=obj.SOSlist(scaled,recalc=NameValueArgs.recalc);
            end

            out=sum([list.value]);
        end

        function T=runfminsearch(obj,NameValueArgs)
            arguments
                obj Experiment
                NameValueArgs.options struct=optimset('MaxFunEvals',5000,'MaxIter',5000,'TolX',10^(-4),'TolFun',10^(-4))  % optimset for fminsearch
                NameValueArgs.groups cell ={"n/a"}
                NameValueArgs.rbot double=-2
                NameValueArgs.rtop double=4
                NameValueArgs.iterations double=1
                NameValueArgs.weight % group to have SOS multiplied by 10
                NameValueArgs.ogtable logical =false % if true, just return a table with the current settings
            end

            options=NameValueArgs.options;
            groups=NameValueArgs.groups;
            rbot=NameValueArgs.rbot;
            rtop=NameValueArgs.rtop;
            iterations=NameValueArgs.iterations;
        

            vals={obj.opts.k_cap,obj.opts.k_del,obj.opts.r_cap,obj.opts.r_del,obj.opts.k_rel};
            fval=getSOS(vals{:},false);
            ogtable=maketable();
            if NameValueArgs.ogtable
                T=ogtable;
                return
            end
            T = table('Size',[2*iterations length(ogtable.Properties.VariableNames)],'VariableNames',ogtable.Properties.VariableNames,'VariableTypes',varfun(@class,ogtable,'OutputFormat','cell'));
            index=0;

            for i=1:iterations
                disp(i)
                r=(rtop-rbot)*rand(1,5)+rbot;

                % not using scaled
                [params,fval] = fminsearch(@(z) getSOS(abs(10^(z(1))),abs(10^(z(2))),abs(10^(z(3))),abs(10^(z(4))),abs(10^(z(5))),false),[r(1);r(2);r(3);r(4);r(5)],options); 
                vals=num2cell(abs(10.^(params)));
                T(index+1,:)=maketable();
                index=index+1;

                % using scaled
                [params,fval] = fminsearch(@(z) getSOS(abs(10^(z(1))),abs(10^(z(2))),abs(10^(z(3))),abs(10^(z(4))),abs(10^(z(5))),true),[r(1);r(2);r(3);r(4);r(5)],options); 
                vals=num2cell(abs(10.^(params)));
                T(index+1,:)=maketable();
                index=index+1;
            end

            % return to OG opts
            obj.opts.k_cap=ogtable.k_cap;
            obj.opts.k_del=ogtable.k_del;
            obj.opts.r_cap=ogtable.r_cap;
            obj.opts.r_del=ogtable.r_del;
            obj.opts.k_rel=ogtable.k_rel;
            
            function val=getSOS(k_cap,k_del,r_cap,r_del,k_rel,scaled)
                obj.opts.k_cap=k_cap;
                obj.opts.k_del=k_del;
                obj.opts.r_cap=r_cap;
                obj.opts.r_del=r_del;
                obj.opts.k_rel=k_rel;
                
                if groups{1}=="n/a"
                    val=obj.SOS();
                else
                    val=0;
                    for j=1:length(groups)
                        weight=1;
                        if isfield(NameValueArgs,"weight")
                            if groups{j}==NameValueArgs.weight
                                weight=10;
                            end
                        end
                        if j==1
                            val=val+(weight*obj.SOS(scaled,group=groups{j}));
                        else
                            val=val+(weight*obj.SOS(scaled,group=groups{j},recalc=false));
                        end
                    end
                end

            end

            function tab=maketable()
                tab=table;
                tab.fval=fval;
                [tab.k_cap,tab.k_del,tab.r_cap,tab.r_del,tab.k_rel]=vals{:};
                grouplist=unique([obj.data.groups]);
                grouptitle=strcat(grouplist," scaled");
                tab.all=obj.SOS();
                allscaled=0;
                for j=1:length(grouplist)
                    tab.(grouplist(j))=obj.SOS(group=grouplist(j),recalc=false);
                    tab.(grouptitle(j))=obj.SOS(true,group=grouplist(j),recalc=false);
                    allscaled=allscaled+tab.(grouptitle(j));
                end

                tab.all_scaled=allscaled;
                
                formins=[obj.data.formin];
                forminnames=[formins.name];
                forminname_gating=strcat(forminnames,"-gating");
                forminname_c_PA=strcat(forminnames,"-c_PA");
                for j=1:length(obj.data)
                    tab.(forminnames(j))=obj.SOS(formin=forminnames(j),recalc=false);
                    tab.(forminname_gating(j))=formins(j).gating;
                    tab.(forminname_c_PA(j))=formins(j).c_PA;
                end
                
                optstable=rows2vars(obj.opts.optionstable);
                optstable=removevars(optstable,{'resultsdir','resultsfolder','python_path','k_cap','k_del','r_cap','r_del','k_rel','OriginalVariableNames'});
                tab=[tab,optstable];
            end
        end

        function applytable(obj,row)
            arguments
                obj Experiment
                row table
            end
            vars=row.Properties.VariableNames;
            optsvar=obj.opts;
            optsvar.cleareq;
            for i=1:length(vars)
                value=row.(vars{i});
                if ismember(vars{i},fieldnames(optsvar))
                    if class(optsvar.(vars{i}))=="double" && class(value)=="string"
                        value=str2double(value);
                    end
                    if optsvar.(vars{i})~=value
                        optsvar.(vars{i})=value;
                    end
                    continue
                end
                if ismember(vars{i},strcat([fieldnames(optsvar.equations)],"_eq"))
                    optsvar.set_equation(0,vars{i}(1:end-3),cellfun(@string,strsplit(char(value),{',','(',')'}),"UniformOutput",false))
                    continue
                end
                formins=[obj.ForminList.name];
                splitformin=strsplit(vars{i},{'-'});
                if length(splitformin)==2
                    if ismember(splitformin{1},formins)
                        j=find(formins==splitformin{1});
                        obj.ForminList(j).(splitformin{2})=value;
                    end
                end
            end
        end
    end
end