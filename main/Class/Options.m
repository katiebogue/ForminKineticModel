classdef Options <handle
    %OPTIONS Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetObservable, AbortSet)
        python_path string
        kpoly_type string {mustBeMember(kpoly_type,{'capture','3st','4st'})}='4st'
        equations struct
        lookup Lookuptable

        % reading sequence options (see get_formin_info.py):
        min_lenPRM double % minimum PRM length (including interruptions)
        nInt double % number of allowed interruptions (counts once for each amino acid-- i.e. if the max int len is 2, AA is acceptable but counts as 2 interruptions)
        max_lenInt double % maximum interruption length
        min_nP double % minimum number of Ps 
        NTopt double % FH1 NT definition options 
        CTopt double % FH1 CT definition options
        PRMscanopt double % PRM scan options

        colors (:,1) % list of hexcodes to be used for plotting
        shapes (:,1) % list of shape symbols to be used for plotting
        resultsdir % location of folder to save results folder to
        resultsfolder % results folder name

        equationstext struct % string for each equation 

        notes string % notes to add to output

        % parameters
        k_cap double % rate constant for PRM + profilin-actin binding (capture) | μM^(-1)s^(-1)
        k_del double % rate constant for barbed end + PRM-profilin-actin binding (delivery) | μM^(-1)s^(-1)
        r_cap double % rate constant for PRM + profilin-actin dissociation (reverse capture) | s^(-1)
        r_del double % rate constant for barbed end + PRM-profilin-actin dissociation (reverse delivery) | s^(-1)
        k_rel double % rate constant for PRM + profilin-actin-barbed end dissociation (release) | s^(-1)
    end

    methods
        function obj = Options(lookup,python_path,kpoly_type,resultsdir,k_cap,k_del,r_cap,r_del,k_rel)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            if nargin>0
            obj.lookup=lookup;
            obj.python_path=python_path;
            obj.kpoly_type=kpoly_type;
            obj.update_points;
            obj.resultsdir=resultsdir;
            obj.update_results_folder;
            obj.k_cap = k_cap;
            obj.k_del = k_del;
            obj.r_cap = r_cap;
            obj.r_del = r_del;
            obj.k_rel = k_rel;
            obj.set_FH1();
            end
        
        end

        function update_results_folder(obj)
            time= datetime('now', 'Format','yyyy-MM-dd HH-mm');
            time= string(time);
            obj.resultsfolder= 'RESULTS_' + time;
        end

        function set_FH1(obj,NameValueArgs)
            arguments
                obj Options
                NameValueArgs.min_lenPRM double=4 % minimum PRM length (including interruptions)
                NameValueArgs.nInt double=1 % number of allowed interruptions (counts once for each amino acid-- i.e. if the max int len is 2, AA is acceptable but counts as 2 interruptions)
                NameValueArgs.max_lenInt double=1 % maximum interruption length
                NameValueArgs.min_nP double=1 % minimum number of Ps 
                NameValueArgs.NTopt double {mustBeInRange(NameValueArgs.NTopt,1,5)}=1 % FH1 NT definition options 
                NameValueArgs.CTopt double {mustBeInRange(NameValueArgs.CTopt,1,3)}=1 % FH1 CT definition options
                NameValueArgs.PRMscanopt double {mustBeInRange(NameValueArgs.PRMscanopt,1,2)}=1 % PRM scan options
            end
            obj.min_lenPRM=NameValueArgs.min_lenPRM;
            obj.nInt=NameValueArgs.nInt; 
            obj.max_lenInt=NameValueArgs.max_lenInt; 
            obj.min_nP=NameValueArgs.min_nP; 
            obj.NTopt=NameValueArgs.NTopt;
            obj.CTopt=NameValueArgs.CTopt;
            obj.PRMscanopt=NameValueArgs.PRMscanopt;
        end

        function set_equation(obj,preset,step,vars)
            arguments
                obj
                preset double {mustBeMember(preset,[0,1])}
            end
            arguments (Repeating)
                step string {mustBeMember(step,{'kcap','kdel','rcap','rdel','krel'})}
                vars (1,:) cell 
            end
            
                
            
            if isempty(obj.equations)
                obj.equations=struct;
            end
            if isempty(obj.equationstext)
                obj.equationstext=struct;
            end
            if preset==0
                vars=vars{:};
                if logical(mod(length(vars),2))
                    if vars{end}==""
                        vars=vars(1:end-1);
                    else
                        eid = 'validvars:sizeNotEven';
                        msg = 'An even number of inputs must be provided to vars.';
                        throwAsCaller(MException(eid,msg))
                    end
                end
                valid_vars=obj.lookup.StatNames;
                scale_types={'linear','negexp','exp','1-','amino','1-base'};
                for k=1:2:length(vars)
                    if ~ismember(vars{k},valid_vars) && ~ismember(vars{k},{'dist_FH2','dist_NT', 'size','gating','c_PA','dist_FH2_start'})
                        eid = 'validvars:varNotValid';
                        msg = 'Variables to include in equations must either be in the reference lookuptable or "dist_FH2","dist_NT","size".';
                        throwAsCaller(MException(eid,msg))
                    end
                    if ~ismember(vars{k+1},scale_types)
                        valid_scale_string="";
                        for j=1:length(scale_types)
                            if j==length(scale_types)
                                valid_scale_string=strcat(valid_scale_string,scale_types{k});
                            else
                                valid_scale_string=strcat(valid_scale_string,scale_types{k},", ");
                            end
                        end
    
                        eid = 'validvars:scaleTypeNotValid';
                        msg = strcat('Sacle type must be one of the following: ',valid_scale_string,".");
                        throwAsCaller(MException(eid,msg))
                    end
                end
                invars=vars;
                if invars{end}==""
                    invars=invars(1:end-1);
                end
                obj.equations.(step{:})=makeeq();
                obj.equationstext.(strcat(step{:},"_eq"))=[sprintf("%s(%s),",invars{1:end})];
            elseif preset==1
                invars={"POcclude","1-","c_PA","linear"};
                obj.equations.kcap=makeeq();
                obj.equationstext.kcap_eq=[sprintf("%s(%s),",invars{1:end})];
                invars={"POcclude","1-base","Prvec0","amino","gating","linear"};
                obj.equations.kdel=makeeq();
                obj.equationstext.kdel_eq=[sprintf("%s(%s),",invars{1:end})];
                invars={"size","negexp"};
                obj.equations.rcap=makeeq();
                obj.equationstext.rcap_eq=[sprintf("%s(%s),",invars{1:end})];
            end
            
            steps={"kcap","kdel","rcap","rdel","krel"};
            for j=1:length(steps)
                if ~isfield(obj.equations,steps{j})
                    obj.equations.(steps{j})=@(PRM) 1;
                end
            end
            

            function fxn = makeeq()
                fxn=@(PRM) 1;
                for i=1:2:length(invars)
                    if invars{i+1}=="linear"
                        fxn=@(PRM) fxn(PRM)*PRM.(invars{i});
                    elseif invars{i+1}=="negexp"
                        fxn=@(PRM) fxn(PRM)*exp(-1*PRM.(invars{i}));
                    elseif invars{i+1}=="exp"
                        fxn=@(PRM) fxn(PRM)*exp(PRM.(invars{i}));
                    elseif invars{i+1}=="1-"
                        fxn=@(PRM) fxn(PRM)*(1-PRM.(invars{i}));
                    elseif invars{i+1}=="amino"
                        fxn=@(PRM) fxn(PRM)*(1.0e33*(PRM.(invars{i}))/(27*6.022e23));
                    elseif invars{i+1}=="1-base"
                        fxn=@(PRM) fxn(PRM)*(1-PRM.(strcat(invars{i},"_Base")));
                    end
                end
            end
            
        end

        function cleareq(obj)
            eqfields=fieldnames(obj.equationstext);
            for i=1:length(eqfields)
                obj.equationstext=rmfield(obj.equationstext,eqfields{i});
            end
            eqfields=fieldnames(obj.equations);
            for i=1:length(eqfields)
                obj.equations.(eqfields{i})=@(PRM) 1;
            end
        end

        function update_points(obj,input,overwrite)
            arguments
                obj Options
                input {mustBeA(input,["double","Formin","Experiment"])}=37
                overwrite logical=false % if false, will only change if the new number of colors is less than the old
            end
            if overwrite
                oldmin=0;
            else
                oldmin=length(obj.colors);
            end

            if class(input)=="double"
                min=input;
                if oldmin>min
                    min=oldmin;
                end
            elseif class(input)=="Formin"
                min=input.PRMCount;
                if oldmin>min
                    min=oldmin;
                end
            elseif class(input)=="Experiment"
                min=length(input.ForminList);
                if oldmin>min
                    min=oldmin;
                end
            end
            
            [obj.colors,obj.shapes]=makepoints(min,'w');
        end

        function tab=optionstable(obj)
            
            tabstruct=struct;

            addel(obj);
            addel(obj.equationstext);

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

        function applytable(obj,row)
            arguments
                obj Options
                row table
            end
            vars=row.Properties.VariableNames;
            obj.cleareq;
            for i=1:length(vars)
                value=row.(vars{i});
                if ismember(vars{i},fieldnames(obj))
                    if class(obj.(vars{i}))=="double" && class(value)=="string"
                        value=str2double(value);
                    end
                    if obj.(vars{i})~=value
                        obj.(vars{i})=value;
                    end
                    continue
                end
                if ismember(vars{i},strcat([fieldnames(obj.equations)],"_eq"))
                    obj.set_equation(0,vars{i}(1:end-3),cellfun(@string,strsplit(char(value),{',','(',')'}),"UniformOutput",false))
                    continue
                end
            end
        end

    end

end