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
        pdf_name % name of pdf to save any figures to
    end

    methods
        function obj = Options(python_path,kpoly_type)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.python_path=python_path;
            obj.kpoly_type=kpoly_type;
            obj.update_points;
        end

        function set_equation(obj,preset,step,vars)
            arguments
                obj
                preset double {mustBeMember(preset,[0,1])}
            end
            arguments (Repeating)
                step (1,1) string {mustBeMember(step,{"kcap","kdel","rcap","rdel","krel"})}
                vars cell {validvars(vars)}
            end
            
            if isempty(obj.equations)
                obj.equations=struct;
            end
            if preset==0
                invars=vars;
                obj.equations.(step)=makeeq();
            elseif preset==1
                invars={"POcclude","1-","c_PA","linear"};
                obj.equations.kcap=makeeq();
                invars={"POcclude","1-base","Prvec0","amino"};
                obj.equations.kdel=makeeq();
                invars={"size","negexp","gating","linear"};
                obj.equations.rcap=makeeq();
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

    end

    methods (Access=private)
        function validvars(a)
                % Must have valid variable names and scale types and be
                % even
                if logical(mod(length(a),2))
                    eid = 'validvars:sizeNotEven';
                    msg = 'An even number of inputs must be provided to vars.';
                    throwAsCaller(MException(eid,msg))
                end
                valid_vars=obj.lookup.StatNames;
                scale_types={"linear","negexp","exp","1-","amino","1-base"};
                for i=1:2:length(a)
                    if ~ismember(a{i},valid_vars) && ~ismember(a{i},{"dist_FH2","dist_NT", "size","gating","c_PA","dist_FH2_start"})
                        eid = 'validvars:varNotValid';
                        msg = 'Variables to include in equations must either be in the reference lookuptable or "dist_FH2","dist_NT","size".';
                        throwAsCaller(MException(eid,msg))
                    end
                    if ~ismember(a{i+1},scale_types)
                        valid_scale_string="";
                        for j=1:length(scale_types)
                            if j==length(scale_types)
                                valid_scale_string=strcat(valid_scale_string,scale_types{i});
                            else
                                valid_scale_string=strcat(valid_scale_string,scale_types{i},", ");
                            end
                        end

                        eid = 'validvars:scaleTypeNotValid';
                        msg = strcat('Sacle type must be one of the following: ',valid_scale_string,".");
                        throwAsCaller(MException(eid,msg))
                    end
                end
                
            end
    end

end