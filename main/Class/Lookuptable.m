classdef Lookuptable <handle 
    %LOOKUPTABLE Summary of this class goes here
    %   Detailed explanation goes here

    properties (GetAccess=public, SetAccess=protected)
        single 
        dimer 
        double 
    end

    properties (SetAccess=private)
        StatNames % names of polymer stats/ outputs from simulation
        NNames % list of Ns that have been added
        AddedStats struct % each stat/ output as a Submodel
        AddedNs struct % each N as a Submodel
        NList Submodel % array of N values as Submodel
        max_length Submodel % Max filament length before extrapolation
    end

    methods
        function value = get.StatNames(obj)
            value = string(fieldnames(obj.AddedStats));
        end
        
        function value = get.NNames(obj)
            value = string(fieldnames(obj.AddedNs));
        end

        function value= get.max_length(obj)
            value=Submodel(find_maxN("single"),find_maxN("double"),find_maxN("dimer"));
            function max=find_maxN(type)
                Ns=fieldnames(obj.(type));
                max=0;
                for i=1:length(Ns)
                    N=str2double(erase(Ns{i},"N"));
                    if N>max
                        max=N;
                    end
                end
            end
        end

        function value = get.NList(obj)
            value=Submodel(find_N("single"),find_N("double"),find_N("dimer"));
            function Nlist=find_N(type)
                Ns=fieldnames(obj.(type));
                Nlist=zeros(1,length(Ns));
                for i=1:length(Ns)
                    N=str2double(erase(Ns{i},"N"));
                    Nlist(i)=N;
                end
            end
        end

        function value = get.AddedStats(obj)
            value=struct;
            addstats("single")
            addstats("double")
            addstats("dimer")
            function addstats(type)
                Ns=fieldnames(obj.(type));
                for i=1:length(Ns)
                    vars=fieldnames(obj.(type).(Ns{i}));
                    for j=1:length(vars)
                        var=vars{j};
                        if ~isfield(value,var)
                            value.(var)=Submodel();
                        end
                        value.(var).(type).(Ns{i})=obj.(type).(Ns{i}).(var);
                    end
                end
            end
        end

        function value = get.AddedNs(obj)
            value=struct;
            addNs("single")
            addNs("double")
            addNs("dimer")
            function addNs(type)
                Ns=fieldnames(obj.(type));
                for i=1:length(Ns)
                    N=Ns{i};
                    if ~isfield(value,N)
                        value.(N)=Submodel();
                    end
                    value.(N).(type)=obj.(type).(N);
                end
            end
        end
       
    end

    methods (Access=public)
        function obj = Lookuptable(ltvar)
            arguments
                ltvar struct
            end
            fields=fieldnames(ltvar);
            if nargin == 1
                if length(fields)>3
                    error("Too many fields in input lookuptable variable")
                elseif isempty(fields)
                    error("No fields in input lookuptable variable")
                else
                   for i=1:length(fields)
                       field=fields{i};
                       if ismember(field,{'single','double','dimer'})
                           obj.(field)=ltvar.(field);
                       else
                           error("Input lookuptable variable contains invalid field name: %s",field)
                       end
                   end
                end

            end
        end

        function addN(obj,Nstruct)
            Ntype=Nstruct.type;
            N=num2str(Nstruct.N(1));
            obj.(Ntype).(strcat("N",N))=Nstruct;
        end

        function output = getstat(obj,NameValueArgs)
            arguments
                obj Lookuptable
                NameValueArgs.N double =[]
                NameValueArgs.type string {mustBeMember(NameValueArgs.type,{"single","double","dimer"})}
                NameValueArgs.Stat string
                NameValueArgs.iSite double=[]
                NameValueArgs.Fil string {mustBeMember(NameValueArgs.Fil,{"a","b"})}
            end
            N=NameValueArgs.N;
            if isfield(NameValueArgs,"type")
                type=NameValueArgs.type;
            else
                type=[];
            end
            if isfield(NameValueArgs,"Stat")
                Stat=NameValueArgs.Stat;
            else
                Stat=[];
            end
            iSite=NameValueArgs.iSite;
            if isfield(NameValueArgs,"Fil")
                Fil=NameValueArgs.Fil;
            else
                Fil=[];
            end


            if ~isempty(N)
                if ~isempty(iSite) && iSite>N
                    error("iSite is larger than FH1 length")
                end
                NName=strcat("N",num2str(N));
                if isempty(Stat)
                    if ismember(NName,obj.NNames)
                        tempout=obj.AddedNs.(NName);
                        if ~isempty(type)
                            tempout=tempout.(type);
                            if ~isempty(Fil) && type~="single"
                                % get filaments from the non submodel
                            end
                        else
                            if ~isempty(Fil)
                                % get filaments from the non submodel
                            end
                        end
                    else
                        for i=1:length(obj.StatNames)
                            out1=obj.extrapolate(obj.StatNames(i));
                            if ~isempty(type)
                                out1=out1.(type);
                                if class(out1)=="scatteredInterpolant"
                                    tempout.(obj.StatNames(i))=zeros(1,N);
                                    for j=1:N
                                        tempout.(obj.StatNames(i))(1,j)=out1(N,j);
                                    end
                                elseif class(out1)=="Filament"
                                    if isempty(Fil)
                                        if class(out1.a)=="scatteredInterpolant"
                                            tempout.(obj.StatNames(i))=zeros(2,N);
                                            for j=1:N
                                                tempout.(obj.StatNames(i))(1,j)=out1.a(N,j);
                                                tempout.(obj.StatNames(i))(2,j)=out1.b(N,j);
                                            end
                                        else
                                            tempout.(obj.StatNames(i))=zeros(2,1);
                                            tempout.(obj.StatNames(i))(1,1)=out1.a(N);
                                            tempout.(obj.StatNames(i))(2,1)=out1.b(N);
                                        end
                                    else
                                        if class(out1.a)=="scatteredInterpolant"
                                            tempout.(obj.StatNames(i))=zeros(1,N);
                                            for j=1:N
                                                tempout.(obj.StatNames(i))(1,j)=out1.(Fil)(N,j);
                                            end
                                        else
                                            tempout.(obj.StatNames(i))=out1.(Fil)(N);
                                        end
                                    end
                                else
                                     tempout.(obj.StatNames(i))=out1(N);
                                end
                            else
                                types=["single","double","dimer"];
                                tempout=Submodel;
                                for k=1:3
                                    temptype=types(k);
                                    out1=out1.(temptype);
                                    if class(out1)=="scatteredInterpolant"
                                        tempout.(temptype).(obj.StatNames(i))=zeros(1,N);
                                        for j=1:N
                                            tempout.(temptype).(obj.StatNames(i))(1,j)=out1(N,j);
                                        end
                                    elseif class(out1)=="Filament"
                                        if isempty(Fil)
                                            if class(out1.a)=="scatteredInterpolant"
                                                tempout.(temptype).(obj.StatNames(i))=zeros(2,N);
                                                for j=1:N
                                                    tempout.(temptype).(obj.StatNames(i))(1,j)=out1.a(N,j);
                                                    tempout.(temptype).(obj.StatNames(i))(2,j)=out1.b(N,j);
                                                end
                                            else
                                                tempout.(temptype).(obj.StatNames(i))=zeros(2,1);
                                                tempout.(temptype).(obj.StatNames(i))(1,1)=out1.a(N);
                                                tempout.(temptype).(obj.StatNames(i))(2,1)=out1.b(N);
                                            end
                                        else
                                            if class(out1.a)=="scatteredInterpolant"
                                                tempout.(temptype).(obj.StatNames(i))=zeros(1,N);
                                                for j=1:N
                                                    tempout.(temptype).(obj.StatNames(i))(1,j)=out1.(Fil)(N,j);
                                                end
                                            else
                                                tempout.(temptype).(obj.StatNames(i))=out1.(Fil)(N);
                                            end
                                        end
                                    else
                                         tempout.(temptype).(obj.StatNames(i))=out1(N);
                                    end
                                end
                            end
                        end
                    end
                else
                    if ismember(NName,obj.NNames)
                        tempout=obj.AddedStats.(Stat);
                        if isempty(type)
                            tempout=Submodel(tempout.single.(NName),tempout.double.(NName),tempout.dimer.(NName));
                            if ~isempty(Fil)
                                % get filament here
                            end
                        else
                            tempout=tempout.(type).(NName);
                            if ~isempty(Fil)
                                % get filament here
                            end
                        end
                    else
                        out1=obj.extrapolate(Stat);
                        if ~isempty(type)
                            out1=out1.(type);
                            if class(out1)=="scatteredInterpolant"
                                tempout=zeros(1,N);
                                    for j=1:N
                                        tempout(1,j)=out1(N,j);
                                    end
                            elseif class(out1)=="Filament"
                                if isempty(Fil)
                                    if class(out1.a)=="scatteredInterpolant"
                                        tempout=zeros(2,N);
                                        for j=1:N
                                            tempout(1,j)=out1.a(N,j);
                                            tempout(2,j)=out1.b(N,j);
                                        end
                                    else
                                        tempout=zeros(2,1);
                                        tempout(1,1)=out1.a(N);
                                        tempout(2,1)=out1.b(N);
                                    end
                                else
                                    if class(out1.a)=="scatteredInterpolant"
                                        tempout=zeros(1,N);
                                        for j=1:N
                                            tempout(1,j)=out1.(Fil)(N,j);
                                        end
                                    else
                                        tempout=out1.(Fil)(N);
                                    end
                                end
                            else
                                 tempout=out1(N);
                            end
                        else
                            types=["single","double","dimer"];
                            tempout=Submodel;
                            for k=1:3
                                temptype=types(k);
                                out1=out1.(temptype);
                                if class(out1)=="scatteredInterpolant"
                                    tempout.(temptype)=zeros(1,N);
                                    for j=1:N
                                        tempout.(temptype)(1,j)=out1(N,j);
                                    end
                                elseif class(out1)=="Filament"
                                    if isempty(Fil)
                                        if class(out1.a)=="scatteredInterpolant"
                                            tempout.(temptype)=zeros(2,N);
                                            for j=1:N
                                                tempout.(temptype)(1,j)=out1.a(N,j);
                                                tempout.(temptype)(2,j)=out1.b(N,j);
                                            end
                                        else
                                            tempout.(temptype)=zeros(2,1);
                                            tempout.(temptype)(1,1)=out1.a(N,1);
                                            tempout.(temptype)(2,1)=out1.b(N,1);
                                        end
                                    else
                                        if class(out1.a)=="scatteredInterpolant"
                                            tempout.(temptype)=zeros(1,N);
                                            for j=1:N
                                                tempout.(temptype)(1,j)=out1.(Fil)(N,j);
                                            end
                                        else
                                            tempout.(temptype)=out1.(Fil)(N,1);
                                        end
                                    end
                                else
                                     tempout.(temptype)=out1(N);
                                end
                            end
                        end
                    end
                end
            else
                if isempty(Stat)
                    if isempty(type)
                        tempout=obj;
                        if ~isempty(Fil)
                            % get filament somehow
                        end
                    else
                        if isempty(Fil) || type=="single"
                            tempout=obj.(type);
                        else
                            % get filament somehow
                        end
                    end
                else
                    tempout=obj.AddedStats.(Stat);
                    if ~isempty(type)
                        if isempty(Fil) || type=="single"
                            tempout=tempout.(type);
                        else
                            % get filament somehow
                        end
                    elseif ~isempty(Fil)
                         % get filament somehow
                    end
                end
            end
            
            output=tempout;

        end

        function output=extrapolate(obj,stat)
            statobj=obj.AddedStats.(stat);
            if stat=="type"
                singlefxn= @(n) "single";
                doubleafxn= @(n) "double";
                dimerafxn= @(n) "double";
                doublebfxn= @(n) "dimer";
                dimerbfxn= @(n) "dimer";
            elseif class(statobj.single.(obj.NNames(1)))=="string"
                singlefxn= @(n) "unable to extrapolate";
                doubleafxn= @(n) "unable to extrapolate";
                dimerafxn= @(n) "unable to extrapolate";
                doublebfxn= @(n) "unable to extrapolate";
                dimerbfxn= @(n) "unable to extrapolate";
            else
                singlefxn=getfxn("single",1);
                doubleafxn=getfxn("double",1);
                dimerafxn=getfxn("dimer",1);
                doublebfxn=getfxn("double",2);
                dimerbfxn=getfxn("dimer",2);
            end
            
            output=Submodel(singlefxn,doubleafxn,dimerafxn,doublebfxn,dimerbfxn);

            function fxn=getfxn(type,fil)
                stats=struct2cell(statobj.(type));
                if size(stats{1},1)==1
                    if size(stats{1},2)==1 && size(stats{2},2)==1
                        stats=cell2mat(stats);
                        [stats,sortIn]=sort(stats);
                        fit = polyfit(obj.NList.(type)(sortIn),stats(:)',1);
                        fxn=@(n) fit(1)*n+fit(2);
                    elseif size(stats{1},2)==2 && size(stats{2},2)==2
                        stats=cell2mat(stats);
                        [stats,sortIn]=sort(stats);
                        fit = polyfit(obj.NList.(type)(sortIn(:,fil)),stats(:,fil),1);
                        fxn=@(n) fit(1)*n+fit(2);
                    else
                        numrows=obj.max_length.(type)*size(stats,1);
                        mat=zeros(numrows,3);
                        for i=1:size(stats,1)
                            for j=1:size(stats{i})
                                mat(i+j,1)=obj.NList.(type)(i);
                                mat(i+j,2)=j;
                                mat(i+j,3)=stats{i}(j);
                            end
                        end
                        fxn = scatteredInterpolant(mat(:,1),mat(:,2),mat(:,3),'linear','nearest');
                    end
                else
                    numrows=obj.max_length.(type)*size(stats,1);
                        mat=zeros(numrows,3);
                        for i=1:size(stats,1)
                            for j=1:size(stats{i},2)
                                mat(i+j,1)=obj.NList.(type)(i);
                                mat(i+j,2)=j;
                                mat(i+j,3)=stats{i}(fil,j);
                            end
                        end
                        fxn = scatteredInterpolant(mat(:,1),mat(:,2),mat(:,3),'linear','nearest');
                end
            end
        end
    end
end