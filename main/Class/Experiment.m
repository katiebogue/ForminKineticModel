classdef Experiment 
    %EXPERIMENT Contains a group of formins and accompanying experimental
    % data, plus functions for plotting, calculating, fitting, etc
    %   
    % Construction:
    %   obj = EXPERIMENT(options,file,type,cPA)
    %       Inputs:
    %           options : Options object
    %           file    : (String).txt file with comma separated values 
    %                       of the format 'name','sequence' or 'name','uniprotID'
    %           type    : (String) which type ('seq' or 'uniprot') the
    %                       information in file is in
    %           cPA     : (double) concentration of profilin-actin to
    %                       set each of the input formins to (default is 0.88)
    % 
    % Following general construction, likely need to set gating factors and
    % add experimental data:
    %   EXPERIMENT.set_gating(forminname, gating)
    %   obj = EXPERIMENT.add_data(forminname,value,type,NameValueArgs)
    % 
    % See also FORMIN, OPTIONS, PRM.

    properties
        ForminList (:,:) Formin % Formins to group together in this experiment
        opts Options % Options object
        data (:,:) struct % experimental data for formin constructs, struct with properties type, formin, value, groups, errtop, and errbot
    end

    methods
        function obj = Experiment(options,file,type,cPA)
            %EXPERIMENT Construct an instance of Experiment class
            %   
            % obj = EXPERIMENT(options,file,type,cPA)
            %
            %   Sets formins according to input file.
            % 
            %   Inputs:
            %       options : Options object
            %       file    : (String).txt file with comma separated values 
            %                   of the format 'name','sequence' or 'name','uniprotID'
            %       type    : (String) which type ('seq' or 'uniprot') the
            %                   information in file is in
            %       cPA     : (double) concentration of profilin-actin to
            %                   set each of the input formins to (default is 0.88)
            % 
            % See also EXPERIMENT, FORMIN.
            arguments
                options Options
                file string % .txt file, comma separated
                type string {mustBeMember(type,{'seq','uniprot'})}
                cPA double =0.88 % [profilin-actin] to assign to all formins
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
            %ADD_FORMIN add a formin to obj.ForminList
            %
            %   obj = EXPERIMENT.ADD_FORMIN(formin)
            %
            %   Does not change the existing object, function output must
            %   be set to itself (or a new object)
            %
            % See also FORMIN, EXPERIMENT.
            arguments
                obj Experiment
                formin Formin
            end
            obj.ForminList(1,end+1)=formin;
        end

        function obj= set.opts(obj,input)
            % run syncforminsettings whenever setting the opts property
            obj.opts=input;
            obj.syncforminsettings;
        end

        function syncforminsettings(obj)
            %SYNCFORMINSETTINGS updates the opts property of each formin in
            %   ForminList to be the opts property of the Experiment object
            %
            %   EXPERIMENT.SYNCFORMINSETTINGS
            %
            %   Results in running update_FH1 for each formin since there
            %   are listeners for changes in opts
            %
            % See also FORMIN, EXPERIMENT, FORMIN/UPDATE_FH1.
            if isfield(obj,"ForminList")
                for i=1:length(obj.ForminList)
                    if obj.ForminList(i).opts~=obj.opts
                        obj.ForminList(i).opts=obj.opts;
                    end
                end
            end
        end

        function T=formintable(obj,save)
            %FORMINTABLE create formin property table
            %
            %   tab = EXPERIMENT.FORMINTABLE create table with all of the
            %   properties of the formins in ForminList
            %
            %   tab = EXPERIMENT.FORMINTABLE(1) create and save table with all of the
            %   properties of the formins in ForminList 
            %
            %   Saves "Per Formin Data.csv" to obj.opts.resultsdir,obj.opts.resultsfolder
            %
            %   Does not save the lastkpoly property.
            %   
            %   If the property is a FilType, saves a separate variable for
            %   single, double, dimer, and dimer/double.
            %   
            %   Rows are formin names, columns are properties.
            % 
            % 
            % See also FORMIN, EXPERIMENT, OPTIONS.
            arguments
                obj Experiment
                save logical=false
            end
            fh1names=[obj.ForminList.name];
            T= table('RowNames',fh1names);
            fields=fieldnames([obj.ForminList]);
            for i=1:length(fields)
                varname=fields{i};
                if varname=="lastkpoly"
                    continue
                end
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
            %PRMTABLE create PRM property table
            %
            %   tab = EXPERIMENT.PRMTABLE create table with all of the
            %   properties of the PRMs in PRMList of each ForminList
            %
            %   tab = EXPERIMENT.PRMTABLE(1) create and save table with all of the
            %   properties of the PRMs in PRMList of each ForminList
            %
            %   Saves "Per PRM Data.csv" to obj.opts.resultsdir,obj.opts.resultsfolder
            %
            %   Does not save the lastkpoly property.
            %   
            %   If the property is a FilType, saves a separate variable for
            %   single, double, dimer, and dimer/double.
            %   
            %   If the property is a Filament, saves a separate variable for
            %   a and b.
            %   
            %   Row for every PRM, columns are properties.
            % 
            % See also FORMIN, EXPERIMENT, OPTIONS, PRM, PRM.GETFIELDS, PRM.GETPROP.
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
            %ADD_DATA add data for a specific formin to obj.data
            %
            %   obj = EXPERIMENT.ADD_DATA(forminname,value,type,NameValueArgs)
            %
            %   Inputs:
            %       forminname  : (String) name of the formin, must match
            %                       one and only one name in ForminList
            %       value       : (double) experimental data value
            %       type        : (String) kpoly data type, either 'single','double','dimer', or 'ratio'
            %       errplus     : (double) how much to add to value to get the top of the error bars (NameValueArgs)
            %       errminus    : (double) how much to subtract from value to get the bottom of the error bars (NameValueArgs)
            %       errtop      : (double) the value of the top of the error bars  (NameValueArgs)
            %       errbot      : (double) the value of the bottom of the error bars (NameValueArgs)
            %       errperc     : (double) percent error (must be less than or equal to 1) (NameValueArgs)
            %       groups      : (String) experiment groups, will be plotted together (NameValueArgs, default is empty string)
            %
            %   Can only supply one type of error values, either errperc or
            %   a combination of (errtop or errplus) and (errbot and
            %   errminus)
            %
            %   obj.data is updated to have an additional entry. Each entry
            %   is a struct object with properties type, formin, value,
            %   groups, errtop, and errbot
            % 
            %   Does not change the existing object, function output must
            %   be set to itself (or a new object)
            %
            % See also FORMIN, EXPERIMENT.
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
            %KPOLYPLOT create a scatterplot of polymerization vs. a
            %parameter for either all formins or all PRMs
            %
            %   fig = EXPERIMENT.KPOLYPLOT(type,parameter,xlab) creates
            %   scatterplots of kpoly vs parameter with specified xlabel
            %   xlab
            %
            %   fig = EXPERIMENT.KPOLYPLOT(type,parameter,xlab,lab_limit) creates
            %   scatterplots of kpoly vs parameter with specified xlabel
            %   xlab and, if type = "formin," with points above value 
            %   lab_limit labeled by formin
            %
            %   fig = EXPERIMENT.KPOLYPLOT(type,parameter,xlab,lab_limit,scale) creates
            %   scatterplots of kpoly (scaled) vs parameter with specified xlabel
            %   xlab and, if type = "formin," with points above value 
            %   lab_limit labeled by formin
            %
            %   fig = EXPERIMENT.KPOLYPLOT(type,parameter,xlab,lab_limit,scale,1)
            %   creates and saves scatterplots of kpoly (scaled) vs parameter with specified xlabel
            %   xlab and, if type = "formin," with points above value 
            %   lab_limit labeled by formin
            %
            %   Inputs:
            %       type      : (string) whether to plot based on 'formin'
            %                   (kpolyplot) or 'PRM' (PRMplot)
            %       parameter : (string) x axes parameter, must be member
            %                   of PRM or Formin properties
            %       xlab      : (string) label for x axes
            %       lab_limit : (double) upper limit (in parameter) to
            %                   label points with formin name, only applies
            %                   if type of 'formin' (default is 1) 
            %       scale     : (string) kpoly axis scale (can be none,
            %                   log2, log10, ln) (default is "none")
            %       save      : (logical) whether or not to save figures
            %                   (default is false)
            %
            %   Uses kpolyplot or PRMplot based on type.
            % 
            % See also FORMIN, EXPERIMENT, KPOLYPLOT, PRMPLOT.
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
            %NTDPLOT create a scatterplot of change in polymerization vs. a
            %parameter for either all formins or all PRMs
            %
            %   fig = EXPERIMENT.NTDPLOT(type,parameter,xlab) creates
            %   scatterplots of change in kpoly vs parameter with specified xlabel
            %   xlab
            %
            %   fig = EXPERIMENT.NTDPLOT(type,parameter,xlab,scale) creates
            %   scatterplots of change in kpoly (scaled) vs parameter with specified xlabel
            %   xlab 
            %
            %   fig = EXPERIMENT.NTDPLOT(type,parameter,xlab,scale,1)
            %   creates and saves scatterplots of change in kpoly (scaled) vs parameter with specified xlabel
            %   xlab 
            %
            %   Inputs:
            %       type      : (string) whether to plot based on 'formin'
            %                   (NTDplot) or 'PRM' (PRM_NTD_plot)
            %       parameter : (string) x axes parameter, must be member
            %                   of PRM or Formin properties
            %       xlab      : (string) label for x axes
            %       scale     : (string) kpoly axis scale (can be none,
            %                   log2, log10, ln) (default is "log2")
            %       save      : (logical) whether or not to save figures
            %                   (default is false)
            %
            %   Uses PRM_NTD_plot or NTDplot based on type.
            % 
            % See also FORMIN, EXPERIMENT, PRM_NTD_PLOT, NTDPLOT.
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
            %FORMINBAR creates overview bargraphs for formins in ForminList
            %
            %   fig = EXPERIMENT.FORMINBAR creates overview bargraphs with
            %   kpoly axes not scaled and ratio axes log2 scaled
            %
            %   fig = EXPERIMENT.FORMINBAR(1) creates and saves overview bargraphs with
            %   kpoly axes not scaled and ratio axes log2 scaled
            %
            %   fig = EXPERIMENT.FORMINBAR(save, kpolyscale) creates and 
            %   (if save) saves overview bargraphs with kpoly axes scaled 
            %   as specified and ratio axes log2 scaled
            %
            %   fig = EXPERIMENT.FORMINBAR(save, kpolyscale, ratioscale) creates and 
            %   (if save) saves overview bargraphs with kpoly axes and
            %   ratio axes each scaled as specified 
            %
            %   Inputs:
            %       save      : (logical) whether or not to save figures
            %                   (default is false)
            %       kpolyscale: (string) kpoly axis scale (can be none,
            %                   log2, log10, ln) (default is "none")
            %       kpolyscale: (string) change in kpoly axis scale (can be
            %                   none, log2, log10, ln) (default is "log2")
            %   
            %   Makes 3 bargraphs:
            %       1. Single, double, and dimer polymerization rates side-by-side 
            %          for each formin
            %       2. Change upon NTD for each formin
            %       3. Number of PRMs for each formin (with PRMs)
            %
            %   Uses forminbar.
            % 
            % See also FORMIN, EXPERIMENT, FORMINBAR.
            arguments
                obj Experiment
                save logical=false
                kpolyscale string {mustBeMember(kpolyscale,{'none','log2','log10','ln'})}="none"
                ratioscale string {mustBeMember(ratioscale,{'none','log2','log10','ln'})}="log2"
            end
            fig=forminbar(obj.ForminList,obj.opts,save,kpolyscale=kpolyscale,ratioscale=ratioscale);
        end

        function fig=polymerstat_change_plot(obj,parameter,stat,xlab,ylab,scale,minus,save)
            %POLYMERSTAT_CHANGE_PLOT sreates 2 scatterplots of change in a polymer statistic vs. a parameter.
            %
            %   fig = EXPERIMENT.POLYMERSTAT_CHANGE_PLOT(parameter,stat,xlab,ylab) 
            %   creates 2 scatterplots of change in stat (scaled as log2) vs. parameter. 
            %
            %   fig = EXPERIMENT.POLYMERSTAT_CHANGE_PLOT(parameter,stat,xlab,ylab,scale) 
            %   creates 2 scatterplots of change in stat (scaled as specified) vs. parameter. 
            %
            %   fig = EXPERIMENT.POLYMERSTAT_CHANGE_PLOT(parameter,stat,xlab,ylab,scale,true) 
            %   creates 2 scatterplots of change in 1-stat (scaled as 
            %   specified) vs. parameter. 
            %
            %   fig = EXPERIMENT.POLYMERSTAT_CHANGE_PLOT(parameter,stat,xlab,ylab,scale,false,true) 
            %   creates and saves 2 scatterplots of change in stat (scaled as 
            %   specified) vs. parameter. 
            %
            %   fig = EXPERIMENT.POLYMERSTAT_CHANGE_PLOT(parameter,stat,xlab,ylab,scale,true,true) 
            %   creates and saves 2 scatterplots of change in 1-stat (scaled as 
            %   specified) vs. parameter. 
            %   
            %   Creates 2 scatterplots:
            %       1. labeled by submodel
            %       2. labeled by formin
            %
            %   if saving, it is recomended to set groot 'defaultfigureposition' to 
            %   [400 250 900 750] in order to avoid the figure being cut off when 
            %   saving as a pdf.
            %   
            %   Inputs:
            %       parameter   : (string) name of parameter containing values to plot against
            %                     stat (must be a property in the PRM class)
            %       stat        : (string)  name of polymerstat to plot (must be a property 
            %                     in the Lookuptable class)
            %       xlab        : (string) x-axis (parameter) label 
            %       ylab        : (string) y-axis (stat) label 
            %       scale       : y-axis scale (can be none, log2, log10,
            %                     ln) (default is log2)
            %       minus       : (logical) whether or not to use 1-stat 
            %                     values (deafult is false)
            %       save        : (logical) whether or not to save the plot 
            %                     to the results pdf in opts (deafult is false)
            %
            %   Uses polymerstat_change_plot.
            % 
            % See also FORMIN, EXPERIMENT, POLYMERSTAT_CHANGE_PLOT.
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
            %FORMINGRAPHIC creates figures with formin schematic and kpoly
            % bar graphs for each formin in ForminList
            %
            %   fig = EXPERIMENT.FORMINGRAPHIC creates figures with formin schematic and kpoly
            % bar graphs for each formin in ForminList
            %
            %   fig = EXPERIMENT.FORMINGRAPHIC(1) creates and saves figures with formin schematic and kpoly
            % bar graphs for each formin in ForminList
            %
            %   Uses formin.formingraphic.
            % 
            % See also FORMIN, EXPERIMENT, FORMIN/FORMINGRAPHIC.
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
            %EXPDATABAR creates side by side bargraphs of simulated and experimental 
            % rates.
            %   fig = EXPERIMENT.EXPDATABAR creates side by side
            %   bargraphs of simulated and experimental data in obj.data. Making
            %   all 6 plots.
            %
            %   fig = EXPERIMENT.EXPDATABAR(scale) creates side by side
            %   bargraphs of simulated and experimental data (scaled as specified) 
            %   in obj.data. Making all 6 plots. 
            %
            %   fig = EXPERIMENT.EXPDATABAR(scale,true) creates and 
            %   saves side by side bargraphs of simulated and experimental data 
            %   (scaled as specified) in obj.data. Making all 6 plots.
            %
            %   fig = EXPERIMENT.EXPDATABAR(group='grp') creates side by 
            %   side bargraphs of simulated and experimental data in obj.data 
            %   with the specified group. Makes 2 plots (see below) unless all 
            %   entries have type="ratio".
            %
            %   fig = EXPERIMENT.EXPDATABAR(scale,group='grp') creates 
            %   side by side bargraphs of simulated and experimental data (scaled
            %   as specified) in obj.data with the specified group. Makes 2 plots
            %   (see below) unless all entries have type="ratio".
            %
            %   fig = EXPERIMENT.EXPDATABAR(scale,true,group='grp')
            %   creates and saves side by side bargraphs of simulated and 
            %   experimental data (scaled as specified) in obj.data with the 
            %   specified group. Makes 2 plots (see below) unless all entries have type="ratio".
            %   
            %   Unless a group is specified, will generate plots for the following:
            %       1. All entries in input obj.data
            %       2. All entries in input obj.data with type="single"
            %       3. All entries in input obj.data with type="double"
            %       4. All entries in input obj.data with type="dimer"
            %       5. All entries in input obj.data with type="ratio"
            %       6. All entries in input obj.data with the same group (2 plots
            %          per group as below)
            %
            %   Will generate 2 plots for every "group" (unless the type="ratio")
            %       1. raw values of experiment and simulation
            %       2. raw experimental values and scaled simulation values
            %           (scaled so the smallest experimental value in the group is
            %           equal to the simulated value)
            %
            %   Inputs:
            %       scale      : method of scaling the polymerization values. Can
            %                    be 'none','log2','log10','ln' (default is 'none')
            %       save       : whether or not to save the plot to the results
            %                    pdf in settings (true/false); deafult is false
            %       group      : value for the "groups" property of objects in 
            %                    data to plot (Name-value argument)
            %
            %   Uses expdatabar.
            % 
            % See also FORMIN, EXPDATABAR.
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

        function fig=alldataplot(obj,scaledTF)
            % ALLDATAPLOT creates side by side bargraphs of simulated and
            % experimental rates for the groups "NTD data,""Fig 3 5," "Fig 3," "Fig
            % 4a," and "Fig 4c, and puts them all in one tiled layout
            %
            %   fig = EXPERIMENT.ALLDATAPLOT creates side by side bargraphs
            %   for the groups (does not include scaled plots generated by
            %   expdatabar)
            %
            %   fig = EXPERIMENT.ALLDATAPLOT(1) creates side by side bargraphs
            %   for the groups and includes the scaled plots as well
            % 
            % See also FORMIN, EXPERIMENT, EXPERIMENT/EXPDATABAR, EXPDATABAR.
            arguments
                obj Experiment
                scaledTF=0 %whether or not to include the scaled plots
            end

            fig=figure('units','centimeters','position',[5,5,45,30],'Name','Experimental Fits');hold on;
            if scaledTF
                tiles = tiledlayout(5,2,'TileSpacing','tight','Padding','none');
                p=9;
            else
                tiles = tiledlayout(2,3,'TileSpacing','tight','Padding','none');
                p=5;
            end
            title(tiles,num2str(obj.opts.getconsts))
        
            figNTD=obj.expdatabar("log2",false,group="NTD data"); % uses log2 scale
            ax_temp = copyobj(copyobj(figNTD(1).Children, get(figNTD(1).Children,'Type')), tiles);
            ax_temp(2).Layout.Tile=p;
            close(figNTD)
            p=p-1;

            fig3=obj.expdatabar(group="Fig 3");
            ax_temp = copyobj(copyobj(fig3(1).Children, get(fig3(1).Children,'Type')), tiles);
            ax_temp(2).Layout.Tile=p;
            if scaledTF
                p=p-1;
                ax_temp = copyobj(copyobj(fig3(2).Children, get(fig3(2).Children,'Type')), tiles);
                ax_temp(2).Layout.Tile=p;
            end
            close(fig3)
            p=p-1;
        
            fig35=obj.expdatabar(group="Fig 3 5");
            ax_temp = copyobj(copyobj(fig35(1).Children, get(fig35(1).Children,'Type')), tiles);
            ax_temp(2).Layout.Tile=p;
            if scaledTF
                p=p-1;
                ax_temp = copyobj(copyobj(fig35(2).Children, get(fig35(2).Children,'Type')), tiles);
                ax_temp(2).Layout.Tile=p;
            end
            close(fig35)
            p=p-1;
        
            fig4a=obj.expdatabar(group="Fig 4a");
            ax_temp = copyobj(copyobj(fig4a(1).Children, get(fig4a(1).Children,'Type')), tiles);
            ax_temp(2).Layout.Tile=p;
            if scaledTF
                p=p-1;
                ax_temp = copyobj(copyobj(fig4a(2).Children, get(fig4a(2).Children,'Type')), tiles);
                ax_temp(2).Layout.Tile=p;
            end
            close(fig4a)
            p=p-1;
        
            fig4c=obj.expdatabar(group="Fig 4c");
            ax_temp = copyobj(copyobj(fig4c(1).Children, get(fig4c(1).Children,'Type')), tiles);
            ax_temp(2).Layout.Tile=p;
            if scaledTF
                p=p-1;
                ax_temp = copyobj(copyobj(fig4c(2).Children, get(fig4c(2).Children,'Type')), tiles);
                ax_temp(2).Layout.Tile=p;
            end
            close(fig4c)
        end

        function makeresults(obj)
            % MAKERESULTS make the figures that were generated by the old
            % wholepackage code 
            % 
            % EXPERIMENT.MAKERESULTS
            % 
            % Creates and saves:
            %   - formin schematics for each formin
            %   - overview figures 
            %   - correlation plots per PRM and per formin
            %   - polymer stat plots
            %   - experimental data bargraphs
            %   - per formin and per PRM tables
            % 
            % See also EXPERIMENT, FORMIN, LOOKUPTABLE, EXPERIMENT/FORMINGRAPHIC,
            % EXPERIMENT/FORMINBAR, EXPERIMENT/NTDPLOT, EXPERIMENT/KPOLYPLOT, EXPERIMENT/POLYMERSTAT_CHANGE_PLOT,
            % EXPERIMENT/EXPDATABAR, EXPERIMENT/FORMINTABLE, EXPERIMENT/PRMTABLE.
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
            if ~isempty(Experiment1.data)
                obj.expdatabar("none",true);
            end

            % Make overview tables
            obj.formintable(true);
            obj.PRMtable(true);
        end

        function set_gating_file(obj,file)
            %SET_GATING_FILE set gating factors for formins according to an
            %input file
            %
            %   EXPERIMENT.SET_GATING_FILE(file)
            %
            %   Set gating factors according to .txt file with formin names
            %   followed by gating factors (delimated by whitespace)
            %
            % See also FORMIN, EXPERIMENT, EXPERIMENT/SET_GATING.
            filestruct = importdata(file); 
            forminlist = filestruct.textdata;
            for i = 1:length(forminlist)
                forminname = convertCharsToStrings(forminlist(i));   %takes the names (every other string)
                gating = filestruct.data(i);  
                obj.set_gating(forminname,gating)
            end
        end

        function set_gating(obj,forminname,gating)
            %SET_GATING set gating factor for formin
            %
            %   EXPERIMENT.SET_GATING(forminname, gating)
            %
            %   Inputs:
            %       forminname : (string) name of formin to set gating
            %                   factor, must be a name in ForminList
            %       gating     : (double) gating factor
            %
            % See also FORMIN, EXPERIMENT.
            arguments
                obj Experiment
                forminname string 
                gating double
            end
            formin=-1;
            indx=find(obj.forminnames==forminname);
            if indx
                formin=obj.ForminList(indx);
            else
                error("no formin in this Experiment has the name %s",forminname)
            end
            formin.gating=gating;
        end

        function out=SOSlist(obj,scaled,NameValueArgs)
            %SOSLIST compute individual squared errors for values in data based
            %on simulated kpolys of the corresponding formins
            %
            %   out= EXPERIMENT.SOSLIST compute SOS for all
            %   data points in obj.data
            %
            %   out= EXPERIMENT.SOSLIST(group='grp') compute SOS for all
            %   data points of the group type 'grp'
            %
            %   out= EXPERIMENT.SOSLIST(true, group='grp') compute SOS for all
            %   data points of the group type 'grp' and scales values by the smallest kpoly
            %   (no scaling if the type is ratio)
            %
            %   out= EXPERIMENT.SOSLIST(formin='form') compute the squared
            %   difference for the formin with name 'form'
            %
            %   out= EXPERIMENT.SOSLIST(...,recalc=false) compute SOS as
            %   specified using the last calculated kpoly instead of recalculating.
            %
            %   Note: the SOS value is set to 0 if its within error bars
            %   
            %   Inputs:
            %       scaled : (logical) whether to scale all points (only
            %               applies if group is specified, default is false)
            %       formin : (string) formin to compute SOS (NameValueArgs)
            %       group  : (string) group to compute SOS (NameValueArgs)
            %       recalc : (Formin) whether or not to reculate kpoly
            %               values (NameValueArgs, default is true)
            %
            % See also FORMIN, EXPERIMENT.
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
                % compute squared difference between experimental and
                % simulated kpoly, scaling the simulated value if specified
                % and setting val to 0 if within error bars
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
            %SOS compute sum of squared errors for values in data based
            %on simulated kpolys of the corresponding formins (sums output
            %of obj.SOSlist)
            %
            %   out= EXPERIMENT.SOS compute and sum SOS for all
            %   data points in obj.data
            %
            %   out= EXPERIMENT.SOS(group='grp') compute and sum SOS for all
            %   data points of the group type 'grp'
            %
            %   out= EXPERIMENT.SOS(true, group='grp') compute and sum  SOS for all
            %   data points of the group type 'grp' and scales values by the smallest kpoly
            %   (no scaling if the type is ratio)
            %
            %   out= EXPERIMENT.SOS(formin='form') compute the squared
            %   difference for the formin with name 'form'
            %
            %   out= EXPERIMENT.SOS(...,recalc=false) compute SOS as
            %   specified using the last calculated kpoly instead of recalculating.
            %
            %   Note: the SOS value is set to 0 if its within error bars
            %   
            %   Inputs:
            %       scaled : (logical) whether to scale all points (only
            %               applies if group is specified, default is false)
            %       formin : (string) formin to compute SOS (NameValueArgs)
            %       group  : (string) group to compute SOS (NameValueArgs)
            %       recalc : (Formin) whether or not to reculate kpoly
            %               values (NameValueArgs, default is true)
            %
            % See also FORMIN, EXPERIMENT, EXPERIMENT/SOSLIST.
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
            %RUNFMINSEARCH run fminsearch on SOS to fit parameters in obj.opts based
            %on data in obj.data
            %
            %   tab= EXPERIMENT.RUNFMINSEARCH(NameValueArgs)
            %
            %   Note: runs fminsearch on all 6 rate parameters, even if the
            %   3 state model is used (k_cap, k_del, r_cap, r_cap_exp,
            %   r_del, k_rel)
            %   
            %   Inputs:
            %       options     : optimset for fminsearch (default is optimset('MaxFunEvals',5000,'MaxIter',5000,'TolX',10^(-4),'TolFun',10^(-4)))
            %       groups      : groups to use for fit, if none specified, use all data in obj.data
            %       rbot        : lower bound on parameters in log10 space (default is -2)
            %       rtop        : upper bound on parameters in log10 space (default is 4)
            %       iterations  : number of fminsearch runs to run (default is 1)
            %       weight      : group to have SOS values multiplied by 10 (weighted more heavily in the fit)
            %       ogtable     : if true, just return a table with the current settings (default is false)
            %   
            %   Output is a table with entries generated by modifing 
            %   the results of opts.optionstable for the fminsearch 
            %   resulted set of options
            %
            % See also FORMIN, EXPERIMENT, EXPERIMENT/SOS, OPTIONS/OPTIONSTABLE.
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
        

            vals={obj.opts.k_cap,obj.opts.k_del,obj.opts.r_cap,obj.opts.r_cap_exp, obj.opts.r_del,obj.opts.k_rel};
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
                r=(rtop-rbot)*rand(1,6)+rbot;
                r=[vals{1},vals{2},vals{3},vals{4},vals{5},vals{6}];

                % not using scaled
                [params,fval] = fminsearch(@(z) getSOS(abs(10^(z(1))),abs(10^(z(2))),abs(10^(z(3))),abs(10^(z(4))),abs(10^(z(5))),abs(10^(z(6))),false),[r(1);r(2);r(3);r(4);r(5);r(6)],options); 
                vals=num2cell(abs(10.^(params)));
                T(index+1,:)=maketable();
                index=index+1;

                % using scaled
                [params,fval] = fminsearch(@(z) getSOS(abs(10^(z(1))),abs(10^(z(2))),abs(10^(z(3))),abs(10^(z(4))),abs(10^(z(5))),abs(10^(z(6))),true),[r(1);r(2);r(3);r(4);r(5);r(6)],options); 
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
            obj.opts.r_cap_exp=ogtable.r_cap_exp;
            
            function val=getSOS(k_cap,k_del,r_cap,r_cap_exp,r_del,k_rel,scaled)
                obj.opts.k_cap=k_cap;
                obj.opts.k_del=k_del;
                obj.opts.r_cap=r_cap;
                obj.opts.r_cap_exp=r_cap_exp;
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
                [tab.k_cap,tab.k_del,tab.r_cap,tab.r_cap_exp,tab.r_del,tab.k_rel]=vals{:};
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
                optstable=removevars(optstable,{'resultsdir','resultsfolder','python_path','k_cap','k_del','r_cap','r_cap_exp','r_del','k_rel','OriginalVariableNames'});
                tab=[tab,optstable];
            end
        end

        function applytable(obj,row)
            %APPLYTABLE update properties (formin and options) accoridng to
            %the input table
            %
            %   tab= EXPERIMENT.APPLYTABLE(row)
            %
            %   Inputs:
            %       row: table to apply
            %
            %   Searches for table variables with the same name as an
            %   Options property and applies them.
            %   Searches for table variables corresponding to any of the
            %   rate equations and runs obj.set_equation accordingly. These
            %   entries are assumed to have the name "rate_eq" and be of
            %   the format of obj.equationstext entries.
            %   Searches for table variables corresponding to formin
            %   properties and applies them. These entries are assumed to 
            %   have the name "forminname - property"
            %
            % See also FORMIN, EXPERIMENT, OPTIONS.
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

        function names = forminnames(obj)
            %FORMINNAMES returns array of the names of formins in
            %ForminList
            %
            %   names= EXPERIMENT.FORMINNAMES
            %
            % See also FORMIN, EXPERIMENT.
            L=length(obj.ForminList);
            names=[];
            for i=1:L
                names=[names;obj.ForminList(i).name];
            end
        end

        function fracs = fracmultibind(obj)
            %FRACMULTIBIND computes fraction of time with multiple
            %simultaneous binding for each formin in ForminList
            %
            %   fracs= EXPERIMENT.FRACMULTIBIND
            %
            %   Output in an Nx3 array where N is the number of formins. 
            %   A fraction is calculated for the single, double, and dimer
            %   filament types (hence the x3)
            %
            %   Uses pmultibind
            %
            % See also FORMIN, EXPERIMENT, PMULTIBIND.
            L=length(obj.ForminList);
            fracs=zeros(L,3);

            for i=1:L
                fract=pmultibind(obj.ForminList(i));
                fracs(i,1)=fract.single;
                fracs(i,2)=fract.double;
                fracs(i,3)=fract.dimer;
            end      
        end

        function [fig,fracs] = plotmultibind(obj)
            %PLOTMULTIBIND creates bargraph of the fraction of time spent
            %with multiple PRMs bound for single, double, and dimer
            %filaments for each formin in ForminList
            %
            %   fracs= EXPERIMENT.PLOTMULTIBIND
            %
            % See also FORMIN, EXPERIMENT, PMULTIBIND, EXPERIMENT/FRACMULTIBIND.
            fracs=obj.fracmultibind;
            fh1names=obj.forminnames;
            fig=figure;
            fracbind_table = table(fracs, 'RowNames', fh1names);
        
            
            fracbind_bar = bar(fracbind_table{:,:});
            set(gca,'xtick',1:length(fh1names), 'xticklabel',fracbind_table.Properties.RowNames)
            xtickangle(90)
            legend("single","double","dimer")
            xlabel('Formins')
            ylabel("Fraction of time spent with multiple simultaneous binding")
        
            title('Fraction of time spent with multiple simultaneous binding per Formin')
        end
    end
end