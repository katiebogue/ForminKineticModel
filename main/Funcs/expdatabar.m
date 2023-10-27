function fig = expdatabar(datastruct,settings,scale,save,NameValueArgs)
% EXPDATABAR creates side by side bargraphs of simulated and experimental 
% rates.
    %   fig = EXPDATABAR(datastruct,settings) creates side by side
    %   bargraphs of simulated and experimental data in datastruct. Making
    %   all 6 plots.
    %
    %   fig = EXPDATABAR(datastruct,settings,'scale') creates side by side
    %   bargraphs of simulated and experimental data (scaled as specified) 
    %   in datastruct. Making all 6 plots. 
    %
    %   fig = EXPDATABAR(datastruct,settings,'scale',true) creates and 
    %   saves side by side bargraphs of simulated and experimental data 
    %   (scaled as specified) in datastruct. Making all 6 plots.
    %
    %   fig = EXPDATABAR(datastruct,settings,group='grp') creates side by 
    %   side bargraphs of simulated and experimental data in datastruct 
    %   with the specified group. Makes 2 plots (see below) unless all 
    %   entries have type="ratio".
    %
    %   fig = EXPDATABAR(datastruct,settings,'scale',group='grp') creates 
    %   side by side bargraphs of simulated and experimental data (scaled
    %   as specified) in datastruct with the specified group. Makes 2 plots
    %   (see below) unless all entries have type="ratio".
    %
    %   fig = EXPDATABAR(datastruct,settings,'scale',true,group='grp')
    %   creates and saves side by side bargraphs of simulated and 
    %   experimental data (scaled as specified) in datastruct with the 
    %   specified group. Makes 2 plots (see below) unless all entries have type="ratio".
    %   
    %   Unless a group is specified, will generate plots for the following:
    %       1. All entries in input datastruct
    %       2. All entries in input datastruct with type="single"
    %       3. All entries in input datastruct with type="double"
    %       4. All entries in input datastruct with type="dimer"
    %       5. All entries in input datastruct with type="ratio"
    %       6. All entries in input datastruct with the same group (2 plots
    %          per group as below)
    %
    %   Will generate 2 plots for every "group" (unless the type="ratio")
    %       1. raw values of experiment and simulation
    %       2. raw experimental values and scaled simulation values
    %           (scaled so the smallest experimental value in the group is
    %           equal to the simulated value)
    %
    %   Calls figuresave
    %
    %   Inputs:
    %       datastruct : structure array with properties type, formin,
    %                    value, groups, errtop, and errbot
    %       settings   : Options class
    %       scale      : method of scaling the polymerization values. Can
    %                    be 'none','log2','log10','ln' (default is 'none')
    %       save       : whether or not to save the plot to the results
    %                    pdf in settings (true/false); deafult is false
    %       group      : value for the "groups" property of objects in 
    %                    datastruct to plot (Name-value argument)
    %   
    %   See also EXPERIMENT, FORMIN, PRM, OPTIONS, FIGURESAVE.
    arguments
        datastruct struct
        settings Options
        scale string {mustBeMember(scale,{'none','log2','log10','ln'})}="none"
        save logical=false
        NameValueArgs.group string
    end

    if scale=="none"
        yscale=@(x) x;
        ylab="";
    elseif scale=="log2"
        ylab="log_{2}";
        yscale=@(x) log2(x);
    elseif scale=="log10"
        ylab="log_{10}";
        yscale=@(x) log10(x);
    elseif scale=="ln"
        ylab="ln";
        yscale=@(x) log(x);
    end
    
    formins=[datastruct.formin];
    kpolys=[formins.kpoly];
    types=[datastruct.type];

    simdata=zeros(1,length(datastruct));
    
    for i=1:length(datastruct)
        type=types(i);
        simdata(i)=kpolys(i).(type);
    end

    
    fig=[];

    if isfield(NameValueArgs,"group")
        makegroupplot(NameValueArgs.group);
        return
    end
    makeplot("all",datastruct,simdata,1);
    makeplot("single",datastruct,simdata,1);
    makeplot("double",datastruct,simdata,1);
    makeplot("dimer",datastruct,simdata,1);
    makeplot("ratio",datastruct,simdata,1);

    groups=unique([datastruct.groups]);
    for i=1:length(groups)
        makegroupplot(groups(i));
    end

    function makegroupplot(group)
        newdatastruct=datastruct;
        newsimdata=simdata;
        for j=length(newdatastruct):-1:1
            if ~ismember(group,newdatastruct(j).groups)
                newdatastruct(j)=[];
                newsimdata(j)=[];
            end
        end
        expdata=[newdatastruct.value];
        minloc=expdata==min(expdata);
        scaler=newsimdata(minloc)/min(expdata);
        makescaler=true;

        newsave=false;
        if save
            save=false;
            newsave=true;
        end
        makeplot("all",newdatastruct,newsimdata,1);
        if length(unique([newdatastruct.type]))==1
            type=newdatastruct(1).type;
            ylabel(strcat(ylab," K_{poly} ",type));
            if type=="ratio"
                makescaler=false;
            end
        end
        title(group)
        if newsave
            figuresave(gcf,settings,append('experimental data bar--',group,'.fig'));
        end
        
        if makescaler
            makeplot("all",newdatastruct,newsimdata,scaler);
            if length(unique([newdatastruct.type]))==1
                type=newdatastruct(1).type;
                ylabel(strcat(ylab," K_{poly} ",type," scaled"));
            end
            title(strcat(group," scaled"))
            if newsave
                figuresave(gcf,settings,append('experimental data bar--',group,'.fig'));
            end
        end
    end
    function makeplot(type,datastruct,simdata,scaler)
        if type=="all"
            typesimdata=simdata;
            typevalues=datastruct;
        else
            typesimdata=simdata(types==type);
            typevalues=datastruct(types==type);
        end

        if length(typesimdata)==0
            return
        end
        
        newfig=figure;
        
        if length(typesimdata)==1
            kpolybar=bar(1,[yscale([typevalues.value]'),yscale(typesimdata'./scaler)]);
        else
            kpolybar=bar([yscale([typevalues.value]'),yscale(typesimdata'./scaler)]);
        end

        typeformins=[typevalues.formin];
    
        set(gca,'xticklabel',[typeformins.name]);
        hold on
        
        erlocs=(0:(length(typesimdata)-1))+0.85;

        if scale=="none"
            ytop=[typevalues.errtop];
            ybot=[typevalues.errbot];
        else
            ytop=yscale(exp(1))*([typevalues.errtop]./[typevalues.value]);
            ybot=yscale(exp(1))*([typevalues.errbot]./[typevalues.value]);
        end


        er = errorbar(erlocs,yscale([typevalues.value]'), ybot,ytop);   
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';
        
        set(kpolybar(1), 'FaceColor','b');
        set(kpolybar(2), 'FaceColor','r');
        legend( 'Experimental', 'model');
        xlabel('Formins');
        if scaler==1
            ylabel(strcat(ylab," K_{poly} ",type));
        else
            ylabel(strcat(ylab," K_{poly} ",type," scaled"));
        end

        if save
            figuresave(gcf,settings,append('experimental data bar--',type,'.fig'));
        end

        fig=[fig,newfig];
    end
end