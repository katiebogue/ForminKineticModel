function fig = expdatabar(datastruct,settings,scale,save,NameValueArgs)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
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