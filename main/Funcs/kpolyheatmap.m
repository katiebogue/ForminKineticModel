function [fig,h]=kpolyheatmap(options,sweep_type,col,smooth,save)
% KPOLYHEATMAP  ...
    %
    %   KPOLYHEATMAP() ...
    %
    %   Inputs:
    %       options  : 
    %
    %
    %   See also .

arguments
    options Options
    sweep_type string
    col double
    smooth logical=true
    save logical=false
end
fig=figure;

if save
    set(groot,'defaultfigureposition',[400 250 900 750]) % helps prevent cut offs in figs
end


kpoly_lab="log_{2}(k_{poly} N terminal dimerized/k_{poly} double)";

if sweep_type=="NT dist v r_cap"
    param="r_cap";
    x_label="Distance from PRM to N-term";
    scatn=0;
    PRM_loc=400;
    PRM_size=10;

    PRM_lab=strcat("loc: ",num2str(PRM_loc)," size: ",num2str(PRM_size));

    ogopt=options.(param);
    param_vals=10.^([-10:0.5:20]);
    xvals=PRM_loc:600;

    kpolyratios=zeros(length(param_vals)*length(xvals),1);
    paramvals_matrix=kpolyratios;
    xvals_matrix=kpolyratios;

    index=0;

    for i=1:length(xvals)
        FH1_length=xvals(i);
        formin1=Formin("formin1",options,c_PA=1,gating=1,length=FH1_length,PRMloc=PRM_loc,PRMsize=PRM_size);
        for n=1:length(param_vals)
            options.(param)=param_vals(n);
            index=index+1;
            kpolyratios(index,1)=formin1.kpoly.ratio;
            paramvals_matrix(index,1)=param_vals(n);
            xvals_matrix(index,1)=xvals(i)-PRM_loc;
        end
    end
    options.(param)=ogopt;

    if smooth
        for n=1:length(param_vals)
            index=paramvals_matrix==param_vals(n);
            kpolyratios(index)=smoothdata(kpolyratios(index),"lowess","SmoothingFactor",.58);
        end
        kpoly_lab=strcat(kpoly_lab," (smoothed)");
    end

    tbl=table(xvals_matrix,log10(paramvals_matrix),log2(kpolyratios));

    h = heatmap(tbl,'xvals_matrix','Var2','ColorVariable','Var3');
    h.XLabel = x_label;
    h.ColorMethod = 'none';
    h.NodeChildren(3).YDir='normal'; 
    if col==1
        h.Colormap=[jet];
    elseif col==2
        h.Colormap=[parula];
        
        if not(smooth)
            red = [1 0 0];
            h.Colormap=[parula;red];
            h.ColorLimits = [min(log2(kpolyratios)) 0];
        end
    end

    strchar=char(param);
        if strchar(1:2)=="k_" || strchar(1:2)=="r_"
            param=strcat(strchar(1:2),"{",strchar(3:end),"}");
        end
    param_lab=strcat("log_{10}(",param,")");
    h.YLabel = param_lab;

    prvec=getmeanstat(options.lookup,'Prvec0',PRM_loc);
    pocc=getmeanstat(options.lookup,'POcclude',PRM_loc);
    pocc0=getmeanstat(options.lookup,'POcclude',1);

    kcapavg=log10(options.k_cap*(1-pocc));
    kdelavg=log10(options.k_del*(1.0e33*(prvec)/(27*6.022e23))*(1-pocc0));

    lab=strcat("kcap: ",num2str(kcapavg)," kdel: ",num2str(kdelavg));

    h.Title = {kpoly_lab,PRM_lab,lab};

    % Convert each number in the array into a string
    CustomXLabels = string(xvals);
    CustomYLabels = string(log10(param_vals));
    % Replace all but the fifth elements by spaces
    CustomXLabels(mod(xvals,5) ~= 0) = " ";
    CustomYLabels(mod(log10(param_vals),2) ~= 0) = " ";
    % Set the 'XDisplayLabels' property of the heatmap 
    % object 'h' to the custom x-axis tick labels
    h.XDisplayLabels = CustomXLabels;
    h.YDisplayLabels = CustomYLabels;
end

if save
    figuresave(gcf,options,append("heatmapsweep_",sweep_type,'.fig'));
end
end

function value=getmeanstat(lookup,stat,isite)
    a=lookup.getstat(type='double',Stat=stat,iSite=isite,fil='a');
    a=structfun(@mean,a);
    a=mean(a);

    b=lookup.getstat(type='double',Stat=stat,iSite=isite,fil='b');
    b=structfun(@mean,b);
    b=mean(b);

    dob=mean([a,b]);

    a=lookup.getstat(type='dimer',Stat=stat,iSite=isite,fil='a');
    a=structfun(@mean,a);
    a=mean(a);

    b=lookup.getstat(type='dimer',Stat=stat,iSite=isite,fil='b');
    b=structfun(@mean,b);
    b=mean(b);

    dim=mean([a,b]);

    value=mean([dim,dob]);
end