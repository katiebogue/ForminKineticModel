function [fig,h]=kpolyheatmap(options,sweep_type,col,smooth,smoothfac,save)
% KPOLYHEATMAP  creates heatmap of kpoly ratios for a fictitious formin
% construct swept across formin and/or rate parameters as specified by the
% input
    %
    %   [fig,h]=KPOLYHEATMAP(options,sweep_type,col) creates heatmap of sweep
    %   specified by sweep_type with the specified colormap options, uses
    %   smoothing factor of 0.58
    %
    %   [fig,h]=KPOLYHEATMAP(options,sweep_type,col,0) creates heatmap of sweep
    %   specified by sweep_type with the specified colormap options, with
    %   no smoothing
    %
    %   [fig,h]=KPOLYHEATMAP(options,sweep_type,col,1,smoothfac) creates heatmap of
    %   sweep specified by sweep_type with the specified colormap options, 
    %   uses specified smoothing factor 
    %
    %   [fig,h]=KPOLYHEATMAP(options,sweep_type,col,0,x,1) creates and saves
    %   heatmap of sweep specified by sweep_type with the specified
    %   colormap options, with no smoothing
    %
    %   [fig,h]=KPOLYHEATMAP(options,sweep_type,col,1,smoothfac,1) creates and
    %   saves heatmap of sweep specified by sweep_type with the specified
    %   colormap options, uses specified smoothing factor
    %
    %   Inputs:
    %       options    : Options object with parameters and rate equations
    %       sweep_type : (string) the type of sweep to run (currently only
    %                    set to work with "NT dist v r_cap")
    %       col        : (double) colormap type
    %                       1- use jet
    %                       2- use parula, if smooth not true, set all increasing value to red
    %                       3- load 'customcolorbar_red_blue.mat' and set color limits to min values
    %                       [x, x] - use 'customcolorbar_red_blue.mat' and
    %                               set colorlimits to [x, x]
    %       smooth     : (logical) whether or not to use smoothing on the
    %                     heatmap (default is true)
    %       smoothfac  : (double) smoothing factor, applied if smooth=true
    %                    (default is 0.58)
    %       save       : (logical) whether or not to save the resulting
    %                    heatmap (default is false)
    %
    %   Outputs:
    %       fig : Figure with heatmap
    %       h   : heatmap object
    %
    %   Currently only works with sweeps across NT distance and r_cap
    %
    %   Calls figuresave
    %   Loads customcolorbar_red_blue.mat
    % 
    %   See also FORMIN, OPTIONS, PRM, FIGURESAVE, FILAMENTSCHEMATIC, FORMINKPOLYBAR, KPOLYMERIZATION.

arguments
    options Options
    sweep_type string
    col double
    smooth logical=true
    smoothfac double =0.58
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
    PRM_loc=100;
    PRM_size=10;

    PRM_lab=strcat("loc: ",num2str(PRM_loc)," size: ",num2str(PRM_size));

    ogopt=options.(param);
    param_vals=10.^([-4:0.5:16]);
    xvals=PRM_loc:400;

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
            vals=kpolyratios(index);
            smootheddata=smoothdata(vals,"lowess","SmoothingFactor",smoothfac);
            smootheddata(smootheddata<0)=vals(smootheddata<0);
            kpolyratios(index)=smootheddata;
        end
        kpoly_lab=strcat(kpoly_lab," (smoothed:",num2str(smoothfac),")");
    end

    tbl=table(xvals_matrix,log10(paramvals_matrix),log2(kpolyratios));

    h = heatmap(tbl,'xvals_matrix','Var2','ColorVariable','Var3');
    h.XLabel = x_label;
    h.ColorMethod = 'none';
    h.GridVisible="off";
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
    elseif col==3
        load('customcolorbar_red_blue.mat');
        h.Colormap=CustomColormap;
        minn=min(log2(kpolyratios(kpolyratios~=0)));
        h.ColorLimits=[minn abs(minn)];
    elseif length(col)==2
        load('customcolorbar_red_blue.mat');
        h.Colormap=CustomColormap;
        h.ColorLimits=col;
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
    CustomXLabels = string(xvals-PRM_loc);
    CustomYLabels = string(log10(param_vals));
    % Replace all but the fifth elements by spaces
    CustomXLabels(mod(xvals-PRM_loc,20) ~= 0) = " ";
    CustomYLabels(mod(log10(param_vals),4) ~= 0) = " ";
    % Set the 'XDisplayLabels' property of the heatmap 
    % object 'h' to the custom x-axis tick labels
    h.XDisplayLabels = CustomXLabels;
    h.YDisplayLabels = CustomYLabels;
    s = struct(h); 
    s.XAxis.TickLabelRotation = 0;   % horizontal
end

if save
    figuresave(gcf,options,append("heatmapsweep_",sweep_type," ", num2str(options.k_cap)," kdel ", num2str(options.k_del),'.fig'),true);
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