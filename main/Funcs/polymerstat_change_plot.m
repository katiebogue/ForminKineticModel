function fig=polymerstat_change_plot(forminList,parameter,stat,xlab,ylab,scale,settings,minus,save)
% POLYMERSTAT_CHANGE_PLOT  Creates 2 scatterplots of change in a polymer 
% statistic vs. a parameter.
    %
    %   fig = POLYMERSTAT_CHANGE_PLOT(forminList,x,stat,'xlab','ylab',scale,settings) 
    %   creates 2 scatterplots of change in stat (scaled as specified) vs. x. 
    %
    %   fig = POLYMERSTAT_CHANGE_PLOT(forminList,x,stat,'xlab','ylab',scale,settings,true) 
    %   creates 2 scatterplots of change in 1-stat (scaled as specified) vs. 
    %   x. 
    %
    %   fig = POLYMERSTAT_CHANGE_PLOT(forminList,x,stat,'xlab','ylab',scale,settings,true,true) 
    %   creates and saves 2 scatterplots of change in 1-stat (scaled as 
    %   specified) vs. x. 
    %
    %   fig = POLYMERSTAT_CHANGE_PLOT(forminList,x,stat,'xlab','ylab',scale,settings,false,true) 
    %   creates and saves 2 scatterplots of change in stat (scaled as 
    %   specified) vs. x. 
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
    %       forminList  : array of Formins to gather data from
    %       parameter   : name of parameter containing values to plot against
    %                     stat (must be a property in the PRM class)
    %       stat        : name of polymerstat to plot (must be a property 
    %                     in the Lookuptable class)
    %       xlab        : x-axis (parameter) label (string) 
    %       ylab        : y-axis (stat) label (string)
    %       scale       : y-axis scale (can be none, log2, log10, ln)
    %       settings    : Options class
    %       minus       : whether or not to use 1-stat values (true/false); 
    %                     deafult is false
    %       save        : whether or not to save the plot to the results
    %                     pdf in settings (true/false); deafult is false
    %
    %   Calls figuresave
    %
    %   See also FIGURESAVE, FORMIN, PRM, EXPERIMENT, OPTIONS, LOOKUPTABLE.

    arguments
        forminList Formin
        parameter string
        stat string
        xlab string
        ylab string
        scale {mustBeMember(scale,{'none','log2','log10','ln'})}
        settings Options
        minus logical=false
        save logical=false
    end

    titlefig=(append('Change in ',ylab,' with NTD vs. ',xlab,' per individual PRM'));
    if minus
        title_=({append('Change in 1- ',ylab),append(' with NTD vs. ',xlab,' per individual PRM')});
    else
        title_=({append('Change in ',ylab),append(' with NTD vs. ',xlab,' per individual PRM')});
    end

    if scale=="none"
        yscale=@(x) x;
        ylab=strcat(ylab," N terminal dimerized/double");
    elseif scale=="log2"
        ylab=strcat("log_{2}(",ylab," N terminal dimerized/double)");
        yscale=@(x) log2(x);
    elseif scale=="log10"
        ylab=strcat("log_{10}(",ylab," N terminal dimerized/double)");
        yscale=@(x) log10(x);
    elseif scale=="ln"
        ylab=strcat("ln(",ylab," N terminal dimerized/double)");
        yscale=@(x) log(x);
    end

    if minus
        ylab=append("1- ",ylab);
    end

    PRMList=[forminList.PRMList];
    xdata=[PRMList.(parameter)];

    ydata=[PRMList.getprop((stat))];
    if minus
        ydouble=[ydata.double];
        ydimer=[ydata.dimer];
        yadouble=1-[ydouble.a];
        ybdouble=1-[ydouble.b];
        yadimer=1-[ydimer.a];
        ybdimer=1-[ydimer.b];
        ya=yadimer./yadouble;
        yb=ybdimer./ybdouble;
    else
        ydataratio=[ydata.ratio];
        ya=[ydataratio.a];
        yb=[ydataratio.b];
    end
    

    %% Plot based on each filament
    fig(1)=figure;
    scatter_3a_2a = scatter(xdata,yscale(ya),'filled', 'o');
    hold on
    scatter_3b_2b = scatter(xdata,yscale(yb), 'filled', 's');
    hold on

    lne=yline(yscale(1),'LineWidth',2);
    lne.Annotation.LegendInformation.IconDisplayStyle='off';
    
    xlabel(xlab)
    ylabel(ylab)
    legend('filament a', 'filament b');
    
    title(title_)
     
    if save
        figuresave(gcf,settings,append('stats.',titlefig,'per fil.fig'));
    end
    
    %% Plot based on formin
    fig(2)=figure;
    
    for i=1:length(forminList)
        x=[forminList(i).PRMList.(parameter)];
        y=[forminList(i).PRMList.getprop((stat))];
        if minus
            ydouble=[y.double];
            ydimer=[y.dimer];
            yadouble=1-[ydouble.a];
            ybdouble=1-[ydouble.b];
            yadimer=1-[ydimer.a];
            ybdimer=1-[ydimer.b];
            ya=yadimer./yadouble;
            yb=ybdimer./ybdouble;
            y=[ya,yb];
        else
            y=[y.ratio];
            y=[y.a,y.b];
        end
        x=[x,x];
        kpoly3a_2a_individual_scatter = scatter(x,yscale(y),'filled',settings.shapes(i),'MarkerFaceColor', settings.colors(i),'MarkerEdgeColor','k');
        hold on
    end

    lne=yline(yscale(1),'LineWidth',2);
    lne.Annotation.LegendInformation.IconDisplayStyle='off';
    
    xlabel(xlab)
    ylabel(ylab)
    
    legend([forminList.name], 'Location','southoutside','NumColumns',6,'FontSize',8);
    
    title(title_)
     
    if save
        figuresave(gcf,settings,append('stats.',titlefig,'per formin.fig'));
    end
end