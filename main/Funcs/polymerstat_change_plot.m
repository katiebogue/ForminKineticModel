function fig=polymerstat_change_plot(forminList,parameter,stat,xlab,ylab,settings,save)
% POLYMERSTAT_CHANGE_PLOT  Creates and saves a scatterplot of change in polymer
%                  statistics vs. a parameter.
    %
    %   POLYMERSTAT_CHANGE_PLOT(x,'xlab',Data,settings) creates and saves a
    %   scatterplot of change in polymer statistics vs. x. 
    %   
    %   Plots change in 3 polymer stats:
    %       1. pr_0
    %       2. 1-pocc
    %       3. 1-pocc_0
    %
    %   Inputs:
    %       forminList  : array of Formins to gather data from
    %       parameter   : name of parameter containing values to plot against
    %                     polymerstat
    %       stat        : name of polymerstat to plot
    %       xlab        : x-axis (parameter) label (string) 
    %       ylab        : y-axis (stat) label (string)
    %       ettings    : Options class
    %       save        : whether or not to save the plot to the results
    %                     pdf in settings (true/false); deafult is false
    %
    %   Calls figuresave
    %
    %   See also .

    arguments
        forminList Formin
        parameter string
        stat string
        xlab string
        ylab string
        settings Options
        save logical=false
    end

    PRMList=[forminList.PRMList];
    xdata=[PRMList.(parameter)];

    ydata=[PRMList.getprop((stat))];
    ydataratio=[ydata.ratio];
    ya=[ydataratio.a];
    yb=[ydataratio.b];
    
    
    ylab=append(ylab,' Dimer/Double');
    title_=({append('Change in ',ylab),append(' with NTD vs. ',xlab,' per individual PRM')});
    titlefig=(append('Change in ',ylab,' with NTD vs. ',xlab,' per individual PRM'));

    %% Plot based on each filament
    fig(1)=figure;
    scatter_3a_2a = scatter(xdata,ya,'filled', 'o');
    hold on
    scatter_3b_2b = scatter(xdata,yb, 'filled', 's');
    hold on
    
    xlabel(xlab)
    ylabel(ylab)
    legend('filament a', 'filament b');
    
    title(title_)
     
    if save
        figuresave(gcf,settings,append('stats.',titlefig,'.fig'));
    end
    
    %% Plot based on formin
    fig(2)=figure;
    
    for i=1:length(forminList)
        x=[forminList(i).PRMList.(parameter)];
        y=[forminList(i).PRMList.getprop((stat))];
        y=[y.ratio];
        y=[y.a,y.b];
        x=[x,x];
        kpoly3a_2a_individual_scatter = scatter(x,y,'filled',settings.shapes(i),'MarkerFaceColor', settings.colors(i),'MarkerEdgeColor','k');
        hold on
    end
    
    xlabel(xlab)
    ylabel(ylab)
    
    legend([forminList.name], 'Location','southoutside','NumColumns',6,'FontSize',8);
    
    title(title_)
     
    if save
        figuresave(gcf,settings,append('stats.',titlefig,'.fig'));
        savefig(h,fullfile("figures",append('stats.',titlefig,'.fig')));
    end
end