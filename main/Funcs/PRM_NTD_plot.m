function fig = PRM_NTD_plot(forminList,parameter,xlab,settings,save)
% PRM_NTD_plot  Creates and saves 2 scatterplots of change in
% polymerization vs. a parameter for each PRM.
    %
    %   PRM_NTD_PLOT(x,'xlab',Data,settings) creates and saves 2
    %   scatterplots of change in polymerization vs. x for each PRM
    %   
    %   Creates 2 scatterplots:
    %       1. labeled by submodel
    %       2. labeled by formin
    %
    %   Inputs:
    %       forminList  : array of Formins to gather data from
    %       parameter   : name of parameter containing values to plot against
    %                     polymerization
    %       xlab        : x-axis (parameter) label (string)
    %       settings    : Options class
    %       save        : whether or not to save the plot to the results
    %                     pdf in settings (true/false); deafult is false
    %   
    %   Calls figuresave but only saves one .fig file for both figures
    %
    %   See also .

    arguments
        forminList Formin
        parameter string
        xlab string
        settings Options
        save logical=false
    end
    
    title_=append('Change in Polymerization Rates vs. ',xlab,' per individual PRM');
    %% Plot based on each filament
    fig(1)=figure;
    PRMList=[forminList.PRMList];
    xdata=[PRMList.(parameter)];
    kpoly=[PRMList.kpoly];
    ratio=[kpoly.ratio];

    kpoly3a_2a_individual_scatter_length = scatter(xdata,[ratio.a],'filled', 'o');
    hold on
    kpoly3b_2b_individual_scatter_length = scatter(xdata,[ratio.b], 'filled', 's');
    hold on
    
    xlabel(xlab)
    ylabel('(kpoly Dimer/ kpoly Double)')
    legend('filament a', 'filament b');
    
    title(title_)
    
    if save
        figuresave(gcf,settings,append(gca().Title.String,'.fig'));
    end
    
    
    %% Plot based on formin
    fig(2)=figure;

    for i=1:length(forminList)
        x=[forminList(i).PRMList.(parameter)];
        y=[forminList(i).PRMList.kpoly];
        y=[y.ratio];
        y=[y.a,y.b];
        x=[x,x];
        kpoly3a_2a_individual_scatter = scatter(x,y,'filled',settings.shapes(i),'MarkerFaceColor', settings.colors(i),'MarkerEdgeColor','k');
        hold on
    end
    
    xlabel(xlab)
    ylabel('(kpoly Dimer/ kpoly Double)')
    
    legend([forminList.name], 'Location','northeastoutside');
    
    title(title_)
    
    if save
        figuresave(gcf,settings,append(gca().Title.String,'.fig'));
        savefig(fig,fullfile("figures",append(gca().Title.String,'.fig')));
    end
    

end