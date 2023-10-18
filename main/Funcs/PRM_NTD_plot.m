function fig = PRM_NTD_plot(forminList,parameter,xlab,scale,settings,save)
% PRM_NTD_plot  Creates 2 scatterplots of change in polymerization vs. a 
% parameter for each PRM.
    %
    %   fig = PRM_NTD_PLOT(forminList, x,'xlab','scale',settings) creates 2
    %   scatterplots of change in polymerization (scaled as specified) vs. x. 
    %
    %   fig = PRM_NTD_PLOT(forminList, x,'xlab','scale',settings,true) 
    %   creates and saves 2 scatterplots of change in polymerization 
    %   (scaled as specified) vs. x. 
    %   
    %   if saving, it is recomended to set groot 'defaultfigureposition' to 
    %   [400 250 900 750] in order to avoid the figure being cut off when 
    %   saving as a pdf.
    %
    %   Creates 2 scatterplots:
    %       1. labeled by submodel
    %       2. labeled by formin
    %
    %   Inputs:
    %       forminList  : array of Formins to gather data from
    %       parameter   : name of parameter containing values to plot against
    %                     polymerization (must be a property in the PRM
    %                     class)
    %       xlab        : x-axis (parameter) label (string)
    %       scale       : y-axis scale (can be none, log2, log10, ln)
    %       settings    : Options class
    %       save        : whether or not to save the plot to the results
    %                     pdf in settings (true/false); deafult is false
    %   
    %   Calls figuresave
    %
    %   See also FIGURESAVE, FORMIN, PRM, EXPERIMENT, OPTIONS, NTDPLOT, KPOLYMERIZATION.

    arguments
        forminList Formin
        parameter string
        xlab string
        scale string {mustBeMember(scale,{'none','log2','log10','ln'})} % scale for ratio
        settings Options
        save logical=false
    end

    if scale=="none"
        yscale=@(x) x;
        ylab="k_{poly} N terminal dimerized/k_{poly} double";
    elseif scale=="log2"
        ylab="log_{2}(k_{poly} N terminal dimerized/k_{poly} double)";
        yscale=@(x) log2(x);
    elseif scale=="log10"
        ylab="log_{10}(k_{poly} N terminal dimerized/k_{poly} double)";
        yscale=@(x) log10(x);
    elseif scale=="ln"
        ylab="ln(k_{poly} N terminal dimerized/k_{poly} double)";
        yscale=@(x) log(x);
    end
    
    title_=append('Change in Polymerization Rates vs. ',xlab,' per individual PRM');
    %% Plot based on each filament
    fig(1)=figure;
    PRMList=[forminList.PRMList];
    xdata=[PRMList.(parameter)];
    kpoly=[PRMList.kpoly];
    ratio=[kpoly.ratio];

    kpoly3a_2a_individual_scatter_length = scatter(xdata,yscale([ratio.a]),'filled', 'o');
    hold on
    kpoly3b_2b_individual_scatter_length = scatter(xdata,yscale([ratio.b]), 'filled', 's');
    hold on

    lne=yline(yscale(1),'LineWidth',2);
    lne.Annotation.LegendInformation.IconDisplayStyle='off';
    
    xlabel(xlab)
    ylabel(ylab)
    legend('filament a', 'filament b');
    
    title(title_)
    
    if save
        figuresave(gcf,settings,append(gca().Title.String,'per fil.fig'));
    end
    
    
    %% Plot based on formin
    fig(2)=figure;

    for i=1:length(forminList)
        x=[forminList(i).PRMList.(parameter)];
        y=[forminList(i).PRMList.kpoly];
        y=[y.ratio];
        y=[y.a,y.b];
        x=[x,x];
        kpoly3a_2a_individual_scatter = scatter(x,yscale(y),'filled',settings.shapes(i),'MarkerFaceColor', settings.colors(i),'MarkerEdgeColor','k');
        hold on
    end

    lne=yline(yscale(1),'LineWidth',2);
    lne.Annotation.LegendInformation.IconDisplayStyle='off';
    
    xlabel(xlab)
    ylabel(ylab)
    
    legend([forminList.name], 'Location','northeastoutside');
    
    title(title_)
    
    if save
        figuresave(gcf,settings,append(gca().Title.String,'per formin.fig'));
    end
    

end