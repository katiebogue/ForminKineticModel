function fig = kpolyplot(forminList, parameter,xlab,lab_limit,settings,save)
% KPOLYPLOT  Creates and saves a scatterplot of polymerization vs. a
%            parameter. 
    %
    %   KPOLYPLOT(x,'xlab',N,d,settings) creates and saves a
    %   scatterplot of polymerization vs. x, labeling points beyond N
    %   
    %   Inputs:
    %       forminList   : array of Formins to gather data from
    %       parameter   : name of parameter containing values to plot against
    %                     polymerization
    %       xlab        : x-axis (parameter) label (string)
    %       lab_limit   : upper limit (in parameter) to label points with
    %                     formin name 
    %       settings    : Options class
    %       save        : whether or not to save the plot to the results
    %                     pdf in settings (true/false); deafult is false
    %
    %   Calls figuresave
    %
    %   See also FIGURESAVE, FORMIN, OPTIONS, EXPERIMENT.
    arguments
        forminList Formin
        parameter string
        xlab string
        lab_limit double
        settings Options
        save logical=false
    end

    if save
        set(groot,'defaultfigureposition',[400 250 900 750])
    end

    fig=figure;
    xdata=[forminList.(parameter)];
    kpoly=[forminList.kpoly];
    kpoly1_scatter = scatter(xdata,[kpoly.single], 'filled','o');
    hold on
    kpoly2_scatter = scatter(xdata,[kpoly.double], 'filled','s');
    hold on
    kpoly3_scatter = scatter(xdata,[kpoly.dimer], 'filled','p');
    xlabel(xlab)
    ylabel('(kpoly)')
    legend('Single', 'Double', 'N-Dimer', 'Location','northeastoutside');
    
    labels = [forminList.name];
    labelpoints(xdata,[kpoly.single],labels,'N',0.005, 1,'outliers_lim', {[-inf lab_limit; -inf inf]}, 'FontSize', 7)
    labelpoints(xdata,[kpoly.double],labels,'E',0.005, 1,'outliers_lim', {[-inf lab_limit; -inf inf]}, 'FontSize', 7)
    labelpoints(xdata,[kpoly.dimer],labels,'N',0.005, 1,'outliers_lim', {[-inf lab_limit; -inf inf]}, 'FontSize', 7)
    
    title(append('Polymerization Rates vs. ',xlab))
    
    if save
        figuresave(gcf,settings.pdf_name,append(gca().Title.String,'.fig'));
    end
end