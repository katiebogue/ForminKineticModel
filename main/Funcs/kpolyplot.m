function fig = kpolyplot(forminList, parameter,xlab,lab_limit,scale,settings,save)
% KPOLYPLOT  Creates a scatterplot of polymerization vs. a parameter. 
    %
    %   fig = KPOLYPLOT(formins,x,'xlab',N,'scale',settings) creates a
    %   scatterplot of polymerization (scaled as specified) vs. x, labeling 
    %   points with x values beyond N.
    %
    %   fig = KPOLYPLOT(formins,x,'xlab',N,'scale',settings,true) creates 
    %   and saves a scatterplot of polymerization (scaled as specified) 
    %   vs. x, labeling points with x values beyond N.
    %
    %   if saving, it is recomended to set groot 'defaultfigureposition' to 
    %   [400 250 900 750] in order to avoid the figure being cut off when 
    %   saving as a pdf.
    %   
    %   Inputs:
    %       forminList  : array of Formins to gather data from
    %       parameter   : name of parameter containing values to plot against
    %                     polymerization (must be a property in the Formin
    %                     class)
    %       xlab        : x-axis (parameter) label (string)
    %       lab_limit   : upper limit (in parameter) to label points with
    %                     formin name 
    %       scale       : kpoly axis scale (can be none, log2, log10, ln)
    %       settings    : Options class
    %       save        : whether or not to save the plot to the results
    %                     pdf in settings (true/false); deafult is false
    %
    %   Calls figuresave
    %
    %   See also FIGURESAVE, FORMIN, OPTIONS, EXPERIMENT, PRMPLOT, KPOLYMERIZATION.
    arguments
        forminList Formin
        parameter string
        xlab string
        lab_limit double
        scale string {mustBeMember(scale,{'none','log2','log10','ln'})} % scale for kpoly
        settings Options
        save logical=false
    end

    if scale=="none"
        yscale=@(x) x;
        ylab="k_{poly}";
    elseif scale=="log2"
        ylab="log_{2}(k_{poly})";
        yscale=@(x) log2(x);
    elseif scale=="log10"
        ylab="log_{10}(k_{poly})";
        yscale=@(x) log10(x);
    elseif scale=="ln"
        ylab="ln(k_{poly})";
        yscale=@(x) log(x);
    end

    fig=figure;
    xdata=[forminList.(parameter)];
    kpoly=[forminList.kpoly];
    kpoly1_scatter = scatter(xdata,yscale([kpoly.single]), 'filled','o');
    hold on
    kpoly2_scatter = scatter(xdata,yscale([kpoly.double]), 'filled','s');
    hold on
    kpoly3_scatter = scatter(xdata,yscale([kpoly.dimer]), 'filled','p');
    xlabel(xlab)
    ylabel(ylab)
    legend('Single', 'Double', 'N-Dimer', 'Location','northeastoutside');
    
    labels = [forminList.name];
    labelpoints(xdata,yscale([kpoly.single]),labels,'N',0.005, 1,'outliers_lim', {[-inf lab_limit; -inf inf]}, 'FontSize', 7)
    labelpoints(xdata,yscale([kpoly.double]),labels,'E',0.005, 1,'outliers_lim', {[-inf lab_limit; -inf inf]}, 'FontSize', 7)
    labelpoints(xdata,yscale([kpoly.dimer]),labels,'N',0.005, 1,'outliers_lim', {[-inf lab_limit; -inf inf]}, 'FontSize', 7)
    
    title(append('Polymerization Rates vs. ',xlab))
    
    if save
        figuresave(gcf,settings,append(gca().Title.String,'.fig'));
    end
end