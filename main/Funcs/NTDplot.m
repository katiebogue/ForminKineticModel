function fig=NTDplot(forminList,parameter,xlab,scale,settings,save)
% NTDplot  Creates a scatterplot of change in polymerization vs. a
%          parameter.
    %
    %   fig = NTDPLOT(forminList, x,'xlab','scale',settings) creates a
    %   scatterplot of change in polymerization (scaled as specified) vs. x. 
    %
    %   fig = NTDPLOT(forminList, x,'xlab','scale',settings,true) creates 
    %   and saves a scatterplot of change in polymerization (scaled as 
    %   specified) vs. x. 
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
    %       scale       : y-axis scale (can be none, log2, log10, ln)
    %       settings    : Options class
    %       save        : whether or not to save the plot to the results
    %                     pdf in settings (true/false); deafult is false
    %
    %   Calls figuresave
    %
    %   See also FIGURESAVE, FORMIN, EXPERIMENT, OPTIONS, PRM_NTD_PLOT, KPOLYMERIZATION.

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

    fig=figure;
    xdata=[forminList.(parameter)];
    kpoly=[forminList.kpoly];
    ratio=[kpoly.ratio];
    kpoly_diff_scatter = scatter(xdata,yscale(ratio),'filled');
    xlabel(xlab)
    ylabel(ylab)
    labels = [forminList.name];
    labelpoints(xdata,yscale(ratio),labels,'N',0.005,'outliers_lim',{[-inf inf; -0.6 0]}, 'FontSize', 7)
    
    lne=yline(yscale(1),'LineWidth',2);
    lne.Annotation.LegendInformation.IconDisplayStyle='off';

    title_=append('Change in Polymerization Rates vs. ',xlab);
    title(title_)
    if save
        figuresave(gcf,settings,append(title_,'.fig'));
    end
end