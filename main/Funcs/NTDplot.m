function fig=NTDplot(forminList,parameter,xlab,title_,settings,save)
% NTDplot  Creates and saves a scatterplot of change in polymerization vs.
%          a parameter.
    %
    %   NTDPLOT(forminList, x,'xlab','title',settings) creates and saves a
    %   scatterplot of change in polymerization vs. x. 
    %   
    %   Inputs:
    %       forminList  : array of Formins to gather data from
    %       parameter   : name of parameter containing values to plot against
    %                     polymerization
    %       xlab        : x-axis (parameter) label (string)
    %       title_      : title of scatterplot (string)
    %                     all_log_kpoly3_2_nobind, all_fh1_names_nobind 
    %       settings    : Options class
    %       save        : whether or not to save the plot to the results
    %                     pdf in settings (true/false); deafult is false
    %
    %   Calls figuresave
    %
    %   See also FIGURESAVE, FORMIN, EXPERIMENT, OPTIONS.

    arguments
        forminList Formin
        parameter string
        xlab string
        title_ string
        settings Options
        save logical=false
    end
    
    fig=figure;
    xdata=[forminList.(parameter)];
    kpoly=[forminList.kpoly];
    ratio=[kpoly.ratio];
    kpoly_diff_scatter = scatter(xdata,ratio,'filled');
    xlabel(xlab)
    ylabel('(kpoly N terminal dimerized/kpoly double)')
    labels = [forminList.name];
    labelpoints(xdata,ratio,labels,'N',0.005,'outliers_lim',{[-inf inf; -0.6 0]}, 'FontSize', 7)
    
    lne=yline(1,'LineWidth',2);
    lne.Annotation.LegendInformation.IconDisplayStyle='off';

    title(title_)
    if save
        figuresave(gcf,settings,append(title_,'.fig'));
    end
end