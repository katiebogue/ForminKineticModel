function NTDplot(parameter,xlab,title_,Data,settings)
% NTDplot  Creates and saves a scatterplot of change in polymerization vs.
%          a parameter.
    %
    %   NTDPLOT(x,'xlab','title',Data,settings) creates and saves a
    %   scatterplot of change in polymerization vs. x. 
    %   
    %   Inputs:
    %       parameter   : array containing values to plot against change in
    %                     polymerization
    %       xlab        : x-axis (parameter) label (string)
    %       title_      : title of scatterplot (string)
    %       Data        : structure containing data
    %                     must include: 
    %                     all_log_kpoly3_2_nobind, all_fh1_names_nobind 
    %       settings    : structure containing settings variables
    %                     must include: 
    %                     settings_variable, pdf_name
    %
    %   Calls figuresave
    %
    %   See also FIGURESAVE, MAKE_FORMIN_PLOTS, WHOLEPACKAGE.

    kpoly_diff_scatter = scatter(parameter,Data.all_log_kpoly3_2_nobind, 'filled');
    xlabel(xlab)
    ylabel('log_2(kpoly N terminal dimerized/kpoly double)')
    labels = Data.all_fh1_names_nobind;
    labelpoints(parameter,Data.all_log_kpoly3_2_nobind,labels,'N',0.005,'outliers_lim',{[-inf inf; -0.6 0]}, 'FontSize', 7)
    
    title(title_)
    subtitle(settings.settings_variable)
     
    figuresave(gcf,settings.pdf_name,append(title_,'.fig'));
    close all
end