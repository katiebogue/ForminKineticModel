function kpolyplot(parameter,xlab,lab_limit,d,settings)
% KPOLYPLOT  Creates and saves a scatterplot of polymerization vs. a
%            parameter. 
    %
    %   KPOLYPLOT(x,'xlab',N,d,settings) creates and saves a
    %   scatterplot of polymerization vs. x, labeling points beyond N
    %   
    %   Inputs:
    %       parameter   : array containing values to plot against
    %                     polymerization
    %       xlab        : x-axis (parameter) label (string)
    %       lab_limit   : upper limit (in parameter) to label points with
    %                     formin name
    %       d           : structure containing data
    %                     must include: 
    %                     all_kpoly1_nobind, all_kpoly2_nobind,
    %                     all_kpoly3_nobind, all_fh1_names_nobind  
    %       settings    : structure containing settings variables
    %                     must include: 
    %                     settings_variable, pdf_name
    %
    %   Calls figuresave
    %
    %   See also FIGURESAVE, MAKE_FORMIN_PLOTS, WHOLEPACKAGE.

kpoly1_scatter = scatter(parameter,d.all_kpoly1_nobind, 'filled','o');
hold on
kpoly2_scatter = scatter(parameter,d.all_kpoly2_nobind, 'filled','s');
hold on
kpoly3_scatter = scatter(parameter,d.all_kpoly3_nobind, 'filled','p');
xlabel('Length of FH1 domain (1st PRM to FH2)')
ylabel('log_2(kpoly)')
legend('Single', 'Double', 'N-Dimer', 'Location','northeastoutside');

labels = d.all_fh1_names_nobind;
labelpoints(parameter,d.all_kpoly1_nobind,labels,'N',0.005, 1,'outliers_lim', {[-inf lab_limit; 5 9]}, 'FontSize', 7)
labelpoints(parameter,d.all_kpoly2_nobind,labels,'E',0.005, 1,'outliers_lim', {[-inf lab_limit; 5 9]}, 'FontSize', 7)
labelpoints(parameter,d.all_kpoly3_nobind,labels,'N',0.005, 1,'outliers_lim', {[-inf lab_limit; 5 9]}, 'FontSize', 7)

title(append('Polymerization Rates vs. ',xlab))
subtitle(settings.settings_variable)
 
figuresave(gcf,settings.pdf_name,append(gca().Title.String,'.fig'));
close all
end