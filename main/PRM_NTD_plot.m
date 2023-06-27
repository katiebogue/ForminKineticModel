function PRM_NTD_plot(parameter,xlab,d,settings)
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
    %       parameter   : array containing values to plot against change
    %                     in polymerization
    %       xlab        : x-axis (parameter) label (string)
    %       d           : structure containing data
    %                     must include: 
    %                     all_fh1_names_nobind, all_kpoly3a_2a,
    %                     all_kpoly3b_2b, all_iSite_tot
    %       settings    : structure containing settings variables
    %                     must include: 
    %                     settings_variable, pdf_name, shapes, colors
    %   
    %   Calls figuresave but only saves one .fig file for both figures
    %
    %   See also FIGURESAVE, MAKE_FORMIN_PLOTS, WHOLEPACKAGE, MAKEPOINTS.
    
title_=append('Change in Polymerization Rates vs. ',xlab,' per individual PRM');
%% Plot based on each filament
h(1)=figure;
kpoly3a_2a_individual_scatter_length = scatter(parameter,d.all_kpoly3a_2a,'filled', 'o');
hold on
kpoly3b_2b_individual_scatter_length = scatter(parameter,d.all_kpoly3b_2b, 'filled', 's');
hold on

xlabel(xlab)
ylabel('log_2(kpoly Dimer/ kpoly Double)')
legend('filament a', 'filament b');

title(title_)
subtitle(settings.settings_variable)
 
figuresave(gcf,settings.pdf_name,append(gca().Title.String,'.fig'));


%% Plot based on formin
h(2)=figure;
n=1;

while n<= length(d.all_kpoly3a_2a)
    for LOOPY = 1:length(d.all_iSite_tot)
        PRMnum = d.all_iSite_tot(LOOPY);
        endpoint = n + PRMnum-1;
        x=[parameter(n:endpoint);parameter(n:endpoint)];
        y=[d.all_kpoly3a_2a(n:endpoint);d.all_kpoly3b_2b(n:endpoint)];
        kpoly3a_2a_individual_scatter_loc = scatter(x,y, 'filled',settings.shapes(LOOPY), 'MarkerFaceColor', settings.colors(LOOPY),'MarkerEdgeColor','k');
        hold on
        n=n+PRMnum;
    end
end

xlabel(xlab)
ylabel('log_2(kpoly Dimer/ kpoly Double)')

legend(d.all_fh1_names_nobind, 'Location','northeastoutside');

title(title_)
subtitle(settings.settings_variable) 
 
figuresave(gcf,settings.pdf_name,append(gca().Title.String,'.fig'));
savefig(h,fullfile("figures",append(gca().Title.String,'.fig')));

close all

end