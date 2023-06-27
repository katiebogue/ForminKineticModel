function make_bargraph_overviews(data,settings)
% MAKE_BARGRAPH_OVERVIEWS  Makes and saves overview bargraphs for wholepackage.
    %
    %   make_bargraph_overviews(data,settings)
    %
    %   Makes 3 bargraphs based on input structures (data and settings):
    %       1. Single, double, and dimer side-by-side for each formin
    %       2. Change upon NTD for each formin
    %       3. Number of PRMs for each formin (with PRMs)
    %   
    %   Inputs:
    %       data     : structure containing data
    %                  must include: 
    %                  all_kpoly1, all_kpoly2, all_kpoly3, all_fh1_names,
    %                  all_log_kpoly3_2, all_iSite_tot, all_fh1_names_nobind  
    %       settings : structure containing settings variables
    %                  must include: 
    %                  settings_variable, pdf_name
    %
    %   Calls figuresave
    %   
    %   See also FIGURESAVE, MAKE_FORMIN_PLOTS.

    close all 

    % bar graph showing all polymerization rates
    kpoly_table = table(data.all_kpoly1, data.all_kpoly2, data.all_kpoly3, 'RowNames', data.all_fh1_names);
    sorted_kpoly_table = sortrows(kpoly_table);
   
    kpoly_bar = bar(sorted_kpoly_table{:,:});
    set(gca,'xtick',1:length(data.all_fh1_names), 'xticklabel',sorted_kpoly_table.Properties.RowNames);
    xtickangle(90);
    set(kpoly_bar(1), 'FaceColor','b');
    set(kpoly_bar(2), 'FaceColor','r');
    set(kpoly_bar(3), 'FaceColor','g');
    legend( 'single', 'double', 'N-term dimerized');
    xlabel('Formins');
    ylabel('log_2(kpoly)');
    ylim([0 12]);
    
    title('Polymerization Rates per Formin')
    subtitle(settings.settings_variable)
    
    figuresave(gcf,settings.pdf_name,append(gca().Title.String,'.fig'));
    
    close all
    
    % Change in Polymerization Rates w/ Dimerization per Formin
    kpoly_table_ratio = table(data.all_log_kpoly3_2, 'RowNames', data.all_fh1_names);
    sorted_kpoly_table_ratio = sortrows(kpoly_table_ratio);
    
    kpoly_bar_ratio = bar(sorted_kpoly_table_ratio{:,:});
    set(gca,'xtick',1:length(data.all_fh1_names), 'xticklabel',sorted_kpoly_table_ratio.Properties.RowNames)
    xtickangle(90)
    set(kpoly_bar_ratio(1), 'FaceColor','m')
    xlabel('Formins')
    ylabel('log_2(kpoly N terminal dimerized/kpoly double)')
    ylim([-1.7 0.3])
    
    title('Change in Polymerization Rates with NTD per Formin')
    subtitle(settings.settings_variable)
    
    figuresave(gcf,settings.pdf_name,append(gca().Title.String,'.fig'));
    
    close all
    
    % Number of PRMs
    hb = bar(data.all_iSite_tot, 'stacked');
    set(hb(1), 'FaceColor','b');
    hold on
    
    xTick=get(gca,'xtick');
    set(gca,'xtick',1:length(data.all_fh1_names_nobind), 'xticklabel', data.all_fh1_names_nobind)
    xtickangle(90)
    ylim([0 35])
    
    xlabel('Formins')
    ylabel('Binding sites')
    title('Number of Binding Sites per Formin')
    subtitle(settings.settings_variable)
    
    figuresave(gcf,settings.pdf_name,append(gca().Title.String,'.fig'));
    
    close all
end