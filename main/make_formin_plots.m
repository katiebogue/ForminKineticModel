%makes all of the correlation plot figures at the end of wholepackage

%% Overview figures
close all
 kpoly_table = table(all_kpoly1, all_kpoly2, all_kpoly3, 'RowNames', all_fh1_names);
 sorted_kpoly_table = sortrows(kpoly_table);
 
%% bar graph showing all polymerization rates
 kpoly_bar = bar(sorted_kpoly_table{:,:});
 set(gca,'xtick',[1:length(all_fh1_names)], 'xticklabel',sorted_kpoly_table.Properties.RowNames);
 xtickangle(90);
 set(kpoly_bar(1), 'FaceColor','b');
 set(kpoly_bar(2), 'FaceColor','r');
 set(kpoly_bar(3), 'FaceColor','g');
 legend( 'single', 'double', 'N-term dimerized');
 xlabel('Formins');
 ylabel('log_2(kpoly)');
 ylim([0 12]);
 
 title('Polymerization Rates per Formin')
 subtitle(settings_variable)
 
 saveas(gcf, append('temp.pdf'))
 append_pdfs(pdf_name, append('temp.pdf'))

close all

%% Change in Polymerization Rates w/ Dimerization per Formin

 kpoly_table_ratio = table(all_log_kpoly3_2, 'RowNames', all_fh1_names);
 sorted_kpoly_table_ratio = sortrows(kpoly_table_ratio);
 
 kpoly_bar_ratio = bar(sorted_kpoly_table_ratio{:,:});
 set(gca,'xtick',[1:length(all_fh1_names)], 'xticklabel',sorted_kpoly_table_ratio.Properties.RowNames)
 xtickangle(90)
 set(kpoly_bar_ratio(1), 'FaceColor','m')
 xlabel('Formins')
 ylabel('log_2(kpoly N terminal dimerized/kpoly double)')
 ylim([-1.8 0.2])
 
 title('Change in Polymerization Rates w/ Dimerization per Formin')
 subtitle(settings_variable)
 
 saveas(gcf, append('temp.pdf'))
 append_pdfs(pdf_name, append('temp.pdf'))

 
close all

%% Number of PRMs

hb = bar(all_iSite_tot, 'stacked');
set(hb(1), 'FaceColor','b');
%set(hb(2), 'FaceColor','r')
hold on

%legend( '6PI', '6P')
%names = {'Diap1--Human', 'Diap2--Human', 'Diap3--Human', 'Diap1--Mouse', 'Diap2--Mouse', 'Diap3--Mouse', 'Diap1--Rat','Diap3--Rat','DAAM1--Human', 'DAAM2--Human', 'DAAM1--Mouse', 'DAAM2--Mouse', 'CAPU--FruitFly', 'FMN1--Human', 'FMN2--Human', 'FMN1--Mouse', 'FMN2--Mouse', 'INF2--Mouse', 'FHOD1--Human', 'FHOD3--Human', 'FHOD1--Mouse', 'FHOD3--Mouse', 'BNR1--Yeast', 'CDC12P--Yeast', 'BNI1P--Yeast'};
xTick=get(gca,'xtick'); 
% xMax=max(xtick); 
% xMin=min(xtick); 
% newXTick=linspace(xMin,xMax,25); 
set(gca,'xtick',[1:length(all_fh1_names_nobind)], 'xticklabel', all_fh1_names_nobind)
xtickangle(90)
%set(gca, 'XTickLabel', {'Diap1--Human', 'Diap2--Human', 'Diap3--Human', 'Diap1--Mouse', 'Diap2--Mouse', 'Diap3--Mouse', 'Diap1--Rat','Diap3--Rat','DAAM1--Human', 'DAAM2--Human', 'DAAM1--Mouse', 'DAAM2--Mouse', 'CAPU--FruitFly', 'FMN1--Human', 'FMN2--Human', 'FMN1--Mouse', 'FMN2--Mouse', 'INF2--Mouse','})
ylim([0 35])

xlabel('Formins')
ylabel('Binding sites')
title('Number of Binding Sites per Formin')
subtitle(settings_variable)

saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

%save fig with data

fig= figure('Name','All Data');
uit = uitable(fig,'Units','Normalized','Position',[0 0 1 1],'ColumnWidth','auto','ColumnName',{'k_poly Single', 'k_poly Double', 'k_poly N-Dimer','log_2(kpoly ratio)','# Binding Sites','FH1 length', 'Mean PRM size','Mean PRM size x # Binding Sites'},'Data',[all_kpoly1_nobind, all_kpoly2_nobind, all_kpoly3_nobind, all_log_kpoly3_2_nobind, all_iSite_tot, all_fh1_length, all_mean_PP_length, all_PP_length_x_PP_isite_tot]);
uit.RowName = all_fh1_names_nobind

saveas(fig,fig_name);

%set(gcf,'PaperPosition',[0 0 8.5 11])

close all

%% Formin correlation plots

%% Change in Polymerization Rates vs Number of PRMs
kpoly_diff_PRM_scatter = scatter(all_iSite_tot,all_log_kpoly3_2_nobind, 'filled');
 xlabel('Number of PRMs')
 ylabel('log_2(kpoly N terminal dimerized/kpoly double)')
 labels = all_fh1_names_nobind
 labelpoints(all_iSite_tot,all_log_kpoly3_2_nobind,labels,'N',0.005,'outliers_lim',{[-inf inf; -0.6 0]}, 'FontSize', 7)
%  xpos = all_log_kpoly3_2_nobind(1);
%  ypos = all_iSite_tot(1);
%  labels = all_fh1_names_nobind(1);
%  h = labelpoints (xpos, ypos, labels)

 title('Change in Polymerization Rates vs Number of PRMs')
 subtitle(settings_variable)
 
saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

%% Polymerization Rates vs Number of PRMs
kpoly1_PRM_scatter = scatter(all_iSite_tot,all_kpoly1_nobind, 'filled','o');
hold on
kpoly2_PRM_scatter = scatter(all_iSite_tot,all_kpoly2_nobind, 'filled','s');
hold on
kpoly3_PRM_scatter = scatter(all_iSite_tot,all_kpoly3_nobind, 'filled','p');
xlabel('Number of PRMs')
ylabel('log_2(kpoly)')
legend('Single', 'Double', 'N-Dimer', 'Location','northeastoutside');

labels = all_fh1_names_nobind
labelpoints(all_iSite_tot,all_kpoly1_nobind,labels,'N',0.005,'outliers_lim', {[-inf 20; 5 9]}, 'FontSize', 7)
labelpoints(all_iSite_tot,all_kpoly2_nobind,labels,'E',0.005,'outliers_lim', {[-inf 20; 5 9]}, 'FontSize', 7)
labelpoints(all_iSite_tot,all_kpoly3_nobind,labels,'N',0.005,'outliers_lim', {[-inf 20; 5 9]}, 'FontSize', 7)

 title('Polymerization Rates vs Number of PRMs')
 subtitle(settings_variable)

saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

%% Change in Polymerization Rates vs Length of FH1 Domain
kpoly_diff_length_scatter = scatter(all_fh1_length,all_log_kpoly3_2_nobind, 'filled');
 xlabel('Length of FH1 domain (1st PRM to FH2)')
 ylabel('log_2(kpoly N terminal dimerized/kpoly double)')
 
labels = all_fh1_names_nobind
labelpoints(all_fh1_length,all_log_kpoly3_2_nobind,labels,'N',0.005,'outliers_lim', {[-inf 300; -0.6 0]}, 'FontSize', 7)
 
 title('Change in Polymerization Rates vs Length of FH1 Domain')
 subtitle(settings_variable)

saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

%% Polymerization Rates vs Length of FH1 Domain
kpoly1_length_scatter = scatter(all_fh1_length,all_kpoly1_nobind, 'filled','o');
hold on
kpoly2_length_scatter = scatter(all_fh1_length,all_kpoly2_nobind, 'filled','s');
hold on
kpoly3_length_scatter = scatter(all_fh1_length,all_kpoly3_nobind, 'filled','p');
xlabel('Length of FH1 domain (1st PRM to FH2)')
ylabel('log_2(kpoly)')
legend('Single', 'Double', 'N-Dimer', 'Location','northeastoutside');

labels = all_fh1_names_nobind
labelpoints(all_fh1_length,all_kpoly1_nobind,labels,'N',0.005, 1,'outliers_lim', {[-inf 400; 5 9]}, 'FontSize', 7)
labelpoints(all_fh1_length,all_kpoly2_nobind,labels,'E',0.005, 1,'outliers_lim', {[-inf 400; 5 9]}, 'FontSize', 7)
labelpoints(all_fh1_length,all_kpoly3_nobind,labels,'N',0.005, 1,'outliers_lim', {[-inf 400; 5 9]}, 'FontSize', 7)


 title('Polymerization Rates vs Length of FH1 Domain')
 subtitle(settings_variable)

saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

%% Change in Polymerization Rates vs Mean PRM size
kpoly_diff_size_scatter = scatter(all_mean_PP_length,all_log_kpoly3_2_nobind, 'filled');
 xlabel('Mean PRM size')
 ylabel('log_2(kpoly N terminal dimerized/kpoly double)')

labels = all_fh1_names_nobind
labelpoints(all_mean_PP_length,all_log_kpoly3_2_nobind,labels,'N',0.005,'outliers_lim', {[-inf inf; -0.6 0]}, 'FontSize', 7)

 
 title('Change in Polymerization Rates vs Mean PRM size')
 subtitle(settings_variable)
 
saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

%% Change in Polymerization Rates vs Mean PRM size x Number of PRMs
kpoly_diff_size_scatter = scatter(all_PP_length_x_PP_isite_tot,all_log_kpoly3_2_nobind, 'filled');
 xlabel('Mean PRM size x #PRMs')
 ylabel('log_2(kpoly N terminal dimerized/kpoly double)')

labels = all_fh1_names_nobind
labelpoints(all_PP_length_x_PP_isite_tot,all_log_kpoly3_2_nobind,labels,'N',0.005,'outliers_lim', {[-inf inf; -0.6 0]}, 'FontSize', 7)

 
 title('Change in Polymerization Rates vs Mean PRM size x Number of PRMs')
 subtitle(settings_variable)
 
saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

%% Per PRM plots
%creates colors
bg = [1 1 1;0.9255 0.9412 0.9451];
rgbmatrix = distinguishable_colors(length(all_fh1_names_nobind));
colors = [];
for k = 1:length(rgbmatrix)
    hexcode = rgb2hex(rgbmatrix(k,:));
    hexcode = convertCharsToStrings(hexcode);
    colors = [colors, hexcode];
end

%create matrix of plot shapes
shapes = ['o';'s';'d';'^';'v';'<';'p';'h'];
for k = 1:length(all_fh1_names_nobind)
    shapes1 = ['o';'s';'d';'^';'v';'<';'p';'h'];
    shapes = [shapes; shapes1];
end 

% convert to log2 scale

all_log_kp1 = log2(all_kp1);

all_log_kp2a = log2(all_kp2a);
all_log_kp2b = log2(all_kp2b);

all_log_kp3a = log2(all_kp3a);
all_log_kp3b = log2(all_kp3b);

%kplabels = [];
% colors = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F"]
% for colorindex = 1:length(all_iSite_tot)/2
%     colors = [colors colors];
% end

%% Vs. PRM length

% All settings

kpoly1_individual_scatter = scatter(all_PP_length,all_log_kp1, 'filled','o');
hold on
kpoly2a_individual_scatter = scatter(all_PP_length,all_log_kp2a, 'filled','s');
hold on
kpoly2b_individual_scatter = scatter(all_PP_length,all_log_kp2b, 'filled','s');
hold on
kpoly3a_individual_scatter = scatter(all_PP_length,all_log_kp3a, 'filled','p');
hold on
kpoly3b_individual_scatter = scatter(all_PP_length,all_log_kp3b, 'filled','p');

xlabel('Length of PP')
ylabel('log_2(kpoly)')
legend('Single', 'Double-1', 'Double-2', 'N-Dimer-1', 'N-Dimer-2','Location','northeastoutside');

 title('Polymerization Rates per individual PRM')
 subtitle(settings_variable)
 
saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

% single
n=1;

while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY);
        endpoint = n + PRMnum-1;
        kpoly1_individual_scatter_length = scatter(all_PP_length(n:endpoint),all_log_kp1(n:endpoint),'filled',shapes(LOOPY),'MarkerFaceColor', colors(LOOPY),'MarkerEdgeColor','k');
        hold on
        n=n+PRMnum;
    end
end

xlabel('Length of PRM')
ylabel('log_2(kpoly)')

legend(all_fh1_names_nobind, 'Location','northeastoutside');

 title('Polymerization Rates vs. PP length per individual PRM (single)')
 subtitle(settings_variable) 
 
saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

% double
n=1;

while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY);
        endpoint = n + PRMnum-1;
        kpoly2a_individual_scatter_length = scatter(all_PP_length(n:endpoint),all_log_kp2a(n:endpoint), 'filled',shapes(LOOPY), 'MarkerFaceColor', colors(LOOPY),'MarkerEdgeColor','k');
        hold on
        n=n+PRMnum;
    end
end

xlabel('Length of PRM')
ylabel('log_2(kpoly)')

legend(all_fh1_names_nobind, 'Location','northeastoutside');

n=1;
while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY);
        endpoint = n + PRMnum-1;
        kpoly2b_individual_scatter_length = scatter(all_PP_length(n:endpoint),all_log_kp2b(n:endpoint), 'filled',shapes(LOOPY), 'MarkerFaceColor', colors(LOOPY),'MarkerEdgeColor','k');
        hold on
        n=n+PRMnum;
    end
end

 title('Polymerization Rates vs. PP length per individual PRM (double)')
 subtitle(settings_variable) 
 
saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

% dimer
n=1;

while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY);
        endpoint = n + PRMnum-1;
        kpoly3a_individual_scatter_length = scatter(all_PP_length(n:endpoint),all_log_kp3a(n:endpoint), 'filled',shapes(LOOPY), 'MarkerFaceColor', colors(LOOPY),'MarkerEdgeColor','k');
        hold on
        n=n+PRMnum;
    end
end

xlabel('Length of PRM')
ylabel('log_2(kpoly)')

legend(all_fh1_names_nobind, 'Location','northeastoutside');

n=1;
while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY);
        endpoint = n + PRMnum-1;
        kpoly3b_individual_scatter_length = scatter(all_PP_length(n:endpoint),all_log_kp3b(n:endpoint), 'filled',shapes(LOOPY), 'MarkerFaceColor', colors(LOOPY),'MarkerEdgeColor','k');
        hold on
        n=n+PRMnum;
    end
end

 title('Polymerization Rates vs. PP length per individual PRM (dimer)')
 subtitle(settings_variable) 
 
saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all
%% Distance from PRM to FH2

% All settings

kpoly1_individual_scatter_loc = scatter(all_PP_location,all_log_kp1, 'filled','o');
hold on
kpoly2a_individual_scatter_loc = scatter(all_PP_location,all_log_kp2a, 'filled','s');
hold on
kpoly2b_individual_scatter_loc = scatter(all_PP_location,all_log_kp2b, 'filled','s');
hold on
kpoly3a_individual_scatter_loc = scatter(all_PP_location,all_log_kp3a, 'filled','p');
hold on
kpoly3b_individual_scatter_loc = scatter(all_PP_location,all_log_kp3b, 'filled','p');

xlabel('Distance from PRM to FH2')
ylabel('log_2(kpoly)')
legend('Single', 'Double-1', 'Double-2', 'N-Dimer-1', 'N-Dimer-2','Location','northeastoutside');

 title('Polymerization Rates vs. PP dist to FH2 per individual PRM')
 subtitle(settings_variable)
 
saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

% single
n=1;

while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY);
        endpoint = n + PRMnum-1;
        kpoly1_individual_scatter_length_loc = scatter(all_PP_location(n:endpoint),all_log_kp1(n:endpoint),'filled',shapes(LOOPY),'MarkerFaceColor', colors(LOOPY),'MarkerEdgeColor','k');
        hold on
        n=n+PRMnum;
    end
end

xlabel('Distance from PRM to FH2')
ylabel('log_2(kpoly)')

legend(all_fh1_names_nobind, 'Location','northeastoutside');

 title('Polymerization Rates vs. PP dist to FH2 per individual PRM (single)')
 subtitle(settings_variable) 
 
saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

% double
n=1;

while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY);
        endpoint = n + PRMnum-1;
        kpoly2a_individual_scatter_length_loc = scatter(all_PP_location(n:endpoint),all_log_kp2a(n:endpoint), 'filled',shapes(LOOPY), 'MarkerFaceColor', colors(LOOPY),'MarkerEdgeColor','k');
        hold on
        n=n+PRMnum;
    end
end

xlabel('Distance from PRM to FH2')
ylabel('log_2(kpoly)')

legend(all_fh1_names_nobind, 'Location','northeastoutside');

n=1;
while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY);
        endpoint = n + PRMnum-1;
        kpoly2b_individual_scatter_length_loc = scatter(all_PP_location(n:endpoint),all_log_kp2b(n:endpoint), 'filled',shapes(LOOPY), 'MarkerFaceColor', colors(LOOPY),'MarkerEdgeColor','k');
        hold on
        n=n+PRMnum;
    end
end

 title('Polymerization Rates vs. PP dist to FH2 per individual PRM (double)')
 subtitle(settings_variable) 
 
saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

% dimer
n=1;

while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY);
        endpoint = n + PRMnum-1;
        kpoly3a_individual_scatter_length_loc = scatter(all_PP_location(n:endpoint),all_log_kp3a(n:endpoint), 'filled',shapes(LOOPY), 'MarkerFaceColor', colors(LOOPY),'MarkerEdgeColor','k');
        hold on
        n=n+PRMnum;
    end
end

xlabel('Distance from PRM to FH2')
ylabel('log_2(kpoly)')

legend(all_fh1_names_nobind, 'Location','northeastoutside');

n=1;
while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY);
        endpoint = n + PRMnum-1;
        kpoly3b_individual_scatter_length_loc = scatter(all_PP_location(n:endpoint),all_log_kp3b(n:endpoint), 'filled',shapes(LOOPY), 'MarkerFaceColor', colors(LOOPY),'MarkerEdgeColor','k');
        hold on
        n=n+PRMnum;
    end
end

 title('Polymerization Rates vs. PP dist to FH2 per individual PRM (dimer)')
 subtitle(settings_variable) 
 
saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

%% Distance from PRM to end

% All settings

kpoly1_individual_scatter_dis= scatter(all_PP_dist_end,all_log_kp1, 'filled','o');
hold on
kpoly2a_individual_scatter_dis= scatter(all_PP_dist_end,all_log_kp2a, 'filled','s');
hold on
kpoly2b_individual_scatter_dis= scatter(all_PP_dist_end,all_log_kp2b, 'filled','s');
hold on
kpoly3a_individual_scatter_dis= scatter(all_PP_dist_end,all_log_kp3a, 'filled','p');
hold on
kpoly3b_individual_scatter_dis= scatter(all_PP_dist_end,all_log_kp3b, 'filled','p');

xlabel('Distance from PRM to N-term')
ylabel('log_2(kpoly)')
legend('Single', 'Double-1', 'Double-2', 'N-Dimer-1', 'N-Dimer-2','Location','northeastoutside');

 title('Polymerization Rates vs. PP dist to N-term per individual PRM')
 subtitle(settings_variable)
 
saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

% single
n=1;

while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY);
        endpoint = n + PRMnum-1;
        kpoly1_individual_scatter_length_dis= scatter(all_PP_dist_end(n:endpoint),all_log_kp1(n:endpoint),'filled',shapes(LOOPY),'MarkerFaceColor', colors(LOOPY),'MarkerEdgeColor','k');
        hold on
        n=n+PRMnum;
    end
end

xlabel('Distance from PRM to N-term')
ylabel('log_2(kpoly)')

legend(all_fh1_names_nobind, 'Location','northeastoutside');

 title('Polymerization Rates vs. PP dist to N-term per individual PRM (single)')
 subtitle(settings_variable) 
 
saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

% double
n=1;

while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY);
        endpoint = n + PRMnum-1;
        kpoly2a_individual_scatter_length_dis= scatter(all_PP_dist_end(n:endpoint),all_log_kp2a(n:endpoint), 'filled',shapes(LOOPY), 'MarkerFaceColor', colors(LOOPY),'MarkerEdgeColor','k');
        hold on
        n=n+PRMnum;
    end
end

xlabel('Distance from PRM to N-term')
ylabel('log_2(kpoly)')

legend(all_fh1_names_nobind, 'Location','northeastoutside');

n=1;
while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY);
        endpoint = n + PRMnum-1;
        kpoly2b_individual_scatter_length_dis= scatter(all_PP_dist_end(n:endpoint),all_log_kp2b(n:endpoint), 'filled',shapes(LOOPY), 'MarkerFaceColor', colors(LOOPY),'MarkerEdgeColor','k');
        hold on
        n=n+PRMnum;
    end
end

 title('Polymerization Rates vs. PP dist to N-term per individual PRM (double)')
 subtitle(settings_variable) 
 
saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

% dimer
n=1;

while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY);
        endpoint = n + PRMnum-1;
        kpoly3a_individual_scatter_length_dis= scatter(all_PP_dist_end(n:endpoint),all_log_kp3a(n:endpoint), 'filled',shapes(LOOPY), 'MarkerFaceColor', colors(LOOPY),'MarkerEdgeColor','k');
        hold on
        n=n+PRMnum;
    end
end

xlabel('Distance from PRM to N-term')
ylabel('log_2(kpoly)')

legend(all_fh1_names_nobind, 'Location','northeastoutside');

n=1;
while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY);
        endpoint = n + PRMnum-1;
        kpoly3b_individual_scatter_length_dis= scatter(all_PP_dist_end(n:endpoint),all_log_kp3b(n:endpoint), 'filled',shapes(LOOPY), 'MarkerFaceColor', colors(LOOPY),'MarkerEdgeColor','k');
        hold on
        n=n+PRMnum;
    end
end

 title('Polymerization Rates vs. PP dist to N-term per individual PRM (dimer)')
 subtitle(settings_variable) 
 
saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

%% Change in Polymerization Rates vs. PP length per individual PRM

% two filaments

kpoly3a_2a_individual_scatter_length = scatter(all_PP_length,all_kpoly3a_2a,'filled', 'o');
hold on
kpoly3b_2b_individual_scatter_length = scatter(all_PP_length,all_kpoly3b_2b, 'filled', 's');
hold on

xlabel('Length of PRM')
ylabel('log_2(kpoly Dimer/ kpoly Double)')
legend('filament a', 'filament b');

 title('Change in Polymerization Rates vs. PP length per individual PRM')
 subtitle(settings_variable)
 
saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

% formin legend
n=1;

while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY);
        endpoint = n + PRMnum-1;
        kpoly3a_2a_individual_scatter_length = scatter(all_PP_length(n:endpoint),all_kpoly3a_2a(n:endpoint), 'filled',shapes(LOOPY), 'MarkerFaceColor', colors(LOOPY),'MarkerEdgeColor','k');
        hold on
        n=n+PRMnum;
    end
end

xlabel('Length of PRM')
ylabel('log_2(kpoly Dimer/ kpoly Double)')

legend(all_fh1_names_nobind, 'Location','northeastoutside');

n=1;
while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY);
        endpoint = n + PRMnum-1;
        kpoly3b_2b_individual_scatter_length = scatter(all_PP_length(n:endpoint),all_kpoly3b_2b(n:endpoint), 'filled',shapes(LOOPY), 'MarkerFaceColor', colors(LOOPY),'MarkerEdgeColor','k');
        hold on
        n=n+PRMnum;
    end
end

 title('Change in Polymerization Rates vs. PP length per individual PRM')
 subtitle(settings_variable) 
 
saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

%% Change in Polymerization Rates vs. PP dist to FH2 per individual PRM

% two filaments

kpoly3a_2a_individual_scatter_loc = scatter(all_PP_location,all_kpoly3a_2a,'filled', 'o');
hold on
kpoly3b_2b_individual_scatter_loc = scatter(all_PP_location,all_kpoly3b_2b, 'filled', 's');
hold on

xlabel('Distance from PRM to FH2')
ylabel('log_2(kpoly Dimer/ kpoly Double)')
legend('filament a', 'filament b');

 title('Change in Polymerization Rates vs. PP dist to FH2 per individual PRM')
 subtitle(settings_variable)
 
saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

% formin legend
n=1;

while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY);
        endpoint = n + PRMnum-1;
        kpoly3a_2a_individual_scatter_loc = scatter(all_PP_location(n:endpoint),all_kpoly3a_2a(n:endpoint), 'filled',shapes(LOOPY), 'MarkerFaceColor', colors(LOOPY),'MarkerEdgeColor','k');
        hold on
        n=n+PRMnum;
    end
end

xlabel('Distance from PRM to FH2')
ylabel('log_2(kpoly Dimer/ kpoly Double)')

legend(all_fh1_names_nobind, 'Location','northeastoutside');

n=1;
while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY);
        endpoint = n + PRMnum-1;
        kpoly3b_2b_individual_scatter_loc = scatter(all_PP_location(n:endpoint),all_kpoly3b_2b(n:endpoint), 'filled',shapes(LOOPY), 'MarkerFaceColor', colors(LOOPY),'MarkerEdgeColor','k');
        hold on
        n=n+PRMnum;
    end
end

 title('Change in Polymerization Rates vs. PP dist to FH2 per individual PRM')
 subtitle(settings_variable) 
 
saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

%% Change in Polymerization Rates vs. PP dist to FH2 per individual PRM

% two filaments

kpoly3a_2a_individual_scatter_dis = scatter(all_PP_dist_end,all_kpoly3a_2a,'filled', 'o');
hold on
kpoly3b_2b_individual_scatter_dis = scatter(all_PP_dist_end,all_kpoly3b_2b, 'filled', 's');
hold on

xlabel('Distance from PRM to N-term')
ylabel('log_2(kpoly Dimer/ kpoly Double)')
legend('filament a', 'filament b');

 title('Change in Polymerization Rates vs. PP dist to N-term per individual PRM')
 subtitle(settings_variable)
 
saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

% formin legend
n=1;

while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY);
        endpoint = n + PRMnum-1;
        kpoly3a_2a_individual_scatter_dis = scatter(all_PP_dist_end(n:endpoint),all_kpoly3a_2a(n:endpoint), 'filled',shapes(LOOPY), 'MarkerFaceColor', colors(LOOPY),'MarkerEdgeColor','k');
        hold on
        n=n+PRMnum;
    end
end

xlabel('Distance from PRM to N-term')
ylabel('log_2(kpoly Dimer/ kpoly Double)')

legend(all_fh1_names_nobind, 'Location','northeastoutside');

n=1;
while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY);
        endpoint = n + PRMnum-1;
        kpoly3b_2b_individual_scatter_dis = scatter(all_PP_dist_end(n:endpoint),all_kpoly3b_2b(n:endpoint), 'filled',shapes(LOOPY), 'MarkerFaceColor', colors(LOOPY),'MarkerEdgeColor','k');
        hold on
        n=n+PRMnum;
    end
end

 title('Change in Polymerization Rates vs. PP dist to N-term per individual PRM')
 subtitle(settings_variable) 
 
saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

%% Add notes page
fig_notes= figure('Name','Notes');
uit = uitable(fig_notes,'Units','Normalized','Position',[0 0 1 1],'ColumnWidth','auto','Data',{min_PP_length;interruptions;delivery;del_state;ownseq;ownseqfile;c_PA;k_paf;k_pab;k_paf_rev;r_PF_rev;r_paf_rev;notes});
uit.RowName = {'Min PRM length','Interruptions','Delivery?','States','Using sequence file?','sequence file','Profilin-actin concentration (μM)','capture const. k_paf (μM^(-1)s^(-1))','delivery const, k_pab (μM^(-1)s^(-1))','reverse capture rate, k_paf_rev (s^(-1))','ring opening rate, r_PF_rev (s^(-1))', 'reverse delivery rate, r_paf_rev (s^(-1))','Notes'}

saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all


