clear
clc
close all
%% (0) INTRO

% saves pdf containing the following info:
% for each formin listed in ForminTypes.txt output:
% schematic displaying filament length and binding site location/strength
% bar graph for calculated polymerization rate for single/double/dimerized
% also shows bar graph for # binding sites for all formins and bar graph
% for all polymerization rates

% calls the following other scripts: 
    % gather_info.py %PP_interruption.py % make_filament_schematic.m
    % find pocc.m % calculate_kpoly.m
  % the following 4 from export_fig:
    % user_string.m % using_hg2.m
    % ghostscript.m % append_pdfs.m

% imports the following look-up tables
    % ForminTypes.txt
    % single_300.txt %double_200.txt %dimer_122.txt

%% (0.1) formatting options:    
%test with less formins:
testing = 'N';

%Determining Binding Sites:
% if interruptions = Y, binding sites are calculated allowing for non proline interruptions
% if interruptions = N, binding sites are calculated without interruptions
interruptions = 'N';

min_PP_length = 4; %defines the minimum length of polyproline region to be considered a binding site

%Defining FH1 domain 
% if standard = Y, all FH1 domains will be defined as starting at first P of a series of at least four Ps with a max of 1 interruption & ending at FH2 domain
% if standard = N, fh1 domains will be defined by uniprot and, if uniprot does not have one defined, the fh1 domain will be defined as above
standard = 'Y';

% if opt = 0, graphs across all formin types have same axes
% if opt = 1, each formin will have individually scaled graphs
% if opt = 2, graphs have same axes unless FH1 length is over 300 amino acids
opt1 = 1;

% if opt2 = 1, graphs are titled with fh1 name and species
% if opt2 = 2, graphs are titled with fh1 name, species, length, and number binding sites
opt2 = 1;

% if opt3 = 0, minimal extrapolation; excludes fh1 with length > 200
% extrapolation always occurs for dimerized > 122 (due to simulation run time)
% if opt3 = 1, all 25 fh1s included; extrapolation occurs in addition for:
    % all filiments length > 300
    % double (and dimer) with length > 200 (due to large number errors)
opt3 = 1;

% if opt4 = 0, saves pdf with each formin on a different page
% if opt4 = 1, creates (but not saves) matlab figures with 3 fh1 per figure
    % if opt3 = 0, grid will have gaps in place of fh1 with length > 200
opt4 = 0;

% location of look-up tables
path = '../../PolymerData/' ;

%location of saved pdf
%filepath = '/Users/Katiebogue/MATLAB/GitHub/ForminKineticModel/main/Results/';

% name of outputed pdf       % must end with '.pdf'
time= datestr(now, 'yyyy-mm-dd HH-MM');
time= convertCharsToStrings(time);
if interruptions == 'Y'
     int_var= 'with_int'
     settings_variable = 'Minimum PP length of ' + string(min_PP_length) + ' ' + 'with interruptions';
end

if interruptions == 'N'
   int_var= 'without_int'
   settings_variable = 'Minimum PP length of ' + string(min_PP_length) + ' ' + 'without interruptions';
end

pdf_name = 'RESULTS_' + time + '_' + 'PPlnth-' + string(min_PP_length) + '_' + int_var + '.pdf';

fig_name = 'RESULTS_' + time + '_' + 'PPlnth-' + string(min_PP_length) + '_' + int_var + '.fig';


%% (1) read output files and extract all values of p_occ

m1 = dlmread('single1_300.txt');        %lookup tables from C code
m2 = dlmread('double_200.txt');
m3 = dlmread('dimer_122.txt');

find_pocc   

%% (2) Further Initialization

% retrives all names and uniprot queries of formin types into seperate strings
if testing == 'Y'
    Name_Query = char(importdata('ForminTypestest.txt')); 
end
if testing == 'N'
    Name_Query = char(importdata('ForminTypes.txt')); 
end 
Name_Query = strsplit(Name_Query);

% adds current matlab path to python paths if necessary
if count(py.sys.path,'') == 0
   insert(py.sys.path,int32(0),'');
end

all_fh1_names = [];
all_fh1_names_nobind = [];
all_kpoly1 = [];
all_kpoly2 = [];
all_kpoly3 = [];

all_log_kpoly3_2 = [];

all_kpoly1_nobind = [];
all_kpoly2_nobind = [];
all_kpoly3_nobind = [];

all_log_kpoly3_2_nobind = [];

all_iSite_tot = [];

all_fh1_length = [];

all_mean_PP_length = [];

all_PP_length = [];

all_PP_location = [];

all_fh1_length_PP = [];

all_kp1 = [];

all_kp2a = [];
all_kp2b = [];

all_kp3a = [];
all_kp3b = [];

all_fh1_names_perPRM = [];

cd '/Users/Katiebogue/MATLAB/GitHub/ForminKineticModel/main/Results';

for LOOP = 1:length(Name_Query)/2
    
fh1_name = convertCharsToStrings(Name_Query(2*LOOP -1));   %takes the names (every other string)
query = convertCharsToStrings(Name_Query(2*LOOP));         %takes the lookup values (every other string)

%     fh1_name = 'fhod3-human'
%     query = 'Q2V2M9 ' % use thise lines (and comment out for loop) to test for 1 specific formin

%% (3) Calls Python and UNIPROT

%calls py script to get polymer stats from UNIPROT
py.importlib.import_module('bioservices')

if interruptions == 'Y'
    py.importlib.import_module('PP_interruption')
    x = py.PP_interruption.gathering_int(query,min_PP_length,standard);   %uses python function defined in PP_interruption that outputs [length of FH1 domain, array of location of PP (middle), array of length of PP] 
end

if interruptions == 'N'
    py.importlib.import_module('gather_info')
    x = py.gather_info.gathering(query,min_PP_length,standard);   %uses python function defined in gather_info that outputs [length of FH1 domain, array of location of PP (middle), array of length of PP] 
end


%changes variable format to matlab doubles
fh1_length = cell2mat(x(1));  %length of FH1 domain from gathering
all_fh1_length = [all_fh1_length; fh1_length];
x2 = {cell(x(2))}; %pp_index_vec  
x2 = x2{1}{1};  %contains locations of PP domains (middle) as a list
pp_index_vec = [];
for i = 1:length(x2)   %for each PP
    pp_index_vec = [pp_index_vec x2{i}]; %locations of binding sites along a single filament
    all_PP_location = [all_PP_location; x2{i}];
    all_fh1_length_PP = [all_fh1_length_PP; fh1_length];
    all_fh1_names_perPRM = [all_fh1_names_perPRM; fh1_name];

end  %makes x2 cells instead of list
x5 = {cell(x(5))}; %pp_index_vec  
x5 = x5{1}{1};  %contains locations of PP domains (middle) as a list
pp_index_vec_start = [];
for i = 1:length(x5)   %for each PP
    pp_index_vec_start = [pp_index_vec_start x5{i}]; %locations of binding sites along a single filament
end  %makes x2 cells instead of list
%does the same thing but for length of PP track:
x3 = {cell(x(3))}; %pp_length_vec
x3 = x3{1}{1};
pp_length_vec = [];
for i = 1:length(x3)
    pp_length_vec = [pp_length_vec x3{i}]; %number of polyprorlines at each binding site
    all_PP_length = [all_PP_length; x3{i}];

end

x4 = {cell(x(4))}; %pp_length_vec
x4 = x4{1}{1};
P_index_vec = [];
for i = 1:length(x4)
    P_index_vec = [P_index_vec x4{i}]; %location of every P

end

mean_PP_length = mean(pp_length_vec);
all_mean_PP_length = [all_mean_PP_length; mean_PP_length];


%sets more variables and constants
iSite_tot = length(pp_index_vec); %number of total binding sites on one filament

all_iSite_tot = [all_iSite_tot; iSite_tot];


% skips fh1 with length over 200 based on chosen opt3
if opt3 == 0
 if fh1_length > 200
     continue
 end
end

%% (4) interpolating and graphing
% interpolates values of p_occ using inputed files from (1)
% makes filament schematic
% calculates/graphs polymerization rate

fh1_length_vec = fh1_length*ones(length(pp_index_vec),1); %makes collumn with each cell = length of FH1
pp_percentage_vec = transpose(1/fh1_length*pp_index_vec); %location of PP over length of FH1


p_occ1 = [];
p_occ2a = [];
p_occ2b = [];
p_occ3a = [];
p_occ3b = [];

% extracts specific probabilities for current fh1 from arrays set in find_pocc 
% see (intro) regarding which probabilities are extrapolated
% if not extrapolated, exact simulated probabilities used

all_fh1_names_nobind = [all_fh1_names_nobind; fh1_name];

if iSite_tot == 0 %skips if no binding sites
    all_kpoly1_nobind = [all_kpoly1_nobind; 0];
    all_kpoly2_nobind = [all_kpoly2_nobind; 0];
    all_kpoly3_nobind = [all_kpoly3_nobind; 0];

    all_log_kpoly3_2_nobind = [all_log_kpoly3_2_nobind; 0];
    
    continue
end

if fh1_length <= 300
    [row,~] = find(X1(:,1) == fh1_length & X1(:,2) == pp_index_vec); %num of row corresponding to fH1 length AND location of PP
    p_occ1 = X1(row,3); %all probabilities for that length at each PP location
    
    if fh1_length <= 200
        [row,~] = find(X2a(:,1) == fh1_length & X2a(:,2) == pp_index_vec);
        p_occ2a = X2a(row,3);
        [row,~] = find(X2b(:,1) == fh1_length & X2b(:,2) == pp_index_vec);
        p_occ2b = X2b(row,3);
    else
        p_occ2a = F2a(fh1_length_vec,pp_index_vec');
        p_occ2b = F2b(fh1_length_vec,pp_index_vec');
    end
    
    if fh1_length <= 121
         [row,~] = find(X3a(:,1) == fh1_length & X3a(:,2) == pp_index_vec);
         p_occ3a = X3a(row,3);
         [row,~] = find(X3b(:,1) == fh1_length & X3b(:,2) == pp_index_vec);
         p_occ3b = X3b(row,3);
    else
        p_occ3a = F3a(fh1_length_vec,pp_index_vec');
        p_occ3b = F3b(fh1_length_vec,pp_index_vec');
    end
else 
    p_occ1 = F1(fh1_length_vec,pp_index_vec');
    p_occ2a = F2a(fh1_length_vec,pp_index_vec');
    p_occ2b = F2b(fh1_length_vec,pp_index_vec'); 
    p_occ3a = F3a(fh1_length_vec,pp_index_vec');
    p_occ3b = F3b(fh1_length_vec,pp_index_vec');
end

% sets title containing name, length, and number of binding sites
% used based on opt2
lenstr = num2str(fh1_length);
iSitestr = num2str(iSite_tot);
fh1_info = strcat(fh1_name,' \\ Length = ',lenstr,' \\ ',iSitestr,' binding sites');

% creates all graphs
make_filament_schematic

calculate_kpoly

all_kp1 = [all_kp1; kp1x];

all_kp2a = [all_kp2a; kp2a];
all_kp2b = [all_kp2b; kp2b];

all_kp3a = [all_kp3a; kp3a];
all_kp3b = [all_kp3b; kp3b];
%% add kpoly combined chart
log_kpoly1= log2(k_poly1);
log_kpoly2= log2(k_poly2);
log_kpoly3= log2(k_poly3);

%Kpoly Ndimer comparison
log_kpoly3_2= log2(k_poly3./k_poly2);


%y = [k_poly1, k_poly2, k_poly3];
% hb = bar(LOOP, y);
% set(hb(1), 'FaceColor','b')
% set(hb(2), 'FaceColor','r')
% set(hb(3), 'FaceColor','g')

hold on
all_fh1_names = [all_fh1_names; fh1_name];
all_kpoly1 = [all_kpoly1; log_kpoly1];
all_kpoly2 = [all_kpoly2; log_kpoly2];
all_kpoly3 = [all_kpoly3; log_kpoly3];

all_log_kpoly3_2 = [all_log_kpoly3_2; log_kpoly3_2];


all_kpoly1_nobind = [all_kpoly1_nobind; log_kpoly1];
all_kpoly2_nobind = [all_kpoly2_nobind; log_kpoly2];
all_kpoly3_nobind = [all_kpoly3_nobind; log_kpoly3];

all_log_kpoly3_2_nobind = [all_log_kpoly3_2_nobind; log_kpoly3_2];

%pp_index_vec
%fh1_length_vec
%row

%% (5) create pdfs

%appends all pdfs into one document called pdf_name

if opt4 == 0
    if LOOP == 1
        saveas(gcf, pdf_name) %gcf= current figure handle
    end

    if LOOP > 1
        saveas(gcf, append('temp.pdf'))
        append_pdfs(pdf_name, append('temp.pdf')) %add to pdf created with the first one
    end
else
    if LOOP == 3
        saveas(gcf, pdf_name) %saves after 3 figures added
    end
    if LOOP > 3 && rem(LOOP,3) == 0  %subsequent 3 figures added
        saveas(gcf, append('temp.pdf'))
        append_pdfs(pdf_name, append('temp.pdf')) %adds the rest to the one pdf
    end
end
    


end
 
if opt4 == 1
    if rem(LOOP,3) == 0 %if the total number of formins is divisible by 3
        
    else
        saveas(gcf, append('temp.pdf'))
        append_pdfs(pdf_name, append('temp.pdf'))
    end
end

all_PP_length_x_PP_isite_tot = all_iSite_tot.*all_mean_PP_length;

all_PP_dist_end = all_fh1_length_PP-all_PP_location;

all_kpoly3a_2a = all_kp3a./all_kp2a;
all_kpoly3b_2b = all_kp3b./all_kp2b;


data_table_all = table(all_kpoly1_nobind, all_kpoly2_nobind, all_kpoly3_nobind, all_log_kpoly3_2_nobind, all_iSite_tot, all_fh1_length, all_mean_PP_length, all_PP_length_x_PP_isite_tot,'RowNames', all_fh1_names_nobind);

%writetable(data_table_all,spreadsheet_name,'Sheet',1);


close all
 kpoly_table = table(all_kpoly1, all_kpoly2, all_kpoly3, 'RowNames', all_fh1_names);
 sorted_kpoly_table = sortrows(kpoly_table);
 
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

%add page with data in numbers to pdf

fig= figure('Name','All Data');
uit = uitable(fig,'Units','Normalized','Position',[0 0 1 1],'ColumnWidth','auto','ColumnName',{'k_poly Single', 'k_poly Double', 'k_poly N-Dimer','log_2(kpoly ratio)','# Binding Sites','FH1 length', 'Mean PRM size','Mean PRM size x # Binding Sites'},'Data',[all_kpoly1_nobind, all_kpoly2_nobind, all_kpoly3_nobind, all_log_kpoly3_2_nobind, all_iSite_tot, all_fh1_length, all_mean_PP_length, all_PP_length_x_PP_isite_tot]);
uit.RowName = all_fh1_names_nobind

saveas(fig,fig_name);

%set(gcf,'PaperPosition',[0 0 8.5 11])

close all

kpoly_diff_PRM_scatter = scatter(all_iSite_tot,all_log_kpoly3_2_nobind, 'filled');
 xlabel('Number of PRMs')
 ylabel('log_2(kpoly N terminal dimerized/kpoly double)')
 labels = all_fh1_names_nobind
 labelpoints(all_iSite_tot,all_log_kpoly3_2_nobind,labels,'N',0.15, 1,'outliers_SD', 1)
%  xpos = all_log_kpoly3_2_nobind(1);
%  ypos = all_iSite_tot(1);
%  labels = all_fh1_names_nobind(1);
%  h = labelpoints (xpos, ypos, labels)

 title('Change in Polymerization Rates vs Number of PRMs')
 subtitle(settings_variable)
 
saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

kpoly1_PRM_scatter = scatter(all_iSite_tot,all_kpoly1_nobind, 'filled');
hold on
kpoly2_PRM_scatter = scatter(all_iSite_tot,all_kpoly2_nobind, 'filled');
hold on
kpoly3_PRM_scatter = scatter(all_iSite_tot,all_kpoly3_nobind, 'filled');
xlabel('Number of PRMs')
ylabel('kpoly')
legend('Single', 'Double', 'N-Dimer', 'Location','northeastoutside');

labels = all_fh1_names_nobind
labelpoints(all_iSite_tot,all_kpoly1_nobind,labels,'N',0.15, 1,'outliers_SD', 1)
labelpoints(all_iSite_tot,all_kpoly2_nobind,labels,'E',0.15, 1,'outliers_SD', 1)
labelpoints(all_iSite_tot,all_kpoly3_nobind,labels,'N',0.15, 1,'outliers_SD', 1)

 title('Polymerization Rates vs Number of PRMs')
 subtitle(settings_variable)

saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

kpoly_diff_length_scatter = scatter(all_fh1_length,all_log_kpoly3_2_nobind, 'filled');
 xlabel('Length of FH1 domain (1st PRM to FH2)')
 ylabel('log_2(kpoly N terminal dimerized/kpoly double)')
 
labels = all_fh1_names_nobind
labelpoints(all_fh1_length,all_log_kpoly3_2_nobind,labels,'N',0.15, 1,'outliers_SD', 1)
 
 title('Change in Polymerization Rates vs Length of FH1 Domain')
 subtitle(settings_variable)

saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

kpoly1_length_scatter = scatter(all_fh1_length,all_kpoly1_nobind, 'filled');
hold on
kpoly2_length_scatter = scatter(all_fh1_length,all_kpoly2_nobind, 'filled');
hold on
kpoly3_length_scatter = scatter(all_fh1_length,all_kpoly3_nobind, 'filled');
xlabel('Length of FH1 domain (1st PRM to FH2)')
ylabel('kpoly')
legend('Single', 'Double', 'N-Dimer', 'Location','northeastoutside');

labels = all_fh1_names_nobind
labelpoints(all_fh1_length,all_kpoly1_nobind,labels,'N',0.15, 1,'outliers_SD', 1)
labelpoints(all_fh1_length,all_kpoly2_nobind,labels,'E',0.15, 1,'outliers_SD', 1)
labelpoints(all_fh1_length,all_kpoly3_nobind,labels,'N',0.15, 1,'outliers_SD', 1)


 title('Polymerization Rates vs Length of FH1 Domain')
 subtitle(settings_variable)

saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

kpoly_diff_size_scatter = scatter(all_mean_PP_length,all_log_kpoly3_2_nobind, 'filled');
 xlabel('Mean PRM size')
 ylabel('log_2(kpoly N terminal dimerized/kpoly double)')

labels = all_fh1_names_nobind
labelpoints(all_mean_PP_length,all_log_kpoly3_2_nobind,labels,'N',0.15, 1,'outliers_SD', 1)

 
 title('Change in Polymerization Rates vs Mean PRM size')
 subtitle(settings_variable)
 
saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

kpoly_diff_size_scatter = scatter(all_PP_length_x_PP_isite_tot,all_log_kpoly3_2_nobind, 'filled');
 xlabel('Mean PRM size x #PRMs')
 ylabel('log_2(kpoly N terminal dimerized/kpoly double)')

labels = all_fh1_names_nobind
labelpoints(all_PP_length_x_PP_isite_tot,all_log_kpoly3_2_nobind,labels,'N',0.15, 1,'outliers_SD', 1)

 
 title('Change in Polymerization Rates vs Mean PRM size x Number of PRMs')
 subtitle(settings_variable)
 
saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all

n=1
kplabels = [];
colors = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F"]
for colorindex = 1:length(all_iSite_tot)/2
    colors = [colors colors];
end
while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY)
        endpoint = n + PRMnum-1
    
        kpoly1_individual_scatter_length = scatter(all_PP_length(n:endpoint),all_kp1(n:endpoint), 'd','filled', 'MarkerFaceColor', colors(LOOPY));
        hold on
        kpoly2a_individual_scatter_length = scatter(all_PP_length(n:endpoint),all_kp2a(n:endpoint),'o', 'filled', 'MarkerFaceColor', colors(LOOPY));
        hold on
        kpoly2b_individual_scatter_length = scatter(all_PP_length(n:endpoint),all_kp2b(n:endpoint),'s', 'filled', 'MarkerFaceColor', colors(LOOPY));
        hold on
        kpoly3a_individual_scatter_length = scatter(all_PP_length(n:endpoint),all_kp3a(n:endpoint),'^', 'filled', 'MarkerFaceColor', colors(LOOPY));
        hold on
        kpoly3b_individual_scatter_length = scatter(all_PP_length(n:endpoint),all_kp3b(n:endpoint),'p', 'filled', 'MarkerFaceColor', colors(LOOPY));
        hold on
        kp1label = all_fh1_names_nobind(LOOPY)+ ' single'
        kp2alabel = all_fh1_names_nobind(LOOPY)+ ' double-1'
        kp2blabel = all_fh1_names_nobind(LOOPY)+ ' double-2'
        kp3alabel = all_fh1_names_nobind(LOOPY)+ ' dimer-1'
        kp3blabel = all_fh1_names_nobind(LOOPY)+ ' dimer-2'
        kplabels = [kplabels kp1label kp2alabel kp2blabel kp3alabel kp3blabel]
        n=n+PRMnum
    end
end

xlabel('Length of PRM')
ylabel('kpoly')

legend(kplabels, 'Location','northeastoutside');

 title('Polymerization Rates vs. PP length per individual PRM')
 subtitle(settings_variable)

saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all
%%%%%%

n=1
kplabels = [];
while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY)
        endpoint = n + PRMnum-1
    
        kpoly1_individual_scatter_loc = scatter(all_PP_location(n:endpoint),all_kp1(n:endpoint), 'd','filled', 'MarkerFaceColor', colors(LOOPY));
        hold on
        kpoly2a_individual_scatter_loc = scatter(all_PP_location(n:endpoint),all_kp2a(n:endpoint),'o', 'filled', 'MarkerFaceColor', colors(LOOPY));
        hold on
        kpoly2b_individual_scatter_loc = scatter(all_PP_location(n:endpoint),all_kp2b(n:endpoint),'s', 'filled', 'MarkerFaceColor', colors(LOOPY));
        hold on
        kpoly3a_individual_scatter_loc = scatter(all_PP_location(n:endpoint),all_kp3a(n:endpoint),'^', 'filled', 'MarkerFaceColor', colors(LOOPY));
        hold on
        kpoly3b_individual_scatter_loc = scatter(all_PP_location(n:endpoint),all_kp3b(n:endpoint),'p', 'filled', 'MarkerFaceColor', colors(LOOPY));
        hold on
        kp1label = all_fh1_names_nobind(LOOPY)+ ' single'
        kp2alabel = all_fh1_names_nobind(LOOPY)+ ' double-1'
        kp2blabel = all_fh1_names_nobind(LOOPY)+ ' double-2'
        kp3alabel = all_fh1_names_nobind(LOOPY)+ ' dimer-1'
        kp3blabel = all_fh1_names_nobind(LOOPY)+ ' dimer-2'
        kplabels = [kplabels kp1label kp2alabel kp2blabel kp3alabel kp3blabel]
        n=n+PRMnum
    end
end

xlabel('Distance from PRM to FH2')
ylabel('kpoly')

legend(kplabels, 'Location','northeastoutside');

 title('Polymerization Rates vs. PP dist to FH2 per individual PRM')
 subtitle(settings_variable)

saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all
%%%%

n=1
kplabels = [];
while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY)
        endpoint = n + PRMnum-1
    
        kpoly1_individual_scatter_dis = scatter(all_PP_dist_end(n:endpoint),all_kp1(n:endpoint), 'd','filled', 'MarkerFaceColor', colors(LOOPY));
        hold on
        kpoly2a_individual_scatter_dis = scatter(all_PP_dist_end(n:endpoint),all_kp2a(n:endpoint),'o', 'filled', 'MarkerFaceColor', colors(LOOPY));
        hold on
        kpoly2b_individual_scatter_dis = scatter(all_PP_dist_end(n:endpoint),all_kp2b(n:endpoint),'s', 'filled', 'MarkerFaceColor', colors(LOOPY));
        hold on
        kpoly3a_individual_scatter_dis = scatter(all_PP_dist_end(n:endpoint),all_kp3a(n:endpoint),'^', 'filled', 'MarkerFaceColor', colors(LOOPY));
        hold on
        kpoly3b_individual_scatter_dis = scatter(all_PP_dist_end(n:endpoint),all_kp3b(n:endpoint),'p', 'filled', 'MarkerFaceColor', colors(LOOPY));
        hold on
        kp1label = all_fh1_names_nobind(LOOPY)+ ' single'
        kp2alabel = all_fh1_names_nobind(LOOPY)+ ' double-1'
        kp2blabel = all_fh1_names_nobind(LOOPY)+ ' double-2'
        kp3alabel = all_fh1_names_nobind(LOOPY)+ ' dimer-1'
        kp3blabel = all_fh1_names_nobind(LOOPY)+ ' dimer-2'
        kplabels = [kplabels kp1label kp2alabel kp2blabel kp3alabel kp3blabel]
        n=n+PRMnum
    end
end

xlabel('Distance from PRM to end')
ylabel('kpoly')

legend(kplabels, 'Location','northeastoutside');

 title('Polymerization Rates vs. PP dist to end per individual PRM')
 subtitle(settings_variable)

saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))


close all
%%%%



n=1
kplabels = [];
while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY)
        endpoint = n + PRMnum-1
    
        kpoly3a_2a_individual_scatter_length = scatter(all_PP_length(n:endpoint),all_kpoly3a_2a(n:endpoint),'d','filled', 'MarkerFaceColor', colors(LOOPY));
        hold on
        kpoly3b_2b_individual_scatter_length = scatter(all_PP_length(n:endpoint),all_kpoly3b_2b(n:endpoint), 'o', 'filled', 'MarkerFaceColor', colors(LOOPY));
        hold on
        kpalabel = all_fh1_names_nobind(LOOPY)+ ' a'
        kpblabel = all_fh1_names_nobind(LOOPY)+ ' b'
        kplabels = [kplabels kpalabel kpblabel]
        n=n+PRMnum
    end
end

xlabel('Length of PRM')
ylabel('kpoly Dimer/ kpoly Double')

legend(kplabels, 'Location','northeastoutside');

 title('Change in Polymerization Rates vs. PP length per individual PRM')
 subtitle(settings_variable)

saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all
%%%%

n=1
kplabels = [];
while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY)
        endpoint = n + PRMnum-1
    
        kpoly3a_2a_individual_scatter_loc = scatter(all_PP_location(n:endpoint),all_kpoly3a_2a(n:endpoint),'d','filled', 'MarkerFaceColor', colors(LOOPY));
        hold on
        kpoly3b_2b_individual_scatter_loc = scatter(all_PP_location(n:endpoint),all_kpoly3b_2b(n:endpoint), 'o', 'filled', 'MarkerFaceColor', colors(LOOPY));
        hold on
        kpalabel = all_fh1_names_nobind(LOOPY)+ ' a'
        kpblabel = all_fh1_names_nobind(LOOPY)+ ' b'
        kplabels = [kplabels kpalabel kpblabel]
        n=n+PRMnum
    end
end

xlabel('Distance from PRM to FH2')
ylabel('kpoly Dimer/ kpoly Double')

legend(kplabels, 'Location','northeastoutside');

 title('Change in Polymerization Rates vs. PP dist to FH2 per individual PRM')
 subtitle(settings_variable)

saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))

close all
%%%%
n=1
kplabels = [];
while n<= length(all_kp1)
    for LOOPY = 1:length(all_iSite_tot)
        PRMnum = all_iSite_tot(LOOPY)
        endpoint = n + PRMnum-1
    
        kpoly3a_2a_individual_scatter_dis = scatter(all_PP_dist_end(n:endpoint),all_kpoly3a_2a(n:endpoint),'d','filled', 'MarkerFaceColor', colors(LOOPY));
        hold on
        kpoly3b_2b_individual_scatter_dis = scatter(all_PP_dist_end(n:endpoint),all_kpoly3b_2b(n:endpoint), 'o', 'filled', 'MarkerFaceColor', colors(LOOPY));
        hold on
        kpalabel = all_fh1_names_nobind(LOOPY)+ ' a'
        kpblabel = all_fh1_names_nobind(LOOPY)+ ' b'
        kplabels = [kplabels kpalabel kpblabel];
        n=n+PRMnum
    end
end

xlabel('Distance from PRM to end')
ylabel('kpoly Dimer/ kpoly Double')

legend(kplabels, 'Location','northeastoutside');

 title('Change in Polymerization Rates vs. PP dist to end per individual PRM')
 subtitle(settings_variable)

saveas(gcf, append('temp.pdf'))
append_pdfs(pdf_name, append('temp.pdf'))


% deletes temporary pdf and closes all matlab figures

close all
delete 'temp.pdf'



