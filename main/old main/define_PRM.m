clc
clear all
pdf_name = 'results11/20.pdf';

% location of look-up tables
path = '../../PolymerData/' ;


%% (1) read output files and extract all values of p_occ

m1 = dlmread(append(path,'single1_300.txt'));
m2 = dlmread(append(path,'double_200.txt'));
m3 = dlmread(append(path,'dimer_122.txt'));

find_pocc

%% (2) Further Initialization

% retrives all names and uniprot queries of formin types into seperate strings
Name_Query = char(importdata('ForminTypes.txt'));
Name_Query = strsplit(Name_Query);

% adds current matlab path to python paths if necessary
% %if count(py.sys.path,'') == 0
%     insert(py.sys.path,int32(0),'');
% end

for LOOP = 1:length(Name_Query)/2
fh1_name = convertCharsToStrings(Name_Query(2*LOOP -1));
query = convertCharsToStrings(Name_Query(2*LOOP));

%     fh1_name = 'fhod3-human'
%     query = 'Q2V2M9 ' % use thise lines (and comment out for loop) to test for 1 specific formin

%% (3) Calls Python and UNIPROT
%M=py.list({'x','y','z'})
%query='Q24120'
%py.help('textwrap')
%calls py script to get polymer stats from UNIPROT
py.importlib.import_module('bioservices')
py.importlib.import_module('PP_interruption')
x = py.PP_interruption.gathering_int(query);

%py.importlib.import_module('gather_info')
%x = py.gather_info.gathering(query);


%changes variable format to matlab doubles
fh1_length = cell2mat(x(1));
x2 = {cell(x(2))}; %pp_index_vec
x2 = x2{1}{1};
% x4 = {cell(x(4))}; %pp_index_vec
% x4= x4{1}{1};
pp_index_vec = [];
pp_index_vec2 = [];
for i = 1:length(x2)
    pp_index_vec = [pp_index_vec x2{i}]; %locations of binding sites along a single filament
end
% for i = 1:length(x4)
%     pp_index_vec2 = [pp_index_vec2 x4{i}];
% end

x3 = {cell(x(3))}; %pp_length_vec
x3 = x3{1}{1};
% x5 = {cell(x(5))}; %pp_length_vec
% x5 = x5{1}{1};
pp_length_vec = [];
pp_length_vec2 = [];
for i = 1:length(x3)
    pp_length_vec = [pp_length_vec x3{i}]; %number of polyprolines at each binding site
end
% for i = 1:length(x5)
%     pp_length_vec2 = [pp_length_vec2 x5{i}];
% end

%sets more variables and constants
iSite_tot = length(pp_index_vec); %number of total binding sites on one filament
iSite_tot2 = length(pp_index_vec2);
%disp(fh1_name);

y=[iSite_tot, iSite_tot2];
hb = bar(LOOP, y, 'stacked');
set(hb(1), 'FaceColor','b')
set(hb(2), 'FaceColor','r')
hold on

% bar(LOOP, iSite_tot,'r')
% hold on
% bar(LOOP, iSite_tot2, 'b')
% hold on
% x=['Diap1--Human', 'Diap2--Human', 'Diap3--Human', 'Diap1--Mouse', 'Diap2--Mouse', 'Diap3--Mouse', 'Diap1--Rat','Diap3--Rat','DAAM1--Human', 'DAAM2--Human', 'DAAM1--Mouse', 'DAAM2--Mouse', 'CAPU--FruitFly', 'FMN1--Human', 'FMN2--Human', 'FMN1--Mouse', 'FMN2--Mouse', 'INF2--Mouse', 'FHOD1--Human', 'FHOD3--Human', 'FHOD1--Mouse', 'FHOD3--Mouse', 'BNR1--Yeast', 'CDC12P--Yeast', 'BNI1P--Yeast']; 
%y =iSite_tot;
end

% for LOOP = 1:length(Name_Query)/2
% bar(LOOP, iSite_tot2, 'b')
% hold on    
% end
legend( '6PI', '6P')
names = {'Diap1--Human', 'Diap2--Human', 'Diap3--Human', 'Diap1--Mouse', 'Diap2--Mouse', 'Diap3--Mouse', 'Diap1--Rat','Diap3--Rat','DAAM1--Human', 'DAAM2--Human', 'DAAM1--Mouse', 'DAAM2--Mouse', 'CAPU--FruitFly', 'FMN1--Human', 'FMN2--Human', 'FMN1--Mouse', 'FMN2--Mouse', 'INF2--Mouse', 'FHOD1--Human', 'FHOD3--Human', 'FHOD1--Mouse', 'FHOD3--Mouse', 'BNR1--Yeast', 'CDC12P--Yeast', 'BNI1P--Yeast'};
xTick=get(gca,'xtick'); 
% xMax=max(xtick); 
% xMin=min(xtick); 
% newXTick=linspace(xMin,xMax,25); 
set(gca,'xtick',[1:25], 'xticklabel', names)
xtickangle(90)
%set(gca, 'XTickLabel', {'Diap1--Human', 'Diap2--Human', 'Diap3--Human', 'Diap1--Mouse', 'Diap2--Mouse', 'Diap3--Mouse', 'Diap1--Rat','Diap3--Rat','DAAM1--Human', 'DAAM2--Human', 'DAAM1--Mouse', 'DAAM2--Mouse', 'CAPU--FruitFly', 'FMN1--Human', 'FMN2--Human', 'FMN1--Mouse', 'FMN2--Mouse', 'INF2--Mouse','})
ylim([0 35])

xlabel('Formins')
ylabel('Binding sites')
savefig('6PI_int_mostrecent.fig')



