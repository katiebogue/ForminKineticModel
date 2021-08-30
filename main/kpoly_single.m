% opt3 = 0
% location of look-up tables
% path = '../../PolymerData/' ;


%% (1) read output files and extract all values of p_occ

m1 = dlmread('single1_300.txt');
m2 = dlmread('double_200.txt');
m3 = dlmread('dimer_122.txt');

find_pocc

%% (2) Further Initialization

% retrives all names and uniprot queries of formin types into seperate strings
Name_Query = char(importdata('ForminTypes.txt'));
Name_Query = strsplit(Name_Query);

% adds current matlab path to python paths if necessary
if count(py.sys.path,'') == 0
     insert(py.sys.path,int32(0),'');
 end

all_fh1_names = [];
all_kpoly1 = [];
all_kpoly2 = [];
all_kpoly3 = [];

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
x = py.PP_interruption.gathering(query);

%changes variable format to matlab doubles
fh1_length = cell2mat(x(1));
x2 = {cell(x(2))}; %pp_index_vec
x2 = x2{1}{1};
% x4 = {cell(x(4))}; %pp_index_vec
% x4= x4{1}{1};
pp_index_vec = [];
% pp_index_vec2 = [];
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
% pp_length_vec2 = [];
for i = 1:length(x3)
    pp_length_vec = [pp_length_vec x3{i}]; %number of polyprorlines at each binding site
end
% for i = 1:length(x5)
%     pp_length_vec2 = [pp_length_vec2 x5{i}];
% end

%sets more variables and constants
iSite_tot = length(pp_index_vec); %number of total binding sites on one filament
% iSite_tot2 = length(pp_index_vec2);

    
        
%% (4) interpolating and graphing
% interpolates values of p_occ using inputed files from (1)
% makes filament schematic
% calculates/graphs polymerization rate

fh1_length_vec = fh1_length*ones(length(pp_index_vec),1);
pp_percentage_vec = transpose(1/fh1_length*pp_index_vec);

p_occ1 = [];
p_occ2a = [];
p_occ2b = [];
p_occ3a = [];
p_occ3b = [];

% extracts specific probabilities for current fh1 from arrays set in find_pocc 
% see (intro) regarding which probabilities are extrapolated
% if not extrapolated, exact simulated probabilities used

if iSite_tot == 0 
    continue
end

% if opt3 == 0
%  if fh1_length > 200 
%      continue
%  end
% end

if fh1_length <= 300
    [row,~] = find(X1(:,1) == fh1_length & X1(:,2) == pp_index_vec);
    p_occ1 = X1(row,3);
    
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
%sets constants
k_paf=10;
c_PA=2.5;


%calculates kpoly for all three states
kp1 = k_paf*c_PA*(1-p_occ1).*pp_length_vec';
k_poly1 = sum(kp1);


kp2a = k_paf*c_PA*(1-p_occ2a).*pp_length_vec';
kp2b = k_paf*c_PA*(1-p_occ2b).*pp_length_vec';
k_poly2 = sum(kp2a) + sum(kp2b);


kp3a = k_paf*c_PA*(1-p_occ3a).*pp_length_vec';
kp3b = k_paf*c_PA*(1-p_occ3b).*pp_length_vec';
k_poly3 = sum(kp3a) + sum(kp3b);


log_kpoly1=log(k_poly1);
log_kpoly2=log(k_poly2);
log_kpoly3= log(k_poly3);


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

pp_index_vec
fh1_length_vec
row
end

kpoly_table = table(all_kpoly1, all_kpoly2, all_kpoly3, 'RowNames', all_fh1_names);
sorted_kpoly_table = sortrows(kpoly_table);

kpoly_bar = bar(sorted_kpoly_table{:,:});
set(gca,'xtick',[1:25], 'xticklabel',sorted_kpoly_table.Properties.RowNames)
xtickangle(90)
set(kpoly_bar(1), 'FaceColor','b')
set(kpoly_bar(2), 'FaceColor','r')
set(kpoly_bar(3), 'FaceColor','g')
legend( 'single', 'double', 'N-term dimerized')
xlabel('Formins')
ylabel('log(kpoly)')
ylim([0 10])

% names = {'Diap1--Human', 'Diap2--Human', 'Diap3--Human', 'Diap1--Mouse', 'Diap2--Mouse', 'Diap3--Mouse', 'Diap1--Rat','Diap3--Rat','DAAM1--Human', 'DAAM2--Human', 'DAAM1--Mouse', 'DAAM2--Mouse', 'CAPU--FruitFly', 'FMN1--Human', 'FMN2--Human', 'FMN1--Mouse', 'FMN2--Mouse', 'INF2--Mouse', 'FHOD1--Human', 'FHOD3--Human', 'FHOD1--Mouse', 'FHOD3--Mouse', 'BNR1--Yeast', 'CDC12P--Yeast', 'BNI1P--Yeast'};
% xTick=get(gca,'xtick'); 
% legend( 'single', 'double', 'N-term dimerized')
% % xMax=max(xtick); 
% % xMin=min(xtick); 
% % newXTick=linspace(xMin,xMax,25); 
% set(gca,'xtick',[1:25], 'xticklabel', names)
% xtickangle(90)
% xlabel('Formins')
% ylabel('kpoly')
% ylim([0 1800])


