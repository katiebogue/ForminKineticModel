% FORMIN_NTD_SWEEPS generates Lookuptable, Options, and Experiment...
    % 
    % 
    %
    % See also FORMIN, EXPERIMENT, LOOKUPTABLE, OPTIONS, EXPDATABAR,
    % KPOLYMERIZATION, PRM, EXPERIMENT/EXPDATABAR.

%% Input files/ paths
ltfile="N600_lookup.mat"; % output file from polymer-c; must be on matlab path
%ltfile="prvec_runs_lookup_N260.mat"; % output file from polymer-c; must be on matlab path

pythonpath="/Users/katiebogue/MATLAB/GitHub/ForminKineticModel/main/python"; % path to python files
resultsloc="/Users/katiebogue/MATLAB/GitHub/Data/ForminKineticmodel_data/Results"; % path to location to save results
forminfile="ForminTypes.txt"; % file containing sequences
%% create lookuptable
lt=(load(ltfile,'lookuptable').lookuptable);
lt=Lookuptable(lt);
%% create options object
% modify this line to change the rate constants:
opts=Options(lt,pythonpath,...
    "4st",...       % kpoly type
    resultsloc,...
    50,...          % k_cap
    4*10^(-1.9),... % k_del
    10^(5.1),...    % r_cap
    2*10^(9),...    % r_del
    8*10^(10));     % k_rel

opts.set_equation(1); % using preset #1 (see Options class)
opts.NTopt=2; 

%% create experiment object
Experiment1=Experiment(opts,forminfile,"uniprot",2.5); % using the input sequence option and concentration of profilin actin of 2.5
%Experiment1.set_gating("BNI1",0.5);

%% Load best fit tables
load("best fits.mat")
fit_3st(:,{'NTopt', 'CTopt'}) = [];
fit_4st(:,{'NTopt', 'CTopt'}) = [];
fit_4st_krel(:,{'NTopt', 'CTopt'}) = [];

%% Make heatmaps
set(groot,'defaultfigureposition',[400 250 900 750]) % helps prevent cut offs in figs
titles="log_2(k_{poly} N terminal dimerized/k_{poly} double)";

Experiment1.applytable(fit_4st)
NTDtable_4st=makeNTDtable(Experiment1);
figure
h1=makeheatmap(NTDtable_4st);
h1.Title = {titles,"4st"};
opts.update_results_folder
opts.resultsfolder=strcat(opts.resultsfolder,"NTDsweep_","4st");
figuresave(gcf,opts,append('NTDsweep_','4st','.fig'),true);


Experiment1.applytable(fit_3st)
NTDtable_3st=makeNTDtable(Experiment1);
figure
h2=makeheatmap(NTDtable_3st);
h2.Title = {titles,"3st"};
opts.update_results_folder
opts.resultsfolder=strcat(opts.resultsfolder,"NTDsweep_","3st");
figuresave(gcf,opts,append('NTDsweep_','3st','.fig'),true);

Experiment1.applytable(fit_4st_krel)
NTDtable_4st_krel=makeNTDtable(Experiment1);
figure
h3=makeheatmap(NTDtable_4st_krel);
h3.Title = {titles,"4st_krel"};
opts.update_results_folder
opts.resultsfolder=strcat(opts.resultsfolder,"NTDsweep_","4st_krel");
figuresave(gcf,opts,append('NTDsweep_','4st_krel','.fig'),true);

function NTDtable=makeNTDtable(exp)
    numformins=length(exp.ForminList);
    NTDtable=table('Size',[numformins*600 5],'VariableTypes',["double","double","double","double","string"],'VariableNames',{'doubles','dimers','ratios','NTD_dists','formin_name'});
    for i=1:numformins
        formini=exp.ForminList(i);
        [doubles,dimers,ratios,NTD_dists]=NTD_predictions(formini);
        x=(600*(i-1)+1);
        y=600*i;
        NTDtable(x:y,1)=array2table(doubles');
        NTDtable(x:y,2)=array2table(dimers');
        NTDtable(x:y,3)=array2table(ratios');
        NTDtable(x:y,4)=array2table(NTD_dists');
        NTDtable(x:y,5)={formini.name};
    end
end
function [doubles,dimers,ratios,NTD_dists]=NTD_predictions(formin)
    doubles=NaN(1,600);
    dimers=NaN(1,600);
    ratios=NaN(1,600);
    NTD_dists=[1:600];
    for i=0:599
        if formin.length<600
            kpoly=formin.kpoly;
            NTD_dist=formin.PRMList(1,formin.PRMCount).dist_NT;
            doubles(NTD_dist)=kpoly.double;
            dimers(NTD_dist)=kpoly.dimer;
            ratios(NTD_dist)=log2(kpoly.ratio);
            formin.add_length(1)
        else
            formin.add_length(-i)
            return
        end
    end
end

function h=makeheatmap(tab)
    h = heatmap(tab,'formin_name','NTD_dists','ColorVariable','ratios');
    h.ColorMethod = 'none';
    h.NodeChildren(3).YDir='normal';
    load('customcolorbar_red_blue_large.mat');
    h.Colormap=CustomColormap;
    maxratio=max(abs(tab.ratios));  
    h.ColorLimits=[-maxratio,maxratio];
    yvals=[1:600];
    CustomYLabels = string(yvals);
    CustomYLabels(mod(yvals,20) ~= 0) = " ";
    h.YDisplayLabels = CustomYLabels;
    h.GridVisible = 'off';

end