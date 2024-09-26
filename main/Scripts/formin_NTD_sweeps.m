% FORMIN_NTD_SWEEPS generates Lookuptable, Options, and Experiment for
% formins and saves heatmaps of simulated kpoly values accross multiple NTD
% locations
    % 
    % Make sure to set the lookuptable file, forminfile, and gatingfile 
    %
    % Uses a set of parameters that can be modified, currently the best fit
    % parameters as of 9/21/24
    %
    % Has code commented that could be used to apply best fit params from
    % tables
    %
    % See also FORMIN, EXPERIMENT, LOOKUPTABLE, OPTIONS, EXPDATABAR,
    % KPOLYMERIZATION, PRM, EXPERIMENT/EXPDATABAR.

    %% Input files/ paths
%ltfile="N600_lookup.mat"; % output file from polymer-c; must be on matlab path
ltfile="prvec_runs_lookup.mat"; % output file from polymer-c; must be on matlab path
Nmax=400;
% ltfile="GST_N300.mat"; % output file from polymer-c; must be on matlab path
% Nmax=300;

pythonpath="/Users/katiebogue/MATLAB/GitHub/ForminKineticModel/main/python"; % path to python files
resultsloc="/Users/katiebogue/MATLAB/GitHub/Data/ForminKineticmodel_data/Results"; % path to location to save results
forminfile="ForminTypes.txt"; % file containing sequences
gatingfile="ForminTypes_gating.txt";
%% create lookuptable
lt=(load(ltfile,'lookuptable').lookuptable);
lt=Lookuptable(lt);
%% create options object
% modify this line to change the rate constants:
opts=Options(lt,pythonpath,...
    "3st",...       % kpoly type
    resultsloc,...
    73.0227,...          % k_cap
    0.022909,... % k_del
    37896.5784,...    % r_cap
    1,...    % r_del
    1);     % k_rel

opts.r_cap_exp=0.86103;
opts.set_equation(1); % using preset #1 (see Options class)
opts.NTopt=2; 

%% create experiment object
Experiment1=Experiment(opts,forminfile,"uniprot",0.88); % using the input sequence option and concentration of profilin actin of 0.88
Experiment1.set_gating_file(gatingfile);
%Experiment1.set_gating("BNI1",0.5);

%% Make heatmaps
set(groot,'defaultfigureposition',[400 250 900 750]) % helps prevent cut offs in figs
titles="log_2(k_{poly} N terminal dimerized/k_{poly} double)";

NTDtable=makeNTDtable(Experiment1,20,Nmax);
figure
h1=makeheatmap(NTDtable,Nmax);
h1.Title = {titles,"3st"};
opts.update_results_folder
opts.resultsfolder=strcat(opts.resultsfolder,"NTDsweep_","3st");
figuresave(gcf,opts,append('NTDsweep_','3st','.fig'),true);
figure
h2=makeheatmap_std(NTDtable,Nmax,'fh1length','rollingstd');
h2.Title = {"rolling standard deviation of log2 ratios","3st"};
figuresave(gcf,opts,append('NTDsweep_stds','3st','.fig'),true);
figure
h2=makeheatmap_std(NTDtable,Nmax,'fh1length','std_dimer');
h2.Title = {"rolling standard deviation of dimer kpoly","3st"};
figuresave(gcf,opts,append('NTDsweep_stds_dimer','3st','.fig'),true);
figure
h2=makeheatmap_std(NTDtable,Nmax,'fh1length','std_double');
h2.Title = {"rolling standard deviation of double kpoly","3st"};
figuresave(gcf,opts,append('NTDsweep_stds_double','3st','.fig'),true);
figure
h2=makeheatmap_std(NTDtable,Nmax,'NTD_dists','rollingstd');
h2.Title = {"rolling standard deviation of log2 ratios","3st"};
figuresave(gcf,opts,append('NTDsweep_stds_NT','3st','.fig'),true);
figure
h2=makeheatmap_std(NTDtable,Nmax,'NTD_dists','std_dimer');
h2.Title = {"rolling standard deviation of dimer kpoly","3st"};
figuresave(gcf,opts,append('NTDsweep_stds_dimer_NT','3st','.fig'),true);
figure
h2=makeheatmap_std(NTDtable,Nmax,'NTD_dists','std_double');
h2.Title = {"rolling standard deviation of double kpoly","3st"};
figuresave(gcf,opts,append('NTDsweep_stds_double_NT','3st','.fig'),true);
figure
NTDmaxmin=max_min_table(NTDtable);
b=make_barplot(NTDmaxmin);
title([titles,"3st"]);
figuresave(gcf,opts,append('NTDsweep_bar_','3st','.fig'),true);

%% look for variance cut off
dimdob_dif=NTDtable.doubles-NTDtable.dimers;
dimer_std=movstd(NTDtable.dimers,5,"omitnan");
double_std=movstd(NTDtable.doubles,5,"omitnan");
ratio_std=movstd(NTDtable.ratios_raw,5,"omitnan");
NTDtable.dimer_check=dimer_std.*0.1>abs(dimdob_dif);
NTDtable.double_check=double_std.*0.1>abs(dimdob_dif);
NTDtable.ratio_check=ratio_std.*0.1>abs(dimdob_dif);

NTDtable.dimer_check_diff=dimer_std.*0.1-abs(dimdob_dif);
NTDtable.double_check_diff=double_std.*0.1-abs(dimdob_dif);
NTDtable.ratio_check_diff=ratio_std.*0.1-abs(dimdob_dif);

figure
h = heatmap(NTDtable,'formin_name','fh1length','ColorVariable','dimer_check');
h.ColorMethod = 'none';
h.Colormap=sky(2);
h.NodeChildren(3).YDir='normal';
yvals=str2double(h.YData);
CustomYLabels = string(yvals);
CustomYLabels(mod(yvals,20) ~= 0) = " ";
h.YDisplayLabels = CustomYLabels;
h.GridVisible = 'off';
figuresave(gcf,opts,append('NTDsweep_dimer_check','3st','.fig'),true);

figure
h = heatmap(NTDtable,'formin_name','fh1length','ColorVariable','double_check');
h.ColorMethod = 'none';
h.Colormap=sky(2);
h.NodeChildren(3).YDir='normal';
yvals=str2double(h.YData);
CustomYLabels = string(yvals);
CustomYLabels(mod(yvals,20) ~= 0) = " ";
h.YDisplayLabels = CustomYLabels;
h.GridVisible = 'off';
figuresave(gcf,opts,append('NTDsweep_double_check','3st','.fig'),true);


%% Load best fit tables
% load("best fits.mat")
% fit_3st(:,{'NTopt', 'CTopt'}) = [];
% fit_4st(:,{'NTopt', 'CTopt'}) = [];
% fit_4st_krel(:,{'NTopt', 'CTopt'}) = [];


% Experiment1.applytable(fit_4st)
% NTDtable_4st=makeNTDtable(Experiment1);
% figure
% h1=makeheatmap(NTDtable_4st);
% h1.Title = {titles,"4st"};
% opts.update_results_folder
% opts.resultsfolder=strcat(opts.resultsfolder,"NTDsweep_","4st");
% figuresave(gcf,opts,append('NTDsweep_','4st','.fig'),true);
% NTDmaxmin_4st=max_min_table(NTDtable_4st);
% figure
% b=make_barplot(NTDmaxmin_4st);
% title([titles,"4st"]);
% figuresave(gcf,opts,append('NTDsweep_bar_','4st','.fig'),true);


% Experiment1.applytable(fit_3st)
% NTDtable_3st=makeNTDtable(Experiment1);
% NTDmaxmin_3st=max_min_table(NTDtable_3st);
% figure
% h2=makeheatmap(NTDtable_3st);
% h2.Title = {titles,"3st"};
% opts.update_results_folder
% opts.resultsfolder=strcat(opts.resultsfolder,"NTDsweep_","3st");
% figuresave(gcf,opts,append('NTDsweep_','3st','.fig'),true);
% NTDmaxmin_3st=max_min_table(NTDtable_3st);
% figure
% b=make_barplot(NTDmaxmin_3st);
% title([titles,"3st"]);
% figuresave(gcf,opts,append('NTDsweep_bar_','3st','.fig'),true);

% Experiment1.applytable(fit_4st_krel)
% NTDtable_4st_krel=makeNTDtable(Experiment1);
% NTDmaxmin_4st_krel=max_min_table(NTDtable_4st_krel);
% figure
% h3=makeheatmap(NTDtable_4st_krel);
% h3.Title = {titles,"4st_krel"};
% opts.update_results_folder
% opts.resultsfolder=strcat(opts.resultsfolder,"NTDsweep_","4st_krel");
% figuresave(gcf,opts,append('NTDsweep_','4st_krel','.fig'),true);
% NTDmaxmin_4st_krel=max_min_table(NTDtable_4st_krel);
% figure
% b=make_barplot(NTDmaxmin_4st_krel);
% title([titles,"4st_krel"]);
% figuresave(gcf,opts,append('NTDsweep_bar_','4st_krel','.fig'),true);

function NTDtable=makeNTDtable(exp,rollingavg,Nmax)
% MAKENTDTABLE runs NTD_predictions for each formin and logs resulting
% kpoly values, ratios, and other stats in a table
arguments
    exp %Experiment with all of the formins to make predictions for
    rollingavg %number of nearby values to use in rolling avarge of ratios
    Nmax %max fh1 length value
end
    numformins=length(exp.ForminList);
    NTDtable=table('Size',[numformins*Nmax 10],'VariableTypes',["double","double","double","double","double","double","double","double","double","string"],'VariableNames',{'doubles','dimers','ratios','ratios_raw','NTD_dists','rollingstd','fh1length','std_dimer','std_double','formin_name'});
    for i=1:numformins
        formini=exp.ForminList(i);
        [doubles,dimers,ratios,NTD_dists,fh1ength]=NTD_predictions(formini,Nmax);
        x=(Nmax*(i-1)+1);
        y=Nmax*i;
        NTDtable(x:y,1)=array2table(doubles');
        NTDtable(x:y,2)=array2table(dimers');
        NTDtable(x:y,3)=array2table(movmean(ratios',rollingavg,"omitnan"));
        NTDtable(x:y,4)=array2table(ratios');
        NTDtable(x:y,5)=array2table(NTD_dists');
        NTDtable(x:y,6)=array2table(movstd(ratios',5,"omitnan"));
        NTDtable(x:y,7)=array2table(fh1ength');
        NTDtable(x:y,8)=array2table((movstd(dimers',5,"omitnan"))./mean(dimers,"omitnan"));
        NTDtable(x:y,9)=array2table((movstd(doubles',5,"omitnan"))./mean(doubles,"omitnan"));
        NTDtable(x:y,10)={formini.name};
    end
end
function [doubles,dimers,ratios,NTD_dists,fh1ength]=NTD_predictions(formin,Nmax)
% NTD_predictions sweeps through all possible NTD locations up to fh1 length=Nmax for
% the input formin and returns kpoly values, corresponding NTD distances
% and FH1 lengths
arguments
    formin %formin to make predicitons for
    Nmax %max fh1 length value
end
    doubles=NaN(1,Nmax);
    dimers=NaN(1,Nmax);
    ratios=NaN(1,Nmax);
    NTD_dists=[1:Nmax];
    fh1ength=NaN(1,Nmax);
    for i=0:(Nmax-1)
        if formin.length<Nmax+1
            kpoly=formin.kpoly;
            NTD_dist=formin.PRMList(1,formin.PRMCount).dist_NT;
            doubles(NTD_dist)=kpoly.double;
            dimers(NTD_dist)=kpoly.dimer;
            ratioval=log2(kpoly.ratio);
            fh1ength(NTD_dist)=formin.length;
            if ratioval==-Inf
                ratioval=NaN;
            end
            ratios(NTD_dist)=ratioval;
            formin.add_length(1)
        else
            formin.add_length(-i)
            return
        end
    end
end

function h=makeheatmap(tab,Nmax)
% MAKEHEATMAP creates heatmap of Kpoly ratios for each formin for each
% possible NTD dist, from a table generated by makeNTDtable
arguments 
    tab %table generated by makeNTDtable
    Nmax %max fh1 length value
end
    ogtab=tab;
    % x=movmean(tab.ratios,10,"omitnan");
    % tab.ratios=x;
    h = heatmap(tab,'formin_name','NTD_dists','ColorVariable','ratios');
    h.ColorMethod = 'none';
    h.NodeChildren(3).YDir='normal';
    load('customcolorbar_red_blue_large.mat');
    h.Colormap=CustomColormap;
    max_min_tab=max_min_table(tab);
    allmaxmin=[max_min_tab.max; abs(max_min_tab.min)];
    allmaxmin=sort(allmaxmin);
    allmaxmin = allmaxmin(~isnan(allmaxmin));
    maxratio=allmaxmin(end);
    if maxratio-allmaxmin(end-1)>5
        maxratio=allmaxmin(end-1);
    end
    h.ColorLimits=[-maxratio,maxratio];
    h.ColorLimits=[-5,5];
    yvals=[1:Nmax];
    CustomYLabels = string(yvals);
    CustomYLabels(mod(yvals,20) ~= 0) = " ";
    h.YDisplayLabels = CustomYLabels;
    h.GridVisible = 'off';

end
function h=makeheatmap_std(tab,Nmax,disttype,stdtype)
% MAKEHEATMAP_STD creates heatmap of standard deviations of kpoly/ratios 
% for each formin for each possible NTD dist, from a table generated by makeNTDtable
arguments 
    tab %table generated by makeNTDtable
    Nmax %max fh1 length value
    disttype % y axis for heatmaps, either NTD distances or total FH1 length
    stdtype % color variable for heatmaps, a column in tab that contains std values
end
    h = heatmap(tab,'formin_name',disttype,'ColorVariable',stdtype);
    h.ColorMethod = 'none';
    h.NodeChildren(3).YDir='normal';
    tab.ratios=tab.(stdtype);
    max_min_tab=max_min_table(tab);
    allmaxmin=[max_min_tab.max; abs(max_min_tab.min)];
    allmaxmin=sort(allmaxmin);
    allmaxmin = allmaxmin(~isnan(allmaxmin));
    maxratio=allmaxmin(end);
    if maxratio-allmaxmin(end-1)>1
        maxratio=allmaxmin(end-1);
    end
    h.ColorLimits=[0,maxratio];
    yvals=str2double(h.YData);
    CustomYLabels = string(yvals);
    CustomYLabels(mod(yvals,20) ~= 0) = " ";
    h.YDisplayLabels = CustomYLabels;
    h.GridVisible = 'off';
end
function mmtab=max_min_table(tab)
% MAX_MIN_TABLE creates table of max and min ratios for each formin in the
% input table
    formins=unique(tab.formin_name);
    numformins=length(formins);
    mmtab=table('Size',[numformins 3],'VariableTypes',["double","double","string"],'VariableNames',{'max','min','formin_name'});
    mmtab.formin_name=formins;
    for i=1:numformins
        formin_name=mmtab.formin_name(i);
        ext=get_max_min(tab,formin_name);
        mmtab{i,'max'}=ext(1);
        mmtab{i,'min'}=ext(2);
    end
end
function ext=get_max_min(tab,formin_name)
% GET_MAX_MIN finds max and min ratios for a given formin found in the
% input table
    rows=matches(tab.formin_name,formin_name);
    vals=tab{rows,'ratios'};
    ext=[max(vals), min(vals)];
end 
function b=make_barplot(mmtab)
% MAKE_BARPLOT creates a barplot of the max and min values in the input
% mmtab, which should be generated by max_min_tab
    b=bar(categorical(mmtab.('formin_name')'),[mmtab.('min')';mmtab.('max')'],'stacked');
    xtickangle(90)
    legend({'minimum','maximum'})
    b(1).FaceColor='blue';
    b(2).FaceColor='red';
end