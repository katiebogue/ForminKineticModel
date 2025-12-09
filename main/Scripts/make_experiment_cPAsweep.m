% MAKE_EXPERIMENT_CPASWEEP generates Lookuptable, Options, and Experiment
% objects from input file.
    % 
    % Specific for Courtemanche and pollard
    % BNI1P single PRM constructs.
    %
    % See also FORMIN, EXPERIMENT, LOOKUPTABLE, OPTIONS, EXPDATABAR,
    % KPOLYMERIZATION, PRM, EXPERIMENT/EXPDATABAR.

%% Input files/ paths
ltfile="prvec_runs_lookup.mat"; % output file from polymer-c; must be on matlab path

pythonpath="/Users/katiebogue/MATLAB/GitHub/ForminKineticModel/main/python"; % path to python files
resultsloc="/Users/katiebogue/MATLAB/GitHub/Data/ForminKineticmodel_data/Results"; % path to location to save results
forminfile="courtemanche_values_cPA_5.csv"; % file containing sequences, data, and cPA values

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
opts.NTopt=5; % for input sequences

%% create experiment object
Experiment1=Experiment(opts,'','null');

tab=readtable(forminfile);
for i=1:height(tab)
    formin1=Formin(tab.formin_name{i},opts,c_PA=tab.cPA(i),sequence=tab.sequence{i},gating=0.5);
    Experiment1=Experiment1.add_formin(formin1);
    Experiment1=Experiment1.add_data(tab.formin_name{i},tab.kpoly(i),'double',errtop=tab.errtop(i),errbot=tab.errbot(i));
end

%% Make plots
set(groot,'defaultfigureposition',[400 250 1500 800]) % helps prevent cut offs in figs
fig=Experiment1.expdatabar;

