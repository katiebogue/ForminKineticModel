% MAKE_EXPERIMENT_OBJ_EXCEL creates experiment object based on input
% sequence, experimental data, and parameters for formins from an input
% excel file
%
    % Can modify the ltfile to change which lookuptable (polymer-c outputs) to use 
    % Can modify the forminfile to change the input excel file (but be sure
    % to use the format in formintest.xlsx
    %
    % See also FORMIN, EXPERIMENT, LOOKUPTABLE, OPTIONS,
    % KPOLYMERIZATION, PRM.

%% Input files/ paths
ltfile="LookupObject355_largeractin.mat"; % output file from polymer-c; must be on matlab path
forminfile="formintest.xlsx";

pythonpath="/Users/katiebogue/MATLAB/GitHub/ForminKineticModel/main/python"; % path to python files
resultsloc="/Users/katiebogue/MATLAB/GitHub/Data/ForminKineticmodel_data/Results"; % path to location to save results


%% load lookuptable
lt=(load(ltfile,'lt').lt);

%% create options object
% modify this line to change the rate constants:
opts=Options(lt,pythonpath,...
    "3st",...       % kpoly type
    resultsloc,...
    10^3.2778,...          % k_cap
    10^-2.8678,... % k_del
    10^3.15,...    % r_cap
    1,...    % r_del
    1);     % k_rel

opts.r_cap_exp=0.57121;
opts.set_equation(1); % using preset #1 (see Options class)
opts.NTopt=5; % for input sequences

%% create experiment object
Experiment1=Experiment(opts,forminfile,"combo");
