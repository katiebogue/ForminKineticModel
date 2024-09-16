% BNI1_PRED 
    % 
    % 
    %
    % See also FORMIN, EXPERIMENT, LOOKUPTABLE, OPTIONS, EXPDATABAR,
    % KPOLYMERIZATION, PRM, EXPERIMENT/EXPDATABAR.


  ltfile="prvec_runs_lookup.mat"; % output file from polymer-c; must be on matlab path  
  ltfile2="GST_N141.mat"; % output file from polymer-c; must be on matlab path

  pythonpath="/Users/katiebogue/MATLAB/GitHub/ForminKineticModel/main/python"; % path to python files
  resultsloc="/Users/katiebogue/MATLAB/GitHub/Data/ForminKineticmodel_data/Results"; % path to location to save results
  forminfile="BNI1 GST.txt"; % file containing sequences
%% create lookuptable
  lt=(load(ltfile,'lookuptable').lookuptable);
  lt=Lookuptable(lt);
  lt2=(load(ltfile2,'lookuptable').lookuptable);
  lt.addN(lt2.dimer.N141);
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
Experiment1=Experiment(opts,forminfile,"seq",0.88); % using the input sequence option and concentration of profilin actin of 2.5
Experiment1.set_gating("BNI1",0.5);

set(groot,'defaultfigureposition',[400 250 900 750]) % helps prevent cut offs in figs

Experiment1.formingraphic(true);
Experiment1.forminbar(true);