% MAKE_EXPERIMENT_OBJS generates Lookuptable, Options, and Experiment
% objects from input files.
    % 
    % Specific for FHODB and Capu NTD data and Courtemanche and pollard
    % BNI1P single PRM constructs.
    %
    % Generates plots of simulated vs. experimental values using the
    % expdatabar function in the Experiment class.
    %
    % See also FORMIN, EXPERIMENT, LOOKUPTABLE, OPTIONS, EXPDATABAR,
    % KPOLYMERIZATION, PRM, EXPERIMENT/EXPDATABAR.

%% Input files/ paths
ltfile="N600_lookup.mat"; % output file from polymer-c; must be on matlab path
ltfile="prvec_runs_lookup.mat"; % output file from polymer-c; must be on matlab path

pythonpath="/Users/katiebogue/MATLAB/GitHub/ForminKineticModel/main/python"; % path to python files
resultsloc="/Users/katiebogue/MATLAB/GitHub/Data/ForminKineticmodel_data/Results"; % path to location to save results
forminfile="Experimental constructs.txt"; % file containing sequences
%forminfile="Quinlan_FHODB variations.txt"; % file containing sequences

%ltfile2="basesep_35.5_output.mat"; % need to use base sep 35.5 runs for capu and fhod dimer
ltfile2="dimerdisttest.mat"; % test with larger dimerdist runs for capu and fhod dimer

% for FhodB and Capu, should we use SD or SEM for the error bars?
ertype='SEM';
%ertype='SD';
%% create lookuptable
lt=(load(ltfile,'lookuptable').lookuptable);
lt=Lookuptable(lt);
 lt2=(load(ltfile2,'lookuptable').lookuptable);
 lt.addN(lt2.dimer.N123);
% lt.addN(lt2.dimer.N255);

%% create options object
% modify this line to change the rate constants:
opts=Options(lt,pythonpath,...
    "3st",...       % kpoly type
    resultsloc,...
    50,...          % k_cap
    4*10^(-1.9),... % k_del
    10^(5.1),...    % r_cap
    2*10^(9),...    % r_del
    8*10^(10));     % k_rel

opts.set_equation(1); % using preset #1 (see Options class)
opts.NTopt=5; % for input sequences

%% create experiment object
Experiment1=Experiment(opts,forminfile,"seq",2.5); % using the input sequence option and concentration of profilin actin of 2.5

%% set gating 
% (gating for Capu is 1, which is the default, so it does not need to be set)
Experiment1.set_gating("FHODB",0.5);
Experiment1.set_gating("BNI1P(pPD18)",0.5);
Experiment1.set_gating("BNI1P(pPD28)",0.5);
Experiment1.set_gating("BNI1P(pPD37)",0.5);
Experiment1.set_gating("BNI1P(pPD65)",0.5);
Experiment1.set_gating("BNI1P(pPB37)",0.5);
Experiment1.set_gating("BNI1P(pPC37)",0.5);
Experiment1.set_gating("BNI1P(pPB18)",0.5);
Experiment1.set_gating("BNI1P(pPC18)",0.5);

%% add experimental data
% format: (formin name, experimental value, type, name-value args)
if ertype=='SEM'
    Experiment1=Experiment1.add_data("FHODB",0.6519,'ratio',errminus=0.04767829547,errplus=0.04767829547,groups=["NTD data"]);
    Experiment1=Experiment1.add_data("CAPU",1,'ratio',errminus=0.04767829547,errplus=0.04767829547,groups=["NTD data"]); % these error bars are a guess!
elseif ertype=='SD' 
    Experiment1=Experiment1.add_data("FHODB",0.6519,'ratio',errminus=0.1168,errplus=0.1168,groups=["NTD data"]);
    Experiment1=Experiment1.add_data("CAPU",1,'ratio',errminus=0.1168,errplus=0.1168,groups=["NTD data"]); % these error bars are a guess!
end
Experiment1=Experiment1.add_data("BNI1P(pPD18)",9.554461,'double',errminus=0.46825234,errplus=0.46829044,groups=["Fig 3","Fig 4c"]);
Experiment1=Experiment1.add_data("BNI1P(pPD28)",7.642275,'double',errminus=0.68292357,errplus=0.72194778,groups=["Fig 3"]);
Experiment1=Experiment1.add_data("BNI1P(pPB37)",9.617512,'double',errminus=1.17073171,errplus=0.60709018,groups=["Fig 4a"]);
Experiment1=Experiment1.add_data("BNI1P(pPC37)",6.885804,'double',errminus=0.88888889,errplus=0.9105691,groups=["Fig 4a"]);
Experiment1=Experiment1.add_data("BNI1P(pPD37)",2.70152246,'double',errminus=0.5203252,errplus=0.49864499,groups=["Fig 3","Fig 4a"]);
Experiment1=Experiment1.add_data("BNI1P(pPD65)",2.471568,'double',errminus=0.33170574,errplus=0.33170573,groups=["Fig 3"]);
%Experiment1=Experiment1.add_data("BNI1P(pPB18)",1.486034,'double',errminus=0.9907,errplus=1.0019,groups=["Fig 4c"]);
%Experiment1=Experiment1.add_data("BNI1P(pPC18)",4.294227,'double',errminus=0.5773,errplus=0.6145,groups=["Fig 4c"]);

Experiment1=Experiment1.add_data("BNI1P(pPD18)",8.17234102,'double',errminus=0.917030689999999,errplus=0.897556700000001,groups=["Fig 3 5"]);
Experiment1=Experiment1.add_data("BNI1P(pPD28)",7.48941745,'double',errminus=1.0926777,errplus=1.11218982,groups=["Fig 3 5"]);
Experiment1=Experiment1.add_data("BNI1P(pPD37)",2.70895242,'double',errminus=0.52678865,errplus=0.50731465,groups=["Fig 3 5"]);
Experiment1=Experiment1.add_data("BNI1P(pPD65)",2.35773458,'double',errminus=0.46829046,errplus=0.56585096,groups=["Fig 3 5"]);


%% make plots
fig3=Experiment1.expdatabar(group="Fig 3");
fig4a=Experiment1.expdatabar(group="Fig 4a");
figNTD=Experiment1.expdatabar("log2",false,group="NTD data"); % uses log2 scale
