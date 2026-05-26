% MAKEFIGSFORPAPER makes figures for 2025 Bogue et. al paper.
%
% Settings and file locations for each section should be updated for each
% run.

temptf=true; % is this a temp one or the final one?
%% make lookup tables 
cd('/Users/katiebogue/MATLAB/GitHub/PolymerData') 
addpath '/Users/katiebogue/MATLAB/GitHub/polymer-c/Analysis'
%makeLookupMat('/Users/katiebogue/MATLAB/GitHub/Data/polymer-c_data/larger_actin_runs/FH1_35_updated',2,"largeractin_35_dimdob_temp.mat")
makeLookupMat('/Users/katiebogue/MATLAB/GitHub/Data/polymer-c_data/larger_actin_runs/FH1_1_updated',2,"largeractin_1_dimdob_temp.mat")
makeLookupMat('/Users/katiebogue/MATLAB/GitHub/Data/polymer-c_data/larger_actin_runs/FH1_16_updated',2,"largeractin_16_dimdob_temp.mat")

%% create lookup table files
cd('/Users/katiebogue/MATLAB/GitHub/PolymerData')
load('largeractin_1_dimdob_temp.mat')
lt_1=lookuptable;
load('largeractin_16_dimdob_temp.mat')
lt_16=lookuptable;
load('largeractin_35_dimdob_temp.mat')
lt_35=lookuptable;
load('largeractin_355.mat')
lt_1.single=lookuptable.single;
lt_16.single=lookuptable.single;
lt_35.single=lookuptable.single;
clear lookuptable
lookup_1=lt_1;
lookup_16=lt_16;
lookup_35=lt_35;
lt_1=Lookuptable(lookup_1);
lt_16=Lookuptable(lookup_16);
lt_35=Lookuptable(lookup_35);
cd('/Users/katiebogue/MATLAB/GitHub/MatLab_Workspaces/')
if temptf
    save("largeractinLookuptabs_temp.mat")
else
    save("largeractinLookuptabs.mat")
end
%% make polymerheatmaps
saveTF=false;
time= datetime('now', 'Format','yyyy-MM-dd HH-mm');
time= string(time);
resultsfolder= 'RESULTS_' + time;
resultsfolder= strcat("/Users/katiebogue/MATLAB/GitHub/Data/files_for_comp_paper/",resultsfolder);
mkdir(fullfile(resultsfolder,"polymerheatmaps"))
savefigfolder=fullfile(resultsfolder,"polymerheatmaps");
polymerstatheatmap("POcclude",lt_16,"16.667","ratio",false,saveTF,savefigfolder,[-1 1])
polymerstatheatmap("POcclude",lt_35,"35.5","ratio",false,saveTF,savefigfolder,[-1 1])
polymerstatheatmap("Prvec0",lt_35,"35.5","ratio",false,saveTF,savefigfolder,[-3 3])
polymerstatheatmap("POcclude",lt_16,"16.667","ratio",false,saveTF,savefigfolder,[-1 1])
polymerstatheatmap("Prvec0",lt_16,"16.667","ratio",false,saveTF,savefigfolder,[-3 3])
polymerstatheatmap("POcclude",lt_1,"1","ratio",false,saveTF,savefigfolder,[-1 1])
polymerstatheatmap("Prvec0",lt_1,"1","ratio",false,saveTF,savefigfolder,[-3 3])
load("forminNTDsweepexp.mat")
polymerstatheatmap("POcclude",lt_35,"35.5","ratio",false,saveTF,savefigfolder,[-1 1],Experiment1)
polymerstatheatmap("Prvec0",lt_35,"35.5","ratio",false,saveTF,savefigfolder,[-3 3],Experiment1)
%% make combo sweeps (note-- need to update this function if the mcmc best fit changes)
makecombosweeps(temptf,resultsfolder)

%% update forminexperimentobjs.mat
cd('/Users/katiebogue/MATLAB/GitHub/MatLab_Workspaces/')
load("forminexperimentobjs.mat")
opts.lookup=lt_35;
save("forminexperimentobjs.mat","opts","Experiment_FHOD","Experiment_BNI1","Experiment_Capu")

%% update BNI1 histograms in mcmc folder
makeBNI1NTDmat('/Users/katiebogue/MATLAB/GitHub/Data/files_for_comp_paper/MCMC_RESULTS_2025-09-10 15-06_3st_nondim0_prcalc0_prfit0_errtype3_fitrexp1_BNI1fit')

%% make formin prediction figures
makeforminntdsweeps(lt_35,resultsfolder)