load("largeractinLookuptabs.mat")

saveTF=true;
savefigfolder="/Users/katiebogue/MATLAB/GitHub/Data/files_for_comp_paper/polymerheatmaps_temp";

polymerstatheatmap("POcclude",lt_35,"35.5","ratio",false,saveTF,savefigfolder,[-1 1])
polymerstatheatmap("Prvec0",lt_35,"35.5","ratio",false,saveTF,savefigfolder,[-3 3])

polymerstatheatmap("POcclude",lt_16,"16.667","ratio",false,saveTF,savefigfolder,[-1 1])
polymerstatheatmap("Prvec0",lt_16,"16.667","ratio",false,saveTF,savefigfolder,[-3 3])

polymerstatheatmap("POcclude",lt_1,"1","ratio",false,saveTF,savefigfolder,[-1 1])
polymerstatheatmap("Prvec0",lt_1,"1","ratio",false,saveTF,savefigfolder,[-3 3])

load("forminNTDsweepexp.mat")
polymerstatheatmap("POcclude",lt_35,"35.5","ratio",false,saveTF,savefigfolder,[-1 1],Experiment1)
polymerstatheatmap("Prvec0",lt_35,"35.5","ratio",false,saveTF,savefigfolder,[-3 3],Experiment1)