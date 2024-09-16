load('Experiments.mat')

y=MCMCParamfit(Experiment1,"4st",1,0,1000,50,1e8,0.01);

y=MCMCParamfit(Experiment1,"3st",1,0,1000,50,1e8,0.01);