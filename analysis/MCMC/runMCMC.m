Exp=Experiment1;
exptype=1; % 1= is an experiment, 2= is a struct with rates, data, and resultsfolder and resultsdir, and fh1sizes and prmlocs
type='3st';
errtype=3; % 1= use separate sigma values, 2= divide each by SEM, 3=rmse
matfileTF=0; % whether to save to a matfile or to keep things in memory
NTCHECK = 3000;
NTADAPT =100;
NTMAX =3*10^5;
KSCRITICAL =0.02; %0.02
nondim= 1; % whether to use nondimensionality
prcalc= 0; 
prfit = 0;% whether or not to fit x and y delivery location as
xloc = 0;
yloc = 0;
titleadd = "";
fitrexp = 1;
plotsameexp = 1;

MCMCParamfit(Exp,exptype,type,errtype,matfileTF, NTCHECK,NTADAPT,NTMAX,KSCRITICAL,nondim,prcalc, prfit, xloc, yloc, titleadd, fitrexp,plotsameexp)

