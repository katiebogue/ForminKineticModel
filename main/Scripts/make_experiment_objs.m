% this is a script that...

ltfile="N600_lookup.mat"; % output file from polymer-c; must be on matlab path
pythonpath="/Users/katiebogue/MATLAB/GitHub/ForminKineticModel/main/python"; % path to python files
resultsloc="/Users/katiebogue/MATLAB/GitHub/Data/ForminKineticmodel_data/Results"; % path to location to save results
forminfile="Experimental constructs.txt"; % file containing sequences

ltfile2="basesep_35.5_output.mat"; % need to use base sep 35.5 runs for capu and fhod dimer

% create lookuptable
lt=(load(ltfile,'lookuptable').lookuptable);
lt=Lookuptable(lt);
lt2=(load(ltfile2,'lookuptable').lookuptable);
lt.addN(lt2.dimer.N123);
lt.addN(lt2.dimer.N255);

% create options
opts=Options(lt,pythonpath,"4st",resultsloc,50,4*10^(-1.9),10^(5.1),2*10^(9),8*10^(10));
opts.set_equation(1);
opts.NTopt=5; % for input sequences

% create experiment
Experiment1=Experiment(opts,forminfile,"seq",2.5);

% set gating
Experiment1.set_gating("FHODB",0.2);
Experiment1.set_gating("BNI1P(pPD18)",0.5);
Experiment1.set_gating("BNI1P(pPD28)",0.5);
Experiment1.set_gating("BNI1P(pPD37)",0.5);
Experiment1.set_gating("BNI1P(pPD65)",0.5);
Experiment1.set_gating("BNI1P(pPB37)",0.5);
Experiment1.set_gating("BNI1P(pPC37)",0.5);
Experiment1.set_gating("BNI1P(pPB18)",0.5);
Experiment1.set_gating("BNI1P(pPC18)",0.5);

% add experimental data
Experiment1=Experiment1.add_data("FHODB",0.6519,'ratio',errminus=0.1168,errplus=0.1168,groups=["NTD data"]);
Experiment1=Experiment1.add_data("CAPU",1,'ratio',errminus=0.1168,errplus=0.1168,groups=["NTD data"]); % these error bars are a guess!
Experiment1=Experiment1.add_data("BNI1P(pPD18)",9.554461,'double',errminus=0.46825234,errplus=0.46829044,groups=["Fig 3","Fig 4c"]);
Experiment1=Experiment1.add_data("BNI1P(pPD28)",7.642275,'double',errminus=0.68292357,errplus=0.72194778,groups=["Fig 3"]);
Experiment1=Experiment1.add_data("BNI1P(pPB37)",9.617512,'double',errminus=1.17073171,errplus=0.60709018,groups=["Fig 4a"]);
Experiment1=Experiment1.add_data("BNI1P(pPC37)",6.885804,'double',errminus=0.88888889,errplus=0.9105691,groups=["Fig 4a"]);
Experiment1=Experiment1.add_data("BNI1P(pPD37)",2.70152246,'double',errminus=0.5203252,errplus=0.49864499,groups=["Fig 3","Fig 4a"]);
Experiment1=Experiment1.add_data("BNI1P(pPD65)",2.471568,'double',errminus=0.33170574,errplus=0.33170573,groups=["Fig 3"]);
Experiment1=Experiment1.add_data("BNI1P(pPB18)",1.486034,'double',errminus=0.9907,errplus=1.0019,groups=["Fig 4c"]);
Experiment1=Experiment1.add_data("BNI1P(pPC18)",4.294227,'double',errminus=0.5773,errplus=0.6145,groups=["Fig 4c"]);


% make plots

fig3=Experiment1.expdatabar(group="Fig 3");
fig4a=Experiment1.expdatabar(group="Fig 4a");
figNTD=Experiment1.expdatabar("log2",false,group="NTD data");
