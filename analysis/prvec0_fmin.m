% PRVEC0_FMIN Run fminsearch 1000 times for each prvec0 option.
%
    %   Runs Experiment.runfminsearch 1000 times for each prvec0 option
    %   specified in the script (this can be modified).
    %
    %   Be sure to set an experiment object with appropriate lookuptable
    %   and data added
    %
    %   Adds the outputs to tables_capudimerdist_3st and
    %   tables_capudimerdist_4st_krelscale, table variables that should
    %   already exist before running the script.
    %
    % See also FORMIN, EXPERIMENT, LOOKUPTABLE, OPTIONS, EXPDATABAR,
    % KPOLYMERIZATION, PRM, EXPERIMENT/RUNFMINSEARCH.

%% Set Experiment
Experiment1=Experiment1;

%% Run fminsearch 1000 times for each Prvec0 option

tables.Prvec0=Experiment1.runfminsearch("iterations",1000,"groups",{"NTD data","Fig 3","Fig 4a"},"weight","NTD data");

%prvec_opts={"Prvec0_halfup","Prvec0_halfup_op","Prvec0_op","Prvec0_up","Prvec0_up_op","Prvec_cen","Prvec_cen_halfup","Prvec_cen_up","Prvec_offcen","Prvec_offcen_halfup","Prvec_offcen_halfup_op","Prvec_offcen_op","Prvec_offcen_up","Prvec_offcen_up_op"};

prvec_opts={"Prvec_offcen","Prvec0_halfup","Prvec_offcen_halfup"};

Experiment1.opts.set_equation(0,"krel",{"size","negexp"})
Experiment1.opts.kpoly_type="4st";
Experiment1.opts.set_equation(0,"kdel",{"POcclude","1-base","Prvec0","amino","gating","linear"})
tables_capudimerdist_3st.Prvec0=Experiment1.runfminsearch("iterations",1000,"groups",{"NTD data","Fig 3","Fig 4a"},"weight","NTD data");

for i=1:length(prvec_opts)
    disp(prvec_opts{i})
    Experiment1.opts.set_equation(0,"kdel",{"POcclude","1-base",prvec_opts{i},"amino","gating","linear"})
    tables_capudimerdist_4st_krelscale.(prvec_opts{i})=Experiment1.runfminsearch("iterations",1000,"groups",{"NTD data","Fig 3","Fig 4a"},"weight","NTD data");
end

