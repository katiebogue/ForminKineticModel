% BNI1_SWEEPS generates Lookuptable, Options, and Experiment...
    % 
    % 
    %
    % See also FORMIN, EXPERIMENT, LOOKUPTABLE, OPTIONS, EXPDATABAR,
    % KPOLYMERIZATION, PRM, EXPERIMENT/EXPDATABAR.

%% Input files/ paths
%ltfile="N600_lookup.mat"; % output file from polymer-c; must be on matlab path
ltfile="prvec_runs_lookup_N260.mat"; % output file from polymer-c; must be on matlab path

pythonpath="/Users/katiebogue/MATLAB/GitHub/ForminKineticModel/main/python"; % path to python files
resultsloc="/Users/katiebogue/MATLAB/GitHub/Data/ForminKineticmodel_data/Results"; % path to location to save results
forminfile="BNI1_FH2 to PRM.txt"; % file containing sequences
%% create lookuptable
lt=(load(ltfile,'lookuptable').lookuptable);
lt=Lookuptable(lt);
%% create options object
% modify this line to change the rate constants:
opts=Options(lt,pythonpath,...
    "4st",...       % kpoly type
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
Experiment1.set_gating("BNI1",0.5);

%% Load best fit tables
load("best fits.mat")

Experiment1.applytable(fit_4st)
Bni1=Experiment1.ForminList;
[doubles_4st,dimers_4st,ratios_4st,NTD_dists_4st]=NTD_predictions(Bni1,"4st");

Experiment1.applytable(fit_3st)
Bni1=Experiment1.ForminList;
[doubles_3st,dimers_3st,ratios_3st,NTD_dists_3st]=NTD_predictions(Bni1,"3st");

Experiment1.applytable(fit_4st_krel)
Bni1=Experiment1.ForminList;
[doubles_4st_krel,dimers_4st_krel,ratios_4st_krel,NTD_dists_4st_krel]=NTD_predictions(Bni1,"4st_krel");

function [doubles,dimers,ratios,NTD_dists]=NTD_predictions(formin,fit)
    doubles=zeros(1,489);
    dimers=zeros(1,489);
    ratios=zeros(1,489);
    NTD_dists=zeros(1,489);
    for i=0:488
        formin.add_length(1)
        kpoly=formin.kpoly;
        doubles(i+1)=kpoly.double;
        dimers(i+1)=kpoly.dimer;
        ratios(i+1)=kpoly.ratio;
        NTD_dists(i+1)=formin.PRMList(1,4).dist_NT;
    end
    formin.add_length(-488)

    figure
    index=ratios<1;
    scatter(NTD_dists(index),log2(ratios(index)),"cyan","filled")
    hold on
    index=ratios==1;
    scatter(NTD_dists(index),log2(ratios(index)),"yellow","filled")
    hold on
    index=ratios>1;
    scatter(NTD_dists(index),log2(ratios(index)),"magenta","filled")
    yline(0)
    xlabel("Dist from last PRM to NT")
    ylabel("log_{2}(Kpoly dimer/double")
    title("BNI1 NTD predictions")
    subtitle(strcat("best fit ",fit))
end


