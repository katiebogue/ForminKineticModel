% MAKE_HEATMAPS generates and saves kpoly heatmap figures
    % 
    % ... 
    %
    % See also .

%% Set variables
ltfile="prvec_runs_lookup.mat"; % output file from polymer-c; must be on matlab path
%ltfile="N600_lookup.mat";

pythonpath="/Users/katiebogue/MATLAB/GitHub/ForminKineticModel/main/python"; % path to python files
resultsloc="/Users/katiebogue/MATLAB/GitHub/Data/ForminKineticmodel_data/Results"; % path to location to save results


%% create lookuptable
lt=(load(ltfile,'lookuptable').lookuptable);
lt=Lookuptable(lt);

%% create options object
% modify this line to change the rate constants:
opts=Options(lt,pythonpath,...
    "3st",...       % kpoly type
    resultsloc,...
    1,...          % k_cap
    1,... % k_del
    1,...    % r_cap
    1,...    % r_del
    1);     % k_rel
opts.r_cap_exp=1;

opts.set_equation(1); % using preset #1 (see Options class)

%% Make heatmaps
opts.update_results_folder
opts.resultsfolder=strcat(opts.resultsfolder,"heatmap");
newheatmap(opts,100000,0.01,1,true,0.58)
newheatmap(opts,1,0.01,1,true,0.58)
newheatmap(opts,100,0.01,1,true,0.58)
newheatmap(opts,1000,0.01,1,true,0.58)
newheatmap(opts,10000,0.01,1,true,0.58)

function newheatmap(opts,kcap,kdel,rcap,smooth,smoothfac)
    %opts.update_results_folder
    %opts.resultsfolder=strcat(opts.resultsfolder,"heatmap",num2str(kcap));
    opts.k_cap=kcap;
    opts.k_del=kdel;
    opts.r_cap=rcap;
    kpolyheatmap(opts,"NT dist v r_cap",[-2 2],smooth,smoothfac,true);
end