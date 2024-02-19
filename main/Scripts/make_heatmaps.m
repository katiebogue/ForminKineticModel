% MAKE_HEATMAPS generates and saves kpoly heatmap figures
    % 
    % ... 
    %
    % See also .

%% Set variables
ltfile="N600_lookup.mat"; % output file from polymer-c; must be on matlab path

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

opts.set_equation(1); % using preset #1 (see Options class)

%% Make heatmaps
newheatmap(opts,1,10,1,true,0.58)
newheatmap(opts,100,0.01,1,true,0.58)
newheatmap(opts,1000,0.01,1,true,0.58)

function newheatmap(opts,kcap,kdel,rcap,smooth,smoothfac)
    opts.update_results_folder
    opts.resultsfolder=strcat(opts.resultsfolder,"heatmap",num2str(kcap));
    opts.k_cap=kcap;
    opts.k_del=kdel;
    opts.r_cap=rcap;
    kpolyheatmap(opts,"NT dist v r_cap",1,smooth,smoothfac,true);
end