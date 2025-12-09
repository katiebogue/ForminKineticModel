function makecombosweeps(tempTF,resultsfolder)
%MAKECOMBOSWEEPS generates kpoly ratio heatmaps based on existing .mat file
%lookuptables for 3 different parameter regiemes
%
%   MAKECOMBOSWEEPS(tempTF,resultsfolder) 
%
%   Inputs:
%       tempTF        : (Bool) whether to use the lookup table
%       largeractinLookuptabs_temp.mat (true) or largeractinLookuptabs.mat
%       (false)
%       resultsfolder : (String) location to save figures to
% 
%  The lookuptable file is hardcoded and must be manually changed.
%  The 3 parameter regimes are hardcoded and must be manually changed.
% 
% See also OPTIONS, KPOLYHEATMAP.


    pythonpath="/Users/katiebogue/MATLAB/GitHub/ForminKineticModel/main/python"; % path to python files
    resultsloc=resultsfolder;
    if tempTF
        load("largeractinLookuptabs_temp.mat")
    else
        load("largeractinLookuptabs.mat")
    end
    
    opts_1=Options(lt_1,pythonpath,...
        "3st",...       % kpoly type
        resultsloc,...
        10^3.2778,...          % k_cap
        10^-2.8678,... % k_del
        10^3.15,...    % r_cap
        1,...    % r_del
        1);     % k_rel
    
    opts_1.r_cap_exp=0.57121;
    opts_1.set_equation(1); % using preset #1 (see Options class)
    opts_1.NTopt=2; 
    
    opts_16=Options(lt_16,pythonpath,...
        "3st",...       % kpoly type
        resultsloc,...
        10^3.2778,...          % k_cap
        10^-2.8678,... % k_del
        10^3.15,...    % r_cap
        1,...    % r_del
        1);     % k_rel
    
    opts_16.r_cap_exp=0.57121;
    opts_16.set_equation(1); % using preset #1 (see Options class)
    opts_16.NTopt=2;
    
    opts_35=Options(lt_35,pythonpath,...
        "3st",...       % kpoly type
        resultsloc,...
        10^3.2778,...          % k_cap
        10^-2.8678,... % k_del
        10^3.15,...    % r_cap
        1,...    % r_del
        1);     % k_rel
    
    opts_35.r_cap_exp=0.57121;
    opts_35.set_equation(1); % using preset #1 (see Options class)
    opts_35.NTopt=2;
    
    opts_35.resultsfolder=strcat(opts_35.resultsfolder,"combosweep_35");
    opts_16.resultsfolder=strcat(opts_16.resultsfolder,"combosweep_16");
    opts_1.resultsfolder=strcat(opts_1.resultsfolder,"combosweep_1");
    
    
    opts_1.k_cap=1;
    opts_1.k_del=1;
    opts_1.r_cap=10^4;
    opts_1.r_cap_exp=0.8;
    
    
    opts_16.k_cap=1;
    opts_16.k_del=1;
    opts_16.r_cap=10^4;
    opts_16.r_cap_exp=0.8;
    
    opts_35.k_cap=1;
    opts_35.k_del=1;
    opts_35.r_cap=10^4;
    opts_35.r_cap_exp=0.8;
    
    opts_1.k_cap=15000;
    opts_16.k_cap=15000;
    opts_35.k_cap=15000;
    
    [fig,h]=kpolyheatmap(opts_1,"NT dist v CT dist",[-3 3]);
    [fig,h]=kpolyheatmap(opts_16,"NT dist v CT dist",[-3 3]);
    [fig,h]=kpolyheatmap(opts_35,"NT dist v CT dist",[-3 3]);
    
    opts_1.k_cap=10^9;
    opts_16.k_cap=10^9;
    opts_35.k_cap=10^9;
    
    [fig,h]=kpolyheatmap(opts_1,"NT dist v CT dist",[-3 3]);
    [fig,h]=kpolyheatmap(opts_16,"NT dist v CT dist",[-3 3]);
    [fig,h]=kpolyheatmap(opts_35,"NT dist v CT dist",[-3 3]);
    
    
    opts_1.k_cap=1;
    opts_16.k_cap=1;
    opts_35.k_cap=1;
    
    
    opts_1.k_del=10^5;
    opts_16.k_del=10^5;
    opts_35.k_del=10^5;
    
    
    [fig,h]=kpolyheatmap(opts_1,"NT dist v CT dist",[-3 3]);
    [fig,h]=kpolyheatmap(opts_16,"NT dist v CT dist",[-3 3]);
    [fig,h]=kpolyheatmap(opts_35,"NT dist v CT dist",[-3 3]);
end

