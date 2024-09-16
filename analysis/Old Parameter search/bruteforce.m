% HPC3 slurm function script to brute force search parameter space

function bruteforce(outfile, outloc)
% must run from the directory with bruteforce.mat

% Output, one line separated by comma:
% 1- regular SOS (no error bars, simulated prvec0)
% 2- SOS with error bars, simulated prvec0
% 3- SOS with error bars, calculated prvec0
% 4- NTD SOS (no error bars, simulated prvec0)
% 5- NTD SOS error bars, simulated prvec0
% 6- NTD SOS error bars, calculated prvec0
% 7- fig3 regular SOS (no error bars, simulated prvec0)
% 8- fig3 SOS error bars, simulated prvec0
% 9- fig3 SOS error bars, calculated prvec0
% 10- fig4 regular SOS (no error bars, simulated prvec0)
% 11- fig4 SOS error bars, simulated prvec0
% 12- fig4 SOS error bars, calculated prvec0
% 13- k_paf
% 14- k_pab
% 15- k_paf_rev
% 16- r_PF_rev
% 17- r_paf_rev
settings
load("bruteforce.mat") % file with workspace
rng shuffle

r=8*rand(1,5)+-3; %random numbers between -3 and 8
vals={10^r(1),10^r(2),10^r(3),10^r(4),10^r(5)};
outputs=[SOS(vals{:});SOS_err(vals{:});SOS_errc(vals{:});SOSNTD(vals{:});SOSNTD_err(vals{:});SOSNTD_err_cr(vals{:}) ;SOS3(vals{:});SOS3_err(vals{:});SOS3_errc(vals{:}) ;SOS4a(vals{:}) ;SOS4a_err(vals{:});SOS4a_errc(vals{:});vals{1} ;vals{2} ;vals{3};vals{4};vals{5}];

cd(outloc);
file=fopen(outfile,'w');
fprintf(file,'%f',outputs(1));
fprintf(file,',%f',outputs(2:end));

end