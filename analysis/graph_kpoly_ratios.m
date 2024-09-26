% GRAPH_KPOLY_RATIOS creates bargraphs of dimer/double Kpoly ratios
% simulated for a list of formins
% 
% Reads in data from fh1_info_u200.mat; kp_ratiou200.mat; all_names.mat;
% kp_ratio_all.mat

%% for all with length under 200
% minimal extrapolation

%kpoly ratios are dimer / double
load('fh1_info_u200.mat');
load('kp_ratio_u200.mat');

X1 = categorical(fh1_info_u200);
c = parula(21); %?

figure()
hold on

for n = 1:21
    if kp_ratio_u200(n) < 1
        bar(X1(n), log10(kp_ratio_u200(n)),'b')
    else
        bar(X1(n), log10(kp_ratio_u200(n)),'r')
    end
end

ylim([-0.1,0.1])

% this should be changed but ah
ylabel('double faster    --  no change --     dimer faster')

title('ratio of dimer:double polymerization rates for formins under length 200')
hold off

%% all the formins

load('all_names.mat');
load('kp_ratio_all.mat');

X2 = categorical(all_names);
c = parula(25); %?

figure()
hold on

for n = 1:25
    if kp_ratio_all(n) < 1
        bar(X2(n), log10(kp_ratio_all(n)),'b')
    else
        bar(X2(n), log10(kp_ratio_all(n)),'r')
    end
end

ylim([-0.5,0.5])

ylabel('double faster    --  no change --     dimer faster')

title('ratio of dimer:double polymerization rates for all formins')
hold off