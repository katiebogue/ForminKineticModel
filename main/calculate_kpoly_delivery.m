%% 
% calculates polymerization rate and outputs bar charts for all 3 states
% for double and dimer shows individual k_poly for each filiment and overall/combined
% uses variables: fh1_length, pp_length_vec, p_occ1/2ab/3ab, opt1/2/3/4

%% Calculations

%calculates kcapture for all three states
kc1 = k_paf*c_PA*(1-p_occ1).*pp_length_vec';
k_cap1 = sum(kc1);
kc1x = kc1;

kc2a = k_paf*c_PA*(1-p_occ2a).*pp_length_vec';
kc2b = k_paf*c_PA*(1-p_occ2b).*pp_length_vec';
k_cap2 = sum(kc2a) + sum(kc2b);

kc3a = k_paf*c_PA*(1-p_occ3a).*pp_length_vec';
kc3b = k_paf*c_PA*(1-p_occ3b).*pp_length_vec';
k_cap3 = sum(kc3a) + sum(kc3b);

%calculates kdelivery for all three states
kd1 = k_pab*(1-p_occ1_0).*(1.0e33*p_r1/(27*6.022e23));
k_del1 = sum(kd1);
kd1x = kd1;

kd2a = k_pab*(1-p_occ2a_0).*(1.0e33*p_r2a/(27*6.022e23));
kd2b = k_pab*(1-p_occ2b_0).*(1.0e33*p_r2b/(27*6.022e23));
k_del2 = sum(kd2a) + sum(kd2b);

kd3a = k_pab*(1-p_occ3a_0).*(1.0e33*p_r3a/(27*6.022e23));
kd3b = k_pab*(1-p_occ3b_0).*(1.0e33*p_r3b/(27*6.022e23));
k_del3 = sum(kd3a) + sum(kd3b);

%calculates kpoly for all three states
% for 3 states:
% kpoly= (1/kd + (kd + kc_rev)/(kd * kc))^(-1)
% for 4 states:
% kpoly= 1/r_PF_rev + (kd_rev + r_PF_rev)/(kd * r_PF_rev) + (kc_rev * kd_rev + kc_rev * r_PF_rev + kd * r_PF_rev)/(kc * kd * r_PF_rev)

%single
if del_state == 3
    kp1 = (1./kd1) + ((kd1 + k_paf_rev)./(kd1.*kc1)); %find kpoly for each individual PRM
end
if del_state == 4
    kp1 = (1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./(kd1 .* r_PF_rev)) + (((k_paf_rev .* r_paf_rev) + (k_paf_rev .* r_PF_rev) + (kd1 .* r_PF_rev))./(kc1 .* kd1 .* r_PF_rev)); %find kpoly for each individual PRM
end
kp1 = 1./kp1; %inverse above
k_poly1 = sum(kp1); %sum up the kpolys for all the PRMs
kp1x = kp1;

%double
if del_state == 3
    kp2a = (1./kd2a) + ((kd2a + k_paf_rev)./(kd2a.*kc2a));
    kp2b = (1./kd2b) + ((kd2b + k_paf_rev)./(kd2b.*kc2b));
end
if del_state == 4
    kp2a = (1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./(kd2a .* r_PF_rev)) + (((k_paf_rev .* r_paf_rev) + (k_paf_rev .* r_PF_rev) + (kd2a .* r_PF_rev))./(kc2a .* kd2a .* r_PF_rev)); %find kpoly for each individual PRM
    kp2b = (1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./(kd2b .* r_PF_rev)) + (((k_paf_rev .* r_paf_rev) + (k_paf_rev .* r_PF_rev) + (kd2b .* r_PF_rev))./(kc2a .* kd2b .* r_PF_rev)); %find kpoly for each individual PRM
end
kp2a = 1./kp2a; %inverse above
kp2b = 1./kp2b; %inverse above
k_poly2 = sum(kp2a) + sum(kp2b);  %sum up the kpolys for all the PRMs

%dimer
if del_state == 3
    kp3a = (1./kd3a) + ((kd3a + k_paf_rev)./(kd3a.*kc3a));
    kp3b = (1./kd3b) + ((kd3b + k_paf_rev)./(kd3b.*kc3b));
end
if del_state == 4
    kp3a = (1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./(kd3a .* r_PF_rev)) + (((k_paf_rev .* r_paf_rev) + (k_paf_rev .* r_PF_rev) + (kd3a .* r_PF_rev))./(kc3a .* kd3a .* r_PF_rev)); %find kpoly for each individual PRM
    kp3b = (1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./(kd3b .* r_PF_rev)) + (((k_paf_rev .* r_paf_rev) + (k_paf_rev .* r_PF_rev) + (kd3b .* r_PF_rev))./(kc3a .* kd3b .* r_PF_rev)); %find kpoly for each individual PRM
end
kp3a = 1./kp3a; %inverse above
kp3b = 1./kp3b; %inverse above
k_poly3 = sum(kp3a) + sum(kp3b);  %sum up the kpolys for all the PRMs

% CURRENTLY UNUSED
% calculates ratio of kpoly dimer to kpoly double
kp_ratio = k_poly3 / k_poly2;


%% Graphing
% creates polymerization bar charts for all three states

max_kp = max([k_cap1, k_cap2, k_cap3]) + 5;

    
%single
kplot1 = [kc1'; kd1'; kp1'];
fil=[1 2 3];

%formatting options based on opt4
if opt4 == 1
    if LOOP == 1
        subplot(3,4,2) %3x4 grid w/ axis on 2nd cell
    else
        if rem(LOOP,3) == 1
            subplot(3,4,2) %3x4 grid w/ axis on 2nd cell
        elseif rem(LOOP,3) == 2
            subplot(3,4,6) %3x4 grid w/ axis on 6th cell
        elseif rem(LOOP,3) == 0
            subplot(3,4,10) %3x4 grid w/ axis on 10th cell
        end
    end
else
    subplot(1,4,2) %1x4 grid w/ axis on 2nd cell
end

single_bar = bar(fil,kplot1,0.5, 'stacked')
title('Single');
%xlim([0.5 1.5]);
xlabel('Filaments');
ylabel('k');
xticklabels({'k_{cap}', 'k_{del}' , 'k_{poly}'});

for i = 1:length(kp1)
    single_bar(i).FaceColor = colors(i);
end

%formatting options based on opt1
if opt1 == 1
    ylim([0,max_kp])
elseif opt1 ==2
    if fh1_length <= 300
        ylim([1,275])
    else
        ylim([1,max_kp])
    end 
else
    ylim([0,580])
end

%double
kplot2 = [(kc2a+kc2b)'; (kd2a+kd2b)'; (kp2a+kp2b)'];
fil=[1 2 3];

% formatting options based on opt4
if opt4 == 1
    if LOOP == 1
        subplot(3,4,3) %3x4 grid w/ axis on 3rd cell
    else
        if rem(LOOP,3) == 1
            subplot(3,4,3) %3x4 grid w/ axis on 3rd cell
        elseif rem(LOOP,3) == 2
            subplot(3,4,7) %3x4 grid w/ axis on 7th cell
        elseif rem(LOOP,3) == 0
            subplot(3,4,11) %3x4 grid w/ axis on 11th cell
        end
    end
else
    subplot(1,4,3) %1x4 grid w/ axis on 3rd cell
end

double_bar = bar(fil, kplot2 ,0.5, 'stacked')
title('Double-Filament');
xlabel('Filaments');
ylabel('k');
xticklabels({'k_{cap}', 'k_{del}' , 'k_{poly}'});

for i = 1:length(kp1)
    double_bar(i).FaceColor = colors(i);
end

%formatting options based on opt1
if opt1 == 1
    ylim([0,max_kp])
elseif opt1 ==2
    if fh1_length <= 300
        ylim([1,275])
    else
        ylim([1,max_kp])
    end 
else
    ylim([0,580])
end

%dimer
kplot3 = [(kc3a+kc3b)'; (kd3a+kd3b)'; (kp3a+kp3b)'];
fil=[1, 2, 3];

%formatting options based on opt4
if opt4 == 1
    if LOOP == 1
        subplot(3,4,4) %3x4 grid w/ axis on 4th cell
    else
        if rem(LOOP,3) == 1
            subplot(3,4,4) %3x4 grid w/ axis on 4th cell
        elseif rem(LOOP,3) == 2
            subplot(3,4,8) %3x4 grid w/ axis on 8th cell
        elseif rem(LOOP,3) == 0
            subplot(3,4,12) %3x4 grid w/ axis on 12th cell
        end
    end
else
    subplot(1,4,4) %1x4 grid w/ axis on 4th cell
end

dimer_bar = bar(fil, kplot3 ,0.5, 'stacked')
title('N-Dimerized');
xlabel('Filaments');
ylabel('k');
xticklabels({'k_{cap}', 'k_{del}' , 'k_{poly}'});

for i = 1:length(kp1)
    dimer_bar(i).FaceColor = colors(i);
end

%formatting options based on opt1
if opt1 == 1
    ylim([0,max_kp])
elseif opt1 ==2
    if fh1_length <= 300
        ylim([1,275])
    else
        ylim([1,max_kp])
    end 
else
    ylim([0,580])
end
