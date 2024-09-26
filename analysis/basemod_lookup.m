% BASEMOD_LOOKUP collects prvec0 values from lookuptables for base
% separation distances 10 and 20, and makes scatterplots for the prvec0
% ratio (dimer/double) vs. FH1 length for each table
    % 
    % Uses the tables generated from the old format. The table files used
    % are double_basesep10_1,10,20-250.txt,
    % dimer_basesep10_1,10,20-250.txt, double_basesep20_1,10,20-250.txt,
    % and dimer_basesep20_1,10,20-250.txt 

%% base sep 10
numErr = 0; 
m2 = dlmread('double_basesep10_1,10,20-250.txt');
m3 = dlmread('dimer_basesep10_1,10,20-250.txt');
% double

NFil = 2;
p_occ2a_bs10 = [];
p_occ2b_bs10 = [];
p_r2a_bs10 = [];
p_r2b_bs10 = [];
length_array2 = 10:10:250;
length_array2 = [1,length_array2];
length_vec2_bs10 = [];
iSite_vec2_bs10 = [];


%POcc1(iN, isite)=[occlusion probability of the iNth N value at isite]
i=0;
for iN=length_array2
    i=i+1;
    N_All =NFil*iN;
    for iy =1:iN
        p_a = m2(i, 16 + 2*(N_All +1) + 7*(iy - 1)); %lookup table 2
        p_b = m2(i, 16 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil));
        r_a = m2(i, 19 + 2*(N_All +1) + 7*(iy - 1)); %lookup table 2
        r_b = m2(i, 19 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil));
        if p_a > 1 || p_b >1
            numErr = numErr+1;
        end
        length_vec2_bs10 = [length_vec2_bs10 iN]; %list of filament lengths
        iSite_vec2_bs10 = [iSite_vec2_bs10 iy]; %all iSites used
        p_occ2a_bs10 = [p_occ2a_bs10 p_a];
        p_occ2b_bs10 = [p_occ2b_bs10 p_b];
        p_r2a_bs10 = [p_r2a_bs10 r_a];
        p_r2b_bs10 = [p_r2b_bs10 r_b];
    end
end

X2a_bs10 = [length_vec2_bs10; iSite_vec2_bs10; p_occ2a_bs10; p_r2a_bs10]' ;

X2b_bs10 = [length_vec2_bs10; iSite_vec2_bs10; p_occ2b_bs10; p_r2b_bs10]';

% dimer

NFil=2;
p_occ3a_bs10 = [];
p_occ3b_bs10 = [];
p_r3a_bs10 = [];
p_r3b_bs10 = [];
length_array3 = 10:10:250;
length_array3 = [1,length_array3];
length_vec3_bs10 = [];
iSite_vec3_bs10 = [];

%POcc1(iN, isite)=[occlusion probability of the iNth N value at isite]
i=0;
for iN=length_array3
    i=i+1;
    N_All =NFil*iN;
    for iy =1:iN
        length_vec3_bs10 = [length_vec3_bs10 iN]; %list of filament lengths
        iSite_vec3_bs10 = [iSite_vec3_bs10 iy]; %all iSites used in percent of total filament length
        p_occ3a_bs10 = [p_occ3a_bs10 m3(i, 16 + 2*(N_All +1) + 7*(iy - 1))];
        p_occ3b_bs10 = [p_occ3b_bs10 m3(i, 16 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil))];
        p_r3a_bs10 = [p_r3a_bs10 m3(i, 19 + 2*(N_All +1) + 7*(iy - 1))];
        p_r3b_bs10 = [p_r3b_bs10 m3(i, 19 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil))];
    end
end

X3a_bs10 = [length_vec3_bs10; iSite_vec3_bs10; p_occ3a_bs10; p_r3a_bs10]' ;

X3b_bs10 = [length_vec3_bs10; iSite_vec3_bs10; p_occ3b_bs10; p_r3b_bs10]';

if numErr > 0
    disp('simulation probability error')
end

%% base sep 20
numErr = 0; 
m2 = dlmread('double_basesep20_1,10,20-250.txt');
m3 = dlmread('dimer_basesep20_1,10,20-250.txt');
% double

NFil = 2;
p_occ2a_bs20 = [];
p_occ2b_bs20 = [];
p_r2a_bs20 = [];
p_r2b_bs20 = [];
length_array2 = 10:10:250;
length_array2 = [1,length_array2];
length_vec2_bs20 = [];
iSite_vec2_bs20 = [];


%POcc1(iN, isite)=[occlusion probability of the iNth N value at isite]
i=0;
for iN=length_array2
    i=i+1;
    N_All =NFil*iN;
    for iy =1:iN
        p_a = m2(i, 16 + 2*(N_All +1) + 7*(iy - 1)); %lookup table 2
        p_b = m2(i, 16 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil));
        r_a = m2(i, 19 + 2*(N_All +1) + 7*(iy - 1)); %lookup table 2
        r_b = m2(i, 19 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil));
        if p_a > 1 || p_b >1
            numErr = numErr+1;
        end
        length_vec2_bs20 = [length_vec2_bs20 iN]; %list of filament lengths
        iSite_vec2_bs20 = [iSite_vec2_bs20 iy]; %all iSites used
        p_occ2a_bs20 = [p_occ2a_bs20 p_a];
        p_occ2b_bs20 = [p_occ2b_bs20 p_b];
        p_r2a_bs20 = [p_r2a_bs20 r_a];
        p_r2b_bs20 = [p_r2b_bs20 r_b];
    end
end

X2a_bs20 = [length_vec2_bs20; iSite_vec2_bs20; p_occ2a_bs20; p_r2a_bs20]' ;

X2b_bs20 = [length_vec2_bs20; iSite_vec2_bs20; p_occ2b_bs20; p_r2b_bs20]';

% dimer

NFil=2;
p_occ3a_bs20 = [];
p_occ3b_bs20 = [];
p_r3a_bs20 = [];
p_r3b_bs20 = [];
length_array3 = 10:10:250;
length_array3 = [1,length_array3];
length_vec3_bs20 = [];
iSite_vec3_bs20 = [];

%POcc1(iN, isite)=[occlusion probability of the iNth N value at isite]
i=0;
for iN=length_array3
    i=i+1;
    N_All =NFil*iN;
    for iy =1:iN
        length_vec3_bs20 = [length_vec3_bs20 iN]; %list of filament lengths
        iSite_vec3_bs20 = [iSite_vec3_bs20 iy]; %all iSites used in percent of total filament length
        p_occ3a_bs20 = [p_occ3a_bs20 m3(i, 16 + 2*(N_All +1) + 7*(iy - 1))];
        p_occ3b_bs20 = [p_occ3b_bs20 m3(i, 16 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil))];
        p_r3a_bs20 = [p_r3a_bs20 m3(i, 19 + 2*(N_All +1) + 7*(iy - 1))];
        p_r3b_bs20 = [p_r3b_bs20 m3(i, 19 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil))];
    end
end

X3a_bs20 = [length_vec3_bs20; iSite_vec3_bs20; p_occ3a_bs20; p_r3a_bs20]' ;

X3b_bs20 = [length_vec3_bs20; iSite_vec3_bs20; p_occ3b_bs20; p_r3b_bs20]';


if numErr > 0
    disp('simulation probability error')
end

%% Plots

% base sep 10
figure
scatter(length_vec2_bs10,(p_r3a_bs10./p_r2a_bs10),"filled")
hold on
scatter(length_vec2_bs10,(p_r3b_bs10./p_r2b_bs10),"filled")
xlabel('FH1 Length')
ylabel('prvec0 dimer/double')
title('basesepdist=10')


% base sep 20
figure
scatter(length_vec2_bs20,(p_r3a_bs20./p_r2a_bs20),"filled")
hold on
scatter(length_vec2_bs20,(p_r3b_bs20./p_r2b_bs20),"filled")
xlabel('FH1 Length')
ylabel('prvec0 dimer/double')
title('basesepdist=20')