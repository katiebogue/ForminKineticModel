%%
% finds probability of occlusion for all three states from simulation
% generates:
    % 3 column array with length / binding site location / probability of
    % occlusion for all possible binding sites along simulated fh1 lengths
    % scattered interpolant object

% the array (X1, X2a/b, X3a/b) contains exact probabilities
% the interpolants (F1, F2a/b, F3a/b) contain interpolated and extrapolated
% probabilities
% both are used later to calculate polymerization rate
% see (wholepackage (intro)) for info about which one is used

% numErr gives number of errors found in simulation
% i.e. some x > 1 is recorded as being probability of occlusion

%% single

numErr = 0;

%initialization
NFil1=1;
p_occ1_all = [];
p_r1_all = [];
length_array1 = 1:300; %length of specific filament from sweep
length_vec1 = [];
iSite_vec1 = [];

%POcc1(iN, isite)=[occlusion probability of the iNth N value at isite]
for iN=length_array1
    N_All = NFil1*iN;
    for iy =1:iN
        length_vec1 = [length_vec1 iN]; %list of filament lengths
        iSite_vec1 = [iSite_vec1 iy]; %all iSites used
        p_occ1_all = [p_occ1_all m1(iN, 16 + 2*(N_All +1) + 7*(iy - 1))]; %lookup table 1
        p_r1_all = [p_r1_all m1(iN, 19 + 2*(N_All +1) + 7*(iy - 1))]; %lookup table 1
        
        if p_occ1_all(end) > 1
            numErr = numErr+1;
        end
    end
end

X1 = [length_vec1' iSite_vec1' p_occ1_all' p_r1_all'];
F1 = scatteredInterpolant(X1(:,1),X1(:,2),X1(:,3),'linear','nearest');
FF1 = scatteredInterpolant(X1(:,1),X1(:,2),X1(:,4),'linear','nearest');

%% double

NFil = 2;
p_occ2a_all = [];
p_occ2b_all = [];
p_r2a_all = [];
p_r2b_all = [];
length_array2 = 1:200;
length_vec2 = [];
iSite_vec2 = [];


%POcc1(iN, isite)=[occlusion probability of the iNth N value at isite]
for iN=length_array2
    N_All =NFil*iN;
    for iy =1:iN
        p_a = m2(iN, 16 + 2*(N_All +1) + 7*(iy - 1)); %lookup table 2
        p_b = m2(iN, 16 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil));
        r_a = m2(iN, 19 + 2*(N_All +1) + 7*(iy - 1)); %lookup table 2
        r_b = m2(iN, 19 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil));
        if p_a > 1 || p_b >1
            numErr = numErr+1;
        end
        length_vec2 = [length_vec2 iN]; %list of filament lengths
        iSite_vec2 = [iSite_vec2 iy]; %all iSites used
        p_occ2a_all = [p_occ2a_all p_a];
        p_occ2b_all = [p_occ2b_all p_b];
        p_r2a_all = [p_r2a_all r_a];
        p_r2b_all = [p_r2b_all r_b];
    end
end

X2a = [length_vec2; iSite_vec2; p_occ2a_all; p_r2a_all]' ;
F2a = scatteredInterpolant(X2a(:,1),X2a(:,2),X2a(:,3),'linear','nearest');
FF2a = scatteredInterpolant(X2a(:,1),X2a(:,2),X2a(:,4),'linear','nearest');


X2b = [length_vec2; iSite_vec2; p_occ2b_all; p_r2b_all]';
F2b = scatteredInterpolant(X2b(:,1),X2b(:,2),X2b(:,3),'linear','nearest');
FF2b = scatteredInterpolant(X2b(:,1),X2b(:,2),X2b(:,4),'linear','nearest');

%% dimer

NFil=2;
p_occ3a_all = [];
p_occ3b_all = [];
p_r3a_all = [];
p_r3b_all = [];
length_array3 = 1:122;
length_vec3 = [];
iSite_vec3 = [];

%POcc1(iN, isite)=[occlusion probability of the iNth N value at isite]
for iN=length_array3
    N_All =NFil*iN;
    for iy =1:iN
        length_vec3 = [length_vec3 iN]; %list of filament lengths
        iSite_vec3 = [iSite_vec3 iy]; %all iSites used in percent of total filament length
        p_occ3a_all = [p_occ3a_all m3(iN, 16 + 2*(N_All +1) + 7*(iy - 1))];
        p_occ3b_all = [p_occ3b_all m3(iN, 16 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil))];
        p_r3a_all = [p_r3a_all m3(iN, 19 + 2*(N_All +1) + 7*(iy - 1))];
        p_r3b_all = [p_r3b_all m3(iN, 19 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil))];
    end
end


X3a = [length_vec3' iSite_vec3' p_occ3a_all' p_r3a_all'];
X3b = [length_vec3' iSite_vec3' p_occ3b_all' p_r3b_all'];

F3a = scatteredInterpolant(X3a(:,1),X3a(:,2),X3a(:,3),'linear','nearest');
F3b = scatteredInterpolant(X3b(:,1),X3b(:,2),X3b(:,3),'linear','nearest');

FF3a = scatteredInterpolant(X3a(:,1),X3a(:,2),X3a(:,4),'linear','nearest');
FF3b = scatteredInterpolant(X3b(:,1),X3b(:,2),X3b(:,4),'linear','nearest');


if numErr > 0
    disp('simulation probability error')
end
