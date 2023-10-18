%% base sep 35.5
m3 = dlmread('dimer_basesep35.5_123_255.txt');
% dimer

NFil=2;
p_occ3a_bs35 = [];
p_occ3b_bs35 = [];
p_r3a_bs35 = [];
p_r3b_bs35 = [];
length_array3 = [123,255];
length_vec3_bs35 = [];
iSite_vec3_bs35 = [];

%POcc1(iN, isite)=[occlusion probability of the iNth N value at isite]
i=0;
for iN=length_array3
    i=i+1;
    N_All =NFil*iN;
    for iy =1:iN
        length_vec3_bs35 = [length_vec3_bs35 iN]; %list of filament lengths
        iSite_vec3_bs35 = [iSite_vec3_bs35 iy]; %all iSites used in percent of total filament length
        p_occ3a_bs35 = [p_occ3a_bs35 m3(i, 16 + 2*(N_All +1) + 7*(iy - 1))];
        p_occ3b_bs35 = [p_occ3b_bs35 m3(i, 16 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil))];
        p_r3a_bs35 = [p_r3a_bs35 m3(i, 19 + 2*(N_All +1) + 7*(iy - 1))];
        p_r3b_bs35 = [p_r3b_bs35 m3(i, 19 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil))];
    end
end

X3a_bs35 = [length_vec3_bs35; iSite_vec3_bs35; p_occ3a_bs35; p_r3a_bs35]' ;

X3b_bs35 = [length_vec3_bs35; iSite_vec3_bs35; p_occ3b_bs35; p_r3b_bs35]';
