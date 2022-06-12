% compares effects of different dimerization forces 
% for length = 32
% all simulated filiments have 2 strands

% calls dmr_force_explr

err = 0;
NFil = 2; N = 32;

% simulation txts saved as .mat
% n1 = dlmread('N32_2_downKSCRIT.txt');
% n2 = dlmread('dmr_F10_5.txt');
% n3 = dlmread('dmr_F10_4.txt');
% n4 = dlmread('dmr_F10_3.txt');
% n5 = dlmread('dmr_F10_2.txt');
% n6 = dlmread('dmr_F05.txt');
% n7 = dlmread('dmr_F1.txt');
% n8 = dlmread('dmr_F5.txt');
% n9 = dlmread('N32_3_downKSCRIT.txt');
% n10 = dlmread('dmr_F15.txt');

%save('dmr_force_explr.mat', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'n10');
load('dmr_force_explr.mat');

%n-dimerization force used in each simulation
forces = [n1(6) n2(6) n3(6) n4(6) n5(6) n6(6) n7(6) n8(6) n9(6) n10(6)];

% extract all probabilities of occlusion

% the following contain as rows all probabilities of occlusion 
p_occ2a_all = [];
p_occ2b_all = [];
 
p_occ2a = []; %temporary variables
p_occ2b = [];
for iN=32
    N_All =NFil*iN;
    for iy =1:iN
        p_a = n1(16 + 2*(N_All +1) + 7*(iy - 1));
        p_b = n1(16 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil));
        if p_a > 1 || p_b >1
            err = err+1;
        end
        p_occ2a = [p_occ2a p_a];
        p_occ2b = [p_occ2b p_b];
    end
end

p_occ2a_all = [p_occ2a_all; p_occ2a];
p_occ2b_all = [p_occ2b_all; p_occ2b]; 

p_occ2a = [];
p_occ2b = [];
for iN=32
    N_All =NFil*iN;
    for iy =1:iN
        p_a = n2(16 + 2*(N_All +1) + 7*(iy - 1));
        p_b = n2(16 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil));
        if p_a > 1 || p_b >1
            err = err+1;
        end
        p_occ2a = [p_occ2a p_a];
        p_occ2b = [p_occ2b p_b];
    end
end

p_occ2a_all = [p_occ2a_all; p_occ2a];
p_occ2b_all = [p_occ2b_all; p_occ2b]; 

p_occ2a = [];
p_occ2b = [];
for iN=32
    N_All =NFil*iN;
    for iy =1:iN
        p_a = n3(16 + 2*(N_All +1) + 7*(iy - 1));
        p_b = n3(16 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil));
        if p_a > 1 || p_b >1
            err = err+1;
        end
        p_occ2a = [p_occ2a p_a];
        p_occ2b = [p_occ2b p_b];
    end
end

p_occ2a_all = [p_occ2a_all; p_occ2a];
p_occ2b_all = [p_occ2b_all; p_occ2b]; 

p_occ2a = [];
p_occ2b = [];
for iN=32
    N_All =NFil*iN;
    for iy =1:iN
        p_a = n4(16 + 2*(N_All +1) + 7*(iy - 1));
        p_b = n4(16 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil));
        if p_a > 1 || p_b >1
            err = err+1;
        end
        p_occ2a = [p_occ2a p_a];
        p_occ2b = [p_occ2b p_b];
    end
end

p_occ2a_all = [p_occ2a_all; p_occ2a];
p_occ2b_all = [p_occ2b_all; p_occ2b]; 


p_occ2a = [];
p_occ2b = [];
for iN=32
    N_All =NFil*iN;
    for iy =1:iN
        p_a = n5(16 + 2*(N_All +1) + 7*(iy - 1));
        p_b = n5(16 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil));
        if p_a > 1 || p_b >1
            err = err+1;
        end
        p_occ2a = [p_occ2a p_a];
        p_occ2b = [p_occ2b p_b];
    end
end

p_occ2a_all = [p_occ2a_all; p_occ2a];
p_occ2b_all = [p_occ2b_all; p_occ2b]; 


p_occ2a = [];
p_occ2b = [];
for iN=32
    N_All =NFil*iN;
    for iy =1:iN
        p_a = n6(16 + 2*(N_All +1) + 7*(iy - 1));
        p_b = n6(16 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil));
        if p_a > 1 || p_b >1
            err = err+1;
        end
        p_occ2a = [p_occ2a p_a];
        p_occ2b = [p_occ2b p_b];
    end
end

p_occ2a_all = [p_occ2a_all; p_occ2a];
p_occ2b_all = [p_occ2b_all; p_occ2b]; 

p_occ2a = [];
p_occ2b = [];
for iN=32
    N_All =NFil*iN;
    for iy =1:iN
        p_a = n7(16 + 2*(N_All +1) + 7*(iy - 1));
        p_b = n7(16 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil));
        if p_a > 1 || p_b >1
            err = err+1;
        end
        p_occ2a = [p_occ2a p_a];
        p_occ2b = [p_occ2b p_b];
    end
end

p_occ2a_all = [p_occ2a_all; p_occ2a];
p_occ2b_all = [p_occ2b_all; p_occ2b]; 

p_occ2a = [];
p_occ2b = [];
for iN=32
    N_All =NFil*iN;
    for iy =1:iN
        p_a = n8(16 + 2*(N_All +1) + 7*(iy - 1));
        p_b = n8(16 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil));
        if p_a > 1 || p_b >1
            err = err+1;
        end
        p_occ2a = [p_occ2a p_a];
        p_occ2b = [p_occ2b p_b];
    end
end

p_occ2a_all = [p_occ2a_all; p_occ2a];
p_occ2b_all = [p_occ2b_all; p_occ2b]; 



p_occ2a = [];
p_occ2b = [];
for iN=32
    N_All =NFil*iN;
    for iy =1:iN
        p_a = n9(16 + 2*(N_All +1) + 7*(iy - 1));
        p_b = n9(16 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil));
        if p_a > 1 || p_b >1
            err = err+1;
        end
        p_occ2a = [p_occ2a p_a];
        p_occ2b = [p_occ2b p_b];
    end
end

p_occ2a_all = [p_occ2a_all; p_occ2a];
p_occ2b_all = [p_occ2b_all; p_occ2b]; 


p_occ2a = [];
p_occ2b = [];
for iN=32
    N_All =NFil*iN;
    for iy =1:iN
        p_a = n10(16 + 2*(N_All +1) + 7*(iy - 1));
        p_b = n10(16 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil));
        if p_a > 1 || p_b >1
            err = err+1;
        end
        p_occ2a = [p_occ2a p_a];
        p_occ2b = [p_occ2b p_b];
    end
end

p_occ2a_all = [p_occ2a_all; p_occ2a];
p_occ2b_all = [p_occ2b_all; p_occ2b]; 


%x axis: iSites

iSites = 1:32;

% Set color vector
c = parula(10);

% Visualize the result
figure()
hold on
for kk = 1:10
  plot(iSites,p_occ2a_all(kk,:),'Color',c(kk,:))
end
forces_name = string(forces(1:10));
xlabel('iSite_location');
ylabel('probability of occlusion');

hleg = legend(forces_name);
htitle = get(hleg,'Title');
set(htitle,'String','Dimer Force')


hold off
