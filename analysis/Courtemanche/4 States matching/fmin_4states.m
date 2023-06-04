%% Data
fig3k=[ 9.554461, 7.642275, 3.000000, 2.471568];
fig3n=[25,36,47,75];
fig3p=[5,5,5,5];
fig3pro=[4,4,4,4];

fig4k=[ 9.617512, 6.885804, 3.000000, 1.486034, 4.294227, 9.554461];
fig4n=[48, 45, 47, 26, 23, 25];
fig4p=[14,  7,  5, 14,  7,  5];
fig4pro=[13,  6,  4, 13,  6,  4];

allk=[9.554461, 7.642275, 3.000000, 2.471568, 9.617512, 6.885804, 1.486034, 4.294227];
alln=[25, 36, 47, 75, 48, 45, 26, 23];
allp=[5,  5,  5,  5, 14,  7, 14,  7];
allpro=[4,  4,  4,  4, 13,  6, 13,  6];

%% Equations that do not change
kc= @(k_cap,c_PA,p_occ) k_cap*c_PA*(1-p_occ);  
kd= @(k_pab,p_occ_0,p_r) k_pab*(1-p_occ_0)*p_r;
p_r=@(n) 2.6*10^7*n.^(-3/2);

kp= @(k_cap,c_PA,p_occ,k_pab,p_occ_0,n,k_paf_rev,r_PF_rev,r_paf_rev) (1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./(kd(k_pab,p_occ_0,p_r(n)) .* r_PF_rev)) + (((k_paf_rev .* r_paf_rev) + (k_paf_rev .* r_PF_rev) + (kd(k_pab,p_occ_0,p_r(n)) .* r_PF_rev))./(kc(k_cap,c_PA,p_occ) .* kd(k_pab,p_occ_0,p_r(n)) .* r_PF_rev));

% (520000*2.5*exp(-L)*kpab*kpaf*R)/(n^(3/2)*(((2*10^6*exp(-L)*kpab*R)/n^(3/2))+exp(-L)*M*(exp(-L)*R*rpafrev)+0.1*2.5*kpaf*(((2*10^6*exp(-L)*kpab*R)/n^(3/2))+exp(-L)*R+rpafrev)))

%% Estimate inputs
p_occ=0.9;
p_occ_0=0.8;
c_PA=2.5;

%% Equations that could change

% Model 1
r_PF_rev1= @(R,L) R*exp(-L);

k_paf_rev1=@(M,L) M*exp(-L);

k_cap1=@(k_paf,L) k_paf; 

% Model 2
r_PF_rev2= @(R,L) R*exp(-L);

k_paf_rev2=@(M,L) M*exp(-L);

k_cap2=@(k_paf,L) k_paf*L;

% Model 3
r_PF_rev3= @(R,L) R*exp(-L);

k_paf_rev3=@(M,L) M; 

k_cap3=@(k_paf,L) k_paf*L;

% Model 4
r_PF_rev4= @(R,L) R;

k_paf_rev4=@(M,L) M;

k_cap4=@(k_paf,L) k_paf*L;


%% Final Kpoly eq

Kpoly1= @(n,L,k_paf,k_pab,M,R,r_paf_rev) kp(k_cap1(k_paf,L),c_PA,p_occ,k_pab,p_occ_0,n,k_paf_rev1(M,L),r_PF_rev1(R,L),r_paf_rev);
Kpoly2= @(n,L,k_paf,k_pab,M,R,r_paf_rev) kp(k_cap2(k_paf,L),c_PA,p_occ,k_pab,p_occ_0,n,k_paf_rev2(M,L),r_PF_rev2(R,L),r_paf_rev);
Kpoly3= @(n,L,k_paf,k_pab,M,R,r_paf_rev) kp(k_cap3(k_paf,L),c_PA,p_occ,k_pab,p_occ_0,n,k_paf_rev3(M,L),r_PF_rev3(R,L),r_paf_rev);
Kpoly4= @(n,L,k_paf,k_pab,M,R,r_paf_rev) kp(k_cap4(k_paf,L),c_PA,p_occ,k_pab,p_occ_0,n,k_paf_rev4(M,L),r_PF_rev4(R,L),r_paf_rev);

%% Sum of squares

SOS1= @(kpol,n,L,k_paf,k_pab,M,R,r_paf_rev) sum(abs((Kpoly1(n,L,abs(k_paf),abs(k_pab),abs(M),abs(R),abs(r_paf_rev))-kpol).^2));
SOS2= @(kpol,n,L,k_paf,k_pab,M,R,r_paf_rev) sum(abs((Kpoly2(n,L,abs(k_paf),abs(k_pab),abs(M),abs(R),abs(r_paf_rev))-kpol).^2));
SOS3= @(kpol,n,L,k_paf,k_pab,M,R,r_paf_rev) sum(abs((Kpoly3(n,L,abs(k_paf),abs(k_pab),abs(M),abs(R),abs(r_paf_rev))-kpol).^2));
SOS4= @(kpol,n,L,k_paf,k_pab,M,R,r_paf_rev) sum(abs((Kpoly4(n,L,abs(k_paf),abs(k_pab),abs(M),abs(R),abs(r_paf_rev))-kpol).^2));


%% Table
% sz=[0 8]; 
% varTypes = ["double","string", "string","double","double","double","double","double"];
% varNames = ["SOS","Model","L","k_paf","k_pab","M","R","r_paf_rev"];
% 
% T=table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

%% fminsearch

options = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolX',10^(-1),'TolFun',10^(-3)); 

[params_AA,fval_AA] = fminsearch(@(z) SOS2(allk,alln,allp,z(1),z(2),z(3),z(4),z(5)),[500;10;10^8;10^7;10^7],options);
[params_pro,fval_pro] = fminsearch(@(z) SOS1(allk,alln,allpro,z(1),z(2),z(3),z(4),z(5)),[500;10;10^8;10^7;10^7],options);

%% rand fminsearch loops
j=i;
for i=j:100000
    disp(i)
    r=15*rand(1,5)+-2;
    % Model 1
        modelnum="Model 1"; 
        
        [params_pro,fval_pro] = fminsearch(@(z) SOS1(allk,alln,allpro,z(1),z(2),z(3),z(4),z(5)),[10^r(1);10^r(2);10^r(3);10^r(4);10^r(5)],options); 
        rowData={fval_pro,modelnum,"pro",params_pro(1),params_pro(2),params_pro(3),params_pro(4),params_pro(5)};
        Tnew=[T;rowData];
        T=Tnew;
        
        [params_AA,fval_AA] = fminsearch(@(z) SOS1(allk,alln,allp,z(1),z(2),z(3),z(4),z(5)),[10^r(1);10^r(2);10^r(3);10^r(4);10^r(5)],options);
        rowData={fval_AA,modelnum,"AA",params_AA(1),params_AA(2),params_AA(3),params_AA(4),params_AA(5)};
        Tnew=[T;rowData];
        T=Tnew;
    % Model 2
        modelnum="Model 2"; 
        
        [params_pro,fval_pro] = fminsearch(@(z) SOS2(allk,alln,allpro,z(1),z(2),z(3),z(4),z(5)),[10^r(1);10^r(2);10^r(3);10^r(4);10^r(5)],options); 
        rowData={fval_pro,modelnum,"pro",params_pro(1),params_pro(2),params_pro(3),params_pro(4),params_pro(5)};
        Tnew=[T;rowData];
        T=Tnew;
        
        [params_AA,fval_AA] = fminsearch(@(z) SOS2(allk,alln,allp,z(1),z(2),z(3),z(4),z(5)),[10^r(1);10^r(2);10^r(3);10^r(4);10^r(5)],options);
        rowData={fval_AA,modelnum,"AA",params_AA(1),params_AA(2),params_AA(3),params_AA(4),params_AA(5)};
        Tnew=[T;rowData];
        T=Tnew;
    % Model 3
        modelnum="Model 3"; 
        
        [params_pro,fval_pro] = fminsearch(@(z) SOS3(allk,alln,allpro,z(1),z(2),z(3),z(4),z(5)),[10^r(1);10^r(2);10^r(3);10^r(4);10^r(5)],options); 
        rowData={fval_pro,modelnum,"pro",params_pro(1),params_pro(2),params_pro(3),params_pro(4),params_pro(5)};
        Tnew=[T;rowData];
        T=Tnew;
        
        [params_AA,fval_AA] = fminsearch(@(z) SOS3(allk,alln,allp,z(1),z(2),z(3),z(4),z(5)),[10^r(1);10^r(2);10^r(3);10^r(4);10^r(5)],options);
        rowData={fval_AA,modelnum,"AA",params_AA(1),params_AA(2),params_AA(3),params_AA(4),params_AA(5)};
        Tnew=[T;rowData];
        T=Tnew;
    % Model 4
        modelnum="Model 4"; 
        
        [params_pro,fval_pro] = fminsearch(@(z) SOS4(allk,alln,allpro,z(1),z(2),z(3),z(4),z(5)),[10^r(1);10^r(2);10^r(3);10^r(4);10^r(5)],options); 
        rowData={fval_pro,modelnum,"pro",params_pro(1),params_pro(2),params_pro(3),params_pro(4),params_pro(5)};
        Tnew=[T;rowData];
        T=Tnew;
        
        [params_AA,fval_AA] = fminsearch(@(z) SOS4(allk,alln,allp,z(1),z(2),z(3),z(4),z(5)),[10^r(1);10^r(2);10^r(3);10^r(4);10^r(5)],options);
        rowData={fval_AA,modelnum,"AA",params_AA(1),params_AA(2),params_AA(3),params_AA(4),params_AA(5)};
        Tnew=[T;rowData];
        T=Tnew;
    % Save table to workspace
        %save("fminsearch_rand.mat") 
end

%% fminsearch loop
% 
% %modelnum="Model 1"; 
% %modelnum="Model 2"; 
% %modelnum="Model 3"; 
% modelnum="Model 4"; 
% 
% 
% 
% % using pro
% for a= 50:200:850
%     disp(a)
%     for b=50:200:850
%         for c=1:1:4
%             for d=1:1:4
%                 for e=1:1:4
%                     [params_pro,fval_pro] = fminsearch(@(z) SOS(allk,alln,allpro,z(1),z(2),z(3),z(4),z(5)),[a;b;10^c;10^d;10^e],options);
%                     rowData={fval_pro,modelnum,"pro",params_pro(1),params_pro(2),params_pro(3),params_pro(4),params_pro(5)};
%                     Tnew=[T;rowData];
%                     T=Tnew; 
%                     %disp(rowData)
%                 end
%             end
%         end
%     end
% end
% 
% % using AA
% for a= 50:200:850
%     disp(a)
%     for b=50:200:850
%         for c=1:1:4
%             for d=1:1:4
%                 for e=1:1:4
%                     [params_AA,fval_AA] = fminsearch(@(z) SOS(allk,alln,allp,z(1),z(2),z(3),z(4),z(5)),[a;b;10^c;10^d;10^e],options);
%                     rowData={fval_AA,modelnum,"AA",params_AA(1),params_AA(2),params_AA(3),params_AA(4),params_AA(5)};
%                     Tnew=[T;rowData];
%                     T=Tnew; 
%                     %disp(rowData)
%                 end
%             end
%         end
%     end
% end
%                     
