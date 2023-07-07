k_paf=@(k_paf) k_paf;

c_PA=2.5; % μM %concentration of profilin-actin

k_pab=@(k_pab) k_pab;

k_paf_rev=@(k_paf_rev) k_paf_rev;

r_PF_rev=@(r_PF_rev) r_PF_rev;

r_paf_rev=@(r_paf_rev) r_paf_rev;

kdb = @(o,pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) sum(1./((1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./((k_pab.*(1-o.p_occ2a_0(:)).*(1.0e33.*o.p_r2a(:)/(27*6.022e23))) .* r_PF_rev)) + ((((k_paf_rev.*exp(-1.*pp_length_vec)) .* r_paf_rev) + ((k_paf_rev.*exp(-1.*pp_length_vec)) .* r_PF_rev) + ((k_pab.*(1-o.p_occ2a_0(:)).*(1.0e33.*o.p_r2a(:)/(27*6.022e23))).* r_PF_rev))./((k_paf.*c_PA.*(1-o.p_occ2a(:))).* (k_pab.*(1-o.p_occ2a_0(:)).*(1.0e33.*o.p_r2a(:)/(27*6.022e23))) .* r_PF_rev))))...
    + sum(1./((1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./((k_pab.*(1-o.p_occ2b_0(:)).*(1.0e33.*o.p_r2b(:)/(27*6.022e23))) .* r_PF_rev)) + ((((k_paf_rev.*exp(-1.*pp_length_vec)) .* r_paf_rev) + ((k_paf_rev.*exp(-1.*pp_length_vec)) .* r_PF_rev) + ((k_pab.*(1-o.p_occ2b_0(:)).*(1.0e33.*o.p_r2b(:)/(27*6.022e23))).* r_PF_rev))./((k_paf.*c_PA.*(1-o.p_occ2b(:))).* (k_pab.*(1-o.p_occ2b_0(:)).*(1.0e33.*o.p_r2b(:)/(27*6.022e23))) .* r_PF_rev))));

kdim = @(o,pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) sum(1./((1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./((k_pab.*(1-o.p_occ3a_0).*(1.0e33.*o.p_r3a/(27*6.022e23))) .* r_PF_rev)) + ((((k_paf_rev.*exp(-1.*pp_length_vec)) .* r_paf_rev) + ((k_paf_rev.*exp(-1.*pp_length_vec)) .* r_PF_rev) + ((k_pab.*(1-o.p_occ3a_0).*(1.0e33.*o.p_r3a/(27*6.022e23))).* r_PF_rev))./((k_paf.*c_PA.*(1-o.p_occ3a)).* (k_pab.*(1-o.p_occ3a_0).*(1.0e33.*o.p_r3a/(27*6.022e23))) .* r_PF_rev))))...
    + sum(1./((1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./((k_pab.*(1-o.p_occ3b_0).*(1.0e33.*o.p_r3b/(27*6.022e23))) .* r_PF_rev)) + ((((k_paf_rev.*exp(-1.*pp_length_vec)) .* r_paf_rev) + ((k_paf_rev.*exp(-1.*pp_length_vec)) .* r_PF_rev) + ((k_pab.*(1-o.p_occ3b_0).*(1.0e33.*o.p_r3b/(27*6.022e23))).* r_PF_rev))./((k_paf.*c_PA.*(1-o.p_occ3b)).* (k_pab.*(1-o.p_occ3b_0).*(1.0e33.*o.p_r3b/(27*6.022e23))) .* r_PF_rev))));

%%
load('lookup_ltvar.mat');

ratios=[-0.6,0];

FHOD.name='FHODB';
FHOD.pp_length_vec=[6; 15; 4];
FHOD.pp_index_vec=[47,218,239];
FHOD.fh1_length=255;
FHOD.o=polymerstats(FHOD,lt);
FHOD.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(FHOD.o,FHOD.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD.kdim=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim(FHOD.o,FHOD.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

Capu.name='Capu';
Capu.pp_length_vec=[6;8;9;14;9];
Capu.pp_index_vec=[33,47,61,78,96];
Capu.fh1_length=123;
Capu.o=polymerstats(Capu,lt);
Capu.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(Capu.o,Capu.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
Capu.kdim=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim(Capu.o,Capu.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

%%
Dlength=5;
Blength=14;
Clength=7;

D18.name='BNI1P(pPD18)';
D18.pp_length_vec=Dlength;
D18.pp_index_vec=25;
D18.fh1_length=88;
D18.o=polymerstats(D18,lt);
D18.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(D18.o,D18.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

D28.name='BNI1P(pPD28)';
D28.pp_length_vec=Dlength;
D28.pp_index_vec=36;
D28.fh1_length=88;
D28.o=polymerstats(D28,lt);
D28.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(D28.o,D28.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

D37.name='BNI1P(pPD37)';
D37.pp_length_vec=Dlength;
D37.pp_index_vec=47;
D37.fh1_length=88;
D37.o=polymerstats(D37,lt);
D37.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(D37.o,D37.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

D65.name='BNI1P(pPD65)';
D65.pp_length_vec=Dlength;
D65.pp_index_vec=75;
D65.fh1_length=88;
D65.o=polymerstats(D65,lt);
D65.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(D65.o,D65.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

B37.name='BNI1P(pPB37)';
B37.pp_length_vec=Blength;
B37.pp_index_vec=48;
B37.fh1_length=94;
B37.o=polymerstats(B37,lt);
B37.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(B37.o,B37.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

C37.name='BNI1P(pPC37)';
C37.pp_length_vec=Clength;
C37.pp_index_vec=45;
C37.fh1_length=87;
C37.o=polymerstats(C37,lt);
C37.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(C37.o,C37.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

B18.name='BNI1P(pPB18)';
B18.pp_length_vec=Blength;
B18.pp_index_vec=26;
B18.fh1_length=94;
B18.o=polymerstats(B18,lt);
B18.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(B18.o,B18.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

C18.name='BNI1P(pPC18)';
C18.pp_length_vec=Clength;
C18.pp_index_vec=23;
C18.fh1_length=87;
C18.o=polymerstats(C18,lt);
C18.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(C18.o,C18.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

%% NTD eq
SOSNTD= @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) 10*(sum(abs(([...
        (log2(FHOD.kdim(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(FHOD.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)))),...
        (log2(Capu.kdim(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(Capu.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))))...
    ]-ratios).^2)));

SOSFHOD=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) (sum(abs(([...
        (log2(FHOD.kdim(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(FHOD.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))))...
    ]-ratios(1)).^2)));

SODFHOD=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) (((([...
        (log2(FHOD.kdim(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(FHOD.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))))...
    ]-ratios(1)))));

SOSCapu=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) (sum(abs(([...
        (log2(Capu.kdim(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(Capu.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))))...
    ]-ratios(2)).^2)));

%[params_NTD,fval_NTD] = fminsearch(@(z) SOSNTD(abs(z(1)),abs(z(2)),abs(z(3)),abs(z(4)),abs(z(5))),[457,10,32704292,1250000,50000],options);
%vals=num2cell(abs(params_NTD));

% for i=1:1000
%     disp(i)
%     r=15*rand(1,5)+-3;
%         [params,fval] = fminsearch(@(z) SOSFHOD(abs(z(1)),abs(z(2)),abs(z(3)),abs(z(4)),abs(z(5))),[10^r(1);10^r(2);10^r(3);10^r(4);10^r(5)],options); 
%         vals=num2cell(abs(params));
%         FHOD_T=[FHOD_T;[{fval},vals(:)']];
% end

for i=1:1000
    disp(i)
    r=15*rand(1,5)+-3;
        [params,fval] = fminsearch(@(z) SOSCapu(abs(z(1)),abs(z(2)),abs(z(3)),abs(z(4)),abs(z(5))),[10^r(1);10^r(2);10^r(3);10^r(4);10^r(5)],options); 
        vals=num2cell(abs(params));
        Capu_T=[Capu_T;[{fval},vals(:)']];
end
%% Fig 3 eq
fig3k=[ 9.554461, 7.642275, 3.000000, 2.471568];
fig3sim=[D18,D28,D37,D65];
fig3ksim={fig3sim.kdob};
fig3ksimf = @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) cellfun(@(f) f(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev), fig3ksim(:));

fig3ksimf_r = @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) fig3ksimf(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(D37.kdob(1,1,1,1,1)/3);

SOS3= @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) sum(abs(([...
        fig3ksimf(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ]'-fig3k).^2));

SOS3_r= @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) sum(abs(([...
        fig3ksimf_r(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ]'-fig3k).^2));

%% Fig 4 eq
fig4ak=[ 9.617512, 6.885804, 3.000000];
fig4asim=[B37,C37,D37];
fig4aksim={fig4asim.kdob};
fig4aksimf = @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) cellfun(@(f) f(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev), fig4aksim(:));

fig4aksimf_r = @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) fig4aksimf(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(D37.kdob(1,1,1,1,1)/3);


SOS4a= @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) sum(abs(([...
        fig4aksimf(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ]'-fig4ak).^2));

SOS4a_r= @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) sum(abs(([...
        fig4aksimf_r(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ]'-fig4ak).^2));

%% Full sos
%SOSs={SOS4a,SOS3,SOSNTD};
SOS=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) SOSNTD(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)+...
    SOS3(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)+...
    SOS4a(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

SOS_r=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) SOSNTD(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)+...
    SOS3_r(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)+...
    SOS4a_r(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolX',10^(-10),'TolFun',10^(-10)); 
%options = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolX',10^(-15),'TolFun',10^(-15),'PlotFcns','optimplotfval','Display','iter','TypicalX',[25,10^(1),10^(15),10^(10),10^(10)]); 
%options=optimoptions('fsolve', 'Display','iter','ScaleProblem','none','InitDamping',0.001,'FunctionTolerance',10^(-9),'StepTolerance',10^(-10));

%[params,fval] = fminsearch(@(z) SOS(abs(z(1)),abs(z(2)),abs(z(3)),abs(z(4)),abs(z(5))),[25,1,32704292,1250000,50000],options);
%[params,fval] = fminsearch(@(z) SOS_r(abs(z(1)),abs(z(2)),abs(z(3)),abs(z(4)),abs(z(5))),[457,10,32704292,1250000,50000],options);
[params_FHOD,fval_FHOD] = fminsearch(@(z) SOSFHOD(abs(z(1)),abs(z(2)),abs(z(3)),abs(z(4)),abs(z(5))),[25,1,32704292,1250000,50000],options);
[params,fval]=fminsearch(@(z) SOSFHOD(abs(z(1)),abs(z(2)),abs(z(3)),abs(z(4)),abs(z(5))),[25,1,32704292,1250000,50000],options)


for i=1:1000
    disp(i)
    r=15*rand(1,5)+-3;
    % Not ratio
        [params,fval] = fminsearch(@(z) SOS(abs(z(1)),abs(z(2)),abs(z(3)),abs(z(4)),abs(z(5))),[10^r(1);10^r(2);10^r(3);10^r(4);10^r(5)],options); 
        vals=num2cell(abs(params));
        T=[T;[{fval},vals(:)']];
    % ratio
        [params,fval] = fminsearch(@(z) SOS_r(abs(z(1)),abs(z(2)),abs(z(3)),abs(z(4)),abs(z(5))),[10^r(1);10^r(2);10^r(3);10^r(4);10^r(5)],options); 
        vals=num2cell(abs(params));
        T=[T;[{fval},vals(:)']];
end

%% NTD Bar
%x=[log2(FHOD.kdim/FHOD.kdob); log2(Capu.kdim/Capu.kdob)];
figure
x=[...
        (log2(FHOD.kdim(vals{:})./(FHOD.kdob(vals{:})))),...
        (log2(Capu.kdim(vals{:})./(Capu.kdob(vals{:}))))...
    ];
y={FHOD.name; Capu.name};
kpoly_table_ratio = table(x', 'RowNames', y);
kpoly_bar_ratio = bar(kpoly_table_ratio{:,:});
set(gca,'xtick',1:2, 'xticklabel',kpoly_table_ratio.Properties.RowNames)
xtickangle(90)
set(kpoly_bar_ratio(1), 'FaceColor','m')
xlabel('Formins')
ylabel('log_2(kpoly N terminal dimerized/kpoly double)')
ylim([-0.62 0.3])

title('Change in Polymerization Rates with NTD per Formin')
%% Fig 3
figure
kpolybar=bar([fig3k',fig3ksimf(vals{:})]);
set(gca,'xticklabel',{fig3sim.name});

    set(kpolybar(1), 'FaceColor','b');
    set(kpolybar(2), 'FaceColor','r');
    legend( 'Experimental', 'model');
    xlabel('Formins');
    ylabel('kpoly');
    title("Fig 3")

%% Fig 4
figure
kpolybar=bar([fig4ak',fig4aksimf(vals{:})]);
set(gca,'xticklabel',{fig4asim.name});

set(kpolybar(1), 'FaceColor','b');
set(kpolybar(2), 'FaceColor','r');
legend( 'Experimental', 'model');
xlabel('Formins');
ylabel('kpoly');
title("Fig 4a")