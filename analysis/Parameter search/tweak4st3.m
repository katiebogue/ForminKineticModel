k_paf=@(k_paf) k_paf;

c_PA=2.5; % μM %concentration of profilin-actin

k_pab=@(k_pab) k_pab;

k_paf_rev=@(k_paf_rev) k_paf_rev;

r_PF_rev=@(r_PF_rev) r_PF_rev;

r_paf_rev=@(r_paf_rev) r_paf_rev;

pr_db=@(pp_index_vec) (2*pi.*pp_index_vec./3).^(-3/2); % where b=1 (since we adjust into AA in kdb)

pr_dim=@(pp_index_vec,fh1_length) (2*pi.*pp_index_vec.*(2*fh1_length-pp_index_vec)./(3*2*fh1_length)).^(-3/2); % where b=1 (since we adjust into AA in kdim)

kdb = @(o,pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) sum(1./((1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./((k_pab.*(1-o.p_occ2a_0(:)).*(1.0e33.*o.p_r2a(:)/(27*6.022e23))) .* r_PF_rev)) + ((((k_paf_rev.*exp(-1.*pp_length_vec)) .* r_paf_rev) + ((k_paf_rev.*exp(-1.*pp_length_vec)) .* r_PF_rev) + ((k_pab.*(1-o.p_occ2a_0(:)).*(1.0e33.*o.p_r2a(:)/(27*6.022e23))).* r_PF_rev))./((k_paf.*c_PA.*(1-o.p_occ2a(:))).* (k_pab.*(1-o.p_occ2a_0(:)).*(1.0e33.*o.p_r2a(:)/(27*6.022e23))) .* r_PF_rev))))...
    + sum(1./((1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./((k_pab.*(1-o.p_occ2b_0(:)).*(1.0e33.*o.p_r2b(:)/(27*6.022e23))) .* r_PF_rev)) + ((((k_paf_rev.*exp(-1.*pp_length_vec)) .* r_paf_rev) + ((k_paf_rev.*exp(-1.*pp_length_vec)) .* r_PF_rev) + ((k_pab.*(1-o.p_occ2b_0(:)).*(1.0e33.*o.p_r2b(:)/(27*6.022e23))).* r_PF_rev))./((k_paf.*c_PA.*(1-o.p_occ2b(:))).* (k_pab.*(1-o.p_occ2b_0(:)).*(1.0e33.*o.p_r2b(:)/(27*6.022e23))) .* r_PF_rev))));

kdim = @(o,pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) sum(1./((1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./((k_pab.*(1-o.p_occ3a_0(:)).*(1.0e33.*o.p_r3a(:)/(27*6.022e23))) .* r_PF_rev)) + ((((k_paf_rev.*exp(-1.*pp_length_vec)) .* r_paf_rev) + ((k_paf_rev.*exp(-1.*pp_length_vec)) .* r_PF_rev) + ((k_pab.*(1-o.p_occ3a_0(:)).*(1.0e33.*o.p_r3a(:)/(27*6.022e23))).* r_PF_rev))./((k_paf.*c_PA.*(1-o.p_occ3a(:))).* (k_pab.*(1-o.p_occ3a_0(:)).*(1.0e33.*o.p_r3a(:)/(27*6.022e23))) .* r_PF_rev))))...
    + sum(1./((1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./((k_pab.*(1-o.p_occ3b_0(:)).*(1.0e33.*o.p_r3b(:)/(27*6.022e23))) .* r_PF_rev)) + ((((k_paf_rev.*exp(-1.*pp_length_vec)) .* r_paf_rev) + ((k_paf_rev.*exp(-1.*pp_length_vec)) .* r_PF_rev) + ((k_pab.*(1-o.p_occ3b_0(:)).*(1.0e33.*o.p_r3b(:)/(27*6.022e23))).* r_PF_rev))./((k_paf.*c_PA.*(1-o.p_occ3b(:))).* (k_pab.*(1-o.p_occ3b_0(:)).*(1.0e33.*o.p_r3b(:)/(27*6.022e23))) .* r_PF_rev))));

kdb_calc = @(o,pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) sum(1./((1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./((k_pab.*(1-o.p_occ2a_0(:)).*(1.0e33.*o.pr_2_calc(:)/(27*6.022e23))) .* r_PF_rev)) + ((((k_paf_rev.*exp(-1.*pp_length_vec)) .* r_paf_rev) + ((k_paf_rev.*exp(-1.*pp_length_vec)) .* r_PF_rev) + ((k_pab.*(1-o.p_occ2a_0(:)).*(1.0e33.*o.pr_2_calc(:)/(27*6.022e23))).* r_PF_rev))./((k_paf.*c_PA.*(1-o.p_occ2a(:))).* (k_pab.*(1-o.p_occ2a_0(:)).*(1.0e33.*o.pr_2_calc(:)/(27*6.022e23))) .* r_PF_rev))))...
    + sum(1./((1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./((k_pab.*(1-o.p_occ2b_0(:)).*(1.0e33.*o.pr_2_calc(:)/(27*6.022e23))) .* r_PF_rev)) + ((((k_paf_rev.*exp(-1.*pp_length_vec)) .* r_paf_rev) + ((k_paf_rev.*exp(-1.*pp_length_vec)) .* r_PF_rev) + ((k_pab.*(1-o.p_occ2b_0(:)).*(1.0e33.*o.pr_2_calc(:)/(27*6.022e23))).* r_PF_rev))./((k_paf.*c_PA.*(1-o.p_occ2b(:))).* (k_pab.*(1-o.p_occ2b_0(:)).*(1.0e33.*o.pr_2_calc(:)/(27*6.022e23))) .* r_PF_rev))));

kdim_calc = @(o,pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) sum(1./((1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./((k_pab.*(1-o.p_occ3a_0(:)).*(1.0e33.*o.pr_3_calc(:)/(27*6.022e23))) .* r_PF_rev)) + ((((k_paf_rev.*exp(-1.*pp_length_vec)) .* r_paf_rev) + ((k_paf_rev.*exp(-1.*pp_length_vec)) .* r_PF_rev) + ((k_pab.*(1-o.p_occ3a_0(:)).*(1.0e33.*o.pr_3_calc(:)/(27*6.022e23))).* r_PF_rev))./((k_paf.*c_PA.*(1-o.p_occ3a(:))).* (k_pab.*(1-o.p_occ3a_0(:)).*(1.0e33.*o.pr_3_calc(:)/(27*6.022e23))) .* r_PF_rev))))...
    + sum(1./((1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./((k_pab.*(1-o.p_occ3b_0(:)).*(1.0e33.*o.pr_3_calc(:)/(27*6.022e23))) .* r_PF_rev)) + ((((k_paf_rev.*exp(-1.*pp_length_vec)) .* r_paf_rev) + ((k_paf_rev.*exp(-1.*pp_length_vec)) .* r_PF_rev) + ((k_pab.*(1-o.p_occ3b_0(:)).*(1.0e33.*o.pr_3_calc(:)/(27*6.022e23))).* r_PF_rev))./((k_paf.*c_PA.*(1-o.p_occ3b(:))).* (k_pab.*(1-o.p_occ3b_0(:)).*(1.0e33.*o.pr_3_calc(:)/(27*6.022e23))) .* r_PF_rev))));

kdb2_calc = @(o,pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) sum(1./((1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./((k_pab.*(1-o.p_occ2a_0(:)).*(1.0e33.*o.pr_2_calc(:)/(27*6.022e23))) .* r_PF_rev)) + ((((k_paf_rev) .* r_paf_rev) + ((k_paf_rev) .* r_PF_rev) + ((k_pab.*(1-o.p_occ2a_0(:)).*(1.0e33.*o.pr_2_calc(:)/(27*6.022e23))).* r_PF_rev))./((k_paf.*c_PA.*(1-o.p_occ2a(:))).* (k_pab.*(1-o.p_occ2a_0(:)).*(1.0e33.*o.pr_2_calc(:)/(27*6.022e23))) .* r_PF_rev))))...
    + sum(1./((1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./((k_pab.*(1-o.p_occ2b_0(:)).*(1.0e33.*o.pr_2_calc(:)/(27*6.022e23))) .* r_PF_rev)) + ((((k_paf_rev) .* r_paf_rev) + ((k_paf_rev) .* r_PF_rev) + ((k_pab.*(1-o.p_occ2b_0(:)).*(1.0e33.*o.pr_2_calc(:)/(27*6.022e23))).* r_PF_rev))./((k_paf.*c_PA.*(1-o.p_occ2b(:))).* (k_pab.*(1-o.p_occ2b_0(:)).*(1.0e33.*o.pr_2_calc(:)/(27*6.022e23))) .* r_PF_rev))));

kdim2_calc = @(o,pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) sum(1./((1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./((k_pab.*(1-o.p_occ3a_0(:)).*(1.0e33.*o.pr_3_calc(:)/(27*6.022e23))) .* r_PF_rev)) + ((((k_paf_rev) .* r_paf_rev) + ((k_paf_rev) .* r_PF_rev) + ((k_pab.*(1-o.p_occ3a_0(:)).*(1.0e33.*o.pr_3_calc(:)/(27*6.022e23))).* r_PF_rev))./((k_paf.*c_PA.*(1-o.p_occ3a(:))).* (k_pab.*(1-o.p_occ3a_0(:)).*(1.0e33.*o.pr_3_calc(:)/(27*6.022e23))) .* r_PF_rev))))...
    + sum(1./((1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./((k_pab.*(1-o.p_occ3b_0(:)).*(1.0e33.*o.pr_3_calc(:)/(27*6.022e23))) .* r_PF_rev)) + ((((k_paf_rev) .* r_paf_rev) + ((k_paf_rev) .* r_PF_rev) + ((k_pab.*(1-o.p_occ3b_0(:)).*(1.0e33.*o.pr_3_calc(:)/(27*6.022e23))).* r_PF_rev))./((k_paf.*c_PA.*(1-o.p_occ3b(:))).* (k_pab.*(1-o.p_occ3b_0(:)).*(1.0e33.*o.pr_3_calc(:)/(27*6.022e23))) .* r_PF_rev))));

kdb3_calc = @(o,pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) sum(1./((1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./((k_pab.*(1-o.p_occ2a_0(:)).*(1.0e33.*o.pr_2_calc(:)/(27*6.022e23))) .* r_PF_rev)) + ((((k_paf_rev) .* r_paf_rev) + ((k_paf_rev) .* r_PF_rev) + ((k_pab.*(1-o.p_occ2a_0(:)).*(1.0e33.*o.pr_2_calc(:)/(27*6.022e23))).* r_PF_rev))./((k_paf.*pp_length_vec.*c_PA.*(1-o.p_occ2a(:))).* (k_pab.*(1-o.p_occ2a_0(:)).*(1.0e33.*o.pr_2_calc(:)/(27*6.022e23))) .* r_PF_rev))))...
    + sum(1./((1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./((k_pab.*(1-o.p_occ2b_0(:)).*(1.0e33.*o.pr_2_calc(:)/(27*6.022e23))) .* r_PF_rev)) + ((((k_paf_rev) .* r_paf_rev) + ((k_paf_rev) .* r_PF_rev) + ((k_pab.*(1-o.p_occ2b_0(:)).*(1.0e33.*o.pr_2_calc(:)/(27*6.022e23))).* r_PF_rev))./((k_paf.*pp_length_vec.*c_PA.*(1-o.p_occ2b(:))).* (k_pab.*(1-o.p_occ2b_0(:)).*(1.0e33.*o.pr_2_calc(:)/(27*6.022e23))) .* r_PF_rev))));

kdim3_calc = @(o,pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) sum(1./((1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./((k_pab.*(1-o.p_occ3a_0(:)).*(1.0e33.*o.pr_3_calc(:)/(27*6.022e23))) .* r_PF_rev)) + ((((k_paf_rev) .* r_paf_rev) + ((k_paf_rev) .* r_PF_rev) + ((k_pab.*(1-o.p_occ3a_0(:)).*(1.0e33.*o.pr_3_calc(:)/(27*6.022e23))).* r_PF_rev))./((k_paf.*pp_length_vec.*c_PA.*(1-o.p_occ3a(:))).* (k_pab.*(1-o.p_occ3a_0(:)).*(1.0e33.*o.pr_3_calc(:)/(27*6.022e23))) .* r_PF_rev))))...
    + sum(1./((1./r_PF_rev) + ((r_paf_rev + r_PF_rev)./((k_pab.*(1-o.p_occ3b_0(:)).*(1.0e33.*o.pr_3_calc(:)/(27*6.022e23))) .* r_PF_rev)) + ((((k_paf_rev) .* r_paf_rev) + ((k_paf_rev) .* r_PF_rev) + ((k_pab.*(1-o.p_occ3b_0(:)).*(1.0e33.*o.pr_3_calc(:)/(27*6.022e23))).* r_PF_rev))./((k_paf.*pp_length_vec.*c_PA.*(1-o.p_occ3b(:))).* (k_pab.*(1-o.p_occ3b_0(:)).*(1.0e33.*o.pr_3_calc(:)/(27*6.022e23))) .* r_PF_rev))));

% the gating factors scale k_pab (delivery)
Capu_gating=1; % from (Zweifel & Courtemanche, 2020)
FHOD_gating=0.2; % cannot find specifically for fhodB but FHOD is listed as 0.2 in (Zweifel & Courtemanche, 2020)
bni1_gating=0.5; % from (Zweifel & Courtemanche, 2020)
%%
load('lookup_ltvar.mat');

%ratios=[-0.6,0];
ratios=[-0.6172640408,0];
ratios_err=[0.2584555583 0.2584555583]; % +/- this value

FHOD.name='FHODB';
FHOD.pp_length_vec=[6; 15; 4];
FHOD.p_length_vec=[6; 15; 3];
FHOD.pp_length_vec_max=[6; 10; 4]; %max at 10
FHOD.pp_index_vec=[47,218,239];
FHOD.fh1_length=255;
FHOD.o=polymerstats(FHOD,lt);
FHOD.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(FHOD.o,FHOD.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD.kdim=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim(FHOD.o,FHOD.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(FHOD.o,FHOD.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD.kdimc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim_calc(FHOD.o,FHOD.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

FHODbs20.name='FHODB';
FHODbs20.pp_length_vec=[6; 15; 4];
FHODbs20.p_length_vec=[6; 15; 3];
FHODbs20.pp_index_vec=[47,218,239];
FHODbs20.fh1_length=250;
FHODbs20.o=polymerstats(FHODbs20,lt_bs20);
FHODbs20.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(FHODbs20.o,FHODbs20.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHODbs20.kdim=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim(FHODbs20.o,FHODbs20.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHODbs20.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(FHODbs20.o,FHODbs20.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHODbs20.kdimc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim_calc(FHODbs20.o,FHODbs20.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

FHODbs35.name='FHODB';
FHODbs35.pp_length_vec=[6; 15; 4];
FHODbs35.p_length_vec=[6; 15; 3];
FHODbs35.pp_index_vec=[47,218,239];
FHODbs35.fh1_length=255;
FHODbs35.o=polymerstats(FHODbs35,lt_bs35);
FHODbs35.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(FHODbs35.o,FHODbs35.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHODbs35.kdim=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim(FHODbs35.o,FHODbs35.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHODbs35.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(FHODbs35.o,FHODbs35.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHODbs35.kdimc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim_calc(FHODbs35.o,FHODbs35.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);


FHOD2.name='FHODB'; % allowing 2 interuptions
FHOD2.pp_length_vec=[6; 15; 6];
FHOD2.pp_index_vec=[47,218,238];
FHOD2.fh1_length=255;
FHOD2.o=polymerstats(FHOD2,lt);
FHOD2.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(FHOD2.o,FHOD2.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD2.kdim=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim(FHOD2.o,FHOD2.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD2.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(FHOD2.o,FHOD2.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD2.kdimc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim_calc(FHOD2.o,FHOD2.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);


FHOD3.name='FHODB'; % allowing 2 interuptions and broke 15 into 2
FHOD3.pp_length_vec=[6; 7; 8; 6];
FHOD3.pp_index_vec=[47,214,221,238];
FHOD3.fh1_length=255;
FHOD3.o=polymerstats(FHOD3,lt);
FHOD3.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(FHOD3.o,FHOD3.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD3.kdim=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim(FHOD3.o,FHOD3.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD3.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(FHOD3.o,FHOD3.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD3.kdimc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim_calc(FHOD3.o,FHOD3.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);



FHOD4.name='FHODB'; % broke 15 into 2
FHOD4.pp_length_vec=[6; 7; 8; 4];
FHOD4.pp_index_vec=[47,214,221,239];
FHOD4.fh1_length=255;
FHOD4.o=polymerstats(FHOD4,lt);
FHOD4.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(FHOD4.o,FHOD4.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD4.kdim=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim(FHOD4.o,FHOD4.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD4.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(FHOD4.o,FHOD4.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD4.kdimc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim_calc(FHOD4.o,FHOD4.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);


FHOD5.name='FHODB'; %removed last PRM
FHOD5.pp_length_vec=[6; 15];
FHOD5.pp_index_vec=[47,218];
FHOD5.fh1_length=255;
FHOD5.o=polymerstats(FHOD5,lt);
FHOD5.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(FHOD5.o,FHOD5.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD5.kdim=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim(FHOD5.o,FHOD5.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD5.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(FHOD5.o,FHOD5.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD5.kdimc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim_calc(FHOD5.o,FHOD5.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);


FHOD6.name='FHODB'; % removed first PRM
FHOD6.pp_length_vec=[15; 4];
FHOD6.pp_index_vec=[218,239];
FHOD6.fh1_length=255;
FHOD6.o=polymerstats(FHOD6,lt);
FHOD6.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(FHOD6.o,FHOD6.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD6.kdim=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim(FHOD6.o,FHOD6.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD6.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(FHOD6.o,FHOD6.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD6.kdimc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim_calc(FHOD6.o,FHOD6.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);


FHOD7.name='FHODB'; % only last PRM
FHOD7.pp_length_vec=[4];
FHOD7.pp_index_vec=[239];
FHOD7.fh1_length=255;
FHOD7.o=polymerstats(FHOD7,lt);
FHOD7.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(FHOD7.o,FHOD7.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD7.kdim=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim(FHOD7.o,FHOD7.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD7.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(FHOD7.o,FHOD7.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD7.kdimc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim_calc(FHOD7.o,FHOD7.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);


Capu.name='Capu';
Capu.pp_length_vec=[6;8;9;14;9];
Capu.p_length_vec=[6;7;9;14;9];
Capu.pp_length_vec_max=[6;8;9;10;9];
Capu.pp_index_vec=[33,47,61,78,96];
Capu.fh1_length=123;
Capu.o=polymerstats(Capu,lt);
Capu.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(Capu.o,Capu.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
Capu.kdim=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim(Capu.o,Capu.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
Capu.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(Capu.o,Capu.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
Capu.kdimc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim_calc(Capu.o,Capu.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

Capubs20.name='Capu';
Capubs20.pp_length_vec=[6;8;9;14;9];
Capubs20.p_length_vec=[6;7;9;14;9];
Capubs20.pp_index_vec=[33,47,61,78,96];
Capubs20.fh1_length=120;
Capubs20.o=polymerstats(Capubs20,lt_bs20);
Capubs20.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(Capubs20.o,Capubs20.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
Capubs20.kdim=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim(Capubs20.o,Capubs20.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
Capubs20.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(Capubs20.o,Capubs20.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
Capubs20.kdimc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim_calc(Capubs20.o,Capubs20.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

Capubs35.name='Capu';
Capubs35.pp_length_vec=[6;8;9;14;9];
Capubs35.p_length_vec=[6;7;9;14;9];
Capubs35.pp_index_vec=[33,47,61,78,96];
Capubs35.fh1_length=123;
Capubs35.o=polymerstats(Capubs35,lt_bs35);
Capubs35.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(Capubs35.o,Capubs35.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
Capubs35.kdim=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim(Capubs35.o,Capubs35.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
Capubs35.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(Capubs35.o,Capubs35.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
Capubs35.kdimc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim_calc(Capubs35.o,Capubs35.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);


%%
Dlength=5;
Blength=14;
Clength=7;

D18.name='BNI1P(pPD18)';
D18.pp_length_vec=Dlength;
D18.p_length_vec=Dlength-1;
D18.pp_index_vec=25;
D18.fh1_length=88;
D18.o=polymerstats(D18,lt);
D18.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(D18.o,D18.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
D18.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(D18.o,D18.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);


D28.name='BNI1P(pPD28)';
D28.pp_length_vec=Dlength;
D28.p_length_vec=Dlength-1;
D28.pp_index_vec=36;
D28.fh1_length=88;
D28.o=polymerstats(D28,lt);
D28.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(D28.o,D28.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
D28.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(D28.o,D28.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);


D37.name='BNI1P(pPD37)';
D37.pp_length_vec=Dlength;
D37.p_length_vec=Dlength-1;
D37.pp_index_vec=47;
D37.fh1_length=88;
D37.o=polymerstats(D37,lt);
D37.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(D37.o,D37.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
D37.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(D37.o,D37.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);


D65.name='BNI1P(pPD65)';
D65.pp_length_vec=Dlength;
D65.p_length_vec=Dlength-1;
D65.pp_index_vec=75;
D65.fh1_length=88;
D65.o=polymerstats(D65,lt);
D65.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(D65.o,D65.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
D65.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(D65.o,D65.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);


B37.name='BNI1P(pPB37)';
B37.pp_length_vec=Blength;
B37.pp_length_vec_max=10;
B37.p_length_vec=Blength-1;
B37.pp_index_vec=48;
B37.fh1_length=94;
B37.o=polymerstats(B37,lt);
B37.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(B37.o,B37.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
B37.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(B37.o,B37.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);


C37.name='BNI1P(pPC37)';
C37.pp_length_vec=Clength;
C37.p_length_vec=Clength-1;
C37.pp_index_vec=45;
C37.fh1_length=87;
C37.o=polymerstats(C37,lt);
C37.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(C37.o,C37.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
C37.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(C37.o,C37.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);


B18.name='BNI1P(pPB18)';
B18.pp_length_vec=Blength;
B18.pp_length_vec_max=10;
B18.p_length_vec=Blength-1;
B18.pp_index_vec=26;
B18.fh1_length=94;
B18.o=polymerstats(B18,lt);
B18.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(B18.o,B18.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
B18.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(B18.o,B18.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);


C18.name='BNI1P(pPC18)';
C18.pp_length_vec=Clength;
C18.p_length_vec=Clength-1;
C18.pp_index_vec=23;
C18.fh1_length=87;
C18.o=polymerstats(C18,lt);
C18.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(C18.o,C18.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
C18.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(C18.o,C18.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

%%
D18bs20.name='BNI1P(pPD18)';
D18bs20.pp_length_vec=Dlength;
D18bs20.p_length_vec=Dlength-1;
D18bs20.pp_index_vec=25;
D18bs20.fh1_length=90;
D18bs20.o=polymerstats(D18bs20,lt_bs20);
D18bs20.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(D18bs20.o,D18bs20.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
D18bs20.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(D18bs20.o,D18bs20.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);


D28bs20.name='BNI1P(pPD28)';
D28bs20.pp_length_vec=Dlength;
D28bs20.p_length_vec=Dlength-1;
D28bs20.pp_index_vec=36;
D28bs20.fh1_length=90;
D28bs20.o=polymerstats(D28bs20,lt_bs20);
D28bs20.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(D28bs20.o,D28bs20.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
D28bs20.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(D28bs20.o,D28bs20.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);


D37bs20.name='BNI1P(pPD37)';
D37bs20.pp_length_vec=Dlength;
D37bs20.p_length_vec=Dlength-1;
D37bs20.pp_index_vec=47;
D37bs20.fh1_length=90;
D37bs20.o=polymerstats(D37bs20,lt_bs20);
D37bs20.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(D37bs20.o,D37bs20.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
D37bs20.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(D37bs20.o,D37bs20.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);


D65bs20.name='BNI1P(pPD65)';
D65bs20.pp_length_vec=Dlength;
D65bs20.p_length_vec=Dlength-1;
D65bs20.pp_index_vec=75;
D65bs20.fh1_length=90;
D65bs20.o=polymerstats(D65bs20,lt_bs20);
D65bs20.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(D65bs20.o,D65bs20.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
D65bs20.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(D65bs20.o,D65bs20.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);


B37bs20.name='BNI1P(pPB37)';
B37bs20.pp_length_vec=Blength;
B37bs20.p_length_vec=Blength-1;
B37bs20.pp_index_vec=48;
B37bs20.fh1_length=90;
B37bs20.o=polymerstats(B37bs20,lt_bs20);
B37bs20.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(B37bs20.o,B37bs20.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
B37bs20.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(B37bs20.o,B37bs20.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);


C37bs20.name='BNI1P(pPC37)';
C37bs20.pp_length_vec=Clength;
C37bs20.p_length_vec=Clength-1;
C37bs20.pp_index_vec=45;
C37bs20.fh1_length=90;
C37bs20.o=polymerstats(C37bs20,lt_bs20);
C37bs20.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(C37bs20.o,C37bs20.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
C37bs20.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(C37bs20.o,C37bs20.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);


B18bs20.name='BNI1P(pPB18)';
B18bs20.pp_length_vec=Blength;
B18bs20.p_length_vec=Blength-1;
B18bs20.pp_index_vec=26;
B18bs20.fh1_length=90;
B18bs20.o=polymerstats(B18bs20,lt);
B18bs20.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(B18bs20.o,B18bs20.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
B18bs20.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(B18bs20.o,B18bs20.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);


C18bs20.name='BNI1P(pPC18)';
C18bs20.pp_length_vec=Clength;
C18bs20.p_length_vec=Clength-1;
C18bs20.pp_index_vec=23;
C18bs20.fh1_length=90;
C18bs20.o=polymerstats(C18bs20,lt);
C18bs20.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(C18bs20.o,C18bs20.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
C18bs20.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(C18bs20.o,C18bs20.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
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

SOSNTD_err=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        (log2(FHOD.kdim(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(FHOD.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)))),...
        (log2(Capu.kdim(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(Capu.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))))...
    ],ratios,ratios_err,ratios_err);

SOSNTD_err_bs20=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        (log2(FHODbs20.kdim(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(FHODbs20.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)))),...
        (log2(Capubs20.kdim(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(Capubs20.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))))...
    ],ratios,ratios_err,ratios_err);

SOSNTD_err_gate_bs35=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        (log2(FHODbs35.kdim(k_paf,k_pab*FHOD_gating,k_paf_rev,r_PF_rev,r_paf_rev)./(FHODbs35.kdob(k_paf,k_pab*FHOD_gating,k_paf_rev,r_PF_rev,r_paf_rev)))),...
        (log2(Capubs35.kdim(k_paf,k_pab*Capu_gating,k_paf_rev,r_PF_rev,r_paf_rev)./(Capubs35.kdob(k_paf,k_pab*Capu_gating,k_paf_rev,r_PF_rev,r_paf_rev))))...
    ],ratios,ratios_err,ratios_err);

SOSNTD_err_bs35=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        (log2(FHODbs35.kdim(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(FHODbs35.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)))),...
        (log2(Capubs35.kdim(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(Capubs35.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))))...
    ],ratios,ratios_err,ratios_err);

SOSNTD_err2=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        (log2(FHOD2.kdim(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(FHOD2.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)))),...
        (log2(Capu.kdim(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(Capu.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))))...
    ],ratios,ratios_err,ratios_err);

SOSNTD_err3=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        (log2(FHOD3.kdim(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(FHOD3.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)))),...
        (log2(Capu.kdim(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(Capu.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))))...
    ],ratios,ratios_err,ratios_err);

SOSNTD_err4=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        (log2(FHOD4.kdim(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(FHOD4.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)))),...
        (log2(Capu.kdim(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(Capu.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))))...
    ],ratios,ratios_err,ratios_err);

SOSNTD_err5=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        (log2(FHOD5.kdim(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(FHOD5.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)))),...
        (log2(Capu.kdim(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(Capu.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))))...
    ],ratios,ratios_err,ratios_err);

SOSNTD_err6=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        (log2(FHOD6.kdim(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(FHOD6.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)))),...
        (log2(Capu.kdim(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(Capu.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))))...
    ],ratios,ratios_err,ratios_err);

SOSNTD_err_cr=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        (log2(FHOD.kdimc(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(FHOD.kdobc(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)))),...
        (log2(Capu.kdimc(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(Capu.kdobc(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))))...
    ],ratios,ratios_err,ratios_err);
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
fig3_err_top=[0.46829044 0.72194778 0.8 0.33170573];		
fig3_err_bot=[0.46825234 0.68292357 0.8 0.33170574];

fig3sim=[D18,D28,D37,D65];
fig3ksim={fig3sim.kdob};
fig3ksimc={fig3sim.kdobc};

fig3sim_bs20=[D18bs20,D28bs20,D37bs20,D65bs20];
fig3ksim_bs20={fig3sim_bs20.kdob};
fig3ksimc_bs20={fig3sim_bs20.kdobc};
fig3ksimf_bs20 = @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) cellfun(@(f) f(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev), fig3ksim_bs20(:));
fig3ksimfc_bs20 = @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) cellfun(@(f) f(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev), fig3ksimc_bs20(:));
fig3ksimf_r_bs20 = @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) fig3ksimf_bs20(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(D65bs20.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)/2.471568);
fig3ksimf_noD37_bs20 = @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) cellfun(@(f) f(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev), fig3ksim_bs20([1 2 4]));
fig3ksimf_r_noD37_bs20 = @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) fig3ksimf_noD37_bs20(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(D65bs20.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)/2.471568);


fig3ksimfc = @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) cellfun(@(f) f(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev), fig3ksimc(:));

fig3ksimf = @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) cellfun(@(f) f(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev), fig3ksim(:));

fig3ksimf_noD37 = @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) cellfun(@(f) f(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev), fig3ksim([1 2 4]));

fig3ksimf_r = @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) fig3ksimf(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(D65.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)/2.471568);

SOS3= @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) sum(abs(([...
        fig3ksimf(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ]'-fig3k).^2));

SOS3_r= @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) sum(abs(([...
        fig3ksimf_r(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ]'-fig3k).^2));

SOS3_err=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        fig3ksimf(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ],fig3k,fig3_err_top,fig3_err_bot);

SOS3_err_bs20=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        fig3ksimf_bs20(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ],fig3k,fig3_err_top,fig3_err_bot);

SOS3_err_r_bs20=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        fig3ksimf_r_bs20(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ],fig3k,fig3_err_top,fig3_err_bot);

SOS3_err_noD37_bs20=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        fig3ksimf_noD37_bs20(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ],fig3k([1 2 4]),fig3_err_top([1 2 4]),fig3_err_bot([1 2 4]));

SOS3_err_r_noD37_bs20=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        fig3ksimf_r_noD37_bs20(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ],fig3k([1 2 4]),fig3_err_top([1 2 4]),fig3_err_bot([1 2 4]));

SOS3_err_r_noD37_gate_bs20=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        fig3ksimf_r_noD37_bs20(k_paf,k_pab*bni1_gating,k_paf_rev,r_PF_rev,r_paf_rev)...
    ],fig3k([1 2 4]),fig3_err_top([1 2 4]),fig3_err_bot([1 2 4]));


SOS3_err_noD37=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        fig3ksimf_noD37(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ],fig3k([1 2 4]),fig3_err_top([1 2 4]),fig3_err_bot([1 2 4]));

SOS3_errc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        fig3ksimfc(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ],fig3k,fig3_err_top,fig3_err_bot);

%% Fig 4 eq
fig4ak=[ 9.617512, 6.885804, 2.70152246];
fig4_err_top=[0.60709018 0.9105691 0.49864499];		
fig4_err_bot=[1.17073171 0.88888889 0.5203252];
fig4asim=[B37,C37,D37];
fig4aksim={fig4asim.kdob};
fig4aksimc={fig4asim.kdobc};

fig4aksimfc = @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) cellfun(@(f) f(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev), fig4aksimc(:));

fig4aksimf = @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) cellfun(@(f) f(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev), fig4aksim(:));

fig4aksimf_r = @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) fig4aksimf(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(D37.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)/3);

fig4asim_bs20=[B37bs20,C37bs20,D37bs20];
fig4aksim_bs20={fig4asim_bs20.kdob};
fig4aksimc_bs20={fig4asim_bs20.kdobc};
fig4aksimfc_bs20 = @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) cellfun(@(f) f(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev), fig4aksimc_bs20(:));
fig4aksimf_bs20 = @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) cellfun(@(f) f(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev), fig4aksim_bs20(:));
fig4aksimf_r_bs20 = @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) fig4aksimf_bs20(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(D37bs20.kdob(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)/3);


SOS4a= @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) sum(abs(([...
        fig4aksimf(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ]'-fig4ak).^2));

SOS4a_r= @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) sum(abs(([...
        fig4aksimf_r(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ]'-fig4ak).^2));

SOS4a_err=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        fig4aksimf(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ],fig4ak,fig4_err_top,fig4_err_bot);

SOS4a_err_bs20=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        fig4aksimf_bs20(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ],fig4ak,fig4_err_top,fig4_err_bot);

SOS4a_err_r_bs20=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        fig4aksimf_r_bs20(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ],fig4ak,fig4_err_top,fig4_err_bot);

SOS4a_err_r_gate_bs20=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        fig4aksimf_r_bs20(k_paf,k_pab*bni1_gating,k_paf_rev,r_PF_rev,r_paf_rev)...
    ],fig4ak,fig4_err_top,fig4_err_bot);

SOS4a_errc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        fig4aksimfc(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ],fig4ak,fig4_err_top,fig4_err_bot);

%% Full sos
%SOSs={SOS4a,SOS3,SOSNTD};
SOS=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) SOSNTD(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)+...
    SOS3(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)+...
    SOS4a(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

SOS_r=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) SOSNTD(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)+...
    SOS3_r(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)+...
    SOS4a_r(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

SOS_err=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) (10*SOSNTD_err(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))+...
    SOS3_err(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)+...
    SOS4a_err(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

SOS_err2=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) (10*SOSNTD_err2(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))+...
    SOS3_err(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)+...
    SOS4a_err(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

SOS_err_noF3D37=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) (10*SOSNTD_err(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))+...
    SOS3_err_noD37(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)+...
    SOS4a_err(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

SOS_errc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) (10*SOSNTD_err_cr(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))+...
    SOS3_errc(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)+...
    SOS4a_errc(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

SOS_err_bs20=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) (10*SOSNTD_err_bs20(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))+...
    SOS3_err_bs20(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)+...
    SOS4a_err_bs20(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

SOS_err_bs35=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) (10*SOSNTD_err_bs35(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))+...
    SOS3_err_bs20(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)+...
    SOS4a_err_bs20(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

SOS_err_r_bs35=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) (100*SOSNTD_err_bs35(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))+...
    SOS3_err_r_bs20(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)+...
    SOS4a_err_r_bs20(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

SOS_err_r_noF3D37_bs35=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) (100*SOSNTD_err_bs35(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))+...
    SOS3_err_r_noD37_bs20(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)+...
    SOS4a_err_r_bs20(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

SOS_err_r_noF3D37_gate_bs35=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) (100*SOSNTD_err_gate_bs35(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))+...
    SOS3_err_r_noD37_bs20(k_paf,k_pab*bni1_gating,k_paf_rev,r_PF_rev,r_paf_rev)+...
    SOS4a_err_r_bs20(k_paf,k_pab*bni1_gating,k_paf_rev,r_PF_rev,r_paf_rev);

SOS_err_noF3D37_gate_bs35=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) (100*SOSNTD_err_gate_bs35(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))+...
    SOS3_err_noD37_bs20(k_paf,k_pab*bni1_gating,k_paf_rev,r_PF_rev,r_paf_rev)+...
    SOS4a_err_bs20(k_paf,k_pab*bni1_gating,k_paf_rev,r_PF_rev,r_paf_rev);

SOS_err_noF3D37_bs35=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) (100*SOSNTD_err_bs35(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))+...
    SOS3_err_noD37_bs20(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)+...
    SOS4a_err_r_bs20(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

SOS_err_r_bs20=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) (10*SOSNTD_err_bs20(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev))+...
    SOS3_err_r_bs20(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)+...
    SOS4a_err_r_bs20(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

options = optimset('MaxFunEvals',500000,'MaxIter',50000,'TolX',10^(-10),'TolFun',10^(-10),'Display','Iter'); 
options = optimset('MaxFunEvals',60000,'MaxIter',60000,'TolX',10^(-6),'TolFun',10^(-6)); 
%options = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolX',10^(-15),'TolFun',10^(-15),'PlotFcns','optimplotfval','Display','iter','TypicalX',[25,10^(1),10^(15),10^(10),10^(10)]); 
%options=optimoptions('fsolve', 'Display','iter','ScaleProblem','none','InitDamping',0.001,'FunctionTolerance',10^(-9),'StepTolerance',10^(-10));

%[params,fval] = fminsearch(@(z) SOS(abs(z(1)),abs(z(2)),abs(z(3)),abs(z(4)),abs(z(5))),[25,1,32704292,1250000,50000],options);
%[params,fval] = fminsearch(@(z) SOS_r(abs(z(1)),abs(z(2)),abs(z(3)),abs(z(4)),abs(z(5))),[457,10,32704292,1250000,50000],options);
[params_FHOD,fval_FHOD] = fminsearch(@(z) SOSFHOD(abs(z(1)),abs(z(2)),abs(z(3)),abs(z(4)),abs(z(5))),[25,1,32704292,1250000,50000],options);
[params,fval]=fminsearch(@(z) SOSFHOD(abs(z(1)),abs(z(2)),abs(z(3)),abs(z(4)),abs(z(5))),[25,1,32704292,1250000,50000],options)
[params,fval]=fminsearch(@(z) SOS_err(abs(z(1)),abs(z(2)),abs(z(3)),abs(z(4)),abs(z(5))),[25,1,32704292,1250000,50000],options)


for i=1:1
    disp(i)
    r=15*rand(1,5)+-3;
    % Not ratio
        [params,fval] = fminsearch(@(z) SOS(abs(z(1)),abs(z(2)),abs(z(3)),abs(z(4)),abs(z(5))),[10^r(1);10^r(2);10^r(3);10^r(4);10^r(5)],options); 
        vals=num2cell(abs(params));
        T=[T;[{fval},vals(:)']];
    % ratio
        [params,fval] = fminsearch(@(z) SOS_r(abs(z(1)),abs(z(2)),abs(z(3)),abs(z(4)),abs(z(5))),[10^r(1);10^r(2);10^r(3);10^r(4);10^r(5)],options); 
        vals=num2cell(abs(params));
        T=[T;[{fval}, SOS(vals{:}),SOSNTD_err(vals{:}),SOS3_err(vals{:}),SOS4a_err(vals{:}), vals(:)']];
end

parfor i=1:1000
    r=6*rand(1,5)+-2;
        [params,fval] = fminsearch(@(z) SOS_err_r_noF3D37_gate_bs35(abs(z(1)),abs(z(2)),abs(z(3)),abs(z(4)),abs(z(5))),[10^r(1);10^r(2);10^r(3);10^r(4);10^r(5)],options);
        disp(fval)
        vals=num2cell(abs(params));
        T_err_r_gate_F3D37_bs35(i,:)=[{fval}, SOS_err_noF3D37_gate_bs35(vals{:}),SOSNTD_err_gate_bs35(vals{:}),SOS3_err_r_noD37_gate_bs20(vals{:}),SOS4a_err_r_gate_bs20(vals{:}), vals(:)'];
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