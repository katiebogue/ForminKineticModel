% generates workspace to be used in bruteforce.m
load("nofunctions.mat") % file with workspace

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


%%
ratios=[-0.6172640408,0];
ratios_err=[0.2584555583 0.2584555583];

FHOD.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(FHOD.o,FHOD.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD.kdim=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim(FHOD.o,FHOD.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(FHOD.o,FHOD.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD.kdimc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim_calc(FHOD.o,FHOD.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);


FHOD2.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(FHOD2.o,FHOD2.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD2.kdim=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim(FHOD2.o,FHOD2.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD2.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(FHOD2.o,FHOD2.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD2.kdimc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim_calc(FHOD2.o,FHOD2.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);



FHOD3.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(FHOD3.o,FHOD3.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD3.kdim=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim(FHOD3.o,FHOD3.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD3.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(FHOD3.o,FHOD3.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD3.kdimc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim_calc(FHOD3.o,FHOD3.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);




FHOD4.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(FHOD4.o,FHOD4.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD4.kdim=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim(FHOD4.o,FHOD4.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD4.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(FHOD4.o,FHOD4.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD4.kdimc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim_calc(FHOD4.o,FHOD4.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);



FHOD5.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(FHOD5.o,FHOD5.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD5.kdim=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim(FHOD5.o,FHOD5.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD5.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(FHOD5.o,FHOD5.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD5.kdimc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim_calc(FHOD5.o,FHOD5.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);



FHOD6.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(FHOD6.o,FHOD6.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD6.kdim=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim(FHOD6.o,FHOD6.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD6.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(FHOD6.o,FHOD6.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD6.kdimc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim_calc(FHOD6.o,FHOD6.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);



FHOD7.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(FHOD7.o,FHOD7.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD7.kdim=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim(FHOD7.o,FHOD7.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD7.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(FHOD7.o,FHOD7.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
FHOD7.kdimc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim_calc(FHOD7.o,FHOD7.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);



Capu.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(Capu.o,Capu.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
Capu.kdim=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim(Capu.o,Capu.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
Capu.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(Capu.o,Capu.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
Capu.kdimc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdim_calc(Capu.o,Capu.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

%%
D18.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(D18.o,D18.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
D18.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(D18.o,D18.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);



D28.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(D28.o,D28.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
D28.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(D28.o,D28.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);



D37.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(D37.o,D37.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
D37.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(D37.o,D37.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);



D65.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(D65.o,D65.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
D65.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(D65.o,D65.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);



B37.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(B37.o,B37.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
B37.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(B37.o,B37.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);



C37.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(C37.o,C37.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
C37.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(C37.o,C37.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);



B18.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(B18.o,B18.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
B18.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(B18.o,B18.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);



C18.kdob=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb(C18.o,C18.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);
C18.kdobc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) kdb_calc(C18.o,C18.pp_length_vec,k_paf,c_PA,k_pab,k_paf_rev,r_PF_rev,r_paf_rev);

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

%% Fig 3 eq
fig3k=[ 9.554461, 7.642275, 3.000000, 2.471568];
fig3_err_top=[0.46829044 0.72194778 0.8 0.33170573];		
fig3_err_bot=[0.46825234 0.68292357 0.8 0.33170574];
fig3sim=[D18,D28,D37,D65];
fig3ksim={fig3sim.kdob};
fig3ksimc={fig3sim.kdobc};

fig3ksimfc = @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) cellfun(@(f) f(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev), fig3ksimc(:));

fig3ksimf = @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) cellfun(@(f) f(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev), fig3ksim(:));

fig3ksimf_noD37 = @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) cellfun(@(f) f(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev), fig3ksim([1 2 4]));

fig3ksimf_r = @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) fig3ksimf(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(D37.kdob(1,1,1,1,1)/3);

SOS3= @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) sum(abs(([...
        fig3ksimf(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ]'-fig3k).^2));

SOS3_r= @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) sum(abs(([...
        fig3ksimf_r(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ]'-fig3k).^2));

SOS3_err=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        fig3ksimf(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ],fig3k,fig3_err_top,fig3_err_bot);

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

fig4aksimf_r = @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) fig4aksimf(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)./(D37.kdob(1,1,1,1,1)/3);


SOS4a= @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) sum(abs(([...
        fig4aksimf(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ]'-fig4ak).^2));

SOS4a_r= @(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) sum(abs(([...
        fig4aksimf_r(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ]'-fig4ak).^2));

SOS4a_err=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        fig4aksimf(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ],fig4ak,fig4_err_top,fig4_err_bot);

SOS4a_errc=@(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev) errSOS([...
        fig4aksimfc(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)...
    ],fig4ak,fig4_err_top,fig4_err_bot);

%% Full sos
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

save("bruteforce.mat")
