kpoly=function(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev){
  kcapa=k_paf*MData[,'c_PA']*(1-MData[,'all_p_occ2a'])
  kcapb=k_paf*MData[,'c_PA']*(1-MData[,'all_p_occ2b'])
  kdela=k_pab*(1-MData[,'all_p_occ2a_0'])*(1.0e33*MData[,'all_p_r2a']/(27*6.022e23))
  kdelb=k_pab*(1-MData[,'all_p_occ2b_0'])*(1.0e33*MData[,'all_p_r2b']/(27*6.022e23))
  kpa=(1/r_PF_rev) + ((r_paf_rev + r_PF_rev)/(kdela * r_PF_rev)) + (((k_paf_rev * r_paf_rev) + (k_paf_rev * r_PF_rev) + (kdela * r_PF_rev))/(kcapa * kdela * r_PF_rev))
  kpb=(1/r_PF_rev) + ((r_paf_rev + r_PF_rev)/(kdelb * r_PF_rev)) + (((k_paf_rev * r_paf_rev) + (k_paf_rev * r_PF_rev) + (kdelb * r_PF_rev))/(kcapb * kdelb * r_PF_rev))
  kpa=1/kpa
  kpb=1/kpb
  value=kpa+kpb
  return(value)
}
