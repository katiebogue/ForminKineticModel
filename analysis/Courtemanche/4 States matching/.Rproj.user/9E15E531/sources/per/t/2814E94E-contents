addkpoly=function(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev,des){
  n=mnum+1
  mod=paste("Model",n,sep=" ")
  
  new_kpoly=kpoly(k_paf,k_pab,k_paf_rev,r_PF_rev,r_paf_rev)
  
  MData<<-addcol(new_kpoly,mod)
  print(mod)
  
  plot3(mod)
  
  plot4(mod)
  
  x=as.data.frame(des)
  row.names(x)=mod
  names(x)=names(Models)
  Models<<-rbind(Models,x)
  
  mnum<<-n
  
}

