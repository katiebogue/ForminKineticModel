plot3=function(x){
myplot=barplot(t(cbind(MData[c('BNI1P(pPD18)','BNI1P(pPD28)','BNI1P(pPD37)','BNI1P(pPD65)'),'courtvals'],MData[c('BNI1P(pPD18)','BNI1P(pPD28)','BNI1P(pPD37)','BNI1P(pPD65)'),names(MData)==x]/(MData['BNI1P(pPD65)',names(MData)==x]/MData['BNI1P(pPD65)','courtvals']))),beside = TRUE,ylim=c(0,11),col=c("aquamarine3","coral"),ylab="Scaled Kpoly",names.arg = c("pPD18","pPD28","pPD37","pPD65"),main=x,sub=wrap_strings(Models[row.names(Models)==x,'Description'],100),cex.sub=0.7)
arrows(x0=myplot[1,],y0=MData[c('BNI1P(pPD18)','BNI1P(pPD28)','BNI1P(pPD37)','BNI1P(pPD65)'),'courtvals']+(Fig3Er*MData[c('BNI1P(pPD18)','BNI1P(pPD28)','BNI1P(pPD37)','BNI1P(pPD65)'),'courtvals']),y1=MData[c('BNI1P(pPD18)','BNI1P(pPD28)','BNI1P(pPD37)','BNI1P(pPD65)'),'courtvals']-(Fig3Er*MData[c('BNI1P(pPD18)','BNI1P(pPD28)','BNI1P(pPD37)','BNI1P(pPD65)'),'courtvals']),angle = 90,code = 3,length = 0.1)
legend("topright",legend=c("Experiment","Model"),col=c("aquamarine3","coral"),pch=15)
}

plot4=function(i){
  prms=c("BNI1P(pPB37)","BNI1P(pPC37)","BNI1P(pPD37)","BNI1P(pPB18)","BNI1P(pPC18)","BNI1P(pPD18)")
  x37=(MData[prms[3],names(MData)==i])/MData[prms[3],'courtvals']
  x18=(MData[prms[4],names(MData)==i])/MData[prms[4],'courtvals']
  myplot=barplot(
    t(cbind(
      MData[prms,'courtvals'],
      c((MData[prms[1:3],names(MData)==i]/x37),
        (MData[prms[4:6],names(MData)==i]/x18))
      )
      ),
    beside = TRUE,
    col=c("aquamarine3","coral"),
    ylab="Scaled Kpoly",main=i,
    ylim=c(0,13),
    names.arg=c("pPB37","pPC37","pPD37","pPB18","pPC18","pPD18"),
    axes=FALSE,
    sub=wrap_strings(Models[row.names(Models)==i,'Description'],100),cex.sub=0.7)
  arrows(x0=myplot[1,],y0=MData[prms,'courtvals']+(Fig4Er*MData[prms,'courtvals']),y1=MData[prms,'courtvals']-(Fig4Er*MData[prms,'courtvals']),angle = 90,code = 3,length = 0.1)
  legend("topright",legend=c("Experiment","Model"),col=c("aquamarine3","coral"),pch=15,inset=0.02)
  axis(2, seq(0, 12, by=2), labels=TRUE, tck=-.03)
  print(x37)
  print(x18)
}
