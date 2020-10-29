#Threshold GAM function
th.gam<-
  function(formula_tgam,formula_base,th.name,a,b,n_incr,data,th.criteria='aic',...){
    data=data
    threshold.variable<-eval(parse(text=paste("data", th.name,sep = "$")))
    low<-quantile(threshold.variable,a)
    high<-quantile(threshold.variable,b)
    th.explore<-seq(low,high,length.out=n_incr)
    aic.all<-NA*(th.explore)
    gcv.all<-NA*(th.explore)
    gam.all<-as.list(1:(length(th.explore)))
    for(i in 1:length(th.explore)){
      data$th<-ifelse(threshold.variable<=th.explore[i],'below_th','above_th')
      gam.all[[i]]<-gam(formula_tgam,data=data,...)
      aic.all[i]<-AIC(gam.all[[i]])
      gcv.all[i]<-gam.all[[i]]$gcv.ubre
    }
    gam.base<-gam(formula_base,data=data,...)
    index<-ifelse(th.criteria=='aic',order(aic.all)[1],order(gcv.all)[1])
    th.value<-th.explore[index]
    best.model<-gam.all[[index]]
    best.aic<-sort(aic.all)[1]
    aic.base<-AIC(gam.base)
    gcv.base<-gam.base$gcv.ubre
    invisible(list(gam.all=gam.all,aic.all=aic.all,
                   gcv.all=gcv.all,th.explore=th.explore,th.value=th.value,best.aic=best.aic,
                   best.model=best.model,gam.base=gam.base,aic.base=aic.base,gcv.base=gcv.base))
  }
